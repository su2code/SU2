/*!
 * \file CFlameletSolver.cpp
 * \brief Main subroutines for the flamelet model solver.
 * \author D. Mayer, T. Economon
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

#include "../../include/solvers/CFlameletSolver.hpp"
#include "../../include/variables/CFlameletVariable.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"

CFlameletSolver::CFlameletSolver(void) : CScalarSolver() {
  
  Inlet_ScalarVars = NULL;
  
}

CFlameletSolver::CFlameletSolver(CGeometry *geometry,
                                 CConfig *config,
                                 unsigned short iMesh)
: CScalarSolver(geometry, config) {
  
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf;
  
  bool turbulent = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));

  bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));


  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarSolver
   will have variable numbers of equations. ---*/
 
  nVar     = config->GetNScalars();
  nPrimVar = config->GetNScalars();

  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  
  /*--- Fluid model pointer initialization ---*/
  
  FluidModel = NULL;
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual    [iVar] = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i  [iVar] = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j  [iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Res_Conv     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv    [iVar] = 0.0;
  Res_Visc     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc    [iVar] = 0.0;
  
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliar vector related with the solution ---*/
  
  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
  
  /*--- Define some auxiliar vector related with the geometry ---*/
  
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Define some auxiliar vector related with the flow solution ---*/
  
  FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of the structure of the whole Jacobian ---*/
  
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Combustion Scalar). MG level: " << iMesh <<"." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*transpose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = -1.0e15;
  upperlimit[0] =  1.0e15;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  
  /*--- Set up fluid model for the diffusivity ---*/
  
  su2double Diffusivity_Ref = 1.0;
  
  su2double DiffusivityND = config->GetDiffusivity_Constant()/Diffusivity_Ref;
  
  config->SetDiffusivity_ConstantND(DiffusivityND);
  
  /*--- Scalar variable state at the far-field. ---*/
  
  Scalar_Inf = new su2double[nVar];
  
  //FIXME daniel: scalar_inf should be set depending on inlet temperature
  //              can be done in the initial condition, also set Inlet_Scalar_vars!
  for (iVar = 0; iVar < nVar; iVar++){
    Scalar_Inf[iVar] = config->GetScalar_Init(iVar);
  }
  /*--- Initialize the solution to the far-field state everywhere. ---*/
  
  nodes = new CFlameletVariable(Scalar_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Initialize quantities for SlidingMesh Interface ---*/
  
  unsigned long iMarker;
  
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
  
  /*-- Allocation of inlets has to happen in derived classes
   (not CScalarSolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  
  Inlet_Position = nDim*2+2;
  if (turbulent) {
    if (turb_SA) Inlet_Position += 1;
    else if (turb_SST) Inlet_Position += 2;
  }
  
  Inlet_ScalarVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_ScalarVars[iMarker] = new su2double*[nVertex[iMarker]];
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_ScalarVars[iMarker][iVertex] = new su2double[nVar];
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        Inlet_ScalarVars[iMarker][iVertex][0] = Scalar_Inf[iVar];
    }
  }
  
}

CFlameletSolver::~CFlameletSolver(void) {
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if ( SlidingState != NULL ) {
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
  
}

void CFlameletSolver::Preprocessing(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CConfig *config,
                                            unsigned short iMesh,
                                            unsigned short iRKStep,
                                            unsigned short RunTime_EqSystem,
                                            bool Output) {
  
  unsigned long ErrorCounter = 0;
  unsigned long InnerIter = config->GetInnerIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
                           (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_scalar   = ((config->GetKind_SlopeLimit_Scalar() != NO_LIMITER) &&
                           (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (!disc_adjoint) Jacobian.SetValZero();
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Upwind second order reconstruction ---*/
  
  if (limiter_scalar) SetSolution_Limiter(geometry, config);
  
  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
  
}

void CFlameletSolver::Postprocessing(CGeometry *geometry,
                                             CSolver **solver_container,
                                             CConfig *config,
                                             unsigned short iMesh) { }

unsigned long CFlameletSolver::SetPrimitive_Variables(CSolver **solver_container,
                                                              CConfig *config,
                                                              bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double Density, Temperature, lam_visc = 0.0, eddy_visc = 0.0, Cp;
  unsigned short turb_model = config->GetKind_Turb_Model();
  unsigned long n_not_in_domain = 0;
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the density, temperature, Cp, and laminar viscosity. ---*/
    
    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    Cp          = solver_container[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint);
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    lam_visc    = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
    }

    /*--- Compute and store the mass diffusivity. ---*/

    CFluidModel * fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
    su2double * scalars = nodes->GetSolution(iPoint);

    n_not_in_domain += fluid_model_local->SetTDState_T(Temperature,scalars);

    for(int i_scalar=0; i_scalar < config->GetNScalars(); ++i_scalar){
      nodes->SetDiffusivity(iPoint, fluid_model_local->GetMassDiffusivity(), i_scalar);
    }

    

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  SetNTableMisses(n_not_in_domain);

  return ErrorCounter;
  
}

void CFlameletSolver::SetInitialCondition(CGeometry **geometry,
                                                  CSolver ***solver_container,
                                                  CConfig *config,
                                                  unsigned long ExtIter) {
  su2double *coords;
  bool Restart   = (config->GetRestart() || config->GetRestart_Flow());
  
  
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE)
      cout << "Initializing progress variable and temperature (initial condition)." << endl;

    su2double *scalar_init  = new su2double[nVar];
    su2double *flame_offset = config->GetFlameOffset();
    su2double *flame_normal = config->GetFlameNormal();

    su2double prog_unburnt    = 0.0;
    su2double prog_burnt      = 0.2; // 2.09379237374982e-01;
    su2double flame_thickness = config->GetFlameThickness();
    su2double burnt_thickness = config->GetBurntThickness();
    su2double flamenorm       = sqrt( flame_normal[0]*flame_normal[0]
                                     +flame_normal[1]*flame_normal[1]
                                     +flame_normal[2]*flame_normal[2]);

    su2double temp_inlet = 300.;
    su2double prog_inlet = 0.0;

    su2double enth_inlet;
    su2double point_loc;
    su2double dist;
    unsigned long n_not_iterated  = 0;
    unsigned long n_not_in_domain = 0;

    CFluidModel *fluid_model_local;

    vector<string>     look_up_tags;
    vector<su2double*> look_up_data;
    string name_enth = config->GetScalarName(I_ENTHALPY);
    string name_prog = config->GetScalarName(I_PROG_VAR);
    

    for (unsigned long i_mesh = 0; i_mesh <= config->GetnMGLevels(); i_mesh++) {

      fluid_model_local = solver_container[i_mesh][FLOW_SOL]->GetFluidModel();

      prog_burnt = fluid_model_local->GetLookUpTable()->GetTableLimitsProg().second;

      for (unsigned long i_point = 0; i_point < geometry[i_mesh]->GetnPoint(); i_point++) {
        
        for (unsigned long i_var = 0; i_var < nVar; i_var++)
          Solution[i_var] = 0.0;
        
        coords = geometry[i_mesh]->nodes->GetCoord(i_point);

        /* determine if our location is above or below the plane, assuming the normal 
           is pointing towards the burned region*/ 
        point_loc = flame_normal[0]*(coords[0]-flame_offset[0]) 
                  + flame_normal[1]*(coords[1]-flame_offset[1]) 
                  + flame_normal[2]*(coords[2]-flame_offset[2]);

        /* compute the exact distance from point to plane */          
        point_loc = point_loc/flamenorm;

        if (point_loc <= 0){ /* unburnt region */
          scalar_init[I_PROG_VAR] = prog_unburnt;

        } else if ( (point_loc > 0) && (point_loc <= flame_thickness) ){ /* flame zone */
          scalar_init[I_PROG_VAR] = 0.5*(prog_unburnt + prog_burnt);

        } else if ( (point_loc > flame_thickness) && (point_loc <= flame_thickness + burnt_thickness) ){ /* burnt region */
          scalar_init[I_PROG_VAR] = prog_burnt;

        } else { /* unburnt region */
          scalar_init[I_PROG_VAR] = prog_unburnt;
        }

        n_not_iterated         += fluid_model_local->GetEnthFromTemp(&enth_inlet, prog_inlet,temp_inlet);
        scalar_init[I_ENTHALPY] = enth_inlet;

        n_not_in_domain        += fluid_model_local->GetLookUpTable()->LookUp_ProgEnth(look_up_tags, look_up_data, scalar_init[I_PROG_VAR], scalar_init[I_ENTHALPY],name_prog,name_enth);

        // nijso: 
        // skip progress variable and enthalpy
        // we can make an init based on the lookup table. 
        for(int i_scalar = 0; i_scalar < config->GetNScalars(); ++i_scalar){
          if ( (i_scalar != I_ENTHALPY) && (i_scalar != I_PROG_VAR) )
            scalar_init[i_scalar] = config->GetScalar_Init(i_scalar);
        }

        
        
        solver_container[i_mesh][SCALAR_SOL]->GetNodes()->SetSolution(i_point, scalar_init);

      }

      solver_container[i_mesh][SCALAR_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][SCALAR_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][FLOW_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->Preprocessing( geometry[i_mesh], solver_container[i_mesh], config, i_mesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      
    }

    delete[] scalar_init;

    if (rank == MASTER_NODE && (n_not_in_domain > 0 || n_not_iterated > 0))
      cout << endl;
    
    if (rank == MASTER_NODE && n_not_in_domain > 0)
      cout << " !!! Initial condition: Number of points outside of table domain: " << n_not_in_domain << " !!!" << endl;

    if (rank == MASTER_NODE && n_not_iterated > 0)
      cout << " !!! Initial condition: Number of points in which enthalpy could not be iterated: " << n_not_iterated << " !!!" << endl;

    if (rank == MASTER_NODE && (n_not_in_domain > 0 || n_not_iterated > 0))
      cout << endl;

  }
}

void CFlameletSolver::SetPreconditioner(CGeometry *geometry, CSolver **solver_container,  CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  su2double  BetaInc2, Density, dRhodT, dRhodC, Temperature, Cp, Delta;
  
  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Access the primitive variables at this node. ---*/
    
    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    BetaInc2    = solver_container[FLOW_SOL]->GetNodes()->GetBetaInc2(iPoint);
    Cp          = solver_container[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint);
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    
    /*--- We need the derivative of the equation of state to build the
     preconditioning matrix. For now, the only option is the ideal gas
     law, but in the future, dRhodT should be in the fluid model. ---*/
    
    if (variable_density) {
      dRhodT = -Density/Temperature;
    } else {
      dRhodT = 0.0;
    }
    
    /*--- Passive scalars have no impact on the density. ---*/
    
    dRhodC = 0.0;
    
    /*--- Modify matrix diagonal with term including volume and time step. ---*/
    
    su2double Vol = geometry->nodes->GetVolume(iPoint);
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*
                   solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    
    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/
    
    if (implicit) {
      
      for (int i_var = 0; i_var < nVar; i_var++) {
        for (int j_var = 0; j_var < nVar; j_var++) {
          Jacobian_i[i_var][j_var] = 0.0;
        }
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        total_index = iPoint*nVar+iVar;
        
        su2double scalar = nodes->GetSolution(iPoint, iVar);
        
        /*--- Add the extra Jacobian term to the scalar system. ---*/
        
        su2double Jaccomp = scalar * dRhodC + Density;
        
        Jacobian_i[iVar][iVar] = Jaccomp*Delta;
        
      }
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
}

void CFlameletSolver::Source_Residual(CGeometry       *geometry,
                                      CSolver        **solver_container,
                                      CNumerics      **numerics_container,
                                      CConfig         *config,
                                      unsigned short   iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool implicit            = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool axisymmetric        = config->GetAxisymmetric();
  bool viscous             = config->GetViscous();
  unsigned short n_scalars = config->GetNScalars();
  unsigned short n_lookups = config->GetNLookups();

  /*--- Combustion source term implementation. (active by default)  ---*/
  su2double delta_enth_lut;
  su2double delta_temp_lut;
  su2double delta_source_prog_lut;
  su2double delta_source_energy_lut;
  
  su2double  temperature_dummy = 73;
  su2double* scalars   = new su2double[n_scalars];
  su2double* sources   = new su2double[n_scalars];

  //pair<su2double, su2double> limits_lut_prog;
  
  CNumerics *second_numerics = numerics_container[SOURCE_SECOND_TERM];

  CNumerics *first_numerics = numerics_container[SOURCE_FIRST_TERM];

  CFluidModel *fluid_model_local;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Get volume of the dual cell. ---*/
    
    su2double Volume = geometry->nodes->GetVolume(iPoint);
    
    /*--- Compute the production term. ---*/
    for (iVar = 0; iVar < nVar; iVar ++)
      scalars[iVar] = nodes->GetSolution(iPoint, iVar);
    
    fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
    
    /*--- clip progress variable to table limits ---*/
    
    // limits_lut_prog  = fluid_model_local->GetTableLimitsProg();
    fluid_model_local->SetScalarSources(scalars);

    /*--- Implicit part for production term (to do). ---*/

    for(int i_lookup = 0; i_lookup < n_lookups; ++i_lookup){
      su2double next_lookup = fluid_model_local->GetLookupScalar(i_lookup);
      nodes->SetLookupScalar(iPoint, next_lookup, i_lookup);
    }

    for(int i_scalar = 0; i_scalar < n_scalars; ++i_scalar){
      su2double next_source = fluid_model_local->GetSourceScalar(i_scalar);
      nodes->SetSourceScalar(iPoint, next_source, i_scalar);
      Residual[i_scalar] = next_source * Volume;
    }



//  /*--- Implicit part for production term (to do). ---*/
//     for (int i_var = 0; i_var < nVar; i_var++) {
//       for (int j_var = 0; j_var < nVar; j_var++) {
//         Jacobian_i[i_var][j_var] = 0.0;
//       }
//     }
//     Jacobian_i[0][0] = Volume*fluid_model_local->GetdSourcePVdPV();
    
    
//     /*--- Add Residual ---*/
    
//     LinSysRes.SubtractBlock(iPoint, Residual);
    
//     /*--- Implicit part ---*/
    
//     if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);



    first_numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);

    auto residual = first_numerics->ComputeResidual(config);

    /*--- Add Residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, residual);
    
    /*--- Implicit part ---*/
    
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
    
  }
  
  /*--- Axisymmetry source term for the scalar equation. ---*/
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);
      
      /*--- Scalar variables w/o reconstruction ---*/
      
      second_numerics->SetScalarVar(nodes->GetSolution(iPoint), NULL);
      
      /*--- Mass diffusivity coefficients. ---*/
      
      second_numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                                         NULL);
      
      /*--- Set control volume ---*/
      
      second_numerics->SetVolume(geometry->nodes->GetVolume(iPoint));
      
      /*--- Set y coordinate ---*/
      
      second_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                NULL);
      
      /*--- If viscous, we need gradients for extra terms. ---*/
      
      if (viscous) {
        
        /*--- Gradient of the scalar variables ---*/
        
        second_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), NULL);
        
      }
      
      /*--- Compute Source term Residual ---*/
      
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }

  delete[] scalars;
  delete[] sources;
  
}

void CFlameletSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config,
                               unsigned short val_marker) {

  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint, total_index, not_used;
  su2double *Coords;
  su2double enth_inlet;

  bool        grid_movement = config->GetGrid_Movement      (          );
  string      Marker_Tag    = config->GetMarker_All_TagBound(val_marker);
  su2double   temp_inlet    = config->GetInlet_Ttotal       (Marker_Tag);
  su2double  *inlet_scalar  = config->GetInlet_ScalarVal    (Marker_Tag);
 
  CFluidModel  *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();

  not_used                 = fluid_model_local->GetEnthFromTemp(&enth_inlet, inlet_scalar[I_PROG_VAR], temp_inlet);
  inlet_scalar[I_ENTHALPY] = enth_inlet;

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

   /* Dirichlet boundary condition at the inlet for scalars */

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

       if (geometry->nodes->GetDomain(iPoint)) {
    
        nodes->SetSolution_Old(iPoint, inlet_scalar);

        LinSysRes.SetBlock_Zero(iPoint);

        for (iVar = 0; iVar < nVar; iVar++) {
          nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
        }

        /*--- Includes 1 in the diagonal ---*/
         for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }      
    }
  }
}

void CFlameletSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker) {

  unsigned long i_point, i_vertex, total_index;
  su2double *normal;
  su2double *v_outlet, *v_domain;

  bool grid_movement  = config->GetGrid_Movement();
  bool implicit       = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  

  CFluidModel *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();  

  normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (i_vertex = 0; i_vertex < geometry->nVertex[val_marker]; i_vertex++) {
    
    i_point = geometry->vertex[val_marker][i_vertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(i_point)) {

        /*--- get the normal vector for i_vertex ---*/
        geometry->vertex[val_marker][i_vertex]->GetNormal(normal);

        /* negate sign for outward facing normal */
        for (int i_dim = 0; i_dim < nDim; ++i_dim)
          normal[i_dim] *= -1;

        /*--- set the normal ---*/
        conv_numerics->SetNormal(normal);

        /*--- get solution at the farfield boundary node ---*/
        v_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(i_point);

        /*--- get solution on the outlet ---*/
        v_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, i_vertex);

        /* --- Neumann condition for velocity ---*/
        for (int i_dim = 0; i_dim < nDim; ++i_dim)
            v_outlet[i_dim+1] = v_domain[i_dim+1];

        /* --- Neumann condition for density ---*/
        v_outlet[nDim+2] = v_domain[nDim+2];

        conv_numerics->SetPrimitive(v_domain, v_outlet);

        conv_numerics->SetScalarVar(nodes->GetSolution(i_point),
                                    nodes->GetSolution(i_point));

        if (dynamic_grid)
          conv_numerics->SetGridVel(geometry->nodes->GetGridVel(i_point),
                                    geometry->nodes->GetGridVel(i_point));

        /*--- compute the residual using an upwind scheme ---*/
        auto residual = conv_numerics->ComputeResidual(config);

        /*--- Update residual value ---*/
        LinSysRes.AddBlock(i_point, residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.AddBlock2Diag(i_point, residual.jacobian_i);
    }
  }
  delete [] normal;
}

void CFlameletSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics *conv_numerics,
                                         CNumerics *visc_numerics,
                                         CConfig *config,
                                         unsigned short val_marker) {

  unsigned short iVar, jVar, iDim;
  unsigned long iVertex, iPoint, total_index;

  bool implicit                   = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  string Marker_Tag               = config->GetMarker_All_TagBound(val_marker);
  su2double temp_wall             = config->GetIsothermal_Temperature(Marker_Tag);
  CFluidModel *fluid_model_local  = solver_container[FLOW_SOL]->GetFluidModel();    
  su2double enth_wall, prog_wall;
  unsigned long n_not_iterated    = 0;

  bool use_weak_bc                = config->GetUseWeakScalarBC();
  su2double *normal;
  su2double *coord_i, *coord_j;
  unsigned long point_normal;
  su2double area;
  su2double dist_ij;
  su2double dEnth_dn;
  su2double dT_dn;
  su2double mass_diffusivity;

  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      if (use_weak_bc){

        point_normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        normal = geometry->vertex[val_marker][iVertex]->GetNormal();

        area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) area += normal[iDim]*normal[iDim];
        area = sqrt(area);

        coord_i = geometry->nodes->GetCoord(iPoint);
        coord_j = geometry->nodes->GetCoord(point_normal);
        dist_ij = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          dist_ij += (coord_j[iDim]-coord_i[iDim])*(coord_j[iDim]-coord_i[iDim]);
        dist_ij = sqrt(dist_ij);

        prog_wall       = solver_container[SCALAR_SOL]->GetNodes()->GetSolution(iPoint)[I_PROG_VAR];
        n_not_iterated += fluid_model_local->GetEnthFromTemp(&enth_wall, prog_wall, temp_wall);

        dT_dn    = -(fluid_model_local->GetTemperature() - temp_wall)/dist_ij;

        dEnth_dn = -(solver_container[SCALAR_SOL]->GetNodes()->GetSolution(point_normal)[I_ENTHALPY]
                     - enth_wall)/dist_ij;

        mass_diffusivity = fluid_model_local->GetMassDiffusivity();

        Res_Visc[I_ENTHALPY] = mass_diffusivity*dEnth_dn*area;

        LinSysRes.SubtractBlock(iPoint, Res_Visc);

        if(implicit) {

          Jacobian_i[I_ENTHALPY][I_ENTHALPY] = -mass_diffusivity * area;
          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

        }

      } else {

        /*--- Set enthalpy on the wall ---*/

        prog_wall       = solver_container[SCALAR_SOL]->GetNodes()->GetSolution(iPoint)[I_PROG_VAR];
        n_not_iterated += fluid_model_local->GetEnthFromTemp(&enth_wall, prog_wall, temp_wall);

        /*--- Impose the value of the enthalpy as a strong boundary
        condition (Dirichlet) and remove any 
        contribution to the residual at this node. ---*/

        nodes->SetSolution(iPoint, I_ENTHALPY, enth_wall);
        nodes->SetSolution_Old(iPoint, I_ENTHALPY, enth_wall);

        LinSysRes(iPoint, I_ENTHALPY) = 0.0;

        nodes->SetVal_ResTruncError_Zero(iPoint, I_ENTHALPY);

        if (implicit) {

          total_index = iPoint*nVar+I_ENTHALPY;

          Jacobian.DeleteValsRowi(total_index);
        }
      }
    }

  }

  if (rank == MASTER_NODE && n_not_iterated > 0){
    cout << " !!! Isothermal wall bc ("  << Marker_Tag << "): Number of points in which enthalpy could not be iterated: " << n_not_iterated << " !!!" << endl;
  }

}
