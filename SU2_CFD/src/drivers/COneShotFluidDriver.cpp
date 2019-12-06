/*!
 * \file COneShotFluid.cpp
 * \brief The main subroutines for driving one-shot problems.
 * \author B. Munguía
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

#include "../../include/drivers/COneShotFluidDriver.hpp"

COneShotFluidDriver::COneShotFluidDriver(char* confFile,
                                         unsigned short val_nZone,
                                         SU2_Comm MPICommunicator) : CDiscAdjSinglezoneDriver(confFile, val_nZone, MPICommunicator) {
  unsigned short iDV, jDV, iDV_Value;

  /*---------- One-shot works on all design variables - find total number of design variables ---------*/
  nDV_Total = 0;
  for (iDV = 0; iDV  < config->GetnDV(); iDV++){
    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++){
      nDV_Total++;
    }
  }

  nConstr = config->GetnConstr();

  Gradient = new su2double[nDV_Total];
  Gradient_Old = new su2double[nDV_Total];

  ShiftLagGrad = new su2double[nDV_Total];
  ShiftLagGrad_Old = new su2double[nDV_Total];

  AugLagGrad = new su2double[nDV_Total];
  AugLagGradAlpha = new su2double[nDV_Total];
  AugLagGradBeta = new su2double[nDV_Total];
  AugLagGradGamma = new su2double*[nDV_Total];
  AugLagGrad_Old = new su2double[nDV_Total];

  DesignVarUpdate = new su2double[nDV_Total];
  DesignVar = new su2double[nDV_Total];
  SearchDirection = new su2double[nDV_Total];
  ActiveSetDV = new bool[nDV_Total];

  BFGS_Inv = new su2double*[nDV_Total];

  if(nConstr > 0){
    ConstrFunc = new su2double[nConstr];
    ConstrFunc_Store = new su2double[nConstr];
    Lambda = new su2double[nConstr];
    Lambda_Old = new su2double[nConstr];
    Lambda_Store = new su2double[nConstr];
    Lambda_Store_Old = new su2double[nConstr];
    AugLagLamGrad = new su2double[nConstr];
    BCheck_Inv = new su2double*[nConstr];
  }

  BFGS_Init = config->GetBFGSInitValue();

  for (iDV = 0; iDV  < nDV_Total; iDV++){
    Gradient[iDV] = 0.0;
    Gradient_Old[iDV] = 0.0;
    ShiftLagGrad[iDV] = 0.0;
    ShiftLagGrad_Old[iDV] = 0.0;
    AugLagGrad[iDV] = 0.0;
    AugLagGradAlpha[iDV] = 0.0;
    AugLagGradBeta[iDV] = 0.0;
    AugLagGradGamma[iDV] = new su2double[nConstr];
    for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++) AugLagGradGamma[iDV][iConstr] = 0.0;
    AugLagGrad_Old[iDV] = 0.0;
    DesignVarUpdate[iDV] = 0.0;
    DesignVar[iDV] = 0.0;
    SearchDirection[iDV] = 0.0;
    ActiveSetDV[iDV]=false;
    BFGS_Inv[iDV] = new su2double[nDV_Total];
    for (jDV = 0; jDV < nDV_Total; jDV++){
      BFGS_Inv[iDV][jDV] = 0.0;
      if (iDV==jDV) BFGS_Inv[iDV][jDV] = BFGS_Init;
    }
  }

  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    ConstrFunc[iConstr] = 0.0;
    ConstrFunc_Store[iConstr] = 0.0;
    Lambda[iConstr] = config->GetMultiplierStart(iConstr);
    Lambda_Old[iConstr] = config->GetMultiplierStart(iConstr);
    Lambda_Store[iConstr] = config->GetMultiplierStart(iConstr);
    Lambda_Store_Old[iConstr] = config->GetMultiplierStart(iConstr);
    AugLagLamGrad[iConstr] = 0.0;
    BCheck_Inv[iConstr] = new su2double[nConstr];
    for (unsigned short jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck_Inv[iConstr][jConstr] = 0.0;
    }
    BCheck_Inv[iConstr][iConstr] = config->GetOneShotGamma(iConstr);
  }
  BCheck_Norm = 1./config->GetOneShotGamma(0);

  /*----- calculate values for bound projection algorithm -------*/
  lb=-config->GetBound()*config->GetDesignScale();
  ub=config->GetBound()*config->GetDesignScale();
  epsilon=(ub-lb)/2.0;

  /*---- calculate line search parameter ----*/
  CWolfeOne= 1E-4*config->GetDesignScale();

  grid_movement[ZONE_0][INST_0] = new CVolumetricMovement(geometry, config);
  surface_movement[ZONE_0]      = new CSurfaceMovement();

}

COneShotFluidDriver::~COneShotFluidDriver(void){

  /*----- free allocated memory -------*/
  unsigned short iDV, iConstr;
  for (iDV = 0; iDV  < nDV_Total; iDV++){
    delete [] BFGS_Inv[iDV];
  }
  delete [] BFGS_Inv;
  delete [] Gradient;
  delete [] Gradient_Old;
  delete [] ShiftLagGrad;
  delete [] ShiftLagGrad_Old;

  delete [] AugLagGrad;
  delete [] AugLagGradAlpha;
  delete [] AugLagGradBeta;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    delete [] AugLagGradGamma[iDV];
  }
  delete [] AugLagGradGamma;
  delete [] AugLagGrad_Old;

  delete [] DesignVarUpdate;
  delete [] DesignVar;
  delete [] SearchDirection;
  delete [] ActiveSetDV;

  if(nConstr > 0){
    delete [] ConstrFunc;
    delete [] Lambda;
    delete [] Lambda_Old;
    delete [] Lambda_Store;
    delete [] Lambda_Store_Old;
    delete [] ConstrFunc_Store;
    delete [] AugLagLamGrad;
  }

}

void COneShotFluidDriver::Preprocess(unsigned long TimeIter) {

  config->SetTimeIter(TimeIter);

}


void COneShotFluidDriver::Run(){

  unsigned long OneShotIter, nOneShotIter = config->GetnInner_Iter();;

  for(OneShotIter = 0; OneShotIter < nOneShotIter; OneShotIter++) {

    config->SetInnerIter(OneShotIter);

    /*--- Preprocess the one-shot iteration ---*/

    iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
   
    /*--- Run an iteration of the one-shot solver ---*/
    RunOneShot();

    /*--- Screen output ---*/
    StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    if(StopCalc) {
      for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
        if(config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR && ConstrFunc[iConstr] > 0.) {
          StopCalc = false;
          break;
        }
      }
    }

    if(StopCalc) break;

  }

}

void COneShotFluidDriver::RunOneShot(){

  su2double stepsize = 1.0, stepsize_tmp, tol = config->GetOneShotSearchTol();
  unsigned short ArmijoIter = 0, nArmijoIter = config->GetOneShotSearchIter();
  unsigned long InnerIter = config->GetInnerIter();
  bool bool_tol = false;
  unsigned short ALPHA_TERM = 0, BETA_TERM = 1, GAMMA_TERM = 2, TOTAL_AUGMENTED = 3, TOTAL_AUGMENTED_OLD = 4;

  /*--- Store the old solution and the old design for line search ---*/
  solver[ADJFLOW_SOL]->SetStoreSolution();
  solver[ADJFLOW_SOL]->SetMeshPointsOld(config, geometry);

  /*--- This is the line search loop that is only called once, if no update is performed ---*/
  do {

    if(InnerIter > config->GetOneShotStart() && InnerIter < config->GetOneShotStop()){      

      if(ArmijoIter > 0){
        /*--- Parabolic backtracking ---*/
        stepsize_tmp = UpdateStepSizeQuadratic();
        stepsize     = UpdateStepSizeBound(stepsize_tmp, stepsize/10., stepsize/2.);
        if(stepsize < tol) {
          stepsize = tol;
          bool_tol = true;
        }

        /*---Load the old design for line search---*/
        solver[ADJFLOW_SOL]->LoadMeshPointsOld(config, geometry);

      }
      // else{
      //   /*--- Store gradient of augmented Lagrangian wrt multiplier ---*/
      //   StoreLambdaGrad();
      // }

      // /*--- Compute and store GradL dot p ---*/
      // StoreGradDotDir(true);

      // /*--- Update constraint multiplier ---*/
      // LoadOldLambda();
      // UpdateLambda(stepsize);

      /*--- Load the old solution for line search (either y_k or y_k-1) ---*/
      solver[ADJFLOW_SOL]->LoadSolution();

      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      if ((ArmijoIter != nArmijoIter-1 && !bool_tol) || (!config->GetZeroStep())) {
        ComputeDesignVarUpdate(stepsize);
        config->SetKind_SU2(SU2_DEF); // set SU2_DEF as the solver
        SurfaceDeformation(geometry, config, surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0]);
        config->SetKind_SU2(SU2_CFD); // set SU2_CFD as the solver

        /*--- Evaluate the objective at the old solution, new design ---*/
      
        solver[FLOW_SOL]->Pressure_Forces(geometry, config);
        solver[FLOW_SOL]->Momentum_Forces(geometry, config);
        solver[FLOW_SOL]->Friction_Forces(geometry, config);
                  
        if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
            solver[FLOW_SOL]->Buffet_Monitoring(geometry, config);
        }
        SetObjFunction();
        ObjFunc_Store = ObjFunc;
      }
      else {
        stepsize = 0.0;
        grid_movement[ZONE_0][INST_0]->UpdateDualGrid(geometry, config);
        ComputeDesignVarUpdate(0.0);
      }

    }

    /*--- Do a primal and adjoint update ---*/
    PrimalDualStep();
    solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);

    if(InnerIter > config->GetOneShotStart() && InnerIter < config->GetOneShotStop()){
      // StoreLambdaGrad();
      /*--- Update constraint multiplier ---*/
      LoadOldLambda();
      UpdateLambda(stepsize);
      // UpdateLambda(1.0);

      /*--- Compute and store GradL dot p ---*/
      StoreGradDotDir(true);
    }

    /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
    if((Armijo_Iter != nArmijoIter-1) || (!config->GetZeroStep())) CalculateLagrangian();

    ArmijoIter++;

  } while((InnerIter > config->GetOneShotStart()) && 
          (InnerIter < config->GetOneShotStop())  &&
          (!CheckFirstWolfe(true)) && (ArmijoIter < nArmijoIter) && (!bool_tol));

  /*--- Store number of search iterations ---*/
  solver[ADJFLOW_SOL]->SetArmijoIter(ArmijoIter);

  /*--- Store FFD info in file ---*/
  if (((config->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) ||
       (config->GetDesign_Variable(0) == FFD_CONTROL_POINT))   &&
       InnerIter > config->GetOneShotStart()                   && 
       InnerIter < config->GetOneShotStop()                    &&
       (!config->GetZeroStep() || CheckFirstWolfe(true))) {
    surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container[ZONE_0][INST_0], config_container, false);
    config->SetMesh_FileName(config->GetMesh_Out_FileName());
  }

  /*--- Compute alpha, beta, gamma at first one-shot iteration, or recompute if line search failed ---*/
  if(InnerIter > 0) solver[ADJFLOW_SOL]->CalculateRhoTheta(config);
  if(InnerIter == config->GetOneShotStart()) {
    solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    if((nConstr > 0) && (!config->GetConstPrecond())) ComputePreconditioner();
    solver[ADJFLOW_SOL]->CalculateGamma(config, BCheck_Norm, ConstrFunc);

    /*--- Recalculate Lagrangian with new Alpha, Beta, and Gamma ---*/
    CalculateLagrangian();
    SetAugLagGrad(TOTAL_AUGMENTED_OLD);
  }
  else if(InnerIter > config->GetOneShotStart() && 
          InnerIter < config->GetOneShotStop()  && 
          ((!CheckFirstWolfe(true)) || (ArmijoIter > nArmijoIter-1) || (bool_tol))){
    /*--- Perform new line search on just multiplier ---*/
    if(nConstr > 0 && config->GetZeroStep()) {
      su2double stepsize_mu = 1.0;
      ArmijoIter = 0;
      bool_tol = false;
      do {
        if(ArmijoIter > 0){
          /*--- Parabolic backtracking ---*/
          stepsize_tmp = UpdateStepSizeQuadratic();
          stepsize_mu  = UpdateStepSizeBound(stepsize_tmp, stepsize_mu/10., stepsize_mu/2.);
          if(stepsize_mu < tol) {
            stepsize_mu  = 0.;
            bool_tol     = true;
          }
        }
        /*--- Compute and store GradL dot p ---*/
        StoreGradDotDir(false);

        /*--- Update constraint multiplier ---*/
        LoadOldLambda();
        UpdateLambda(stepsize_mu);

        /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
        CalculateLagrangian();

        ArmijoIter++;

      } while((!CheckFirstWolfe(false)) && (ArmijoIter < nArmijoIter) && (!bool_tol));
    }
    // solver[ADJFLOW_SOL]->LoadStepSolution(1.0);
    // solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);
    // LoadOldLambda();
    // UpdateLambda(1.0);
    solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    solver[ADJFLOW_SOL]->CalculateGamma(config, BCheck_Norm, ConstrFunc);

    /*--- Recalculate Lagrangian with new Alpha, Beta, and Gamma ---*/
    CalculateLagrangian();
    SetAugLagGrad(TOTAL_AUGMENTED_OLD);
  }

  if(InnerIter >= config->GetOneShotStart() && 
     InnerIter < config->GetOneShotStop()   && 
     InnerIter > 0) {
    /*--- Calculate Lagrangian with new Alpha, Beta, and Gamma ---*/
    // solver[ADJFLOW_SOL]->CalculateRhoTheta(config);
    // // solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    // // solver[ADJFLOW_SOL]->CalculateGamma(config, BCheck_Norm);
    /*--- Store the multiplier and constraint function, then recalculate Lagrangian for next iteration ---*/
    StoreOldLambda();
    StoreConstrFunction();
    // CalculateLagrangian();
    // SetAugLagGrad(TOTAL_AUGMENTED_OLD);
  }

  /*--- Store Deltay and DeltaBary ---*/
  solver[ADJFLOW_SOL]->SetStoreSolutionDelta();

  if(InnerIter >= config->GetOneShotStart() && InnerIter < config->GetOneShotStop()){

    /*--- Update design variable ---*/
    UpdateDesignVar();

    /*--- N_u ---*/
    solver[ADJFLOW_SOL]->SetSensitivityShiftedLagrangian(geometry);
    solver[ADJFLOW_SOL]->SetSaveSolution();
    solver[ADJFLOW_SOL]->LoadSolution();
    solver[ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry);
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, 1.0);

    if((nConstr > 0) && (!config->GetConstPrecond())) ComputePreconditioner();

    /*--- Gamma*h^T*h_u ---*/
    if(nConstr > 0) {
      ComputeGammaTerm();
      solver[ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry);
      solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, 1.0);
      solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry); //Lagrangian
      ProjectMeshSensitivities();
      SetAugLagGrad(GAMMA_TERM);
    }
    solver[ADJFLOW_SOL]->LoadSolution();

    /*--- Alpha*Deltay^T*G_u ---*/
    ComputeAlphaTerm();
    solver[ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry);
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, 1.0);
    solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry); //Lagrangian
    ProjectMeshSensitivities();
    SetAugLagGrad(ALPHA_TERM);
    solver[ADJFLOW_SOL]->LoadSolution();

    /*--- Beta*DeltaBary^T*N_yu ---*/
    solver[ADJFLOW_SOL]->UpdateStateVariable(config);
    ComputeBetaTerm();
    solver[ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry, config);
    solver[ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry);
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, 1.0);
    solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry); //Lagrangian
    ProjectMeshSensitivities();
    SetAugLagGrad(BETA_TERM);
    solver[ADJFLOW_SOL]->LoadSaveSolution();

    /*--- Projection of the gradient N_u---*/
    solver[ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry);
    ProjectMeshSensitivities();
    SetShiftLagGrad();

    /*--- Projection of the gradient L_u---*/
    SetAugLagGrad(TOTAL_AUGMENTED);

    /*--- Use N_u to compute the active set (bound constraints) ---*/
    ComputeActiveSet(stepsize);

    /*--- Do a BFGS update to approximate the inverse preconditioner ---*/
    if(InnerIter > config->GetOneShotStart()) BFGSUpdate(config);

    /*--- Compute the search direction for the line search procedure ---*/
    ComputeSearchDirection();

    StoreLagrangianInformation();
  }
}

void COneShotFluidDriver::PrimalDualStep(){

  /*--- Note: Unsteady cases not applicable to the one-shot method yet! ---*/

  SetRecording(COMBINED);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);


  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  SetAdj_ConstrFunction(Lambda);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->Iterate(output_container[ZONE_0], integration_container, geometry_container,
                     solver_container, numerics_container, config_container,
                     surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Extract the computed sensitivity values. ---*/
  solver[ADJFLOW_SOL]->SetSensitivity(geometry,solver,config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

}

void COneShotFluidDriver::SetRecording(unsigned short kind_recording){
  unsigned long InnerIter = config->GetInnerIter();

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the flow iteration (conservative variables and node coordinates) --- */

  if (kind_recording != NONE){

    AD::StartRecording();

    if (rank == MASTER_NODE && kind_recording == MainVariables && (InnerIter == 0)) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store the primal computational graph." << endl;
      cout << "Combined recording of flow and design variables." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }

    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  }

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0, INST_0, kind_recording);

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun(kind_recording);

  RecordingState = kind_recording;

  iteration->RegisterOutput(solver_container, geometry_container, config_container, output_container[ZONE_0], ZONE_0, INST_0);

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  SetConstrFunction();

  AD::StopRecording();

}

void COneShotFluidDriver::SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double* Gradient){

  su2double DV_Value, *VarCoord, Sensitivity, my_Gradient, localGradient;
  unsigned short iDV_Value = 0, iMarker, nMarker, iDim, nDim, iDV, nDV, nDV_Value;
  unsigned long iVertex, nVertex, iPoint;
  unsigned long nDV_Count = 0;

  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  VarCoord = NULL;

  AD::Reset();

  /*--- Discrete adjoint gradient computation ---*/

  /*--- Start recording of operations ---*/

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  for (iDV = 0; iDV < nDV; iDV++){

    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization with su2double resets the index ---*/

      DV_Value = 0.0;

      AD::RegisterInput(DV_Value);

      config->SetDV_Value(iDV, iDV_Value, DV_Value);

    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Stop the recording --- */

  AD::StopRecording();

  /*--- Create a structure to identify points that have been already visited.
   * We need that to make sure to set the sensitivity of surface points only once
   *  (Markers share points, so we would visit them more than once in the loop over the markers below) ---*/

  bool* visited = new bool[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
    visited[iPoint] = false;
  }

  /*--- Initialize the derivatives of the output of the surface deformation routine
   * with the discrete adjoints from the CFD solution ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint      = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!visited[iPoint]){
          VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();

          for (iDim = 0; iDim < nDim; iDim++){
            Sensitivity = geometry->GetSensitivity(iPoint, iDim);
            SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
          }
          visited[iPoint] = true;
        }
      }
    }
  }

  delete [] visited;

  /*--- Compute derivatives and extract gradient ---*/

  AD::ComputeAdjoint();

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      DV_Value = config->GetDV_Value(iDV, iDV_Value);
      my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      localGradient = my_Gradient;
#endif

      Gradient[nDV_Count] = localGradient;
      nDV_Count++;
    }
  }
  
  AD::ClearAdjoints();

  AD::Reset();
}

void COneShotFluidDriver::SurfaceDeformation(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement){

  unsigned short iMarker, iDV, iDV_Value, nDV_Value;
  bool allmoving=true;
  unsigned long nDV_Count = 0;

  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      config->SetDV_Value(iDV,iDV_Value, DesignVarUpdate[nDV_Count]/config->GetDesignScale());
      nDV_Count++;
    }
  }

  /*--- Surface grid deformation using design variables ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- For scale, translation and rotation if all boundaries are moving they are set via volume method
   * Otherwise, the surface deformation has been set already in SetSurface_Deformation.  --- */
  /*--- Loop over markers, set flag to false if any are not moving ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    if (config->GetMarker_All_DV(iMarker) == NO)
      allmoving = false;
  }

  /*--- Volumetric grid deformation/transformations ---*/

  if (config->GetDesign_Variable(0) == SCALE && allmoving) {

    grid_movement->SetVolume_Scaling(geometry, config, false);

  } else if (config->GetDesign_Variable(0) == TRANSLATION && allmoving) {

    grid_movement->SetVolume_Translation(geometry, config, false);

  } else if (config->GetDesign_Variable(0) == ROTATION && allmoving) {

    grid_movement->SetVolume_Rotation(geometry, config, false);

  } else if (config->GetDesign_Variable(0) != FFD_SETTING) {

    grid_movement->SetVolume_Deformation(geometry, config, true, false);

  }

}

void COneShotFluidDriver::BFGSUpdate(CConfig *config){
  unsigned long iDV, jDV, kDV, lDV;

  su2double *yk, *sk;
  su2double vk=0;
  su2double normyk=0;
  su2double normsk=0;

  yk=new su2double[nDV_Total];
  sk=new su2double[nDV_Total];
  for (iDV = 0; iDV < nDV_Total; iDV++){
    yk[iDV]=ProjectionSet(iDV, AugLagGrad[iDV]-AugLagGrad_Old[iDV], false);
    sk[iDV]=ProjectionSet(iDV, DesignVarUpdate[iDV], false);
    vk+=yk[iDV]*sk[iDV];
    normyk+=yk[iDV]*yk[iDV];
    normsk+=sk[iDV]*sk[iDV];
  }

  if (vk>0){
    su2double** MatA;
    MatA=new su2double*[nDV_Total];
    for (iDV=0;iDV<nDV_Total;iDV++){
      MatA[iDV]=new su2double[nDV_Total];
      for (jDV=0;jDV<nDV_Total;jDV++){
          MatA[iDV][jDV]=0.0;
      }
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      for (jDV=0;jDV<nDV_Total;jDV++){
        MatA[iDV][jDV]=ProjectionPAP(iDV,jDV,BFGS_Inv[iDV][jDV],false)+(1.0/vk)*sk[iDV]*sk[jDV];
        for (kDV=0; kDV<nDV_Total; kDV++){
          MatA[iDV][jDV]+=-(1.0/vk)*sk[iDV]*ProjectionPAP(kDV,jDV,BFGS_Inv[kDV][jDV],false)*yk[kDV]-(1.0/vk)*sk[jDV]*ProjectionPAP(iDV,kDV,BFGS_Inv[iDV][kDV],false)*yk[kDV];
          for (lDV=0; lDV<nDV_Total; lDV++){
            MatA[iDV][jDV]+=(1.0/vk)*(1.0/vk)*sk[iDV]*sk[jDV]*yk[lDV]*ProjectionPAP(lDV,kDV,BFGS_Inv[lDV][kDV],false)*yk[kDV];
          }
        }
      }
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      for (jDV=0;jDV<nDV_Total;jDV++){
        BFGS_Inv[iDV][jDV]=MatA[iDV][jDV];
      }
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      delete [] MatA[iDV];
    }
    delete [] MatA;
    if(config->GetBFGSInit()){
      for (iDV=0;iDV<nDV_Total;iDV++){
        BFGS_Init = vk/normyk;
      }
    }

  }else{
    /*--- Calculate new alpha, beta, gamma, and reset BFGS update if needed ---*/
    // unsigned short TOTAL_AUGMENTED = 3;
    // solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    // solver[ADJFLOW_SOL]->CalculateGamma(config, BCheck_Norm);
    // CalculateLagrangian();
    // SetAugLagGrad(TOTAL_AUGMENTED);
    if(config->GetBoolBFGSReset()){
      for (iDV = 0; iDV < nDV_Total; iDV++){
        for (jDV = 0; jDV < nDV_Total; jDV++){
          BFGS_Inv[iDV][jDV]=0.0;
          if(iDV==jDV){ BFGS_Inv[iDV][jDV]=ProjectionSet(iDV,BFGS_Init,false); }
        }
      }
    }
  }
  delete [] yk;
  delete [] sk;
}

bool COneShotFluidDriver::CheckFirstWolfe(bool design_update){
  unsigned short iDV;
  su2double admissible_step = 0.0;

  if(design_update) {
    for (iDV = 0; iDV < nDV_Total; iDV++){
      // /*--- ShiftLagGrad is the gradient at the old iterate. ---*/
      // admissible_step += DesignVarUpdate[iDV]*ShiftLagGrad[iDV];
      /*--- AugLagGrad is the gradient at the old iterate. ---*/
      admissible_step += DesignVarUpdate[iDV]*AugLagGrad[iDV];
    }
  }
  if (nConstr > 0) {
    unsigned short iConstr;
    for (iConstr = 0; iConstr < nConstr; iConstr++) {
      // admissible_step += (Lambda[iConstr]-Lambda_Old[iConstr])*AugLagLamGrad[iConstr];
      const su2double gamma = config->GetOneShotGamma(iConstr);
      const su2double dh = ConstrFunc[iConstr]-ConstrFunc_Store[iConstr];
      const su2double hdh = ConstrFunc_Store[iConstr]*dh;
      // const bool active = (ConstrFunc_Store[iConstr] - Lambda_Old[iConstr]/gamma > 0.);
      const bool active = (ConstrFunc_Store[iConstr] > 0.);
      // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (hdh <= 0.)) || 
      //    ((active) && (dh <= 0.))) {
      if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
        // admissible_step += (Lambda[iConstr]-Lambda_Old[iConstr])*ConstrFunc_Store[iConstr];
        admissible_step += (Lambda[iConstr]-Lambda_Old[iConstr])*AugLagLamGrad[iConstr];
      }
    }
  }
  admissible_step *= CWolfeOne;

  return (Lagrangian <= Lagrangian_Old + admissible_step);
}

void COneShotFluidDriver::StoreGradDotDir(bool design_update){
  unsigned short iDV;
  GradDotDir = 0.0;

  if(design_update) {
    for (iDV = 0; iDV < nDV_Total; iDV++){
      // /*--- ShiftLagGrad is the gradient at the old iterate. ---*/
      // GradDotDir += DesignVarUpdate[iDV]*ShiftLagGrad[iDV];
      /*--- AugLagGrad is the gradient at the old iterate. ---*/
      GradDotDir += DesignVarUpdate[iDV]*AugLagGrad[iDV];
    }
  }
  if (nConstr > 0) {
    unsigned short iConstr;
    for (iConstr = 0; iConstr < nConstr; iConstr++) {
      // GradDotDir += (Lambda[iConstr]-Lambda_Old[iConstr])*AugLagLamGrad[iConstr];
      const su2double gamma = config->GetOneShotGamma(iConstr);
      const su2double dh = ConstrFunc[iConstr]-ConstrFunc_Store[iConstr];
      const su2double hdh = ConstrFunc_Store[iConstr]*dh;
      // const bool active = (ConstrFunc_Store[iConstr] - Lambda_Old[iConstr]/gamma > 0.);
      const bool active = (ConstrFunc_Store[iConstr] > 0.);
      // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (hdh <= 0.)) || 
      //    ((active) && (dh <= 0.))) {
      if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
        // GradDotDir += (Lambda[iConstr]-Lambda_Old[iConstr])*ConstrFunc_Store[iConstr];
        GradDotDir += (Lambda[iConstr]-Lambda_Old[iConstr])*AugLagLamGrad[iConstr];
      }
    }
  }
}

su2double COneShotFluidDriver::UpdateStepSizeQuadratic(){
  return -GradDotDir/(2.*(Lagrangian - Lagrangian_Old - GradDotDir));
}

su2double COneShotFluidDriver::UpdateStepSizeBound(su2double stepsize, su2double a, su2double b){
  if(stepsize < a) {return a;}
  else if(stepsize > b) {return b;}
  else {return stepsize;}
}


void COneShotFluidDriver::ComputeDesignVarUpdate(su2double stepsize){
  unsigned short iDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    DesignVarUpdate[iDV]=BoundProjection(DesignVar[iDV]+stepsize*SearchDirection[iDV]*config->GetDesignScale())-DesignVar[iDV];
  }
}

void COneShotFluidDriver::ComputeSearchDirection(){
  unsigned short iDV, jDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    SearchDirection[iDV]=0.0;
    for (jDV=0;jDV<nDV_Total;jDV++){
      SearchDirection[iDV]+= BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,-ShiftLagGrad[jDV],false);
    }
    SearchDirection[iDV]=-ProjectionSet(iDV, ShiftLagGrad[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false);
  }
}

void COneShotFluidDriver::StoreLagrangianInformation(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    AugLagGrad_Old[iDV] = AugLagGrad[iDV];
  }
  Lagrangian_Old = Lagrangian;
}

void COneShotFluidDriver::UpdateDesignVar(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    DesignVar[iDV] += DesignVarUpdate[iDV];
  }
}

void COneShotFluidDriver::CalculateLagrangian(){
    
  Lagrangian = 0.0;
  Lagrangian += ObjFunc_Store; //TODO use for BFGS either only objective function or normal Lagrangian

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma(iConstr);
    const su2double helper = ConstrFunc_Store[iConstr] + Lambda[iConstr]/gamma;
    // const bool active = (ConstrFunc_Store[iConstr] - Lambda_Old[iConstr]/gamma > 0.);
    const bool active = (ConstrFunc_Store[iConstr] > 0.);
    /*--- Lagrangian += gamma/2 ||h + mu/gamma - P_I(h+mu/gamma)||^2 ---*/
    if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
      Lagrangian += gamma/2.*helper*helper - 1./(2.*gamma)*Lambda[iConstr]*Lambda[iConstr];
      // Lagrangian += gamma/2.*helper*helper;
    }
  }

  Lagrangian += solver[ADJFLOW_SOL]->CalculateLagrangian(config);

}

su2double COneShotFluidDriver::ProjectionSet(unsigned short iDV, su2double value, bool active){
  if (active) {
      if(!ActiveSetDV[iDV]) value = 0.0;
  } else {
      if(ActiveSetDV[iDV]) value = 0.0;
  }
  return value;
}

su2double COneShotFluidDriver::ProjectionPAP(unsigned short iDV, unsigned short jDV, su2double value, bool active){
  //returns for a Matrix entry a_iDV,jDV of a Matrix A the resulting entry of P*A*P (active or inactive set)
  if (active) {
      if(!ActiveSetDV[iDV] || !ActiveSetDV[jDV]) value = 0.0;
  } else {
      if(ActiveSetDV[iDV] || ActiveSetDV[jDV]) value = 0.0;
  }
  return value;
}

su2double COneShotFluidDriver::BoundProjection(su2double value){
  if(value <= lb) value = lb;
  if(value >= ub) value = ub;
  return value;
}

void COneShotFluidDriver::ComputeActiveSet(su2double stepsize){
  //Compute ||x-P(x-gradx)||
  unsigned short iDV;
  su2double norm = 0.0;

  for (iDV = 0; iDV < nDV_Total; iDV++) {
    norm += (DesignVar[iDV]-BoundProjection(DesignVar[iDV]-stepsize*ShiftLagGrad[iDV]))
          * (DesignVar[iDV]-BoundProjection(DesignVar[iDV]-stepsize*ShiftLagGrad[iDV]));
  }
  norm = sqrt(norm);
  epsilon = min(norm, (ub-lb)/2.0);
  unsigned short nActive = nDV_Total;

  for (iDV = 0; iDV < nDV_Total; iDV++) {
    ActiveSetDV[iDV] = false;
    if(ub-DesignVar[iDV] <= epsilon)      ActiveSetDV[iDV] = true;
    else if(DesignVar[iDV]-lb <= epsilon) ActiveSetDV[iDV] = true;
    else                                       nActive--;
  }
  solver[ADJFLOW_SOL]->SetnActiveDV(nActive);
}

void COneShotFluidDriver::SetShiftLagGrad(){
  unsigned short iDV;
  su2double norm = 0.;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    ShiftLagGrad[iDV] = Gradient[iDV];
    norm += Gradient[iDV]*Gradient[iDV];
  }
  solver[ADJFLOW_SOL]->SetShiftedLagGradNorm(sqrt(norm));
}

void COneShotFluidDriver::SetAugLagGrad(unsigned short kind){
  unsigned short iDV, iConstr;
  unsigned short ALPHA_TERM = 0, BETA_TERM = 1, GAMMA_TERM = 2, TOTAL_AUGMENTED = 3, TOTAL_AUGMENTED_OLD = 4;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    if(kind == ALPHA_TERM) {
      AugLagGradAlpha[iDV] = Gradient[iDV];
    }
    else if(kind == BETA_TERM) {
      AugLagGradBeta[iDV] = Gradient[iDV];
    }
    else if(kind == GAMMA_TERM) {
      for(iConstr = 0; iConstr < nConstr; iConstr++) {
        AugLagGradGamma[iDV][iConstr] = Gradient[iDV];
      }
    }
    else if(kind == TOTAL_AUGMENTED) {
      AugLagGrad[iDV] = ShiftLagGrad[iDV]
                      + AugLagGradAlpha[iDV]*config->GetOneShotAlpha()
                      + AugLagGradBeta[iDV]*config->GetOneShotBeta();
      for(iConstr = 0; iConstr < nConstr; iConstr++) {
        AugLagGrad[iDV] += AugLagGradGamma[iDV][iConstr]*config->GetOneShotGamma(iConstr);
      }
    }   
    else if(kind == TOTAL_AUGMENTED_OLD) {
      AugLagGrad_Old[iDV] = ShiftLagGrad[iDV]
                          + AugLagGradAlpha[iDV]*config->GetOneShotAlpha()
                          + AugLagGradBeta[iDV]*config->GetOneShotBeta();
      for(iConstr = 0; iConstr < nConstr; iConstr++) {
        AugLagGrad_Old[iDV] += AugLagGradGamma[iDV][iConstr]*config->GetOneShotGamma(iConstr);
      }
    }   
  }
}

void COneShotFluidDriver::ComputeGammaTerm(){

  /*--- Note: Not applicable for unsteady code ---*/

  /*--- Initialize the adjoint of the output variables of the iteration with difference of the solution and the solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint_Zero(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 0.0. ---*/

  SetAdj_ObjFunction_Zero();
  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma(iConstr);
    // const bool active = (ConstrFunc_Store[iConstr] - Lambda_Old[iConstr]/gamma > 0.);
    const bool active = (ConstrFunc_Store[iConstr] > 0.);
    if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
      seeding[iConstr] = ConstrFunc[iConstr];
    }
    else {
      seeding[iConstr] = 0.;
    }
  }
  SetAdj_ConstrFunction(seeding);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                 solver_container, numerics_container, config_container,
                                 surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Extract the computed sensitivity values. ---*/
  solver[ADJFLOW_SOL]->SetSensitivity(geometry,solver,config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  delete [] seeding;
}

void COneShotFluidDriver::ComputeAlphaTerm(){

  /*--- Note: Not applicable for unsteady code ---*/

  /*--- Initialize the adjoint of the output variables of the iteration with difference of the solution and the solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint_Update(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 0.0. ---*/

  SetAdj_ObjFunction_Zero();
  SetAdj_ConstrFunction_Zero();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                 solver_container, numerics_container, config_container,
                                 surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Extract the computed sensitivity values. ---*/
  solver[ADJFLOW_SOL]->SetSensitivity(geometry,solver,config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  AD::Reset();
}

void COneShotFluidDriver::ComputeBetaTerm(){

  /*--- Note: Not applicable for unsteady code ---*/

  /*--- For the one shot iteration we have to record for every steady state iteration. ---*/

  /*--- Store the computational graph of one direct iteration with the conservative variables and the mesh coordinates as input. ---*/

  SetRecording(COMBINED);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  SetAdj_ConstrFunction(Lambda);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                 solver_container, numerics_container, config_container,
                                 surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Extract the computed sensitivity values. ---*/
  solver[ADJFLOW_SOL]->SetSensitivity(geometry, solver, config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  AD::Reset();

}

void COneShotFluidDriver::ComputePreconditioner(){

  unsigned short INST_0 = 0;

  unsigned short iConstr, jConstr;

  su2double* seeding = new su2double[nConstr];
  for (iConstr = 0; iConstr < nConstr; iConstr++){
    seeding[iConstr] = 0.0;
  }
  su2double **BCheck = new su2double*[nConstr];
  for (iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr] = new su2double[nConstr];
    for (jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck[iConstr][jConstr] = 0.0;
    }
  }

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    seeding[iConstr] = 1.0;


    iteration->InitializeAdjoint_Zero(solver_container, geometry_container, config_container, ZONE_0, INST_0);


    /*--- Initialize the adjoint of the objective function with 0.0. ---*/

    SetAdj_ObjFunction_Zero();
    SetAdj_ConstrFunction(seeding);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                   solver_container, numerics_container, config_container,
                                   surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    solver[ADJFLOW_SOL]->SetConstrDerivative(iConstr);


    AD::ClearAdjoints();

    solver[ADJFLOW_SOL]->LoadSolution();

    seeding[iConstr]=0.0;

  }

  su2double bcheck=0;
  for (iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr][iConstr] = 1./config->GetOneShotGamma(iConstr);
    for (jConstr = 0; jConstr < nConstr; jConstr++){
      BCheck[iConstr][jConstr] += config->GetOneShotBeta()*solver[ADJFLOW_SOL]->MultiplyConstrDerivative(iConstr,jConstr);
    }
  }
  if (nConstr == 1){
      BCheck_Norm = BCheck[0][0] - 1./config->GetOneShotGamma(0);
      BCheck_Inv[0][0] = 1./BCheck[0][0];
  } else {
    bcheck=1./(BCheck[0][0]*BCheck[1][1]*BCheck[2][2]+BCheck[1][0]*BCheck[2][1]*BCheck[0][2]+BCheck[2][0]*BCheck[0][1]*BCheck[1][2]-BCheck[0][0]*BCheck[2][1]*BCheck[1][2]-BCheck[2][0]*BCheck[1][1]*BCheck[0][2]-BCheck[1][0]*BCheck[0][1]*BCheck[2][2]);
    BCheck_Inv[0][0]=bcheck*(BCheck[1][1]*BCheck[2][2]-BCheck[1][2]*BCheck[2][1]);
    BCheck_Inv[0][1]=bcheck*(BCheck[0][2]*BCheck[2][1]-BCheck[0][1]*BCheck[2][2]);
    BCheck_Inv[0][2]=bcheck*(BCheck[0][1]*BCheck[1][2]-BCheck[0][2]*BCheck[1][1]);
    BCheck_Inv[1][0]=bcheck*(BCheck[1][2]*BCheck[2][0]-BCheck[1][0]*BCheck[2][2]);
    BCheck_Inv[1][1]=bcheck*(BCheck[0][0]*BCheck[2][2]-BCheck[0][2]*BCheck[2][0]);
    BCheck_Inv[1][2]=bcheck*(BCheck[0][2]*BCheck[1][0]-BCheck[0][0]*BCheck[1][2]);
    BCheck_Inv[2][0]=bcheck*(BCheck[1][0]*BCheck[2][1]-BCheck[1][1]*BCheck[2][0]);
    BCheck_Inv[2][1]=bcheck*(BCheck[0][1]*BCheck[2][0]-BCheck[0][0]*BCheck[2][1]);
    BCheck_Inv[2][2]=bcheck*(BCheck[0][0]*BCheck[1][1]-BCheck[0][1]*BCheck[1][0]);
  }

  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    delete [] BCheck[iConstr];
  }
  delete [] BCheck;
  delete [] seeding;
}

void COneShotFluidDriver::SetAdj_ObjFunction_Zero(){
  SU2_TYPE::SetDerivative(ObjFunc, 0.0);
}

void COneShotFluidDriver::ProjectMeshSensitivities(){
  config->SetKind_SU2(SU2_DOT); // set SU2_DOT as solver
  // get the dependency of the volumetric grid movement
  grid_movement[ZONE_0][INST_0]->SetVolume_Deformation(geometry, config, false, true);
  surface_movement[ZONE_0]->CopyBoundary(geometry, config);
  // project sensitivities (surface) on design variables
  SetProjection_AD(geometry, config, surface_movement[ZONE_0] , Gradient);
  config->SetKind_SU2(SU2_CFD); // set SU2_CFD as solver
}

void COneShotFluidDriver::SetAdj_ConstrFunction(su2double *seeding){

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    if (rank == MASTER_NODE){
      SU2_TYPE::SetDerivative(ConstrFunc[iConstr], SU2_TYPE::GetValue(seeding[iConstr]));
    } else {
      SU2_TYPE::SetDerivative(ConstrFunc[iConstr], 0.0);
    }
  }

}

void COneShotFluidDriver::SetAdj_ConstrFunction_Zero(){
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    SU2_TYPE::SetDerivative(ConstrFunc[iConstr], 0.0);
  }
}

void COneShotFluidDriver::SetConstrFunction(){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  su2double FunctionValue = 0.0;

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFunc[iConstr] = 0.0;

    FunctionValue = solver[FLOW_SOL]->Evaluate_ConstrFunc(config, iConstr);

    /*--- Flip sign on GEQ constraint (we want descent on (Target-FuncVal)) ---*/
    if(config->GetKind_ConstrFuncType(iConstr) == GEQ_CONSTR) {
      ConstrFunc[iConstr] = config->GetConstraintScale(iConstr)*(config->GetConstraintTarget(iConstr) - FunctionValue);
    }
    else {
      ConstrFunc[iConstr] = config->GetConstraintScale(iConstr)*(FunctionValue - config->GetConstraintTarget(iConstr));
    }

    if (rank == MASTER_NODE){
      AD::RegisterOutput(ConstrFunc[iConstr]);
    }
  }
}

void COneShotFluidDriver::StoreOldLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda_Old[iConstr] = Lambda[iConstr];
    Lambda_Store_Old[iConstr] = Lambda_Store[iConstr];
  }
}

void COneShotFluidDriver::LoadOldLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda[iConstr] = Lambda_Old[iConstr];
    Lambda_Store[iConstr] = Lambda_Store_Old[iConstr];
  }
}

void COneShotFluidDriver::UpdateLambda(su2double stepsize){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    su2double helper = 0.0;
    const su2double gamma = config->GetOneShotGamma(iConstr);
    const su2double dh = ConstrFunc[iConstr]-ConstrFunc_Store[iConstr];
    const su2double hdh = ConstrFunc_Store[iConstr]*dh;
    // const bool active = (ConstrFunc_Store[iConstr] - Lambda_Old[iConstr]/gamma > 0.);
    const bool active = (ConstrFunc_Store[iConstr] > 0.);

    // /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/
    // helper = 0.0;
    // for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
    //   if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (ConstrFunc_Store[iConstr] + Lambda_Old[iConstr]/gamma > 0.)) {
    //     helper += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
    //   }
    //   else {
    //     helper -= BCheck_Inv[iConstr][jConstr]*Lambda_Old[jConstr]/gamma;
    //     // helper -= Lambda_Old[iConstr]/stepsize;
    //     // break;
    //   }
    // }
    // Lambda[iConstr] = Lambda_Old[iConstr] + helper*stepsize*config->GetMultiplierScale(iConstr);

    // /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/
    // helper = 0.0;
    // for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
    //   helper += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
    // }
    // Lambda[iConstr] = Lambda_Store[iConstr];
    // if(config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR && ConstrFunc_Store[iConstr] + Lambda_Old[iConstr]/gamma <= 0.) {
    //   Lambda[iConstr] = 0.;
    // }
    // /*--- Only update if constraint violation improves ---*/
    // // else if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (ConstrFunc_Store[iConstr]*(ConstrFunc[iConstr]-ConstrFunc_Store[iConstr]) < 0.)) ||
    // //         (ConstrFunc[iConstr] - ConstrFunc_Store[iConstr] < 0.)) {
    // else if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (ConstrFunc_Store[iConstr]*(ConstrFunc[iConstr]-ConstrFunc_Store[iConstr]) < 0.)) ||
    //         (ConstrFunc[iConstr] - ConstrFunc_Store[iConstr] < 0.)) {
    //   Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    //   // Lambda[iConstr] = Lambda_Old[iConstr] + helper*stepsize*config->GetMultiplierScale(iConstr);
    //   // Lambda_Store[iConstr] = Lambda[iConstr];
    // }
    // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (ConstrFunc_Store[iConstr]*(ConstrFunc[iConstr]-ConstrFunc_Store[iConstr]) < 0.)) ||
    //    (ConstrFunc[iConstr] - ConstrFunc_Store[iConstr] < 0.)) {
    //   Lambda_Store[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    // }
    // // Lambda_Store[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);

    /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/
    for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
      helper += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
    }
    // Lambda[iConstr] = Lambda_Store[iConstr];
    /*--- Only update if constraint violation improves ---*/
    // if((config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) && (!active) && (dh <= 0.)) {
    if((config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) && (!active)) {
      Lambda[iConstr] -= stepsize*Lambda_Old[iConstr]*config->GetMultiplierScale(iConstr);
      Lambda_Store[iConstr] -= stepsize*Lambda_Store[iConstr];
      // Lambda_Store[iConstr] -= stepsize*Lambda_Store[iConstr];
      // Lambda_Store[iConstr] -= stepsize*Lambda_Store[iConstr]*config->GetMultiplierScale(iConstr);
    }
    // else if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (hdh <= 0.)) || (dh <= 0.)) {
    else {
      Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
      // Lambda_Store[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
      // Lambda_Store[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
      // Lambda[iConstr] = Lambda_Old[iConstr] + helper*stepsize*config->GetMultiplierScale(iConstr);
      // Lambda_Store[iConstr] = Lambda[iConstr];
    }
    // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (hdh <= 0.)) || (dh <= 0.)) {
    //   Lambda_Store[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    // }
    // Lambda_Store[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);

    // /*--- gamma*(h-P_I(h+mu/gamma)) ---*/
    // if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (ConstrFunc_Store[iConstr] + Lambda_Old[iConstr]/gamma > 0.)) {
    //   Lambda[iConstr] = Lambda_Old[iConstr] + stepsize*gamma*ConstrFunc_Store[iConstr];
    // }
    // else {
    //   Lambda[iConstr] = 0.;
    // }
    // Lambda_Store[iConstr] += stepsize*gamma*ConstrFunc_Store[iConstr];

    if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) {
      if(Lambda[iConstr]*ConstrFunc_Store[iConstr] < 0.) {
        Lambda[iConstr] = ConstrFunc_Store[iConstr];
      }
    }
    // else if(ConstrFunc_Store[iConstr] + Lambda[iConstr]/gamma <= 0.){
    //   // Lambda[iConstr] = 0.;
    //   helper = 0.0;
    //   for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
    //     helper += BCheck_Inv[iConstr][jConstr]*Lambda_Old[jConstr]/gamma;
    //   }
    //   Lambda[iConstr] = max(Lambda_Old[iConstr] - helper*stepsize*config->GetMultiplierScale(iConstr), 0.0);
    // }
    // else {
    //   Lambda_Store[iConstr] = Lambda[iConstr];
    // }
    else {
      Lambda[iConstr] = max(Lambda[iConstr], 0.);
      Lambda_Store[iConstr] = max(Lambda_Store[iConstr], 0.);
    }
  }
}

void COneShotFluidDriver::StoreLambdaGrad() {
  if(nConstr > 0) {
    unsigned short iConstr, iVar, nVar = solver[ADJFLOW_SOL]->GetnVar();
    unsigned long iPoint, nPointDomain = geometry->GetnPointDomain();
    const su2double beta = config->GetOneShotBeta();
    for (iConstr = 0; iConstr < nConstr; iConstr++) {
      const su2double gamma = config->GetOneShotGamma(iConstr);
      // const bool active = (ConstrFunc_Store[iConstr] - Lambda_Old[iConstr]/gamma > 0.);
      const bool active = (ConstrFunc_Store[iConstr] > 0.);
      su2double my_Gradient = 0.;
      // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (ConstrFunc_Store[iConstr]*(ConstrFunc[iConstr]-ConstrFunc_Store[iConstr]) < 0.)) || 
      //    ((ConstrFunc[iConstr] - ConstrFunc_Store[iConstr] < 0.) && (ConstrFunc[iConstr] + Lambda_Old[iConstr]/gamma > 0.))) {
      if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
        // my_Gradient += ConstrFunc[iConstr] + Lambda[iConstr]/gamma;
        my_Gradient += ConstrFunc[iConstr];
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          for (iVar = 0; iVar < nVar; iVar++) {
            my_Gradient += beta
                * solver[ADJFLOW_SOL]->GetConstrDerivative(iConstr, iPoint, iVar)
                * solver[ADJFLOW_SOL]->GetNodes()->GetSolution_Delta(iPoint,iVar);
          }
        }
      }
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_Gradient, &AugLagLamGrad[iConstr], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  AugLagLamGrad[iConstr] = my_Gradient;
#endif
    }
  }
}

void COneShotFluidDriver::CheckLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) {
      if(Lambda[iConstr]*ConstrFunc_Store[iConstr] < 0.) {
        Lambda[iConstr] = 0.;
        for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++) {
          Lambda[iConstr] += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
        }
      }
    }
    else {
      Lambda[iConstr] = max(Lambda[iConstr], 0.);
    }
  }
}

void COneShotFluidDriver::StoreConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFunc_Store[iConstr] = ConstrFunc[iConstr];
  }
}
