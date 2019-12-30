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
  /*---------- One-shot works on all design variables - find total number of design variables ---------*/
  OneShotIter = 0;
  nDV_Total = 0;
  for (unsigned short iDV = 0; iDV  < config->GetnDV(); iDV++){
    for (unsigned short iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++){
      nDV_Total++;
    }
  }

  nConstr = config->GetnConstr();

  Gradient = new su2double[nDV_Total];
  ShiftLagGrad = new su2double[nDV_Total];
  ShiftLagGradOld = new su2double[nDV_Total];

  AugLagGrad = new su2double[nDV_Total];
  AugLagGradAlpha = new su2double[nDV_Total];
  AugLagGradBeta = new su2double[nDV_Total];
  AugLagGradGamma = new su2double[nDV_Total];
  AugLagGradOld = new su2double[nDV_Total];

  DesignVarUpdate = new su2double[nDV_Total];
  DesignVar = new su2double[nDV_Total];
  SearchDirection = new su2double[nDV_Total];
  ActiveSetDV = new bool[nDV_Total];

  BFGSInv = new su2double*[nDV_Total];

  if(nConstr > 0){
    ConstrFunc = new su2double[nConstr];
    ConstrFuncStore = new su2double[nConstr];
    ConstrFuncOld = new su2double[nConstr];
    Lambda = new su2double[nConstr];
    LambdaOld = new su2double[nConstr];
    LambdaStore = new su2double[nConstr];
    LambdaTilde = new su2double[nConstr];
    LambdaTildeOld = new su2double[nConstr];
    LambdaTildeStore = new su2double[nConstr];
    AugLagLamGrad = new su2double[nConstr];
    BCheckInv = new su2double*[nConstr];
  }

  for (unsigned short iDV = 0; iDV  < nDV_Total; iDV++){
    Gradient[iDV] = 0.0;
    ShiftLagGrad[iDV] = 0.0;
    ShiftLagGradOld[iDV] = 0.0;
    AugLagGrad[iDV] = 0.0;
    AugLagGradAlpha[iDV] = 0.0;
    AugLagGradBeta[iDV] = 0.0;
    AugLagGradGamma[iDV] = 0.0;
    AugLagGradOld[iDV] = 0.0;
    DesignVarUpdate[iDV] = 0.0;
    DesignVar[iDV] = 0.0;
    SearchDirection[iDV] = 0.0;
    ActiveSetDV[iDV]=false;
    BFGSInv[iDV] = new su2double[nDV_Total];
    for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
      BFGSInv[iDV][jDV] = 0.0;
    }
    BFGSInv[iDV][iDV] = 1.0;
  }

  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    ConstrFunc[iConstr] = -1.0E-16;
    ConstrFuncStore[iConstr] = -1.0E-16;
    ConstrFuncOld[iConstr] = -1.0E-16;
    Lambda[iConstr] = 0.0;
    LambdaOld[iConstr] = 0.0;
    LambdaStore[iConstr] = 0.0;
    LambdaTilde[iConstr] = 0.0;
    LambdaTildeOld[iConstr] = 0.0;
    LambdaTildeStore[iConstr] = 0.0;
    AugLagLamGrad[iConstr] = 0.0;
    BCheckInv[iConstr] = new su2double[nConstr];
    for (unsigned short jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheckInv[iConstr][jConstr] = 0.0;
    }
    BCheckInv[iConstr][iConstr] = config->GetOneShotGamma();
  }
  BCheckNorm = 1./config->GetOneShotGamma();

  /*----- calculate values for bound projection algorithm -------*/
  lb=-config->GetBound()*config->GetDesignScale();
  ub=config->GetBound()*config->GetDesignScale();
  epsilon=(ub-lb)/2.0;
  stepsize0 = 1.0;

  /*---- calculate line search parameter ----*/
  CWolfeOne= 1E-4;
  // CWolfeTwo= 1.0-CWolfeOne;
  CWolfeTwo= 9.0E-1;

  grid_movement[ZONE_0][INST_0] = new CVolumetricMovement(geometry, config);
  surface_movement[ZONE_0]      = new CSurfaceMovement();

  /*--- Store some initial values ---*/
  ObjFunc        = 1.0;
  ObjFuncStore  = 1.0;
  Lagrangian     = 1.0;
  LagrangianOld = -1.0;
  GradDotDir     = 1.0;
  GradDotDirOld  = 1.0;

}

COneShotFluidDriver::~COneShotFluidDriver(void){

  /*----- free allocated memory -------*/
  for (unsigned short iDV = 0; iDV  < nDV_Total; iDV++){
    delete [] BFGSInv[iDV];
  }
  delete [] BFGSInv;
  delete [] Gradient;
  delete [] ShiftLagGrad;
  delete [] ShiftLagGradOld;

  delete [] AugLagGrad;
  delete [] AugLagGradAlpha;
  delete [] AugLagGradBeta;
  delete [] AugLagGradGamma;
  delete [] AugLagGradOld;

  delete [] DesignVarUpdate;
  delete [] DesignVar;
  delete [] SearchDirection;
  delete [] ActiveSetDV;

  if(nConstr > 0){
    delete [] ConstrFunc;
    delete [] ConstrFuncStore;
    delete [] ConstrFuncOld;
    delete [] Lambda;
    delete [] LambdaOld;
    delete [] LambdaStore;
    delete [] LambdaTilde;
    delete [] LambdaTildeOld;
    delete [] LambdaTildeStore;
    delete [] AugLagLamGrad;
  }

}

void COneShotFluidDriver::Preprocess(unsigned long TimeIter) {

  config->SetTimeIter(TimeIter);

}


void COneShotFluidDriver::Run(){

  unsigned long nOneShotIter = config->GetnInner_Iter();;

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
        if((OneShotIter <= config->GetOneShotStart()) || 
           ((config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) && 
            ((ConstrFunc[iConstr] > 0.) || (Lambda[iConstr] < 0.)))) {
          StopCalc = false;
          break;
        }
      }
    }

    if(StopCalc) break;

  }

}

void COneShotFluidDriver::RunOneShot(){

  su2double stepsize = stepsize0, stepsizel = 0., stepsizer = 1., tol = config->GetOneShotSearchTol();
  unsigned short ArmijoIter = 0, nArmijoIter = config->GetOneShotSearchIter();
  unsigned short ArmijoFlag = 1;
  bool bool_tol = false;

  /*--- Store the old solution and the old design for line search ---*/
  // solver[ADJFLOW_SOL]->SetOldStoreSolution();
  solver[ADJFLOW_SOL]->SetStoreSolution();
  solver[ADJFLOW_SOL]->SetMeshPointsOld(config, geometry);

  /*--- Perform line search on the design ---*/
  do {

    if((OneShotIter > config->GetOneShotStart()) && (OneShotIter < config->GetOneShotStop())){      

      if(ArmijoIter > 0){
        /*--- Parabolic backtracking ---*/
        su2double stepsize_tmp = UpdateStepSizeQuadratic();
        if(ArmijoFlag == 1) {
          // stepsizer = stepsize;
          stepsize = UpdateStepSizeBound(stepsize_tmp, stepsize/10., stepsize/2.);
          // stepsize  = 0.5*(stepsizel+stepsize);
        }
        // else if(ArmijoFlag == 2) {
        //   // stepsize = min(UpdateStepSizeBound(stepsize_tmp, stepsize*1.5, stepsize*7.5), 1.0);
        //   if(ArmijoIter == 1) {
        //     ArmijoFlag = 0;
        //     break;
        //   }
        //   else {
        //     stepsizel = stepsize;
        //     stepsize  = 0.5*(stepsize+stepsizer);
        //   }
        // }
        if(stepsize < tol) {
          stepsize = tol;
          bool_tol = true;
        }

        /*---Load the old design and solution for line search---*/
        solver[ADJFLOW_SOL]->LoadMeshPointsOld(config, geometry);
        solver[ADJFLOW_SOL]->LoadSolution();

        /*--- Preprocess to recompute primitive variables ---*/
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
      }
      // else {
      //   StoreLambdaGrad();
      // }
      // else{
      //   // UpdateLambda(1.0);
      //   // UpdateLambda(stepsize);
      //   ComputeDesignVarUpdate(stepsize);
      //   StoreGradDotDir();
      //   if(GradDotDir >= 0) {
      //     stepsize = 0.0;
      //     bool_tol = true;
      //     ComputeDesignVarUpdate(0.0);
      //     PrimalDualStep();
      //     solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);
      //     StoreObjFunction();
      //     StoreConstrFunction();
      //     // UpdateLambda(1.0);
      //     ArmijoIter = 1;
      //     break;
      //   }
      // }

      // solver[ADJFLOW_SOL]->LoadSaveSolution();
      // --- Preprocess to recompute primitive variables ---
      // solver[FLOW_SOL]->Preprocessing(geometry, solver, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      if ((ArmijoIter < nArmijoIter-1) && (!bool_tol)) {
        ComputeDesignVarUpdate(stepsize);
        config->SetKind_SU2(SU2_DEF); // set SU2_DEF as the solver
        SurfaceDeformation(surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0]);
        config->SetKind_SU2(SU2_CFD); // set SU2_CFD as the solver

        /*--- Evaluate the objective at the old solution, new design ---*/
        ComputeFunctionals();
        StoreObjFunction();
        StoreConstrFunction();

        // LoadOldLambda();
        // UpdateLambda(stepsize);
      }
      else {
        stepsize = 0.0;
        grid_movement[ZONE_0][INST_0]->UpdateDualGrid(geometry, config);
        ComputeDesignVarUpdate(0.0);

        /*--- Evaluate the objective at the old solution, new design ---*/
        ComputeFunctionals();
        StoreObjFunction();
        StoreConstrFunction();

        // LoadOldLambda();
        // UpdateLambda(1.0);
      }

      LoadOldLambda();
      // UpdateLambda(1.0);
      UpdateLambda(stepsize);

    }

    /*--- Do a primal and adjoint update ---*/
    PrimalDualStep();
    solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);

    /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
    if ((OneShotIter > config->GetOneShotStart()) && 
        (OneShotIter < config->GetOneShotStop())  &&
        ((ArmijoIter < nArmijoIter-1) && (!bool_tol))) {
      /*--- Compute and store GradL dot p ---*/
      StoreLambdaGrad();
      StoreGradDotDir(true);
      CalculateLagrangian();
      ArmijoFlag = CheckArmijo(true);
    }

    ArmijoIter++;

  } while((OneShotIter > config->GetOneShotStart()) && 
          (OneShotIter < config->GetOneShotStop())  &&
          (ArmijoFlag != 0) && (ArmijoIter < nArmijoIter) && (!bool_tol));

  /*--- Save solution ---*/
  solver[ADJFLOW_SOL]->SetSaveSolution();

  /*--- Store number of search iterations ---*/
  solver[ADJFLOW_SOL]->SetArmijoIter(ArmijoIter);

  /*--- Perform line search on the multiplier ---*/
  // if((OneShotIter > config->GetOneShotStart()) && 
  //    (OneShotIter < config->GetOneShotStop())  && 
  //    (ArmijoFlag != 0)) {

  //   bool bool_tol_feas = false;
  //   unsigned short ArmijoIterFeas = 0, ArmijoFlagFeas = 1;
  //   su2double stepsizefeas = 1.0, stepsizel = 0.0, stepsizer = 1.0;

  //   /*--- Store GradL ---*/
  //   StoreLambdaGrad();

  //   do {

  //     if(ArmijoIterFeas > 0){
  //       /*--- Parabolic backtracking ---*/
  //       su2double stepsize_tmp = UpdateStepSizeQuadratic();
  //       if(ArmijoFlagFeas == 1) {
  //         // stepsizer    = stepsizefeas;
  //         stepsizefeas = UpdateStepSizeBound(stepsize_tmp, stepsizefeas/10., stepsizefeas/2.);
  //         // stepsizefeas = 0.5*(stepsizel+stepsizefeas);
  //       }
  //       // else if(ArmijoFlagFeas == 2) {
  //       //   // stepsizefeas = min(UpdateStepSizeBound(stepsize_tmp, stepsizefeas*1.5, stepsizefeas*7.5), 1.0);
  //       //   if(ArmijoIter == 1) {
  //       //     ArmijoFlag = 0;
  //       //     break;
  //       //   }
  //       //   else {
  //       //     stepsizel = stepsize;
  //       //     stepsize  = 0.5*(stepsize+stepsizer);
  //       //   }
  //       // }
  //       if(stepsizefeas < tol) {
  //         stepsizefeas  = 0.0;
  //         bool_tol_feas = true;
  //       }

  //     }

  //     LoadOldLambda();
  //     // UpdateLambda(1.0);
  //     UpdateLambda(stepsizefeas);

  //     /*--- Compute and store GradL dot p ---*/
  //     StoreGradDotDir(false);

  //     /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
  //     if ((ArmijoIterFeas < nArmijoIter-1) && (!bool_tol_feas)) {
  //       CalculateLagrangian();
  //       ArmijoFlagFeas = CheckArmijo(false);
  //     }

  //     ArmijoIterFeas++;

  //   } while((ArmijoFlagFeas != 0) && (ArmijoIterFeas < nArmijoIter) && (!bool_tol_feas));
  // }

  /*--- Store FFD info in file ---*/
  if (((config->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) ||
       (config->GetDesign_Variable(0) == FFD_CONTROL_POINT))   &&
      (OneShotIter > config->GetOneShotStart())                && 
      (OneShotIter < config->GetOneShotStop())                 &&
       ((!config->GetZeroStep()) || (ArmijoFlag == 0))) {
    surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container[ZONE_0][INST_0], config_container, false);
    config->SetMesh_FileName(config->GetMesh_Out_FileName());
  }

  /*--- Compute alpha, beta, gamma at first one-shot iteration, or recompute if line search failed ---*/
  if(OneShotIter > 0) solver[ADJFLOW_SOL]->CalculateRhoTheta(config);
  if(OneShotIter == config->GetOneShotStart()) {
    solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    if(nConstr > 0) ComputePreconditioner();
    solver[ADJFLOW_SOL]->CalculateGamma(config, BCheckNorm, ConstrFunc, Lambda);
  }
  else if((OneShotIter > config->GetOneShotStart()) && 
          // (OneShotIter < config->GetOneShotStop())){
          (OneShotIter < config->GetOneShotStop()) &&
          (ArmijoFlag != 0)){
    solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    solver[ADJFLOW_SOL]->CalculateGamma(config, BCheckNorm, ConstrFunc, Lambda);

    /*--- Recalculate Lagrangian and gradient with new alpha/beta/gamma ---*/
    SetAugLagGrad(TOTAL_AUGMENTED_OLD);

  }

  /*--- Store the multiplier and constraint function, then recalculate Lagrangian for next iteration ---*/
  StoreObjFunction();
  StoreConstrFunction();
  CheckLambda();
  CalculateLagrangian();
  StoreOldLambda();
  StoreOldConstrFunction();

  /*--- Store Deltay and DeltaBary ---*/
  solver[ADJFLOW_SOL]->SetStoreSolutionDelta();

  if((OneShotIter >= config->GetOneShotStart()) && (OneShotIter < config->GetOneShotStop())){

    /*--- Update design variable ---*/
    UpdateDesignVar();

    /*--- N_u ---*/
    solver[ADJFLOW_SOL]->SetSensitivityShiftedLagrangian(geometry);
    // solver[ADJFLOW_SOL]->SetSaveSolution();
    // solver[ADJFLOW_SOL]->LoadSolution();

    if(nConstr > 0) ComputePreconditioner();

    /*--- Gamma*h^T*h_u ---*/
    if(nConstr > 0) {
      ComputeGammaTerm();
      solver[ADJFLOW_SOL]->SetSensitivityLagrangian(geometry, GAMMA_TERM);
      // solver[ADJFLOW_SOL]->LoadSolution();
    }

    /*--- Alpha*Deltay^T*G_u ---*/
    ComputeAlphaTerm();
    solver[ADJFLOW_SOL]->SetSensitivityLagrangian(geometry, ALPHA_TERM);
    // solver[ADJFLOW_SOL]->LoadSolution();

    /*--- Beta*DeltaBary^T*N_yu ---*/
    ComputeBetaTerm();
    // solver[ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry, config);
    solver[ADJFLOW_SOL]->SetSensitivityLagrangian(geometry, BETA_TERM);
    solver[ADJFLOW_SOL]->LoadSaveSolution();

    /*--- Projection of the gradient N_u ---*/
    ProjectMeshSensitivities();

    /*--- Projection of the gradient L_u---*/
    SetAugLagGrad(TOTAL_AUGMENTED);

    /*--- Use N_u to compute the active set (bound constraints) ---*/
    ComputeActiveSet(stepsize);

    /*--- Do a BFGS update to approximate the inverse preconditioner ---*/
    if(OneShotIter > config->GetOneShotStart()) BFGSUpdate(config);

    /*--- Compute the search direction for the line search procedure ---*/
    ComputeSearchDirection();

    StoreLagrangianInformation();
  }

  /*--- Modifiy initial line search guess based on success of line search ---*/
  if(OneShotIter > config->GetOneShotStart()) {
    // if((!bool_tol) && (ArmijoIter < nArmijoIter) && (stepsize < stepsize0/2.0)) {
    //   stepsize0 = max(10.0*tol, stepsize0/2.0);
    //   // stepsize0 = stepsize;
    // }
    // else {
    //   stepsize0 = min(1.0, stepsize0*2.0);
    //   // stepsize0 = stepsize;
    // }
    // else if(((!bool_tol) && (ArmijoIter < nArmijoIter)) || (ArmijoFlag == 2)) {
    // // else {
    //   stepsize0 = min(1.0, stepsize0*2.0);
    // }

    // if((!bool_tol) && (ArmijoIter < nArmijoIter)) {
    //   StoreOldGradDotDir();
    //   ComputeDesignVarUpdate(1.0);
    //   StoreGradDotDir();
    //   if(GradDotDirOld < 0 && GradDotDir < 0) {
    //     stepsize0 = max(10.0*tol, min(1.0, 1.01*GradDotDirOld/GradDotDir));
    //   }
    //   else{
    //     // stepsize0 = min(1.0, 2.0*stepsize0);
    //     stepsize0 = 1.0;
    //   }
    // }
    // else{
    //   // stepsize0 = min(1.0, 2.0*stepsize0);
    //   stepsize0 = 1.0;
    // }
  }

}

void COneShotFluidDriver::PrimalDualStep(){

  /*--- Note: Unsteady cases not applicable to the one-shot method yet! ---*/

  // SetRecording(NONE);
  // solver[ADJFLOW_SOL]->LoadSolution();
  SetRecording(COMBINED);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);


  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma();
    // const bool active = (ConstrFuncStore[iConstr] + LambdaOld[iConstr]/gamma > 0.);
    const bool active = (ConstrFuncStore[iConstr] > 0.);
    const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);
    if((eqconstr) || (active)) {
      seeding[iConstr] = Lambda[iConstr];
    }
    else {
      seeding[iConstr] = 0.;
    }
  }
  SetAdj_ConstrFunction(seeding);

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

  delete [] seeding;

}

void COneShotFluidDriver::SetRecording(unsigned short kind_recording){
  unsigned long InnerIter = config->GetInnerIter();

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the flow iteration (conservative variables and node coordinates) --- */

  if (kind_recording != NONE){

    AD::StartRecording();

    if ((rank == MASTER_NODE) && (kind_recording == MainVariables) && (InnerIter == 0)) {
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

void COneShotFluidDriver::ComputeFunctionals(){
  solver[FLOW_SOL]->Pressure_Forces(geometry, config);
  solver[FLOW_SOL]->Momentum_Forces(geometry, config);
  solver[FLOW_SOL]->Friction_Forces(geometry, config);
                  
  if((config->GetBuffet_Monitoring()) || (config->GetKind_ObjFunc() == BUFFET_SENSOR)){
    solver[FLOW_SOL]->Buffet_Monitoring(geometry, config);
  }

  ObjFunc = 0.0;
  
  direct_output->SetHistory_Output(geometry, solver, config,
                                   config->GetTimeIter(),
                                   config->GetOuterIter(),
                                   config->GetInnerIter());

  /*--- Specific scalar objective functions ---*/

  solver[FLOW_SOL]->SetTotal_ComboObj(0.0);

  /*--- Surface based obj. function ---*/

  solver[FLOW_SOL]->Evaluate_ObjFunc(config);
  ObjFunc += solver[FLOW_SOL]->GetTotal_ComboObj();

  su2double FunctionValue = 0.0;

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const bool geqconstr = (config->GetKind_ConstrFuncType(iConstr) == GEQ_CONSTR);
    ConstrFunc[iConstr] = 0.0;

    FunctionValue = solver[FLOW_SOL]->Evaluate_ConstrFunc(config, iConstr);

    /*--- Flip sign on GEQ constraint (we want descent on (Target-FuncVal)) ---*/
    if(geqconstr) {
      ConstrFunc[iConstr] = config->GetConstraintScale(iConstr)*(config->GetConstraintTarget(iConstr) - FunctionValue);
    }
    else {
      ConstrFunc[iConstr] = config->GetConstraintScale(iConstr)*(FunctionValue - config->GetConstraintTarget(iConstr));
    }
  }
}

void COneShotFluidDriver::SetProjection_AD(CSurfaceMovement *surface_movement){

  su2double DV_Value, *VarCoord, Sensitivity, my_Gradient, localGradient;
  unsigned short nMarker, nDim, nDV, nDV_Value;
  unsigned long nVertex, jPoint;

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

  for (unsigned short iDV = 0; iDV < nDV; iDV++){

    nDV_Value =  config->GetnDV_Value(iDV);

    for (unsigned short iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

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

  for(short kind_gradient = 4; kind_gradient >= 0; kind_gradient--) {

    if(kind_gradient < 3) solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry, kind_gradient);
    else                  solver[ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry);

    for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
      visited[iPoint] = false;
    }

    /*--- Initialize the derivatives of the output of the surface deformation routine
     * with the discrete adjoints from the CFD solution ---*/

    for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES) {
        nVertex = geometry->nVertex[iMarker];
        for (unsigned long iVertex = 0; iVertex <nVertex; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (!visited[jPoint]){
            VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();

            for (unsigned short iDim = 0; iDim < nDim; iDim++){
              Sensitivity = geometry->GetSensitivity(jPoint, iDim);
              SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
            }
            visited[jPoint] = true;
          }
        }
      }
    }

    /*--- Compute derivatives and extract gradient ---*/

    AD::ComputeAdjoint();

    unsigned long nDV_Count = 0;
    for (unsigned short iDV = 0; iDV  < nDV; iDV++){
      nDV_Value =  config->GetnDV_Value(iDV);

      for (unsigned short iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
        DV_Value = config->GetDV_Value(iDV, iDV_Value);
        my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
        AD::ResetInput(DV_Value);

  #ifdef HAVE_MPI
      SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #else
        localGradient = my_Gradient;
  #endif

        Gradient[nDV_Count] = localGradient;
        nDV_Count++;
      }
    }

    if(kind_gradient < 3) SetAugLagGrad(kind_gradient);
    else                  SetShiftLagGrad();
    
    AD::ClearAdjoints();
  }

  delete [] visited;

  AD::Reset();
}

void COneShotFluidDriver::SurfaceDeformation(CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement){

  unsigned short nDV_Value;
  bool allmoving=true;
  unsigned long nDV_Count = 0;

  for (unsigned short iDV = 0; iDV < config->GetnDV(); iDV++) {
    nDV_Value =  config->GetnDV_Value(iDV);

    for (unsigned short iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      config->SetDV_Value(iDV,iDV_Value, DesignVarUpdate[nDV_Count]/config->GetDesignScale());
      nDV_Count++;
    }
  }

  /*--- Surface grid deformation using design variables ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- For scale, translation and rotation if all boundaries are moving they are set via volume method
   * Otherwise, the surface deformation has been set already in SetSurface_Deformation.  --- */
  /*--- Loop over markers, set flag to false if any are not moving ---*/
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
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

  su2double vk = 0;
  su2double normyk = 0;
  su2double normsk = 0;

  su2double* yk = new su2double[nDV_Total];
  su2double* sk = new su2double[nDV_Total];
  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
    yk[iDV] = ProjectionSet(iDV, AugLagGrad[iDV]-AugLagGradOld[iDV], false);
    sk[iDV] = ProjectionSet(iDV, DesignVarUpdate[iDV], false);
    vk += yk[iDV]*sk[iDV];
    normyk += yk[iDV]*yk[iDV];
    normsk += sk[iDV]*sk[iDV];
  }

  // if ((vk > 0) && (GradDotDir < 0)){
  if (vk > 0){
    su2double** MatA = new su2double*[nDV_Total];
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      MatA[iDV] = new su2double[nDV_Total];
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
          MatA[iDV][jDV] = 0.0;
      }
    }
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
        MatA[iDV][jDV] = ProjectionPAP(iDV,jDV,BFGSInv[iDV][jDV],false)+(1.0/vk)*sk[iDV]*sk[jDV];
        for (unsigned short kDV = 0; kDV < nDV_Total; kDV++){
          MatA[iDV][jDV] += -(1.0/vk)*sk[iDV]*ProjectionPAP(kDV,jDV,BFGSInv[kDV][jDV],false)*yk[kDV]
                            -(1.0/vk)*sk[jDV]*ProjectionPAP(iDV,kDV,BFGSInv[iDV][kDV],false)*yk[kDV];
          for (unsigned short lDV = 0; lDV < nDV_Total; lDV++){
            MatA[iDV][jDV] += (1.0/vk)*(1.0/vk)*sk[iDV]*sk[jDV]*yk[lDV]*ProjectionPAP(lDV,kDV,BFGSInv[lDV][kDV],false)*yk[kDV];
          }
        }
      }
    }
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
        BFGSInv[iDV][jDV] = MatA[iDV][jDV];
      }
    }
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      delete [] MatA[iDV];
    }
    delete [] MatA;
  }else{
    // solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    // solver[ADJFLOW_SOL]->CalculateGamma(config, BCheckNorm, ConstrFunc, Lambda);

    // /*--- Recalculate gradient with new alpha/beta/gamma ---*/
    // SetAugLagGrad(TOTAL_AUGMENTED);
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
        BFGSInv[iDV][jDV] = 0.0;
      }
      BFGSInv[iDV][iDV] = 1.0;
    }
  }
  delete [] yk;
  delete [] sk;
}

unsigned short COneShotFluidDriver::CheckArmijo(bool designing){
  su2double admissible_step = 0.0, admissible_step_new = 0.0;

  if(designing) {
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      /*--- ShiftLagGrad is the gradient at the old iterate. ---*/
      // admissible_step += DesignVarUpdate[iDV]*ShiftLagGradOld[iDV];
      /*--- AugLagGrad is the gradient at the old iterate. ---*/
      admissible_step += DesignVarUpdate[iDV]*AugLagGrad[iDV];
    }
  }
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    // const bool active = (ConstrFuncStore[iConstr] > 0.);
    // const su2double gamma = config->GetOneShotGamma();
    // if(active) {
    //   admissible_step += (Lambda[iConstr]-LambdaOld[iConstr])*ConstrFuncOld[iConstr];
    // }
    // else {
    //   admissible_step -= (Lambda[iConstr]-LambdaOld[iConstr])*Lambda[iConstr]/gamma;
    // }
    /*--- AugLagLamGrad is the gradient at the old iterate. ---*/
    admissible_step -= (Lambda[iConstr]-LambdaOld[iConstr])*AugLagLamGrad[iConstr];
  }
  
  /*--- Return 0 if satisfied, 1 if 1st condition not satisfied, 2 if 2nd condition not satisfied ---*/
  if (Lagrangian > LagrangianOld - CWolfeOne*abs(admissible_step)) {
  // if (Lagrangian > LagrangianOld + CWolfeOne*admissible_step) {
    return 1;
  }
  // else if (abs(admissible_step_new) > CWolfeTwo*abs(admissible_step)) {
  // else if (Lagrangian < LagrangianOld - CWolfeTwo*abs(admissible_step)) {
  //   return 2;
  // }
  else {
    return 0;
  }
}

void COneShotFluidDriver::StoreGradDotDir(bool designing){

  GradDotDir = 0.0;

  if(designing) {
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      /*--- ShiftLagGrad is the gradient at the old iterate. ---*/
      // GradDotDir += DesignVarUpdate[iDV]*ShiftLagGradOld[iDV];
      /*--- AugLagGrad is the gradient at the old iterate. ---*/
      GradDotDir += DesignVarUpdate[iDV]*AugLagGrad[iDV];
    }
  }
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    // const bool active = (ConstrFuncStore[iConstr] > 0.);
    // const su2double gamma = config->GetOneShotGamma();
    // if(active) {
    //   GradDotDir += (Lambda[iConstr]-LambdaOld[iConstr])*ConstrFuncOld[iConstr];
    // }
    // else {
    //   GradDotDir -= (Lambda[iConstr]-LambdaOld[iConstr])*Lambda[iConstr]/gamma;
    // }
    /*--- AugLagLamGrad is the gradient at the old iterate. ---*/
    GradDotDir -= (Lambda[iConstr]-LambdaOld[iConstr])*AugLagLamGrad[iConstr];
  }
  GradDotDir = -abs(GradDotDir);
}

void COneShotFluidDriver::StoreOldGradDotDir(){
  GradDotDirOld = GradDotDir;
}

su2double COneShotFluidDriver::UpdateStepSizeQuadratic(){
  su2double step = -GradDotDir/(2.*(Lagrangian - LagrangianOld - GradDotDir));
  return step;
}

su2double COneShotFluidDriver::UpdateStepSizeBound(su2double stepsize, su2double a, su2double b){
  if(stepsize < a) {return a;}
  else if(stepsize > b) {return b;}
  else {return stepsize;}
}


void COneShotFluidDriver::ComputeDesignVarUpdate(su2double stepsize){
  for (unsigned short iDV=0;iDV<nDV_Total;iDV++){
    DesignVarUpdate[iDV]=BoundProjection(DesignVar[iDV]+stepsize*SearchDirection[iDV]*config->GetDesignScale())-DesignVar[iDV];
  }
}

void COneShotFluidDriver::ComputeSearchDirection(){
  for (unsigned short iDV=0;iDV<nDV_Total;iDV++){
    SearchDirection[iDV]=0.0;
    for (unsigned short jDV=0;jDV<nDV_Total;jDV++){
      SearchDirection[iDV]+= BFGSInv[iDV][jDV]*ProjectionSet(jDV,-ShiftLagGrad[jDV],false);
    }
    SearchDirection[iDV]=-ProjectionSet(iDV, ShiftLagGrad[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false);
  }
}

void COneShotFluidDriver::StoreLagrangianInformation(){
  for (unsigned short iDV=0; iDV<nDV_Total; iDV++){
    AugLagGradOld[iDV] = AugLagGrad[iDV];
    ShiftLagGradOld[iDV] = ShiftLagGrad[iDV];
  }
  LagrangianOld = Lagrangian;
}

void COneShotFluidDriver::UpdateDesignVar(){
  for (unsigned short iDV=0; iDV<nDV_Total; iDV++){
    DesignVar[iDV] += DesignVarUpdate[iDV];
  }
}

void COneShotFluidDriver::CalculateLagrangian(){
    
  Lagrangian = 0.0;
  Lagrangian += ObjFuncStore; //TODO use for BFGS either only objective function or normal Lagrangian

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma();
    const su2double helper = ConstrFuncStore[iConstr] + Lambda[iConstr]/gamma;
    // const bool active = (ConstrFuncStore[iConstr] + LambdaOld[iConstr]/gamma > 0.);
    const bool active = (ConstrFuncStore[iConstr] > 0.);
    const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);
    // const bool active = (LambdaTilde[iConstr] > 0.);
    /*--- Lagrangian += gamma/2 ||h + mu/gamma - P_I(h+mu/gamma)||^2 ---*/
    if((eqconstr) || (active)) {
      Lagrangian += gamma/2.*helper*helper - 1./(2.*gamma)*Lambda[iConstr]*Lambda[iConstr];
      // Lagrangian += gamma/2.*helper*helper;
    }
    else {
      Lagrangian -= 1./(2.*gamma)*Lambda[iConstr]*Lambda[iConstr];
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
  su2double norm = 0.0;

  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++) {
    norm += (DesignVar[iDV]-BoundProjection(DesignVar[iDV]-stepsize*ShiftLagGrad[iDV]))
          * (DesignVar[iDV]-BoundProjection(DesignVar[iDV]-stepsize*ShiftLagGrad[iDV]));
  }
  norm = sqrt(norm);
  epsilon = min(norm, (ub-lb)/2.0);
  unsigned short nActive = nDV_Total;

  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++) {
    ActiveSetDV[iDV] = false;
    if(ub-DesignVar[iDV] <= epsilon)      ActiveSetDV[iDV] = true;
    else if(DesignVar[iDV]-lb <= epsilon) ActiveSetDV[iDV] = true;
    else                                  nActive--;
  }
  solver[ADJFLOW_SOL]->SetnActiveDV(nActive);
}

void COneShotFluidDriver::SetShiftLagGrad(){
  su2double norm = 0.;
  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
    ShiftLagGrad[iDV] = Gradient[iDV];
    norm += Gradient[iDV]*Gradient[iDV];
  }
  solver[ADJFLOW_SOL]->SetShiftedLagGradNorm(sqrt(norm));
}

void COneShotFluidDriver::SetAugLagGrad(unsigned short kind){
  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
    if(kind == ALPHA_TERM) {
      AugLagGradAlpha[iDV] = Gradient[iDV];
    }
    else if(kind == BETA_TERM) {
      AugLagGradBeta[iDV] = (Gradient[iDV]-ShiftLagGrad[iDV])/config->GetFDStep();
    }
    else if(kind == GAMMA_TERM) {
      AugLagGradGamma[iDV] = Gradient[iDV];
    }
    else if(kind == TOTAL_AUGMENTED) {
      AugLagGrad[iDV] = ShiftLagGrad[iDV]
                      + AugLagGradAlpha[iDV]*config->GetOneShotAlpha()
                      + AugLagGradBeta[iDV]*config->GetOneShotBeta()
                      + AugLagGradGamma[iDV]*config->GetOneShotGamma();
    }   
    else if(kind == TOTAL_AUGMENTED_OLD) {
      AugLagGradOld[iDV] = ShiftLagGradOld[iDV]
                          + AugLagGradAlpha[iDV]*config->GetOneShotAlpha()
                          + AugLagGradBeta[iDV]*config->GetOneShotBeta()
                          + AugLagGradGamma[iDV]*config->GetOneShotGamma();
    }   
  }
}

void COneShotFluidDriver::ComputeGammaTerm(){

  /*--- Note: Not applicable for unsteady code ---*/

  /*--- Initialize the adjoint of the output variables of the iteration with difference of the solution and the solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  solver[ADJFLOW_SOL]->ResetInputs(geometry, config);
  iteration->InitializeAdjoint_Zero(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 0.0. ---*/

  SetAdj_ObjFunction_Zero();
  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma();
    // const bool active = (ConstrFuncStore[iConstr] + Lambda[iConstr]/gamma > 0.);
    const bool active = (ConstrFuncStore[iConstr] > 0.);
    const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);
    // const bool active = (LambdaTilde[iConstr] > 0.);
    if((eqconstr) || (active)) {
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

  solver[ADJFLOW_SOL]->ResetInputs(geometry, config);
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

  // SetRecording(NONE);
  solver[ADJFLOW_SOL]->ResetInputs(geometry, config);
  solver[ADJFLOW_SOL]->UpdateStateVariable(config);
  SetRecording(COMBINED);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma();
    // const bool active = (ConstrFuncStore[iConstr] + Lambda[iConstr]/gamma > 0.);
    const bool active = (ConstrFuncStore[iConstr] > 0.);
    const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);
    // const bool active = (LambdaTilde[iConstr] > 0.);
    if((eqconstr) || (active)) {
      seeding[iConstr] = Lambda[iConstr];
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
  solver[ADJFLOW_SOL]->SetSensitivity(geometry, solver, config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  AD::Reset();

  delete [] seeding;

}

void COneShotFluidDriver::ComputePreconditioner(){

  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    seeding[iConstr] = 0.0;
  }
  su2double **BCheck = new su2double*[nConstr];
  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr] = new su2double[nConstr];
    for (unsigned short jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck[iConstr][jConstr] = 0.0;
    }
  }

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    // const bool active = (ConstrFunc[iConstr] > 0.);
    // if(active) {
      seeding[iConstr] = 1.0;

      solver[ADJFLOW_SOL]->ResetInputs(geometry, config);
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

      // solver[ADJFLOW_SOL]->LoadSolution();

      seeding[iConstr]=0.0;
    // }
  }

  su2double bcheck=0;
  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    // BCheck[iConstr][iConstr] = 1./config->GetOneShotGamma();
    for (unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
      BCheck[iConstr][jConstr] += config->GetOneShotBeta()*solver[ADJFLOW_SOL]->MultiplyConstrDerivative(iConstr,jConstr);
    }
  }
  if (nConstr == 1){
      const su2double gamma = config->GetOneShotGamma();
      // const bool active = (ConstrFuncStore[0] + Lambda[0]/gamma > 0.);
      const bool active = (ConstrFunc[0] > 0.);
      // const bool active = (LambdaTilde[0] > 0.);
      if(active) {
        BCheckNorm = BCheck[0][0];
        solver[ADJFLOW_SOL]->CalculateGamma(config, BCheckNorm, ConstrFunc, Lambda);
        SetAugLagGrad(TOTAL_AUGMENTED_OLD);
        // BCheckNorm = BCheck[0][0] - 1./gamma;
        BCheckInv[0][0] = 1./(BCheck[0][0]+1./config->GetOneShotGamma());
      }
      else {
        BCheckNorm = 1.01/gamma;
        BCheckInv[0][0] = gamma;
      }
  } else {
    bcheck=1./(BCheck[0][0]*BCheck[1][1]*BCheck[2][2]+BCheck[1][0]*BCheck[2][1]*BCheck[0][2]+BCheck[2][0]*BCheck[0][1]*BCheck[1][2]-BCheck[0][0]*BCheck[2][1]*BCheck[1][2]-BCheck[2][0]*BCheck[1][1]*BCheck[0][2]-BCheck[1][0]*BCheck[0][1]*BCheck[2][2]);
    BCheckInv[0][0]=bcheck*(BCheck[1][1]*BCheck[2][2]-BCheck[1][2]*BCheck[2][1]);
    BCheckInv[0][1]=bcheck*(BCheck[0][2]*BCheck[2][1]-BCheck[0][1]*BCheck[2][2]);
    BCheckInv[0][2]=bcheck*(BCheck[0][1]*BCheck[1][2]-BCheck[0][2]*BCheck[1][1]);
    BCheckInv[1][0]=bcheck*(BCheck[1][2]*BCheck[2][0]-BCheck[1][0]*BCheck[2][2]);
    BCheckInv[1][1]=bcheck*(BCheck[0][0]*BCheck[2][2]-BCheck[0][2]*BCheck[2][0]);
    BCheckInv[1][2]=bcheck*(BCheck[0][2]*BCheck[1][0]-BCheck[0][0]*BCheck[1][2]);
    BCheckInv[2][0]=bcheck*(BCheck[1][0]*BCheck[2][1]-BCheck[1][1]*BCheck[2][0]);
    BCheckInv[2][1]=bcheck*(BCheck[0][1]*BCheck[2][0]-BCheck[0][0]*BCheck[2][1]);
    BCheckInv[2][2]=bcheck*(BCheck[0][0]*BCheck[1][1]-BCheck[0][1]*BCheck[1][0]);
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
  SetProjection_AD(surface_movement[ZONE_0]);
  config->SetKind_SU2(SU2_CFD); // set SU2_CFD as solver
}

void COneShotFluidDriver::SetAdj_ConstrFunction(su2double *seeding){

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

  su2double FunctionValue = 0.0;

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const bool geqconstr = (config->GetKind_ConstrFuncType(iConstr) == GEQ_CONSTR);
    ConstrFunc[iConstr] = 0.0;

    FunctionValue = solver[FLOW_SOL]->Evaluate_ConstrFunc(config, iConstr);

    /*--- Flip sign on GEQ constraint (we want descent on (Target-FuncVal)) ---*/
    if(geqconstr) {
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

void COneShotFluidDriver::StoreLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    LambdaStore[iConstr] = Lambda[iConstr];
    LambdaTildeStore[iConstr] = LambdaTilde[iConstr];
  }
}

void COneShotFluidDriver::LoadLambdaStore(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda[iConstr] = LambdaStore[iConstr];
    LambdaTilde[iConstr] = LambdaTildeStore[iConstr];
  }
}

void COneShotFluidDriver::StoreOldLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    LambdaOld[iConstr] = Lambda[iConstr];
    LambdaTildeOld[iConstr] = LambdaTilde[iConstr];
    solver[ADJFLOW_SOL]->SetLambdaValue(iConstr, Lambda[iConstr]);
  }
}

void COneShotFluidDriver::LoadOldLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda[iConstr] = LambdaOld[iConstr];
    LambdaTilde[iConstr] = LambdaTildeOld[iConstr];
  }
}

void COneShotFluidDriver::UpdateLambda(su2double stepsize){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    su2double helper = 0.0;
    const su2double gamma = config->GetOneShotGamma();
    const bool active = (ConstrFuncStore[iConstr] > 0.);
    const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);

    /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/
    if((active) || (eqconstr)) {
      for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
        helper += BCheckInv[iConstr][jConstr]*ConstrFuncOld[jConstr];
      }
      Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    }
    else {
      Lambda[iConstr] -= stepsize*LambdaOld[iConstr];
    }

    /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/

    // if(active) Lambda[iConstr] = LambdaTilde[iConstr];

    // for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
    //   helper += BCheckInv[iConstr][jConstr]*ConstrFuncStore[jConstr];
    // }

    // Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);

    // Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);

    // /*--- gamma*(h-P_I(h+mu/gamma)) ---*/
    // if(active) Lambda[iConstr] = LambdaTilde[iConstr];
    // if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
    //   Lambda[iConstr] += stepsize*gamma*ConstrFuncOld[iConstr];
    // }
    // else {
    //   Lambda[iConstr] = 0.;
    // }
    // LambdaTilde[iConstr] += stepsize*gamma*ConstrFuncOld[iConstr];


    // if(config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) {
    //   Lambda[iConstr] = max(Lambda[iConstr], 0.);
    //   // LambdaTilde[iConstr] = max(LambdaTilde[iConstr], 0.);
    // }
  }
}

void COneShotFluidDriver::CheckLambda() {
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    su2double helper = 0.0;
    const su2double gamma = config->GetOneShotGamma();
    // const bool active = (ConstrFuncStore[iConstr] + Lambda[iConstr]/gamma > 0.);
    const bool active = (ConstrFuncStore[iConstr] > 0.);
    const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);
    // const bool active = (LambdaTildeOld[iConstr] > 0.);
    if((!eqconstr) && (active)) {
      // if(Lambda[iConstr] < 0.0) {
        // // InitializeLambdaTilde(iConstr);
        // Lambda[iConstr] = max(LambdaTilde[iConstr], 0.0);
      // }

      Lambda[iConstr] = max(Lambda[iConstr], 0.0);
    }
  }
}

void COneShotFluidDriver::StoreLambdaGrad() {
  if(nConstr > 0) {
    unsigned short nVar = solver[ADJFLOW_SOL]->GetnVar();
    unsigned long nPointDomain = geometry->GetnPointDomain();
    const su2double beta = config->GetOneShotBeta();
    for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
      const su2double gamma = config->GetOneShotGamma();
      // const bool active = (ConstrFuncOld[iConstr] + LambdaOld[iConstr]/gamma > 0.);
      const bool active = (ConstrFuncStore[iConstr] > 0.);
      const bool eqconstr = (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR);
      su2double my_Gradient = 0.;
      if((eqconstr) || (active)) {
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
          for (unsigned short iVar = 0; iVar < nVar; iVar++) {
            my_Gradient += beta
                * solver[ADJFLOW_SOL]->GetConstrDerivative(iConstr, iPoint, iVar)
                * solver[ADJFLOW_SOL]->GetNodes()->GetSolution_Delta(iPoint,iVar);
          }
        }
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_Gradient, &AugLagLamGrad[iConstr], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  AugLagLamGrad[iConstr] = my_Gradient;
#endif
        // AugLagLamGrad[iConstr] += ConstrFuncOld[iConstr] + LambdaOld[iConstr]/gamma;
        AugLagLamGrad[iConstr] += ConstrFuncOld[iConstr];
      }
      else {
        AugLagLamGrad[iConstr] = -LambdaOld[iConstr]/gamma;
      }
    }
  }
}

void COneShotFluidDriver::InitializeLambdaTilde(unsigned short iConstr) {
  unsigned short nVar = solver[ADJFLOW_SOL]->GetnVar();
  unsigned long nPointDomain = geometry->GetnPointDomain();
  const su2double beta = config->GetOneShotBeta();
  const su2double gamma = config->GetOneShotGamma();
  su2double my_Lambda = 0., Lambda_Init;
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      my_Lambda -= beta
          * solver[ADJFLOW_SOL]->GetConstrDerivative(iConstr, iPoint, iVar)
          * solver[ADJFLOW_SOL]->GetNodes()->GetSolution_DeltaStore(iPoint,iVar);
    }
  }
#ifdef HAVE_MPI
SU2_MPI::Allreduce(&my_Lambda, &Lambda_Init, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
Lambda_Init = my_Lambda;
#endif
  Lambda_Init -= ConstrFuncOld[iConstr];
  LambdaTilde[iConstr] = gamma*Lambda_Init;
  // if(config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) LambdaTilde[iConstr] = max(LambdaTilde[iConstr], 0.0);
}

void COneShotFluidDriver::StoreObjFunction(){
  ObjFuncStore = ObjFunc;
}

void COneShotFluidDriver::StoreConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFuncStore[iConstr] = ConstrFunc[iConstr];
  }
}

void COneShotFluidDriver::StoreOldConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFuncOld[iConstr] = ConstrFunc[iConstr];
  }
}
