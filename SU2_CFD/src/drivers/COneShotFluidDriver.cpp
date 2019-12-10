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
  Gradient_Old = new su2double[nDV_Total];

  ShiftLagGrad = new su2double[nDV_Total];

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
    ConstrFunc_Old = new su2double[nConstr];
    Lambda = new su2double[nConstr];
    Lambda_Old = new su2double[nConstr];
    Lambda_Store = new su2double[nConstr];
    Lambda_Tilde = new su2double[nConstr];
    Lambda_Tilde_Old = new su2double[nConstr];
    Lambda_Tilde_Store = new su2double[nConstr];
    AugLagLamGrad = new su2double[nConstr];
    BCheck_Inv = new su2double*[nConstr];
  }

  BFGS_Init = config->GetBFGSInitValue();

  for (unsigned short iDV = 0; iDV  < nDV_Total; iDV++){
    Gradient[iDV] = 0.0;
    Gradient_Old[iDV] = 0.0;
    ShiftLagGrad[iDV] = 0.0;
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
    for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
      BFGS_Inv[iDV][jDV] = 0.0;
      if (iDV==jDV) BFGS_Inv[iDV][jDV] = BFGS_Init;
    }
  }

  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    ConstrFunc[iConstr] = 0.0;
    ConstrFunc_Store[iConstr] = 0.0;
    ConstrFunc_Old[iConstr] = 0.0;
    Lambda[iConstr] = 0.0;
    Lambda_Old[iConstr] = 0.0;
    Lambda_Store[iConstr] = 0.0;
    Lambda_Tilde[iConstr] = 0.0;
    Lambda_Tilde_Old[iConstr] = 0.0;
    Lambda_Tilde_Store[iConstr] = 0.0;
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
  for (unsigned short iDV = 0; iDV  < nDV_Total; iDV++){
    delete [] BFGS_Inv[iDV];
  }
  delete [] BFGS_Inv;
  delete [] Gradient;
  delete [] Gradient_Old;
  delete [] ShiftLagGrad;

  delete [] AugLagGrad;
  delete [] AugLagGradAlpha;
  delete [] AugLagGradBeta;
  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
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
    delete [] ConstrFunc_Store;
    delete [] ConstrFunc_Old;
    delete [] Lambda;
    delete [] Lambda_Old;
    delete [] Lambda_Store;
    delete [] Lambda_Tilde;
    delete [] Lambda_Tilde_Old;
    delete [] Lambda_Tilde_Store;
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
            ((ConstrFunc[iConstr] > 0.) || (Lambda[iConstr] < 1.0E-16)))) {
          StopCalc = false;
          break;
        }
      }
    }

    if(StopCalc) break;

  }

}

void COneShotFluidDriver::RunOneShot(){

  su2double stepsize = 1.0, tol = config->GetOneShotSearchTol();
  unsigned short ArmijoIter = 0, nArmijoIter = config->GetOneShotSearchIter();
  bool bool_tol = false;
  unsigned short ALPHA_TERM = 0, BETA_TERM = 1, GAMMA_TERM = 2, TOTAL_AUGMENTED = 3, TOTAL_AUGMENTED_OLD = 4;

  /*--- Store the old solution and the old design for line search ---*/
  solver[ADJFLOW_SOL]->SetStoreSolution();
  solver[ADJFLOW_SOL]->SetMeshPointsOld(config, geometry);

  // /*--- Perform line search on just multiplier ---*/
  // if(nConstr > 0 && OneShotIter > config->GetOneShotStart() && OneShotIter < config->GetOneShotStop()) {
  //   StoreLambdaGrad();

  //   StoreObjFunction();
  //   StoreConstrFunction();

  //   /*--- Do a primal and adjoint update ---*/
  //   PrimalDualStep();
  //   solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);

  //   stepsize = 1.0;
  //   ArmijoIter = 0;
  //   bool_tol = false;
  //   do {
  //     if(ArmijoIter > 0){
  //       /*--- Parabolic backtracking ---*/
  //       su2double stepsize_tmp = UpdateStepSizeQuadratic();
  //       stepsize  = UpdateStepSizeBound(stepsize_tmp, stepsize/10., stepsize/2.);
  //       if(stepsize < tol) {
  //         stepsize = 0.;
  //         bool_tol = true;
  //       }
  //     }
  //     /*--- Compute and store GradL dot p ---*/
  //     StoreGradDotDir(false);

  //     /*--- Update constraint multiplier ---*/
  //     LoadOldLambda();
  //     UpdateLambda(stepsize);

  //     /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
  //     CalculateLagrangian();

  //     ArmijoIter++;

  //   } while((!CheckFirstWolfe(false)) && (ArmijoIter < nArmijoIter) && (!bool_tol));
  //   StoreLambda();
  //   LoadOldLambda();
  //   solver[ADJFLOW_SOL]->LoadSolution();
  // }

  /*--- Perform line search on the design ---*/
  stepsize = 1.0;
  ArmijoIter = 0;
  bool_tol = false;
  do {

    if(OneShotIter > config->GetOneShotStart() && OneShotIter < config->GetOneShotStop()){      

      if(ArmijoIter > 0){
        /*--- Parabolic backtracking ---*/
        su2double stepsize_tmp = UpdateStepSizeQuadratic();
        stepsize     = UpdateStepSizeBound(stepsize_tmp, stepsize/10., stepsize/2.);
        if(stepsize < tol) {
          stepsize = tol;
          bool_tol = true;
        }

        /*---Load the old design and solution for line search---*/
        solver[ADJFLOW_SOL]->LoadMeshPointsOld(config, geometry);
        solver[ADJFLOW_SOL]->LoadSolution();

      }
      // else {
      //   // LoadOldLambda();
      //   // UpdateLambda(1.0);
      // }

      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      if (((ArmijoIter != nArmijoIter-1) && (!bool_tol)) || (!config->GetZeroStep())) {
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
        SetObjFunction(false);
        StoreObjFunction();
        SetConstrFunction(false);
        StoreConstrFunction();

        // /*--- Update constraint multiplier ---*/
        // LoadOldLambda();
        // UpdateLambda(1.0);
        // // UpdateLambda(stepsize);
      }
      else {
        stepsize = 0.0;
        grid_movement[ZONE_0][INST_0]->UpdateDualGrid(geometry, config);
        ComputeDesignVarUpdate(0.0);

        /*--- Evaluate the objective at the old solution, new design ---*/
      
        solver[FLOW_SOL]->Pressure_Forces(geometry, config);
        solver[FLOW_SOL]->Momentum_Forces(geometry, config);
        solver[FLOW_SOL]->Friction_Forces(geometry, config);
                  
        if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
          solver[FLOW_SOL]->Buffet_Monitoring(geometry, config);
        }
        SetObjFunction(false);
        StoreObjFunction();
        SetConstrFunction(false);
        StoreConstrFunction();

        // /*--- Update constraint multiplier ---*/
        // LoadOldLambda();
        // UpdateLambda(1.0);
      }

      /*--- Update constraint multiplier ---*/
      LoadOldLambda();
      UpdateLambda(1.0);

      /*--- Compute and store GradL dot p ---*/
      StoreLambdaGrad();
      StoreGradDotDir(true);

    }

    /*--- Do a primal and adjoint update ---*/
    PrimalDualStep();
    solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);

    /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
    if (((ArmijoIter != nArmijoIter-1) && (!bool_tol)) || (!config->GetZeroStep())) CalculateLagrangian();

    ArmijoIter++;

  } while((OneShotIter > config->GetOneShotStart()) && 
          (OneShotIter < config->GetOneShotStop())  &&
          (!CheckFirstWolfe(true)) && (ArmijoIter < nArmijoIter) && (!bool_tol));

  /*--- Store number of search iterations ---*/
  solver[ADJFLOW_SOL]->SetArmijoIter(ArmijoIter);

  /*--- Load multipliers from first line search ---*/
  LoadLambdaStore();
  // UpdateLambda(1.0);

  /*--- Store FFD info in file ---*/
  if (((config->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) ||
       (config->GetDesign_Variable(0) == FFD_CONTROL_POINT))   &&
       OneShotIter > config->GetOneShotStart()                   && 
       OneShotIter < config->GetOneShotStop()                    &&
       (!config->GetZeroStep() || CheckFirstWolfe(true))) {
    surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container[ZONE_0][INST_0], config_container, false);
    config->SetMesh_FileName(config->GetMesh_Out_FileName());
  }

  /*--- Compute alpha, beta, gamma at first one-shot iteration, or recompute if line search failed ---*/
  if(OneShotIter > 0) solver[ADJFLOW_SOL]->CalculateRhoTheta(config);
  if(OneShotIter == config->GetOneShotStart()) {
    solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    solver[ADJFLOW_SOL]->SetSaveSolution();
    solver[ADJFLOW_SOL]->LoadSolution();
    if((nConstr > 0) && (!config->GetConstPrecond())) ComputePreconditioner();
    solver[ADJFLOW_SOL]->LoadSaveSolution();
    solver[ADJFLOW_SOL]->CalculateGamma(config, BCheck_Norm, ConstrFunc, Lambda);
  }
  else if(OneShotIter > config->GetOneShotStart() && 
          OneShotIter < config->GetOneShotStop()  && 
          ((!CheckFirstWolfe(true)) || (ArmijoIter > nArmijoIter-1) || (bool_tol))){
    // /*--- Perform line search on just multiplier ---*/
    // if(nConstr > 0) {
    //   StoreLambdaGrad();

    //   su2double stepsize_mu = 1.0;
    //   ArmijoIter = 0;
    //   bool_tol = false;
    //   do {
    //     if(ArmijoIter > 0){
    //       /*--- Parabolic backtracking ---*/
    //       su2double stepsize_tmp = UpdateStepSizeQuadratic();
    //       stepsize_mu  = UpdateStepSizeBound(stepsize_tmp, stepsize_mu/10., stepsize_mu/2.);
    //       if(stepsize_mu < tol) {
    //         stepsize_mu = 0.;
    //         bool_tol = true;
    //       }
    //     }
    //     /*--- Compute and store GradL dot p ---*/
    //     StoreGradDotDir(false);

    //     /*--- Update constraint multiplier ---*/
    //     LoadOldLambda();
    //     UpdateLambda(stepsize_mu);

    //     /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
    //     CalculateLagrangian();

    //     ArmijoIter++;

    //   } while((!CheckFirstWolfe(false)) && (ArmijoIter < nArmijoIter) && (!bool_tol));
    // }

    solver[ADJFLOW_SOL]->CalculateAlphaBeta(config);
    solver[ADJFLOW_SOL]->CalculateGamma(config, BCheck_Norm, ConstrFunc, Lambda);

    /*--- Recalculate Lagrangian and gradient with new Alpha, Beta, Gamma, and Lambda ---*/
    CalculateLagrangian();
    SetAugLagGrad(TOTAL_AUGMENTED_OLD);
  }

  // else if(OneShotIter > config->GetOneShotStart() && 
  //         OneShotIter < config->GetOneShotStop()){
  //   /*--- Perform line search on just multiplier ---*/
  //   if(nConstr > 0) {
  //     StoreLambdaGrad();

  //     su2double stepsize_mu = 1.0;
  //     ArmijoIter = 0;
  //     bool_tol = false;
  //     do {
  //       if(ArmijoIter > 0){
  //         /*--- Parabolic backtracking ---*/
  //         su2double stepsize_tmp = UpdateStepSizeQuadratic();
  //         stepsize_mu  = UpdateStepSizeBound(stepsize_tmp, stepsize_mu/10., stepsize_mu/2.);
  //         if(stepsize_mu < tol) {
  //           stepsize_mu = 0.;
  //           bool_tol = true;
  //         }
  //       }
  //       /*--- Compute and store GradL dot p ---*/
  //       StoreGradDotDir(false);

  //       /*--- Update constraint multiplier ---*/
  //       LoadOldLambda();
  //       UpdateLambda(stepsize_mu);

  //       /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
  //       CalculateLagrangian();

  //       ArmijoIter++;

  //     } while((!CheckFirstWolfe(false)) && (ArmijoIter < nArmijoIter) && (!bool_tol));
  //   }

  //   /*--- Recalculate Lagrangian with new Lambda ---*/
  //   CalculateLagrangian();
  // }
 
  /*--- Store the multiplier and constraint function, then recalculate Lagrangian for next iteration ---*/
  StoreOldLambda();
  StoreOldConstrFunction();
  StoreObjFunction();
  StoreConstrFunction();

  /*--- Store Deltay and DeltaBary ---*/
  solver[ADJFLOW_SOL]->SetStoreSolutionDelta();

  if(OneShotIter >= config->GetOneShotStart() && OneShotIter < config->GetOneShotStop()){

    /*--- Update design variable ---*/
    UpdateDesignVar();

    /*--- N_u ---*/
    solver[ADJFLOW_SOL]->SetSensitivityShiftedLagrangian(geometry);
    solver[ADJFLOW_SOL]->SetSaveSolution();
    solver[ADJFLOW_SOL]->LoadSolution();

    if((nConstr > 0) && (!config->GetConstPrecond())) ComputePreconditioner();

    /*--- Gamma*h^T*h_u ---*/
    if(nConstr > 0) {
      ComputeGammaTerm();
      solver[ADJFLOW_SOL]->SetSensitivityLagrangian(geometry, GAMMA_TERM);
      solver[ADJFLOW_SOL]->LoadSolution();
    }

    /*--- Alpha*Deltay^T*G_u ---*/
    ComputeAlphaTerm();
    solver[ADJFLOW_SOL]->SetSensitivityLagrangian(geometry, ALPHA_TERM);
    solver[ADJFLOW_SOL]->LoadSolution();

    /*--- Beta*DeltaBary^T*N_yu ---*/
    solver[ADJFLOW_SOL]->UpdateStateVariable(config);
    ComputeBetaTerm();
    solver[ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry, config);
    solver[ADJFLOW_SOL]->SetSensitivityLagrangian(geometry, BETA_TERM);
    solver[ADJFLOW_SOL]->LoadSaveSolution();

    /*--- Projection of the gradient N_u---*/
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

  // if(OneShotIter == config->GetOneShotStart()) {
  //   for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
  //     InitializeLambdaTilde(iConstr);
  //     // Lambda_Tilde[iConstr] = 0.0;
  //     Lambda[iConstr] = Lambda_Tilde[iConstr];
  //   }
  //   StoreLambda();
  //   StoreOldLambda();
  // }

}

void COneShotFluidDriver::PrimalDualStep(){

  /*--- Note: Unsteady cases not applicable to the one-shot method yet! ---*/

  SetRecording(COMBINED);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);


  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma(iConstr);
    const bool active = (ConstrFunc_Old[iConstr] + Lambda_Old[iConstr]/gamma > 0.);
    // const bool active = (ConstrFunc[iConstr] > 0.);
    if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
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

  SetObjFunction(true);

  SetConstrFunction(true);

  AD::StopRecording();

}

void COneShotFluidDriver::SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double* Gradient){

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

  for(unsigned short kind_gradient = 0; kind_gradient < 4; kind_gradient++) {

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

void COneShotFluidDriver::SurfaceDeformation(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement){

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
    yk[iDV] = ProjectionSet(iDV, AugLagGrad[iDV]-AugLagGrad_Old[iDV], false);
    sk[iDV] = ProjectionSet(iDV, DesignVarUpdate[iDV], false);
    vk += yk[iDV]*sk[iDV];
    normyk += yk[iDV]*yk[iDV];
    normsk += sk[iDV]*sk[iDV];
  }

  if (vk>0){
    su2double** MatA = new su2double*[nDV_Total];
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      MatA[iDV] = new su2double[nDV_Total];
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
          MatA[iDV][jDV] = 0.0;
      }
    }
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
        MatA[iDV][jDV] = ProjectionPAP(iDV,jDV,BFGS_Inv[iDV][jDV],false)+(1.0/vk)*sk[iDV]*sk[jDV];
        for (unsigned short kDV = 0; kDV < nDV_Total; kDV++){
          MatA[iDV][jDV] += -(1.0/vk)*sk[iDV]*ProjectionPAP(kDV,jDV,BFGS_Inv[kDV][jDV],false)*yk[kDV]
                            -(1.0/vk)*sk[jDV]*ProjectionPAP(iDV,kDV,BFGS_Inv[iDV][kDV],false)*yk[kDV];
          for (unsigned short lDV = 0; lDV < nDV_Total; lDV++){
            MatA[iDV][jDV] += (1.0/vk)*(1.0/vk)*sk[iDV]*sk[jDV]*yk[lDV]*ProjectionPAP(lDV,kDV,BFGS_Inv[lDV][kDV],false)*yk[kDV];
          }
        }
      }
    }
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
        BFGS_Inv[iDV][jDV] = MatA[iDV][jDV];
      }
    }
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      delete [] MatA[iDV];
    }
    delete [] MatA;
    if(config->GetBFGSInit()){
      for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
        BFGS_Init = vk/normyk;
      }
    }

  }else{
    /*--- Calculate new alpha, beta, gamma, and reset BFGS update if needed ---*/
    if(config->GetBoolBFGSReset()){
      for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
        for (unsigned short jDV = 0; jDV < nDV_Total; jDV++){
          BFGS_Inv[iDV][jDV] = 0.0;
          if(iDV==jDV){ BFGS_Inv[iDV][jDV] = ProjectionSet(iDV,BFGS_Init,false); }
        }
      }
    }
  }
  delete [] yk;
  delete [] sk;
}

bool COneShotFluidDriver::CheckFirstWolfe(bool design_update){
  su2double admissible_step = 0.0;

  if(design_update) {
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      /*--- ShiftLagGrad is the gradient at the old iterate. ---*/
      // admissible_step += DesignVarUpdate[iDV]*ShiftLagGrad[iDV];
      /*--- AugLagGrad is the gradient at the old iterate. ---*/
      admissible_step += DesignVarUpdate[iDV]*AugLagGrad[iDV];
    }
  }
  else {
    if (nConstr > 0) {
      for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
        const su2double gamma = config->GetOneShotGamma(iConstr);
        const su2double dh = ConstrFunc_Store[iConstr]-ConstrFunc_Old[iConstr];
        const su2double hdh = ConstrFunc_Old[iConstr]*dh;
        // const bool active = (ConstrFunc_Old[iConstr] + Lambda_Old[iConstr]/gamma > 0.);
        const bool active = (ConstrFunc_Store[iConstr] > 0.);
        // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (hdh <= 0.)) || 
           // ((active) && (dh <= 0.))) {
        // if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active && dh < 0.)) {
          // admissible_step -= (Lambda[iConstr]-Lambda_Old[iConstr])*ConstrFunc_Old[iConstr];
          admissible_step -= (Lambda[iConstr]-Lambda_Old[iConstr])*AugLagLamGrad[iConstr];
        // }
      }
    }
  }
  admissible_step *= CWolfeOne;

  return (Lagrangian <= Lagrangian_Old + admissible_step);
}

void COneShotFluidDriver::StoreGradDotDir(bool design_update){

  GradDotDir = 0.0;

  if(design_update) {
    for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
      /*--- ShiftLagGrad is the gradient at the old iterate. ---*/
      // GradDotDir += DesignVarUpdate[iDV]*ShiftLagGrad[iDV];
      /*--- AugLagGrad is the gradient at the old iterate. ---*/
      GradDotDir += DesignVarUpdate[iDV]*AugLagGrad[iDV];
    }
  }
  else {
    if (nConstr > 0) {
      for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
        const su2double gamma = config->GetOneShotGamma(iConstr);
        const su2double dh = ConstrFunc_Store[iConstr]-ConstrFunc_Old[iConstr];
        const su2double hdh = ConstrFunc_Old[iConstr]*dh;
        // const bool active = (ConstrFunc_Old[iConstr] + Lambda_Old[iConstr]/gamma > 0.);
        const bool active = (ConstrFunc_Store[iConstr] > 0.);
        // if(((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) && (hdh <= 0.)) || 
           // ((active) && (dh <= 0.))) {
        // if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active && dh < 0.)) {
          // GradDotDir -= (Lambda[iConstr]-Lambda_Old[iConstr])*ConstrFunc_Old[iConstr];
          GradDotDir -= (Lambda[iConstr]-Lambda_Old[iConstr])*AugLagLamGrad[iConstr];
        // }
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
  for (unsigned short iDV=0;iDV<nDV_Total;iDV++){
    DesignVarUpdate[iDV]=BoundProjection(DesignVar[iDV]+stepsize*SearchDirection[iDV]*config->GetDesignScale())-DesignVar[iDV];
  }
}

void COneShotFluidDriver::ComputeSearchDirection(){
  for (unsigned short iDV=0;iDV<nDV_Total;iDV++){
    SearchDirection[iDV]=0.0;
    for (unsigned short jDV=0;jDV<nDV_Total;jDV++){
      SearchDirection[iDV]+= BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,-ShiftLagGrad[jDV],false);
    }
    SearchDirection[iDV]=-ProjectionSet(iDV, ShiftLagGrad[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false);
  }
}

void COneShotFluidDriver::StoreLagrangianInformation(){
  for (unsigned short iDV=0; iDV<nDV_Total; iDV++){
    AugLagGrad_Old[iDV] = AugLagGrad[iDV];
  }
  Lagrangian_Old = Lagrangian;
}

void COneShotFluidDriver::UpdateDesignVar(){
  for (unsigned short iDV=0; iDV<nDV_Total; iDV++){
    DesignVar[iDV] += DesignVarUpdate[iDV];
  }
}

void COneShotFluidDriver::CalculateLagrangian(){
    
  Lagrangian = 0.0;
  Lagrangian += ObjFunc_Store; //TODO use for BFGS either only objective function or normal Lagrangian

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma(iConstr);
    const su2double helper = ConstrFunc_Store[iConstr] + Lambda[iConstr]/gamma;
    const bool active = (ConstrFunc_Store[iConstr] + Lambda[iConstr]/gamma > 0.);
    // const bool active = (ConstrFunc_Store[iConstr] > 0.);
    /*--- Lagrangian += gamma/2 ||h + mu/gamma - P_I(h+mu/gamma)||^2 ---*/
    if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
      // Lagrangian += gamma/2.*helper*helper - 1./(2.*gamma)*Lambda[iConstr]*Lambda[iConstr];
      Lagrangian += gamma/2.*helper*helper;
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
  unsigned short ALPHA_TERM = 0, BETA_TERM = 1, GAMMA_TERM = 2, TOTAL_AUGMENTED = 3, TOTAL_AUGMENTED_OLD = 4;
  for (unsigned short iDV = 0; iDV < nDV_Total; iDV++){
    if(kind == ALPHA_TERM) {
      AugLagGradAlpha[iDV] = Gradient[iDV];
    }
    else if(kind == BETA_TERM) {
      AugLagGradBeta[iDV] = Gradient[iDV];
    }
    else if(kind == GAMMA_TERM) {
      for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
        AugLagGradGamma[iDV][iConstr] = Gradient[iDV];
      }
    }
    else if(kind == TOTAL_AUGMENTED) {
      AugLagGrad[iDV] = ShiftLagGrad[iDV]
                      + AugLagGradAlpha[iDV]*config->GetOneShotAlpha()
                      + AugLagGradBeta[iDV]*config->GetOneShotBeta();
      for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
        AugLagGrad[iDV] += AugLagGradGamma[iDV][iConstr]*config->GetOneShotGamma(iConstr);
      }
    }   
    else if(kind == TOTAL_AUGMENTED_OLD) {
      AugLagGrad_Old[iDV] = ShiftLagGrad[iDV]
                          + AugLagGradAlpha[iDV]*config->GetOneShotAlpha()
                          + AugLagGradBeta[iDV]*config->GetOneShotBeta();
      for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
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
    const bool active = (ConstrFunc[iConstr] + Lambda[iConstr]/gamma > 0.);
    // const bool active = (ConstrFunc[iConstr] > 0.);
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
  su2double* seeding = new su2double[nConstr];
  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    const su2double gamma = config->GetOneShotGamma(iConstr);
    const bool active = (ConstrFunc[iConstr] + Lambda[iConstr]/gamma > 0.);
    // const bool active = (ConstrFunc[iConstr] > 0.);
    if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
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
  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr][iConstr] = 1./config->GetOneShotGamma(iConstr);
    for (unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
      BCheck[iConstr][jConstr] += config->GetOneShotBeta()*solver[ADJFLOW_SOL]->MultiplyConstrDerivative(iConstr,jConstr);
    }
  }
  if (nConstr == 1){
      const bool active = (ConstrFunc[0] + Lambda[0]/gamma > 0.);
      if(active) {
        BCheck_Norm = BCheck[0][0] - 1./config->GetOneShotGamma(0);
        BCheck_Inv[0][0] = 1./BCheck[0][0];
      }
      else {
        BCheck_Norm = 2./config->GetOneShotGamma(0);
        BCheck_Inv[0][0] = config->GetOneShotGamma(0);
      }
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

void COneShotFluidDriver::SetObjFunction(bool registering){

  bool heat         = (config->GetWeakly_Coupled_Heat());
  bool turbo        = (config->GetBoolTurbomachinery());

  ObjFunc = 0.0;
  
  direct_output->SetHistory_Output(geometry, solver, config,
                                   config->GetTimeIter(),
                                   config->GetOuterIter(),
                                   config->GetInnerIter());

  /*--- Specific scalar objective functions ---*/

  switch (config->GetKind_Solver()) {
  case DISC_ADJ_INC_EULER:       case DISC_ADJ_INC_NAVIER_STOKES:      case DISC_ADJ_INC_RANS:    
  case DISC_ADJ_EULER:           case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
  case DISC_ADJ_FEM_EULER:       case DISC_ADJ_FEM_NS:                 case DISC_ADJ_FEM_RANS:
  case ONE_SHOT_EULER:           case ONE_SHOT_NAVIER_STOKES:          case ONE_SHOT_RANS:

    solver[FLOW_SOL]->SetTotal_ComboObj(0.0);

//    if (config->GetnMarker_Analyze() != 0)
//      output->SpecialOutput_AnalyzeSurface(solver[FLOW_SOL], geometry, config, false);

//    if ((config->GetnMarker_Analyze() != 0) && compressible)
//      output->SpecialOutput_Distortion(solver[FLOW_SOL], geometry, config, false);

//    if (config->GetnMarker_NearFieldBound() != 0)
//      output->SpecialOutput_SonicBoom(solver[FLOW_SOL], geometry, config, false);

//    if (config->GetPlot_Section_Forces())
//      output->SpecialOutput_SpanLoad(solver[FLOW_SOL], geometry, config, false);

    /*--- Surface based obj. function ---*/

    solver[FLOW_SOL]->Evaluate_ObjFunc(config);
    ObjFunc += solver[FLOW_SOL]->GetTotal_ComboObj();
    if (heat){
      if (config->GetKind_ObjFunc() == TOTAL_HEATFLUX) {
        ObjFunc += solver[HEAT_SOL]->GetTotal_HeatFlux();
      }
      else if (config->GetKind_ObjFunc() == TOTAL_AVG_TEMPERATURE) {
        ObjFunc += solver[HEAT_SOL]->GetTotal_AvgTemperature();
      }
    }

    /*--- This calls to be moved to a generic framework at a next stage         ---*/
    /*--- Some things that are currently hacked into output must be reorganized ---*/
    if (turbo){

      solver[FLOW_SOL]->SetTotal_ComboObj(0.0);
      output_legacy->ComputeTurboPerformance(solver[FLOW_SOL], geometry, config);

      switch (config_container[ZONE_0]->GetKind_ObjFunc()){
      case ENTROPY_GENERATION:
        solver[FLOW_SOL]->AddTotal_ComboObj(output_legacy->GetEntropyGen(config->GetnMarker_TurboPerformance() - 1, config->GetnSpanWiseSections()));
        break;
      case FLOW_ANGLE_OUT:
        solver[FLOW_SOL]->AddTotal_ComboObj(output_legacy->GetFlowAngleOut(config->GetnMarker_TurboPerformance() - 1, config->GetnSpanWiseSections()));
        break;
      case MASS_FLOW_IN:
        solver[FLOW_SOL]->AddTotal_ComboObj(output_legacy->GetMassFlowIn(config->GetnMarker_TurboPerformance() - 1, config->GetnSpanWiseSections()));
        break;
      default:
        break;
      }

      ObjFunc = solver[FLOW_SOL]->GetTotal_ComboObj();

    }

    break;
  case DISC_ADJ_FEM:
    switch (config->GetKind_ObjFunc()){
    case REFERENCE_GEOMETRY:
        ObjFunc = solver[FEA_SOL]->GetTotal_OFRefGeom();
        break;
    case REFERENCE_NODE:
        ObjFunc = solver[FEA_SOL]->GetTotal_OFRefNode();
        break;
    case VOLUME_FRACTION:
        ObjFunc = solver[FEA_SOL]->GetTotal_OFVolFrac();
        break;
    default:
        ObjFunc = 0.0;  // If the objective function is computed in a different physical problem
        break;
    }
    break;
  }

  /*--- Scale objective for one-shot ---*/
  switch (config->GetKind_Solver()) {
    case ONE_SHOT_EULER: case ONE_SHOT_NAVIER_STOKES: case ONE_SHOT_RANS:
      ObjFunc *= config->GetObjScale();
      break;
  }

  if (rank == MASTER_NODE && registering){
    AD::RegisterOutput(ObjFunc);
  }

}

void COneShotFluidDriver::SetConstrFunction(bool registering){

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

    if (rank == MASTER_NODE && registering){
      AD::RegisterOutput(ConstrFunc[iConstr]);
    }
  }
}

void COneShotFluidDriver::StoreLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda_Store[iConstr] = Lambda[iConstr];
    Lambda_Tilde_Store[iConstr] = Lambda_Tilde[iConstr];
  }
}

void COneShotFluidDriver::LoadLambdaStore(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda[iConstr] = Lambda_Store[iConstr];
    Lambda_Tilde[iConstr] = Lambda_Tilde_Store[iConstr];
  }
}

void COneShotFluidDriver::StoreOldLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda_Old[iConstr] = Lambda[iConstr];
    Lambda_Tilde_Old[iConstr] = Lambda_Tilde[iConstr];
  }
}

void COneShotFluidDriver::LoadOldLambda(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lambda[iConstr] = Lambda_Old[iConstr];
    Lambda_Tilde[iConstr] = Lambda_Tilde_Old[iConstr];
  }
}

void COneShotFluidDriver::UpdateLambda(su2double stepsize){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    su2double helper = 0.0;
    const su2double gamma = config->GetOneShotGamma(iConstr);
    const su2double dh = ConstrFunc[iConstr]-ConstrFunc_Old[iConstr];
    const su2double hdh = ConstrFunc[iConstr]*dh;
    const bool active = (ConstrFunc_Old[iConstr] + Lambda_Old[iConstr]/gamma > 0.);
    // const bool active = (ConstrFunc_Old[iConstr] > 0.);

    // /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/
    // for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
    //   helper += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Old[jConstr];
    // }
    // if(active) Lambda[iConstr] = Lambda_Tilde[iConstr];

    // if((config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) && (!active)) {
    //   Lambda[iConstr] = 0.0;
    //   // InitializeLambdaTilde(iConstr);
    //   // Lambda_Old[iConstr] = Lambda_Tilde[iConstr];
    //   // Lambda_Tilde_Old[iConstr] = Lambda_Tilde[iConstr];
    //   // Lambda_Tilde[iConstr] = -gamma*ConstrFunc_Old[iConstr];
    // }
    // // else if ((active && dh < 0.) || (config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR && hdh < 0.)){
    // else {
    //   // InitializeLambdaTilde(iConstr);
    //   // Lambda[iConstr] = Lambda_Tilde[iConstr];
    //   Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    //   // Lambda_Tilde[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    // }
    // Lambda_Tilde[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);

    /*--- BCheck^(-1)*(h-P_I(h+mu/gamma)) ---*/

    if((config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) && (!active)) {
      for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
        helper -= BCheck_Inv[iConstr][jConstr]*Lambda_Old[jConstr]/gamma;
      }
      // helper = -gamma*ConstrFunc_Old[iConstr]-Lambda_Old[iConstr];
      Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    }
    else {
      for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
        helper += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Old[jConstr];
      }
      Lambda[iConstr] += helper*stepsize*config->GetMultiplierScale(iConstr);
    }

    // /*--- gamma*(h-P_I(h+mu/gamma)) ---*/
    // if(active) Lambda[iConstr] = Lambda_Tilde[iConstr];
    // if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
    //   Lambda[iConstr] += stepsize*gamma*ConstrFunc_Old[iConstr];
    // }
    // else {
    //   Lambda[iConstr] = 0.;
    // }
    // Lambda_Tilde[iConstr] += stepsize*gamma*ConstrFunc_Old[iConstr];


    // if(config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) {
    //   Lambda[iConstr] = max(Lambda[iConstr], 0.);
    //   // Lambda_Tilde[iConstr] = max(Lambda_Tilde[iConstr], 0.);
    // }
  }
}

void COneShotFluidDriver::StoreLambdaGrad() {
  if(nConstr > 0) {
    unsigned short nVar = solver[ADJFLOW_SOL]->GetnVar();
    unsigned long nPointDomain = geometry->GetnPointDomain();
    const su2double beta = config->GetOneShotBeta();
    for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++) {
      const su2double gamma = config->GetOneShotGamma(iConstr);
      // const bool active = (ConstrFunc_Old[iConstr] + Lambda_Old[iConstr]/gamma > 0.);
      const bool active = (ConstrFunc_Store[iConstr] > 0.);
      su2double my_Gradient = 0.;
      // if((config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) || (active)) {
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
          for (unsigned short iVar = 0; iVar < nVar; iVar++) {
            my_Gradient += beta
                * solver[ADJFLOW_SOL]->GetConstrDerivative(iConstr, iPoint, iVar)
                * solver[ADJFLOW_SOL]->GetNodes()->GetSolution_DeltaStore(iPoint,iVar);
          }
        }
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_Gradient, &AugLagLamGrad[iConstr], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  AugLagLamGrad[iConstr] = my_Gradient;
#endif
        AugLagLamGrad[iConstr] += ConstrFunc_Old[iConstr] + Lambda_Old[iConstr]/gamma;
      // }
    }
  }
}

void COneShotFluidDriver::InitializeLambdaTilde(unsigned short iConstr) {
  unsigned short nVar = solver[ADJFLOW_SOL]->GetnVar();
  unsigned long nPointDomain = geometry->GetnPointDomain();
  const su2double beta = config->GetOneShotBeta();
  const su2double gamma = config->GetOneShotGamma(iConstr);
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
  Lambda_Init -= ConstrFunc_Old[iConstr];
  Lambda_Tilde[iConstr] = gamma*Lambda_Init;
  // if(config->GetKind_ConstrFuncType(iConstr) != EQ_CONSTR) Lambda_Tilde[iConstr] = max(Lambda_Tilde[iConstr], 0.0);
}

void COneShotFluidDriver::StoreObjFunction(){
  ObjFunc_Store = ObjFunc;
}

void COneShotFluidDriver::StoreConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFunc_Store[iConstr] = ConstrFunc[iConstr];
  }
}

void COneShotFluidDriver::StoreOldConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFunc_Old[iConstr] = ConstrFunc[iConstr];
  }
}
