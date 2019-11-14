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

  ShiftedLagrangianGradient = new su2double[nDV_Total];
  ShiftedLagrangianGradient_Old = new su2double[nDV_Total];

  AugmentedLagrangianGradient = new su2double[nDV_Total];
  AugmentedLagrangianGradient_Old = new su2double[nDV_Total];

  DesignVarUpdate = new su2double[nDV_Total];
  DesignVariable = new su2double[nDV_Total];
  SearchDirection = new su2double[nDV_Total];
  activeset = new bool[nDV_Total];

  BFGS_Inv = new su2double*[nDV_Total];

  if(nConstr > 0){
    ConstrFunc = new su2double[nConstr];
    ConstrFunc_Store = new su2double[nConstr];
    Multiplier = new su2double[nConstr];
    Multiplier_Old = new su2double[nConstr];
    AugmentedLagrangianMultiplierGradient = new su2double[nConstr];
    BCheck_Inv = new su2double*[nConstr];
  }

  BFGS_Init = config->GetBFGSInitValue();

  for (iDV = 0; iDV  < nDV_Total; iDV++){
    Gradient[iDV] = 0.0;
    Gradient_Old[iDV] = 0.0;
    ShiftedLagrangianGradient[iDV] = 0.0;
    ShiftedLagrangianGradient_Old[iDV] = 0.0;
    AugmentedLagrangianGradient[iDV] = 0.0;
    AugmentedLagrangianGradient_Old[iDV] = 0.0;
    DesignVarUpdate[iDV] = 0.0;
    DesignVariable[iDV] = 0.0;
    SearchDirection[iDV] = 0.0;
    activeset[iDV]=false;
    BFGS_Inv[iDV] = new su2double[nDV_Total];
    for (jDV = 0; jDV < nDV_Total; jDV++){
      BFGS_Inv[iDV][jDV] = 0.0;
      if (iDV==jDV) BFGS_Inv[iDV][jDV] = BFGS_Init;
    }
  }

  for (unsigned short iConstr = 0; iConstr  < nConstr; iConstr++){
    ConstrFunc[iConstr] = 0.0;
    ConstrFunc_Store[iConstr] = 0.0;
    Multiplier[iConstr] = 0.0;
    Multiplier_Old[iConstr] = 0.0;
    AugmentedLagrangianMultiplierGradient[iConstr] = 0.0;
    BCheck_Inv[iConstr] = new su2double[nConstr];
    for (unsigned short jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck_Inv[iConstr][jConstr] = 0.0;
    }
    BCheck_Inv[iConstr][iConstr] = config->GetOneShotGamma();
  }
  BCheck_Norm = pow(config->GetBCheckEpsilon(), su2double(nConstr));

  /*----- calculate values for bound projection algorithm -------*/
  lb=-config->GetBound()*config->GetDesignScale();
  ub=config->GetBound()*config->GetDesignScale();
  epsilon=(ub-lb)/2.0;

  /*---- calculate line search parameter ----*/
  cwolfeone= 1E-4*config->GetDesignScale();

  grid_movement[ZONE_0][INST_0] = new CVolumetricMovement(geometry, config);
  surface_movement[ZONE_0]      = new CSurfaceMovement();

}

COneShotFluidDriver::~COneShotFluidDriver(void){

  /*----- free allocated memory -------*/
  unsigned short iDV;
  for (iDV = 0; iDV  < nDV_Total; iDV++){
    delete [] BFGS_Inv[iDV];
  }
  delete [] BFGS_Inv;
  delete [] Gradient;
  delete [] Gradient_Old;
  delete [] ShiftedLagrangianGradient;
  delete [] ShiftedLagrangianGradient_Old;

  delete [] AugmentedLagrangianGradient;
  delete [] AugmentedLagrangianGradient_Old;

  delete [] DesignVarUpdate;
  delete [] DesignVariable;
  delete [] SearchDirection;
  delete [] activeset;

  if(nConstr > 0){
    delete [] ConstrFunc;
    delete [] Multiplier;
    delete [] Multiplier_Old;
    delete [] ConstrFunc_Store;
    delete [] AugmentedLagrangianMultiplierGradient;
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
    const bool StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                             solver_container, numerics_container, config_container,
                                             surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    if(StopCalc) break;

  }

}

void COneShotFluidDriver::RunOneShot(){

  su2double stepsize = 1.0, stepsize_tmp, tol = config->GetOneShotSearchTol();
  unsigned short ArmijoIter = 0, nArmijoIter = config->GetOneShotSearchIter();
  unsigned long InnerIter = config->GetInnerIter();
  bool bool_tol = false;

  /*--- Store the old solution and the old design for line search ---*/
  solver[ADJFLOW_SOL]->SetStoreSolution();
  solver[ADJFLOW_SOL]->SetMeshPointsOld(config, geometry);

  /*--- This is the line search loop that is only called once, if no update is performed ---*/
  do {

    if(InnerIter > config->GetOneShotStart() && InnerIter < config->GetOneShotStop()){      

      // if(ArmijoIter > 1) {
      //   /*--- Cubic backtracking ---*/
      //   stepsize_tmp = UpdateStepSizeCubic(stepsize, stepsize_p);
      //   Lagrangian_p = Lagrangian;
      //   stepsize_p   = stepsize;
      //   stepsize     = UpdateStepSizeBound(stepsize_tmp, stepsize/10., stepsize/2.);
      //   if(stepsize < tol) {
      //     stepsize = tol;
      //     bool_tol = true;
      //   }

      //   /*---Load the old design for line search---*/
      //   solver[ADJFLOW_SOL]->LoadMeshPointsOld(config, geometry);
      //   LoadMultiplier();
      //   UpdateMultiplier(stepsize);
      // }
      // else if(ArmijoIter > 0){
      if(ArmijoIter > 0){
        /*--- Parabolic backtracking ---*/
        stepsize_tmp = UpdateStepSizeQuadratic();
        // Lagrangian_p = Lagrangian;
        // stepsize_p   = stepsize;
        stepsize     = UpdateStepSizeBound(stepsize_tmp, stepsize/10., stepsize/2.);
        if(stepsize < tol) {
          stepsize = tol;
          bool_tol = true;
        }

        /*---Load the old design for line search---*/
        solver[ADJFLOW_SOL]->LoadMeshPointsOld(config, geometry);
        LoadMultiplier();
        // UpdateMultiplier(stepsize);
      }
      else{
        /*--- Store and update constraint multiplier ---*/
        StoreMultiplier();
        StoreMultiplierGrad();
        // UpdateMultiplier(1.0);
      }

      /*--- Compute and store GradL dot p ---*/
      StoreGradDotDir();

      /*--- Update multiplier ---*/
      UpdateMultiplier(stepsize);

      /*--- Load the old solution for line search (either y_k or y_k-1) ---*/
      solver[ADJFLOW_SOL]->LoadSolution();

      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      if ((ArmijoIter != nArmijoIter && !bool_tol) || (!config->GetZeroStep())) {
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
        LoadMultiplier();
        ComputeDesignVarUpdate(0.0);
      }

    }

    /*--- Do a primal and adjoint update ---*/
    PrimalDualStep();
    solver[ADJFLOW_SOL]->SetSolutionDelta(geometry);

    /*--- Update multiplier ---*/
    // if(ArmijoIter == 0) UpdateMultiplier(1.0);

    /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
    CalculateLagrangian();

    ArmijoIter++;

  } while(InnerIter > config->GetOneShotStart() && 
          InnerIter < config->GetOneShotStop()  &&
          !CheckFirstWolfe()                    && 
          ArmijoIter < nArmijoIter+1            &&
          !bool_tol);

  /*--- Store number of search iterations ---*/
  solver[ADJFLOW_SOL]->SetArmijoIter(ArmijoIter);

  /*--- Store FFD info in file ---*/
  if (((config->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) ||
       (config->GetDesign_Variable(0) == FFD_CONTROL_POINT))   &&
       InnerIter > config->GetOneShotStart()                   && 
       InnerIter < config->GetOneShotStop()                    &&
       (!config->GetZeroStep() || CheckFirstWolfe())) {
    surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container[ZONE_0][INST_0], config_container, false);
    config->SetMesh_FileName(config->GetMesh_Out_FileName());
  }

  if(InnerIter >= config->GetOneShotStart() && 
     InnerIter < config->GetOneShotStop()   && 
     InnerIter > 1) {
    /*--- Calculate Lagrangian with new Alpha, Beta, and Gamma ---*/
    solver[ADJFLOW_SOL]->CalculateRhoTheta(config);
    solver[ADJFLOW_SOL]->CalculateAlphaBetaGamma(config, BCheck_Norm);
    /*--- Store the constraint function, and set the multiplier to 0 if the sign is opposite ---*/
    StoreConstrFunction();
    // CheckMultiplier();
    CalculateLagrangian();
  }

  /*--- Store Deltay and DeltaBary ---*/
  solver[ADJFLOW_SOL]->SetStoreSolutionDelta();

  if(InnerIter >= config->GetOneShotStart() && InnerIter < config->GetOneShotStop()){

    /*--- Update design variable ---*/
    UpdateDesignVariable();

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
      solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, config->GetOneShotGamma());
    }
    solver[ADJFLOW_SOL]->LoadSolution();

    /*--- Alpha*Deltay^T*G_u ---*/
    ComputeAlphaTerm();
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, config->GetOneShotAlpha());
    solver[ADJFLOW_SOL]->LoadSolution();

    /*--- Beta*DeltaBary^T*N_yu ---*/
    solver[ADJFLOW_SOL]->UpdateStateVariable(config);
    ComputeBetaTerm();
    solver[ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry, config);
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry, config->GetOneShotBeta());
    solver[ADJFLOW_SOL]->LoadSaveSolution();

    /*--- Projection of the gradient N_u---*/
    solver[ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry);
    ProjectMeshSensitivities();
    SetShiftedLagrangianGradient();

    /*--- Projection of the gradient L_u---*/
    solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry); //Lagrangian
    ProjectMeshSensitivities();
    SetAugmentedLagrangianGradient();

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
  SetAdj_ConstrFunction(Multiplier);

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
    yk[iDV]=ProjectionSet(iDV, AugmentedLagrangianGradient[iDV]-AugmentedLagrangianGradient_Old[iDV], false);
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
    // solver[ADJFLOW_SOL]->CalculateAlphaBetaGamma(config, BCheck_Norm);
    // CalculateLagrangian();
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

bool COneShotFluidDriver::CheckFirstWolfe(){
  unsigned short iDV;
  su2double admissible_step = 0.0;

  for (iDV = 0; iDV < nDV_Total; iDV++){
    /*--- AugmentedLagrangianGradient is the gradient at the old iterate. ---*/
    admissible_step += DesignVarUpdate[iDV]*AugmentedLagrangianGradient[iDV];
    // admissible_step += DesignVarUpdate[iDV]*ShiftedLagrangianGradient[iDV];
  }
  if (nConstr > 0) {
    unsigned short iConstr;
    for (iConstr = 0; iConstr < nConstr; iConstr++) {
      admissible_step += (Multiplier[iConstr]-Multiplier_Old[iConstr])*AugmentedLagrangianMultiplierGradient[iConstr];
    }
  }
  admissible_step *= cwolfeone;

  return (Lagrangian <= Lagrangian_Old + admissible_step);
}

void COneShotFluidDriver::StoreGradDotDir(){
  unsigned short iDV;
  GradDotDir = 0.0;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    /*--- AugmentedLagrangianGradient is the gradient at the old iterate. ---*/
    GradDotDir += DesignVarUpdate[iDV]*AugmentedLagrangianGradient[iDV];
    // GradDotDir += DesignVarUpdate[iDV]*ShiftedLagrangianGradient[iDV];
  }
  if (nConstr > 0) {
    unsigned short iConstr;
    for (iConstr = 0; iConstr < nConstr; iConstr++) {
      GradDotDir += (Multiplier[iConstr]-Multiplier_Old[iConstr])*AugmentedLagrangianMultiplierGradient[iConstr];
    }
  }
}

su2double COneShotFluidDriver::UpdateStepSizeQuadratic(){
  return -GradDotDir/(2.*(Lagrangian - Lagrangian_Old - GradDotDir));
}

su2double COneShotFluidDriver::UpdateStepSizeCubic(su2double stepsize, su2double stepsize_p){
  const su2double tmp = 1./(pow(stepsize,2.)*pow(stepsize_p,2.)*(stepsize-stepsize_p));

  const su2double vec1 = Lagrangian   - Lagrangian_Old - stepsize*GradDotDir,
                  vec2 = Lagrangian_p - Lagrangian_Old - stepsize_p*GradDotDir;

  const su2double a = tmp*(pow(stepsize_p,2.)*vec1  - pow(stepsize,2.)*vec2),
                  b = tmp*(-pow(stepsize_p,3.)*vec1 + pow(stepsize,3.)*vec2);

  if(fabs(a) < 1.0E-16) {
    return -GradDotDir/(2.*b); // Function is quadratic
  }
  else {
    const su2double d = pow(b,2.) - 3.*a*GradDotDir; // Discriminant
    return (-b + sqrt(d))/(3.*a); // Function is cubic
  }
}

su2double COneShotFluidDriver::UpdateStepSizeBound(su2double stepsize, su2double a, su2double b){
  if(stepsize < a) {return a;}
  else if(stepsize > b) {return b;}
  else {return stepsize;}
}


void COneShotFluidDriver::ComputeDesignVarUpdate(su2double stepsize){
  unsigned short iDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    DesignVarUpdate[iDV]=BoundProjection(DesignVariable[iDV]+stepsize*SearchDirection[iDV]*config->GetDesignScale())-DesignVariable[iDV];
  }
}

void COneShotFluidDriver::ComputeSearchDirection(){
  unsigned short iDV, jDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    SearchDirection[iDV]=0.0;
    for (jDV=0;jDV<nDV_Total;jDV++){
      SearchDirection[iDV]+= BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,-ShiftedLagrangianGradient[jDV],false);
    }
    SearchDirection[iDV]=-ProjectionSet(iDV, ShiftedLagrangianGradient[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false);
  }
}

void COneShotFluidDriver::ComputeNegativeSearchDirection(){
  unsigned short iDV, jDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    SearchDirection[iDV]=0.0;
    for (jDV=0;jDV<nDV_Total;jDV++){
      SearchDirection[iDV]+=BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,ShiftedLagrangianGradient[jDV],false);
    }
    SearchDirection[iDV]=ProjectionSet(iDV, ShiftedLagrangianGradient[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false);
  }
}

bool COneShotFluidDriver::CheckDescent(){
  unsigned short iDV;
  su2double product = 0.0;
  for (iDV=0;iDV<nDV_Total;iDV++){
    product+=SearchDirection[iDV]*ProjectionSet(iDV, AugmentedLagrangianGradient_Old[iDV], false);
  }
  return (product<=0.0);
}

void COneShotFluidDriver::StoreLagrangianInformation(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    AugmentedLagrangianGradient_Old[iDV] = AugmentedLagrangianGradient[iDV];
  }
  Lagrangian_Old = Lagrangian;
}

void COneShotFluidDriver::UpdateDesignVariable(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    DesignVariable[iDV] += DesignVarUpdate[iDV];
  }
}

void COneShotFluidDriver::CalculateLagrangian(){
    
  Lagrangian = 0.0;
  Lagrangian += ObjFunc_Store; //TODO use for BFGS either only objective function or normal Lagrangian

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    // Lagrangian += config->GetOneShotGamma()/2.*ConstrFunc_Store[iConstr]*ConstrFunc_Store[iConstr]
    //             + ConstrFunc_Store[iConstr]*Multiplier[iConstr]
    //             + config->GetBCheckEpsilon()*Multiplier[iConstr]*Multiplier[iConstr];
    su2double helper = ConstrFunc_Store[iConstr] + Multiplier[iConstr]/config->GetOneShotGamma();
    /*--- Lagrangian += gamma/2 ||(h + mu/gamma)_+||^2 ---*/
    if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) {
      Lagrangian += config->GetOneShotGamma()/2.*helper*helper;
    }
    else {
      Lagrangian += config->GetOneShotGamma()/2.*max(helper,0.)*max(helper,0.);
    }
  }

  Lagrangian += solver[ADJFLOW_SOL]->CalculateLagrangian(config);

}

su2double COneShotFluidDriver::ProjectionSet(unsigned short iDV, su2double value, bool active){
  if (active) {
      if(!activeset[iDV]) value = 0.0;
  } else {
      if(activeset[iDV]) value = 0.0;
  }
  return value;
}

su2double COneShotFluidDriver::ProjectionPAP(unsigned short iDV, unsigned short jDV, su2double value, bool active){
  //returns for a Matrix entry a_iDV,jDV of a Matrix A the resulting entry of P*A*P (active or inactive set)
  if (active) {
      if(!activeset[iDV] || !activeset[jDV]) value = 0.0;
  } else {
      if(activeset[iDV] || activeset[jDV]) value = 0.0;
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
    norm += (DesignVariable[iDV]-BoundProjection(DesignVariable[iDV]-stepsize*ShiftedLagrangianGradient[iDV]))
          * (DesignVariable[iDV]-BoundProjection(DesignVariable[iDV]-stepsize*ShiftedLagrangianGradient[iDV]));
  }
  norm = sqrt(norm);
  epsilon = min(norm, (ub-lb)/2.0);
  unsigned short nActive = nDV_Total;

  for (iDV = 0; iDV < nDV_Total; iDV++) {
    activeset[iDV] = false;
    if(ub-DesignVariable[iDV] <= epsilon)      activeset[iDV] = true;
    else if(DesignVariable[iDV]-lb <= epsilon) activeset[iDV] = true;
    else                                       nActive--;
  }
  solver[ADJFLOW_SOL]->SetnActiveDV(nActive);
}

void COneShotFluidDriver::SetShiftedLagrangianGradient(){
  unsigned short iDV;
  su2double norm = 0.;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    ShiftedLagrangianGradient[iDV] = Gradient[iDV];
    norm += Gradient[iDV]*Gradient[iDV];
  }
  solver[ADJFLOW_SOL]->SetShiftedLagGradNorm(sqrt(norm));
}

void COneShotFluidDriver::SetAugmentedLagrangianGradient(){
  unsigned short iDV;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    AugmentedLagrangianGradient[iDV] = Gradient[iDV];
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
      if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR || 
         ConstrFunc[iConstr] + Multiplier[iConstr]/config->GetOneShotGamma() > 0.) {
        seeding[iConstr] = ConstrFunc[iConstr];
      }
      else {
        seeding[iConstr] = 0.;
      }
    }
    SetAdj_ConstrFunction(ConstrFunc);

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
    SetAdj_ConstrFunction(Multiplier);

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
    BCheck[iConstr][iConstr] = 1./config->GetOneShotGamma();
    for (jConstr = 0; jConstr < nConstr; jConstr++){
      BCheck[iConstr][jConstr] += config->GetOneShotBeta()*solver[ADJFLOW_SOL]->MultiplyConstrDerivative(iConstr,jConstr);
    }
  }
  if (nConstr == 1){
      BCheck_Norm = BCheck[0][0] - 1./config->GetOneShotGamma();
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

void COneShotFluidDriver::StoreMultiplier(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Multiplier_Old[iConstr] = Multiplier[iConstr];
  }
}

void COneShotFluidDriver::LoadMultiplier(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Multiplier[iConstr] = Multiplier_Old[iConstr];
  }
}

void COneShotFluidDriver::UpdateMultiplier(su2double stepsize){
  su2double helper;
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    /*--- BCheck^(-1)*h ---*/
    helper = 0.0;
    for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++){
       helper += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
    }
    Multiplier[iConstr] += helper*stepsize*ConstrFunc_Store[iConstr]*config->GetMultiplierScale(iConstr);
    // /*--- gamma*h ---*/
    // Multiplier[iConstr] += config->GetOneShotGamma()*stepsize*ConstrFunc_Store[iConstr]*config->GetMultiplierScale(iConstr);
    if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) {
      if(Multiplier[iConstr]*ConstrFunc_Store[iConstr] < 0.) {
        Multiplier[iConstr] = ConstrFunc_Store[iConstr];
      }
    }
    else {
      Multiplier[iConstr] = max(Multiplier[iConstr], 0.);
    }
  }
}

void COneShotFluidDriver::StoreMultiplierGrad() {
  if(nConstr > 0) {
    unsigned short iConstr, iVar, nVar = solver[ADJFLOW_SOL]->GetnVar();
    unsigned long iPoint, nPointDomain = geometry->GetnPointDomain();
    su2double beta = config->GetOneShotBeta(), gamma = config->GetOneShotGamma();
    for (iConstr = 0; iConstr < nConstr; iConstr++) {
      su2double my_Gradient = 0.;
      if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR || 
         ConstrFunc[iConstr] + Multiplier[iConstr]/gamma > 0.) {
        my_Gradient = ConstrFunc[iConstr] + 1./gamma*Multiplier[iConstr];
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          for (iVar = 0; iVar < nVar; iVar++) {
            my_Gradient += beta
                * solver[ADJFLOW_SOL]->GetConstrDerivative(iConstr, iPoint, iVar)
                * solver[ADJFLOW_SOL]->GetNodes()->GetSolution_Delta(iPoint,iVar);
          }
        }
      }
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_Gradient, &AugmentedLagrangianMultiplierGradient[iConstr], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  AugmentedLagrangianMultiplierGradient[iConstr] = my_Gradient;
#endif
    }
  }
}

void COneShotFluidDriver::CheckMultiplier(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    if(config->GetKind_ConstrFuncType(iConstr) == EQ_CONSTR) {
      if(Multiplier[iConstr]*ConstrFunc_Store[iConstr] < 0.) {
        Multiplier[iConstr] = 0.;
        for(unsigned short jConstr = 0; jConstr < nConstr; jConstr++) {
          Multiplier[iConstr] += BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
        }
      }
    }
    else {
      Multiplier[iConstr] = max(Multiplier[iConstr], 0.);
    }
  }
}

void COneShotFluidDriver::StoreConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFunc_Store[iConstr] = ConstrFunc[iConstr];
  }
}
