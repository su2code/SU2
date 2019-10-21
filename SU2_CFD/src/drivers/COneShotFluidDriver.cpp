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

#include "../../include/drivers/COneShotFluidOutput.hpp"

COneShotFluidDriver::COneShotFluidDriver(char* confFile,
                                         unsigned short val_nZone,
                                         SU2_Comm MPICommunicator) : CDiscAdjSinglezoneDriver(confFile, val_nZone, MPICommunicator) {
  unsigned short iDV, jDV, iDV_Value;

  /*---------- One-shot works on all design variables - find total number of design variables ---------*/
  nDV_Total = 0;
  for (iDV = 0; iDV  < config_container[ZONE_0]->GetnDV(); iDV++){
    for (iDV_Value = 0; iDV_Value < config_container[ZONE_0]->GetnDV_Value(iDV); iDV_Value++){
      nDV_Total++;
    }
  }

  nConstr = config_container[ZONE_0]->GetnConstr();

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
    Multiplier = new su2double[nConstr];
    Multiplier_Old = new su2double[nConstr];
    ConstrFunc_Store = new su2double[nConstr];
    BCheck_Inv = new su2double*[nConstr];
  }

  nBFGSmax = config_container[ZONE_0]->GetLimitedMemoryIter();
  nBFGS = 0;
  ykvec = new su2double*[nBFGSmax];
  skvec = new su2double*[nBFGSmax];

  BFGS_Init = config_container[ZONE_0]->GetBFGSInitValue();

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
    BCheck_Inv[iConstr] = new su2double[nConstr];
    for (unsigned short jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck_Inv[iConstr][jConstr] = 0.0;
    }
    BCheck_Inv[iConstr][iConstr] = config_container[ZONE_0]->GetBCheckEpsilon();
  }
  BCheck_Norm = sqrt(su2double(nConstr))*config_container[ZONE_0]->GetBCheckEpsilon();

  for (unsigned short iBFGS = 0; iBFGS < nBFGSmax; iBFGS++){
    ykvec[iBFGS] = new su2double[nDV_Total];
    skvec[iBFGS] = new su2double[nDV_Total];
    for (jDV = 0; jDV < nDV_Total; jDV++){
      ykvec[iBFGS][jDV] = 0.0;
      skvec[iBFGS][jDV] = 0.0;
    }
  }

  /*----- calculate values for bound projection algorithm -------*/
  lb=-config_container[ZONE_0]->GetBound()*config_container[ZONE_0]->GetDesignScale();
  ub=config_container[ZONE_0]->GetBound()*config_container[ZONE_0]->GetDesignScale();
  epsilon=(ub-lb)/2.0;

  /*---- calculate line search parameter ----*/
  cwolfeone= 1E-4*config_container[ZONE_0]->GetDesignScale();

  for (unsigned short iZone = 0; iZone < nZone; iZone++){
    grid_movement[iZone][INST_0] = new CVolumetricMovement(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
    surface_movement[iZone] = new CSurfaceMovement();
  }

  update = false;

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
  }

  for (unsigned short iBFGS = 0; iBFGS < nBFGSmax; iBFGS++){
    delete [] ykvec[iBFGS];
    delete [] skvec[iBFGS];
  }
  delete [] ykvec;
  delete [] skvec;

}

void COneShotFluidDriver::Preprocess(unsigned long TimeIter) {

  unsigned short iZone;

  for(iZone = 0; iZone < nZone; iZone++) {

    /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
    /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/

    config_container[iZone]->SetExtIter(TimeIter);

    /*--- NOTE: Inv Design Routines moved to CDiscAdjFluidIteration::Preprocess ---*/

    /*--- Preprocess the adjoint iteration ---*/

    iteration_container[iZone][INST_0]->Preprocess(output, integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox, iZone, INST_0);

  }

}


void COneShotFluidDriver::Run(){

  config_container[ZONE_0]->SetIntIter(TimeIter);
 
  /*--- Run an iteration of the one-shot solver ---*/
  RunOneShot();

  /*--- Screen output ---*/
  bool StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                     solver_container, numerics_container, config_container,
                                     surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

}

void COneShotFluidDriver::RunOneShot(){

  su2double stepsize = config_container[ZONE_0]->GetStepSize();
  unsigned short maxcounter = config_container[ZONE_0]->GetOneShotMaxCounter();
  unsigned short whilecounter = 0;

  /*--- Store the old solution and the old design for line search ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->StoreSolution();
    solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->StoreMeshPoints(config_container[iZone], geometry_container[iZone][INST_0][MESH_0]);
  }

  /*--- This is the line search loop that is only called once, if no update is performed ---*/
  do {

    if(TimeIter>config_container[ZONE_0]->GetOneShotStart() && TimeIter<config_container[ZONE_0]->GetOneShotStop()){

      if(whilecounter > 0){
        /*--- Armijo line search (halve step) ---*/
        stepsize=stepsize*0.5;

        /*---Load the old design for line search---*/
        for (iZone = 0; iZone < nZone; iZone++){
          solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->LoadMeshPoints(config_container[iZone], geometry_container[iZone][INST_0][MESH_0]);
        }
        LoadMultiplier();
        UpdateMultiplier(stepsize);
      }
      else{
        /*--- Update constraint multiplier ---*/
        StoreMultiplier();
        UpdateMultiplier(1.0);
      }

      /*--- Load the old solution for line search (either y_k or y_k-1) ---*/
      for (iZone = 0; iZone < nZone; iZone++){
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->LoadSolution();
      }

      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      if (whilecounter != maxcounter || (!config_container[ZONE_0]->GetZeroStep())) {
        ComputeDesignVarUpdate(stepsize);
      }
      else {
        LoadMultiplier();
        ComputeDesignVarUpdate(0.0);
      }
      
      for (iZone = 0; iZone < nZone; iZone++){
        config_container[iZone]->SetKind_SU2(SU2_DEF); // set SU2_DEF as the solver
        SurfaceDeformation(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], surface_movement[iZone], grid_movement[iZone][INST_0]);
        config_container[iZone]->SetKind_SU2(SU2_CFD); // set SU2_CFD as the solver
      }

    }

    /*--- Do a primal and adjoint update ---*/
    PrimalDualStep();

    /*--- Estimate Alpha, Beta, and Gamma ---*/
    if(TimeIter >= config_container[ZONE_0]->GetOneShotStart() && TimeIter > 1) solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->CalculateAlphaBetaGamma(config_container[ZONE_0]);

    /*--- Calculate Lagrangian with old Alpha, Beta, and Gamma ---*/
    CalculateLagrangian(true);

    whilecounter++;

  } while(TimeIter > config_container[ZONE_0]->GetOneShotStart() && 
          TimeIter < config_container[ZONE_0]->GetOneShotStop() &&
          (!CheckFirstWolfe()) && whilecounter<maxcounter+1);

  /*--- Store FFD info in file ---*/
  if (((config_container[ZONE_0]->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) ||
       (config_container[ZONE_0]->GetDesign_Variable(0) == FFD_CONTROL_POINT))   &&
       update) {
    surface_movement[ZONE_0]->WriteFFDInfo(surface_movement, geometry_container[ZONE_0][INST_0], config_container, false);
    config_container[ZONE_0]->SetMesh_FileName(config_container[ZONE_0]->GetMesh_Out_FileName());
  }

  if (!CheckFirstWolfe() && config_container[ZONE_0]->GetZeroStep()) stepsize = 0.0;

  /*--- Calculate Lagrangian with new Alpha, Beta, and Gamma ---*/
  if(TimeIter > 1) solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->SetAlphaBetaGamma(config_container[ZONE_0], BCheck_Norm);
  CalculateLagrangian(true);

  if(TimeIter >= config_container[ZONE_0]->GetOneShotStart() && TimeIter < config_container[ZONE_0]->GetOneShotStop()){

    /*--- Update design variable ---*/
    UpdateDesignVariable();

    /*--- Store the constraint function, and set the multiplier to 0 if the sign is opposite ---*/
    StoreConstrFunction();
    CheckMultiplier();

    /*--- N_u ---*/
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SaveSensitivity(geometry_container[iZone][INST_0][MESH_0]);
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->StoreSaveSolution();
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry_container[iZone][INST_0][MESH_0]);
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry_container[iZone][INST_0][MESH_0],1.0);
    }   

    if((nConstr > 0) && (!config_container[ZONE_0]->GetConstPrecond())) ComputePreconditioner();

    /*--- Gamma*h^T*h_u ---*/
    if(nConstr > 0) {
      ComputeGammaTerm();
      for (iZone = 0; iZone < nZone; iZone++){
        solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]->GetOneShotGamma());
      }
    }

    /*--- Alpha*Deltay^T*G_u ---*/
    ComputeAlphaTerm();
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]->GetOneShotAlpha());
    }

    /*--- Beta*DeltaBary^T*N_yu ---*/
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->LoadSolution();
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->UpdateStateVariable(config_container[iZone]);
    }
    ComputeBetaTerm();
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]->GetOneShotBeta());
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
    }

    /*--- Projection of the gradient N_u---*/
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry_container[iZone][INST_0][MESH_0]);
    }
    ProjectMeshSensitivities();
    SetShiftedLagrangianGradient();

    /*--- Projection of the gradient L_u---*/
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry_container[iZone][INST_0][MESH_0]); //Lagrangian
    }
    ProjectMeshSensitivities();
    SetAugmentedLagrangianGradient();

    /*--- Use N_u to compute the active set (bound constraints) ---*/
    ComputeActiveSet(stepsize);

    /*--- Do a BFGS update to approximate the inverse preconditioner ---*/
    if(TimeIter > config_container[ZONE_0]->GetOneShotStart()) BFGSUpdate(config_container[ZONE_0]);

    /*--- Compute the search direction for the line search procedure ---*/
    ComputeSearchDirection();

    StoreLagrangianInformation();
  }
}

void COneShotFluidDriver::PrimalDualStep(){

  unsigned short iZone = 0;
  unsigned short iInst = 0;

  /*--- Note: Unsteady cases not applicable to the one-shot method yet! ---*/

  SetRecording(COMBINED);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

    config_container[iZone]->SetIntIter(0);

    iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone, iInst);

  }

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  SetAdj_ConstrFunction(Multiplier);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Iterate(output, integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, iZone, iInst);
  }

  /*--- Extract the computed sensitivity values. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();
}

void COneShotFluidDriver::SetRecording(unsigned short kind_recording){
  unsigned short iZone;
  unsigned short iInst = 0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->SetRecording(solver_container, geometry_container, config_container, iZone, INST_0, kind_recording);
  }

  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if (kind_recording != NONE){

    AD::StartRecording();

    if (rank == MASTER_NODE && kind_recording == MainVariables && (TimeIter == 0)) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store the primal computational graph." << endl;
      cout << "Combined recording of flow and design variables." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }

    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->RegisterInput(solver_container, geometry_container, config_container, iZone, iInst, kind_recording);
    }

  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->SetDependencies(solver_container, geometry_container, numerics_container, config_container, iZone, iInst, kind_recording);
  }

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun(kind_recording);

  RecordingState = kind_recording;

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->RegisterOutput(solver_container, geometry_container, config_container, output, iZone, iInst);
  }

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
  AD::Reset();

}

void COneShotFluidDriver::SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double* Gradient){
  unsigned short iDV, nDV, iFFDBox, nDV_Value, iMarker, iDim;
  unsigned long iVertex, iPoint;
  su2double delta_eps, my_Gradient, localGradient, *Normal, dS, *VarCoord, Sensitivity,
  dalpha[3], deps[3], dalpha_deps;
  bool *UpdatePoint, MoveSurface, Local_MoveSurface;
  CFreeFormDefBox **FFDBox;
    unsigned long nDV_Count = 0;

  int rank = SU2_MPI::GetRank();

  nDV = config->GetnDV();

  /*--- Boolean controlling points to be updated ---*/

  UpdatePoint = new bool[geometry->GetnPoint()];

  /*--- Definition of the FFD deformation class ---*/

  unsigned short nFFDBox = MAX_NUMBER_FFD;
  FFDBox = new CFreeFormDefBox*[nFFDBox];
  for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) FFDBox[iFFDBox] = NULL;

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (nDV_Value != 1){
      SU2_MPI::Error("The projection using finite differences currently only supports a fixed direction of movement for FFD points.", CURRENT_FUNCTION);
    }
  }

  /*--- Continuous adjoint gradient computation ---*/

  for (iDV = 0; iDV < nDV; iDV++) {
    config_container[ZONE_0]->SetDV_Value(iDV,0, 1E-4);
    MoveSurface = true;
    Local_MoveSurface = true;

    /*--- Free Form deformation based ---*/

    if ((config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
        (config->GetDesign_Variable(iDV) == FFD_NACELLE) ||
        (config->GetDesign_Variable(iDV) == FFD_GULL) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST) ||
        (config->GetDesign_Variable(iDV) == FFD_ROTATION) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS) ||
        (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK)) {

      /*--- Read the FFD information in the first iteration ---*/

      if (iDV == 0) {


        /*--- Read the FFD information from the grid file ---*/

        surface_movement->ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());

        /*--- If the FFDBox was not defined in the input file ---*/
        if (!surface_movement->GetFFDBoxDefinition()) {
          SU2_MPI::Error("The input grid doesn't have the entire FFD information!", CURRENT_FUNCTION);
        }

        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

          surface_movement->CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);

          surface_movement->CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);

        }

      }

      /*--- Apply the control point change ---*/

      MoveSurface = false;

      for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

        /*--- Reset FFD box ---*/
        switch (config->GetDesign_Variable(iDV) ) {
          case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER_2D :        Local_MoveSurface = surface_movement->SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS_2D :     Local_MoveSurface = surface_movement->SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST_2D :         Local_MoveSurface = surface_movement->SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_NACELLE :          Local_MoveSurface = surface_movement->SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_GULL :             Local_MoveSurface = surface_movement->SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST :            Local_MoveSurface = surface_movement->SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ROTATION :         Local_MoveSurface = surface_movement->SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER :           Local_MoveSurface = surface_movement->SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS :        Local_MoveSurface = surface_movement->SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_SURFACE :  Local_MoveSurface = surface_movement->SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
        }

        /*--- Recompute cartesian coordinates using the new control points position ---*/

        if (Local_MoveSurface) {
          MoveSurface = true;
          surface_movement->SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, true);
        }

      }

    }

    /*--- Hicks Henne design variable ---*/

    else if (config->GetDesign_Variable(iDV) == HICKS_HENNE) {
      surface_movement->SetHicksHenne(geometry, config, iDV, true);
    }

    /*--- Surface bump design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SURFACE_BUMP) {
      surface_movement->SetSurface_Bump(geometry, config, iDV, true);
    }

    /*--- Kulfan (CST) design variable ---*/

    else if (config->GetDesign_Variable(iDV) == CST) {
      surface_movement->SetCST(geometry, config, iDV, true);
    }

    /*--- Displacement design variable ---*/

    else if (config->GetDesign_Variable(iDV) == TRANSLATION) {
      surface_movement->SetTranslation(geometry, config, iDV, true);
    }

    /*--- Angle of Attack design variable ---*/

 /*   else if (config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) {
      Gradient[iDV][0] = config->GetAoA_Sens();
    }*/

    /*--- Scale design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SCALE) {
      surface_movement->SetScale(geometry, config, iDV, true);
    }

    /*--- Rotation design variable ---*/

    else if (config->GetDesign_Variable(iDV) == ROTATION) {
      surface_movement->SetRotation(geometry, config, iDV, true);
    }

    /*--- NACA_4Digits design variable ---*/

    else if (config->GetDesign_Variable(iDV) == NACA_4DIGITS) {
      surface_movement->SetNACA_4Digits(geometry, config);
    }

    /*--- Parabolic design variable ---*/

    else if (config->GetDesign_Variable(iDV) == PARABOLIC) {
      surface_movement->SetParabolic(geometry, config);
    }

    /*--- Design variable not implement ---*/

    else {
      if (rank == MASTER_NODE)
        cout << "Design Variable not implemented yet." << endl;
    }

    /*--- Load the delta change in the design variable (finite difference step). ---*/

    if ((config->GetDesign_Variable(iDV) != ANGLE_OF_ATTACK) &&
        (config->GetDesign_Variable(iDV) != FFD_ANGLE_OF_ATTACK)) {

      /*--- If the Angle of attack is not involved, reset the value of the gradient ---*/

      my_Gradient = 0.0; Gradient[nDV_Count] = 0.0;

      if (MoveSurface) {

        delta_eps = config->GetDV_Value(iDV);

        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;

        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry->GetnPointDomain()) && UpdatePoint[iPoint]) {

                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                su2double Prod = 0.0;
                VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();

                dS = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  Sensitivity = geometry->GetSensitivity(iPoint,iDim);
                  Prod+=Normal[iDim]*Sensitivity;
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);
                Sensitivity = -Prod/dS;

                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }

                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
        }

      }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    localGradient = my_Gradient;
#endif
    Gradient[nDV_Count] = localGradient;
    nDV_Count++;
    }
  }

  /*--- Delete memory for parameterization. ---*/

  if (FFDBox != NULL) {
    for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) {
      if (FFDBox[iFFDBox] != NULL) {
        delete FFDBox[iFFDBox];
      }
    }
    delete [] FFDBox;
  }

  delete [] UpdatePoint;
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

    grid_movement->SetVolume_Deformation(geometry, config, false, false);

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

void COneShotFluidDriver::LBFGSUpdate(CConfig *config){
  unsigned long iDV;
  su2double vk=0.0;

  if (nBFGS == nBFGSmax) {
    nBFGS--;
    for (unsigned short iBFGS=0; iBFGS<nBFGSmax-1; iBFGS++){
      for (iDV = 0; iDV < nDV_Total; iDV++){
        ykvec[iBFGS][iDV]= ykvec[iBFGS+1][iDV];
        skvec[iBFGS][iDV]= skvec[iBFGS+1][iDV];
      }
    }
  }
  for (iDV = 0; iDV < nDV_Total; iDV++){
    ykvec[nBFGS][iDV]=ProjectionSet(iDV, AugmentedLagrangianGradient[iDV]-AugmentedLagrangianGradient_Old[iDV], false);
    skvec[nBFGS][iDV]=ProjectionSet(iDV, DesignVarUpdate[iDV], false);
    SearchDirection[iDV]=-ShiftedLagrangianGradient[iDV];
    vk+=ykvec[nBFGS][iDV]*skvec[nBFGS][iDV];
  }
  if(vk>0){
    nBFGS = nBFGS + 1;
    if(config->GetBFGSInit()){
      su2double helper = 0.0;
      for (iDV=0;iDV<nDV_Total;iDV++){
        helper+=ykvec[nBFGS-1][iDV]*ykvec[nBFGS-1][iDV];
      }
      for (iDV=0;iDV<nDV_Total;iDV++){
        for (unsigned short jDV=0;jDV<nDV_Total;jDV++){
          BFGS_Inv[iDV][jDV]= (ykvec[nBFGS-1][iDV]*skvec[nBFGS-1][jDV])/helper;
        }
      }
    }
  }else{
    nBFGS = 0;
  }

  LBFGSUpdateRecursive(config, nBFGS);
}

void COneShotFluidDriver::LBFGSUpdateRecursive(CConfig *config, unsigned short nCounter){
  unsigned long iDV, jDV;
  su2double *helper = new su2double [nDV_Total];
  su2double alpha = 0.0;
  su2double alphahelpone = 0.0;
  su2double alphahelptwo = 0.0;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    SearchDirection[iDV]=ProjectionSet(iDV,SearchDirection[iDV],false);
  }
  if(nCounter == 0){
    for (iDV=0;iDV<nDV_Total;iDV++){
      helper[iDV]=0.0;
      for (jDV=0;jDV<nDV_Total;jDV++){
        helper[iDV]+=BFGS_Inv[iDV][jDV]*SearchDirection[jDV];
      }
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      SearchDirection[iDV] = helper[iDV];
    }
  }
  else{
    for (iDV=0;iDV<nDV_Total;iDV++){
      ykvec[nCounter-1][iDV] = ProjectionSet(iDV, ykvec[nCounter-1][iDV], false);
      skvec[nCounter-1][iDV] = ProjectionSet(iDV, skvec[nCounter-1][iDV], false);
      alphahelpone+=skvec[nCounter-1][iDV]*SearchDirection[iDV];
      alphahelptwo+=ykvec[nCounter-1][iDV]*skvec[nCounter-1][iDV];
    }
    alpha = alphahelpone/alphahelptwo;
    for (iDV=0;iDV<nDV_Total;iDV++){
      SearchDirection[iDV] -= alpha*ykvec[nCounter-1][iDV];
    }
    LBFGSUpdateRecursive(config, nCounter-1);
    alphahelpone = 0.0;
    for (iDV=0;iDV<nDV_Total;iDV++){
      alphahelpone+=ykvec[nCounter-1][iDV]*SearchDirection[iDV];
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      SearchDirection[iDV] += (alpha - alphahelpone/alphahelptwo)*skvec[nCounter-1][iDV];
      SearchDirection[iDV] = ProjectionSet(iDV,SearchDirection[iDV],false);
    }
  }
  delete [] helper;
}

bool COneShotFluidDriver::CheckFirstWolfe(){
  unsigned short iDV;
  su2double admissible_step = 0.0;

  for (iDV=0;iDV<nDV_Total;iDV++){
    /*--- ShiftedLagrangianGradient is the gradient at the old iterate (for One_Shot it is N_u and not L_u) ---*/
    admissible_step += DesignVarUpdate[iDV]*ShiftedLagrangianGradient[iDV];
  }
  admissible_step *= cwolfeone;

  return (Lagrangian<=Lagrangian_Old+admissible_step);
}

void COneShotFluidDriver::ComputeDesignVarUpdate(su2double stepsize){
  unsigned short iDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    DesignVarUpdate[iDV]=BoundProjection(DesignVariable[iDV]+stepsize*SearchDirection[iDV]*config_container[ZONE_0]->GetDesignScale())-DesignVariable[iDV];
  }
}

void COneShotFluidDriver::ComputeSearchDirection(){
  unsigned short iDV, jDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    if(!config_container[ZONE_0]->GetLimitedMemory()){
      SearchDirection[iDV]=0.0;
      for (jDV=0;jDV<nDV_Total;jDV++){
        SearchDirection[iDV]+= BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,-ShiftedLagrangianGradient[jDV],false);
      }
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
  update = true;
}

void COneShotFluidDriver::CalculateLagrangian(bool augmented){
  
  unsigned short iZone;
  
  Lagrangian = 0.0;
  Lagrangian += ObjFunc; //TODO use for BFGS either only objective function or normal Lagrangian

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    Lagrangian += ConstrFunc[iConstr]*Multiplier[iConstr];
  }

  if(augmented){
    for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
      Lagrangian += config_container[ZONE_0]->GetOneShotGamma()/2.*ConstrFunc[iConstr]*ConstrFunc[iConstr];
    }
    for (iZone = 0; iZone < nZone; iZone++) {
      Lagrangian += solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->CalculateLagrangianPart(config_container[iZone], augmented);
    }
  }

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
      if(!activeset[iDV]||!activeset[jDV]) value = 0.0;
  } else {
      if(activeset[iDV]||activeset[jDV]) value = 0.0;
  }
  return value;
}

su2double COneShotFluidDriver::BoundProjection(su2double value){
  if(value<=lb) value = lb;
  if(value>=ub) value = ub;
  return value;
}

void COneShotFluidDriver::ComputeActiveSet(su2double stepsize){
  //Compute ||x-P(x-gradx)||
  unsigned short iDV;
  su2double norm = 0.0;
  for (iDV=0; iDV<nDV_Total; iDV++) {
    norm+=(DesignVariable[iDV]-BoundProjection(DesignVariable[iDV]-stepsize*ShiftedLagrangianGradient[iDV]))*(DesignVariable[iDV]-BoundProjection(DesignVariable[iDV]-stepsize*ShiftedLagrangianGradient[iDV]));
  }
  norm=sqrt(norm);
  if(norm<(ub-lb)/2.0) {
    epsilon=norm;
  }
  for (iDV=0; iDV<nDV_Total; iDV++) {
    activeset[iDV]=false;
    if(ub-DesignVariable[iDV]<=epsilon) activeset[iDV]=true;
    if(DesignVariable[iDV]-lb<=epsilon) activeset[iDV]=true;
  }
}

void COneShotFluidDriver::SetShiftedLagrangianGradient(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    ShiftedLagrangianGradient[iDV] = Gradient[iDV];
  }
}

void COneShotFluidDriver::SetAugmentedLagrangianGradient(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    AugmentedLagrangianGradient[iDV] = Gradient[iDV];
  }
}

void COneShotFluidDriver::ComputeGammaTerm(){
    unsigned short iZone;
    unsigned short iInst = 0;

    /*--- Note: Not applicable for unsteady code ---*/

    /*--- Initialize the adjoint of the output variables of the iteration with difference of the solution and the solution
     *    of the previous iteration. The values are passed to the AD tool. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      config_container[iZone]->SetIntIter(0);
      iteration_container[iZone][INST_0]->InitializeAdjoint_Zero(solver_container, geometry_container, config_container, iZone, iInst);

    }

    /*--- Initialize the adjoint of the objective function with 0.0. ---*/

    SetAdj_ObjFunction_Zero();
    SetAdj_ConstrFunction(ConstrFunc);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->Iterate_No_Residual(output, integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, iZone, iInst);
    }

    /*--- Extract the computed sensitivity values. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();
}

void COneShotFluidDriver::ComputeAlphaTerm(){
    unsigned short iZone;
    unsigned short iInst = 0;

    /*--- Note: Not applicable for unsteady code ---*/

    /*--- Initialize the adjoint of the output variables of the iteration with difference of the solution and the solution
     *    of the previous iteration. The values are passed to the AD tool. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      config_container[iZone]->SetIntIter(0);
      iteration_container[iZone][INST_0]->InitializeAdjoint_Update(solver_container, geometry_container, config_container, iZone, iInst);

    }

    /*--- Initialize the adjoint of the objective function with 0.0. ---*/

    SetAdj_ObjFunction_Zero();
    SetAdj_ConstrFunction_Zero();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->Iterate_No_Residual(output, integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, iZone, iInst);
    }

    /*--- Extract the computed sensitivity values. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    AD::Reset();
}

void COneShotFluidDriver::ComputeBetaTerm(){
    unsigned short iZone = 0;
    unsigned short iInst = 0;

    /*--- Note: Not applicable for unsteady code ---*/

    /*--- For the one shot iteration we have to record for every steady state iteration. ---*/

    /*--- Store the computational graph of one direct iteration with the conservative variables and the mesh coordinates as input. ---*/

    SetRecording(COMBINED);

      /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *    of the previous iteration. The values are passed to the AD tool. ---*/

    for (iZone = 0; iZone < nZone; iZone++) {

      config_container[iZone]->SetIntIter(0);
      iteration_container[iZone][INST_0]->InitializeAdjoint(solver_container, geometry_container, config_container, iZone, iInst);

    }
    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();
    SetAdj_ConstrFunction(Multiplier);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->Iterate_No_Residual(output, integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, iZone, iInst);
    }

    /*--- Extract the computed sensitivity values. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    AD::Reset();
}

void COneShotFluidDriver::ComputePreconditioner(){

  unsigned short iInst = 0;

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

    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetIntIter(0);
      iteration_container[iZone][INST_0]->InitializeAdjoint_Zero(solver_container, geometry_container, config_container, iZone, iInst);
    }

    /*--- Initialize the adjoint of the objective function with 0.0. ---*/

    SetAdj_ObjFunction_Zero();
    SetAdj_ConstrFunction(seeding);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      iteration_container[iZone][INST_0]->Iterate_No_Residual(output, integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, iZone, iInst);
    }

    for (iZone = 0; iZone < nZone; iZone++) {
      solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->SetConstrDerivative(iConstr);
    }


    AD::ClearAdjoints();

    seeding[iConstr]=0.0;

  }

  su2double bcheck=0;
  for (iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr][iConstr] = config_container[ZONE_0]->GetBCheckEpsilon();
    for (jConstr = 0; jConstr  < nConstr; jConstr++){
      for (iZone = 0; iZone < nZone; iZone++) {
        BCheck[iConstr][jConstr] += config_container[ZONE_0]->GetOneShotBeta()*solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->MultiplyConstrDerivative(iConstr,jConstr);
      }
    }
  }
  if (nConstr == 1){
    if(BCheck[0][0]-config_container[ZONE_0]->GetBCheckEpsilon() > 0.) {
      BCheck_Norm = BCheck[0][0] - config_container[ZONE_0]->GetBCheckEpsilon()/config_container[ZONE_0]->GetMultiplierScale(0);
      BCheck_Inv[0][0] = 1./BCheck[0][0];
    }
    else{
      if (rank == MASTER_NODE) cout << "BCheck not positive definite!!!" << endl;
      BCheck_Norm = config_container[ZONE_0]->GetBCheckEpsilon();
      BCheck_Inv[0][0] = 1./config_container[ZONE_0]->GetBCheckEpsilon();
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
  for (iZone = 0; iZone < nZone; iZone++){
    config_container[iZone]->SetKind_SU2(SU2_DOT); // set SU2_DOT as solver
    // get the dependency of the volumetric grid movement
    grid_movement[iZone][INST_0]->SetVolume_Deformation(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], false, true);
  }
  for (iZone = 0; iZone < nZone; iZone++){
    surface_movement[iZone]->CopyBoundary(geometry_container[iZone][INST_0][MESH_0], config_container[iZone]);
    // project sensitivities (surface) on design variables
    if (config_container[iZone]->GetProjectionAD()){
      SetProjection_AD(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], surface_movement[iZone] , Gradient);
    }else{
      SetProjection_FD(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], surface_movement[iZone] , Gradient);
    }

    config_container[iZone]->SetKind_SU2(SU2_CFD); // set SU2_CFD as solver
  }
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

    FunctionValue = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Evaluate_ConstrFunc(config_container[ZONE_0], iConstr);

    ConstrFunc[iConstr] = config_container[ZONE_0]->GetConstraintScale(iConstr)*(FunctionValue - config_container[ZONE_0]->GetConstraintTarget(iConstr));

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
       helper+= BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
    }
    Multiplier[iConstr] = Multiplier[iConstr] + helper*stepsize*config_container[ZONE_0]->GetMultiplierScale(iConstr);
  }
}

void COneShotFluidDriver::CheckMultiplier(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    if(Multiplier[iConstr]*ConstrFunc_Store[iConstr] < 0.) Multiplier[iConstr] = 0.;
  }
}

void COneShotFluidDriver::StoreConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    ConstrFunc_Store[iConstr] = ConstrFunc[iConstr];
  }
}
