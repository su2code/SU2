/*!
 * \file COneShotDriver.cpp
 * \brief The main subroutines for one-shot problems.
 * \author T.Dick
 * \version 7.0.6 "Blackbird"
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

#include "../../include/drivers/COneShotDriver.hpp"

COneShotFluidDriver::COneShotFluidDriver(char* confFile,
                                         unsigned short val_nZone,
                                         SU2_Comm MPICommunicator) : CDiscAdjSinglezoneDriver(confFile, val_nZone, MPICommunicator) {
  unsigned short iDV;

  /*--- Store the pointers ---*/
  /* Since this class inherits from CDiscAdjSinglezoneDriver relevant pointers to
   * config, solver, iteration, geometry, numerics in the ZONE_0 should already be available
   */

  /*---------- One-shot works on all design variables - get the total number of design variables from Config ---------*/
  nDV_Total = config->GetnDV_total();

  nConstr = config[ZONE_0]->GetnConstr();

  Gradient = new su2double[nDV_Total];
  Gradient_Old = new su2double[nDV_Total];

  DesignVarUpdate = new su2double[nDV_Total];
  DesignVariable = new su2double[nDV_Total];
  SearchDirection = new su2double[nDV_Total];

  for (iDV = 0; iDV  < nDV_Total; iDV++){
    Gradient[iDV] = 0.0;
    Gradient_Old[iDV] = 0.0;
    DesignVarUpdate[iDV] = 0.0;
    DesignVariable[iDV] = 0.0;
    SearchDirection[iDV] = 0.0;
  }

  /*----- calculate values for bound projection algorithm -------*/
  lb=-config->GetBound()*config->GetDesignScale();
  ub=config->GetBound()*config->GetDesignScale();

  grid_movement[ZONE_0][INST_0] = new CVolumetricMovement(geometry, config);
  surface_movement[ZONE_0] = new CSurfaceMovement();

  update = false;

}

COneShotFluidDriver::~COneShotFluidDriver(void){

  /*----- free allocated memory -------*/
  delete [] Gradient;
  delete [] Gradient_Old;

  delete [] DesignVarUpdate;
  delete [] DesignVariable;
  delete [] SearchDirection;

}

void COneShotFluidDriver::Preprocess(unsigned long TimeIter) {

  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- NOTE: Inv Design Routines moved to CDiscAdjFluidIteration::Preprocess ---*/

  /*--- Preprocess the adjoint iteration ---*/

  iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                        solver_container, numerics_container, config_container,
                        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Check what is a good initialization for this recording ---*/
  /*
  if (RecordingState != MainVariables){
    MainRecording();
  }
  */

}


void COneShotFluidDriver::Run(){

  config->SetIntIter(TimeIter);

  /*--- Run an iteration of the one-shot solver ---*/
  RunOneShot();

  /*--- Screen output ---*/
  bool steady = (config->GetUnsteady_Simulation() == STEADY);
  bool output_history = false;

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime - StartTime;

  /*--- If convergence was reached --*/
  StopCalc = integration[ADJFLOW_SOL]->GetConvergence();

  /*--- Store objective and constraint for screen/history output ---*/
  solver[ADJFLOW_SOL]->SetObjFunc_Value(ObjFunc);

  /*
  if(nConstr > 0) {
    solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->SetConFunc_Value(0.0);
    for(unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
      solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->AddConFunc_Value(ConstrFunc[iConstr]);
    }
  }
  */

  /*--- Write the convergence history for the fluid (only screen output) ---*/

  /*--- Write an output ----*/
  // look into new output structure and fix this.
  output_history = (steady && !((config->GetnInner_Iter()==1)));
  //if (output_history) output->SetConvHistory_Body(NULL, geometry_container, solver_container,
  //                                                config_container, integration_container, false, UsedTime, ZONE_0, INST_0);

}

void COneShotFluidDriver::RunOneShot(){

  su2double stepsize = config->GetStepSize();
  unsigned short maxcounter = config->GetOneShotMaxCounter();
  unsigned short whilecounter = 0;

  /*--- Store the old solution and the old design for line search ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[ADJFLOW_SOL]->StoreSolution();
    solver_container[ADJFLOW_SOL]->StoreMeshPoints(config, geometry);
  }

  /*--- This is the main One Shot loop ---*/
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

  /*--- Calculate Lagrangian function ---*/
  // CalculateLagrangian(true);

  if(TimeIter >= config->GetOneShotStart() && TimeIter < config->GetOneShotStop()){

    /*--- Update design variable ---*/
    UpdateDesignVariable();

    ComputePreconditioner();

    /*--- Compute the search direction for the line search procedure ---*/
    ComputeSearchDirection();
  }
}

void COneShotFluidDriver::PrimalDualStep(){

  unsigned short iZone = 0;
  unsigned short iInst = 0;

  /*--- Note: Unsteady cases not applicable to the one-shot method yet! ---*/

  SetRecording(SOLUTION_AND_MESH);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  config->SetIntIter(0);
  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, iZone, iInst);


  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  // SetAdj_ConstrFunction(Multiplier);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->Iterate(output, integration_container, geometry_container,
                     solver_container, numerics_container, config_container,
                     surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);


  /*--- Extract the computed sensitivity values. ---*/
  solver[ADJFLOW_SOL]->SetSensitivity(geometry, solver, config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();
}

void COneShotFluidDriver::SetRecording(unsigned short kind_recording){

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if (kind_recording == SOLUTION_AND_MESH){

    AD::StartRecording();

    if (rank == MASTER_NODE && kind_recording == MainVariables && (TimeIter == 0)) {
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

  iteration->RegisterOutput(solver_container, geometry_container, config_container, output, ZONE_0, INST_0);

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  // Enable this option later
  //SetConstrFunction();

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

void COneShotFluidDriver::ComputeSearchDirection(){
 // should be something like preconditioner * gradient + constraint projection
}

bool COneShotFluidDriver::CheckDescent(){
  // to be filled
}

void COneShotFluidDriver::UpdateDesignVariable(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    DesignVariable[iDV] += DesignVarUpdate[iDV];
  }
  update = true;
}

void COneShotFluidDriver::ComputePreconditioner(){
  // here we will use the Sobolev Hessian approximation
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

