/*!
 * \file COneShotSinglezoneDriver.cpp
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

#include "../../include/drivers/COneShotSinglezoneDriver.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutputLegacy.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"

COneShotSinglezoneDriver::COneShotSinglezoneDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator) :
    CDiscAdjSinglezoneDriver(confFile, val_nZone, MPICommunicator), StopNext(false) {

  /*--- Store the pointers ---*/
  /* Since this class inherits from CDiscAdjSinglezoneDriver relevant pointers to
   * config, solver, iteration, geometry, numerics in the ZONE_0 should already be available
   */

  /*---------- One-shot works on all design variables - get the total number of design variables from Config ---------*/
  nDV_Total = config->GetnDV_Total();

  /*--- Preprocess the additional COutput class for the flow output ---*/
  flowoutput = COutputFactory::CreateOutput(EULER, config, nDim);
  flowoutput->PreprocessHistoryOutput(config, false);
  flowoutput->PreprocessVolumeOutput(config);

  /*--- Initialize the members of the driver class. ---*/
  OneShotIter=0;

  for (auto iConstr=0; iConstr<config->GetnConstr(); iConstr++) {
    ConstrFunc.push_back(0.0);
    multiplier.push_back(config->GetInitialMultiplier(iConstr));
  }

  for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
    gradient.push_back(0.0);
    design.push_back(0.0);
  }

}

COneShotSinglezoneDriver::~COneShotSinglezoneDriver(void){

}

void COneShotSinglezoneDriver::Preprocess(unsigned long TimeIter) {

  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- NOTE: Inv Design Routines moved to CDiscAdjFluidIteration::Preprocess ---*/

  /*--- Preprocess the adjoint iteration ---*/

  iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                        solver_container, numerics_container, config_container,
                        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

}

void COneShotSinglezoneDriver::Run(){
 if(config->GetOneShotMode() == PIGGYBACK) {
   PiggyBack();
 } else if (config->GetOneShotMode() == ONESHOT) {
   RunOneShot();
 } else if (config->GetOneShotMode() == NO_MODE) {
   if (rank == MASTER_NODE)  {
     cout << "Error: OneShot driver called without a mode!" << endl;
   }
 }
}

void COneShotSinglezoneDriver::Postprocess() {

  if(config->GetOneShotMode() == PIGGYBACK) {

    /*--- For Piggyback smooth the calculated geometry sensitivities if required. ---*/
    if (config->GetSmoothGradient()) {
      DerivativeTreatment();
    }
  }

}

void COneShotSinglezoneDriver::RunOneShot(){

  su2double CombinedFunc=0.0;
  bool isconverged=false;

  /// main loop for the optimization method
  while (OneShotIter<config->GetOneShotIter()) {

    /// call to Piggyback
    PiggyBack();

    /// get the function value, gradient, hessian
    CombinedFunc = ObjFunc;
    for (auto iConstr=0; iConstr < config->GetnConstr(); iConstr++){
      CombinedFunc+=multiplier[iConstr]*ConstrFunc[iConstr];
    }
    ComputePreconditioner();

    /// do a linesearch and design update
    Linesearch(CombinedFunc, gradient, config);

    /// check for convergence and update for next iteration

    UpdateDesignVariable();
    if (isconverged) { break; }
    OneShotIter++;
  }

}

void COneShotSinglezoneDriver::PiggyBack() {

  // note: use inherited iteration count from the adjoint driver for now.
  nPiggyIter = nAdjoint_Iter;

  // in the first and last loop we might want to do more iterations to be more accurate.
  if (OneShotIter==0 || OneShotIter==(config->GetOneShotIter()-1)) {
    nPiggyIter += config->GetAddInnerIter();
  }

  /*--- Main loop for the Piggyback iteration ---*/
  for (auto PiggyIter = 0ul; PiggyIter < nPiggyIter; PiggyIter++) {

    /*--- Run the PrimalDualStep for Piggyback iteration --*/

    PrimalDualStep(PiggyIter);

    /*--- If maximum number of iteration is reached or the iteration converged, force flow output. ---*/

    if (PiggyIter == nPiggyIter-1 || StopCalc) {
      flowoutput->SetResult_Files(geometry, config, solver, PiggyIter, true);
    }

    /*--- If the convergence detection in PrimalDualStep succeded, then stop ---*/

    if (StopCalc) break;

  }

}

void COneShotSinglezoneDriver::PrimalDualStep(unsigned long iPiggyIter){

  /*--- In the last Piggyback step we need to calculate the mesh sensitivities as well ---*/

  unsigned short kind_recording;
  if (iPiggyIter == nPiggyIter-1 || StopNext) {
    kind_recording = SOLUTION_AND_MESH;
  } else {
    if (config->GetOneShotMode()==PIGGYBACK) {
      kind_recording = SOLUTION_VARIABLES;
    } else {
      kind_recording = SOLUTION_AND_MESH;
    }
  }

  SetRecording(kind_recording);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  config->SetInnerIter(iPiggyIter);

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);


  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  SetAdj_ConstrFunction(multiplier);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/

  iteration->Iterate(output_container[ZONE_0], integration_container, geometry_container,
                     solver_container, numerics_container, config_container,
                     surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Extract the computed sensitivity values. ---*/

  if (RecordingState == SOLUTION_AND_MESH) {
    solver[ADJFLOW_SOL]->SetSensitivity(geometry, solver, config);
  }

  /*--- Monitor the pseudo-time ---*/
  /* if we are converged we need one more step to record and tape everything correct before we stop!
   * so we use an intermediate flag to do one more iteration */

  StopCalc = StopNext;
  StopNext = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                solver_container, numerics_container, config_container,
                                surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

  /*--- Output files for steady state simulations. ---*/

  if (!config->GetTime_Domain()) {
    iteration->Output(output_container[ZONE_0], geometry_container, solver_container,
                      config_container, iPiggyIter, false, ZONE_0, INST_0);
    // for OneShot set flow output aswell
    flowoutput->SetResult_Files(geometry, config, solver, iPiggyIter);
  }

}

void COneShotSinglezoneDriver::SetRecording(unsigned short kind_recording){

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if (kind_recording != NONE){

    AD::StartRecording();

    if (rank == MASTER_NODE && (TimeIter == 0)) {
      cout << "Direct iteration to store the primal computational graph." << endl;
    }
    if (rank == MASTER_NODE && kind_recording == SOLUTION_AND_MESH) {
      cout << "Combined recording of flow and design variables." << endl;
    }

    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  }

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0, INST_0, kind_recording);

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun(kind_recording);

  RecordingState = kind_recording;

  iteration->RegisterOutput(solver_container, geometry_container, config_container, output_container[ZONE_0], ZONE_0, INST_0);

  /*--- Evaluate the objective function and store it --- */

  SetObjFunction();
  // evaluate constraints as well if there are some.
  SetConstrFunction();

  AD::StopRecording();

}

void COneShotSinglezoneDriver::SetConstrFunction(){

  unsigned short Kind_ConstrFunc, nConstr=config->GetnConstr();
  su2double FunctionValue = 0.0;

  if (rank == MASTER_NODE && nConstr != 0) { std::cout << "Constraint Function Value: "; }

  for (auto iConstr=0; iConstr < nConstr; iConstr++){

    ConstrFunc[iConstr] = 0.0;
    Kind_ConstrFunc = config->GetKind_ConstrFunc(iConstr);
    FunctionValue = solver[FLOW_SOL]->Evaluate_ConstrFunc(config, iConstr);
    // The sign in this equation is a matter of cenvention, just ensure, that you choose the multiplier accordingly.
    ConstrFunc[iConstr] = config->GetConstrFuncScale(iConstr)*(FunctionValue - config->GetConstrFuncTarget(iConstr));
    if (rank == MASTER_NODE){
      std::cout << FunctionValue << ", " << config->GetConstrFuncScale(iConstr) << ", " << multiplier[iConstr] << "; ";
      AD::RegisterOutput(ConstrFunc[iConstr]);
    }
  }

  if (rank == MASTER_NODE && nConstr != 0) { std::cout<<std::endl; }
}

void COneShotSinglezoneDriver::SetAdj_ConstrFunction(vector<su2double> seeding){

  for (auto iConstr=0; iConstr < config->GetnConstr(); iConstr++){
    if (rank == MASTER_NODE){
      SU2_TYPE::SetDerivative(ConstrFunc[iConstr], SU2_TYPE::GetValue(seeding[iConstr]));
    }
  }

}

void COneShotSinglezoneDriver::DeformGeometry(vector<su2double>& delta_design, CConfig *config) {

  config->SetKind_SU2(SU2_DEF);

  unsigned short iDV, iDV_Value, nDV_Value;
  unsigned long iDVtotal=0;

  for (iDV = 0; iDV < config->GetnDV_Total(); iDV++) {
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      config->SetDV_Value(iDV,iDV_Value, delta_design[iDVtotal]);
      iDVtotal++;
    }
  }

  /*--- Surface grid deformation using design variables ---*/
  surface_movement[ZONE_0]->SetSurface_Deformation(geometry, config);

  /*--- Volumetric grid deformation/transformations ---*/
  if (config->GetDesign_Variable(0) != FFD_SETTING) {
    grid_movement[ZONE_0][INST_0]->SetVolume_Deformation(geometry, config, true);
  }

  config->SetKind_SU2(SU2_CFD);

}

void COneShotSinglezoneDriver::Linesearch(su2double funcValue, vector<su2double> funcGrad, CConfig *config) {

  su2double length, maxLength=config->GetMaxOneShotStepsize();
  vector<su2double> delta_design(funcGrad.size(), 0.0);
  bool isDescent=false;

  /// Reduce to max length and change the sign (gradient descend!)
  length = L2Norm(funcGrad);
  if (length > maxLength) {
    for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
      delta_design[iDV] = -(maxLength/length)*funcGrad[iDV];
    }
  }

  /// Store geometry and solution

  /// Deform the geometry
  DeformGeometry(delta_design, config);

  /// Update the function

  /// Check descend
  isDescent = CheckDescent();

  /// Reset geometry and solution

  /// Set return gradient
  for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
    gradient[iDV] = delta_design[iDV];
  }

}

bool COneShotSinglezoneDriver::CheckDescent(){
  return true;
}

void COneShotSinglezoneDriver::UpdateDesignVariable() {
  for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
    design[iDV] -= gradient[iDV];
  }
}

void COneShotSinglezoneDriver::ComputePreconditioner() {

  /*--- get the sensitivities from the adjoint solver to work with ---*/
  solver[GRADIENT_SMOOTHING]->SetSensitivity(geometry,solver,config);
  /*--- precondition gradient and extract thew result. ---*/
  solver[GRADIENT_SMOOTHING]->ApplyGradientSmoothingDV(geometry, solver[ADJFLOW_SOL], numerics[GRADIENT_SMOOTHING], surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0], config);
  gradient = solver[GRADIENT_SMOOTHING]->GetDeltaP();
}

su2double COneShotSinglezoneDriver::L2Norm(vector<su2double>& vector) {
  su2double norm=0.0;
  for (auto i=0; i<vector.size(); i++) {
    norm+=vector[i]*vector[i];
  }
  return sqrt(norm);
}
