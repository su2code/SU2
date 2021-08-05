/*!
 * \file COneShotSinglezoneDriver.cpp
 * \brief The main subroutines for one-shot problems.
 * \author T.Dick
 * \version 7.1.1 "Blackbird"
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

#include <chrono>

#include "../../SU2_CFD/include/output/CMeshOutput.hpp"

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
    delta_design.push_back(0.0);
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

  /// prepare the output
  if (rank==MASTER_NODE) {
    WriteOneShotHistoryHeader();
  }

  /// main loop for the optimization method
  while (OneShotIter<config->GetOneShotIter()) {

    /// reset the convergence flags and reset the recording before running piggyback again
    StopCalc=false;
    StopNext=false;

    /* We need to record with NONE as input to ensure that all auxiliary variable dependencies
     * are removed for future runs. */
    SetRecording(RECORDING::CLEAR_INDICES);

    /// call to Piggyback
    PiggyBack();

    /// get the function value, gradient, hessian
    CombinedFunc = ComputeCombFunction();
    ComputePreconditioner();

    /// do a linesearch and design update
    Linesearch(CombinedFunc, gradient, config);

    /// update for next iteration
    UpdateDesignVariable();

    if (rank==MASTER_NODE) {
      WriteOneShotHistory(CombinedFunc);
    }

    /// check for convergence
    if (isconverged) { break; }
    OneShotIter++;
  }

  /// after the loop finishes store the final design
  OutputDesign("finaldesign.txt");

  /// for debugging purposes output the final mesh
  COutput* mesh_output = new CMeshOutput(config, geometry->GetnDim());
  mesh_output->PreprocessVolumeOutput(config);
  mesh_output->PreprocessHistoryOutput(config, false);
  mesh_output->Load_Data(geometry, config, nullptr);
  mesh_output->WriteToFile(config, geometry, MESH, "mesh_oneshot_sol");

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

  RECORDING kind_recording;
  if (iPiggyIter == nPiggyIter-1 || StopNext) {
    kind_recording = RECORDING::SOLUTION_AND_MESH;
  } else {
    kind_recording = RECORDING::SOLUTION_VARIABLES;
  }

  /*--- Store the computational graph of one direct iteration. ---*/

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

  if (RecordingState == RECORDING::SOLUTION_AND_MESH) {
    solver[MainSolver]->SetSensitivity(geometry, config);
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

void COneShotSinglezoneDriver::SetRecording(RECORDING kind_recording){

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if (kind_recording != RECORDING::CLEAR_INDICES){

    AD::StartRecording();

    if (rank == MASTER_NODE && (TimeIter == 0)) {
      cout << "Direct iteration to store the primal computational graph." << endl;
    }
    if (rank == MASTER_NODE && kind_recording == RECORDING::SOLUTION_AND_MESH) {
      cout << "Combined recording of flow and design variables." << endl;
    }

    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  }

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0, INST_0, kind_recording);

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun(kind_recording);

  RecordingState = kind_recording;

  iteration->RegisterOutput(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Evaluate the objective function and store it --- */

  SetObjFunction();
  // evaluate constraints as well if there are some.
  SetConstrFunction();

  AD::StopRecording();

}

void COneShotSinglezoneDriver::SetConstrFunction(){

  unsigned short Kind_ConstrFunc, nConstr=config->GetnConstr();
  su2double FunctionValue = 0.0;
  unsigned int iPlane=0;

  for (auto iConstr=0; iConstr < nConstr; iConstr++){

    ConstrFunc[iConstr] = 0.0;
    Kind_ConstrFunc = config->GetKind_ConstrFunc(iConstr);
    if (Kind_ConstrFunc != THICKNESS_CONSTRAINT) {
      FunctionValue = solver[FLOW_SOL]->Evaluate_ConstrFunc(config, iConstr);
      // The sign in this equation is a matter of convention, just ensure, that you choose the multiplier accordingly.
      ConstrFunc[iConstr] = config->GetConstrFuncScale(iConstr)*(FunctionValue - config->GetConstrFuncTarget(iConstr));
      if (rank == MASTER_NODE) {
        AD::RegisterOutput(ConstrFunc[iConstr]);
      }
    } else if (Kind_ConstrFunc == THICKNESS_CONSTRAINT) {
      FunctionValue = solver[MainSolver]->EvaluateGeometryFunction(geometry, config, iPlane);
      ConstrFunc[iConstr] = config->GetConstrFuncScale(iConstr)*(FunctionValue - config->GetConstrFuncTarget(iConstr));
      iPlane++;
    }
  }

}

void COneShotSinglezoneDriver::SetAdj_ConstrFunction(vector<su2double> seeding){

  for (auto iConstr=0; iConstr < config->GetnConstr(); iConstr++){
    // For geometric constraint the gradient is obtained using FD
    if (config->GetKind_ConstrFunc(iConstr) != THICKNESS_CONSTRAINT) {
      if (rank == MASTER_NODE){
        SU2_TYPE::SetDerivative(ConstrFunc[iConstr], SU2_TYPE::GetValue(seeding[iConstr]));
      }
    }
  }

}

void COneShotSinglezoneDriver::DeformGeometry(vector<su2double>& deltaVector, CConfig *config) {

  config->SetKind_SU2(SU2_COMPONENT::SU2_DEF);

  unsigned short iDV, iDV_Value, nDV_Value;
  unsigned long iDVtotal=0;

  for (iDV = 0; iDV < config->GetnDV_Total(); iDV++) {
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      config->SetDV_Value(iDV,iDV_Value, deltaVector[iDVtotal]);
      iDVtotal++;
    }
  }

  /*--- Surface grid deformation using design variables ---*/
  surface_movement[ZONE_0]->SetSurface_Deformation(geometry, config);

  /*--- Volumetric grid deformation/transformations ---*/
  if (config->GetDesign_Variable(0) != FFD_SETTING) {
    grid_movement[ZONE_0][INST_0]->SetVolume_Deformation(geometry, config, true);
  }

  config->SetKind_SU2(SU2_COMPONENT::SU2_CFD);

}

void COneShotSinglezoneDriver::Linesearch(su2double funcValue, vector<su2double> funcGrad, CConfig *config) {

  su2double length, maxLength=config->GetMaxOneShotStepsize();
  bool isDescend=false;
  su2double factor=-1.0;

  /*! Just scale the design step to avoid to large steps. */
  if (!config->GetCheckAndReset()) {

    /// Reduce to max length and change the sign (gradient descend!)
    length = L2Norm(funcGrad);
    if (length > maxLength) {
      factor = -(maxLength/length);
    }
    for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
      delta_design[iDV] = factor*funcGrad[iDV];
    }

    /// Deform the geometry
    DeformGeometry(delta_design, config);

  /*! Check wether or not we have descend. */
  } else {

    /// get original design change
    for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
      delta_design[iDV] = factor*funcGrad[iDV];
    }

    /// Store geometry and solution
    StoreSolutionAndMesh();

    /// backtracking loop
    while (!isDescend) {

      /// update the design
      DeformGeometry(delta_design, config);

      /* We need to record with RECORDING::CLEAR_INDICES as input to ensure that all auxiliary variable dependencies
       * are removed for future runs. */
      SetRecording(RECORDING::CLEAR_INDICES);

      /// call to Piggyback
      PiggyBack();

      /// Check descend
      isDescend = CheckDescent();
      if (isDescend) { break; }

      /// reduce size of design change
      for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
        delta_design[iDV] = delta_design[iDV]/10;
      }
    }

  }
}

bool COneShotSinglezoneDriver::CheckDescent(){
  return true;
}

void COneShotSinglezoneDriver::UpdateDesignVariable() {
  for (auto iDV=0; iDV<config->GetnDV_Total(); iDV++) {
    design[iDV] = BoundProjection(design[iDV] + delta_design[iDV]);
  }
}

su2double COneShotSinglezoneDriver::ComputeCombFunction() {

  su2double CombinedFunc = ObjFunc;
  for (auto iConstr=0; iConstr < config->GetnConstr(); iConstr++){
    // adapt the multiplier sign
    if ((ConstrFunc[iConstr]>=0 && multiplier[iConstr]<0) ||
         (ConstrFunc[iConstr]<0 && multiplier[iConstr]>=0)) {
      multiplier[iConstr] = -multiplier[iConstr];
    }
    CombinedFunc+=multiplier[iConstr]*ConstrFunc[iConstr];
  }
  return CombinedFunc;
}

void COneShotSinglezoneDriver::ComputePreconditioner() {

   vector<su2double> geoGrad;
   vector<su2double> storage(config->GetnDV_Total(), 0.0);
   unsigned int iDV, iDV_Value, iDV_index, nDV_Value;
   unsigned int nDV = config->GetnDV();
   su2double** DV_Value = config->GetDV_Pointer();

   /*--- Add the geometric gradients if needed ---*/
   for (auto iConstr=0; iConstr < config->GetnConstr(); iConstr++){
     if (config->GetKind_ConstrFunc(iConstr)==THICKNESS_CONSTRAINT) {

       iDV_index=0;
       for (iDV = 0; iDV  < nDV; iDV++){
         nDV_Value =  config->GetnDV_Value(iDV);
         for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
           storage[iDV_index] = DV_Value[iDV][iDV_Value];
           DV_Value[iDV][iDV_Value] = 0.001;
           iDV_index++;
         }
       }

       geoGrad = solver[MainSolver]->EvaluateGeometryGradient(geometry, surface_movement[ZONE_0], config);

       iDV_index=0;
       for (iDV = 0; iDV  < nDV; iDV++){
         nDV_Value =  config->GetnDV_Value(iDV);
         for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
           DV_Value[iDV][iDV_Value] = storage[iDV_index];
           iDV_index++;
         }
       }

       break;
     }
   }

  /*--- get the sensitivities from the adjoint solver to work with ---*/
  solver[GRADIENT_SMOOTHING]->SetSensitivity(geometry,config,solver[ADJFLOW_SOL]);

  /*--- precondition gradient and extract the result. ---*/
  if (config->GetSobMode()==ONLY_GRAD) {
    solver[GRADIENT_SMOOTHING]->RecordTapeAndCalculateOriginalGradient(geometry, surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0], config);
    for (auto iDV = 0; iDV < geoGrad.size(); iDV++) {
      gradient[iDV] += geoGrad[iDV];
    }
  } else {
    solver[GRADIENT_SMOOTHING]->ApplyGradientSmoothingDV(geometry, solver[MainSolver], numerics[GRADIENT_SMOOTHING], surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0], config, geoGrad);
  }

  // get the treated gradient back
  gradient = solver[GRADIENT_SMOOTHING]->GetDeltaP();

    ofstream outGrad ("geometric_gradient" + std::to_string(rank) + ".csv");
    outGrad.precision(17);
    for (iDV = 0; iDV < geoGrad.size(); iDV++) {
      outGrad << geoGrad[iDV] << ",";
    }
    outGrad.close();


}

void COneShotSinglezoneDriver::WriteOneShotHistoryHeader() {

  std::ofstream outfile;
  outfile.open("historyoneshot.txt", std::ios_base::out); // overwrite earlier history
  outfile << "comb function value; " << "comb function gradient; " << "current design; "  << "delta design; " << "obj function value; " << "constraint value; " << "multiplier; " << endl;
  outfile.close();

}

void COneShotSinglezoneDriver::WriteOneShotHistory(su2double &combFuncValue) {

  std::ofstream outfile;
  outfile.open("historyoneshot.txt", std::ios_base::app); // append instead of overwrite
  outfile.precision(15);
  outfile << std::scientific;

  outfile << combFuncValue << "; ";

  for (auto iDV = 0; iDV < gradient.size(); iDV++) {
    outfile << gradient[iDV] << ", ";
  }
  outfile << "; ";

  for (auto iDV = 0; iDV < design.size(); iDV++) {
    outfile << design[iDV] << ", ";
  }
  outfile << "; ";

  for (auto iDV = 0; iDV < delta_design.size(); iDV++) {
    outfile << delta_design[iDV] << ", ";
  }
  outfile << "; ";

  outfile << ObjFunc << "; ";

  for (auto iConstr=0; iConstr<config->GetnConstr(); iConstr++) {
    outfile << ConstrFunc[iConstr] << ", ";
  }
  outfile << "; ";

  for (auto iConstr=0; iConstr<config->GetnConstr(); iConstr++) {
    outfile << multiplier[iConstr] << ", ";
  }
  outfile << "; ";

  outfile << endl;
  outfile.close();

}

void COneShotSinglezoneDriver::OutputDesign(string out_file) {
  unsigned iDV;
  if (rank == MASTER_NODE) {
    ofstream dsgout (out_file);
    dsgout.precision(17);
    for (iDV = 0; iDV < design.size(); iDV++) {
      dsgout << design[iDV] << ",";
    }
    dsgout.close();
  }
}

su2double COneShotSinglezoneDriver::L2Norm(vector<su2double>& vector) {
  su2double norm=0.0;
  for (auto i=0; i<vector.size(); i++) {
    norm+=vector[i]*vector[i];
  }
  return sqrt(norm);
}

su2double COneShotSinglezoneDriver::BoundProjection(su2double value) {
  if(value<-config->GetDVBound()) value = -config->GetDVBound();
  if(value> config->GetDVBound()) value =  config->GetDVBound();
  return value;
}

void COneShotSinglezoneDriver::StoreSolutionAndMesh() {
  solver[FLOW_SOL]->GetNodes()->SetSolution_Store();
  if (config->GetKind_Turb_Model() != NONE){
    solver[TURB_SOL]->GetNodes()->SetSolution_Store();
  }
  solver[MainSolver]->GetNodes()->SetSolution_Store();
  if (config->GetKind_Turb_Model() != NONE){
    solver[ADJTURB_SOL]->GetNodes()->SetSolution_Store();
  }
  solver[MainSolver]->StoreMeshPoints(geometry, config);
}

void COneShotSinglezoneDriver::LoadSolutionAndMesh() {
  solver[FLOW_SOL]->GetNodes()->GetSolution_Store();
  if (config->GetKind_Turb_Model() != NONE){
    solver[TURB_SOL]->GetNodes()->GetSolution_Store();
  }
  solver[MainSolver]->GetNodes()->GetSolution_Store();
  if (config->GetKind_Turb_Model() != NONE){
    solver[ADJTURB_SOL]->GetNodes()->GetSolution_Store();
  }
  solver[MainSolver]->LoadMeshPoints(geometry, config);
  solver[MainSolver]->UpdateAuxiliaryGeometryVariables(geometry_container[ZONE_0][INST_0],grid_movement[ZONE_0][INST_0],config);
}
