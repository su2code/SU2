/*!
 * \file driver_adjoint_singlezone.cpp
 * \brief The main subroutines for driving adjoint single-zone problems.
 * \author R. Sanchez
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

#include "../../include/drivers/CDiscAdjSinglezoneDriver.hpp"

CDiscAdjSinglezoneDriver::CDiscAdjSinglezoneDriver(char* confFile,
                                                   unsigned short val_nZone,
                                                   SU2_Comm MPICommunicator) : CSinglezoneDriver(confFile,
                                                                                                 val_nZone,
                                                                                                 MPICommunicator) {


  /*--- Store the number of internal iterations that will be run by the adjoint solver ---*/
  nAdjoint_Iter = config_container[ZONE_0]->GetnInner_Iter();
  

  /*--- Store the pointers ---*/
  config      = config_container[ZONE_0];
  iteration   = iteration_container[ZONE_0][INST_0];
  solver      = solver_container[ZONE_0][INST_0][MESH_0];
  numerics    = numerics_container[ZONE_0][INST_0][MESH_0];
  geometry    = geometry_container[ZONE_0][INST_0][MESH_0];
  integration = integration_container[ZONE_0][INST_0];

  /*--- Store the recording state ---*/
  RecordingState = NONE;

  /*--- Determine if the problem is a turbomachinery problem ---*/
  bool turbo = config->GetBoolTurbomachinery();
  
  bool compressible = config->GetKind_Regime() == COMPRESSIBLE;

  /*--- Determine if the problem has a mesh deformation solver ---*/
  bool mesh_def = config->GetDeform_Mesh();

  /*--- Initialize the direct iteration ---*/

  switch (config->GetKind_Solver()) {

  case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;
    if (turbo) {
      direct_iteration = new CTurboIteration(config);
      output_legacy = new COutputLegacy(config_container[ZONE_0]);
    }
    else       direct_iteration = new CFluidIteration(config);
    if (compressible) direct_output = new CFlowCompOutput(config, nDim);
    else direct_output = new CFlowIncOutput(config, nDim);
    MainVariables = FLOW_CONS_VARS;
    if (mesh_def) SecondaryVariables = MESH_DEFORM;
    else          SecondaryVariables = MESH_COORDS;
    break;

  case DISC_ADJ_FEM_EULER : case DISC_ADJ_FEM_NS : case DISC_ADJ_FEM_RANS :
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;
    direct_iteration = new CFEMFluidIteration(config);
    direct_output = new CFlowCompFEMOutput(config, nDim);
    MainVariables = FLOW_CONS_VARS;
    SecondaryVariables = MESH_COORDS;
    break;

  case DISC_ADJ_FEM:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: elasticity equation." << endl;
    direct_iteration = new CFEAIteration(config);
    direct_output = new CElasticityOutput(config, nDim);
    MainVariables = FEA_DISP_VARS;
    SecondaryVariables = NONE;
    break;

  case DISC_ADJ_HEAT:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: heat equation." << endl;
    direct_iteration = new CHeatIteration(config);
    direct_output = new CHeatOutput(config, nDim);    
    MainVariables = FLOW_CONS_VARS;
    SecondaryVariables = MESH_COORDS;
    break;

  }
  
 direct_output->PreprocessHistoryOutput(config, false);

}

CDiscAdjSinglezoneDriver::~CDiscAdjSinglezoneDriver(void) {

}

void CDiscAdjSinglezoneDriver::Preprocess(unsigned long TimeIter) {
  
  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- NOTE: Inv Design Routines moved to CDiscAdjFluidIteration::Preprocess ---*/

  /*--- Preprocess the adjoint iteration ---*/

  iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                        solver_container, numerics_container, config_container,
                        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- For the adjoint iteration we need the derivatives of the iteration function with
   *--- respect to the conservative variables. Since these derivatives do not change in the steady state case
   *--- we only have to record if the current recording is different from the main variables. ---*/

  if (RecordingState != MainVariables){

    MainRecording();

  }

}

void CDiscAdjSinglezoneDriver::Run() {

  unsigned long Adjoint_Iter;

  for (Adjoint_Iter = 0; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) {

    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *--- of the previous iteration. The values are passed to the AD tool.
     *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

    config->SetInnerIter(Adjoint_Iter);

    iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/

    iteration->Iterate(output_container[ZONE_0], integration_container, geometry_container,
                         solver_container, numerics_container, config_container,
                         surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    /*--- Monitor the pseudo-time ---*/
    StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();
    
    if (StopCalc) break;

  }

}

void CDiscAdjSinglezoneDriver::Postprocess() {

  switch(config->GetKind_Solver())
  {
    case DISC_ADJ_EULER :     case DISC_ADJ_NAVIER_STOKES :     case DISC_ADJ_RANS :
    case DISC_ADJ_INC_EULER : case DISC_ADJ_INC_NAVIER_STOKES : case DISC_ADJ_INC_RANS :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();
      break;

    case DISC_ADJ_FEM :

      /*--- Apply the boundary condition to clamped nodes ---*/
      iteration->Postprocess(output_container[ZONE_0],integration_container,geometry_container,solver_container,numerics_container,
                             config_container,surface_movement,grid_movement,FFDBox,ZONE_0,INST_0);

      RecordingState = NONE;
      break;
  }//switch

  if (config->GetError_Estimate()) {
    /*--- Compute metric for anisotropic mesh adaptation ---*/
    ComputeMetric();

    output->PreprocessVolumeOutput(config);
    
    /*--- Load the data --- */
    
    direct_output->Load_Data(geometry, config, solver);
    
    /*--- Set the filenames ---*/
    
    direct_output->SetVolume_Filename(config->GetVolume_FileName());
    
    direct_output->SetSurface_Filename(config->GetSurfCoeff_FileName());
    
    for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++){
      unsigned short* FileFormat = config->GetVolumeOutputFiles();
      direct_output->WriteToFile(config, geometry, FileFormat[iFile]);

      /*--- For now we need ASCII restarts for AMGIO ---*/
      direct_output->WriteToFile(config, geometry, RESTART_ASCII);
    }
  }

}

void CDiscAdjSinglezoneDriver::SetRecording(unsigned short kind_recording){

  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/

  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the iteration --- */

  if (kind_recording != NONE){

    AD::StartRecording();

    if (rank == MASTER_NODE && kind_recording == MainVariables) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store the primal computational graph." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
    }
    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  }

  /*--- Set the dependencies of the iteration ---*/

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0,
                             INST_0, kind_recording);

  /*--- Do one iteration of the direct solver ---*/

  DirectRun(kind_recording);

  // NOTE: The inverse design calls were moved to DirectRun() - postprocess

  /*--- Store the recording state ---*/

  RecordingState = kind_recording;

  /*--- Register Output of the iteration ---*/

  iteration->RegisterOutput(solver_container, geometry_container, config_container, output_container[ZONE_0], ZONE_0, INST_0);

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  AD::StopRecording();

}

void CDiscAdjSinglezoneDriver::SetAdj_ObjFunction(){

  bool time_stepping = config->GetTime_Marching() != STEADY;
  unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
  su2double seeding = 1.0;

  if (time_stepping){
    if (TimeIter < IterAvg_Obj){
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else{
      seeding = 0.0;
    }
  }

  if (rank == MASTER_NODE){
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }

}

void CDiscAdjSinglezoneDriver::SetObjFunction(){

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

  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc);
  }

}

void CDiscAdjSinglezoneDriver::DirectRun(unsigned short kind_recording){

  /*--- Mesh movement ---*/

  direct_iteration->SetMesh_Deformation(geometry_container[ZONE_0][INST_0], solver, numerics, config, kind_recording);

  /*--- Zone preprocessing ---*/

  direct_iteration->Preprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Iterate the direct solver ---*/

  direct_iteration->Iterate(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Postprocess the direct solver ---*/

  direct_iteration->Postprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Print the direct residual to screen ---*/

  Print_DirectResidual(kind_recording);

}

void CDiscAdjSinglezoneDriver::Print_DirectResidual(unsigned short kind_recording){

  /*--- Print the residuals of the direct iteration that we just recorded ---*/
  /*--- This routine should be moved to the output, once the new structure is in place ---*/
  if ((rank == MASTER_NODE) && (kind_recording == MainVariables)){

    switch (config->GetKind_Solver()) {

    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
    case DISC_ADJ_FEM_EULER : case DISC_ADJ_FEM_NS : case DISC_ADJ_FEM_RANS :
      cout << "log10[U(0)]: "   << log10(solver[FLOW_SOL]->GetRes_RMS(0))
           << ", log10[U(1)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(1))
           << ", log10[U(2)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(2)) << "." << endl;
      cout << "log10[U(3)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(3));
      if (geometry->GetnDim() == 3) cout << ", log10[U(4)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(4));
      cout << "." << endl;
      if ( config->GetKind_Turb_Model() != NONE && !config->GetFrozen_Visc_Disc()) {
        cout << "log10[Turb(0)]: "   << log10(solver[TURB_SOL]->GetRes_RMS(0));
        if (solver[TURB_SOL]->GetnVar() > 1) cout << ", log10[Turb(1)]: " << log10(solver[TURB_SOL]->GetRes_RMS(1));
        cout << "." << endl;
      }
      if (config->GetWeakly_Coupled_Heat()){
        cout << "log10[Heat(0)]: "   << log10(solver[HEAT_SOL]->GetRes_RMS(0)) << "." << endl;
      }
      break;

    case DISC_ADJ_FEM:

      if (config->GetGeometricConditions() == LARGE_DEFORMATIONS){
        cout << "UTOL-A: "   << log10(solver[FEA_SOL]->GetRes_FEM(0))
             << ", RTOL-A: " << log10(solver[FEA_SOL]->GetRes_FEM(1))
             << ", ETOL-A: " << log10(solver[FEA_SOL]->GetRes_FEM(2)) << "." << endl;
      }
      else{
        if (geometry->GetnDim() == 2){
          cout << "log10[RMS Ux]: "   << log10(solver[FEA_SOL]->GetRes_RMS(0))
               << ", log10[RMS Uy]: " << log10(solver[FEA_SOL]->GetRes_RMS(1)) << "." << endl;
        }
        else{
          cout << "log10[RMS Ux]: "   << log10(solver[FEA_SOL]->GetRes_RMS(0))
               << ", log10[RMS Uy]: " << log10(solver[FEA_SOL]->GetRes_RMS(1))
               << ", log10[RMS Uz]: " << log10(solver[FEA_SOL]->GetRes_RMS(2))<< "." << endl;
        }
      }

      break;

    case DISC_ADJ_HEAT:
      cout << "log10[Cons(0)]: "   << log10(solver[HEAT_SOL]->GetRes_RMS(0)) << "." << endl;
      break;

    }

    cout << "-------------------------------------------------------------------------" << endl << endl;
  }
  else if ((rank == MASTER_NODE) && (kind_recording == SecondaryVariables) && (SecondaryVariables != NONE)){
    cout << endl << "Recording the computational graph with respect to the ";
    switch (SecondaryVariables){
      case MESH_COORDS: cout << "mesh coordinates." << endl;    break;
      default:          cout << "secondary variables." << endl; break;
     }
  }

}

void CDiscAdjSinglezoneDriver::MainRecording(){

  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
   *    as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(NONE);

  /*--- Store the computational graph of one direct iteration with the conservative variables as input. ---*/

  SetRecording(MainVariables);

}

void CDiscAdjSinglezoneDriver::SecondaryRecording(){

  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
   * as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(NONE);

  /*--- Store the computational graph of one direct iteration with the secondary variables as input. ---*/

  SetRecording(SecondaryVariables);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the current iteration. The values are passed to the AD tool. ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed sensitivity values. ---*/
  switch(SecondaryVariables){
  case MESH_COORDS:
    solver[ADJFLOW_SOL]->SetSensitivity(geometry, solver, config);
    break;
  case MESH_DEFORM:
    solver[ADJMESH_SOL]->SetSensitivity(geometry, solver, config);
    break;
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

}

void CDiscAdjSinglezoneDriver::ComputeMetric() {

  CSolver *solver_flow    = solver[FLOW_SOL],
          *solver_turb    = solver[TURB_SOL],
          *solver_adjflow = solver[ADJFLOW_SOL],
          *solver_adjturb = solver[ADJTURB_SOL];

  if (rank == MASTER_NODE)
    cout << endl <<"----------------------------- Compute Metric ----------------------------" << endl;

  //--- 2D
  if(nDim == 2){
    //--- Volume flow grad
    if(rank == MASTER_NODE) cout << "Computing flow volume gradient via L2 Projection." << endl;
    solver_flow->SetGradient_L2Proj2(geometry, config);

    //--- Volume flow Hess
    if(rank == MASTER_NODE) cout << "Computing flow volume Hessian via L2 Projection." << endl;
    solver_flow->SetHessian_L2Proj2(geometry, config);

    //--- Volume adj grad
    if(rank == MASTER_NODE) cout << "Computing adjoint volume gradient via L2 Projection." << endl;
    solver_adjflow->SetGradient_L2Proj2(geometry, config);

    if(config->GetViscous()) {
      //--- Volume turb grad
      if(rank == MASTER_NODE) cout << "Computing turbulent volume gradient via L2 Projection." << endl;
      solver_turb->SetTurbGradient_L2Proj2(geometry, config, solver_flow);

      //--- Volume turb Hess
      if(rank == MASTER_NODE) cout << "Computing turbulent volume Hessian via L2 Projection." << endl;
      solver_turb->SetHessian_L2Proj2(geometry, config);

      //--- Volume adj turb grad
      if(rank == MASTER_NODE) cout << "Computing turbulent adjoint volume gradient via L2 Projection." << endl;
      solver_adjturb->SetGradient_L2Proj2(geometry, config);
    }

    //--- Metric
    if(rank == MASTER_NODE) cout << "Computing goal-oriented metric tensor." << endl;
    SumWeightedHessian2(solver_flow, solver_turb, solver_adjflow, solver_adjturb, geometry);
  }

  //--- 3D
  else{
    //--- Volume flow grad
    if(rank == MASTER_NODE) cout << "Computing flow volume gradient via L2 Projection." << endl;
    solver_flow->SetGradient_L2Proj3(geometry, config);

    //--- Volume flow Hess
    if(rank == MASTER_NODE) cout << "Computing flow volume Hessian via L2 Projection." << endl;
    solver_flow->SetHessian_L2Proj3(geometry, config);

    //--- Volume adj grad
    if(rank == MASTER_NODE) cout << "Computing adjoint volume gradient via L2 Projection." << endl;
    solver_adjflow->SetGradient_L2Proj3(geometry, config);

    if(config->GetViscous()) {
      //--- Volume turb grad
      if(rank == MASTER_NODE) cout << "Computing turbulent volume gradient via L2 Projection." << endl;
      solver_turb->SetTurbGradient_L2Proj3(geometry, config, solver_flow);

      //--- Volume turb Hess
      if(rank == MASTER_NODE) cout << "Computing turbulent volume Hessian via L2 Projection." << endl;
      solver_turb->SetHessian_L2Proj3(geometry, config);

      //--- Volume adj turb grad
      if(rank == MASTER_NODE) cout << "Computing turbulent adjoint volume gradient via L2 Projection." << endl;
      solver_adjturb->SetGradient_L2Proj3(geometry, config);
    }

    //--- Metric
    if(rank == MASTER_NODE) cout << "Computing goal-oriented metric tensor." << endl;
    SumWeightedHessian3(solver_flow, solver_turb, solver_adjflow, solver_adjturb, geometry);
  }
}

void CDiscAdjSinglezoneDriver::SumWeightedHessian2(CSolver   *solver_flow,
                                                 CSolver   *solver_turb,
                                                 CSolver   *solver_adjflow,
                                                 CSolver   *solver_adjturb,
                                                 CGeometry *geometry) {

  unsigned long iPoint, 
                nPointDomain = geometry->GetnPointDomain();
  unsigned short nVarMetr = solver_flow->GetnVar(), 
                 nFluxMetr = 2;  //--- TODO: adjust size for goal (for different solvers, currently Euler) vs. feature
  unsigned short nMetr = 3*(nDim-1);

  su2double localScale = 0.0,
            globalScale = 0.0,
            p = config->GetAdap_Norm(),
            eigmax = 1./(config->GetMesh_Hmin()*config->GetMesh_Hmin()),
            eigmin = 1./(config->GetMesh_Hmax()*config->GetMesh_Hmax()),
            outComplex = su2double(config->GetMesh_Complexity());  // Constraint mesh complexity

  su2double localMinDensity = 1.E16, localMaxDensity = 0., localTotComplex = 0.;
  su2double globalMinDensity = 1.E16, globalMaxDensity = 0., globalTotComplex = 0.;

  su2double **A      = new su2double*[nDim],
            **EigVec = new su2double*[nDim], 
            *EigVal  = new su2double[nDim];

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    A[iDim]      = new su2double[nDim];
    EigVec[iDim] = new su2double[nDim];
  }

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    //--- perform summation of weighted mean flow Hessians
    for (unsigned short iVar = 0; iVar < nVarMetr; ++iVar) {

      for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
        const unsigned short ig = iVar*nDim + iFlux;
        const su2double grad = solver_adjflow->GetNodes()->GetAnisoGrad(iPoint, ig);

        for (unsigned short im = 0; im < nMetr; ++im) {
          const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
          const su2double hess = solver_flow->GetNodes()->GetAnisoHess(iPoint, ih);
          const su2double part = abs(grad)*hess;
          solver_flow->GetNodes()->AddAnisoMetr(iPoint, im,part);
        }
      }
    }

    //--- add viscous Hessian terms
    if(config->GetViscous()) {
      //--- viscous mass flux is 0, so start with momentum
      for (unsigned short iVar = 1; iVar < nVarMetr; ++iVar) {

        for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
          const unsigned short ig = iVar*nDim + iFlux;
          const su2double grad = solver_adjflow->GetNodes()->GetAnisoGrad(iPoint, ig);

          for (unsigned short im = 0; im < nMetr; ++im) {
            const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
            const su2double hess = solver_flow->GetNodes()->GetAnisoViscHess(iPoint, ih);
            const su2double part = abs(grad)*hess;
            solver_flow->GetNodes()->AddAnisoMetr(iPoint, im,part);
          }
        }
      }

      //--- add turbulent terms
      const unsigned short nVarTurbMetr = solver_turb->GetnVar();
      for (unsigned short iVar = 0; iVar < nVarTurbMetr; ++iVar) {

        for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
          const unsigned short ig = iVar*nDim + iFlux;
          const su2double grad = solver_adjturb->GetNodes()->GetAnisoGrad(iPoint, ig);

          for (unsigned short im = 0; im < nMetr; ++im) {
            const unsigned short ih = iFlux*nVarTurbMetr*nMetr + iVar*nMetr + im;  
            const su2double hess = solver_turb->GetNodes()->GetAnisoHess(iPoint, ih) 
                                 + solver_turb->GetNodes()->GetAnisoViscHess(iPoint, ih);
            const su2double part = abs(grad)*hess;
            solver_flow->GetNodes()->AddAnisoMetr(iPoint, im,part);
          }
        }
      }
    }
  }

  //--- communicate the metric values via MPI, then correct boundary terms
  
  if(rank == MASTER_NODE) cout << "Correcting metric." << endl;
  solver_flow->InitiateComms(geometry, config, ANISO_METRIC);
  solver_flow->CompleteComms(geometry, config, ANISO_METRIC);
  solver_flow->CorrectBoundAnisoMetr(geometry, config);

  //--- set tolerance and obtain global scaling
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const su2double a = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 0);
    const su2double b = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 1);
    const su2double c = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 2);
    
    A[0][0] = a; A[0][1] = b;
    A[1][0] = b; A[1][1] = c;

    CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

    const su2double Vol = geometry->node[iPoint]->GetVolume();

    localScale += pow(abs(EigVal[0]*EigVal[1]),p/(2.*p+nDim))*Vol;
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&localScale, &globalScale, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  globalScale = localScale;
#endif

  //--- normalize to achieve Lp metric for constraint complexity, then truncate size
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const su2double a = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 0);
    const su2double b = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 1);
    const su2double c = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 2);
    
    A[0][0] = a; A[0][1] = b;
    A[1][0] = b; A[1][1] = c;

    CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

    const su2double factor = pow(outComplex/globalScale, 2./nDim) * pow(abs(EigVal[0]*EigVal[1]), -1./(2.*p+nDim));

    for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = min(max(abs(factor*EigVal[iDim]),eigmin),eigmax);

    CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 0, A[0][0]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 1, A[0][1]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 2, A[1][1]);

    //--- compute min, max, total complexity
    const su2double Vol = geometry->node[iPoint]->GetVolume();
    const su2double density = sqrt(abs(EigVal[0]*EigVal[1]));

    localMinDensity = min(localMinDensity, density);
    localMaxDensity = max(localMaxDensity, density);
    localTotComplex += density*Vol;
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&localMinDensity, &globalMinDensity, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&localMaxDensity, &globalMaxDensity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&localTotComplex, &globalTotComplex, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  globalMinDensity = localMinDensity;
  globalMaxDensity = localMaxDensity;
  globalTotComplex = localTotComplex;
#endif

  if(rank == MASTER_NODE) {
    cout << "Minimum density: " << globalMinDensity << "." << endl;
    cout << "Maximum density: " << globalMaxDensity << "." << endl;
    cout << "Mesh complexity: " << globalTotComplex << "." << endl;
  }

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    delete [] A[iDim];
    delete [] EigVec[iDim];
  }
  delete [] A;
  delete [] EigVec;
  delete [] EigVal;
}

void CDiscAdjSinglezoneDriver::SumWeightedHessian3(CSolver   *solver_flow,
                                                 CSolver   *solver_turb,
                                                 CSolver   *solver_adjflow,
                                                 CSolver   *solver_adjturb,
                                                 CGeometry *geometry) {

  unsigned long iPoint, nPointDomain = geometry->GetnPointDomain();
  unsigned short nVarMetr = solver_flow->GetnVar(), 
                 nFluxMetr = 3;  //--- TODO: adjust size for goal (for different solvers, currently Euler) vs. feature
  unsigned short nMetr = 3*(nDim-1);

  su2double localScale = 0.0,
            globalScale = 0.0,
            p = config->GetAdap_Norm(),
            eigmax = 1./(config->GetMesh_Hmin()*config->GetMesh_Hmin()),
            eigmin = 1./(config->GetMesh_Hmax()*config->GetMesh_Hmax()),
            outComplex = su2double(config->GetMesh_Complexity());  // Constraint mesh complexity

  su2double localMinDensity = 1.E16, localMaxDensity = 0., localTotComplex = 0.;
  su2double globalMinDensity = 1.E16, globalMaxDensity = 0., globalTotComplex = 0.;

  su2double **A      = new su2double*[nDim],
            **EigVec = new su2double*[nDim], 
            *EigVal  = new su2double[nDim];

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    A[iDim]      = new su2double[nDim];
    EigVec[iDim] = new su2double[nDim];
  }

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    //--- perform summation of weighted mean flow Hessians
    for (unsigned short iVar = 0; iVar < nVarMetr; ++iVar) {

      for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
        const unsigned short ig = iVar*nDim + iFlux;
        const su2double grad = solver_adjflow->GetNodes()->GetAnisoGrad(iPoint, ig);

        for (unsigned short im = 0; im < nMetr; ++im) {
          const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
          const su2double hess = solver_flow->GetNodes()->GetAnisoHess(iPoint, ih);
          const su2double part = abs(grad)*hess;
          solver_flow->GetNodes()->AddAnisoMetr(iPoint, im,part);
        }
      }
    }

    //--- add viscous Hessian terms
    if(config->GetViscous()) {
      //--- viscous mass flux is 0, so start with momentum
      for (unsigned short iVar = 1; iVar < nVarMetr; ++iVar) {

        for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
          const unsigned short ig = iVar*nDim + iFlux;
          const su2double grad = solver_adjflow->GetNodes()->GetAnisoGrad(iPoint, ig);

          for (unsigned short im = 0; im < nMetr; ++im) {
            const unsigned short ih = iFlux*nVarMetr*nMetr + iVar*nMetr + im;  
            const su2double hess = solver_flow->GetNodes()->GetAnisoViscHess(iPoint, ih);
            const su2double part = abs(grad)*hess;
            solver_flow->GetNodes()->AddAnisoMetr(iPoint, im,part);
          }
        }
      }

      //--- add turbulent terms
      const unsigned short nVarTurbMetr = solver_turb->GetnVar();
      for (unsigned short iVar = 0; iVar < nVarTurbMetr; ++iVar) {

        for (unsigned short iFlux = 0; iFlux < nFluxMetr; ++iFlux) {
          const unsigned short ig = iVar*nDim + iFlux;
          const su2double grad = solver_adjturb->GetNodes()->GetAnisoGrad(iPoint, ig);

          for (unsigned short im = 0; im < nMetr; ++im) {
            const unsigned short ih = iFlux*nVarTurbMetr*nMetr + iVar*nMetr + im;  
            const su2double hess = solver_turb->GetNodes()->GetAnisoHess(iPoint, ih) 
                                 + solver_turb->GetNodes()->GetAnisoViscHess(iPoint, ih);
            const su2double part = abs(grad)*hess;
            solver_flow->GetNodes()->AddAnisoMetr(iPoint, im,part);
          }
        }
      }
    }
  }

  //--- communicate the metric values via MPI, then correct boundary terms
  
  if(rank == MASTER_NODE) cout << "Correcting metric." << endl;
  solver_flow->InitiateComms(geometry, config, ANISO_METRIC);
  solver_flow->CompleteComms(geometry, config, ANISO_METRIC);
  solver_flow->CorrectBoundAnisoMetr(geometry, config);

  //--- set tolerance and obtain global scaling
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const su2double a = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 0);
    const su2double b = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 1);
    const su2double c = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 2);
    const su2double d = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 3);
    const su2double e = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 4);
    const su2double f = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 5);

    A[0][0] = a; A[0][1] = b; A[0][2] = c;
    A[1][0] = b; A[1][1] = d; A[1][2] = e;
    A[2][0] = c; A[2][1] = e; A[2][2] = f;

    CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

    const su2double Vol = geometry->node[iPoint]->GetVolume();

    localScale += pow(abs(EigVal[0]*EigVal[1]*EigVal[2]),p/(2.*p+nDim))*Vol;
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&localScale, &globalScale, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  globalScale = localScale;
#endif

  //--- normalize to achieve Lp metric for constraint complexity, then truncate size
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const su2double a = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 0);
    const su2double b = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 1);
    const su2double c = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 2);
    const su2double d = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 3);
    const su2double e = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 4);
    const su2double f = solver_flow->GetNodes()->GetAnisoMetr(iPoint, 5);

    A[0][0] = a; A[0][1] = b; A[0][2] = c;
    A[1][0] = b; A[1][1] = d; A[1][2] = e;
    A[2][0] = c; A[2][1] = e; A[2][2] = f;

    CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

    const su2double factor = pow(outComplex/globalScale, 2./nDim) * pow(abs(EigVal[0]*EigVal[1]*EigVal[2]), -1./(2.*p+nDim));

    for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = min(max(abs(factor*EigVal[iDim]),eigmin),eigmax);

    CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

    //--- store lower triangle to be consistent with AMG
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 0, A[0][0]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 1, A[1][0]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 2, A[1][1]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 3, A[2][0]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 4, A[2][1]);
    solver_flow->GetNodes()->SetAnisoMetr(iPoint, 5, A[2][2]);

    //--- compute min, max, total complexity
    const su2double Vol = geometry->node[iPoint]->GetVolume();
    const su2double density = sqrt(abs(EigVal[0]*EigVal[1]*EigVal[2]));

    localMinDensity = min(localMinDensity, density);
    localMaxDensity = max(localMaxDensity, density);
    localTotComplex += density*Vol;
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&localMinDensity, &globalMinDensity, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&localMaxDensity, &globalMaxDensity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&localTotComplex, &globalTotComplex, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  globalMinDensity = localMinDensity;
  globalMaxDensity = localMaxDensity;
  globalTotComplex = localTotComplex;
#endif

  if(rank == MASTER_NODE) {
    cout << "Minimum density: " << globalMinDensity << "." << endl;
    cout << "Maximum density: " << globalMaxDensity << "." << endl;
    cout << "Mesh complexity: " << globalTotComplex << "." << endl;
  }

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    delete [] A[iDim];
    delete [] EigVec[iDim];
  }
  delete [] A;
  delete [] EigVec;
  delete [] EigVal;
}
