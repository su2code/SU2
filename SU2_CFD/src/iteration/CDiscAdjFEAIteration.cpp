/*!
 * \file CDiscAdjFEAIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/iteration/CDiscAdjFEAIteration.hpp"
#include "../../include/iteration/CFEAIteration.hpp"
#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/output/COutput.hpp"

CDiscAdjFEAIteration::CDiscAdjFEAIteration(const CConfig *config) : CIteration(config), CurrentRecording(NONE) {
  fem_iteration = new CFEAIteration(config);

  // TEMPORARY output only for standalone structural problems
  if ((!config->GetFSI_Simulation()) && (rank == MASTER_NODE)) {
    bool de_effects = config->GetDE_Effects();
    unsigned short iVar;

    /*--- Header of the temporary output file ---*/
    ofstream myfile_res;
    myfile_res.open("Results_Reverse_Adjoint.txt");

    myfile_res << "Obj_Func"
               << " ";
    for (iVar = 0; iVar < config->GetnElasticityMod(); iVar++) myfile_res << "Sens_E_" << iVar << "\t";

    for (iVar = 0; iVar < config->GetnPoissonRatio(); iVar++) myfile_res << "Sens_Nu_" << iVar << "\t";

    if (config->GetTime_Domain()) {
      for (iVar = 0; iVar < config->GetnMaterialDensity(); iVar++) myfile_res << "Sens_Rho_" << iVar << "\t";
    }

    if (de_effects) {
      for (iVar = 0; iVar < config->GetnElectric_Field(); iVar++) myfile_res << "Sens_EField_" << iVar << "\t";
    }

    myfile_res << endl;

    myfile_res.close();
  }
}

CDiscAdjFEAIteration::~CDiscAdjFEAIteration(void) {}
void CDiscAdjFEAIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                      CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                      CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                      CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned long iPoint;
  unsigned short TimeIter = config[val_iZone]->GetTimeIter();
  bool dynamic = (config[val_iZone]->GetTime_Domain());

  int Direct_Iter;

  /*--- For the dynamic adjoint, load direct solutions from restart files. ---*/

  if (dynamic) {
    Direct_Iter = SU2_TYPE::Int(config[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(TimeIter) - 1;

    /*--- We want to load the already converged solution at timesteps n and n-1 ---*/

    /*--- Load solution at timestep n-1 ---*/

    LoadDynamic_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter - 1);

    /*--- Push solution back to correct array ---*/

    solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->Set_Solution_time_n();

    /*--- Push solution back to correct array ---*/

    solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->SetSolution_Accel_time_n();

    /*--- Push solution back to correct array ---*/

    solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->SetSolution_Vel_time_n();

    /*--- Load solution timestep n ---*/

    LoadDynamic_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter);

    /*--- Store FEA solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
      solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Direct(
          iPoint, solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->GetSolution(iPoint));
    }

    for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
      solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Accel_Direct(
          iPoint, solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->GetSolution_Accel(iPoint));
    }

    for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
      solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Vel_Direct(
          iPoint, solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->GetSolution_Vel(iPoint));
    }

  } else {
    /*--- Store FEA solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
      solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Direct(
          iPoint, solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->GetSolution(iPoint));
    }
  }

  solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->Preprocessing(
      geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, 0,
      RUNTIME_ADJFEA_SYS, false);
}

void CDiscAdjFEAIteration::LoadDynamic_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                                unsigned short val_iZone, unsigned short val_iInst,
                                                int val_DirectIter) {
  unsigned short iVar;
  unsigned long iPoint;
  bool update_geo = false;  // TODO: check

  if (val_DirectIter >= 0) {
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Loading FEA solution from direct iteration " << val_DirectIter << "." << endl;
    solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->LoadRestart(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], val_DirectIter, update_geo);
  } else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Setting static conditions at direct iteration " << val_DirectIter << "." << endl;
    /*--- Push solution back to correct array ---*/
    for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetnVar(); iVar++) {
        solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->SetSolution(iPoint, iVar, 0.0);
        solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->SetSolution_Accel(iPoint, iVar, 0.0);
        solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetNodes()->SetSolution_Vel(iPoint, iVar, 0.0);
      }
    }
  }
}

void CDiscAdjFEAIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                   CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                   CSurfaceMovement** surface_movement, CVolumetricMovement*** volume_grid_movement,
                                   CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  bool dynamic = (config[val_iZone]->GetTime_Domain());

  /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

  solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Solution(geometry[val_iZone][val_iInst][MESH_0],
                                                                            config[val_iZone]);

  solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Variables(geometry[val_iZone][val_iInst][MESH_0],
                                                                             config[val_iZone]);
  if (dynamic) {
    integration[val_iZone][val_iInst][ADJFEA_SOL]->SetConvergence(false);
  }
}

void CDiscAdjFEAIteration::SetRecording(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                        CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                        CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                        CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst,
                                        unsigned short kind_recording) {
  unsigned long InnerIter = config[ZONE_0]->GetInnerIter();
  unsigned long TimeIter = config[val_iZone]->GetTimeIter(), DirectTimeIter;
  bool dynamic = (config[val_iZone]->GetTime_Domain());

  DirectTimeIter = 0;
  if (dynamic) {
    DirectTimeIter = SU2_TYPE::Int(config[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(TimeIter) - 1;
  }

  /*--- Reset the tape ---*/

  AD::Reset();

  /*--- We only need to reset the indices if the current recording is different from the recording we want to have ---*/

  if (CurrentRecording != kind_recording && (CurrentRecording != NONE)) {
    solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->SetRecording(geometry[val_iZone][val_iInst][MESH_0],
                                                                   config[val_iZone]);

    /*--- Clear indices of coupling variables ---*/

    SetDependencies(solver, geometry, numerics, config, val_iZone, val_iInst, SOLUTION_AND_MESH);

    /*--- Run one iteration while tape is passive - this clears all indices ---*/

    fem_iteration->Iterate(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement,
                           FFDBox, val_iZone, val_iInst);
  }

  /*--- Prepare for recording ---*/

  solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->SetRecording(geometry[val_iZone][val_iInst][MESH_0],
                                                                 config[val_iZone]);

  /*--- Start the recording of all operations ---*/

  AD::StartRecording();

  /*--- Register FEA variables ---*/

  RegisterInput(solver, geometry, config, val_iZone, val_iInst, kind_recording);

  /*--- Compute coupling or update the geometry ---*/

  SetDependencies(solver, geometry, numerics, config, val_iZone, val_iInst, kind_recording);

  /*--- Set the correct direct iteration number ---*/

  if (dynamic) {
    config[val_iZone]->SetTimeIter(DirectTimeIter);
  }

  /*--- Run the direct iteration ---*/

  fem_iteration->Iterate(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement,
                         FFDBox, val_iZone, val_iInst);

  config[val_iZone]->SetTimeIter(TimeIter);

  /*--- Register structural variables and objective function as output ---*/

  RegisterOutput(solver, geometry, config, val_iZone, val_iInst);

  /*--- Stop the recording ---*/

  AD::StopRecording();

  /*--- Set the recording status ---*/

  CurrentRecording = kind_recording;

  /* --- Reset the number of the internal iterations---*/

  config[ZONE_0]->SetInnerIter(InnerIter);
}

void CDiscAdjFEAIteration::SetRecording(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                        unsigned short val_iZone, unsigned short val_iInst,
                                        unsigned short kind_recording) {
  /*--- Prepare for recording by resetting the solution to the initial converged solution ---*/

  solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->SetRecording(geometry[val_iZone][val_iInst][MESH_0],
                                                                 config[val_iZone]);
}

void CDiscAdjFEAIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                         unsigned short iZone, unsigned short iInst, unsigned short kind_recording) {
  if (kind_recording != MESH_COORDS) {
    /*--- Register structural displacements as input ---*/

    solver[iZone][iInst][MESH_0][ADJFEA_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Register variables as input ---*/

    solver[iZone][iInst][MESH_0][ADJFEA_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  } else {
    /*--- Register topology optimization densities (note direct solver) ---*/

    solver[iZone][iInst][MESH_0][FEA_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Register mesh coordinates for geometric sensitivities ---*/

    geometry[iZone][iInst][MESH_0]->RegisterCoordinates(config[iZone]);
  }
}

void CDiscAdjFEAIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                           CConfig** config, unsigned short iZone, unsigned short iInst,
                                           unsigned short kind_recording) {
  auto dir_solver = solver[iZone][iInst][MESH_0][FEA_SOL];
  auto adj_solver = solver[iZone][iInst][MESH_0][ADJFEA_SOL];
  auto structural_geometry = geometry[iZone][iInst][MESH_0];
  auto structural_numerics = numerics[iZone][iInst][MESH_0][FEA_SOL];

  /*--- Some numerics are only instanciated under these conditions ---*/
  bool fsi = config[iZone]->GetFSI_Simulation() || config[iZone]->GetMultizone_Problem();
  bool nonlinear = config[iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS;
  bool de_effects = config[iZone]->GetDE_Effects() && nonlinear;
  bool element_based = dir_solver->IsElementBased() && nonlinear;

  for (unsigned short iProp = 0; iProp < config[iZone]->GetnElasticityMod(); iProp++) {
    su2double E = adj_solver->GetVal_Young(iProp);
    su2double nu = adj_solver->GetVal_Poisson(iProp);
    su2double rho = adj_solver->GetVal_Rho(iProp);
    su2double rhoDL = adj_solver->GetVal_Rho_DL(iProp);

    /*--- Add dependencies for E and Nu ---*/

    structural_numerics[FEA_TERM]->SetMaterial_Properties(iProp, E, nu);

    /*--- Add dependencies for Rho and Rho_DL ---*/

    structural_numerics[FEA_TERM]->SetMaterial_Density(iProp, rho, rhoDL);

    /*--- Add dependencies for element-based simulations. ---*/

    if (element_based) {
      /*--- Neo Hookean Compressible ---*/
      structural_numerics[MAT_NHCOMP]->SetMaterial_Properties(iProp, E, nu);
      structural_numerics[MAT_NHCOMP]->SetMaterial_Density(iProp, rho, rhoDL);

      /*--- Ideal DE ---*/
      structural_numerics[MAT_IDEALDE]->SetMaterial_Properties(iProp, E, nu);
      structural_numerics[MAT_IDEALDE]->SetMaterial_Density(iProp, rho, rhoDL);

      /*--- Knowles ---*/
      structural_numerics[MAT_KNOWLES]->SetMaterial_Properties(iProp, E, nu);
      structural_numerics[MAT_KNOWLES]->SetMaterial_Density(iProp, rho, rhoDL);
    }
  }

  if (de_effects) {
    for (unsigned short iEField = 0; iEField < adj_solver->GetnEField(); iEField++) {
      structural_numerics[FEA_TERM]->Set_ElectricField(iEField, adj_solver->GetVal_EField(iEField));
      structural_numerics[DE_TERM]->Set_ElectricField(iEField, adj_solver->GetVal_EField(iEField));
    }
  }

  /*--- Add dependencies for element-based simulations. ---*/

  switch (config[iZone]->GetDV_FEA()) {
    case YOUNG_MODULUS:
    case POISSON_RATIO:
    case DENSITY_VAL:
    case DEAD_WEIGHT:
    case ELECTRIC_FIELD:

      for (unsigned short iDV = 0; iDV < adj_solver->GetnDVFEA(); iDV++) {
        su2double dvfea = adj_solver->GetVal_DVFEA(iDV);

        structural_numerics[FEA_TERM]->Set_DV_Val(iDV, dvfea);

        if (de_effects) structural_numerics[DE_TERM]->Set_DV_Val(iDV, dvfea);

        if (element_based) {
          structural_numerics[MAT_NHCOMP]->Set_DV_Val(iDV, dvfea);
          structural_numerics[MAT_IDEALDE]->Set_DV_Val(iDV, dvfea);
          structural_numerics[MAT_KNOWLES]->Set_DV_Val(iDV, dvfea);
        }
      }
      break;
  }

  /*--- MPI dependencies. ---*/

  dir_solver->InitiateComms(structural_geometry, config[iZone], SOLUTION_FEA);
  dir_solver->CompleteComms(structural_geometry, config[iZone], SOLUTION_FEA);

  if (kind_recording == MESH_COORDS) {
    structural_geometry->InitiateComms(structural_geometry, config[iZone], COORDINATES);
    structural_geometry->CompleteComms(structural_geometry, config[iZone], COORDINATES);
  }

  /*--- FSI specific dependencies. ---*/
  if (fsi) {
    /*--- Set relation between solution and predicted displacements, which are the transferred ones. ---*/
    dir_solver->PredictStruct_Displacement(structural_geometry, config[iZone]);
  }

  /*--- Topology optimization dependencies. ---*/

  /*--- We only differentiate wrt to this variable in the adjoint secondary recording. ---*/
  if (config[iZone]->GetTopology_Optimization() && (kind_recording == MESH_COORDS)) {
    /*--- The filter may require the volumes of the elements. ---*/
    structural_geometry->SetElemVolume();
    /// TODO: Ideally there would be a way to capture this dependency without the `static_cast`, but
    ///       making it a virtual method of CSolver does not feel "right" as its purpose could be confused.
    static_cast<CFEASolver*>(dir_solver)->FilterElementDensities(structural_geometry, config[iZone]);
  }
}

void CDiscAdjFEAIteration::RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                          unsigned short iZone, unsigned short iInst) {
  /*--- Register conservative variables as output of the iteration ---*/

  solver[iZone][iInst][MESH_0][ADJFEA_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
}

void CDiscAdjFEAIteration::InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                             unsigned short iZone, unsigned short iInst) {
  /*--- Initialize the adjoints the conservative variables ---*/

  solver[iZone][iInst][MESH_0][ADJFEA_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
}

void CDiscAdjFEAIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                  CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                  CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                  CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {}

bool CDiscAdjFEAIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                   CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                   CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                   CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  /*--- Write the convergence history (only screen output) ---*/

  output->SetHistory_Output(geometry[val_iZone][INST_0][MESH_0], solver[val_iZone][INST_0][MESH_0], config[val_iZone],
                            config[val_iZone]->GetTimeIter(), config[val_iZone]->GetOuterIter(),
                            config[val_iZone]->GetInnerIter());

  return output->GetConvergence();
}
void CDiscAdjFEAIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                       CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                       CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                       CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  bool dynamic = (config[val_iZone]->GetTime_Domain());

  // TEMPORARY output only for standalone structural problems
  if ((!config[val_iZone]->GetFSI_Simulation()) && (rank == MASTER_NODE)) {
    unsigned short iVar;

    bool de_effects = config[val_iZone]->GetDE_Effects();

    /*--- Header of the temporary output file ---*/
    ofstream myfile_res;
    myfile_res.open("Results_Reverse_Adjoint.txt", ios::app);

    myfile_res.precision(15);

    myfile_res << config[val_iZone]->GetTimeIter() << "\t";

    solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->Evaluate_ObjFunc(config[val_iZone]);
    myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][FEA_SOL]->GetTotal_ComboObj() << "\t";

    for (iVar = 0; iVar < config[val_iZone]->GetnElasticityMod(); iVar++)
      myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetTotal_Sens_E(iVar) << "\t";
    for (iVar = 0; iVar < config[val_iZone]->GetnPoissonRatio(); iVar++)
      myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetTotal_Sens_Nu(iVar) << "\t";
    if (dynamic) {
      for (iVar = 0; iVar < config[val_iZone]->GetnMaterialDensity(); iVar++)
        myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetTotal_Sens_Rho(iVar) << "\t";
    }
    if (de_effects) {
      for (iVar = 0; iVar < config[val_iZone]->GetnElectric_Field(); iVar++)
        myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetTotal_Sens_EField(iVar)
                   << "\t";
    }
    for (iVar = 0; iVar < solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetnDVFEA(); iVar++) {
      myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetTotal_Sens_DVFEA(iVar) << "\t";
    }

    myfile_res << endl;

    myfile_res.close();
  }

  // TEST: for implementation of python framework in standalone structural problems
  if ((!config[val_iZone]->GetFSI_Simulation()) && (rank == MASTER_NODE)) {
    /*--- Header of the temporary output file ---*/
    ofstream myfile_res;
    bool outputDVFEA = false;

    switch (config[val_iZone]->GetDV_FEA()) {
      case YOUNG_MODULUS:
        myfile_res.open("grad_young.opt");
        outputDVFEA = true;
        break;
      case POISSON_RATIO:
        myfile_res.open("grad_poisson.opt");
        outputDVFEA = true;
        break;
      case DENSITY_VAL:
      case DEAD_WEIGHT:
        myfile_res.open("grad_density.opt");
        outputDVFEA = true;
        break;
      case ELECTRIC_FIELD:
        myfile_res.open("grad_efield.opt");
        outputDVFEA = true;
        break;
      default:
        outputDVFEA = false;
        break;
    }

    if (outputDVFEA) {
      unsigned short iDV;
      unsigned short nDV = solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetnDVFEA();

      myfile_res << "INDEX"
                 << "\t"
                 << "GRAD" << endl;

      myfile_res.precision(15);

      for (iDV = 0; iDV < nDV; iDV++) {
        myfile_res << iDV;
        myfile_res << "\t";
        myfile_res << scientific << solver[val_iZone][val_iInst][MESH_0][ADJFEA_SOL]->GetTotal_Sens_DVFEA(iDV);
        myfile_res << endl;
      }

      myfile_res.close();
    }
  }
}
