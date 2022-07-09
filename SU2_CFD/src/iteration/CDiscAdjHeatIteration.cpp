/*!
 * \file CDiscAdjHeatIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/iteration/CDiscAdjHeatIteration.hpp"
#include "../../include/output/COutput.hpp"

void CDiscAdjHeatIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                       CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                       CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                       CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  const auto TimeIter = config[val_iZone]->GetTimeIter();
  const bool dual_time_1st = (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool dual_time_2nd = (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool dual_time = (dual_time_1st || dual_time_2nd);

  auto solvers = solver[val_iZone][val_iInst];

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config[val_iZone]->GetTime_Marching() != TIME_MARCHING::STEADY) {
    const int Direct_Iter = static_cast<int>(config[val_iZone]->GetUnst_AdjointIter()) - static_cast<int>(TimeIter) - 2 + dual_time;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (TimeIter == 0) {
      if (dual_time_2nd) {
        /*--- Load solution at timestep n-2 ---*/

        LoadUnsteady_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter - 2);

        /*--- Push solution back to correct array ---*/

        for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          solvers[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n();
          solvers[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n1();
        }
      }
      if (dual_time) {
        /*--- Load solution at timestep n-1 ---*/

        LoadUnsteady_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter - 1);

        /*--- Push solution back to correct array ---*/

        for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          solvers[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n();
        }
      }

      /*--- Load solution timestep n ---*/

      LoadUnsteady_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter);
    }
    if ((TimeIter > 0) && dual_time) {
      /*--- Load solution timestep n - 2 ---*/

      LoadUnsteady_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter - 2);

      /*--- Temporarily store the loaded solution in the Solution_Old array ---*/

      for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
        solvers[iMesh][HEAT_SOL]->Set_OldSolution();

      /*--- Set Solution at timestep n to solution at n-1 ---*/

      for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
        for (auto iPoint = 0ul; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
          solvers[iMesh][HEAT_SOL]->GetNodes()->SetSolution(
              iPoint, solvers[iMesh][HEAT_SOL]->GetNodes()->GetSolution_time_n(iPoint));
        }
      }
      if (dual_time_1st) {
        /*--- Set Solution at timestep n-1 to the previously loaded solution ---*/
        for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          for (auto iPoint = 0ul; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
            solvers[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n(
                iPoint, solvers[iMesh][HEAT_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
          }
        }
      }
      if (dual_time_2nd) {
        /*--- Set Solution at timestep n-1 to solution at n-2 ---*/
        for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          for (auto iPoint = 0ul; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
            solvers[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n(
                iPoint, solvers[iMesh][HEAT_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
          }
        }
        /*--- Set Solution at timestep n-2 to the previously loaded solution ---*/
        for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          for (auto iPoint = 0ul; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
            solvers[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n1(
                iPoint, solvers[iMesh][HEAT_SOL]->GetNodes()->GetSolution_Old(iPoint));
          }
        }
      }
    }
  }

  /*--- Store flow solution also in the adjoint solver in order to be able to reset it later ---*/

  if (TimeIter == 0 || dual_time) {
    for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      for (auto iPoint = 0ul; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
        solvers[iMesh][ADJHEAT_SOL]->GetNodes()->SetSolution_Direct(
            iPoint, solvers[iMesh][HEAT_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
  }

  solvers[MESH_0][ADJHEAT_SOL]->Preprocessing(
      geometry[val_iZone][val_iInst][MESH_0], solvers[MESH_0], config[val_iZone], MESH_0, 0, RUNTIME_ADJHEAT_SYS, false);
}

void CDiscAdjHeatIteration::LoadUnsteady_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                                  unsigned short val_iZone, unsigned short val_iInst,
                                                  int val_DirectIter) {

  if (val_DirectIter >= 0) {
    if (rank == MASTER_NODE)
      cout << " Loading heat solution from direct iteration " << val_DirectIter << " for zone " << val_iZone << "." << endl;

    solver[val_iZone][val_iInst][MESH_0][HEAT_SOL]->LoadRestart(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], val_DirectIter, false);
  }
  else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE)
      cout << " Setting freestream conditions at direct iteration " << val_DirectIter << " for zone " << val_iZone << "." << endl;

    for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->SetFreeStream_Solution(config[val_iZone]);
      solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->Postprocessing(
          geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh], config[val_iZone], iMesh);
    }
  }
}

void CDiscAdjHeatIteration::IterateDiscAdj(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                           unsigned short val_iZone, unsigned short val_iInst, bool CrossTerm) {

  solver[val_iZone][val_iInst][MESH_0][ADJHEAT_SOL]->ExtractAdjoint_Solution(geometry[val_iZone][val_iInst][MESH_0],
                                                                             config[val_iZone], CrossTerm);
}

void CDiscAdjHeatIteration::InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                              unsigned short iZone, unsigned short iInst) {

  /*--- Initialize the adjoints the solution variables ---*/
  solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
}

void CDiscAdjHeatIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                          unsigned short iZone, unsigned short iInst, RECORDING kind_recording) {

  if (kind_recording == RECORDING::SOLUTION_VARIABLES || kind_recording == RECORDING::SOLUTION_AND_MESH) {
    /*--- Register flow and turbulent variables as input ---*/

    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  else if (kind_recording == RECORDING::MESH_COORDS) {
    /*--- Register node coordinates as input ---*/

    geometry[iZone][iInst][MESH_0]->RegisterCoordinates();
  }
  else if (kind_recording == RECORDING::MESH_DEFORM) {
    /*--- Register the variables of the mesh deformation ---*/
    /*--- Undeformed mesh coordinates ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Boundary displacements ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

void CDiscAdjHeatIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                            CConfig** config, unsigned short iZone, unsigned short iInst,
                                            RECORDING kind_recording) {

  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometries = geometry[iZone][iInst];

  if ((kind_recording == RECORDING::MESH_COORDS) ||
      (kind_recording == RECORDING::CLEAR_INDICES) ||
      (kind_recording == RECORDING::SOLUTION_AND_MESH)) {
    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/

    CGeometry::UpdateGeometry(geometries, config[iZone]);

    CGeometry::ComputeWallDistance(config, geometry);
  }

  solvers0[HEAT_SOL]->Set_Heatflux_Areas(geometries[MESH_0], config[iZone]);

  solvers0[HEAT_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone], MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, true);
  solvers0[HEAT_SOL]->Postprocessing(geometries[MESH_0], solvers0, config[iZone], MESH_0);

  solvers0[HEAT_SOL]->InitiateComms(geometries[MESH_0], config[iZone], SOLUTION);
  solvers0[HEAT_SOL]->CompleteComms(geometries[MESH_0], config[iZone], SOLUTION);
}

void CDiscAdjHeatIteration::RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                           unsigned short iZone, unsigned short iInst) {

  solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
}

void CDiscAdjHeatIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                   CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                   CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                   CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Dual time stepping strategy ---*/

  if ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
      (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {
    for (auto iMesh = 0u; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      integration[val_iZone][val_iInst][ADJHEAT_SOL]->SetConvergence(false);
    }
  }
}

bool CDiscAdjHeatIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                    CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                    CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                    CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  output->SetHistory_Output(geometry[val_iZone][INST_0][MESH_0], solver[val_iZone][INST_0][MESH_0], config[val_iZone],
                            config[val_iZone]->GetTimeIter(), config[val_iZone]->GetOuterIter(),
                            config[val_iZone]->GetInnerIter());

  return output->GetConvergence();
}
