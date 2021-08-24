/*!
 * \file CDiscAdjFluidIteration.cpp
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

#include "../../include/iteration/CDiscAdjFluidIteration.hpp"
#include "../../include/output/COutput.hpp"

void CDiscAdjFluidIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                        CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                        CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                        CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  StartTime = SU2_MPI::Wtime();

  unsigned long iPoint;
  unsigned short TimeIter = config[val_iZone]->GetTimeIter();
  bool dual_time_1st = (config[val_iZone]->GetTime_Marching() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config[val_iZone]->GetTime_Marching() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  unsigned short iMesh;
  int Direct_Iter;
  bool heat = config[val_iZone]->GetWeakly_Coupled_Heat();
  bool grid_IsMoving = config[val_iZone]->GetGrid_Movement();

  //  /*--- Read the target pressure for inverse design. ---------------------------------------------*/
  //  if (config[val_iZone]->GetInvDesign_Cp() == YES)
  //    output->SetCp_InverseDesign(solver[val_iZone][val_iInst][MESH_0][FLOW_SOL],
  //    geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], ExtIter);

  //  /*--- Read the target heat flux ----------------------------------------------------------------*/
  //  if (config[ZONE_0]->GetInvDesign_HeatFlux() == YES)
  //    output->SetHeatFlux_InverseDesign(solver[val_iZone][val_iInst][MESH_0][FLOW_SOL],
  //    geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], ExtIter);

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config[val_iZone]->GetTime_Marching()) {
    Direct_Iter = SU2_TYPE::Int(config[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(TimeIter) - 2;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (dual_time) {
      Direct_Iter += 1;
    }

    if (TimeIter == 0) {
      if (dual_time_2nd) {
        /*--- Load solution at timestep n-2 ---*/
        LoadUnsteady_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter - 2);

        /*--- Push solution back to correct array ---*/

        for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
          solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1();
          if (turbulent) {
            solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
            solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1();
          }
          if (heat) {
            solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n();
            solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n1();
          }
          if (grid_IsMoving) {
            geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_n();
            geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_n1();
          }
        }
      }
      if (dual_time) {
        /*--- Load solution at timestep n-1 ---*/
        LoadUnsteady_Solution(geometry, solver, config, val_iZone, val_iInst, Direct_Iter - 1);

        /*--- Push solution back to correct array ---*/

        for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
          if (turbulent) {
            solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
          }
          if (heat) {
            solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n();
          }
          if (grid_IsMoving) {
            geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_n();
          }
        }
      }

      /*--- Load solution timestep n ---*/

      LoadUnsteady_Solution(geometry, solver, config, val_iInst, val_iZone, Direct_Iter);

      if (config[val_iZone]->GetDeform_Mesh()) {
        solver[val_iZone][val_iInst][MESH_0][MESH_SOL]->LoadRestart(
            geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], Direct_Iter, true);
      }

    } else if ((TimeIter > 0) && dual_time) {
      /*---
      Here the primal solutions (only working variables) are loaded and put in the correct order
      into containers. For ALE the mesh coordinates have to be put into the
      correct containers as well, i.e. follow the same logic for the solution.
      Afterwards the GridVelocity is computed based on the Coordinates.
      ---*/

      if (config[val_iZone]->GetDeform_Mesh()) {
        solver[val_iZone][val_iInst][MESH_0][MESH_SOL]->LoadRestart(
            geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], Direct_Iter, true);
      }

      /*--- Load solution timestep n-1 | n-2 for DualTimestepping 1st | 2nd order ---*/
      if (dual_time_1st) {
        LoadUnsteady_Solution(geometry, solver, config, val_iInst, val_iZone, Direct_Iter - 1);
      } else {
        LoadUnsteady_Solution(geometry, solver, config, val_iInst, val_iZone, Direct_Iter - 2);
      }

      /*--- Temporarily store the loaded solution in the Solution_Old array ---*/

      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
        solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->Set_OldSolution();
        if (turbulent) {
          solver[val_iZone][val_iInst][iMesh][TURB_SOL]->Set_OldSolution();
        }
        if (heat) {
          solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->Set_OldSolution();
        }
        if (grid_IsMoving) {
          geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_Old();
        }
      }

      /*--- Set Solution at timestep n to solution at n-1 ---*/

      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
          solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->SetSolution(
              iPoint, solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint));

          if (grid_IsMoving) {
            geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord(
                iPoint, geometry[val_iZone][val_iInst][iMesh]->nodes->GetCoord_n(iPoint));
          }
          if (turbulent) {
            solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->SetSolution(
                iPoint, solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          }
          if (heat) {
            solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->SetSolution(
                iPoint, solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          }
        }
      }
      if (dual_time_1st) {
        /*--- Set Solution at timestep n-1 to the previously loaded solution ---*/
        for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
            solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n(
                iPoint, solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint));

            if (grid_IsMoving) {
              geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_n(
                  iPoint, geometry[val_iZone][val_iInst][iMesh]->nodes->GetCoord_Old(iPoint));
            }
            if (turbulent) {
              solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
            if (heat) {
              solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
          }
        }
      }
      if (dual_time_2nd) {
        /*--- Set Solution at timestep n-1 to solution at n-2 ---*/
        for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
            solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n(
                iPoint, solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint));

            if (grid_IsMoving) {
              geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_n(
                  iPoint, geometry[val_iZone][val_iInst][iMesh]->nodes->GetCoord_n1(iPoint));
            }
            if (turbulent) {
              solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
            if (heat) {
              solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
          }
        }
        /*--- Set Solution at timestep n-2 to the previously loaded solution ---*/
        for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
          for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
            solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1(
                iPoint, solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint));

            if (grid_IsMoving) {
              geometry[val_iZone][val_iInst][iMesh]->nodes->SetCoord_n1(
                  iPoint, geometry[val_iZone][val_iInst][iMesh]->nodes->GetCoord_Old(iPoint));
            }
            if (turbulent) {
              solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1(
                  iPoint, solver[val_iZone][val_iInst][iMesh][TURB_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
            if (heat) {
              solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n1(
                  iPoint, solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
          }
        }
      }

    }  // else if Ext_Iter > 0

    /*--- Compute & set Grid Velocity via finite differences of the Coordinates. ---*/
    if (grid_IsMoving)
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
        geometry[val_iZone][val_iInst][iMesh]->SetGridVelocity(config[val_iZone], TimeIter);

  }  // if unsteady

  /*--- Store flow solution also in the adjoint solver in order to be able to reset it later ---*/

  if (TimeIter == 0 || dual_time) {
    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][iMesh]->GetnPoint(); iPoint++) {
        solver[val_iZone][val_iInst][iMesh][ADJFLOW_SOL]->GetNodes()->SetSolution_Direct(
            iPoint, solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
    if (turbulent && !config[val_iZone]->GetFrozen_Visc_Disc()) {
      for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
        solver[val_iZone][val_iInst][MESH_0][ADJTURB_SOL]->GetNodes()->SetSolution_Direct(
            iPoint, solver[val_iZone][val_iInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
    if (heat) {
      for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
        solver[val_iZone][val_iInst][MESH_0][ADJHEAT_SOL]->GetNodes()->SetSolution_Direct(
            iPoint, solver[val_iZone][val_iInst][MESH_0][HEAT_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
    if (config[val_iZone]->AddRadiation()) {
      for (iPoint = 0; iPoint < geometry[val_iZone][val_iInst][MESH_0]->GetnPoint(); iPoint++) {
        solver[val_iZone][val_iInst][MESH_0][ADJRAD_SOL]->GetNodes()->SetSolution_Direct(
            iPoint, solver[val_iZone][val_iInst][MESH_0][RAD_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
  }

  solver[val_iZone][val_iInst][MESH_0][ADJFLOW_SOL]->Preprocessing(
      geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, 0,
      RUNTIME_ADJFLOW_SYS, false);
  if (turbulent && !config[val_iZone]->GetFrozen_Visc_Disc()) {
    solver[val_iZone][val_iInst][MESH_0][ADJTURB_SOL]->Preprocessing(
        geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, 0,
        RUNTIME_ADJTURB_SYS, false);
  }
  if (heat) {
    solver[val_iZone][val_iInst][MESH_0][ADJHEAT_SOL]->Preprocessing(
        geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, 0,
        RUNTIME_ADJHEAT_SYS, false);
  }
  if (config[val_iZone]->AddRadiation()) {
    solver[val_iZone][val_iInst][MESH_0][ADJRAD_SOL]->Preprocessing(
        geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, 0,
        RUNTIME_ADJRAD_SYS, false);
  }
}

void CDiscAdjFluidIteration::LoadUnsteady_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                                   unsigned short val_iZone, unsigned short val_iInst,
                                                   int val_DirectIter) {
  unsigned short iMesh;
  bool heat = config[val_iZone]->GetWeakly_Coupled_Heat();

  if (val_DirectIter >= 0) {
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Loading flow solution from direct iteration " << val_DirectIter << "." << endl;
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->LoadRestart(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], val_DirectIter, true);
    if (turbulent) {
      solver[val_iZone][val_iInst][MESH_0][TURB_SOL]->LoadRestart(
          geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], val_DirectIter, false);
    }
    if (heat) {
      solver[val_iZone][val_iInst][MESH_0][HEAT_SOL]->LoadRestart(
          geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], val_DirectIter, false);
    }
  } else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Setting freestream conditions at direct iteration " << val_DirectIter << "." << endl;
    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->SetFreeStream_Solution(config[val_iZone]);
      solver[val_iZone][val_iInst][iMesh][FLOW_SOL]->Preprocessing(
          geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh], config[val_iZone], iMesh,
          val_DirectIter, RUNTIME_FLOW_SYS, false);
      if (turbulent) {
        solver[val_iZone][val_iInst][iMesh][TURB_SOL]->SetFreeStream_Solution(config[val_iZone]);
        solver[val_iZone][val_iInst][iMesh][TURB_SOL]->Postprocessing(
            geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh], config[val_iZone], iMesh);
      }
      if (heat) {
        solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->SetFreeStream_Solution(config[val_iZone]);
        solver[val_iZone][val_iInst][iMesh][HEAT_SOL]->Postprocessing(
            geometry[val_iZone][val_iInst][iMesh], solver[val_iZone][val_iInst][iMesh], config[val_iZone], iMesh);
      }
    }
  }
}

void CDiscAdjFluidIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                     CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                     CSurfaceMovement** surface_movement, CVolumetricMovement*** volume_grid_movement,
                                     CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone]);

    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (heat) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone]);

    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

void CDiscAdjFluidIteration::InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                               unsigned short iZone, unsigned short iInst) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  /*--- Initialize the adjoints the conservative variables ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (heat) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][FLOW_SOL]->SetVertexTractionsAdjoint(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

void CDiscAdjFluidIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                           unsigned short iZone, unsigned short iInst, unsigned short kind_recording) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  if (kind_recording == SOLUTION_VARIABLES || kind_recording == SOLUTION_AND_MESH) {
    /*--- Register flow and turbulent variables as input ---*/

    if (config[iZone]->GetFluidProblem()) {
      solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

      solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
    }

    if (turbulent && !frozen_visc) {
      solver[iZone][iInst][MESH_0][ADJTURB_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (heat) {
      solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (config[iZone]->AddRadiation()) {
      solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

      solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
  }

  if (kind_recording == MESH_COORDS) {
    /*--- Register node coordinates as input ---*/

    geometry[iZone][iInst][MESH_0]->RegisterCoordinates(config[iZone]);
  }

  /*--- Register the variables of the mesh deformation ---*/
  if (kind_recording == MESH_DEFORM) {
    /*--- Undeformed mesh coordinates ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Boundary displacements ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

void CDiscAdjFluidIteration::SetRecording(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                          unsigned short iZone, unsigned short iInst, unsigned short kind_recording) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();

  /*--- Prepare for recording by resetting the solution to the initial converged solution ---*/

  if (solver[iZone][iInst][MESH_0][ADJFEA_SOL]) {
    solver[iZone][iInst][MESH_0][ADJFEA_SOL]->SetRecording(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
    solver[iZone][iInst][iMesh][ADJFLOW_SOL]->SetRecording(geometry[iZone][iInst][iMesh], config[iZone]);
  }
  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->SetRecording(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->SetRecording(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][INST_0][MESH_0][ADJRAD_SOL]->SetRecording(geometry[iZone][INST_0][MESH_0], config[iZone]);
  }
}

void CDiscAdjFluidIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                             CConfig** config, unsigned short iZone, unsigned short iInst,
                                             unsigned short kind_recording) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  if ((kind_recording == MESH_COORDS) || (kind_recording == NONE) || (kind_recording == SOLUTION_AND_MESH)) {
    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/

    geometry[iZone][iInst][MESH_0]->UpdateGeometry(geometry[iZone][iInst], config[iZone]);

    CGeometry::ComputeWallDistance(config, geometry);
  }

  /*--- Compute coupling between flow and turbulent equations ---*/
  solver[iZone][iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                        config[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
  solver[iZone][iInst][MESH_0][FLOW_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  solver[iZone][iInst][MESH_0][FLOW_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);

  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                           config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][TURB_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][TURB_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }

  if (heat) {
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Set_Heatflux_Areas(geometry[iZone][iInst][MESH_0], config[iZone]);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                          config[iZone], MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, true);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                           config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][RAD_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                          config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][RAD_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][RAD_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }
}

void CDiscAdjFluidIteration::RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                            COutput* output, unsigned short iZone, unsigned short iInst) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  /*--- Register conservative variables as output of the iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (heat) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][FLOW_SOL]->RegisterVertexTractions(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

void CDiscAdjFluidIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                    CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                    CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                    CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned short iMesh;

  /*--- Dual time stepping strategy ---*/

  if ((config[val_iZone]->GetTime_Marching() == DT_STEPPING_1ST) ||
      (config[val_iZone]->GetTime_Marching() == DT_STEPPING_2ND)) {
    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      integration[val_iZone][val_iInst][ADJFLOW_SOL]->SetConvergence(false);
    }
  }
}

bool CDiscAdjFluidIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                     CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                     CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                     CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime - StartTime;

  /*--- Write the convergence history for the fluid (only screen output) ---*/

  output->SetHistory_Output(geometry[val_iZone][INST_0][MESH_0], solver[val_iZone][INST_0][MESH_0], config[val_iZone],
                            config[val_iZone]->GetTimeIter(), config[val_iZone]->GetOuterIter(),
                            config[val_iZone]->GetInnerIter());

  return output->GetConvergence();
}
void CDiscAdjFluidIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                         CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                         CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                         CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                                         unsigned short val_iInst) {}
