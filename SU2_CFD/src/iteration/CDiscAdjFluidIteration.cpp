/*!
 * \file CDiscAdjFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
                                        CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  StartTime = SU2_MPI::Wtime();

  const auto TimeIter = config[iZone]->GetTimeIter();
  const bool dual_time_1st = (config[iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool dual_time_2nd = (config[iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool dual_time = (dual_time_1st || dual_time_2nd);
  const bool grid_IsMoving = config[iZone]->GetGrid_Movement();
  const bool species = config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE;
  const bool heat = config[iZone]->GetWeakly_Coupled_Heat();
  const bool radiation = config[iZone]->AddRadiation();

  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometries = geometry[iZone][iInst];

  unsigned long nSolvers = 1;
  std::array<int, 5> solversToProcess{{FLOW_SOL}};
  if (turbulent) solversToProcess[nSolvers++] = TURB_SOL;
  if (species) solversToProcess[nSolvers++] = SPECIES_SOL;
  if (heat) solversToProcess[nSolvers++] = HEAT_SOL;
  if (radiation) solversToProcess[nSolvers++] = RAD_SOL;

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config[iZone]->GetTime_Marching() != TIME_MARCHING::STEADY) {
    const int Direct_Iter = static_cast<int>(config[iZone]->GetUnst_AdjointIter()) -
                            static_cast<int>(TimeIter) - 2 + dual_time;

    /*--- For dual-time stepping we want to load the already converged solution at previous timesteps.
     * In general we only load one file and shift the previously loaded solutions, on the first we
     * load one or two more (depending on dual time order). ---*/

    if (dual_time_2nd) {
      LoadUnsteady_Solution(geometry, solver, config, iZone, iInst, Direct_Iter - 2);
    } else if (dual_time_1st) {
      LoadUnsteady_Solution(geometry, solver, config, iZone, iInst, Direct_Iter - 1);
    }

    if (TimeIter == 0) {
      /*--- Push solution back one level. ---*/

      if (dual_time) {
        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          for (auto iSol = 0ul; iSol < nSolvers; ++iSol) {
            solver[iZone][iInst][iMesh][solversToProcess[iSol]]->GetNodes()->Set_Solution_time_n();
          }
          if (grid_IsMoving) {
            geometries[iMesh]->nodes->SetCoord_n();
          }
          if (config[iZone]->GetDynamic_Grid()) {
            geometries[iMesh]->nodes->SetVolume_n();
          }
        }
      }

      /*--- If required load another time step. Push the previous time step to n-1 and the
       loaded time step to n. ---*/

      if (dual_time_2nd) {
        LoadUnsteady_Solution(geometry, solver, config, iZone, iInst, Direct_Iter - 1);

        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          for (auto iSol = 0ul; iSol < nSolvers; ++iSol) {
            solver[iZone][iInst][iMesh][solversToProcess[iSol]]->GetNodes()->Set_Solution_time_n1();
            solver[iZone][iInst][iMesh][solversToProcess[iSol]]->GetNodes()->Set_Solution_time_n();
          }
          if (grid_IsMoving) {
            geometries[iMesh]->nodes->SetCoord_n1();
            geometries[iMesh]->nodes->SetCoord_n();
          }
          if (config[iZone]->GetDynamic_Grid()) {
            geometries[iMesh]->nodes->SetVolume_nM1();
            geometries[iMesh]->nodes->SetVolume_n();
          }
        }
      }

      /*--- Load solution timestep n ---*/

      LoadUnsteady_Solution(geometry, solver, config, iZone, iInst, Direct_Iter);

      if (config[iZone]->GetDeform_Mesh()) {
        solvers0[MESH_SOL]->LoadRestart(geometries, solver[iZone][iInst], config[iZone], Direct_Iter, true);
      }

    }
    if ((TimeIter > 0) && dual_time) {
      /*--- Here the primal solutions (only working variables) are loaded and put in the correct order
       into containers. For ALE the mesh coordinates have to be put into the correct containers as well,
       i.e. follow the same logic for the solution. Afterwards the GridVelocity is computed based on
       the Coordinates. ---*/

      /*--- Temporarily store the loaded volumes into old containers. ---*/
      if (config[iZone]->GetDynamic_Grid()) {
        for (auto iMesh=0; iMesh<=config[iZone]->GetnMGLevels();iMesh++) {
          geometries[iMesh]->nodes->SetVolume_Old();
          geometries[iMesh]->nodes->SetVolume_n_Old();
          geometries[iMesh]->nodes->SetVolume_nM1_Old();
        }
      }

      /*--- Load mesh solver. ---*/
      if (config[iZone]->GetDeform_Mesh()) {
        solvers0[MESH_SOL]->LoadRestart(geometries, solver[iZone][iInst], config[iZone], Direct_Iter, true);
      }

      /*--- Set volumes into correct containers ---*/
      if (config[iZone]->GetDynamic_Grid()) {
        for (auto iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          /*--- If negative iteration number, set default. ---*/
          if (Direct_Iter - 1 - dual_time_2nd < 0) {
            for(auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
              geometries[iMesh]->nodes->SetVolume(iPoint, 0.0);
            }
          }

          /*--- Set currently loaded volume to Volume_nM1 ---*/
          geometries[iMesh]->nodes->SetVolume_n();
          geometries[iMesh]->nodes->SetVolume_nM1();

          /*--- Set Volume_n and Volume from old containers ---*/
          geometries[iMesh]->nodes->SetVolume_n_from_OldnM1();
          geometries[iMesh]->nodes->SetVolume_from_Oldn();
        }
      }

      /*--- Temporarily store the loaded solution and coordinates in the "Old" arrays. ---*/

      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        for (auto iSol = 0ul; iSol < nSolvers; ++iSol) {
          solver[iZone][iInst][iMesh][solversToProcess[iSol]]->Set_OldSolution();
        }
        if (grid_IsMoving) {
          geometries[iMesh]->nodes->SetCoord_Old();
        }
      }

      /*--- Move timestep n to current solution. ---*/

      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        for (auto iSol = 0ul; iSol < nSolvers; ++iSol) {
          auto* s = solver[iZone][iInst][iMesh][solversToProcess[iSol]];
          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            s->GetNodes()->SetSolution(iPoint, s->GetNodes()->GetSolution_time_n(iPoint));
          }
        }
        if (grid_IsMoving) {
          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            geometries[iMesh]->nodes->SetCoord(iPoint, geometries[iMesh]->nodes->GetCoord_n(iPoint));
          }
        }
      }

      /*--- Finally, place the loaded solution in the correct place (n or n-1 depending on order). ---*/

      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        for (auto iSol = 0ul; iSol < nSolvers; ++iSol) {
          auto* s = solver[iZone][iInst][iMesh][solversToProcess[iSol]];
          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            if (dual_time_2nd) {
              /*--- If required also move timestep n-1 to timestep n. ---*/
              s->GetNodes()->Set_Solution_time_n(iPoint, s->GetNodes()->GetSolution_time_n1(iPoint));
              s->GetNodes()->Set_Solution_time_n1(iPoint, s->GetNodes()->GetSolution_Old(iPoint));
            } else {
              s->GetNodes()->Set_Solution_time_n(iPoint, s->GetNodes()->GetSolution_Old(iPoint));
            }
          }
        }
        if (grid_IsMoving) {
          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            if (dual_time_2nd) {
              geometries[iMesh]->nodes->SetCoord_n(iPoint, geometries[iMesh]->nodes->GetCoord_n1(iPoint));
              geometries[iMesh]->nodes->SetCoord_n1(iPoint, geometries[iMesh]->nodes->GetCoord_Old(iPoint));
            } else {
              geometries[iMesh]->nodes->SetCoord_n(iPoint, geometries[iMesh]->nodes->GetCoord_Old(iPoint));
            }
          }
        }
      }

    }  // else if TimeIter > 0

    /*--- Compute & set Grid Velocity via finite differences of the Coordinates. ---*/
    if (grid_IsMoving) {
      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++)
        geometries[iMesh]->SetGridVelocity(config[iZone]);
    }

  }  // if unsteady

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Store flow solution also in the adjoint solver in order to be able to reset it later ---*/

  if (TimeIter == 0 || dual_time) {
    auto SetSolutionDirect = [&](CSolver** solvers, int adj, int primal, unsigned long nPoint) {
      SU2_OMP_FOR_STAT(1024)
      for (auto iPoint = 0ul; iPoint < nPoint; iPoint++)
        solvers[adj]->GetNodes()->SetSolution_Direct(iPoint, solvers[primal]->GetNodes()->GetSolution(iPoint));
      END_SU2_OMP_FOR
    };

    for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
      auto solvers = solver[iZone][iInst][iMesh];
      SetSolutionDirect(solvers, ADJFLOW_SOL, FLOW_SOL, geometries[iMesh]->GetnPoint());
    }
    if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
      SetSolutionDirect(solvers0, ADJTURB_SOL, TURB_SOL, geometries[MESH_0]->GetnPoint());
    }
    if (species) {
      SetSolutionDirect(solvers0, ADJSPECIES_SOL, SPECIES_SOL, geometries[MESH_0]->GetnPoint());
    }
    if (heat) {
      SetSolutionDirect(solvers0, ADJHEAT_SOL, HEAT_SOL, geometries[MESH_0]->GetnPoint());
    }
    if (radiation) {
      SetSolutionDirect(solvers0, ADJRAD_SOL, RAD_SOL, geometries[MESH_0]->GetnPoint());
    }
  }

  solvers0[ADJFLOW_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
                                       MESH_0, 0, RUNTIME_ADJFLOW_SYS, false);

  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solvers0[ADJTURB_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
                                         MESH_0, 0, RUNTIME_ADJTURB_SYS, false);
  }
  if (species) {
    solvers0[ADJSPECIES_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
                                            MESH_0, 0, RUNTIME_ADJSPECIES_SYS, false);
  }
  if (heat) {
    solvers0[ADJHEAT_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
                                         MESH_0, 0, RUNTIME_ADJHEAT_SYS, false);
  }
  if (radiation) {
    solvers0[ADJRAD_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
                                        MESH_0, 0, RUNTIME_ADJRAD_SYS, false);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::LoadUnsteady_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                                   unsigned short iZone, unsigned short iInst, int DirectIter) {

  auto solvers = solver[iZone][iInst];
  auto geometries = geometry[iZone][iInst];
  const bool species = config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE;

  if (DirectIter >= 0) {
    if (rank == MASTER_NODE)
      cout << " Loading flow solution from direct iteration " << DirectIter << " for zone " << iZone << "." << endl;

    solvers[MESH_0][FLOW_SOL]->LoadRestart(geometries, solvers, config[iZone], DirectIter, true);

    if (turbulent) {
      solvers[MESH_0][TURB_SOL]->LoadRestart(geometries, solvers, config[iZone], DirectIter, false);
    }
    if (species) {
      solvers[MESH_0][SPECIES_SOL]->LoadRestart(geometries, solvers, config[iZone], DirectIter, false);
    }
    if (config[iZone]->GetWeakly_Coupled_Heat()) {
      solvers[MESH_0][HEAT_SOL]->LoadRestart(geometries, solvers, config[iZone], DirectIter, false);
    }
  } else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE)
      cout << " Setting freestream conditions at direct iteration " << DirectIter << " for zone " << iZone << "." << endl;

    for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
      solvers[iMesh][FLOW_SOL]->SetFreeStream_Solution(config[iZone]);
      solvers[iMesh][FLOW_SOL]->Preprocessing(geometries[iMesh], solvers[iMesh], config[iZone], iMesh,
                                              DirectIter, RUNTIME_FLOW_SYS, false);
      if (turbulent) {
        solvers[iMesh][TURB_SOL]->SetFreeStream_Solution(config[iZone]);
        solvers[iMesh][TURB_SOL]->Postprocessing(geometries[iMesh], solvers[iMesh], config[iZone], iMesh);
      }
      if (species) {
        solvers[iMesh][SPECIES_SOL]->SetFreeStream_Solution(config[iZone]);
        solvers[iMesh][SPECIES_SOL]->Postprocessing(geometries[iMesh], solvers[iMesh], config[iZone], iMesh);
      }
      if (config[iZone]->GetWeakly_Coupled_Heat()) {
        solvers[iMesh][HEAT_SOL]->SetFreeStream_Solution(config[iZone]);
        solvers[iMesh][HEAT_SOL]->Postprocessing(geometries[iMesh], solvers[iMesh], config[iZone], iMesh);
      }
    }
  }
}

void CDiscAdjFluidIteration::IterateDiscAdj(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                            unsigned short iZone, unsigned short iInst, bool CrossTerm) {
  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometry0 = geometry[iZone][iInst][MESH_0];

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solvers0[ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry0, config[iZone], CrossTerm);

    solvers0[ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry0, config[iZone]);
  }
  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solvers0[ADJTURB_SOL]->ExtractAdjoint_Solution(geometry0, config[iZone], CrossTerm);
  }
  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solvers0[ADJSPECIES_SOL]->ExtractAdjoint_Solution(geometry0, config[iZone], CrossTerm);
  }
  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solvers0[ADJHEAT_SOL]->ExtractAdjoint_Solution(geometry0, config[iZone], CrossTerm);
  }
  if (config[iZone]->AddRadiation()) {
    solvers0[ADJRAD_SOL]->ExtractAdjoint_Solution(geometry0, config[iZone], CrossTerm);

    solvers0[ADJRAD_SOL]->ExtractAdjoint_Variables(geometry0, config[iZone]);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                               unsigned short iZone, unsigned short iInst) {
  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometry0 = geometry[iZone][iInst][MESH_0];

  AD::ResizeAdjoints();

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Initialize the adjoints the conservative variables ---*/

  if (config[iZone]->GetFluidProblem()) {
    solvers0[ADJFLOW_SOL]->SetAdjoint_Output(geometry0, config[iZone]);
  }

  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solvers0[ADJTURB_SOL]->SetAdjoint_Output(geometry0, config[iZone]);
  }

  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solvers0[ADJSPECIES_SOL]->SetAdjoint_Output(geometry0, config[iZone]);
  }

  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solvers0[ADJHEAT_SOL]->SetAdjoint_Output(geometry0, config[iZone]);
  }

  if (config[iZone]->AddRadiation()) {
    solvers0[ADJRAD_SOL]->SetAdjoint_Output(geometry0, config[iZone]);
  }

  if (config[iZone]->GetFluidProblem() && !config[iZone]->GetMultizone_Problem()) {
    solvers0[FLOW_SOL]->SetVertexTractionsAdjoint(geometry0, config[iZone]);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                           unsigned short iZone, unsigned short iInst, RECORDING kind_recording) {
  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometry0 = geometry[iZone][iInst][MESH_0];

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  if (kind_recording == RECORDING::SOLUTION_VARIABLES || kind_recording == RECORDING::SOLUTION_AND_MESH) {
    /*--- Register flow and turbulent variables as input ---*/

    if (config[iZone]->GetFluidProblem()) {
      solvers0[ADJFLOW_SOL]->RegisterSolution(geometry0, config[iZone]);

      solvers0[ADJFLOW_SOL]->RegisterVariables(geometry0, config[iZone]);
    }

    if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
      solvers0[ADJTURB_SOL]->RegisterSolution(geometry0, config[iZone]);
    }
    if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
      solvers0[ADJSPECIES_SOL]->RegisterSolution(geometry0, config[iZone]);
    }
    if (config[iZone]->GetWeakly_Coupled_Heat()) {
      solvers0[ADJHEAT_SOL]->RegisterSolution(geometry0, config[iZone]);
    }
    if (config[iZone]->AddRadiation()) {
      solvers0[ADJRAD_SOL]->RegisterSolution(geometry0, config[iZone]);

      solvers0[ADJRAD_SOL]->RegisterVariables(geometry0, config[iZone]);
    }
  }

  if (kind_recording == RECORDING::MESH_COORDS || kind_recording == RECORDING::SOLUTION_AND_MESH) {
    /*--- Register node coordinates as input ---*/

    geometry0->RegisterCoordinates();
  }

  if (kind_recording == RECORDING::MESH_DEFORM) {
    /*--- Undeformed mesh coordinates ---*/
    solvers0[ADJMESH_SOL]->RegisterSolution(geometry0, config[iZone]);

    /*--- Boundary displacements ---*/
    solvers0[ADJMESH_SOL]->RegisterVariables(geometry0, config[iZone]);
  }
  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                             CConfig** config, unsigned short iZone, unsigned short iInst,
                                             RECORDING kind_recording) {
  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometry0 = geometry[iZone][iInst][MESH_0];

  if ((kind_recording == RECORDING::MESH_COORDS) ||
      (kind_recording == RECORDING::CLEAR_INDICES) ||
      (kind_recording == RECORDING::SOLUTION_AND_MESH)) {
    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/

    CGeometry::UpdateGeometry(geometry[iZone][iInst], config[iZone]);

    CGeometry::ComputeWallDistance(config, geometry);
  }

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Compute coupling between flow, turbulent and species equations ---*/
  solvers0[FLOW_SOL]->Preprocessing(geometry0, solvers0, config[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
  solvers0[FLOW_SOL]->InitiateComms(geometry0, config[iZone], SOLUTION);
  solvers0[FLOW_SOL]->CompleteComms(geometry0, config[iZone], SOLUTION);

  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solvers0[TURB_SOL]->Postprocessing(geometry0, solvers0,
                                                           config[iZone], MESH_0);
    solvers0[TURB_SOL]->InitiateComms(geometry0, config[iZone], SOLUTION);
    solvers0[TURB_SOL]->CompleteComms(geometry0, config[iZone], SOLUTION);
  }
  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solvers0[SPECIES_SOL]->Preprocessing(geometry0, solvers0, config[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
    solvers0[SPECIES_SOL]->InitiateComms(geometry0, config[iZone], SOLUTION);
    solvers0[SPECIES_SOL]->CompleteComms(geometry0, config[iZone], SOLUTION);
  }

  }
  END_SU2_OMP_PARALLEL

  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solvers0[HEAT_SOL]->Set_Heatflux_Areas(geometry0, config[iZone]);
    solvers0[HEAT_SOL]->Preprocessing(geometry0, solvers0, config[iZone], MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, true);
    solvers0[HEAT_SOL]->Postprocessing(geometry0, solvers0, config[iZone], MESH_0);
    solvers0[HEAT_SOL]->InitiateComms(geometry0, config[iZone], SOLUTION);
    solvers0[HEAT_SOL]->CompleteComms(geometry0, config[iZone], SOLUTION);
  }
  if (config[iZone]->AddRadiation()) {
    solvers0[RAD_SOL]->Postprocessing(geometry0, solvers0, config[iZone], MESH_0);
    solvers0[RAD_SOL]->InitiateComms(geometry0, config[iZone], SOLUTION);
    solvers0[RAD_SOL]->CompleteComms(geometry0, config[iZone], SOLUTION);
  }
}

void CDiscAdjFluidIteration::RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                            unsigned short iZone, unsigned short iInst) {
  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometry0 = geometry[iZone][iInst][MESH_0];

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Register conservative variables as output of the iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solvers0[ADJFLOW_SOL]->RegisterOutput(geometry0, config[iZone]);
  }
  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solvers0[ADJTURB_SOL]->RegisterOutput(geometry0, config[iZone]);
  }
  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solvers0[ADJSPECIES_SOL]->RegisterOutput(geometry0, config[iZone]);
  }
  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solvers0[ADJHEAT_SOL]->RegisterOutput(geometry0, config[iZone]);
  }
  if (config[iZone]->AddRadiation()) {
    solvers0[ADJRAD_SOL]->RegisterOutput(geometry0, config[iZone]);
  }
  if (config[iZone]->GetFluidProblem() && !config[iZone]->GetMultizone_Problem()) {
    solvers0[FLOW_SOL]->RegisterVertexTractions(geometry0, config[iZone]);
  }

  }
  END_SU2_OMP_PARALLEL
}

bool CDiscAdjFluidIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                     CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                     CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                     CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime - StartTime;

  /*--- Write the convergence history for the fluid (only screen output) ---*/

  output->SetHistoryOutput(geometry[iZone][INST_0][MESH_0], solver[iZone][INST_0][MESH_0], config[iZone],
                            config[iZone]->GetTimeIter(), config[iZone]->GetOuterIter(),
                            config[iZone]->GetInnerIter());

  return output->GetConvergence();
}
