/*!
 * \file CDiscAdjFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.3.0 "Blackbird"
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
  //bool scalar = (config[iZone]->GetKind_Scalar_Model() != NO_SCALAR_MODEL);

  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometries = geometry[iZone][iInst];

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config[iZone]->GetTime_Marching() != TIME_MARCHING::STEADY) {
    const int Direct_Iter = static_cast<int>(config[iZone]->GetUnst_AdjointIter()) - static_cast<int>(TimeIter) - 2 + dual_time;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (TimeIter == 0) {
      if (dual_time_2nd) {
        /*--- Load solution at timestep n-2 ---*/
        LoadUnsteady_Solution(geometry, solver, config, iZone, iInst, Direct_Iter - 2);

        /*--- Push solution back to correct array ---*/

        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          auto solvers = solver[iZone][iInst][iMesh];

          solvers[FLOW_SOL]->GetNodes()->Set_Solution_time_n();
          solvers[FLOW_SOL]->GetNodes()->Set_Solution_time_n1();
          if (turbulent) {
            solvers[TURB_SOL]->GetNodes()->Set_Solution_time_n();
            solvers[TURB_SOL]->GetNodes()->Set_Solution_time_n1();
          }
          if (species) {
            solvers[SPECIES_SOL]->GetNodes()->Set_Solution_time_n();
            solvers[SPECIES_SOL]->GetNodes()->Set_Solution_time_n1();
          }
          if (heat) {
            solvers[HEAT_SOL]->GetNodes()->Set_Solution_time_n();
            solvers[HEAT_SOL]->GetNodes()->Set_Solution_time_n1();
          }
          if (grid_IsMoving) {
            geometries[iMesh]->nodes->SetCoord_n();
            geometries[iMesh]->nodes->SetCoord_n1();
          }
          if (config[iZone]->GetDynamic_Grid()) {
            geometries[iMesh]->nodes->SetVolume_n();
            geometries[iMesh]->nodes->SetVolume_nM1();
          }
        }
      }
      if (dual_time) {
        /*--- Load solution at timestep n-1 ---*/
        LoadUnsteady_Solution(geometry, solver, config, iZone, iInst, Direct_Iter - 1);

        /*--- Push solution back to correct array ---*/

        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          auto solvers = solver[iZone][iInst][iMesh];

          solvers[FLOW_SOL]->GetNodes()->Set_Solution_time_n();
          if (turbulent) {
            solvers[TURB_SOL]->GetNodes()->Set_Solution_time_n();
          }
          if (species) {
            solvers[SPECIES_SOL]->GetNodes()->Set_Solution_time_n();
          }
          if (heat) {
            solvers[HEAT_SOL]->GetNodes()->Set_Solution_time_n();
          }
          if (grid_IsMoving) {
            geometries[iMesh]->nodes->SetCoord_n();
          }
          if (config[iZone]->GetDynamic_Grid()) {
            geometries[iMesh]->nodes->SetVolume_n();
          }
        }
      }

      /*--- Load solution timestep n ---*/

      LoadUnsteady_Solution(geometry, solver, config, iInst, iZone, Direct_Iter);

      if (config[iZone]->GetDeform_Mesh()) {
        solvers0[MESH_SOL]->LoadRestart(geometries, solver[iZone][iInst], config[iZone], Direct_Iter, true);
      }

    }
    else if ((TimeIter > 0) && dual_time) {
      /*---
      Here the primal solutions (only working variables) are loaded and put in the correct order
      into containers. For ALE the mesh coordinates have to be put into the
      correct containers as well, i.e. follow the same logic for the solution.
      Afterwards the GridVelocity is computed based on the Coordinates.
      ---*/

      /*--- Temporarily store the loaded volumes into old containers ---*/
      if (config[iZone]->GetDynamic_Grid()) {
        for (auto iMesh=0; iMesh<=config[iZone]->GetnMGLevels();iMesh++) {
          geometries[iMesh]->nodes->SetVolume_Old();
          geometries[iMesh]->nodes->SetVolume_n_Old();
          geometries[iMesh]->nodes->SetVolume_nM1_Old();
        }
      }

      /*-- Load mesh solver ---*/
      if (config[iZone]->GetDeform_Mesh()) {
        solvers0[MESH_SOL]->LoadRestart(geometries, solver[iZone][iInst], config[iZone], Direct_Iter, true);
      }

      /*--- Load solution timestep n-1 | n-2 for DualTimestepping 1st | 2nd order ---*/
      if (dual_time_1st) {
        LoadUnsteady_Solution(geometry, solver, config, iInst, iZone, Direct_Iter - 1);
      } else {
        LoadUnsteady_Solution(geometry, solver, config, iInst, iZone, Direct_Iter - 2);

        /*--- Set volumes into correct containers ---*/
        if (config[iZone]->GetDynamic_Grid()) {
          for (auto iMesh=0; iMesh<=config[iZone]->GetnMGLevels();iMesh++) {
            /*--- If negative iteration number, set default ---*/
            if (Direct_Iter - 2 < 0) {
              for(auto iPoint=0ul; iPoint<geometries[iMesh]->GetnPoint();iPoint++) {
                geometries[iMesh]->nodes->SetVolume(iPoint,0.0);
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
      }

      /*--- Temporarily store the loaded solution in the Solution_Old array ---*/

      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        auto solvers = solver[iZone][iInst][iMesh];

        solvers[FLOW_SOL]->Set_OldSolution();
        if (turbulent) {
          solvers[TURB_SOL]->Set_OldSolution();
        }
        if (species) {
          solvers[SPECIES_SOL]->Set_OldSolution();
        }
        if (heat) {
          solvers[HEAT_SOL]->Set_OldSolution();
        }
        if (grid_IsMoving) {
          geometries[iMesh]->nodes->SetCoord_Old();
        }
      }

      /*--- Set Solution at timestep n to solution at n-1 ---*/

      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        auto solvers = solver[iZone][iInst][iMesh];

        for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
          solvers[FLOW_SOL]->GetNodes()->SetSolution(iPoint, solvers[FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint));

          if (grid_IsMoving) {
            geometries[iMesh]->nodes->SetCoord(iPoint, geometries[iMesh]->nodes->GetCoord_n(iPoint));
          }
          if (turbulent) {
            solvers[TURB_SOL]->GetNodes()->SetSolution(iPoint, solvers[TURB_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          }
          if (species) {
            solvers[SPECIES_SOL]->GetNodes()->SetSolution(iPoint, solvers[SPECIES_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          }
          if (heat) {
            solvers[HEAT_SOL]->GetNodes()->SetSolution(iPoint, solvers[HEAT_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          }
        }
      }
      if (dual_time_1st) {
        /*--- Set Solution at timestep n-1 to the previously loaded solution ---*/
        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          auto solvers = solver[iZone][iInst][iMesh];

          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            solvers[FLOW_SOL]->GetNodes()->Set_Solution_time_n(
                iPoint, solvers[FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint));

            if (grid_IsMoving) {
              geometries[iMesh]->nodes->SetCoord_n(iPoint, geometries[iMesh]->nodes->GetCoord_Old(iPoint));
            }
            if (turbulent) {
              solvers[TURB_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solvers[TURB_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
            if (species) {
              solvers[SPECIES_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solvers[SPECIES_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
            if (heat) {
              solvers[HEAT_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solvers[HEAT_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
          }
        }
      }
      if (dual_time_2nd) {
        /*--- Set Solution at timestep n-1 to solution at n-2 ---*/
        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          auto solvers = solver[iZone][iInst][iMesh];

          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            solvers[FLOW_SOL]->GetNodes()->Set_Solution_time_n(
                iPoint, solvers[FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint));

            if (grid_IsMoving) {
              geometries[iMesh]->nodes->SetCoord_n(iPoint, geometries[iMesh]->nodes->GetCoord_n1(iPoint));
            }
            if (turbulent) {
              solvers[TURB_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solvers[TURB_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
            if (species) {
              solvers[SPECIES_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solvers[SPECIES_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
            if (heat) {
              solvers[HEAT_SOL]->GetNodes()->Set_Solution_time_n(
                  iPoint, solvers[HEAT_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
          }
        }
        /*--- Set Solution at timestep n-2 to the previously loaded solution ---*/
        for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
          auto solvers = solver[iZone][iInst][iMesh];

          for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++) {
            solvers[FLOW_SOL]->GetNodes()->Set_Solution_time_n1(
                iPoint, solvers[FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint));

            if (grid_IsMoving) {
              geometries[iMesh]->nodes->SetCoord_n1(iPoint, geometries[iMesh]->nodes->GetCoord_Old(iPoint));
            }
            if (turbulent) {
              solvers[TURB_SOL]->GetNodes()->Set_Solution_time_n1(
                  iPoint, solvers[TURB_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
            if (species) {
              solvers[SPECIES_SOL]->GetNodes()->Set_Solution_time_n1(
                  iPoint, solvers[SPECIES_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
            if (heat) {
              solvers[HEAT_SOL]->GetNodes()->Set_Solution_time_n1(
                  iPoint, solvers[HEAT_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
          }
        }
      }

    }  // else if TimeIter > 0

    /*--- Compute & set Grid Velocity via finite differences of the Coordinates. ---*/
    if (grid_IsMoving)
      for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++)
        geometries[iMesh]->SetGridVelocity(config[iZone]);

  }  // if unsteady

  SU2_OMP_PARALLEL_(if(solvers0[ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Store flow solution also in the adjoint solver in order to be able to reset it later ---*/

  if (TimeIter == 0 || dual_time) {
    for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
      auto solvers = solver[iZone][iInst][iMesh];
      SU2_OMP_FOR_STAT(1024)
      for (auto iPoint = 0ul; iPoint < geometries[iMesh]->GetnPoint(); iPoint++)
        solvers[ADJFLOW_SOL]->GetNodes()->SetSolution_Direct(iPoint, solvers[FLOW_SOL]->GetNodes()->GetSolution(iPoint));
      END_SU2_OMP_FOR
    }
    if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
      SU2_OMP_FOR_STAT(1024)
      for (auto iPoint = 0ul; iPoint < geometries[MESH_0]->GetnPoint(); iPoint++)
        solvers0[ADJTURB_SOL]->GetNodes()->SetSolution_Direct(iPoint, solvers0[TURB_SOL]->GetNodes()->GetSolution(iPoint));
      END_SU2_OMP_FOR
    }
    if (species) {
      SU2_OMP_FOR_STAT(1024)
      for (auto iPoint = 0ul; iPoint < geometries[MESH_0]->GetnPoint(); iPoint++)
        solvers0[ADJSPECIES_SOL]->GetNodes()->SetSolution_Direct(iPoint, solvers0[SPECIES_SOL]->GetNodes()->GetSolution(iPoint));
      END_SU2_OMP_FOR
    }
    if (heat) {
      SU2_OMP_FOR_STAT(1024)
      for (auto iPoint = 0ul; iPoint < geometries[MESH_0]->GetnPoint(); iPoint++)
        solvers0[ADJHEAT_SOL]->GetNodes()->SetSolution_Direct(iPoint, solvers0[HEAT_SOL]->GetNodes()->GetSolution(iPoint));
      END_SU2_OMP_FOR
    }
    //if (scalar) {
    //  SU2_OMP_FOR_STAT(1024)
    //  for (auto iPoint = 0ul; iPoint < geometries[MESH_0]->GetnPoint(); iPoint++) {
    //    solvers0[ADJSCALAR_SOL]->GetNodes()->SetSolution_Direct(iPoint, solvers0[SCALAR_SOL]->GetNodes()->GetSolution(iPoint));
    //  END_SU2_OMP_FOR
    //  }
    }
    if (config[iZone]->AddRadiation()) {
      SU2_OMP_FOR_STAT(1024)
      for (auto iPoint = 0ul; iPoint < geometries[MESH_0]->GetnPoint(); iPoint++)
        solvers0[ADJRAD_SOL]->GetNodes()->SetSolution_Direct(iPoint, solvers0[RAD_SOL]->GetNodes()->GetSolution(iPoint));
      END_SU2_OMP_FOR
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
  //if (scalar) {
  //  solvers0[ADJSCALAR_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
  //                                         MESH_0, 0, RUNTIME_ADJSCALAR_SYS, false);
  }
  if (config[iZone]->AddRadiation()) {
    solvers0[ADJRAD_SOL]->Preprocessing(geometries[MESH_0], solvers0, config[iZone],
                                        MESH_0, 0, RUNTIME_ADJRAD_SYS, false);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::LoadUnsteady_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                                   unsigned short iZone, unsigned short iInst, int DirectIter) {

  auto solvers = solver[iZone][iInst];
  const bool species = config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE;

  if (DirectIter >= 0) {
    if (rank == MASTER_NODE)
      cout << " Loading flow solution from direct iteration " << DirectIter << " for zone " << iZone << "." << endl;

    solvers[MESH_0][FLOW_SOL]->LoadRestart(geometry[iZone][iInst], solvers, config[iZone], DirectIter, true);

    if (turbulent) {
      solvers[MESH_0][TURB_SOL]->LoadRestart(geometry[iZone][iInst], solvers, config[iZone], DirectIter, false);
    }
    if (species) {
      solvers[MESH_0][SPECIES_SOL]->LoadRestart(geometry[iZone][iInst], solvers, config[iZone], DirectIter, false);
    }
    if (config[iZone]->GetWeakly_Coupled_Heat()) {
      solvers[MESH_0][HEAT_SOL]->LoadRestart(geometry[iZone][iInst], solvers, config[iZone], DirectIter, false);
    }
  } else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE)
      cout << " Setting freestream conditions at direct iteration " << DirectIter << " for zone " << iZone << "." << endl;

    for (auto iMesh = 0u; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
      solvers[iMesh][FLOW_SOL]->SetFreeStream_Solution(config[iZone]);
      solvers[iMesh][FLOW_SOL]->Preprocessing(geometry[iZone][iInst][iMesh], solvers[iMesh], config[iZone], iMesh,
                                              DirectIter, RUNTIME_FLOW_SYS, false);
      if (turbulent) {
        solvers[iMesh][TURB_SOL]->SetFreeStream_Solution(config[iZone]);
        solvers[iMesh][TURB_SOL]->Postprocessing(geometry[iZone][iInst][iMesh], solvers[iMesh], config[iZone], iMesh);
      }
      if (species) {
        solvers[iMesh][SPECIES_SOL]->SetFreeStream_Solution(config[iZone]);
        solvers[iMesh][SPECIES_SOL]->Postprocessing(geometry[iZone][iInst][iMesh], solvers[iMesh], config[iZone], iMesh);
      }
      if (config[iZone]->GetWeakly_Coupled_Heat()) {
        solvers[iMesh][HEAT_SOL]->SetFreeStream_Solution(config[iZone]);
        solvers[iMesh][HEAT_SOL]->Postprocessing(geometry[iZone][iInst][iMesh], solvers[iMesh], config[iZone], iMesh);
      }
    }
  }
}

void CDiscAdjFluidIteration::IterateDiscAdj(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                            unsigned short iZone, unsigned short iInst, bool CrossTerm) {

  SU2_OMP_PARALLEL_(if(solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->GetHasHybridParallel())) {

  //bool scalar      = (config[iZone]->GetKind_Scalar_Model() != NO_SCALAR_MODEL);
  
  /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone], CrossTerm);

    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone], CrossTerm);
  }
  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solver[iZone][iInst][MESH_0][ADJSPECIES_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone], CrossTerm);
  }
  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone], CrossTerm);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone], CrossTerm);

    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  //if (scalar) {
   // solver[iZone][iInst][MESH_0][ADJSCALAR_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone], CrossTerm);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                               unsigned short iZone, unsigned short iInst) {

  SU2_OMP_PARALLEL_(if(solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->GetHasHybridParallel())) {

  //const bool scalar = (config[iZone]->GetKind_Scalar_Model() != NO_SCALAR_MODEL);

  /*--- Initialize the adjoints the conservative variables ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solver[iZone][iInst][MESH_0][ADJSPECIES_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  //if (scalar) {
  //  solver[iZone][iInst][MESH_0][ADJSCALAR_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->GetFluidProblem() && config[iZone]->GetSinglezone_Driver()) {
    solver[iZone][iInst][MESH_0][FLOW_SOL]->SetVertexTractionsAdjoint(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                           unsigned short iZone, unsigned short iInst, RECORDING kind_recording) {

  SU2_OMP_PARALLEL_(if(solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->GetHasHybridParallel())) {

  //const bool scalar      = (config[iZone]->GetKind_Scalar_Model() != NO_SCALAR_MODEL);

  if (kind_recording == RECORDING::SOLUTION_VARIABLES || kind_recording == RECORDING::SOLUTION_AND_MESH) {
    /*--- Register flow and turbulent variables as input ---*/

    if (config[iZone]->GetFluidProblem()) {
      solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

      solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
    }

    if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
      solver[iZone][iInst][MESH_0][ADJTURB_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
      solver[iZone][iInst][MESH_0][ADJSPECIES_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (config[iZone]->GetWeakly_Coupled_Heat()) {
      solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    //if (scalar) {
    //  solver[iZone][iInst][MESH_0][ADJSCALAR_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (config[iZone]->AddRadiation()) {
      solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

      solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
  }

  if (kind_recording == RECORDING::MESH_COORDS || kind_recording == RECORDING::SOLUTION_AND_MESH) {
    /*--- Register node coordinates as input ---*/

    geometry[iZone][iInst][MESH_0]->RegisterCoordinates();
  }

  if (kind_recording == RECORDING::MESH_DEFORM) {
    /*--- Undeformed mesh coordinates ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Boundary displacements ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                             CConfig** config, unsigned short iZone, unsigned short iInst,
                                             RECORDING kind_recording) {

  //const bool scalar      = (config[iZone]->GetKind_Scalar_Model() != NO_SCALAR_MODEL);

  if ((kind_recording == RECORDING::MESH_COORDS) ||
      (kind_recording == RECORDING::CLEAR_INDICES) ||
      (kind_recording == RECORDING::SOLUTION_AND_MESH)) {
    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/

    CGeometry::UpdateGeometry(geometry[iZone][iInst], config[iZone]);

    CGeometry::ComputeWallDistance(config, geometry);
  }

  SU2_OMP_PARALLEL_(if(solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->GetHasHybridParallel())) {

  /*--- Compute coupling between flow, turbulent and species equations ---*/
  solver[iZone][iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                        config[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
  solver[iZone][iInst][MESH_0][FLOW_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  solver[iZone][iInst][MESH_0][FLOW_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);

  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solver[iZone][iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                           config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][TURB_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][TURB_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }
  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solver[iZone][iInst][MESH_0][SPECIES_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                             config[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
    solver[iZone][iInst][MESH_0][SPECIES_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][SPECIES_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }

  }
  END_SU2_OMP_PARALLEL

  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Set_Heatflux_Areas(geometry[iZone][iInst][MESH_0], config[iZone]);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                          config[iZone], MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, true);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                           config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }

  //if (scalar) {
  //  solver[iZone][iInst][MESH_0][SCALAR_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
  //                                                          config[iZone], MESH_0, NO_RK_ITER, RUNTIME_SCALAR_SYS, true);
  //  solver[iZone][iInst][MESH_0][SCALAR_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
  //                                                           config[iZone], MESH_0);
  //  solver[iZone][iInst][MESH_0][SCALAR_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  //  solver[iZone][iInst][MESH_0][SCALAR_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  //}

  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][RAD_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                          config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][RAD_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][RAD_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }
}

void CDiscAdjFluidIteration::RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                            unsigned short iZone, unsigned short iInst) {

  SU2_OMP_PARALLEL_(if(solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->GetHasHybridParallel())) {

  //const bool scalar      = (config[iZone]->GetKind_Scalar_Model() != NO_SCALAR_MODEL);

  /*--- Register conservative variables as output of the iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (turbulent && !config[iZone]->GetFrozen_Visc_Disc()) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    solver[iZone][iInst][MESH_0][ADJSPECIES_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->GetWeakly_Coupled_Heat()) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  //}
  //if (scalar) {
  //  solver[iZone][iInst][MESH_0][ADJSCALAR_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->GetFluidProblem() && config[iZone]->GetSinglezone_Driver()) {
    solver[iZone][iInst][MESH_0][FLOW_SOL]->RegisterVertexTractions(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  }
  END_SU2_OMP_PARALLEL
}

void CDiscAdjFluidIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                    CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                    CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                    CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  /*--- Dual time stepping strategy ---*/

  if ((config[iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
      (config[iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {
    for (unsigned short iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
      integration[iZone][iInst][ADJFLOW_SOL]->SetConvergence(false);
    }
  }
}

bool CDiscAdjFluidIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                     CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                     CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                     CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime - StartTime;

  /*--- Write the convergence history for the fluid (only screen output) ---*/

  output->SetHistory_Output(geometry[iZone][INST_0][MESH_0], solver[iZone][INST_0][MESH_0], config[iZone],
                            config[iZone]->GetTimeIter(), config[iZone]->GetOuterIter(),
                            config[iZone]->GetInnerIter());

  return output->GetConvergence();
}
