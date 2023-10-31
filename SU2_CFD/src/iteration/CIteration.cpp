/*!
 * \file iteration_structure.cpp
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

#include "../../include/iteration/CIteration.hpp"

#include "../../include/output/COutput.hpp"
#include "../../include/solvers/CFEASolver.hpp"

void CIteration::SetGrid_Movement(CGeometry** geometry, CSurfaceMovement* surface_movement,
                                  CVolumetricMovement* grid_movement, CSolver*** solver, CConfig* config,
                                  unsigned long IntIter, unsigned long TimeIter) {
  unsigned short Kind_Grid_Movement = config->GetKind_GridMovement();
  bool adjoint = config->GetContinuous_Adjoint();

  unsigned short val_iZone = config->GetiZone();

  /*--- Perform mesh movement depending on specified type ---*/
  switch (Kind_Grid_Movement) {
    case RIGID_MOTION:

      if (rank == MASTER_NODE) cout << endl << " Performing rigid mesh transformation." << endl;

      /*--- Move each node in the volume mesh using the specified type
       of rigid mesh motion. These routines also compute analytic grid
       velocities for the fine mesh. ---*/

      grid_movement->Rigid_Translation(geometry[MESH_0], config, val_iZone, TimeIter);
      grid_movement->Rigid_Plunging(geometry[MESH_0], config, val_iZone, TimeIter);
      grid_movement->Rigid_Pitching(geometry[MESH_0], config, val_iZone, TimeIter);
      grid_movement->Rigid_Rotation(geometry[MESH_0], config, val_iZone, TimeIter);

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement->UpdateMultiGrid(geometry, config);

      break;

    case STEADY_TRANSLATION:
      /*--- Set or update the translating frame mesh movement with the current translation rates,
       * which might be altered via the python interface. ---*/

      if (rank == MASTER_NODE) cout << "\n Setting translational grid velocities." << endl;

      /*--- Set the translational velocity on all grid levels. ---*/

      for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        geometry[iMGlevel]->SetTranslationalVelocity(config, iMGlevel == 0);

      break;

    case ROTATING_FRAME:
      /*--- Set or update the rotating frame mesh movement with the current translation and rotation
       * rates, which might be altered via the python interface. ---*/

      if (rank == MASTER_NODE) cout << "\n Setting rotating frame grid velocities." << endl;

      /*--- Set the grid velocities on all multigrid levels for a steadily
         rotating reference frame. ---*/

      for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++){
        geometry[iMGlevel]->SetRotationalVelocity(config, iMGlevel == 0);
        geometry[iMGlevel]->SetShroudVelocity(config);
      }

      break;
  }

  if (config->GetSurface_Movement(AEROELASTIC) || config->GetSurface_Movement(AEROELASTIC_RIGID_MOTION)) {
    /*--- Apply rigid mesh transformation to entire grid first, if necessary ---*/
    if (IntIter == 0) {
      if (Kind_Grid_Movement == AEROELASTIC_RIGID_MOTION) {
        if (rank == MASTER_NODE) cout << endl << " Performing rigid mesh transformation." << endl;

        /*--- Move each node in the volume mesh using the specified type
         of rigid mesh motion. These routines also compute analytic grid
         velocities for the fine mesh. ---*/

        grid_movement->Rigid_Translation(geometry[MESH_0], config, val_iZone, TimeIter);
        grid_movement->Rigid_Plunging(geometry[MESH_0], config, val_iZone, TimeIter);
        grid_movement->Rigid_Pitching(geometry[MESH_0], config, val_iZone, TimeIter);
        grid_movement->Rigid_Rotation(geometry[MESH_0], config, val_iZone, TimeIter);

        /*--- Update the multigrid structure after moving the finest grid,
         including computing the grid velocities on the coarser levels. ---*/

        grid_movement->UpdateMultiGrid(geometry, config);
      }

    }

    /*--- Use the if statement to move the grid only at selected dual time step iterations. ---*/
    else if (IntIter % config->GetAeroelasticIter() == 0) {
      if (rank == MASTER_NODE) cout << endl << " Solving aeroelastic equations and updating surface positions." << endl;

      /*--- Solve the aeroelastic equations for the new node locations of the moving markers(surfaces) ---*/

      solver[MESH_0][FLOW_SOL]->Aeroelastic(surface_movement, geometry[MESH_0], config, TimeIter);

      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE) cout << " Deforming the volume grid due to the aeroelastic movement." << endl;
      grid_movement->SetVolume_Deformation(geometry[MESH_0], config, true);

      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/

      if (rank == MASTER_NODE) cout << " Computing grid velocities by finite differencing." << endl;
      geometry[MESH_0]->SetGridVelocity(config);

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement->UpdateMultiGrid(geometry, config);
    }
  }

  if (config->GetSurface_Movement(EXTERNAL) || config->GetSurface_Movement(EXTERNAL_ROTATION)) {
    /*--- Apply rigid rotation to entire grid first, if necessary ---*/

    if (Kind_Grid_Movement == EXTERNAL_ROTATION) {
      if (rank == MASTER_NODE) cout << " Updating node locations by rigid rotation." << endl;
      grid_movement->Rigid_Rotation(geometry[MESH_0], config, val_iZone, TimeIter);
    }

    /*--- Load new surface node locations from external files ---*/

    if (rank == MASTER_NODE) cout << " Updating surface locations from file." << endl;
    surface_movement->SetExternal_Deformation(geometry[MESH_0], config, val_iZone, TimeIter);

    /*--- Deform the volume grid around the new boundary locations ---*/

    if (rank == MASTER_NODE) cout << " Deforming the volume grid." << endl;
    grid_movement->SetVolume_Deformation(geometry[MESH_0], config, true);

    /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/

    if (!adjoint) {
      if (rank == MASTER_NODE) cout << " Computing grid velocities by finite differencing." << endl;
      geometry[MESH_0]->SetGridVelocity(config);
    }

    /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

    grid_movement->UpdateMultiGrid(geometry, config);
  }
}

void CIteration::SetMesh_Deformation(CGeometry** geometry, CSolver** solver, CNumerics*** numerics, CConfig* config,
                                     RECORDING kind_recording) {
  if (!config->GetDeform_Mesh()) return;

  /*--- Perform the elasticity mesh movement ---*/

  bool wasActive = false;
  if ((kind_recording != RECORDING::MESH_DEFORM) && !config->GetMultizone_Problem()) {
    /*--- In a primal run, AD::TapeActive returns a false ---*/
    /*--- In any other recordings, the tape is passive during the deformation. ---*/
    wasActive = AD::BeginPassive();
  }

  /*--- Set the stiffness of each element mesh into the mesh numerics ---*/

  solver[MESH_SOL]->SetMesh_Stiffness(numerics[MESH_SOL], config);

  /*--- Deform the volume grid around the new boundary locations ---*/

  solver[MESH_SOL]->DeformMesh(geometry, numerics[MESH_SOL], config);

  /*--- Continue recording. ---*/
  AD::EndPassive(wasActive);
}

void CIteration::Output(COutput* output, CGeometry**** geometry, CSolver***** solver, CConfig** config,
                        unsigned long InnerIter, bool StopCalc, unsigned short val_iZone, unsigned short val_iInst) {
  output->SetResultFiles(geometry[val_iZone][INST_0][MESH_0], config[val_iZone], solver[val_iZone][INST_0][MESH_0],
                          InnerIter);
}
