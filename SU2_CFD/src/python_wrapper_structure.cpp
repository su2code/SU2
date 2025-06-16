/*!
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../include/drivers/CDriver.hpp"
#include "../include/drivers/CSinglezoneDriver.hpp"

void CDriver::PreprocessPythonInterface(CConfig** config, CGeometry**** geometry, CSolver***** solver) {
  int rank = MASTER_NODE;
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  /*--- Initialize boundary conditions customization, this is achieved through the Python wrapper. --- */
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config[iZone]->GetnMarker_PyCustom() > 0) {
      if (rank == MASTER_NODE) {
        cout << "----------------- Python Interface Preprocessing ( Zone " << iZone << " ) -----------------\n";
        cout << "Setting customized boundary conditions for zone " << iZone << endl;
      }
      for (iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        geometry[iZone][INST_0][iMesh]->SetCustomBoundary(config[iZone]);
      }
      geometry[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);

      if (solver[iZone][INST_0][MESH_0][FLOW_SOL]) {
        solver[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], solver[iZone][INST_0], config[iZone]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, etc.)   */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberTimeIter() const { return config_container[selected_zone]->GetnTime_Iter(); }

unsigned long CDriver::GetNumberInnerIter() const { return config_container[selected_zone]->GetnInner_Iter(); }
unsigned long CDriver::GetNumberOuterIter() const { return config_container[selected_zone]->GetnOuter_Iter(); }
void CDriver::SetNumberInnerIter(unsigned long val_iter) { config_container[selected_zone]->SetnInner_Iter(val_iter); }
void CDriver::SetNumberOuterIter(unsigned long val_iter) { config_container[selected_zone]->SetnOuter_Iter(val_iter); }

unsigned long CDriver::GetDensity_FreeStreamND() const {
  return SU2_TYPE::GetValue(config_container[selected_zone]->GetDensity_FreeStreamND());
  }
unsigned long CDriver::GetForce_Ref() const {
  return SU2_TYPE::GetValue(config_container[selected_zone]->GetForce_Ref());
  }

unsigned long CDriver::GetTimeIter() const { return TimeIter; }

passivedouble CDriver::GetUnsteadyTimeStep() const {
  return SU2_TYPE::GetValue(config_container[selected_zone]->GetTime_Step());
}

string CDriver::GetSurfaceFileName() const { return config_container[selected_zone]->GetSurfCoeff_FileName(); }

unsigned long CDriver::GetSolution(unsigned short iSOLVER, unsigned long iPoint, unsigned short iVar) {
  auto solver = solver_container[iZone][INST_0][MESH_0][iSOLVER];
  return SU2_TYPE::GetValue(solver->GetNodes()->GetSolution(iPoint,iVar));
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers                             */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetHeatSourcePosition(passivedouble alpha, passivedouble pos_x, passivedouble pos_y,
                                    passivedouble pos_z) {
  CSolver* solver = solver_container[selected_zone][INST_0][MESH_0][RAD_SOL];

  config_container[selected_zone]->SetHeatSource_Rot_Z(alpha);
  config_container[selected_zone]->SetHeatSource_Center(pos_x, pos_y, pos_z);

  solver->SetVolumetricHeatSource(geometry_container[selected_zone][INST_0][MESH_0], config_container[selected_zone]);
}

void CDriver::SetInletAngle(unsigned short iMarker, passivedouble alpha) {
  const su2double alpha_rad = alpha * PI_NUMBER / 180.0;

  const auto* geometry = geometry_container[selected_zone][INST_0][MESH_0];
  auto* flow_solver = solver_container[selected_zone][INST_0][MESH_0][FLOW_SOL];

  for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; ++iVertex) {
    flow_solver->SetInletFlowDir(iMarker, iVertex, 0, cos(alpha_rad));
    flow_solver->SetInletFlowDir(iMarker, iVertex, 1, sin(alpha_rad));
    if (geometry->GetnDim() == 3) flow_solver->SetInletFlowDir(iMarker, iVertex, 2, 0);
  }
}

void CDriver::SetFarFieldAoA(const passivedouble AoA) {
  config_container[selected_zone]->SetAoA(AoA);
  solver_container[selected_zone][INST_0][MESH_0][FLOW_SOL]->UpdateFarfieldVelocity(config_container[selected_zone]);
}

void CDriver::SetFarFieldAoS(const passivedouble AoS) {
  config_container[selected_zone]->SetAoS(AoS);
  solver_container[selected_zone][INST_0][MESH_0][FLOW_SOL]->UpdateFarfieldVelocity(config_container[selected_zone]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, etc.)   */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CSinglezoneDriver::SetInitialMesh() {
  DynamicMeshUpdate(0);

  SU2_OMP_PARALLEL {
    for (auto iMesh = 0u; iMesh <= main_config->GetnMGLevels(); iMesh++) {
      SU2_OMP_FOR_STAT(roundUpDiv(geometry_container[selected_zone][INST_0][iMesh]->GetnPoint(), omp_get_max_threads()))
      for (auto iPoint = 0ul; iPoint < geometry_container[selected_zone][INST_0][iMesh]->GetnPoint(); iPoint++) {
        /*--- Overwrite fictitious velocities. ---*/
        su2double Grid_Vel[3] = {0.0, 0.0, 0.0};

        /*--- Set the grid velocity for this coarse node. ---*/
        geometry_container[selected_zone][INST_0][iMesh]->nodes->SetGridVel(iPoint, Grid_Vel);
      }
      END_SU2_OMP_FOR
      /*--- Push back the volume. ---*/
      geometry_container[selected_zone][INST_0][iMesh]->nodes->SetVolume_n();
      geometry_container[selected_zone][INST_0][iMesh]->nodes->SetVolume_nM1();
    }
    /*--- Push back the solution so that there is no fictitious velocity at the next step. ---*/
    solver_container[selected_zone][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n();
    solver_container[selected_zone][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n1();
  }
  END_SU2_OMP_PARALLEL
}

void CDriver::BoundaryConditionsUpdate() {
  int rank = MASTER_NODE;

  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  if (rank == MASTER_NODE) cout << "Updating boundary conditions." << endl;
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
  }
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to dynamic mesh */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetTranslationRate(passivedouble xDot, passivedouble yDot, passivedouble zDot) {
  main_config->SetTranslation_Rate(0, xDot);
  main_config->SetTranslation_Rate(1, yDot);
  main_config->SetTranslation_Rate(2, zDot);
}

void CDriver::SetRotationRate(passivedouble rot_x, passivedouble rot_y, passivedouble rot_z) {
  main_config->SetRotation_Rate(0, rot_x);
  main_config->SetRotation_Rate(1, rot_y);
  main_config->SetRotation_Rate(2, rot_z);
}

void CDriver::SetMarkerRotationRate(unsigned short iMarker, passivedouble rot_x, passivedouble rot_y, passivedouble rot_z) {
  config_container[selected_zone]->SetMarkerRotationRate(iMarker, 0, rot_x);
  config_container[selected_zone]->SetMarkerRotationRate(iMarker, 1, rot_y);
  config_container[selected_zone]->SetMarkerRotationRate(iMarker, 2, rot_z);
}

void CDriver::SetMarkerTranslationRate(unsigned short iMarker, passivedouble vel_x, passivedouble vel_y, passivedouble vel_z) {
  config_container[selected_zone]->SetMarkerTranslationRate(iMarker, 0, vel_x);
  config_container[selected_zone]->SetMarkerTranslationRate(iMarker, 1, vel_y);
  config_container[selected_zone]->SetMarkerTranslationRate(iMarker, 2, vel_z);
}

