/*!
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas
 * \version 8.3.0 "Harrier"
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
/* Functions related to farfield flow variables.                                */
//////////////////////////////////////////////////////////////////////////////////
passivedouble CDriver::GetAngleOfAttack() const { return SU2_TYPE::GetValue(main_config->GetAoA()); }

void CDriver::SetAngleOfAttack(const passivedouble AoA) {
  config_container[selected_zone]->SetAoA(AoA);
  solver_container[selected_zone][INST_0][MESH_0][FLOW_SOL]->UpdateFarfieldVelocity(config_container[selected_zone]);
}

passivedouble CDriver::GetAngleOfSideslip() const { return SU2_TYPE::GetValue(main_config->GetAoS()); }

void CDriver::SetAngleOfSideslip(const passivedouble AoS) {
  config_container[selected_zone]->SetAoS(AoS);
  solver_container[selected_zone][INST_0][MESH_0][FLOW_SOL]->UpdateFarfieldVelocity(config_container[selected_zone]);
}

passivedouble CDriver::GetMachNumber() const { return SU2_TYPE::GetValue(main_config->GetMach()); }

void CDriver::SetMachNumber(passivedouble value) {
  main_config->SetMach(value);
  UpdateFarfield();
}

passivedouble CDriver::GetReynoldsNumber() const { return SU2_TYPE::GetValue(main_config->GetReynolds()); }

void CDriver::SetReynoldsNumber(passivedouble value) {
  main_config->SetReynolds(value);
  UpdateFarfield();
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions related to the flow solver solution and variables.                 */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberStateVariables() const {
  if (!main_config->GetFluidProblem()) {
    SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
  }

  return solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
}

unsigned long CDriver::GetNumberPrimitiveVariables() const {
  if (!main_config->GetFluidProblem()) {
    SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
  }

  return solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnPrimVar();
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions related to the adjoint flow solver solution.                       */
//////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetMarkerAdjointForces(unsigned short iMarker, unsigned long iVertex) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetAdjointVertexTractions(iMarker, iVertex, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetCoordinatesCoordinatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dCoordinates_dCoordinates(iPoint, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetMarkerCoordinatesDisplacementsSensitivities(unsigned short iMarker, unsigned long iVertex) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dCoordinates_dDisplacements(iMarker, iVertex, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetObjectiveFarfieldVariablesSensitivities() const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }

  const int nTrim = 2;
  vector<passivedouble> values(nTrim, 0.0);

  su2double mach = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dVariables(0);
  values[0] = SU2_TYPE::GetValue(mach);

  su2double alpha = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dVariables(1);
  values[1] = SU2_TYPE::GetValue(alpha);

  return values;
}

vector<passivedouble> CDriver::GetResidualsFarfieldVariablesSensitivities() const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }

  const int nTrim = 2;
  vector<passivedouble> values(nTrim, 0.0);

  su2double mach = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dVariables(0);
  values[0] = SU2_TYPE::GetValue(mach);

  su2double alpha = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dVariables(1);
  values[1] = SU2_TYPE::GetValue(alpha);

  return values;
}

vector<passivedouble> CDriver::GetObjectiveStatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  const auto nVar = GetNumberStateVariables();
  vector<passivedouble> values(nVar, 0.0);

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dStates(iPoint, iVar);

    values[iVar] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetResidualsStatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  const auto nVar = GetNumberStateVariables();
  vector<passivedouble> values(nVar, 0.0);

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dStates(iPoint, iVar);

    values[iVar] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetForcesStatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  const auto nVar = GetNumberStateVariables();
  vector<passivedouble> values(nVar, 0.0);

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dTractions_dStates(iPoint, iVar);

    values[iVar] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetObjectiveCoordinatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dCoordinates(iPoint, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetResidualsCoordinatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dCoordinates(iPoint, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetForcesCoordinatesSensitivities(unsigned long iPoint) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dTractions_dCoordinates(iPoint, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetMarkerObjectiveDisplacementsSensitivities(unsigned short iMarker, unsigned long iVertex) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dDisplacements(iMarker, iVertex, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetMarkerResidualsDisplacementsSensitivities(unsigned short iMarker, unsigned long iVertex) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dDisplacements(iMarker, iVertex, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

vector<passivedouble> CDriver::GetMarkerForcesDisplacementsSensitivities(unsigned short iMarker, unsigned long iVertex) const {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dTractions_dDisplacements(iMarker, iVertex, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

void CDriver::SetAdjointSourceTerm(vector<passivedouble> values) {
  if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::RESIDUALS) {
    SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
  }

  const auto nPoint = GetNumberNodes();
  const auto nVar = GetNumberStateVariables();

  if (values.size() != nPoint * nVar) {
    SU2_MPI::Error("Size does not match nPoint * nVar!", CURRENT_FUNCTION);
  }

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->SetAdjoint_SourceTerm(iPoint, iVar, values[iPoint * nVar + iVar]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, etc.).  */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberTimeIter() const { return config_container[selected_zone]->GetnTime_Iter(); }

passivedouble CDriver::GetDensityFreeStreamND() const {
  return SU2_TYPE::GetValue(config_container[selected_zone]->GetDensity_FreeStreamND());
  }

passivedouble CDriver::GetForceRef() const {
  return SU2_TYPE::GetValue(config_container[selected_zone]->GetForce_Ref());
  }

unsigned long CDriver::GetTimeIter() const { return TimeIter; }

passivedouble CDriver::GetUnsteadyTimeStep() const {
  return SU2_TYPE::GetValue(config_container[selected_zone]->GetTime_Step());
}

string CDriver::GetSurfaceFileName() const { return config_container[selected_zone]->GetSurfCoeff_FileName(); }

//////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers.                              */
//////////////////////////////////////////////////////////////////////////////////

vector<string> CDriver::GetFluidLoadMarkerTags() const {
  const auto nMarker = main_config->GetnMarker_Fluid_Load();
  vector<string> tags;

  tags.resize(nMarker);

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    tags[iMarker] = main_config->GetMarker_Fluid_Load_TagBound(iMarker);
  }

  return tags;
}

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, etc.). */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void CDriver::UpdateGeometry() {
  geometry_container[ZONE_0][INST_0][MESH_0]->InitiateComms(main_geometry, main_config, MPI_QUANTITIES::COORDINATES);
  geometry_container[ZONE_0][INST_0][MESH_0]->CompleteComms(main_geometry, main_config, MPI_QUANTITIES::COORDINATES);

  geometry_container[ZONE_0][INST_0][MESH_0]->SetControlVolume(main_config, UPDATE);
  geometry_container[ZONE_0][INST_0][MESH_0]->SetBoundControlVolume(main_config, UPDATE);
  geometry_container[ZONE_0][INST_0][MESH_0]->SetMaxLength(main_config);
}

void CDriver::UpdateFarfield() {
  su2double Velocity_Ref = main_config->GetVelocity_Ref();
  su2double Alpha = main_config->GetAoA() * PI_NUMBER / 180.0;
  su2double Beta = main_config->GetAoS() * PI_NUMBER / 180.0;
  su2double Mach = main_config->GetMach();
  su2double Temperature = main_config->GetTemperature_FreeStream();
  su2double Gas_Constant = main_config->GetGas_Constant();
  su2double Gamma = main_config->GetGamma();
  su2double SoundSpeed = sqrt(Gamma * Gas_Constant * Temperature);

  if (nDim == 2) {
    main_config->GetVelocity_FreeStreamND()[0] = cos(Alpha) * Mach * SoundSpeed / Velocity_Ref;
    main_config->GetVelocity_FreeStreamND()[1] = sin(Alpha) * Mach * SoundSpeed / Velocity_Ref;
  }
  if (nDim == 3) {
    main_config->GetVelocity_FreeStreamND()[0] = cos(Alpha) * cos(Beta) * Mach * SoundSpeed / Velocity_Ref;
    main_config->GetVelocity_FreeStreamND()[1] = sin(Beta) * Mach * SoundSpeed / Velocity_Ref;
    main_config->GetVelocity_FreeStreamND()[2] = sin(Alpha) * Mach * SoundSpeed / Velocity_Ref;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions related to adjoint finite element simulations.                     */
//////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetMarkerForceSensitivities(unsigned short iMarker) const {
  if (!main_config->GetStructuralProblem() || !main_config->GetDiscrete_Adjoint()) {
    SU2_MPI::Error("Discrete adjoint structural solver is not defined!", CURRENT_FUNCTION);
  }
  if (main_config->GetKind_DiscreteAdjoint() != ENUM_DISC_ADJ_TYPE::FIXED_POINT) {
    SU2_MPI::Error("Discrete adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
  }

  const auto nVertex = GetNumberMarkerNodes(iMarker);
  vector<passivedouble> values(nVertex * nDim, 0.0);

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      values[iPoint * nDim + iDim] = SU2_TYPE::GetValue(
          solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL]->GetNodes()->GetFlowTractionSensitivity(iPoint, iDim));
    }
  }

  return values;
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions related to dynamic mesh.                                           */
//////////////////////////////////////////////////////////////////////////////////

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
