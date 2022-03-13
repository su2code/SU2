/*!
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas, H. Patel, A. Gastaldi
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


#include "../include/drivers/CDriver.hpp"
#include "../include/drivers/CSinglezoneDriver.hpp"
#include "../../Common/include/toolboxes/geometry_toolbox.hpp"

void CDriver::PythonInterface_Preprocessing(CConfig **config, CGeometry ****geometry, CSolver *****solver){
    
    int rank = MASTER_NODE;
    SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
    
    /* --- Initialize boundary conditions customization, this is achieve through the Python wrapper --- */
    for(iZone=0; iZone < nZone; iZone++) {
        
        if (config[iZone]->GetnMarker_PyCustom() > 0) {
            
            if (rank == MASTER_NODE) cout << endl << "----------------- Python Interface Preprocessing ( Zone "<< iZone <<" ) -----------------" << endl;
            
            if (rank == MASTER_NODE) cout << "Setting customized boundary conditions for zone " << iZone << endl;
            for (iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
                geometry[iZone][INST_0][iMesh]->SetCustomBoundary(config[iZone]);
            }
            geometry[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
            
            if ((config[iZone]->GetKind_Solver() == MAIN_SOLVER::EULER) ||
                (config[iZone]->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES) ||
                (config[iZone]->GetKind_Solver() == MAIN_SOLVER::RANS)) {
                
                solver[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the global performance indices (Lift, Drag, ecc..) */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetDrag(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = config_container[ZONE_0]->GetFinestMesh();
    const auto CDrag = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();
    
    if (coefficient) {
        value = CDrag;
    } else {
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CDrag*scale_factor;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetLift(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = config_container[ZONE_0]->GetFinestMesh();
    const auto CLift = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();
    
    if (coefficient) {
        value = CLift;
    } else {
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CLift*scale_factor;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetRollMoment(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = config_container[ZONE_0]->GetFinestMesh();
    const auto CMx = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();
    
    if (coefficient) {
        value = CMx;
    } else {
        const auto ref_length   = config_container[ZONE_0]->GetRefLength();
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CMx*scale_factor*ref_length;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetPitchMoment(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = config_container[ZONE_0]->GetFinestMesh();
    const auto CMy = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();
    
    if (coefficient) {
        value = CMy;
    } else {
        const auto ref_length   = config_container[ZONE_0]->GetRefLength();
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CMy*scale_factor*ref_length;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetYawMoment(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = config_container[ZONE_0]->GetFinestMesh();
    const auto CMz = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();
    
    if (coefficient) {
        value = CMz;
    } else {
        const auto ref_length   = config_container[ZONE_0]->GetRefLength();
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CMz*scale_factor*ref_length;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetObjective() const {
    CConfig* config = config_container[ZONE_0];
    
    if (config->GetFluidProblem()) {
        CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
        solver->Evaluate_ObjFunc(config);
        
        return SU2_TYPE::GetValue(solver->GetTotal_ComboObj());
    } else {
        if (rank == MASTER_NODE) cout << "Objective function value not available... returning 0." << endl;
        
        return 0.0;
    }
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the far-field flow variables.                      */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetAngleOfAttack() const {
    return SU2_TYPE::GetValue(config_container[ZONE_0]->GetAoA());
}

void CDriver::SetAngleOfAttack(passivedouble value) {
    CConfig* config = config_container[ZONE_0];
    
    config->SetAoA(value);
    
    // Get the angle of attack to the free-stream velocity vector
    su2double velocity_inf_vec[nDim];
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        velocity_inf_vec[iDim] = config->GetVelocity_FreeStreamND()[iDim];
    }
    
    su2double velocity_inf_mag = GeometryToolbox::Norm(nDim, velocity_inf_vec);
    
    su2double alpha_rad = value            * PI_NUMBER/180.0;
    su2double beta_rad  = config->GetAoS() * PI_NUMBER/180.0;
    
    if (nDim == 3) {
        velocity_inf_vec[0] = cos(alpha_rad)*cos(beta_rad)*velocity_inf_mag;
        velocity_inf_vec[1] = sin(beta_rad)*velocity_inf_mag;
        velocity_inf_vec[2] = sin(alpha_rad)*cos(beta_rad)*velocity_inf_mag;
    } else if (nDim == 2) {
        velocity_inf_vec[0] = cos(alpha_rad)*velocity_inf_mag;
        velocity_inf_vec[1] = sin(alpha_rad)*velocity_inf_mag;
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        config->SetVelocity_FreeStreamND(velocity_inf_vec[iDim], iDim);
    }
    
    //  TODO:   Move this code to the solution update (would also be required for side-slip, Mach, and Reynolds)
}

passivedouble CDriver::GetAngleOfSideslip() const {
    return SU2_TYPE::GetValue(config_container[ZONE_0]->GetAoS());
}

void CDriver::SetAngleOfSideslip(passivedouble value) {
    config_container[ZONE_0]->SetAoS(value);
}

passivedouble CDriver::GetMachNumber() const {
    return SU2_TYPE::GetValue(config_container[ZONE_0]->GetMach());
}

void CDriver::SetMachNumber(passivedouble value) {
    config_container[ZONE_0]->SetMach(value);
}

passivedouble CDriver::GetReynoldsNumber() const {
    return SU2_TYPE::GetValue(config_container[ZONE_0]->GetReynolds());
}

void CDriver::SetReynoldsNumber(passivedouble value) {
    config_container[ZONE_0]->SetReynolds(value);
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the flow solver solution and variables.            */
/////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberStateVariables() const {
    return solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
}

unsigned long CDriver::GetNumberPrimitiveVariables() const {
    return solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnPrimVar();
}

vector<passivedouble> CDriver::GetResiduals() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver *solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    vector<passivedouble> values(nPoint*nVar, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->LinSysRes(iPoint, iVar);
            values[iPoint*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerResiduals(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver *solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    const auto nVar    = solver->GetnVar();
    vector<passivedouble> values(nVertex*nVar, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->LinSysRes(iPoint, iVar);
            values[iVertex*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetStates() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    vector<passivedouble> values(nPoint*nVar, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetNodes()->GetSolution(iPoint, iVar);
            values[iPoint*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerStates(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    const auto nVar    = solver->GetnVar();
    vector<passivedouble> values(nVertex*nVar, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetNodes()->GetSolution(iPoint, iVar);
            values[iVertex*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriver::SetStates(vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    if (values.size() != nPoint*nVar) {
        SU2_MPI::Error("Size does not match nPoint * nVar!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            solver->GetNodes()->SetSolution(iPoint, iVar, values[iPoint*nVar + iVar]);
        }
    }
}

void CDriver::SetMarkerStates(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    const auto nVar    = solver->GetnVar();
    if (values.size() != nVertex*nVar) {
        SU2_MPI::Error("Size does not match nVertex * nVar!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            solver->GetNodes()->SetSolution(iPoint, iVar, values[iVertex*nVar + iVar]);
        }
    }
}

vector<passivedouble> CDriver::GetPrimitiveStates() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nPrim  = 5;
    vector<passivedouble> values(nPoint*nPrim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        value = solver->GetNodes()->GetDensity(iPoint);
        values[iPoint*nPrim] = SU2_TYPE::GetValue(value);
        
        value = solver->GetNodes()->GetVelocity(iPoint, 0);
        values[iPoint*nPrim + 1] = SU2_TYPE::GetValue(value);
        value = solver->GetNodes()->GetVelocity(iPoint, 1);
        values[iPoint*nPrim + 2] = SU2_TYPE::GetValue(value);
        value = solver->GetNodes()->GetVelocity(iPoint, 2);
        values[iPoint*nPrim + 3] = SU2_TYPE::GetValue(value);
        
        value = solver->GetNodes()->GetPressure(iPoint);
        values[iPoint*nPrim + 4] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerPrimitiveStates(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    const auto nPrim   = 5;
    vector<passivedouble> values(nVertex*nPrim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        value = solver->GetNodes()->GetDensity(iPoint);
        values[iVertex*nPrim] = SU2_TYPE::GetValue(value);
        
        value = solver->GetNodes()->GetVelocity(iPoint, 0);
        values[iVertex*nPrim + 1] = SU2_TYPE::GetValue(value);
        value = solver->GetNodes()->GetVelocity(iPoint, 1);
        values[iVertex*nPrim + 2] = SU2_TYPE::GetValue(value);
        value = solver->GetNodes()->GetVelocity(iPoint, 2);
        values[iVertex*nPrim + 3] = SU2_TYPE::GetValue(value);
        
        value = solver->GetNodes()->GetPressure(iPoint);
        values[iVertex*nPrim + 4] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<passivedouble> CDriver::GetSpeedOfSound() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        value = solver->GetNodes()->GetSoundSpeed(iPoint);
        values[iPoint] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerSpeedOfSound(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        value = solver->GetNodes()->GetSoundSpeed(iPoint);
        values[iVertex] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerForces(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetVertexTractions(iMarker, iVertex, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the adjoint mesh solver solution.                  */
/////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetAdjointCoordinates() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
    
    if (!config->GetDeform_Mesh() || !config->GetDiscrete_Adjoint()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerAdjointCoordinates(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
    
    if (!config->GetDeform_Mesh() || !config->GetDiscrete_Adjoint()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriver::SetAdjointCoordinates(vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
    
    if (!config->GetDeform_Mesh() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    if (values.size() != nPoint*nDim) {
        SU2_MPI::Error("Size does not match nPoint * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetSolution(iPoint, iDim, values[iPoint*nDim + iDim]);
        }
    }
}

void CDriver::SetMarkerAdjointCoordinates(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
    
    if (!config->GetDeform_Mesh() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnPoint();
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetSolution(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the adjoint flow solver solution.                  */
/////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetAdjointStates() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    vector<passivedouble> values(nPoint*nVar, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetNodes()->GetSolution(iPoint, iVar);
            values[iPoint*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerAdjointStates(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    const auto nVar    = solver->GetnVar();
    vector<passivedouble> values(nVertex*nVar, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetNodes()->GetSolution(iPoint, iVar);
            values[iVertex*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriver::SetAdjointStates(vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    if (values.size() != nPoint*nVar) {
        SU2_MPI::Error("Size does not match nPoint * nVar!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            solver->GetNodes()->SetSolution(iPoint, iVar, values[iPoint*nVar + iVar]);
        }
    }
}

void CDriver::SetMarkerAdjointStates(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    const auto nVar    = solver->GetnVar();
    if (values.size() != nVertex*nVar) {
        SU2_MPI::Error("Size does not match nVertex * nVar!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            solver->GetNodes()->SetSolution(iPoint, iVar, values[iVar]);
        }
    }
}

vector<passivedouble> CDriver::GetMarkerAdjointForces(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetAdjointVertexTractions(iMarker, iVertex, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriver::SetMarkerAdjointForces(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->StoreVertexTractionsAdjoint(iMarker, iVertex, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

vector<passivedouble> CDriver::GetCoordinatesCoordinatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetProd_dCoordinates_dCoordinates(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerCoordinatesDisplacementsSensitivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetProd_dCoordinates_dDisplacements(iMarker, iVertex, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetObjectiveFarfieldVariablesSensitivity() const {
    CConfig* config = config_container[ZONE_0];
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const int nTrim = 2;
    vector<passivedouble> values(nTrim, 0.0);
    su2double value;
    
    value = solver->GetSens_dObjective_dVariables(0);
    values[0] = SU2_TYPE::GetValue(value);
    
    value = solver->GetSens_dObjective_dVariables(1);
    values[1] = SU2_TYPE::GetValue(value);
    
    return values;
}

vector<passivedouble> CDriver::GetResidualsFarfieldVariablesSensitivity() const {
    CConfig* config = config_container[ZONE_0];
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const int nTrim = 2;
    vector<passivedouble> values(nTrim, 0.0);
    su2double value;
    
    value = solver->GetProd_dResiduals_dVariables(0);
    values[0] = SU2_TYPE::GetValue(value);
    
    value = solver->GetProd_dResiduals_dVariables(1);
    values[1] = SU2_TYPE::GetValue(value);
    
    return values;
}

vector<passivedouble> CDriver::GetObjectiveStatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    vector<passivedouble> values(nPoint*nVar, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetSens_dObjective_dStates(iPoint, iVar);
            values[iPoint*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetResidualsStatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    vector<passivedouble> values(nPoint*nVar, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetProd_dResiduals_dStates(iPoint, iVar);
            values[iPoint*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetForcesStatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    const auto nVar   = solver->GetnVar();
    vector<passivedouble> values(nPoint*nVar, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iVar = 0u; iVar < nVar; iVar++) {
            value = solver->GetProd_dTractions_dStates(iPoint, iVar);
            values[iPoint*nVar + iVar] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetObjectiveCoordinatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetSens_dObjective_dCoordinates(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetResidualsCoordinatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetProd_dResiduals_dCoordinates(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetForcesCoordinatesSensitivity() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetProd_dTractions_dCoordinates(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerObjectiveDisplacementsSensitivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetSens_dObjective_dDisplacements(iMarker, iVertex, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerResidualsDisplacementsSensitivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetProd_dResiduals_dDisplacements(iMarker, iVertex, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerForcesDisplacementsSensitivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL];
    
    if (!config->GetFluidProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetProd_dTractions_dDisplacements(iMarker, iVertex, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, etc.).  */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberTimeIterations() const {
    return config_container[ZONE_0]->GetnTime_Iter();
}

unsigned long CDriver::GetCurrentTimeIteration() const{
    return TimeIter;
}

passivedouble CDriver::GetUnsteadyTimeStep() const {
    return SU2_TYPE::GetValue(config_container[ZONE_0]->GetTime_Step());
}

string CDriver::GetSurfaceFileName() const {
    return config_container[ZONE_0]->GetSurfCoeff_FileName();
}

///////////////////////////////////////////////////////////////////////////////
/* Functions related to CHT solver.                                          */
///////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetTemperatures() const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!config->GetFluidProblem()) {
        return {};
    }
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        value = solver->GetNodes()->GetTemperature(iPoint);
        values[iPoint] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerTemperatures(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
    if (!compressible) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        value = solver->GetNodes()->GetTemperature(iPoint);
        values[iVertex] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

void CDriver::SetMarkerTemperatures(unsigned short iMarker, vector<passivedouble> values) {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex  = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex) {
        SU2_MPI::Error("Size does not match nVertex!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        geometry->SetCustomBoundaryTemperature(iMarker, iVertex, values[iVertex]);
    }
}

vector<passivedouble> CDriver::GetMarkerHeatFlux(unsigned short iMarker, bool NormalVector) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
    if (!compressible) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    vector<passivedouble> values_normal(nVertex, 0.0);
    su2double value;
    
    const auto Prandtl_Lam  = config->GetPrandtl_Lam();
    const auto Gas_Constant = config->GetGas_ConstantND();
    const auto Gamma        = config->GetGamma();
    const auto Cp           = (Gamma/(Gamma - 1.0))*Gas_Constant;
    
    if (!NormalVector) {
        for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
            auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            
            auto laminar_viscosity    = solver->GetNodes()->GetLaminarViscosity(iPoint);
            auto thermal_conductivity = Cp*(laminar_viscosity/Prandtl_Lam);
            
            for (auto iDim = 0u; iDim < nDim; iDim++) {
                auto GradT = solver->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
                
                value = -thermal_conductivity*GradT;
                values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
            }
        }
        
        return values;
        
    } else {
        for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
            auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            
            auto laminar_viscosity    = solver->GetNodes()->GetLaminarViscosity(iPoint);
            auto thermal_conductivity = Cp*(laminar_viscosity/Prandtl_Lam);
            
            auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            auto Area   = GeometryToolbox::Norm(nDim, Normal);
            
            su2double dTdn = 0.0;
            
            for (auto iDim = 0u; iDim < nDim; iDim++) {
                auto GradT = solver->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
                dTdn += GradT*(Normal[iDim]/Area);
            }
            
            value = -thermal_conductivity*dTdn;
            values_normal[iVertex] = SU2_TYPE::GetValue(value);
        }
        
        return values_normal;
    }
}

void CDriver::SetMarkerNormalHeatFlux(unsigned short iMarker, vector<passivedouble> values) {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex) {
        SU2_MPI::Error("Size does not match nVertex!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        geometry->SetCustomBoundaryHeatFlux(iMarker, iVertex, values[iVertex]);
    }
}

vector<passivedouble> CDriver::GetMarkerThermalConductivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex, 0.0);
    su2double value;
    
    const auto Prandtl_Lam  = config->GetPrandtl_Lam();
    const auto Gas_Constant = config->GetGas_ConstantND();
    const auto Gamma        = config->GetGamma();
    const auto Cp           = (Gamma/(Gamma - 1.0))*Gas_Constant;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        auto laminar_viscosity = solver->GetNodes()->GetLaminarViscosity(iPoint);
        
        value = Cp*(laminar_viscosity/Prandtl_Lam);
        values[iVertex] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers.                            */
////////////////////////////////////////////////////////////////////////////////

vector<string> CDriver::GetFluidLoadMarkerTags() const {
    CConfig* config = config_container[ZONE_0];
    
    const auto nBoundariesMarker = config->GetnMarker_Fluid_Load();
    vector<string> interfaceBoundariesTagList;
    
    interfaceBoundariesTagList.resize(nBoundariesMarker);
    
    for (auto iMarker = 0u; iMarker < nBoundariesMarker; iMarker++) {
        auto Marker_Tag = config->GetMarker_Fluid_Load_TagBound(iMarker);
        interfaceBoundariesTagList[iMarker] = Marker_Tag;
    }
    
    return interfaceBoundariesTagList;
}

vector<string> CDriver::GetCHTMarkerTags() const {
    CConfig* config = config_container[ZONE_0];
    
    const auto nBoundariesMarker = config->GetnMarker_All();
    vector<string> CHTBoundariesTagList;
    
    for (auto iMarker = 0u; iMarker < nBoundariesMarker; iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) && config->GetMarker_All_PyCustom(iMarker)) {
            auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            CHTBoundariesTagList.push_back(Marker_Tag);
        }
    }
    
    return CHTBoundariesTagList;
}

vector<string> CDriver::GetInletMarkerTags() const {
    CConfig* config = config_container[ZONE_0];
    
    const auto nBoundariesMarker = config->GetnMarker_All();
    vector<string> BoundariesTagList;
    
    for (auto iMarker = 0u; iMarker < nBoundariesMarker; iMarker++) {
        bool isCustomizable = config->GetMarker_All_PyCustom(iMarker);
        bool isInlet = (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
        if (isCustomizable && isInlet) {
            auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            BoundariesTagList.push_back(Marker_Tag);
        }
    }
    
    return BoundariesTagList;
}

void CDriver::SetHeatSourcePosition(passivedouble alpha, passivedouble pos_x, passivedouble pos_y, passivedouble pos_z) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver *solver     = solver_container[ZONE_0][INST_0][MESH_0][RAD_SOL];
    
    config->SetHeatSource_Rot_Z(alpha);
    config->SetHeatSource_Center(pos_x, pos_y, pos_z);
    solver->SetVolumetricHeatSource(geometry, config);
}

void CDriver::SetInletAngle(unsigned short iMarker, passivedouble alpha) {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver *solver     = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    su2double alpha_rad = alpha * PI_NUMBER/180.0;
    
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        solver->SetInlet_FlowDir(iMarker, iVertex, 0, cos(alpha_rad));
        solver->SetInlet_FlowDir(iMarker, iVertex, 1, sin(alpha_rad));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, etc.).  */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CDriver::ResetConvergence() {
    CConfig* config = config_container[ZONE_0];
    
    for(auto iZone = 0u; iZone < nZone; iZone++) {
        switch (config->GetKind_Solver()) {
            case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
            case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
                integration_container[iZone][INST_0][FLOW_SOL]->SetConvergence(false);
                if (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) integration_container[iZone][INST_0][TURB_SOL]->SetConvergence(false);
                if(config_container[iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) integration_container[iZone][INST_0][TRANS_SOL]->SetConvergence(false);
                break;
                
            case MAIN_SOLVER::FEM_ELASTICITY:
                integration_container[iZone][INST_0][FEA_SOL]->SetConvergence(false);
                break;
                
            case MAIN_SOLVER::ADJ_EULER: case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
            case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
                integration_container[iZone][INST_0][ADJFLOW_SOL]->SetConvergence(false);
                if( (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) || (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) )
                    integration_container[iZone][INST_0][ADJTURB_SOL]->SetConvergence(false);
                break;
                
            default:
                break;
        }
    }
}

void CSinglezoneDriver::SetInitialMesh() {
    CConfig* config = config_container[ZONE_0];
    CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    DynamicMeshUpdate(0);
    
    SU2_OMP_PARALLEL {
        // Overwrite fictious velocities
        for (iMesh = 0u; iMesh <= config->GetnMGLevels(); iMesh++) {
            SU2_OMP_FOR_STAT(roundUpDiv(geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(),omp_get_max_threads()))
            for (auto iPoint = 0ul; iPoint < geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(); iPoint++) {
                
                /*--- Overwrite fictitious velocities ---*/
                su2double Grid_Vel[3] = {0.0, 0.0, 0.0};
                
                /*--- Set the grid velocity for this coarse node. ---*/
                geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetGridVel(iPoint, Grid_Vel);
            }
            END_SU2_OMP_FOR
            /*--- Push back the volume. ---*/
            geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetVolume_n();
            geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetVolume_nM1();
        }
        /*--- Push back the solution so that there is no fictious velocity at the next step. ---*/
        solver->GetNodes()->Set_Solution_time_n();
        solver->GetNodes()->Set_Solution_time_n1();
    }
    END_SU2_OMP_PARALLEL
}

void CDriver::UpdateBoundaryConditions() {
    int rank = MASTER_NODE;
    
    SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
    
    if (rank == MASTER_NODE) cout << "Updating boundary conditions." << endl;
    for (auto iZone = 0u; iZone < nZone; iZone++) {
        geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
    }
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to finite elements.                                      */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetMarkerFEAForces(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        su2double NodalForce[3] = {0.0, 0.0, 0.0};
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            NodalForce[iDim] = values[iVertex*nDim + iDim];
        }
        
        solver->GetNodes()->Set_FlowTraction(iPoint, NodalForce);
    }
}

vector<passivedouble> CDriver::GetMarkerFEADisplacements(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}


vector<passivedouble> CDriver::GetMarkerFEAVelocity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetDynamic_Analysis() != DYNAMIC) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution_Vel(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerCurrentFEAVelocity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetDynamic_Analysis() != DYNAMIC) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to adjoint simulations.                                  */
////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetMarkerDisplacementsSensitivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
    
    if (!config->GetDeform_Mesh() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Adjoint mesh solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBoundDisp_Sens(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerForcesSensitivity(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
    
    if (!config->GetStructuralProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriver::SetMarkerAdjointDisplacementSourceTerm(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
    
    if (!config->GetStructuralProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

void CDriver::SetMarkerAdjointVelocitySourceTerm(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
    
    if (!config->GetStructuralProblem() || !config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Adjoint structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}
