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
    
    const auto FinestMesh = main_config->GetFinestMesh();
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
    
    const auto FinestMesh = main_config->GetFinestMesh();
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
    
    const auto FinestMesh = main_config->GetFinestMesh();
    const auto CMx = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();
    
    if (coefficient) {
        value = CMx;
    } else {
        const auto ref_length   = main_config->GetRefLength();
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CMx*scale_factor*ref_length;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetPitchMoment(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = main_config->GetFinestMesh();
    const auto CMy = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();
    
    if (coefficient) {
        value = CMy;
    } else {
        const auto ref_length   = main_config->GetRefLength();
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CMy*scale_factor*ref_length;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetYawMoment(bool coefficient) const {
    su2double value;
    
    const auto FinestMesh = main_config->GetFinestMesh();
    const auto CMz = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();
    
    if (coefficient) {
        value = CMz;
    } else {
        const auto ref_length   = main_config->GetRefLength();
        const auto scale_factor = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
        value = CMz*scale_factor*ref_length;
    }
    
    return SU2_TYPE::GetValue(value);
}

passivedouble CDriver::GetObjective() const {
    if (!main_config->GetFluidProblem()) {
        if (rank == MASTER_NODE) cout << "Objective function value not available... returning 0." << endl;
        
        return 0.0;
    }

    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    solver->Evaluate_ObjFunc(main_config);
    
    return SU2_TYPE::GetValue(solver->GetTotal_ComboObj());
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to design variables.                                  */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetNumberDesignVariables() const {
    return SU2_TYPE::GetValue(main_config->GetnDV());
}

passivedouble CDriver::GetNumberFFDBoxes() const {
    return SU2_TYPE::GetValue(main_config->GetnFFDBox());
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the far-field flow variables.                      */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetAngleOfAttack() const {
    return SU2_TYPE::GetValue(main_config->GetAoA());
}

void CDriver::SetAngleOfAttack(passivedouble value) {
    main_config->SetAoA(value);
    
    // Get the angle of attack to the free-stream velocity vector
    su2double velocity_inf_vec[nDim];
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        velocity_inf_vec[iDim] = main_config->GetVelocity_FreeStreamND()[iDim];
    }
    
    su2double velocity_inf_mag = GeometryToolbox::Norm(nDim, velocity_inf_vec);
    
    su2double alpha_rad = value                * PI_NUMBER/180.0;
    su2double beta_rad  = main_config->GetAoS() * PI_NUMBER/180.0;
    
    if (nDim == 3) {
        velocity_inf_vec[0] = cos(alpha_rad)*cos(beta_rad)*velocity_inf_mag;
        velocity_inf_vec[1] = sin(beta_rad)*velocity_inf_mag;
        velocity_inf_vec[2] = sin(alpha_rad)*cos(beta_rad)*velocity_inf_mag;
    } else if (nDim == 2) {
        velocity_inf_vec[0] = cos(alpha_rad)*velocity_inf_mag;
        velocity_inf_vec[1] = sin(alpha_rad)*velocity_inf_mag;
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        main_config->SetVelocity_FreeStreamND(velocity_inf_vec[iDim], iDim);
    }
    
    //  TODO:   Move this code to the solution update (would also be required for side-slip, Mach, and Reynolds)
}

passivedouble CDriver::GetAngleOfSideslip() const {
    return SU2_TYPE::GetValue(main_config->GetAoS());
}

void CDriver::SetAngleOfSideslip(passivedouble value) {
    main_config->SetAoS(value);
}

passivedouble CDriver::GetMachNumber() const {
    return SU2_TYPE::GetValue(main_config->GetMach());
}

void CDriver::SetMachNumber(passivedouble value) {
    main_config->SetMach(value);
}

passivedouble CDriver::GetReynoldsNumber() const {
    return SU2_TYPE::GetValue(main_config->GetReynolds());
}

void CDriver::SetReynoldsNumber(passivedouble value) {
    main_config->SetReynolds(value);
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the flow solver solution and variables.            */
/////////////////////////////////////////////////////////////////////////////

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

vector<vector<passivedouble>> CDriver::GetResiduals() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetResiduals(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetResiduals(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    const auto nVar = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->LinSysRes(iPoint, iVar);

        values[iVar] = SU2_TYPE::GetValue(value);
    } 
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerResiduals(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerResiduals(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerResiduals(unsigned short iMarker, unsigned long iVertex) const {
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    const auto nVar = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->LinSysRes(iPoint, iVar);

        values[iVar] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetStates() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetStates(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetStates(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    const auto nVar = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);

        values[iVar] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerStates(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerStates(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerStates(unsigned short iMarker, unsigned long iVertex) const {
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    const auto nVar = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);

        values[iVar] = SU2_TYPE::GetValue(value);
    }

    return values;
}

void CDriver::SetStates(vector<vector<passivedouble>> values) {
    const auto nPoint = GetNumberVertices();
    
    if (values.size() != nPoint) {
        SU2_MPI::Error("Invalid number of vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        SetStates(iPoint, values[iPoint]);
    }
}

void CDriver::SetStates(unsigned long iPoint, vector<passivedouble> values) {
    const auto nVar   = GetNumberStateVariables();

    if (values.size() != nVar) {
        SU2_MPI::Error("Invalid number of variables!", CURRENT_FUNCTION);
    }
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->SetSolution(iPoint, iVar, values[iVar]);
    }
}

void CDriver::SetMarkerStates(unsigned short iMarker, vector<vector<passivedouble>> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
    }

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerStates(iMarker, iVertex, values[iVertex]);
    }
}

void CDriver::SetMarkerStates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) { 
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    const auto nVar   = GetNumberStateVariables();
    
    if (values.size() != nVar) {
        SU2_MPI::Error("Invalid number of variables!", CURRENT_FUNCTION);
    }
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->SetSolution(iPoint, iVar, values[iVar]);
    }
}

vector<vector<passivedouble>> CDriver::GetPrimitiveStates() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetPrimitiveStates(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetPrimitiveStates(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    const auto nPrim = nDim + 2;  // TODO: Use nPrimVar ?
    vector<passivedouble> values(nPrim, 0.0);

    CSolver* solver  = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];

    if (nDim == 2) {
        values[0]     = SU2_TYPE::GetValue(solver->GetNodes()->GetDensity(iPoint));
        values[1] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 0));
        values[2] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 1));
        values[3] = SU2_TYPE::GetValue(solver->GetNodes()->GetPressure(iPoint));
    }
    else {
        values[0]     = SU2_TYPE::GetValue(solver->GetNodes()->GetDensity(iPoint));
        values[1] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 0));
        values[2] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 1));
        values[3] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 2));
        values[4] = SU2_TYPE::GetValue(solver->GetNodes()->GetPressure(iPoint));
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerPrimitiveStates(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<vector<passivedouble>> values;

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerPrimitiveStates(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerPrimitiveStates(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    const auto nPrim  = 5;
    vector<passivedouble> values(nPrim, 0.0);

    CSolver* solver  = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];

    if (nDim == 2) {
        values[0]     = SU2_TYPE::GetValue(solver->GetNodes()->GetDensity(iPoint));
        values[1] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 0));
        values[2] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 1));
        values[3] = SU2_TYPE::GetValue(solver->GetNodes()->GetPressure(iPoint));
    }
    else {
        values[0]     = SU2_TYPE::GetValue(solver->GetNodes()->GetDensity(iPoint));
        values[1] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 0));
        values[2] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 1));
        values[3] = SU2_TYPE::GetValue(solver->GetNodes()->GetVelocity(iPoint, 2));
        values[4] = SU2_TYPE::GetValue(solver->GetNodes()->GetPressure(iPoint));
    }

    return values;
}

vector<passivedouble> CDriver::GetSpeedOfSound() const {
    const auto nPoint = GetNumberVertices();
    
    vector<passivedouble> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetSpeedOfSound(iPoint));
    }
    
    return values;
}

passivedouble CDriver::GetSpeedOfSound(unsigned long iPoint) const {
    if (main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
        
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint));
}

vector<passivedouble> CDriver::GetMarkerSpeedOfSound(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<passivedouble> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetSpeedOfSound(iVertex));
    }
    
    return values;
}

passivedouble CDriver::GetMarkerSpeedOfSound(unsigned short iMarker, unsigned long iVertex) const {
    if (main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint));
}

vector<vector<passivedouble>> CDriver::GetMarkerForces(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0u; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerForces(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerForces(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetVertexTractions(iMarker, iVertex, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the adjoint mesh solver solution.                  */
/////////////////////////////////////////////////////////////////////////////

vector<vector<passivedouble>> CDriver::GetAdjointCoordinates() const {
    const auto nPoint = GetNumberVertices();

    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetAdjointCoordinates(iPoint));    
    }
    
    return values;
}

vector<passivedouble> CDriver::GetAdjointCoordinates(unsigned long iPoint) const {
    if (!main_config->GetDeform_Mesh() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }   

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);

    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->GetNodes()->GetSolution(iPoint, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }

    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerAdjointCoordinates(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerAdjointCoordinates(iMarker, iVertex));    
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerAdjointCoordinates(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetDeform_Mesh() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    vector<passivedouble> values(nDim, 0.0);

    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->GetNodes()->GetSolution(iPoint, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }

    return values;
}

void CDriver::SetAdjointCoordinates(vector<vector<passivedouble>> values) {
    const auto nPoint = GetNumberVertices();
    
    if (values.size() != nPoint) {
        SU2_MPI::Error("Invalid number of vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        SetAdjointCoordinates(iPoint, values[iPoint]);
    }
}

void CDriver::SetAdjointCoordinates(unsigned long iPoint, vector<passivedouble> values) {
    if (!main_config->GetDeform_Mesh() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->GetNodes()->SetSolution(iPoint, iDim, values[iDim]);
    }
}

void CDriver::SetMarkerAdjointCoordinates(unsigned short iMarker, vector<vector<passivedouble>> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetAdjointCoordinates(iVertex, values[iVertex]);
    }
}

void CDriver::SetMarkerAdjointCoordinates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    if (!main_config->GetDeform_Mesh() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->GetNodes()->SetSolution(iPoint, iDim, values[iDim]);
    }
}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the adjoint flow solver solution.                  */
/////////////////////////////////////////////////////////////////////////////

vector<vector<passivedouble>> CDriver::GetAdjointStates() const {
    const auto nPoint = GetNumberVertices();

    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetAdjointStates(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetAdjointStates(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    const auto nVar   = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);

        values[iVar] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerAdjointStates(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerAdjointStates(iMarker, iVertex));
    }
    
    return values;
    
}

vector<passivedouble> CDriver::GetMarkerAdjointStates(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    const auto nVar   = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);

    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);

        values[iVar] = SU2_TYPE::GetValue(value);
    }

    return values;
}

void CDriver::SetAdjointStates(vector<vector<passivedouble>> values) {
    const auto nPoint = GetNumberVertices();
    
    if (values.size() != nPoint) {
        SU2_MPI::Error("Invalid number of vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        SetAdjointStates(iPoint, values[iPoint]);
    }
}

void CDriver::SetAdjointStates(unsigned long iPoint, vector<passivedouble> values) {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    const auto nVar = GetNumberStateVariables();

    if (values.size() != nVar) {
        SU2_MPI::Error("Invalid number of variables!", CURRENT_FUNCTION);
    }
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetNodes()->SetSolution(iPoint, iVar, values[iVar]);
    }
}

void CDriver::SetMarkerAdjointStates(unsigned short iMarker, vector<vector<passivedouble>> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerAdjointStates(iMarker, iVertex, values[iVertex]);
    }
}

void CDriver::SetMarkerAdjointStates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex); 
    const auto nVar   = GetNumberStateVariables();
    
    if (values.size() != nVar) {
        SU2_MPI::Error("Invalid number of variables!", CURRENT_FUNCTION);
    }

    for (auto iVar = 0u; iVar < nVar; iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetNodes()->SetSolution(iPoint, iVar, values[iVar]);
    }
}

vector<vector<passivedouble>> CDriver::GetMarkerAdjointForces(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerAdjointForces(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerAdjointForces(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }
    
    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetAdjointVertexTractions(iMarker, iVertex, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

void CDriver::SetMarkerAdjointForces(unsigned short iMarker, vector<vector<passivedouble>> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerAdjointForces(iMarker, iVertex, values[iVertex]);
    }
}

void CDriver::SetMarkerAdjointForces(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->StoreVertexTractionsAdjoint(iMarker, iVertex, iDim, values[iVertex*nDim + iDim]);
    }
}

vector<vector<passivedouble>> CDriver::GetCoordinatesCoordinatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetCoordinatesCoordinatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetCoordinatesCoordinatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dCoordinates_dCoordinates(iPoint, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerCoordinatesDisplacementsSensitivity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerCoordinatesDisplacementsSensitivity(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerCoordinatesDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dCoordinates_dDisplacements(iMarker, iVertex, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<passivedouble> CDriver::GetObjectiveFarfieldVariablesSensitivity() const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
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

vector<passivedouble> CDriver::GetResidualsFarfieldVariablesSensitivity() const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
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

vector<vector<passivedouble>> CDriver::GetObjectiveStatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetObjectiveStatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetObjectiveStatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
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

vector<vector<passivedouble>> CDriver::GetResidualsStatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetResidualsStatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetResidualsStatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    const auto nVar   = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dStates(iPoint, iVar);
        
        values[iVar] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetForcesStatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetForcesStatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetForcesStatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    const auto nVar   = GetNumberStateVariables();
    vector<passivedouble> values(nVar, 0.0);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dTractions_dStates(iPoint, iVar);
        
        values[iVar] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetObjectiveCoordinatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetObjectiveCoordinatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetObjectiveCoordinatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dCoordinates(iPoint, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetResidualsCoordinatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetResidualsCoordinatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetResidualsCoordinatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dCoordinates(iPoint, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetForcesCoordinatesSensitivity() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetForcesCoordinatesSensitivity(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetForcesCoordinatesSensitivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dTractions_dCoordinates(iPoint, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerObjectiveDisplacementsSensitivity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerObjectiveDisplacementsSensitivity(iMarker, iVertex));    
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerObjectiveDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }
    
    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetSens_dObjective_dDisplacements(iMarker, iVertex, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerResidualsDisplacementsSensitivity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerResidualsDisplacementsSensitivity(iMarker, iVertex));    
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerResidualsDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }
    
    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dResiduals_dDisplacements(iMarker, iVertex, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerForcesDisplacementsSensitivity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerForcesDisplacementsSensitivity(iMarker, iVertex));    
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerForcesDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != RESIDUALS) {
        SU2_MPI::Error("Discrete adjoint flow solver does not use residual-based formulation!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }
    
    vector<passivedouble> values(nDim, 0.0);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][ADJFLOW_SOL]->GetProd_dTractions_dDisplacements(iMarker, iVertex, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, etc.).  */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberTimeIterations() const {
    return main_config->GetnTime_Iter();
}

unsigned long CDriver::GetCurrentTimeIteration() const{
    return TimeIter;
}

passivedouble CDriver::GetUnsteadyTimeStep() const {
    return SU2_TYPE::GetValue(main_config->GetTime_Step());
}

string CDriver::GetSurfaceFileName() const {
    return main_config->GetSurfCoeff_FileName();
}

///////////////////////////////////////////////////////////////////////////////
/* Functions related to CHT solver.                                          */
///////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetTemperatures() const {
    const auto nPoint = GetNumberVertices();
    
    vector<passivedouble> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetTemperatures(iPoint));    
    }
    
    return values;
}

passivedouble CDriver::GetTemperatures(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetTemperature(iPoint));
}

vector<passivedouble> CDriver::GetMarkerTemperatures(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<passivedouble> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerTemperatures(iMarker, iVertex));    
    }
    
    return values;
}

passivedouble CDriver::GetMarkerTemperatures(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    if (main_config->GetKind_Regime() != ENUM_REGIME::COMPRESSIBLE) {
        return 0.0;
    }
    else {
        return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetTemperature(iPoint));
    }
}

void CDriver::SetMarkerTemperatures(unsigned short iMarker, vector<passivedouble> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
    }

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerTemperatures(iMarker, iVertex, values[iVertex]);
    }
}

void CDriver::SetMarkerTemperatures(unsigned short iMarker, unsigned long iVertex, passivedouble value) {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }
    
    main_geometry->SetCustomBoundaryTemperature(iMarker, iVertex, value);
}

vector<vector<passivedouble>> CDriver::GetHeatFlux() const {
    const auto nPoint = GetNumberVertices();
    
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetHeatFlux(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetHeatFlux(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values(nDim, 0.0);

    if (main_config->GetKind_Regime() != ENUM_REGIME::COMPRESSIBLE) {
        return values;
    }
    
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    const auto Prandtl_Lam          = main_config->GetPrandtl_Lam();
    const auto Gas_Constant         = main_config->GetGas_ConstantND();
    const auto Gamma                = main_config->GetGamma();
    const auto Cp                   = (Gamma/(Gamma - 1.0))*Gas_Constant;
    const auto laminar_viscosity    = solver->GetNodes()->GetLaminarViscosity(iPoint);
    const auto thermal_conductivity = Cp*(laminar_viscosity/Prandtl_Lam);
        
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        auto GradT = solver->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(-thermal_conductivity*GradT);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetMarkerHeatFlux(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerHeatFlux(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerHeatFlux(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    vector<passivedouble> values(nDim, 0.0);

    if (main_config->GetKind_Regime() != ENUM_REGIME::COMPRESSIBLE) {
        return values;
    }

    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    const auto Prandtl_Lam          = main_config->GetPrandtl_Lam();
    const auto Gas_Constant         = main_config->GetGas_ConstantND();
    const auto Gamma                = main_config->GetGamma();
    const auto Cp                   = (Gamma/(Gamma - 1.0))*Gas_Constant;
    const auto laminar_viscosity    = solver->GetNodes()->GetLaminarViscosity(iPoint);
    const auto thermal_conductivity = Cp*(laminar_viscosity/Prandtl_Lam);

    for (auto iDim = 0u; iDim < nDim; iDim++) {
        auto GradT = solver->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(-thermal_conductivity*GradT);
    }

    return values;
}

vector<passivedouble> CDriver::GetMarkerNormalHeatFlux(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<passivedouble> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerNormalHeatFlux(iMarker, iVertex));
    }
    
    return values;
}

passivedouble CDriver::GetMarkerNormalHeatFlux(unsigned short iMarker, unsigned long iVertex) const {
    vector<passivedouble> values = GetMarkerHeatFlux(iMarker, iVertex);
    passivedouble projected = 0.0;

    const auto Normal = main_geometry->vertex[iMarker][iVertex]->GetNormal();
    const auto Area   = GeometryToolbox::Norm(nDim, Normal);
        
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        projected += values[iDim] * SU2_TYPE::GetValue(Normal[iDim] / Area);
    }

    return projected;
}

void CDriver::SetMarkerNormalHeatFlux(unsigned short iMarker, vector<passivedouble> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
    }

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerNormalHeatFlux(iMarker, iVertex, values[iVertex]);
    }
}

void CDriver::SetMarkerNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble value) {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }

    main_geometry->SetCustomBoundaryHeatFlux(iMarker, iVertex, value);
}

vector<passivedouble> CDriver::GetThermalConductivity() const {
    const auto nPoint = GetNumberVertices();

    vector<passivedouble> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetThermalConductivity(iPoint));
    }
    
    return values;
}

passivedouble CDriver::GetThermalConductivity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }

    const auto Prandtl_Lam       = main_config->GetPrandtl_Lam();
    const auto Gas_Constant      = main_config->GetGas_ConstantND();
    const auto Gamma             = main_config->GetGamma();
    const auto Cp                = (Gamma/(Gamma - 1.0))*Gas_Constant;
    
    const auto laminar_viscosity = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    
    return SU2_TYPE::GetValue(Cp*(laminar_viscosity/Prandtl_Lam));
}

vector<passivedouble> CDriver::GetMarkerThermalConductivity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<passivedouble> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerThermalConductivity(iMarker, iVertex));
    }
    
    return values;
}

passivedouble CDriver::GetMarkerThermalConductivity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }

    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    const auto Prandtl_Lam  = main_config->GetPrandtl_Lam();
    const auto Gas_Constant = main_config->GetGas_ConstantND();
    const auto Gamma        = main_config->GetGamma();
    const auto Cp           = (Gamma/(Gamma - 1.0))*Gas_Constant;
    
    const auto laminar_viscosity = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        
    return SU2_TYPE::GetValue(Cp*(laminar_viscosity/Prandtl_Lam));
}

vector<passivedouble> CDriver::GetLaminarViscosity() const {
    const auto nPoint = GetNumberVertices();

    vector<passivedouble> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetLaminarViscosity(iPoint));
    }
    
    return values;
}

passivedouble CDriver::GetLaminarViscosity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint));
}

vector<passivedouble> CDriver::GetMarkerLaminarViscosity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<passivedouble> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerLaminarViscosity(iMarker, iVertex));
    }
    
    return values;
}

passivedouble CDriver::GetMarkerLaminarViscosity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint));
}

vector<passivedouble> CDriver::GetEddyViscosity() const {
    const auto nPoint = GetNumberVertices();

    vector<passivedouble> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetEddyViscosity(iPoint));
    }
    
    return values;
}

passivedouble CDriver::GetEddyViscosity(unsigned long iPoint) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
    }
    
    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint));
}

vector<passivedouble> CDriver::GetMarkerEddyViscosity(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<passivedouble> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerLaminarViscosity(iMarker, iVertex));
    }
    
    return values;
}

passivedouble CDriver::GetMarkerEddyViscosity(unsigned short iMarker, unsigned long iVertex) const {
    if (!main_config->GetFluidProblem()) {
        SU2_MPI::Error("Flow solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    return SU2_TYPE::GetValue(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint));
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to nonequilibrium flow solver.                           */
////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberNonequilibriumSpecies() const {
    return main_config->GetnSpecies();
}

unsigned long CDriver::GetNumberNonequilibriumStateVariables() const {
    return GetNumberNonequilibriumSpecies() + nDim + 2;
}

unsigned short CDriver::GetNumberNonequilibriumPrimitiveVariables() const {
    unsigned short nPrim;
    
    if (main_config->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES) {
        nPrim = GetNumberNonequilibriumSpecies() + nDim + 10;
    } else {
        nPrim = GetNumberNonequilibriumSpecies() + nDim + 8;
    }
    
    return nPrim;
}

vector<vector<passivedouble>> CDriver::GetNonequilibriumMassFractions() const {
    if (!main_config->GetNEMOProblem()) {
        return {};
    }
    
    const auto nPoint = GetNumberVertices();
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetNonequilibriumMassFractions(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetNonequilibriumMassFractions(unsigned long iPoint) const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!main_config->GetNEMOProblem()) {
        return {};
    }
    
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }
    
    const auto nSpecies = GetNumberNonequilibriumSpecies();
    vector<passivedouble> values;
    
    for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++) {
        auto rho_s = solver->GetNodes()->GetSolution(iPoint, iSpecies);
        auto rho_t = solver->GetNodes()->GetDensity(iPoint);
        su2double value = rho_s/rho_t;
        
        values.push_back(SU2_TYPE::GetValue(value));
    }
    
    return values;
}

vector<vector<passivedouble>> CDriver::GetNonequilibriumStates() const {
    if (!main_config->GetNEMOProblem()) {
        return {};
    }
    
    const auto nPoint = GetNumberVertices();
    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetNonequilibriumStates(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriver::GetNonequilibriumStates(unsigned long iPoint) const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!main_config->GetNEMOProblem()) {
        return {};
    }
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }
    
    const auto nVar = GetNumberNonequilibriumStateVariables();
    vector<passivedouble> values;
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        su2double value = solver->GetNodes()->GetSolution(iPoint, iVar);
        
        values.push_back(SU2_TYPE::GetValue(value));
    }
    
    return values;
}

void CDriver::SetNonequilibriumStates(vector<vector<passivedouble>> values) {
    if (!main_config->GetNEMOProblem()) {
        SU2_MPI::Error("Nonequilibrium low solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nPoint = GetNumberVertices();
    
    if (values.size() != nPoint) {
        SU2_MPI::Error("Invalid number of vertices!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        SetNonequilibriumStates(iPoint, values[iPoint]);
    }
    
}

void CDriver::SetNonequilibriumStates(unsigned long iPoint, vector<passivedouble> values) {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!main_config->GetNEMOProblem()) {
        SU2_MPI::Error("Nonequilibrium flow solver is not defined!", CURRENT_FUNCTION);
    }
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }
    
    const auto nVar = GetNumberNonequilibriumStateVariables();
    
    if (values.size() != nVar) {
        SU2_MPI::Error("Invalid number of state variables!", CURRENT_FUNCTION);
    }
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
        solver->GetNodes()->SetSolution(iPoint, iVar, values[iVar]);
    }
}

vector<passivedouble> CDriver::GetVibrationalTemperatures() const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    if (!main_config->GetNEMOProblem()) {
        return {};
    }
    
    const auto nPoint = GetNumberVertices();
    vector<passivedouble> values(nPoint, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        value = solver->GetNodes()->GetTemperature_ve(iPoint);
        values[iPoint] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers.                            */
////////////////////////////////////////////////////////////////////////////////

vector<string> CDriver::GetFluidLoadMarkerTags() const {
    const auto nMarker = main_config->GetnMarker_Fluid_Load();
    vector<string> tags;
    
    tags.resize(nMarker);
    
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
        tags[iMarker] = main_config->GetMarker_Fluid_Load_TagBound(iMarker);
    }
    
    return tags;
}

vector<string> CDriver::GetCHTMarkerTags() const {
    const auto nMarker = GetNumberMarkers();
    vector<string> tags;
    
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
        if ((main_config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || main_config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) && main_config->GetMarker_All_PyCustom(iMarker)) {
            tags.push_back(main_config->GetMarker_All_TagBound(iMarker));
        }
    }
    
    return tags;
}

vector<string> CDriver::GetInletMarkerTags() const {
    const auto nMarker = GetNumberMarkers();
    vector<string> tags;
    
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
        bool isCustomizable = main_config->GetMarker_All_PyCustom(iMarker);
        bool isInlet = (main_config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
        if (isCustomizable && isInlet) {
            tags.push_back(main_config->GetMarker_All_TagBound(iMarker));
        }
    }
    
    return tags;
}

void CDriver::SetHeatSourcePosition(passivedouble alpha, passivedouble pos_x, passivedouble pos_y, passivedouble pos_z) {
    CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][RAD_SOL];
    
    main_config->SetHeatSource_Rot_Z(alpha);
    main_config->SetHeatSource_Center(pos_x, pos_y, pos_z);
    solver->SetVolumetricHeatSource(main_geometry, main_config);
}

void CDriver::SetInletAngle(unsigned short iMarker, passivedouble alpha) {
    CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
    
    su2double alpha_rad = alpha * PI_NUMBER/180.0;
    
    for (auto iVertex = 0ul; iVertex < main_geometry->nVertex[iMarker]; iVertex++) {
        solver->SetInlet_FlowDir(iMarker, iVertex, 0, cos(alpha_rad));
        solver->SetInlet_FlowDir(iMarker, iVertex, 1, sin(alpha_rad));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, etc.).  */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CDriver::ResetConvergence() {
    for(auto iZone = 0u; iZone < nZone; iZone++) {
        switch (main_config->GetKind_Solver()) {
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
    CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    DynamicMeshUpdate(0);
    
    SU2_OMP_PARALLEL {
        // Overwrite fictious velocities
        for (iMesh = 0u; iMesh <= main_config->GetnMGLevels(); iMesh++) {
            SU2_OMP_FOR_STAT(roundUpDiv(geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(), omp_get_max_threads()))
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
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!main_config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        su2double NodalForce[3] = {0.0, 0.0, 0.0};
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            NodalForce[iDim] = values[iVertex*nDim + iDim];
        }
        
        solver->GetNodes()->Set_FlowTraction(iPoint, NodalForce);
    }
}

vector<passivedouble> CDriver::GetMarkerFEADisplacements(unsigned short iMarker) const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!main_config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}


vector<passivedouble> CDriver::GetMarkerFEAVelocity(unsigned short iMarker) const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!main_config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetDynamic_Analysis() != DYNAMIC) {
        return {};
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetSolution_Vel(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerCurrentFEAVelocity(unsigned short iMarker) const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
    
    if (!main_config->GetStructuralProblem()) {
        SU2_MPI::Error("Structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetDynamic_Analysis() != DYNAMIC) {
        return {};
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
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
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
    
    if (!main_config->GetDeform_Mesh() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint mesh solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Discrete adjoint mesh solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBoundDisp_Sens(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriver::GetMarkerForcesSensitivity(unsigned short iMarker) const {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
    
    if (!main_config->GetStructuralProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Discrete adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriver::SetMarkerAdjointDisplacementSourceTerm(unsigned short iMarker, vector<passivedouble> values) {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
    
    if (!main_config->GetStructuralProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Discrete adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

void CDriver::SetMarkerAdjointVelocitySourceTerm(unsigned short iMarker, vector<passivedouble> values) {
    CSolver* solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
    
    if (!main_config->GetStructuralProblem() || !main_config->GetDiscrete_Adjoint()) {
        SU2_MPI::Error("Discrete adjoint structural solver is not defined!", CURRENT_FUNCTION);
    }
    if (main_config->GetKind_DiscreteAdjoint() != FIXED_POINT) {
        SU2_MPI::Error("Discrete adjoint structural solver does not use fixed-point formulation!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = main_geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}
