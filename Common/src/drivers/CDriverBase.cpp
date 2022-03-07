/*!
 * \file CDriverBase.hpp
 * \brief Base class template for all drivers.
 * \author H. Patel, A. Gastaldi
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser/ General Public
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

#include "../../include/drivers/CDriverBase.hpp"

#include "../../include/geometry/CPhysicalGeometry.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

using namespace std;

CDriverBase::CDriverBase(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator): 
config_file_name(confFile),
StartTime(0.0),
StopTime(0.0),
UsedTime(0.0),
TimeIter(0),
nZone(val_nZone)
{
    
}

CDriverBase::~CDriverBase(void) {
    
}

void CDriverBase::SetContainers_Null() {
    
    /*--- Create pointers to all of the classes that may be used by drivers. In general, the pointers are instantiated down a
     hierarchy over all zones, multigrid levels, equation sets, and equation
     terms as described in the comments below. ---*/
    
    config_container      = nullptr;
    output_container      = nullptr;
    geometry_container    = nullptr;
    solver_container      = nullptr;
    numerics_container    = nullptr;
    
    surface_movement      = nullptr;
    grid_movement         = nullptr;
    FFDBox                = nullptr;
    
    config_container      = new CConfig*[nZone] ();
    output_container      = new COutput*[nZone] ();
    geometry_container    = new CGeometry***[nZone] ();
    solver_container      = new CSolver****[nZone] ();
    numerics_container    = new CNumerics*****[nZone] ();
    surface_movement      = new CSurfaceMovement*[nZone] ();
    grid_movement         = new CVolumetricMovement**[nZone] ();
    FFDBox                = new CFreeFormDefBox**[nZone] ();
    
    nInst                 = new unsigned short[nZone];
    
    for (iZone = 0; iZone < nZone; iZone++) {
        nInst[iZone] = 1;
    }
    
    driver_config         = nullptr;
    driver_output         = nullptr;
}

unsigned short CDriverBase::GetNumberMarkers() const {
    return config_container[ZONE_0]->GetnMarker_All();
}

map<string, unsigned short> CDriverBase::GetMarkerIndices() const {
    const auto nMarker = config_container[ZONE_0]->GetnMarker_All();
    map<string, unsigned short> indexMap;
    
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
        auto tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);

        indexMap[tag] = iMarker;
    }
    
    return indexMap;
}

map<string, string> CDriverBase::GetMarkerTypes() const {
    map<string, string> typeMap;
    string type;
    
    for (auto iMarker = 0u; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
        auto tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
        auto kindBC = config_container[ZONE_0]->GetMarker_All_KindBC(iMarker);
        
        switch(kindBC) {
            case EULER_WALL:
                type = "EULER_WALL";
                break;
            case FAR_FIELD:
                type = "FARFIELD";
                break;
            case ISOTHERMAL:
                type = "ISOTHERMAL";
                break;
            case HEAT_FLUX:
                type = "HEATFLUX";
                break;
            case INLET_FLOW:
                type = "INLET_FLOW";
                break;
            case OUTLET_FLOW:
                type = "OUTLET_FLOW";
                break;
            case SYMMETRY_PLANE:
                type = "SYMMETRY";
                break;
            case SEND_RECEIVE:
                type = "SEND_RECEIVE";
                break;
            default:
                type = "UNKNOWN_TYPE";
        }
        typeMap[tag] = type;
    }
    
    return typeMap;
}

vector<string> CDriverBase::GetMarkerTags() const {
    vector<string> tags;
    const auto nMarker = config_container[ZONE_0]->GetnMarker_All();
    
    for(auto iMarker = 0u; iMarker < nMarker; iMarker++){
        tags.push_back(config_container[ZONE_0]->GetMarker_All_TagBound(iMarker));
    }
    
    return tags;
}

vector<string> CDriverBase::GetDeformableMarkerTags() const {
    vector<string> tags;
    const auto nMarker = config_container[ZONE_0]->GetnMarker_Deform_Mesh();
    
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
        tags.push_back(config_container[ZONE_0]->GetMarker_Deform_Mesh_TagBound(iMarker));
    }
    
    return tags;
}

unsigned long CDriverBase::GetNumberDimensions() const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
}

unsigned long CDriverBase::GetNumberElements() const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnElem();
}

unsigned long CDriverBase::GetNumberMarkerElements(unsigned short iMarker) const {
    if (iMarker >= GetNumberMarkers()) {
        SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnElem_Bound(iMarker);
}

vector<unsigned long> CDriverBase::GetElementIDs() const {
    const auto nElem = GetNumberElements();

    vector<unsigned long> values;

    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        values.push_back(GetElementIDs(iElem));
    }
    
    return values;
}

unsigned long CDriverBase::GetElementIDs(unsigned long iElem) const {
    if (iElem >= GetNumberElements()) {
        SU2_MPI::Error("Element index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->elem[iElem]->GetGlobalIndex();
}

vector<unsigned long> CDriverBase::GetMarkerElementIDs(unsigned short iMarker) const { 
    const auto nBound = GetNumberMarkerElements(iMarker);

    vector<unsigned long> values;
     
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        values.push_back(GetMarkerElementIDs(iMarker, iBound));
    }
    
    return values;
}

unsigned long CDriverBase::GetMarkerElementIDs(unsigned short iMarker, unsigned long iBound) const {
    if (iBound >= GetNumberMarkerElements(iMarker)) {
        SU2_MPI::Error("Marker element index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->bound[iMarker][iBound]->GetGlobalIndex();
}

vector<unsigned long> CDriverBase::GetElementColors() const {
    const auto nElem = GetNumberElements();

    vector<unsigned long> values;

    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        values.push_back(GetElementColors(iElem));
    }
    
    return values;
}

unsigned long CDriverBase::GetElementColors(unsigned long iElem) const {
    if (iElem >= GetNumberElements()) {
        SU2_MPI::Error("Element index exceeds size.", CURRENT_FUNCTION);
    }
    
    return geometry_container[ZONE_0][INST_0][MESH_0]->elem[iElem]->GetColor();
}

vector<unsigned long> CDriverBase::GetMarkerElementColors(unsigned short iMarker) const {
    const auto nBound = GetNumberMarkerElements(iMarker);

    vector<unsigned long> values;
    
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        values.push_back(GetMarkerElementColors(iMarker, iBound));
    }
    
    return values;
}

unsigned long CDriverBase::GetMarkerElementColors(unsigned short iMarker, unsigned long iBound) const {
    if (iBound >= GetNumberMarkerElements(iMarker)) {
        SU2_MPI::Error("Marker element index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->bound[iMarker][iBound]->GetColor();
}

vector<vector<unsigned long>> CDriverBase::GetElementConnectivities() const {
    const auto nElem = GetNumberElements();
    
    vector<vector<unsigned long>> values;
    
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        values.push_back(GetElementConnectivities(iElem));
    }
    
    return values;
}

vector<unsigned long> CDriverBase::GetElementConnectivities(unsigned long iElem) const {
    if (iElem >= GetNumberElements()) {
        SU2_MPI::Error("Element index exceeds size.", CURRENT_FUNCTION);
    }
    
    unsigned short nNode = geometry_container[ZONE_0][INST_0][MESH_0]->elem[iElem]->GetnNodes();
    
    vector<unsigned long> values; 

    for (auto iNode = 0u; iNode < nNode; iNode++) {
        values.push_back(geometry_container[ZONE_0][INST_0][MESH_0]->elem[iElem]->GetNode(iNode));
    }
    
    return values;
}

vector<vector<unsigned long>> CDriverBase::GetMarkerElementConnectivities(unsigned short iMarker) const {
    const auto nBound = GetNumberMarkerElements(iMarker);
    
    vector<vector<unsigned long>> values;

    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        values.push_back(GetMarkerElementConnectivities(iMarker, iBound));
    }
    
    return values;
}

vector<unsigned long> CDriverBase::GetMarkerElementConnectivities(unsigned short iMarker, unsigned long iBound) const {
    if (iBound >= GetNumberMarkerElements(iMarker)) {
        SU2_MPI::Error("Marker element index exceeds size.", CURRENT_FUNCTION);
    }
    
    unsigned short nNode = geometry_container[ZONE_0][INST_0][MESH_0]->bound[iMarker][iBound]->GetnNodes();
    
    vector<unsigned long> values; 

    for (auto iNode = 0u; iNode < nNode; iNode++) {
        values.push_back(geometry_container[ZONE_0][INST_0][MESH_0]->bound[iMarker][iBound]->GetNode(iNode));
    }
    
    return values;
}

unsigned long CDriverBase::GetNumberVertices() const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
}

unsigned long CDriverBase::GetNumberMarkerVertices(unsigned short iMarker) const {
    if (iMarker >= GetNumberMarkers()) {
        SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnVertex(iMarker);
}

unsigned long CDriverBase::GetNumberHaloVertices() const {
    const auto nPoint = GetNumberVertices();
    unsigned long nHalo = 0;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        if (!(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint))) {
            nHalo += 1;
        }
    }
    
    return nHalo;
}

unsigned long CDriverBase::GetNumberMarkerHaloVertices(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    unsigned long nHalo = 0;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
        
        if (!(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint))) {
            nHalo += 1;
        }
    }
    
    return nHalo;
}

vector<unsigned long> CDriverBase::GetMarkerVertexIndex(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<unsigned long> values;

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerVertexIndex(iMarker, iVertex));
    }
    
    return values;
}

unsigned long CDriverBase::GetMarkerVertexIndex(unsigned short iMarker, unsigned long iVertex) const {
    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Marker vertex index exceeds size.", CURRENT_FUNCTION);
    }
    
    return geometry_container[MESH_0][INST_0][ZONE_0]->vertex[iMarker][iVertex]->GetNode();
}

vector<unsigned long> CDriverBase::GetVertexIDs() const {
    const auto nPoint = GetNumberVertices();
    
    vector<unsigned long> values;

    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetVertexIDs(iPoint));
    }
    
    return values;
}

unsigned long CDriverBase::GetVertexIDs(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetGlobalIndex(iPoint);
}

vector<unsigned long> CDriverBase::GetMarkerVertexIDs(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker); 
    
    vector<unsigned long> values;

    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerVertexIDs(iMarker, iVertex));
    }
    
    return values;
}

unsigned long CDriverBase::GetMarkerVertexIDs(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    return geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetGlobalIndex(iPoint);
}

vector<unsigned long> CDriverBase::GetVertexColors() const {
    const auto nPoint = GetNumberVertices();
    
    vector<unsigned long> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetVertexColors(iPoint));
    }
    
    return values;
}

unsigned long CDriverBase::GetVertexColors(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }

    return geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetColor(iPoint);
}

vector<unsigned long> CDriverBase::GetMarkerVertexColors(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);
    
    vector<unsigned long> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerVertexColors(iMarker, iVertex));
    }
    
    return values;
}

unsigned long CDriverBase::GetMarkerVertexColors(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    return geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetColor(iPoint);
}

vector<bool> CDriverBase::GetDomain() const {
    const auto nPoint = GetNumberVertices();

    vector<bool> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetDomain(iPoint));
    }
    
    return values;
}

bool CDriverBase::GetDomain(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }
    
    return geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint);
}

vector<bool> CDriverBase::GetMarkerDomain(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<bool> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerDomain(iMarker, iVertex));
    }
    
    return values;
}

bool CDriverBase::GetMarkerDomain(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);
    
    return geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint);
}

vector<vector<passivedouble>> CDriverBase::GetInitialCoordinates() const {
    const auto nPoint = GetNumberVertices();

    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetInitialCoordinates(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetInitialCoordinates(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values (nDim, 0.0);
    
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        return values;
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord(iPoint, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriverBase::GetMarkerInitialCoordinates(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerInitialCoordinates(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetMarkerInitialCoordinates(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    vector<passivedouble> values(nDim, 0.0);

    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        return values;
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord(iPoint, iDim);

        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

vector<vector<passivedouble>> CDriverBase::GetCoordinates() const {
    const auto nPoint = GetNumberVertices();

    vector<vector<passivedouble>> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(GetCoordinates(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetCoordinates(unsigned long iPoint) const {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values;
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetCoord(iPoint, iDim);
        
        values.push_back(SU2_TYPE::GetValue(value));
    }
    
    return values;
}

vector<vector<passivedouble>> CDriverBase::GetMarkerCoordinates(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerCoordinates(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    vector<passivedouble> values;
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetCoord(iPoint, iDim);
        
        values.push_back(SU2_TYPE::GetValue(value));
    }
    
    return values;
}

void CDriverBase::SetCoordinates(vector<vector<passivedouble>> values) {
    const auto nPoint = GetNumberVertices();
    
    if (values.size() != nPoint) {
        SU2_MPI::Error("Invalid number of vertices !", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        SetCoordinates(iPoint, values[iPoint]);
    }
}

void CDriverBase::SetCoordinates(unsigned long iPoint, vector<passivedouble> values) {
    if (iPoint >= GetNumberVertices()) {
        SU2_MPI::Error("Vertex index exceeds size.", CURRENT_FUNCTION);
    }

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions !", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetCoord(iPoint, iDim, values[iDim]);
    }
}

void CDriverBase::SetMarkerCoordinates(unsigned short iMarker, vector<vector<passivedouble>> values) {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices !", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerCoordinates(iMarker, iVertex, values[iVertex]);
    }
}

void CDriverBase::SetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions !", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetCoord(iPoint, iDim, values[iDim]);
    }
}

vector<vector<passivedouble>> CDriverBase::GetMarkerDisplacements(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerDisplacements(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    vector<passivedouble> values (nDim, 0.0);
    
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        return values;
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetBound_Disp(iPoint, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

void CDriverBase::SetMarkerDisplacements(unsigned short iMarker, vector<vector<passivedouble>> values) {
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }

    const auto nVertex = GetNumberMarkerVertices(iMarker);

    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices !", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerDisplacements(iMarker, iVertex, values[iVertex]);
    }
}

void CDriverBase::SetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }

    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions !", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iPoint, iDim, values[iDim]);
    }
}

vector<vector<passivedouble>> CDriverBase::GetMarkerVelocities(unsigned short iMarker) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerVelocities(iMarker, iVertex));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetMarkerVelocities(unsigned short iMarker, unsigned long iVertex) const {
    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    vector<passivedouble> values (nDim, 0.0);
    
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        return values;
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetBound_Vel(iPoint, iDim);
        
        values[iDim] = SU2_TYPE::GetValue(value);
    }
    
    return values;
}

void CDriverBase::SetMarkerVelocities(unsigned short iMarker, vector<vector<passivedouble>> values) {
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }

    const auto nVertex = GetNumberMarkerVertices(iMarker);

    if (values.size() != nVertex) {
        SU2_MPI::Error("Invalid number of marker vertices !", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        SetMarkerVelocities(iMarker, iVertex, values[iVertex]);
    }
}

void CDriverBase::SetMarkerVelocities(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    if (!config_container[ZONE_0]->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined !", CURRENT_FUNCTION);
    }

    auto iPoint = GetMarkerVertexIndex(iMarker, iVertex);

    if (values.size() != nDim) {
        SU2_MPI::Error("Invalid number of dimensions !", CURRENT_FUNCTION);
    }
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iPoint, iDim, values[iDim]);
    }
}

vector<vector<passivedouble>> CDriverBase::GetMarkerVertexNormals(unsigned short iMarker, bool normalize) const {
    const auto nVertex = GetNumberMarkerVertices(iMarker);

    vector<vector<passivedouble>> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        values.push_back(GetMarkerVertexNormals(iMarker, iVertex, normalize=normalize));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetMarkerVertexNormals(unsigned short iMarker, unsigned long iVertex, bool normalize) const {
    if (iVertex >= GetNumberMarkerVertices(iMarker)) {
        SU2_MPI::Error("Marker vertex index exceeds size.", CURRENT_FUNCTION);
    }

    vector<passivedouble> values;

    auto normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
    auto area   = GeometryToolbox::Norm(nDim, normal);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
        if (normalize) {
            values[iDim] = SU2_TYPE::GetValue(normal[iDim]/area);
        } else {
            values[iDim] = SU2_TYPE::GetValue(normal[iDim]);
        }
    }

    return values;
}


void CDriverBase::CommunicateMeshDisplacements(void) {
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], MESH_DISPLACEMENTS);
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], MESH_DISPLACEMENTS);
}
