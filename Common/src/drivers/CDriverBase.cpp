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

map<string, int> CDriverBase::GetBoundaryMarkerIndices() const {
    CConfig* config = config_container[ZONE_0];
    
    const auto nBoundaryMarkers = config->GetnMarker_All();
    map<string, int>  allBoundariesMap;
    
    for (auto iMarker = 0u; iMarker < nBoundaryMarkers; iMarker++) {
        auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        allBoundariesMap[Marker_Tag] = iMarker;
    }
    
    return allBoundariesMap;
}

map<string, string> CDriverBase::GetBoundaryMarkerTypes() const {
    CConfig* config = config_container[ZONE_0];
    
    map<string, string> allBoundariesTypeMap;
    string Marker_Type;
    
    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
        auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        auto KindBC = config->GetMarker_All_KindBC(iMarker);
        
        switch(KindBC) {
            case EULER_WALL:
                Marker_Type = "EULER_WALL";
                break;
            case FAR_FIELD:
                Marker_Type = "FARFIELD";
                break;
            case ISOTHERMAL:
                Marker_Type = "ISOTHERMAL";
                break;
            case HEAT_FLUX:
                Marker_Type = "HEATFLUX";
                break;
            case INLET_FLOW:
                Marker_Type = "INLET_FLOW";
                break;
            case OUTLET_FLOW:
                Marker_Type = "OUTLET_FLOW";
                break;
            case SYMMETRY_PLANE:
                Marker_Type = "SYMMETRY";
                break;
            case SEND_RECEIVE:
                Marker_Type = "SEND_RECEIVE";
                break;
            default:
                Marker_Type = "UNKNOWN_TYPE";
        }
        allBoundariesTypeMap[Marker_Tag] = Marker_Type;
    }
    
    return allBoundariesTypeMap;
}

vector<string> CDriverBase::GetDeformableMarkerTags() const {
    CConfig* config = config_container[ZONE_0];
    
    const auto nBoundariesMarker = config->GetnMarker_Deform_Mesh();
    vector<string> interfaceBoundariesTagList;
    
    interfaceBoundariesTagList.resize(nBoundariesMarker);
    
    for (auto iMarker = 0u; iMarker < nBoundariesMarker; iMarker++) {
        auto Marker_Tag = config->GetMarker_Deform_Mesh_TagBound(iMarker);
        interfaceBoundariesTagList[iMarker] = Marker_Tag;
    }
    
    return interfaceBoundariesTagList;
}

unsigned long CDriverBase::GetNumberDimensions() const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
}

unsigned long CDriverBase::GetNumberElements() const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnElem();
}

unsigned long CDriverBase::GetNumberElementsMarker(unsigned short iMarker) const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnElem_Bound(iMarker);
}

unsigned long CDriverBase::GetNumberVertices() const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
}

unsigned long CDriverBase::GetNumberVerticesMarker(unsigned short iMarker) const {
    return geometry_container[ZONE_0][INST_0][MESH_0]->GetnVertex(iMarker);
}

unsigned long CDriverBase::GetNumberHaloVertices() const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nPoint = geometry->GetnPoint();
    unsigned long nHaloVertices = 0;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        if (!(geometry->nodes->GetDomain(iPoint))) {
            nHaloVertices += 1;
        }
    }
    
    return nHaloVertices;
}

unsigned long CDriverBase::GetNumberHaloVerticesMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    unsigned long nHaloVertices = 0;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!(geometry->nodes->GetDomain(iPoint))) {
            nHaloVertices += 1;
        }
    }
    
    return nHaloVertices;
}

vector<unsigned long> CDriverBase::GetVertexIDs() const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nPoint = geometry->GetnPoint();
    vector<unsigned long> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(geometry->nodes->GetGlobalIndex(iPoint));
    }
    
    return values;
}

vector<unsigned long> CDriverBase::GetVertexIDsMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<unsigned long> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        values.push_back(geometry->nodes->GetGlobalIndex(iPoint));
    }
    
    return values;
}

vector<unsigned long> CDriverBase::GetElementIDs() const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nElem = geometry->GetnElem();
    vector<unsigned long> values;
    
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        values.push_back(geometry->elem[iElem]->GetGlobalIndex());
    }
    
    return values;
}

vector<unsigned long> CDriverBase::GetElementIDsMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nBound = geometry->GetnElem_Bound(iMarker);
    vector<unsigned long> values;
    
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        values.push_back(geometry->bound[iMarker][iBound]->GetGlobalIndex());
    }
    
    return values;
}

vector<vector<unsigned long>> CDriverBase::GetConnectivity() const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nElem = geometry->GetnElem();
    vector<vector<unsigned long>> values(nElem);
    
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
        unsigned short nNode = geometry->elem[iElem]->GetnNodes();
        
        for (auto iNode = 0u; iNode < nNode; iNode++) {
            values[iElem].push_back(geometry->elem[iElem]->GetNode(iNode));
        }
    }
    
    return values;
}

vector<vector<unsigned long>> CDriverBase::GetConnectivityMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nBound = geometry->GetnElem_Bound(iMarker);
    vector<vector<unsigned long>> values(nBound);
    
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        unsigned short nNode = geometry->bound[iMarker][iBound]->GetnNodes();
        
        for (auto iNode = 0u; iNode < nNode; iNode++) {
            values[iBound].push_back(geometry->bound[iMarker][iBound]->GetNode(iNode));
        }
    }
    
    return values;
}

vector<bool> CDriverBase::GetDomain() const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nPoint = geometry->GetnPoint();
    vector<bool> values;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        values.push_back(geometry->nodes->GetDomain(iPoint));
    }
    
    return values;
}

vector<bool> CDriverBase::GetDomainMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<bool> values;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        values.push_back(geometry->nodes->GetDomain(iPoint));
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetCoordinates() const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nPoint = geometry->GetnPoint();
    vector<passivedouble> values(nPoint*nDim, 0.0);
    su2double value;
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = geometry->nodes->GetCoord(iPoint, iDim);
            values[iPoint*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetCoordinatesMarker(unsigned short iMarker) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = geometry->nodes->GetCoord(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriverBase::SetCoordinates(vector<passivedouble> values) {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nPoint = geometry->GetnPoint();
    if (values.size() != nPoint*nDim) {
        SU2_MPI::Error("Size does not match nPoint * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            geometry->nodes->SetCoord(iPoint, iDim, values[iPoint*nDim + iDim]);
        }
    }
}

void CDriverBase::SetCoordinatesMarker(unsigned short iMarker, vector<passivedouble> values) {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            geometry->nodes->SetCoord(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

vector<passivedouble> CDriverBase::GetDisplacementsMarker(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    if (!config->GetDeform_Mesh()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBound_Disp(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriverBase::SetDisplacementsMarker(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    if (!config->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetBound_Disp(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

vector<passivedouble> CDriverBase::GetVelocitiesMarker(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    if (!config->GetDeform_Mesh()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetBound_Vel(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

void CDriverBase::SetVelocitiesMarker(unsigned short iMarker, vector<passivedouble> values) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    if (!config->GetDeform_Mesh()) {
        SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    if (values.size() != nVertex*nDim) {
        SU2_MPI::Error("Size does not match nVertex * nDim!", CURRENT_FUNCTION);
    }
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            solver->GetNodes()->SetBound_Vel(iPoint, iDim, values[iVertex*nDim + iDim]);
        }
    }
}

vector<passivedouble> CDriverBase::GetInitialCoordinatesMarker(unsigned short iMarker) const {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    if (!config->GetDeform_Mesh()) {
        return {};
    }
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    su2double value;
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            value = solver->GetNodes()->GetMesh_Coord(iPoint, iDim);
            values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(value);
        }
    }
    
    return values;
}

vector<passivedouble> CDriverBase::GetVertexNormalsMarker(unsigned short iMarker, bool UnitNormal) const {
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    
    const auto nVertex = geometry->GetnVertex(iMarker);
    vector<passivedouble> values(nVertex*nDim, 0.0);
    
    for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
        auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        auto Area   = GeometryToolbox::Norm(nDim, Normal);
        
        for (auto iDim = 0u; iDim < nDim; iDim++) {
            if (!UnitNormal) {
                values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(Normal[iDim]);
            } else {
                values[iVertex*nDim + iDim] = SU2_TYPE::GetValue(Normal[iDim]/Area);
            }
        }
    }
    
    return values;
}

void CDriverBase::CommunicateMeshDisplacements(void) {
    CConfig* config     = config_container[ZONE_0];
    CGeometry* geometry = geometry_container[ZONE_0][INST_0][MESH_0];
    CSolver* solver     = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];
    
    solver->InitiateComms(geometry, config, MESH_DISPLACEMENTS);
    solver->CompleteComms(geometry, config, MESH_DISPLACEMENTS);
}
