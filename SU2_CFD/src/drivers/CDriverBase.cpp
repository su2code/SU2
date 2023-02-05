/*!
 * \file CDriverBase.hpp
 * \brief Base class template for all drivers.
 * \author H. Patel, A. Gastaldi
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

using namespace std;

CDriverBase::CDriverBase(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator)
    : config_file_name(confFile), StartTime(0.0), StopTime(0.0), UsedTime(0.0), TimeIter(0), nZone(val_nZone) {}

CDriverBase::~CDriverBase(void) {}

void CDriverBase::SetContainers_Null() {
  /*--- Create pointers to all the classes that may be used by drivers. In general, the pointers are instantiated
   * down a hierarchy over all zones, multi-grid levels, equation sets, and equation terms as described in the comments
   * below. ---*/

  config_container = nullptr;
  output_container = nullptr;
  geometry_container = nullptr;
  solver_container = nullptr;
  numerics_container = nullptr;

  surface_movement = nullptr;
  grid_movement = nullptr;
  FFDBox = nullptr;

  config_container = new CConfig*[nZone]();
  output_container = new COutput*[nZone]();
  geometry_container = new CGeometry***[nZone]();
  solver_container = new CSolver****[nZone]();
  numerics_container = new CNumerics*****[nZone]();
  surface_movement = new CSurfaceMovement*[nZone]();
  grid_movement = new CVolumetricMovement**[nZone]();
  FFDBox = new CFreeFormDefBox**[nZone]();

  nInst = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }

  driver_config = nullptr;
  driver_output = nullptr;

  main_config = nullptr;
  main_geometry = nullptr;
}

unsigned short CDriverBase::GetNumberDesignVariables() const { return main_config->GetnDV(); }

unsigned short CDriverBase::GetNumberFFDBoxes() const { return main_config->GetnFFDBox(); }

unsigned long CDriverBase::GetNumberDimensions() const { return main_geometry->GetnDim(); }

unsigned long CDriverBase::GetNumberElements() const { return main_geometry->GetnElem(); }

unsigned long CDriverBase::GetElementGlobalIndex(unsigned long iElem) const {
  if (iElem >= GetNumberElements()) {
    SU2_MPI::Error("Element index exceeds size.", CURRENT_FUNCTION);
  }
  return main_geometry->elem[iElem]->GetGlobalIndex();
}

vector<unsigned long> CDriverBase::GetElementNodes(unsigned long iElem) const {
  if (iElem >= GetNumberElements()) {
    SU2_MPI::Error("Element index exceeds size.", CURRENT_FUNCTION);
  }
  const auto nNode = main_geometry->elem[iElem]->GetnNodes();
  vector<unsigned long> values(nNode);

  for (auto iNode = 0u; iNode < nNode; iNode++) {
    values[iNode] = main_geometry->elem[iElem]->GetNode(iNode);
  }
  return values;
}

unsigned long CDriverBase::GetNumberNodes() const { return main_geometry->GetnPoint(); }

unsigned long CDriverBase::GetNumberHaloNodes() const {
  return main_geometry->GetnPoint() - main_geometry->GetnPointDomain();
}

unsigned long CDriverBase::GetNodeGlobalIndex(unsigned long iPoint) const {
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
  }
  return main_geometry->nodes->GetGlobalIndex(iPoint);
}

bool CDriverBase::GetNodeDomain(unsigned long iPoint) const {
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
  }
  return main_geometry->nodes->GetDomain(iPoint);
}

vector<passivedouble> CDriverBase::GetInitialCoordinates(unsigned long iPoint) const {
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
  }
  vector<passivedouble> values(nDim, 0.0);

  if (main_config->GetDeform_Mesh()) {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord(iPoint, iDim);
      values[iDim] = SU2_TYPE::GetValue(value);
    }
  }
  return values;
}

vector<passivedouble> CDriverBase::GetCoordinates(unsigned long iPoint) const {
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
  }
  vector<passivedouble> values;

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = main_geometry->nodes->GetCoord(iPoint, iDim);
    values.push_back(SU2_TYPE::GetValue(value));
  }
  return values;
}

void CDriverBase::SetCoordinates(unsigned long iPoint, vector<passivedouble> values) {
  if (iPoint >= GetNumberNodes()) {
    SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
  }
  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
  }
  for (auto iDim = 0u; iDim < nDim; iDim++) {
    main_geometry->nodes->SetCoord(iPoint, iDim, values[iDim]);
  }
}

unsigned short CDriverBase::GetNumberMarkers() const { return main_config->GetnMarker_All(); }

map<string, unsigned short> CDriverBase::GetMarkerIndices() const {
  const auto nMarker = main_config->GetnMarker_All();
  map<string, unsigned short> indexMap;

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    indexMap[main_config->GetMarker_All_TagBound(iMarker)] = iMarker;
  }
  return indexMap;
}

map<string, string> CDriverBase::GetMarkerTypes() const {
  map<string, string> typeMap;
  string type;

  for (auto iMarker = 0u; iMarker < main_config->GetnMarker_All(); iMarker++) {
    auto tag = main_config->GetMarker_All_TagBound(iMarker);
    auto kindBC = main_config->GetMarker_All_KindBC(iMarker);

    switch (kindBC) {
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
        type = "HEAT_FLUX";
        break;
      case HEAT_TRANSFER:
        type = "HEAT_TRANSFER";
        break;
      case INLET_FLOW:
        type = "INLET_FLOW";
        break;
      case OUTLET_FLOW:
        type = "OUTLET_FLOW";
        break;
      case SUPERSONIC_INLET:
        type = "SUPERSONIC_INLET";
        break;
      case SUPERSONIC_OUTLET:
        type = "SUPERSONIC_OUTLET";
        break;
      case RIEMANN_BOUNDARY:
        type = "RIEMANN";
        break;
      case GILES_BOUNDARY:
        type = "GILES";
        break;
      case DISPLACEMENT_BOUNDARY:
        type = "DISPLACEMENT";
        break;
      case LOAD_BOUNDARY:
        type = "LOAD";
        break;
      case PERIODIC_BOUNDARY:
        type = "PERIODIC";
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
  const auto nMarker = main_config->GetnMarker_All();
  vector<string> tags(nMarker);

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    tags[iMarker] = main_config->GetMarker_All_TagBound(iMarker);
  }
  return tags;
}

vector<string> CDriverBase::GetDeformableMarkerTags() const {
  const auto nMarker = main_config->GetnMarker_Deform_Mesh();
  vector<string> tags(nMarker);

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    tags[iMarker] = main_config->GetMarker_Deform_Mesh_TagBound(iMarker);
  }
  return tags;
}

unsigned long CDriverBase::GetNumberMarkerElements(unsigned short iMarker) const {
  if (iMarker >= GetNumberMarkers()) {
    SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
  }
  return main_geometry->GetnElem_Bound(iMarker);
}

unsigned long CDriverBase::GetMarkerElementGlobalIndex(unsigned short iMarker, unsigned long iElem) const {
  if (iElem >= GetNumberMarkerElements(iMarker)) {
    SU2_MPI::Error("Marker element index exceeds size.", CURRENT_FUNCTION);
  }
  return main_geometry->bound[iMarker][iElem]->GetGlobalIndex();
}

vector<unsigned long> CDriverBase::GetMarkerElementNodes(unsigned short iMarker, unsigned long iElem) const {
  if (iElem >= GetNumberMarkerElements(iMarker)) {
    SU2_MPI::Error("Marker element index exceeds size.", CURRENT_FUNCTION);
  }
  unsigned short nNode = main_geometry->bound[iMarker][iElem]->GetnNodes();
  vector<unsigned long> values(nNode);

  for (auto iNode = 0u; iNode < nNode; iNode++) {
    values[iNode] = main_geometry->bound[iMarker][iElem]->GetNode(iNode);
  }
  return values;
}

unsigned long CDriverBase::GetNumberMarkerNodes(unsigned short iMarker) const {
  if (iMarker >= GetNumberMarkers()) {
    SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
  }
  return main_geometry->GetnVertex(iMarker);
}

unsigned long CDriverBase::GetMarkerNode(unsigned short iMarker, unsigned long iVertex) const {
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }
  return geometry_container[MESH_0][INST_0][ZONE_0]->vertex[iMarker][iVertex]->GetNode();
}

vector<passivedouble> CDriverBase::GetMarkerVertexNormals(unsigned short iMarker, unsigned long iVertex,
                                                          bool normalize) const {
  if (iVertex >= GetNumberMarkerNodes(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }
  const auto* normal = main_geometry->vertex[iMarker][iVertex]->GetNormal();
  const su2double area = normalize ? GeometryToolbox::Norm(nDim, normal) : 1.0;

  vector<passivedouble> values(nDim, 0.0);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    values[iDim] = SU2_TYPE::GetValue(normal[iDim] / area);
  }
  return values;
}

vector<passivedouble> CDriverBase::GetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex) const {
  vector<passivedouble> values(nDim, 0.0);

  if (main_config->GetDeform_Mesh()) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetBound_Disp(iPoint, iDim);
      values[iDim] = SU2_TYPE::GetValue(value);
    }
  }
  return values;
}

void CDriverBase::SetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
  if (!main_config->GetDeform_Mesh()) {
    SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
  }
  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
  }
  const auto iPoint = GetMarkerNode(iMarker, iVertex);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iPoint, iDim, values[iDim]);
  }
}

vector<passivedouble> CDriverBase::GetMarkerVelocities(unsigned short iMarker, unsigned long iVertex) const {
  vector<passivedouble> values(nDim, 0.0);

  if (main_config->GetDeform_Mesh()) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetBound_Vel(iPoint, iDim);
      values[iDim] = SU2_TYPE::GetValue(value);
    }
  }
  return values;
}

void CDriverBase::SetMarkerVelocities(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
  if (!main_config->GetDeform_Mesh()) {
    SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
  }
  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
  }
  const auto iPoint = GetMarkerNode(iMarker, iVertex);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Vel(iPoint, iDim, values[iDim]);
  }
}

void CDriverBase::CommunicateMeshDisplacements(void) {
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(main_geometry, main_config, MESH_DISPLACEMENTS);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(main_geometry, main_config, MESH_DISPLACEMENTS);
}
