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
    FFDBox[iZone] = new CFreeFormDefBox*[MAX_NUMBER_FFD];
  }

  driver_config = nullptr;
  driver_output = nullptr;

  main_config = nullptr;
  main_geometry = nullptr;
}

unsigned short CDriverBase::GetNumberDesignVariables() const { return main_config->GetnDV(); }

unsigned short CDriverBase::GetNumberFFDBoxes() const { return main_config->GetnFFDBox(); }

unsigned short CDriverBase::GetNumberFFDBoxCornerPoints(unsigned short iFFDBox) const {
  return FFDBox[ZONE_0][iFFDBox]->GetnCornerPoints();
}

unsigned short CDriverBase::GetNumberFFDBoxControlPoints(unsigned short iFFDBox) const {
  return FFDBox[ZONE_0][iFFDBox]->GetnControlPoints();
}

unsigned long CDriverBase::GetNumberFFDBoxSurfacePoints(unsigned short iFFDBox) const {
  return FFDBox[ZONE_0][iFFDBox]->GetnSurfacePoint();
}

vector<unsigned short> CDriverBase::GetFFDBoxMarkerIDs(unsigned short iFFDBox) const {
  vector<unsigned short> values;

  for (auto iPoint = 0ul; iPoint < GetNumberFFDBoxSurfacePoints(iFFDBox); iPoint++) {
    values.push_back(FFDBox[ZONE_0][iFFDBox]->Get_MarkerIndex(iPoint));
  }
  return values;
}

vector<unsigned long> CDriverBase::GetFFDBoxVertexIDs(unsigned short iFFDBox) const {
  vector<unsigned long> values;

  for (auto iPoint = 0ul; iPoint < GetNumberFFDBoxSurfacePoints(iFFDBox); iPoint++) {
    values.push_back(FFDBox[ZONE_0][iFFDBox]->Get_VertexIndex(iPoint));
  }
  return values;
}

vector<unsigned long> CDriverBase::GetFFDBoxPointIDs(unsigned short iFFDBox) const {
  vector<unsigned long> values;

  for (auto iPoint = 0ul; iPoint < GetNumberFFDBoxSurfacePoints(iFFDBox); iPoint++) {
    values.push_back(FFDBox[ZONE_0][iFFDBox]->Get_PointIndex(iPoint));
  }
  return values;
}

vector<vector<passivedouble>> CDriverBase::GetFFDBoxCornerCoordinates(unsigned short iFFDBox) const {
  vector<vector<passivedouble>> values;

  for (auto iPoint = 0u; iPoint < GetNumberFFDBoxCornerPoints(iFFDBox); iPoint++) {
    vector<passivedouble> coords;

    const su2double* coord = FFDBox[ZONE_0][iFFDBox]->GetCoordCornerPoints(iPoint);

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      coords.push_back(SU2_TYPE::GetValue(coord[iDim]));
    }

    values.push_back(coords);
  }

  return values;
}

vector<vector<passivedouble>> CDriverBase::GetFFDBoxControlPointCoordinates(unsigned short iFFDBox) const {
  vector<vector<passivedouble>> values;

  for (auto iOrder = 0u; iOrder < FFDBox[ZONE_0][iFFDBox]->GetlOrder(); iOrder++) {
    for (auto jOrder = 0u; jOrder < FFDBox[ZONE_0][iFFDBox]->GetmOrder(); jOrder++) {
      for (auto kOrder = 0u; kOrder < FFDBox[ZONE_0][iFFDBox]->GetnOrder(); kOrder++) {
        values.push_back(GetFFDBoxControlPointCoordinates(iFFDBox, iOrder, jOrder, kOrder));
      }
    }
  }

  return values;
}

vector<passivedouble> CDriverBase::GetFFDBoxControlPointCoordinates(unsigned short iFFDBox, unsigned int iOrder,
                                                                    unsigned int jOrder, unsigned int kOrder) const {
  vector<passivedouble> values;

  const su2double* coord = FFDBox[ZONE_0][iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    values.push_back(SU2_TYPE::GetValue(coord[iDim]));
  }

  return values;
}

vector<vector<passivedouble>> CDriverBase::GetFFDBoxSurfaceCoordinates(unsigned short iFFDBox, bool parametric) const {
  vector<vector<passivedouble>> values;

  for (auto iPoint = 0ul; iPoint < GetNumberFFDBoxSurfacePoints(iFFDBox); iPoint++) {
    values.push_back(GetFFDBoxSurfaceCoordinates(iFFDBox, iPoint, parametric));
  }

  return values;
}

vector<passivedouble> CDriverBase::GetFFDBoxSurfaceCoordinates(unsigned short iFFDBox, unsigned long iPoint,
                                                               bool parametric) const {
  vector<passivedouble> values;
  const su2double* coord;

  if (parametric) {
    coord = FFDBox[ZONE_0][iFFDBox]->Get_ParametricCoord(iPoint);
  } else {
    coord = FFDBox[ZONE_0][iFFDBox]->Get_CartesianCoord(iPoint);
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    values.push_back(SU2_TYPE::GetValue(coord[iDim]));
  }

  return values;
}

unsigned short CDriverBase::GetNumberMarkers() const { return main_config->GetnMarker_All(); }

map<string, unsigned short> CDriverBase::GetMarkerIndices() const {
  const auto nMarker = main_config->GetnMarker_All();
  map<string, unsigned short> indexMap;

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    auto tag = main_config->GetMarker_All_TagBound(iMarker);

    indexMap[tag] = iMarker;
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
  const auto nMarker = main_config->GetnMarker_All();

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    tags.push_back(main_config->GetMarker_All_TagBound(iMarker));
  }

  return tags;
}

vector<string> CDriverBase::GetDeformableMarkerTags() const {
  vector<string> tags;
  const auto nMarker = main_config->GetnMarker_Deform_Mesh();

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    tags.push_back(main_config->GetMarker_Deform_Mesh_TagBound(iMarker));
  }

  return tags;
}

unsigned long CDriverBase::GetNumberDimensions() const { return main_geometry->GetnDim(); }

unsigned long CDriverBase::GetNumberElements() const { return main_geometry->GetnElem(); }

unsigned long CDriverBase::GetNumberMarkerElements(unsigned short iMarker) const {
  if (iMarker >= GetNumberMarkers()) {
    SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
  }

  return main_geometry->GetnElem_Bound(iMarker);
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

  return main_geometry->elem[iElem]->GetGlobalIndex();
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

  return main_geometry->bound[iMarker][iBound]->GetGlobalIndex();
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

  unsigned short nNode = main_geometry->elem[iElem]->GetnNodes();

  vector<unsigned long> values;

  for (auto iNode = 0u; iNode < nNode; iNode++) {
    unsigned long iPoint = main_geometry->elem[iElem]->GetNode(iNode);

    values.push_back(main_geometry->nodes->GetGlobalIndex(iPoint));
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

  unsigned short nNode = main_geometry->bound[iMarker][iBound]->GetnNodes();

  vector<unsigned long> values;

  for (auto iNode = 0u; iNode < nNode; iNode++) {
    unsigned long iPoint = main_geometry->bound[iMarker][iBound]->GetNode(iNode);

    values.push_back(main_geometry->nodes->GetGlobalIndex(iPoint));
  }

  return values;
}

unsigned long CDriverBase::GetNumberVertices() const { return main_geometry->GetnPoint(); }

unsigned long CDriverBase::GetNumberMarkerVertices(unsigned short iMarker) const {
  if (iMarker >= GetNumberMarkers()) {
    SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
  }

  return main_geometry->GetnVertex(iMarker);
}

unsigned long CDriverBase::GetNumberHaloVertices() const {
  const auto nPoint = GetNumberVertices();
  unsigned long nHalo = 0;

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    if (!(main_geometry->nodes->GetDomain(iPoint))) {
      nHalo += 1;
    }
  }

  return nHalo;
}

unsigned long CDriverBase::GetNumberMarkerHaloVertices(unsigned short iMarker) const {
  const auto nVertex = GetNumberMarkerVertices(iMarker);
  unsigned long nHalo = 0;

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

    if (!(main_geometry->nodes->GetDomain(iPoint))) {
      nHalo += 1;
    }
  }

  return nHalo;
}

vector<unsigned long> CDriverBase::GetMarkerVertexIndices(unsigned short iMarker) const {
  const auto nVertex = GetNumberMarkerVertices(iMarker);

  vector<unsigned long> values;

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    values.push_back(GetMarkerVertexIndices(iMarker, iVertex));
  }

  return values;
}

unsigned long CDriverBase::GetMarkerVertexIndices(unsigned short iMarker, unsigned long iVertex) const {
  if (iVertex >= GetNumberMarkerVertices(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
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
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  return main_geometry->nodes->GetGlobalIndex(iPoint);
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
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  return main_geometry->nodes->GetGlobalIndex(iPoint);
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
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  return main_geometry->nodes->GetDomain(iPoint);
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
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  return main_geometry->nodes->GetDomain(iPoint);
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
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  if (!main_config->GetDeform_Mesh()) {
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
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  vector<passivedouble> values(nDim, 0.0);

  if (!main_config->GetDeform_Mesh()) {
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
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values;

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = main_geometry->nodes->GetCoord(iPoint, iDim);

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
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  vector<passivedouble> values;

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = main_geometry->nodes->GetCoord(iPoint, iDim);

    values.push_back(SU2_TYPE::GetValue(value));
  }

  return values;
}

void CDriverBase::SetCoordinates(vector<vector<passivedouble>> values) {
  const auto nPoint = GetNumberVertices();

  if (values.size() != nPoint) {
    SU2_MPI::Error("Invalid number of vertices!", CURRENT_FUNCTION);
  }

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    SetCoordinates(iPoint, values[iPoint]);
  }
}

void CDriverBase::SetCoordinates(unsigned long iPoint, vector<passivedouble> values) {
  if (iPoint >= GetNumberVertices()) {
    SU2_MPI::Error("Vertex index exceeds mesh size.", CURRENT_FUNCTION);
  }
  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    main_geometry->nodes->SetCoord(iPoint, iDim, values[iDim]);
  }
}

void CDriverBase::SetMarkerCoordinates(unsigned short iMarker, vector<vector<passivedouble>> values) {
  const auto nVertex = GetNumberMarkerVertices(iMarker);

  if (values.size() != nVertex) {
    SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
  }

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    SetMarkerCoordinates(iMarker, iVertex, values[iVertex]);
  }
}

void CDriverBase::SetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    main_geometry->nodes->SetCoord(iPoint, iDim, values[iDim]);
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
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  vector<passivedouble> values(nDim, 0.0);

  if (!main_config->GetDeform_Mesh()) {
    return values;
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetBound_Disp(iPoint, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

void CDriverBase::SetMarkerDisplacements(unsigned short iMarker, vector<vector<passivedouble>> values) {
  if (!main_config->GetDeform_Mesh()) {
    SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
  }

  const auto nVertex = GetNumberMarkerVertices(iMarker);

  if (values.size() != nVertex) {
    SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
  }

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    SetMarkerDisplacements(iMarker, iVertex, values[iVertex]);
  }
}

void CDriverBase::SetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
  if (!main_config->GetDeform_Mesh()) {
    SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
  }

  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
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
  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  vector<passivedouble> values(nDim, 0.0);

  if (!main_config->GetDeform_Mesh()) {
    return values;
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double value = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetBound_Vel(iPoint, iDim);

    values[iDim] = SU2_TYPE::GetValue(value);
  }

  return values;
}

void CDriverBase::SetMarkerVelocities(unsigned short iMarker, vector<vector<passivedouble>> values) {
  if (!main_config->GetDeform_Mesh()) {
    SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
  }

  const auto nVertex = GetNumberMarkerVertices(iMarker);

  if (values.size() != nVertex) {
    SU2_MPI::Error("Invalid number of marker vertices!", CURRENT_FUNCTION);
  }

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    SetMarkerVelocities(iMarker, iVertex, values[iVertex]);
  }
}

void CDriverBase::SetMarkerVelocities(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
  if (!main_config->GetDeform_Mesh()) {
    SU2_MPI::Error("Mesh solver is not defined!", CURRENT_FUNCTION);
  }

  auto iPoint = GetMarkerVertexIndices(iMarker, iVertex);

  if (values.size() != nDim) {
    SU2_MPI::Error("Invalid number of dimensions!", CURRENT_FUNCTION);
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iPoint, iDim, values[iDim]);
  }
}

vector<vector<passivedouble>> CDriverBase::GetMarkerVertexNormals(unsigned short iMarker, bool normalize) const {
  const auto nVertex = GetNumberMarkerVertices(iMarker);

  vector<vector<passivedouble>> values;

  for (auto iVertex = 0ul; iVertex < nVertex; iVertex++) {
    values.push_back(GetMarkerVertexNormals(iMarker, iVertex, normalize = normalize));
  }

  return values;
}

vector<passivedouble> CDriverBase::GetMarkerVertexNormals(unsigned short iMarker, unsigned long iVertex,
                                                          bool normalize) const {
  if (iVertex >= GetNumberMarkerVertices(iMarker)) {
    SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
  }

  vector<passivedouble> values(nDim, 0.0);

  auto normal = main_geometry->vertex[iMarker][iVertex]->GetNormal();
  auto area = GeometryToolbox::Norm(nDim, normal);

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    if (normalize) {
      values[iDim] = SU2_TYPE::GetValue(normal[iDim] / area);
    } else {
      values[iDim] = SU2_TYPE::GetValue(normal[iDim]);
    }
  }

  return values;
}

void CDriverBase::CommunicateMeshDisplacements(void) {
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(main_geometry, main_config, MESH_DISPLACEMENTS);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(main_geometry, main_config, MESH_DISPLACEMENTS);
}

void CDriverBase::ReadFFDInfo(CGeometry* geometry, CConfig* config, CFreeFormDefBox** FFDBox, string val_mesh_filename) {
  string text_line, iTag;
  ifstream mesh_file;
  su2double CPcoord[3], coord[] = {0, 0, 0};
  unsigned short degree[3], iFFDBox, iCornerPoints, iControlPoints, iMarker, iDegree, jDegree, kDegree, iChar,
      LevelFFDBox, nParentFFDBox, iParentFFDBox, nChildFFDBox, iChildFFDBox, nMarker, *nCornerPoints, *nControlPoints;
  unsigned long iSurfacePoints, iPoint, jPoint, iVertex, nVertex, nPoint, iElem = 0, nElem, my_nSurfPoints, nSurfPoints,
                                                                          *nSurfacePoints;
  su2double XCoord, YCoord;

  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim(), iDim;
  unsigned short SplineOrder[3];
  unsigned short Blending = 0;

  unsigned short nFFDBox = 0;
  unsigned short nLevel = 0;
  bool FFDBoxDefinition = false;

  mesh_file.open(val_mesh_filename);
  if (mesh_file.fail()) {
    SU2_MPI::Error("There is no geometry file (ReadFFDInfo)!!", CURRENT_FUNCTION);
  }

  while (getline(mesh_file, text_line)) {
    /*--- Read the inner elements. ---*/

    string::size_type position = text_line.find("NELEM=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      nElem = atoi(text_line.c_str());
      for (iElem = 0; iElem < nElem; iElem++) {
        getline(mesh_file, text_line);
      }
    }

    /*--- Read the inner points. ---*/

    position = text_line.find("NPOIN=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      nPoint = atoi(text_line.c_str());
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        getline(mesh_file, text_line);
      }
    }

    /*--- Read the boundaries. ---*/

    position = text_line.find("NMARK=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      nMarker = atoi(text_line.c_str());
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        getline(mesh_file, text_line);
        getline(mesh_file, text_line);
        text_line.erase(0, 13);
        nVertex = atoi(text_line.c_str());
        for (iVertex = 0; iVertex < nVertex; iVertex++) {
          getline(mesh_file, text_line);
        }
      }
    }

    /*--- Read the FFDBox information. ---*/

    position = text_line.find("FFD_NBOX=", 0);
    if (position != string::npos) {
      text_line.erase(0, 9);
      nFFDBox = atoi(text_line.c_str());

      if (rank == MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;

      nCornerPoints = new unsigned short[nFFDBox];
      nControlPoints = new unsigned short[nFFDBox];
      nSurfacePoints = new unsigned long[nFFDBox];

      getline(mesh_file, text_line);
      text_line.erase(0, 11);
      nLevel = atoi(text_line.c_str());

      if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;

      for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
        /*--- Read the name of the FFD box. ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 8);

        /*--- Remove extra data from the FFDBox name. ---*/

        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find(" ", 0);
          if (position != string::npos) text_line.erase(position, 1);
          position = text_line.find("\r", 0);
          if (position != string::npos) text_line.erase(position, 1);
          position = text_line.find("\n", 0);
          if (position != string::npos) text_line.erase(position, 1);
        }

        string TagFFDBox = text_line.c_str();

        if (rank == MASTER_NODE) cout << "FFD box tag: " << TagFFDBox << ". ";

        /*--- Read the level of the FFD box. ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 10);
        LevelFFDBox = atoi(text_line.c_str());

        if (rank == MASTER_NODE) cout << "FFD box level: " << LevelFFDBox << ". ";

        /*--- Read the degree of the FFD box. ---*/

        if (nDim == 2) {
          if (polar) {
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[0] = atoi(text_line.c_str());
            degree[1] = 1;
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[2] = atoi(text_line.c_str());
          } else {
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[0] = atoi(text_line.c_str());
            getline(mesh_file, text_line);
            text_line.erase(0, 13);
            degree[1] = atoi(text_line.c_str());
            degree[2] = 1;
          }
        } else {
          getline(mesh_file, text_line);
          text_line.erase(0, 13);
          degree[0] = atoi(text_line.c_str());
          getline(mesh_file, text_line);
          text_line.erase(0, 13);
          degree[1] = atoi(text_line.c_str());
          getline(mesh_file, text_line);
          text_line.erase(0, 13);
          degree[2] = atoi(text_line.c_str());
        }

        if (rank == MASTER_NODE) {
          if (nDim == 2) {
            if (polar)
              cout << "Degrees: " << degree[0] << ", " << degree[2] << "." << endl;
            else
              cout << "Degrees: " << degree[0] << ", " << degree[1] << "." << endl;
          } else
            cout << "Degrees: " << degree[0] << ", " << degree[1] << ", " << degree[2] << "." << endl;
        }

        getline(mesh_file, text_line);
        if (text_line.substr(0, 12) != "FFD_BLENDING") {
          SU2_MPI::Error(
              string("Deprecated FFD information found in mesh file.\n") +
                  string(
                      "FFD information generated with SU2 version <= 4.3 is incompatible with the current version.") +
                  string("Run SU2_DEF again with DV_KIND= FFD_SETTING."),
              CURRENT_FUNCTION);
        }
        text_line.erase(0, 14);
        if (text_line == "BEZIER") {
          Blending = BEZIER;
        }
        if (text_line == "BSPLINE_UNIFORM") {
          Blending = BSPLINE_UNIFORM;
        }

        if (Blending == BSPLINE_UNIFORM) {
          getline(mesh_file, text_line);
          text_line.erase(0, 17);
          SplineOrder[0] = atoi(text_line.c_str());
          getline(mesh_file, text_line);
          text_line.erase(0, 17);
          SplineOrder[1] = atoi(text_line.c_str());
          if (nDim == 3) {
            getline(mesh_file, text_line);
            text_line.erase(0, 17);
            SplineOrder[2] = atoi(text_line.c_str());
          } else {
            SplineOrder[2] = 2;
          }
        }
        if (rank == MASTER_NODE) {
          if (Blending == BSPLINE_UNIFORM) {
            cout << "FFD Blending using B-Splines. ";
            cout << "Order: " << SplineOrder[0] << ", " << SplineOrder[1];
            if (nDim == 3) cout << ", " << SplineOrder[2];
            cout << ". " << endl;
          }
          if (Blending == BEZIER) {
            cout << "FFD Blending using Bezier Curves." << endl;
          }
        }

        FFDBox[iFFDBox] = new CFreeFormDefBox(degree, SplineOrder, Blending);
        FFDBox[iFFDBox]->SetTag(TagFFDBox);
        FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

        /*--- Read the number of parents boxes. ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 12);
        nParentFFDBox = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox << ". ";
        for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
          getline(mesh_file, text_line);

          /*--- Remove extra data from the FFDBox name. ---*/

          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find(" ", 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find("\r", 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find("\n", 0);
            if (position != string::npos) text_line.erase(position, 1);
          }

          string ParentFFDBox = text_line.c_str();
          FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
        }

        /*--- Read the number of children boxes. ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 13);
        nChildFFDBox = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox << "." << endl;

        for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
          getline(mesh_file, text_line);

          /*--- Remove extra data from the FFDBox name. ---*/

          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find(" ", 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find("\r", 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find("\n", 0);
            if (position != string::npos) text_line.erase(position, 1);
          }

          string ChildFFDBox = text_line.c_str();
          FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
        }

        /*--- Read the number of the corner points. ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 18);
        nCornerPoints[iFFDBox] = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Corner points: " << nCornerPoints[iFFDBox] << ". ";
        if (nDim == 2) nCornerPoints[iFFDBox] = nCornerPoints[iFFDBox] * SU2_TYPE::Int(2);

        /*--- Read the coordinates of the corner points. ---*/

        if (nDim == 2) {
          if (polar) {
            getline(mesh_file, text_line);
            istringstream FFDBox_line_1(text_line);
            FFDBox_line_1 >> XCoord;
            FFDBox_line_1 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 4);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 7);

            getline(mesh_file, text_line);
            istringstream FFDBox_line_2(text_line);
            FFDBox_line_2 >> XCoord;
            FFDBox_line_2 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 0);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 3);

            getline(mesh_file, text_line);
            istringstream FFDBox_line_3(text_line);
            FFDBox_line_3 >> XCoord;
            FFDBox_line_3 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 1);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 2);

            getline(mesh_file, text_line);
            istringstream FFDBox_line_4(text_line);
            FFDBox_line_4 >> XCoord;
            FFDBox_line_4 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = -sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 5);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1) * YCoord;
            CPcoord[2] = sin(0.1) * YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 6);

          } else {
            for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
              if (iCornerPoints < nCornerPoints[iFFDBox] / SU2_TYPE::Int(2)) {
                getline(mesh_file, text_line);
                istringstream FFDBox_line(text_line);
                FFDBox_line >> CPcoord[0];
                FFDBox_line >> CPcoord[1];
                CPcoord[2] = -0.5;
              } else {
                CPcoord[0] =
                    FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints - nCornerPoints[iFFDBox] / SU2_TYPE::Int(2));
                CPcoord[1] =
                    FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints - nCornerPoints[iFFDBox] / SU2_TYPE::Int(2));
                CPcoord[2] = 0.5;
              }
              FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, iCornerPoints);
            }
          }

        } else {
          for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
            getline(mesh_file, text_line);
            istringstream FFDBox_line(text_line);
            FFDBox_line >> CPcoord[0];
            FFDBox_line >> CPcoord[1];
            FFDBox_line >> CPcoord[2];
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, iCornerPoints);
          }
        }

        /*--- Read the number of the control points. ---*/

        getline(mesh_file, text_line);
        text_line.erase(0, 19);
        nControlPoints[iFFDBox] = atoi(text_line.c_str());

        if (rank == MASTER_NODE) cout << "Control points: " << nControlPoints[iFFDBox] << ". ";

        /*--- Method to identify if there is a FFDBox definition. ---*/

        if (nControlPoints[iFFDBox] != 0) FFDBoxDefinition = true;

        /*--- Read the coordinates of the control points. ---*/

        for (iControlPoints = 0; iControlPoints < nControlPoints[iFFDBox]; iControlPoints++) {
          getline(mesh_file, text_line);
          istringstream FFDBox_line(text_line);
          FFDBox_line >> iDegree;
          FFDBox_line >> jDegree;
          FFDBox_line >> kDegree;
          FFDBox_line >> CPcoord[0];
          FFDBox_line >> CPcoord[1];
          FFDBox_line >> CPcoord[2];
          FFDBox[iFFDBox]->SetCoordControlPoints(CPcoord, iDegree, jDegree, kDegree);
          FFDBox[iFFDBox]->SetCoordControlPoints_Copy(CPcoord, iDegree, jDegree, kDegree);
        }

        getline(mesh_file, text_line);
        text_line.erase(0, 19);
        nSurfacePoints[iFFDBox] = atoi(text_line.c_str());

        /*--- The surface points parametric coordinates (all the nodes read the FFD
         * information but they only store their part). ---*/

        my_nSurfPoints = 0;
        for (iSurfacePoints = 0; iSurfacePoints < nSurfacePoints[iFFDBox]; iSurfacePoints++) {
          getline(mesh_file, text_line);
          istringstream FFDBox_line(text_line);
          FFDBox_line >> iTag;
          FFDBox_line >> iPoint;

          if (config->GetMarker_All_TagBound(iTag) != -1) {
            iMarker = config->GetMarker_All_TagBound(iTag);
            FFDBox_line >> CPcoord[0];
            FFDBox_line >> CPcoord[1];
            FFDBox_line >> CPcoord[2];

            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if (iPoint == geometry->nodes->GetGlobalIndex(jPoint)) {
                for (iDim = 0; iDim < nDim; iDim++) {
                  coord[iDim] = geometry->nodes->GetCoord(jPoint, iDim);
                }
                FFDBox[iFFDBox]->Set_MarkerIndex(iMarker);
                FFDBox[iFFDBox]->Set_VertexIndex(iVertex);
                FFDBox[iFFDBox]->Set_PointIndex(jPoint);
                FFDBox[iFFDBox]->Set_ParametricCoord(CPcoord);
                FFDBox[iFFDBox]->Set_CartesianCoord(coord);
                my_nSurfPoints++;
              }
            }
          }
        }

        nSurfacePoints[iFFDBox] = my_nSurfPoints;

#ifdef HAVE_MPI
        nSurfPoints = 0;
        SU2_MPI::Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
        if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints << "." << endl;
#else
        nSurfPoints = my_nSurfPoints;
        if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints << "." << endl;
#endif
      }

      delete[] nCornerPoints;
      delete[] nControlPoints;
      delete[] nSurfacePoints;
    }
  }
  mesh_file.close();

  if (nFFDBox == 0) {
    if (rank == MASTER_NODE) cout << "There is no FFD box definition. Just in case, check the .su2 file" << endl;
  }
}
