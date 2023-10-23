/*!
 * \file CDriverBase.hpp
 * \brief Base class template for all drivers.
 * \author H. Patel, A. Gastaldi
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
#include "../../include/variables/CPrimitiveIndices.hpp"

using namespace std;

CDriverBase::CDriverBase(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator)
    : config_file_name(confFile), StartTime(0.0), StopTime(0.0), UsedTime(0.0), TimeIter(0), nZone(val_nZone) {

  /*--- Some initializations are placed here so that they are also seen by the python wrapper. Note that the python
   * wrapper instantiates a driver directly. ---*/

  /*--- MPI is required to be initialized already, e.g, via SU2_MPI::Init, SU2_MPI::Init_thread, or via mpi4py in the
   * python wrapper. ---*/

  /*--- Initialize MeDiPack ---*/
#ifdef HAVE_MPI
#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
  SU2_MPI::Init_AMPI();
#endif
#endif

  /*--- Set up MPI ---*/
  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- OpenMP initialization ---*/
  omp_initialize();

  /*--- Initialize AD ---*/
  AD::Initialize();
}

CDriverBase::~CDriverBase() = default;

void CDriverBase::InitializeContainers() {
  /*--- Create pointers to all the classes that may be used by drivers. In general, the pointers are instantiated
   * down a hierarchy over all zones, multi-grid levels, equation sets, and equation terms as described in the comments
   * below. ---*/

  config_container = new CConfig*[nZone]();
  output_container = new COutput*[nZone]();
  geometry_container = new CGeometry***[nZone]();
  solver_container = new CSolver****[nZone]();
  numerics_container = new CNumerics*****[nZone]();
  surface_movement = new CSurfaceMovement*[nZone]();
  grid_movement = new CVolumetricMovement**[nZone]();

  nInst = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    nInst[iZone] = 1;
  }
}

void CDriverBase::CommonFinalize() {

  if (numerics_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      delete[] numerics_container[iZone];
    }
    delete[] numerics_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;

  if (solver_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      delete[] solver_container[iZone];
    }
    delete[] solver_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSolver container." << endl;

  if (geometry_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (geometry_container[iZone] != nullptr) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++) {
          delete geometry_container[iZone][iInst][MESH_0];
          delete[] geometry_container[iZone][iInst];
        }
        delete[] geometry_container[iZone];
      }
    }
    delete[] geometry_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  if (surface_movement != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      delete surface_movement[iZone];
    }
    delete[] surface_movement;
  }
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

  if (grid_movement != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (grid_movement[iZone] != nullptr) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++) {
          delete grid_movement[iZone][iInst];
        }
        delete[] grid_movement[iZone];
      }
    }
    delete[] grid_movement;
  }
  if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;

  if (config_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      delete config_container[iZone];
    }
    delete[] config_container;
  }
  delete driver_config;
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  if (output_container != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++) {
      delete output_container[iZone];
    }
    delete[] output_container;
  }
  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  delete[] nInst;
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

vector<string> CDriverBase::GetCHTMarkerTags() const {
  vector<string> tags;
  const auto nMarker = main_config->GetnMarker_All();

  // The CHT markers can be identified as the markers that are customizable with a BC type HEAT_FLUX or ISOTHERMAL.
  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    if ((main_config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX ||
         main_config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) &&
        main_config->GetMarker_All_PyCustom(iMarker)) {
      tags.push_back(main_config->GetMarker_All_TagBound(iMarker));
    }
  }
  return tags;
}

vector<string> CDriverBase::GetInletMarkerTags() const {
  vector<string> tags;
  const auto nMarker = main_config->GetnMarker_All();

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    bool isCustomizable = main_config->GetMarker_All_PyCustom(iMarker);
    bool isInlet = (main_config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);

    if (isCustomizable && isInlet) {
      tags.push_back(main_config->GetMarker_All_TagBound(iMarker));
    }
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
  return main_geometry->vertex[iMarker][iVertex]->GetNode();
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

void CDriverBase::CommunicateMeshDisplacements() {
  solver_container[selected_zone][INST_0][MESH_0][MESH_SOL]->InitiateComms(main_geometry, main_config, MESH_DISPLACEMENTS);
  solver_container[selected_zone][INST_0][MESH_0][MESH_SOL]->CompleteComms(main_geometry, main_config, MESH_DISPLACEMENTS);
}

map<string, unsigned short> CDriverBase::GetSolverIndices() const {
  map<string, unsigned short> indexMap;
  for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
    const auto* solver = solver_container[selected_zone][INST_0][MESH_0][iSol];
    if (solver != nullptr) {
      if (solver->GetSolverName().empty()) SU2_MPI::Error("Solver name was not defined.", CURRENT_FUNCTION);
      indexMap[solver->GetSolverName()] = iSol;
    }
  }
  return indexMap;
}

std::map<string, unsigned short> CDriverBase::GetFEASolutionIndices() const {
  if (solver_container[selected_zone][INST_0][MESH_0][FEA_SOL] == nullptr) {
    SU2_MPI::Error("The FEA solver does not exist.", CURRENT_FUNCTION);
  }
  const auto nDim = main_geometry->GetnDim();
  std::map<string, unsigned short> names;
  names["DISPLACEMENT_X"] = 0;
  names["DISPLACEMENT_Y"] = 1;
  if (nDim == 3) names["DISPLACEMENT_Z"] = 2;

  if (main_config->GetTime_Domain()) {
    names["VELOCITY_X"] = nDim;
    names["VELOCITY_Y"] = nDim + 1;
    if (nDim == 3) names["VELOCITY_Z"] = nDim + 2;
    names["ACCELERATION_X"] = 2 * nDim;
    names["ACCELERATION_Y"] = 2 * nDim + 1;
    if (nDim == 3) names["ACCELERATION_Z"] = 2 * nDim + 2;
  }
  return names;
}

map<string, unsigned short> CDriverBase::GetPrimitiveIndices() const {
  if (solver_container[selected_zone][INST_0][MESH_0][FLOW_SOL] == nullptr) {
    SU2_MPI::Error("The flow solver does not exist.", CURRENT_FUNCTION);
  }
  return PrimitiveNameToIndexMap(CPrimitiveIndices<unsigned short>(
      main_config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE,
      main_config->GetNEMOProblem(), nDim, main_config->GetnSpecies()));
}
