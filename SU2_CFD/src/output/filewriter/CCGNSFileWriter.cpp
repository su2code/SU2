/*!
 * \file CCGNSFileWriter.cpp
 * \brief Filewriter class for CGNS format.
 * \author G. Baldan
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/output/filewriter/CCGNSFileWriter.hpp"
#include "../../../include/output/filewriter/CFVMDataSorter.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"

#include <numeric>

const string CCGNSFileWriter::fileExt = ".cgns";

CCGNSFileWriter::CCGNSFileWriter(CParallelDataSorter* valDataSorter, bool isSurf, bool doublePrecision)
    : CFileWriter(valDataSorter, fileExt), isSurface(isSurf), doublePrecisionFields(doublePrecision) {}

void CCGNSFileWriter::WriteData(string val_filename) {

#ifdef HAVE_CGNS

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);
  /*--- Open the CGNS file for writing.  ---*/
  InitializeMeshFile(val_filename);

  if (surfaceMarkers.empty()) {
    WriteZone("Zone");
  } else {
    /*--- One zone per marker, the surface data is sorted again for each of them. ---*/
    for (const auto& marker : surfaceMarkers) {
      dataSorter->SortConnectivity(config, geometry, vector<string>{marker});
      dataSorter->SortOutputData();
      WriteZone(marker);
    }
  }

  /*--- Close the CGNS file. ---*/
  if (rank == MASTER_NODE) CallCGNS(cg_close(cgnsFileID));

#endif
}

void CCGNSFileWriter::SetBoundaryMarkers(CConfig* valConfig, CGeometry* valGeometry,
                                         const CFVMDataSorter* volumeSorter) {
  boundaryMarkers.clear();

  for (unsigned short iMarkerCfg = 0; iMarkerCfg < valConfig->GetnMarker_CfgFile(); iMarkerCfg++) {
    const string tag = valConfig->GetMarker_CfgFile_TagBound(iMarkerCfg);
    const auto kindBC = valConfig->GetMarker_CfgFile_KindBC(tag);
    if (kindBC == SEND_RECEIVE) continue;

    BoundaryMarker marker{tag, kindBC, {}};

    for (unsigned short iMarker = 0; iMarker < valConfig->GetnMarker_All(); iMarker++) {
      if (valConfig->GetMarker_All_TagBound(iMarker) != tag) continue;

      for (unsigned long iElem = 0; iElem < valGeometry->GetnElem_Bound(iMarker); iElem++) {
        const auto* elem = valGeometry->bound[iMarker][iElem];

        /*--- Same rule as for the volume elements: keep the element on the rank where none of its nodes is a halo.
         Also require a node owned by this rank, so that an element whose nodes are all owned by lower ranks is not
         kept again by a higher rank that holds it as a halo element. ---*/
        bool halo = false, owned = false;
        for (unsigned short iNode = 0; iNode < elem->GetnNodes(); iNode++) {
          halo |= volumeSorter->GetHalo(elem->GetNode(iNode));
          owned |= valGeometry->nodes->GetDomain(elem->GetNode(iNode));
        }
        if (halo || !owned) continue;

        marker.conn.push_back(elem->GetVTK_Type());
        for (unsigned short iNode = 0; iNode < elem->GetnNodes(); iNode++)
          marker.conn.push_back(valGeometry->nodes->GetGlobalIndex(elem->GetNode(iNode)) + 1);
      }
    }
    boundaryMarkers.push_back(std::move(marker));
  }
}

void CCGNSFileWriter::SetSurfaceMarkers(CConfig* valConfig, CGeometry* valGeometry) {
  config = valConfig;
  geometry = valGeometry;
  surfaceMarkers.clear();

  for (unsigned short iMarkerCfg = 0; iMarkerCfg < valConfig->GetnMarker_CfgFile(); iMarkerCfg++) {
    const string tag = valConfig->GetMarker_CfgFile_TagBound(iMarkerCfg);
    if (valConfig->GetMarker_CfgFile_Plotting(tag) != YES) continue;

    /*--- Only keep the markers present on at least one rank. ---*/
    int localFound = 0, globalFound = 0;
    for (unsigned short iMarker = 0; iMarker < valConfig->GetnMarker_All(); iMarker++) {
      if (valConfig->GetMarker_All_TagBound(iMarker) == tag && valConfig->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE)
        localFound = 1;
    }
    SU2_MPI::Allreduce(&localFound, &globalFound, 1, MPI_INT, MPI_SUM, SU2_MPI::GetComm());
    if (globalFound > 0) surfaceMarkers.push_back(tag);
  }
}

#ifdef HAVE_CGNS
void CCGNSFileWriter::WriteZone(const string& zoneName) {
  InitializeZone(zoneName);

  /*--- Write point coordinates. ---*/
  WriteField(0, "CoordinateX");
  WriteField(1, "CoordinateY");
  if (nDim == 3) WriteField(2, "CoordinateZ");

  /*--- Write mesh connectivity. ---*/
  if (nDim == 2) {
    WriteConnectivity(LINE, "Lines");
    WriteConnectivity(TRIANGLE, "Triangles");
    WriteConnectivity(QUADRILATERAL, "Quadrilaterals");
  }
  if (nDim == 3) {
    WriteConnectivity(TRIANGLE, "Triangles");
    WriteConnectivity(QUADRILATERAL, "Quadrilaterals");
    WriteConnectivity(TETRAHEDRON, "Tetrahedra");
    WriteConnectivity(PYRAMID, "Pyramids");
    WriteConnectivity(PRISM, "Prisms");
    WriteConnectivity(HEXAHEDRON, "Hexahedra");
  }

  /*--- Write the boundaries of a volume file. ---*/
  if (!isSurface) WriteBoundaries();

  /*--- Initialize and write fields. ---*/
  InitializeFields();

  const auto& fieldNames = dataSorter->GetFieldNames();
  for (unsigned long i = nDim; i < fieldNames.size(); ++i) {
    WriteField(i, fieldNames[i]);
  }
}

void CCGNSFileWriter::InitializeMeshFile(const string& val_filename) {
  nDim = dataSorter->GetnDim();

  /*--- If surface file cell dimension is decreased. ---*/
  const auto nCell = static_cast<int>(nDim - isSurface);

  if (rank == MASTER_NODE) {
    /*--- Remove the previous file if present. ---*/
    remove(val_filename.c_str());

    /*--- Create CGNS file and open in write mode. ---*/
    CallCGNS(cg_open(val_filename.c_str(), CG_MODE_WRITE, &cgnsFileID));

    /*--- Create Base. ---*/
    CallCGNS(cg_base_write(cgnsFileID, "Base", nCell, nDim, &cgnsBase));
  }
}

void CCGNSFileWriter::InitializeZone(const string& zoneName) {
  if (!dataSorter->GetConnectivitySorted()) {
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  nLocalPoints = dataSorter->GetnPoints();
  GlobalElem = static_cast<cgsize_t>(dataSorter->GetnElemGlobal());
  GlobalPoint = static_cast<cgsize_t>(dataSorter->GetnPointsGlobal());
  cumulative = 0;

  if (rank == MASTER_NODE) {
    /*--- Create Zone. The number of cells does not include the boundary elements. ---*/
    array<cgsize_t, 3> zoneData;

    zoneData[0] = GlobalPoint;
    zoneData[1] = GlobalElem;
    zoneData[2] = 0;

    CallCGNS(cg_zone_write(cgnsFileID, cgnsBase, zoneName.substr(0, 32).c_str(), zoneData.data(), Unstructured,
                           &cgnsZone));
  }
}

void CCGNSFileWriter::WriteBoundaries() {
  for (const auto& marker : boundaryMarkers) {
    /*--- Gather the boundary elements of this marker on the master node, in rank order. ---*/
    const unsigned long localSize = marker.conn.size();
    vector<unsigned long> sizes(size);
    SU2_MPI::Allgather(&localSize, 1, MPI_UNSIGNED_LONG, sizes.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

    const auto totalSize = std::accumulate(sizes.begin(), sizes.end(), 0ul);
    if (totalSize == 0) continue;

    if (rank != MASTER_NODE) {
      SendChunked(marker.conn.data(), localSize * sizeof(unsigned long), MASTER_NODE, 2);
      continue;
    }

    vector<unsigned long> conn(totalSize);
    std::copy(marker.conn.begin(), marker.conn.end(), conn.begin());
    auto offset = localSize;
    for (int i = 0; i < size; ++i) {
      if (i == MASTER_NODE) continue;
      RecvChunked(conn.data() + offset, sizes[i] * sizeof(unsigned long), i, 2);
      offset += sizes[i];
    }

    /*--- Convert to the CGNS numbering: with a single element type the node ids of the elements, otherwise a MIXED
     section, i.e. the CGNS element type followed by the node ids of each element, plus the start offsets. ---*/
    vector<cgsize_t> elems, mixed, startOffsets{0};
    bool singleType = true;
    for (size_t pos = 0; pos < conn.size();) {
      const auto type = static_cast<unsigned short>(conn[pos]);
      const auto nNodes = nPointsOfElementType(type);
      singleType &= (type == conn[0]);
      mixed.push_back(GetCGNSType(type));
      for (unsigned short iNode = 1; iNode <= nNodes; ++iNode) {
        elems.push_back(static_cast<cgsize_t>(conn[pos + iNode]));
        mixed.push_back(static_cast<cgsize_t>(conn[pos + iNode]));
      }
      startOffsets.push_back(static_cast<cgsize_t>(mixed.size()));
      pos += nNodes + 1;
    }
    const auto nElem = static_cast<cgsize_t>(startOffsets.size() - 1);

    const string name = marker.name.substr(0, 32);
    cgsize_t range[2] = {cumulative + 1, cumulative + nElem};
    int section;
    if (singleType) {
      CallCGNS(cg_section_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), GetCGNSType(conn[0]), range[0],
                                range[1], 0, elems.data(), &section));
    } else {
      CallCGNS(cg_poly_section_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), MIXED, range[0], range[1], 0,
                                     mixed.data(), startOffsets.data(), &section));
    }
    cumulative += nElem;

    /*--- The BC points to the boundary elements and takes its type from a family with the name of the marker. ---*/
    int bc, family, familyBC;
    CallCGNS(cg_boco_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), FamilySpecified, PointRange, 2, range, &bc));
    CallCGNS(cg_boco_gridlocation_write(cgnsFileID, cgnsBase, cgnsZone, bc, nDim == 3 ? FaceCenter : EdgeCenter));
    CallCGNS(cg_goto(cgnsFileID, cgnsBase, "Zone_t", cgnsZone, "ZoneBC_t", 1, "BC_t", bc, "end"));
    CallCGNS(cg_famname_write(name.c_str()));

    CallCGNS(cg_family_write(cgnsFileID, cgnsBase, name.c_str(), &family));
    CallCGNS(cg_fambc_write(cgnsFileID, cgnsBase, family, "FamBC", GetCGNSBCType(marker.kindBC), &familyBC));
  }
}

BCType_t CCGNSFileWriter::GetCGNSBCType(unsigned short kindBC) {
  switch (kindBC) {
    case EULER_WALL:
      return BCWallInviscid;
    case HEAT_FLUX:
      return BCWallViscousHeatFlux;
    case ISOTHERMAL:
      return BCWallViscousIsothermal;
    case HEAT_TRANSFER:
    case CHT_WALL_INTERFACE:
    case SMOLUCHOWSKI_MAXWELL:
      return BCWallViscous;
    case FAR_FIELD:
      return BCFarfield;
    case SYMMETRY_PLANE:
      return BCSymmetryPlane;
    case INLET_FLOW:
      return BCInflow;
    case OUTLET_FLOW:
      return BCOutflow;
    case SUPERSONIC_INLET:
      return BCInflowSupersonic;
    case SUPERSONIC_OUTLET:
      return BCOutflowSupersonic;
    default:
      return BCTypeUserDefined;
  }
}

void CCGNSFileWriter::WriteField(int iField, const string& FieldName) {
  /*--- The coordinates define the mesh, so they are always written in double precision. Single precision would
   move the points by up to ~1e-7 of the size of the domain, which can be larger than the smallest cells. ---*/

  const bool isCoord = iField < nDim;

  if (isCoord || doublePrecisionFields)
    WriteFieldOfType<double>(iField, FieldName, RealDouble);
  else
    WriteFieldOfType<float>(iField, FieldName, RealSingle);
}

template <class T>
void CCGNSFileWriter::WriteFieldOfType(int iField, const string& FieldName, DataType_t dataType) {
  /*--- Check if field is coordinate. ---*/
  const bool isCoord = iField < nDim;

  /*--- Create send buffer. ---*/
  vector<T> sendBufferField(nLocalPoints);

  for (unsigned long iPoint = 0; iPoint < nLocalPoints; iPoint++) {
    sendBufferField[iPoint] = static_cast<T>(dataSorter->GetData(iField, iPoint));
  }

  if (rank != MASTER_NODE) {
    SendChunked(sendBufferField.data(), nLocalPoints * sizeof(T), MASTER_NODE, 0);
    return;
  }

  vector<T> recvBufferField;

  /*--- Coordinate vector is written in blocks, one for each process. ---*/
  cgsize_t nodeBegin = 1;
  auto nodeEnd = static_cast<cgsize_t>(nLocalPoints);
  if (nLocalPoints > 0) {
    if (isCoord) {
      int CoordinateNumber;
      CallCGNS(cg_coord_partial_write(cgnsFileID, cgnsBase, cgnsZone, dataType, FieldName.c_str(), &nodeBegin, &nodeEnd,
                                      sendBufferField.data(), &CoordinateNumber));
    } else {
      int fieldNumber;
      CallCGNS(cg_field_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, dataType, FieldName.c_str(), &nodeBegin,
                                      &nodeEnd, sendBufferField.data(), &fieldNumber));
    }
  }

  for (int i = 0; i < size; ++i) {
    if (i == MASTER_NODE) continue;
    /*--- In CGNS numbering starts form 1 and ranges are inclusive ---*/
    nodeBegin = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i) + 1);
    nodeEnd = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i + 1));

    const auto recvSize = static_cast<size_t>(nodeEnd - nodeBegin + 1);
    recvBufferField.resize(recvSize);

    RecvChunked(recvBufferField.data(), recvSize * sizeof(T), i, 0);
    if (recvSize == 0) continue;
    if (isCoord) {
      int CoordinateNumber;
      CallCGNS(cg_coord_partial_write(cgnsFileID, cgnsBase, cgnsZone, dataType, FieldName.c_str(), &nodeBegin, &nodeEnd,
                                      recvBufferField.data(), &CoordinateNumber));
    } else {
      int fieldNumber;
      CallCGNS(cg_field_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, dataType, FieldName.c_str(),
                                      &nodeBegin, &nodeEnd, recvBufferField.data(), &fieldNumber));
    }
  }
}

void CCGNSFileWriter::WriteConnectivity(GEO_TYPE type, const string& SectionName) {
  const auto nTotElem = dataSorter->GetnElemGlobal(type);
  if (nTotElem == 0) return;

  /*--- Create new CGNS nodes to store connectivity. Some readers (e.g. the VTK/ParaView CGNS reader)
   store the connectivity size of a section in a 32-bit int, so an element type with more than
   maxSectionEntries connectivity entries is split into several sections with consecutive ranges. ---*/
  const auto elementType = GetCGNSType(type);
  const auto nPointsElem = nPointsOfElementType(type);
  const auto nTotElemCG = static_cast<cgsize_t>(nTotElem);
  const auto maxElemSection = static_cast<cgsize_t>(maxSectionEntries / nPointsElem);
  const auto nSections = (nTotElemCG + maxElemSection - 1) / maxElemSection;

  /*--- First and last element (CGNS numbering starts from 1 and ranges are inclusive) of a section. ---*/
  auto sectionBegin = [&](cgsize_t iSec) { return cumulative + 1 + iSec * maxElemSection; };
  auto sectionEnd = [&](cgsize_t iSec) { return cumulative + std::min(nTotElemCG, (iSec + 1) * maxElemSection); };

  vector<int> cgnsSections(nSections);
  if (rank == MASTER_NODE) {
    for (cgsize_t iSec = 0; iSec < nSections; ++iSec) {
      const string name = nSections == 1 ? SectionName : SectionName + "_" + std::to_string(iSec + 1);
      CallCGNS(cg_section_partial_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), elementType,
                                        sectionBegin(iSec), sectionEnd(iSec), 0, &cgnsSections[iSec]));
    }
  }

  /*--- Write the connectivity of the elements [first, last], which may span more than one section. ---*/
  auto writeBlock = [&](cgsize_t first, cgsize_t last, const cgsize_t* conn) {
    for (cgsize_t iSec = 0; iSec < nSections; ++iSec) {
      const auto lo = std::max(first, sectionBegin(iSec));
      const auto hi = std::min(last, sectionEnd(iSec));
      if (lo > hi) continue;
      CallCGNS(cg_elements_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsSections[iSec], lo, hi,
                                         conn + (lo - first) * nPointsElem));
    }
  };

  /*--- Retrieve element distribution among processes. ---*/
  const auto nLocalElem = dataSorter->GetnElem(type);

  vector<unsigned long> distElem(size);

  SU2_MPI::Allgather(&nLocalElem, 1, MPI_UNSIGNED_LONG, distElem.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  cgsize_t firstElem = cumulative + 1;
  cgsize_t endElem = cumulative + static_cast<cgsize_t>(distElem[rank]);

  /*--- Connectivity is stored in send buffer. ---*/
  sendBufferConnectivity.resize(nLocalElem * nPointsElem);

  for (unsigned long iElem = 0; iElem < nLocalElem; iElem++) {
    for (unsigned long iPoint = 0; iPoint < nPointsElem; iPoint++) {
      sendBufferConnectivity[iPoint + nPointsElem * iElem] =
          static_cast<cgsize_t>(dataSorter->GetElemConnectivity(type, iElem, iPoint));
    }
  }

  if (rank != MASTER_NODE) {
    SendChunked(sendBufferConnectivity.data(), sendBufferConnectivity.size() * sizeof(cgsize_t), MASTER_NODE, 1);
    return;
  }

  /*--- Connectivity vector is written in blocks, one for each process. ---*/
  if (nLocalElem > 0) writeBlock(firstElem, endElem, sendBufferConnectivity.data());

  for (int i = 0; i < size; ++i) {
    if (i == MASTER_NODE) continue;
    /*--- In CGNS numbering starts form 1 and ranges are inclusive ---*/
    firstElem = endElem + 1;
    endElem += static_cast<cgsize_t>(distElem[i]);
    const auto recvSize = static_cast<size_t>(endElem - firstElem + 1) * nPointsElem;
    recvBufferConnectivity.resize(recvSize);

    RecvChunked(recvBufferConnectivity.data(), recvBufferConnectivity.size() * sizeof(cgsize_t), i, 1);

    if (!recvBufferConnectivity.empty()) writeBlock(firstElem, endElem, recvBufferConnectivity.data());
  }
  cumulative += static_cast<cgsize_t>(nTotElem);
}

void CCGNSFileWriter::SendChunked(const void* buf, size_t nBytes, int dest, int tag) {
  const auto* bytes = static_cast<const char*>(buf);
  for (size_t offset = 0; offset < nBytes; offset += maxChunkBytes) {
    const auto count = static_cast<int>(std::min(maxChunkBytes, nBytes - offset));
    SU2_MPI::Send(bytes + offset, count, MPI_CHAR, dest, tag, SU2_MPI::GetComm());
  }
}

void CCGNSFileWriter::RecvChunked(void* buf, size_t nBytes, int source, int tag) {
  auto* bytes = static_cast<char*>(buf);
  for (size_t offset = 0; offset < nBytes; offset += maxChunkBytes) {
    const auto count = static_cast<int>(std::min(maxChunkBytes, nBytes - offset));
    SU2_MPI::Recv(bytes + offset, count, MPI_CHAR, source, tag, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
  }
}

void CCGNSFileWriter::InitializeFields() {
  /*--- Create "Fields" node to store solution. ---*/
  if (rank == MASTER_NODE) CallCGNS(cg_sol_write(cgnsFileID, cgnsBase, cgnsZone, "Fields", Vertex, &cgnsFields));
}
#endif  // HAVE_CGNS
