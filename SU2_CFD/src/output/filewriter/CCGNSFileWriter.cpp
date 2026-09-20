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
#ifdef HAVE_MPI
  CallCGNS(cgp_close(cgnsFileID));
#else
  CallCGNS(cg_close(cgnsFileID));
#endif

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

  /*--- Remove the previous file if present, before any rank opens it. ---*/
  if (rank == MASTER_NODE) remove(val_filename.c_str());

#ifdef HAVE_MPI
  /*--- All ranks open the file and write their own part of the data with the parallel CGNS API. The nodes of the
   file (base, zone, sections, solution, fields, boundary conditions) are metadata and must be created by all
   ranks with the same arguments, only the data itself is written per rank. ---*/

  SU2_MPI::Barrier(SU2_MPI::GetComm());
  CallCGNS(cgp_mpi_comm(SU2_MPI::GetComm()));
  CallCGNS(cgp_pio_mode(CGP_COLLECTIVE));
  CallCGNS(cgp_open(val_filename.c_str(), CG_MODE_WRITE, &cgnsFileID));
#else
  CallCGNS(cg_open(val_filename.c_str(), CG_MODE_WRITE, &cgnsFileID));
#endif

  /*--- Create Base. ---*/
  CallCGNS(cg_base_write(cgnsFileID, "Base", nCell, nDim, &cgnsBase));
}

void CCGNSFileWriter::InitializeZone(const string& zoneName) {
  if (!dataSorter->GetConnectivitySorted()) {
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  nLocalPoints = dataSorter->GetnPoints();
  GlobalElem = static_cast<cgsize_t>(dataSorter->GetnElemGlobal());
  GlobalPoint = static_cast<cgsize_t>(dataSorter->GetnPointsGlobal());
  cumulative = 0;

  /*--- Create Zone. The number of cells does not include the boundary elements. ---*/
  array<cgsize_t, 3> zoneData;

  zoneData[0] = GlobalPoint;
  zoneData[1] = GlobalElem;
  zoneData[2] = 0;

  CallCGNS(
      cg_zone_write(cgnsFileID, cgnsBase, zoneName.substr(0, 32).c_str(), zoneData.data(), Unstructured, &cgnsZone));
}

void CCGNSFileWriter::WriteBoundaries() {
  for (const auto& marker : boundaryMarkers) {
    /*--- Count the elements of this rank and collect the element types it holds. The local connectivity holds the
     VTK type of each element followed by the ids of its nodes. ---*/

    unsigned long nLocalElem = 0, typesMask = 0;
    for (size_t pos = 0; pos < marker.conn.size();) {
      const auto type = static_cast<unsigned short>(marker.conn[pos]);
      typesMask |= 1ul << type;
      nLocalElem++;
      pos += nPointsOfElementType(type) + 1;
    }
    const unsigned long nLocalEntries = marker.conn.size() - nLocalElem;

    /*--- Sizes and offsets of the elements of each rank, which are written as a contiguous range. ---*/

    vector<unsigned long> elemPerRank(size), entriesPerRank(size);
    SU2_MPI::Allgather(&nLocalElem, 1, MPI_UNSIGNED_LONG, elemPerRank.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
    SU2_MPI::Allgather(&nLocalEntries, 1, MPI_UNSIGNED_LONG, entriesPerRank.data(), 1, MPI_UNSIGNED_LONG,
                       SU2_MPI::GetComm());

    const auto nTotElem = std::accumulate(elemPerRank.begin(), elemPerRank.end(), 0ul);
    if (nTotElem == 0) continue;

    auto elemOffset = std::accumulate(elemPerRank.begin(), elemPerRank.begin() + rank, 0ul);
    auto entryOffset = std::accumulate(entriesPerRank.begin(), entriesPerRank.begin() + rank, 0ul);

    /*--- A marker with a single element type is written as a section of that type, one with several types
     (e.g. triangles and quadrilaterals) as a MIXED section. ---*/

    unsigned long globalTypesMask = 0;
    SU2_MPI::Allreduce(&typesMask, &globalTypesMask, 1, MPI_UNSIGNED_LONG, MPI_BOR, SU2_MPI::GetComm());
    const bool singleType = (globalTypesMask & (globalTypesMask - 1)) == 0;

    const string name = marker.name.substr(0, 32);
    const cgsize_t range[2] = {cumulative + 1, cumulative + static_cast<cgsize_t>(nTotElem)};
    const cgsize_t first = range[0] + static_cast<cgsize_t>(elemOffset);
    const cgsize_t last = first + static_cast<cgsize_t>(nLocalElem) - 1;
    int section;

    if (singleType) {
      /*--- Only the ids of the nodes are stored, the type is that of the section. ---*/

      vector<cgsize_t> elems;
      elems.reserve(nLocalEntries);
      for (size_t pos = 0; pos < marker.conn.size();) {
        const auto nNodes = nPointsOfElementType(static_cast<unsigned short>(marker.conn[pos]));
        for (unsigned short iNode = 1; iNode <= nNodes; ++iNode)
          elems.push_back(static_cast<cgsize_t>(marker.conn[pos + iNode]));
        pos += nNodes + 1;
      }

      unsigned short type = 0;
      while ((globalTypesMask >> type) != 1) type++;

      CallCGNS(SectionWrite(name, GetCGNSType(type), range[0], range[1], &section));
      CallCGNS(ElementsWriteData(section, first, last, nLocalElem > 0 ? elems.data() : nullptr));

    } else {
      /*--- The CGNS element type of each element is stored before the ids of its nodes, and the start offset of
       each element in the connectivity array is stored in a second array. ---*/

      const auto nTotEntries = std::accumulate(entriesPerRank.begin(), entriesPerRank.end(), 0ul) + nTotElem;

      vector<cgsize_t> elems, offsets{static_cast<cgsize_t>(entryOffset + elemOffset)};
      elems.reserve(marker.conn.size());
      for (size_t pos = 0; pos < marker.conn.size();) {
        const auto type = static_cast<unsigned short>(marker.conn[pos]);
        const auto nNodes = nPointsOfElementType(type);
        elems.push_back(GetCGNSType(type));
        for (unsigned short iNode = 1; iNode <= nNodes; ++iNode)
          elems.push_back(static_cast<cgsize_t>(marker.conn[pos + iNode]));
        offsets.push_back(offsets.back() + nNodes + 1);
        pos += nNodes + 1;
      }

#ifdef HAVE_MPI
      CallCGNS(cgp_poly_section_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), MIXED, range[0], range[1],
                                      static_cast<cgsize_t>(nTotEntries), 0, &section));
      CallCGNS(cgp_poly_elements_write_data(cgnsFileID, cgnsBase, cgnsZone, section, first, last,
                                            nLocalElem > 0 ? elems.data() : nullptr,
                                            nLocalElem > 0 ? offsets.data() : nullptr));
#else
      CallCGNS(cg_poly_section_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), MIXED, range[0], range[1], 0,
                                     elems.data(), offsets.data(), &section));
#endif
    }
    cumulative += static_cast<cgsize_t>(nTotElem);

    /*--- The BC points to the boundary elements and takes its type from a family with the name of the marker.
     These are metadata nodes, written by all ranks. ---*/

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

  /*--- Each rank writes the points it holds, which are a contiguous range of the points of the file. ---*/
  vector<T> buffer(nLocalPoints);

  for (unsigned long iPoint = 0; iPoint < nLocalPoints; iPoint++) {
    buffer[iPoint] = static_cast<T>(dataSorter->GetData(iField, iPoint));
  }

  cgsize_t nodeBegin = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(rank) + 1);
  cgsize_t nodeEnd = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(rank + 1));

  /*--- A rank without points takes part in the collective call but writes nothing. ---*/
  const T* data = nLocalPoints > 0 ? buffer.data() : nullptr;

  if (isCoord) {
    int coordinateNumber;
#ifdef HAVE_MPI
    CallCGNS(cgp_coord_write(cgnsFileID, cgnsBase, cgnsZone, dataType, FieldName.c_str(), &coordinateNumber));
    CallCGNS(cgp_coord_write_data(cgnsFileID, cgnsBase, cgnsZone, coordinateNumber, &nodeBegin, &nodeEnd, data));
#else
    CallCGNS(cg_coord_partial_write(cgnsFileID, cgnsBase, cgnsZone, dataType, FieldName.c_str(), &nodeBegin, &nodeEnd,
                                    data, &coordinateNumber));
#endif
  } else {
    int fieldNumber;
#ifdef HAVE_MPI
    CallCGNS(cgp_field_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, dataType, FieldName.c_str(), &fieldNumber));
    CallCGNS(
        cgp_field_write_data(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, fieldNumber, &nodeBegin, &nodeEnd, data));
#else
    CallCGNS(cg_field_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, dataType, FieldName.c_str(),
                                    &nodeBegin, &nodeEnd, data, &fieldNumber));
#endif
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

  /*--- The sections are metadata, all ranks create them with the same arguments. ---*/
  vector<int> cgnsSections(nSections);
  for (cgsize_t iSec = 0; iSec < nSections; ++iSec) {
    const string name = nSections == 1 ? SectionName : SectionName + "_" + std::to_string(iSec + 1);
    CallCGNS(SectionWrite(name, elementType, sectionBegin(iSec), sectionEnd(iSec), &cgnsSections[iSec]));
  }

  /*--- Retrieve element distribution among processes, the elements of a rank are a contiguous range. ---*/
  const auto nLocalElem = dataSorter->GetnElem(type);

  vector<unsigned long> distElem(size);
  SU2_MPI::Allgather(&nLocalElem, 1, MPI_UNSIGNED_LONG, distElem.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  cgsize_t firstElem = cumulative + 1;
  for (int i = 0; i < rank; ++i) firstElem += static_cast<cgsize_t>(distElem[i]);
  const cgsize_t endElem = firstElem + static_cast<cgsize_t>(nLocalElem) - 1;

  /*--- Store the connectivity of this rank. ---*/
  vector<cgsize_t> connectivity(nLocalElem * nPointsElem);

  for (unsigned long iElem = 0; iElem < nLocalElem; iElem++) {
    for (unsigned long iPoint = 0; iPoint < nPointsElem; iPoint++) {
      connectivity[iPoint + nPointsElem * iElem] =
          static_cast<cgsize_t>(dataSorter->GetElemConnectivity(type, iElem, iPoint));
    }
  }

  /*--- Write the elements of this rank, which may span more than one section. A rank without elements takes
   part in the collective calls but writes nothing. ---*/
  for (cgsize_t iSec = 0; iSec < nSections; ++iSec) {
    const auto lo = std::max(firstElem, sectionBegin(iSec));
    const auto hi = std::min(endElem, sectionEnd(iSec));
    const bool empty = (nLocalElem == 0) || (lo > hi);
    CallCGNS(ElementsWriteData(cgnsSections[iSec], lo, hi, empty ? nullptr : &connectivity[(lo - firstElem) * nPointsElem]));
  }

  cumulative += nTotElemCG;
}

int CCGNSFileWriter::SectionWrite(const string& name, ElementType_t type, cgsize_t start, cgsize_t end, int* section) {
#ifdef HAVE_MPI
  return cgp_section_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), type, start, end, 0, section);
#else
  return cg_section_partial_write(cgnsFileID, cgnsBase, cgnsZone, name.c_str(), type, start, end, 0, section);
#endif
}

int CCGNSFileWriter::ElementsWriteData(int section, cgsize_t start, cgsize_t end, const cgsize_t* elements) {
#ifdef HAVE_MPI
  return cgp_elements_write_data(cgnsFileID, cgnsBase, cgnsZone, section, start, end, elements);
#else
  if (elements == nullptr) return CG_OK;
  return cg_elements_partial_write(cgnsFileID, cgnsBase, cgnsZone, section, start, end, elements);
#endif
}

void CCGNSFileWriter::InitializeFields() {
  /*--- Create "Fields" node to store solution. ---*/
  CallCGNS(cg_sol_write(cgnsFileID, cgnsBase, cgnsZone, "Fields", Vertex, &cgnsFields));
}
#endif  // HAVE_CGNS
