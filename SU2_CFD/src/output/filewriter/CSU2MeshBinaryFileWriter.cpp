/*!
 * \file CSU2MeshBinaryFileWriter.cpp
 * \brief Filewriter class SU2 binary mesh format.
 * \author E. van der Weide
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
#include "../../../include/output/filewriter/CSU2MeshBinaryFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <vector>

const string CSU2MeshBinaryFileWriter::fileExt = ".su2b";

CSU2MeshBinaryFileWriter::CSU2MeshBinaryFileWriter(CParallelDataSorter *valDataSorter,
                                                   unsigned short valiZone, unsigned short valnZone) :
   CFileWriter(valDataSorter, fileExt), iZone(valiZone), nZone(valnZone) {}

namespace {

/*--- The writer always emits 8-byte (uint64_t) connectivity entries. This keeps
      the implementation simple (no branching on element/point counts) and is
      unconditionally readable by CSU2BinaryMeshReaderBase, which auto-detects
      the connectivity width from the file header. ---*/
using conn_t = uint64_t;
constexpr int32_t SU2B_CONN_TYPE_SIZE = static_cast<int32_t>(sizeof(conn_t));

/*--- Elements are always visited in this fixed order, matching CSU2MeshFileWriter
      (the ASCII writer) so that the two formats produce numerically identical
      meshes for the same input. ---*/
constexpr std::array<GEO_TYPE, 6> ElemTypes = {TRIANGLE, QUADRILATERAL, TETRAHEDRON,
                                                HEXAHEDRON, PRISM, PYRAMID};

FILE* OpenAppend(const string& filename) {
  FILE* f = fopen(filename.c_str(), "ab");
  if (!f) SU2_MPI::Error(string("Unable to open file ") + filename, CURRENT_FUNCTION);
  return f;
}

}  // namespace

void CSU2MeshBinaryFileWriter::WriteData(string val_filename) {

  val_filename.append(fileExt);

  /*--- Write the file-level header (only once, before the first zone) followed
        by the per-zone header (zone_id, n_dim, n_elem). Only rank 0 touches the
        file for this; the implicit synchronization at the first Allreduce below
        (also present in CSU2MeshFileWriter) keeps the other ranks from writing
        before this is done. ---*/

  if (rank == 0) {
    FILE* f = fopen(val_filename.c_str(), (iZone == 0) ? "wb" : "ab");
    if (!f) SU2_MPI::Error(string("Unable to open file ") + val_filename, CURRENT_FUNCTION);

    if (iZone == 0) {
      int32_t size_conn_type = SU2B_CONN_TYPE_SIZE;
      int32_t n_zone = nZone;
      fwrite(&size_conn_type, sizeof(size_conn_type), 1, f);
      fwrite(&n_zone, sizeof(n_zone), 1, f);
    }

    int32_t zone_id = iZone;
    int32_t n_dim = dataSorter->GetnDim();
    conn_t n_elem = dataSorter->GetnElemGlobal();
    fwrite(&zone_id, sizeof(zone_id), 1, f);
    fwrite(&n_dim, sizeof(n_dim), 1, f);
    fwrite(&n_elem, sizeof(n_elem), 1, f);

    fclose(f);
  }

  /*--- Section 1: element offsets. Every rank streams, in turn, the starting
        connectivity-array position of each of its local elements (visited in
        the fixed type order above), starting from the cumulative total left
        behind by the previous ranks. Each rank's accumulator holds only its
        own local contribution (like CSU2MeshFileWriter's nElem/myPoint), so
        summing it across ranks via Allreduce yields the new cumulative total.
        Once all ranks are done, the final total is appended once more as the
        closing sentinel offset[n_elem]. ---*/

  unsigned long connOffset = 0, localConnOffset = 0;

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      FILE* f = OpenAppend(val_filename);
      unsigned long running = connOffset;
      for (auto type : ElemTypes) {
        const conn_t nPointsElem = nPointsOfElementType(type);
        for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
          conn_t value = running;
          fwrite(&value, sizeof(value), 1, f);
          running += nPointsElem + 2;
        }
      }
      fclose(f);
      localConnOffset = running - connOffset;
    }
    SU2_MPI::Allreduce(&localConnOffset, &connOffset, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  }

  if (rank == 0) {
    FILE* f = OpenAppend(val_filename);
    conn_t sentinel = connOffset;
    fwrite(&sentinel, sizeof(sentinel), 1, f);
    fclose(f);
  }

  /*--- Section 2: element connectivity, [VTK_Type, node_0..node_n-1, GlobalIndex]
        per element, in the same fixed type order and the same rank-by-rank
        streaming pattern as CSU2MeshFileWriter uses for the ASCII format
        (including the same "-1" to convert 1-based dataSorter indices to the
        0-based indices used throughout the SU2 mesh formats). ---*/

  unsigned long elemIndexOffset = 0, localElemCount = 0;

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      FILE* f = OpenAppend(val_filename);
      conn_t globalIndex = elemIndexOffset;
      for (auto type : ElemTypes) {
        const auto nPointsElem = nPointsOfElementType(type);
        for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
          conn_t vtkType = type;
          fwrite(&vtkType, sizeof(vtkType), 1, f);
          for (auto iNode = 0u; iNode < nPointsElem; iNode++) {
            conn_t node = dataSorter->GetElemConnectivity(type, iElem, iNode) - 1;
            fwrite(&node, sizeof(node), 1, f);
          }
          fwrite(&globalIndex, sizeof(globalIndex), 1, f);
          globalIndex++;
        }
      }
      fclose(f);
      localElemCount = static_cast<unsigned long>(globalIndex - elemIndexOffset);
    }
    SU2_MPI::Allreduce(&localElemCount, &elemIndexOffset, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  }

  /*--- Section 3: point coordinates and IDs, interleaved. Same rank-by-rank
        streaming pattern as CSU2MeshFileWriter's point section. ---*/

  if (rank == 0) {
    FILE* f = OpenAppend(val_filename);
    conn_t nPointsGlobal = dataSorter->GetnPointsGlobal();
    fwrite(&nPointsGlobal, sizeof(nPointsGlobal), 1, f);
    fclose(f);
  }

  unsigned long myPoint = 0, pointOffset = 0;

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      FILE* f = OpenAppend(val_filename);
      for (auto iPoint = 0ul; iPoint < dataSorter->GetnPoints(); iPoint++) {
        for (auto iDim = 0u; iDim < dataSorter->GetnDim(); iDim++) {
          double coord = dataSorter->GetData(iDim, iPoint);
          fwrite(&coord, sizeof(coord), 1, f);
        }
        conn_t pointID = iPoint + pointOffset;
        fwrite(&pointID, sizeof(pointID), 1, f);
      }
      fclose(f);
      myPoint = dataSorter->GetnPoints();
    }
    SU2_MPI::Allreduce(&myPoint, &pointOffset, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  }

  /*--- Section 4: markers. Mirrors CSU2MeshFileWriter: the marker connectivity
        is not available from the data sorter, so it is read back from the
        "boundary[_iZone].dat" file that CPhysicalGeometry writes (for
        SU2_COMPONENT::SU2_DEF) right after reading the original mesh. Only the
        master rank does this work, exactly as for the ASCII format. ---*/

  if (rank == MASTER_NODE) {
    FILE* f = OpenAppend(val_filename);

    string str = "boundary";
    if (nZone > 1) str += "_" + PrintingToolbox::to_string(iZone);
    str += ".dat";

    ifstream input_file(str);
    if (!input_file.is_open()) SU2_MPI::Error(string("Cannot find ") + str, CURRENT_FUNCTION);

    string text_line;
    while (getline(input_file, text_line)) {

      auto position = text_line.find("NMARK=", 0);
      if (position == string::npos) continue;

      text_line.erase(0, 6);
      const int32_t nMarker_ = atoi(text_line.c_str());
      fwrite(&nMarker_, sizeof(nMarker_), 1, f);

      for (int iMarker = 0; iMarker < nMarker_; iMarker++) {

        getline(input_file, text_line);
        text_line.erase(0, 11);
        for (int iChar = 0; iChar < 20; iChar++) {
          position = text_line.find(' ', 0);
          if (position != string::npos) text_line.erase(position, 1);
          position = text_line.find('\r', 0);
          if (position != string::npos) text_line.erase(position, 1);
          position = text_line.find('\n', 0);
          if (position != string::npos) text_line.erase(position, 1);
        }
        string Marker_Tag = text_line;

        /*--- Write the null-padded marker name. ---*/
        char name_buf[CGNS_STRING_SIZE] = {};
        strncpy(name_buf, Marker_Tag.c_str(), CGNS_STRING_SIZE - 1);
        fwrite(name_buf, sizeof(char), CGNS_STRING_SIZE, f);

        getline(input_file, text_line);
        text_line.erase(0, 13);
        const unsigned long nElem_Bound_ = atoi(text_line.c_str());

        /*--- Consume (but do not use) the SEND_TO= line: periodic/send-receive
              markers are not supported by the binary format yet, matching
              CSU2BinaryMeshReaderBase, which errors out on SEND_RECEIVE. ---*/
        getline(input_file, text_line);

        /*--- Parse this marker's boundary elements into memory, then emit the
              offset array followed by the connectivity array, matching what
              CSU2BinaryMeshReaderBase::ReadSurfaceElementConnectivity expects. ---*/
        vector<unsigned short> vtkTypes(nElem_Bound_);
        vector<std::array<unsigned long, N_POINTS_QUADRILATERAL>> nodes(nElem_Bound_);
        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
          getline(input_file, text_line);
          istringstream bound_line(text_line);
          unsigned short VTK_Type;
          bound_line >> VTK_Type;
          vtkTypes[iElem_Bound] = VTK_Type;
          const auto nPointsElem = nPointsOfElementType(VTK_Type);
          for (unsigned short iNode = 0; iNode < nPointsElem; iNode++) bound_line >> nodes[iElem_Bound][iNode];
        }

        auto nElemBoundConn = static_cast<conn_t>(nElem_Bound_);
        fwrite(&nElemBoundConn, sizeof(nElemBoundConn), 1, f);

        unsigned long running = 0;
        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
          conn_t value = running;
          fwrite(&value, sizeof(value), 1, f);
          running += nPointsOfElementType(vtkTypes[iElem_Bound]) + 1;
        }
        conn_t sentinel = running;
        fwrite(&sentinel, sizeof(sentinel), 1, f);

        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
          conn_t vtkType = vtkTypes[iElem_Bound];
          fwrite(&vtkType, sizeof(vtkType), 1, f);
          const auto nPointsElem = nPointsOfElementType(vtkTypes[iElem_Bound]);
          for (unsigned short iNode = 0; iNode < nPointsElem; iNode++) {
            conn_t node = nodes[iElem_Bound][iNode];
            fwrite(&node, sizeof(node), 1, f);
          }
        }
      }
    }

    input_file.close();
    fclose(f);
  }

  SU2_MPI::Barrier(SU2_MPI::GetComm());
}
