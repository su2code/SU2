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

#include <numeric>
#include <tuple>
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

/*--- fopen() creates new files with permissive default permissions (typically
      0666 before umask, i.e. potentially world-writable). Explicitly restrict
      to owner read/write, group/other read-only, matching the permissions of
      the ASCII mesh files written elsewhere via ofstream. A no-op on Windows,
      which does not use POSIX permission bits. ---*/
void RestrictPermissions(const string& filename) {
#if !defined(_WIN32)
  chmod(filename.c_str(), S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
#endif
}

/*--- Append the bytes of a value to a buffer, the buffers of all ranks are written to the file together. ---*/
template <class T>
void AppendBytes(string& buffer, const T& value) {
  buffer.append(reinterpret_cast<const char*>(&value), sizeof(value));
}

}  // namespace

void CSU2MeshBinaryFileWriter::WriteData(string val_filename) {

  /*--- For multizone cases all zones are written into one file, the zones after the first are appended. ---*/

  OpenMPIFile(val_filename, iZone != 0);

  if (rank == MASTER_NODE) RestrictPermissions(val_filename + fileExt);

  /*--- Write the file-level header (only once, before the first zone) followed by the per-zone header
        (zone_id, n_dim, n_elem). Only the master node writes it. ---*/

  string buffer;

  if (iZone == 0) {
    AppendBytes(buffer, SU2B_CONN_TYPE_SIZE);
    AppendBytes(buffer, static_cast<int32_t>(nZone));
  }

  /*--- Zone IDs are 1-based, matching the "IZONE=" convention of the ASCII format. ---*/
  AppendBytes(buffer, static_cast<int32_t>(iZone + 1));
  AppendBytes(buffer, static_cast<int32_t>(dataSorter->GetnDim()));
  AppendBytes(buffer, static_cast<conn_t>(dataSorter->GetnElemGlobal()));

  WriteMPIString(buffer, MASTER_NODE);

  /*--- Each rank writes the data of its own elements and points, at the position that follows the data of the
        ranks before it. The global offsets and indices of a rank are those of the ranks before it. ---*/

  auto offsetOfRank = [&](unsigned long localCount) {
    vector<unsigned long> counts(size, localCount);
    SU2_MPI::Allgather(&localCount, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
    return std::make_pair(std::accumulate(counts.begin(), counts.begin() + rank, 0ul),
                          std::accumulate(counts.begin(), counts.end(), 0ul));
  };

  /*--- Section 1: element offsets, the starting connectivity-array position of each element, closed by the
        total size as sentinel offset[n_elem]. ---*/

  unsigned long localConnSize = 0, localElemCount = 0;
  for (auto type : ElemTypes) {
    localConnSize += dataSorter->GetnElem(type) * (nPointsOfElementType(type) + 2);
    localElemCount += dataSorter->GetnElem(type);
  }

  unsigned long connOffset, totalConnSize;
  std::tie(connOffset, totalConnSize) = offsetOfRank(localConnSize);

  unsigned long elemOffset, nElemGlobal;
  std::tie(elemOffset, nElemGlobal) = offsetOfRank(localElemCount);

  buffer.clear();
  buffer.reserve(localElemCount * sizeof(conn_t));

  unsigned long running = connOffset;
  for (auto type : ElemTypes) {
    const conn_t nPointsElem = nPointsOfElementType(type);
    for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
      AppendBytes(buffer, static_cast<conn_t>(running));
      running += nPointsElem + 2;
    }
  }
  WriteMPIStringAll(buffer);

  buffer.clear();
  AppendBytes(buffer, static_cast<conn_t>(totalConnSize));
  WriteMPIString(buffer, MASTER_NODE);

  /*--- Section 2: element connectivity, [VTK_Type, node_0..node_n-1, GlobalIndex] per element. ---*/

  buffer.clear();
  buffer.reserve(localConnSize * sizeof(conn_t));

  conn_t globalIndex = elemOffset;
  for (auto type : ElemTypes) {
    const auto nPointsElem = nPointsOfElementType(type);
    for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
      AppendBytes(buffer, static_cast<conn_t>(type));
      for (auto iNode = 0u; iNode < nPointsElem; iNode++)
        AppendBytes(buffer, static_cast<conn_t>(dataSorter->GetElemConnectivity(type, iElem, iNode) - 1));
      AppendBytes(buffer, globalIndex);
      globalIndex++;
    }
  }
  WriteMPIStringAll(buffer);

  /*--- Section 3: point coordinates and IDs, interleaved. ---*/

  buffer.clear();
  AppendBytes(buffer, static_cast<conn_t>(dataSorter->GetnPointsGlobal()));
  WriteMPIString(buffer, MASTER_NODE);

  unsigned long pointOffset, nPointsTotal;
  std::tie(pointOffset, nPointsTotal) = offsetOfRank(dataSorter->GetnPoints());

  buffer.clear();
  buffer.reserve(dataSorter->GetnPoints() * (dataSorter->GetnDim() * sizeof(double) + sizeof(conn_t)));

  for (auto iPoint = 0ul; iPoint < dataSorter->GetnPoints(); iPoint++) {
    for (auto iDim = 0u; iDim < dataSorter->GetnDim(); iDim++)
      AppendBytes(buffer, static_cast<double>(dataSorter->GetData(iDim, iPoint)));
    AppendBytes(buffer, static_cast<conn_t>(iPoint + pointOffset));
  }
  WriteMPIStringAll(buffer);

  /*--- Section 4: markers. Mirrors CSU2MeshFileWriter: the marker connectivity
        is not available from the data sorter, so it is read back from the
        "boundary[_iZone].dat" file that CPhysicalGeometry writes (for
        SU2_COMPONENT::SU2_DEF) right after reading the original mesh. Only the
        master rank does this work, exactly as for the ASCII format. ---*/

  buffer.clear();

  if (rank == MASTER_NODE) {

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
      AppendBytes(buffer, nMarker_);

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

        /*--- Write the null-padded marker name. Uses SU2_BINARY_STRING_SIZE
              (shared with CSU2BinaryMeshReaderBase::SU2_STRING_SIZE) rather
              than CGNS_STRING_SIZE, since the two formats' name field widths
              are independent and must not silently drift apart. ---*/
        char name_buf[SU2_BINARY_STRING_SIZE] = {};
        strncpy(name_buf, Marker_Tag.c_str(), SU2_BINARY_STRING_SIZE - 1);
        buffer.append(name_buf, SU2_BINARY_STRING_SIZE);

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
        AppendBytes(buffer, nElemBoundConn);

        unsigned long boundOffset = 0;
        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
          conn_t value = boundOffset;
          AppendBytes(buffer, value);
          boundOffset += nPointsOfElementType(vtkTypes[iElem_Bound]) + 1;
        }
        conn_t sentinel = boundOffset;
        AppendBytes(buffer, sentinel);

        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
          conn_t vtkType = vtkTypes[iElem_Bound];
          AppendBytes(buffer, vtkType);
          const auto nPointsElem = nPointsOfElementType(vtkTypes[iElem_Bound]);
          for (unsigned short iNode = 0; iNode < nPointsElem; iNode++) {
            conn_t node = nodes[iElem_Bound][iNode];
            AppendBytes(buffer, node);
          }
        }
      }
    }

    input_file.close();
  }

  /*--- Only the master node has the markers, the other ranks write nothing. ---*/

  WriteMPIString(buffer, MASTER_NODE);

  CloseMPIFile();
}
