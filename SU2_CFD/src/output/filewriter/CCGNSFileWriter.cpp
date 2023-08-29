/*!
 * \file CCGNSFileWriter.cpp
 * \brief Filewriter class for CGNS format.
 * \author G. Baldan
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

const string CCGNSFileWriter::fileExt = ".cgns";

CCGNSFileWriter::CCGNSFileWriter(CParallelDataSorter* valDataSorter, bool isSurf)
    : CFileWriter(valDataSorter, fileExt), isSurface(isSurf) {}

void CCGNSFileWriter::WriteData(string val_filename) {

#ifdef HAVE_CGNS

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  /*--- Open the CGNS file for writing.  ---*/
  InitializeMeshFile(val_filename);

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

  /*--- Initialize and write fields. ---*/
  InitializeFields();

  const auto& fieldNames = dataSorter->GetFieldNames();
  for (unsigned long i = nDim; i < fieldNames.size(); ++i) {
    WriteField(i, fieldNames[i]);
  }

  /*--- Close the CGNS file. ---*/
  if (rank == MASTER_NODE) CallCGNS(cg_close(cgnsFileID));

#endif
}

#ifdef HAVE_CGNS
void CCGNSFileWriter::InitializeMeshFile(const string& val_filename) {
  if (!dataSorter->GetConnectivitySorted()) {
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  nLocalPoints = dataSorter->GetnPoints();
  nDim = dataSorter->GetnDim();
  GlobalElem = static_cast<cgsize_t>(dataSorter->GetnElemGlobal());
  GlobalPoint = static_cast<cgsize_t>(dataSorter->GetnPointsGlobal());
  cumulative = 0;

  /*--- If surface file cell dimension is decreased. ---*/
  const auto nCell = static_cast<int>(nDim - isSurface);

  if (rank == MASTER_NODE) {
    /*--- Remove the previous file if present. ---*/
    remove(val_filename.c_str());

    /*--- Create CGNS file and open in write mode. ---*/
    CallCGNS(cg_open(val_filename.c_str(), CG_MODE_WRITE, &cgnsFileID));

    /*--- Create Base. ---*/
    CallCGNS(cg_base_write(cgnsFileID, "Base", nCell, nDim, &cgnsBase));

    /*--- Create Zone. ---*/
    array<cgsize_t, 3> zoneData;

    zoneData[0] = GlobalPoint;
    zoneData[1] = GlobalElem;
    zoneData[2] = 0;

    CallCGNS(cg_zone_write(cgnsFileID, cgnsBase, "Zone", zoneData.data(), Unstructured, &cgnsZone));
  }
}

void CCGNSFileWriter::WriteField(int iField, const string& FieldName) {
  /*--- Check if field is coordinate. ---*/
  const bool isCoord = iField < nDim;

  /*--- Create send buffer. ---*/
  sendBufferField.resize(nLocalPoints);

  for (unsigned long iPoint = 0; iPoint < nLocalPoints; iPoint++) {
    sendBufferField[iPoint] = static_cast<dataPrecision>(dataSorter->GetData(iField, iPoint));
  }

  if (rank != MASTER_NODE) {
    SU2_MPI::Send(sendBufferField.data(), nLocalPoints * sizeof(dataPrecision), MPI_CHAR, MASTER_NODE, 0,
                  SU2_MPI::GetComm());
    return;
  }

  /*--- Coordinate vector is written in blocks, one for each process. ---*/
  cgsize_t nodeBegin = 1;
  auto nodeEnd = static_cast<cgsize_t>(nLocalPoints);

  if (isCoord) {
    int CoordinateNumber;
    CallCGNS(cg_coord_partial_write(cgnsFileID, cgnsBase, cgnsZone, dataType, FieldName.c_str(), &nodeBegin, &nodeEnd,
                                    sendBufferField.data(), &CoordinateNumber));
  } else {
    int fieldNumber;
    CallCGNS(cg_field_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, dataType, FieldName.c_str(), &nodeBegin,
                                    &nodeEnd, sendBufferField.data(), &fieldNumber));
  }

  for (int i = 0; i < size; ++i) {
    if (i == MASTER_NODE) continue;
    /*--- In CGNS numbering starts form 1 and ranges are inclusive ---*/
    nodeBegin = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i) + 1);
    nodeEnd = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i + 1));

    const auto recvSize = static_cast<int>(nodeEnd - nodeBegin + 1);
    recvBufferField.resize(recvSize);

    SU2_MPI::Recv(recvBufferField.data(), recvSize * sizeof(dataPrecision), MPI_CHAR, i, 0, SU2_MPI::GetComm(),
                  MPI_STATUS_IGNORE);
    if (recvSize <= 0) continue;
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

  /*--- Create a new CGNS node to store connectivity. ---*/
  const auto elementType = GetCGNSType(type);

  cgsize_t firstElem = cumulative + 1;
  cgsize_t endElem = cumulative + static_cast<cgsize_t>(nTotElem);

  int cgnsSection;
  if (rank == MASTER_NODE)
    CallCGNS(cg_section_partial_write(cgnsFileID, cgnsBase, cgnsZone, SectionName.c_str(), elementType, firstElem,
                                      endElem, 0, &cgnsSection));

  /*--- Retrieve element distribution among processes. ---*/
  const auto nLocalElem = dataSorter->GetnElem(type);

  vector<unsigned long> distElem(size);

  SU2_MPI::Allgather(&nLocalElem, 1, MPI_UNSIGNED_LONG, distElem.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  firstElem = cumulative + 1;
  endElem = cumulative + static_cast<cgsize_t>(distElem[rank]);

  /*--- Connectivity is stored in send buffer. ---*/
  const auto nPointsElem = nPointsOfElementType(type);
  sendBufferConnectivity.resize(nLocalElem * nPointsElem);

  for (unsigned long iElem = 0; iElem < nLocalElem; iElem++) {
    for (unsigned long iPoint = 0; iPoint < nPointsElem; iPoint++) {
      sendBufferConnectivity[iPoint + nPointsElem * iElem] =
          static_cast<cgsize_t>(dataSorter->GetElemConnectivity(type, iElem, iPoint));
    }
  }

  const auto bufferSize = static_cast<int>(nLocalElem * nPointsElem * sizeof(cgsize_t));
  if (rank != MASTER_NODE) {
    SU2_MPI::Send(sendBufferConnectivity.data(), bufferSize, MPI_CHAR, MASTER_NODE, 1, SU2_MPI::GetComm());
    return;
  }

  /*--- Connectivity vector is written in blocks, one for each process. ---*/
  if (nLocalElem > 0)
    CallCGNS(cg_elements_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsSection, firstElem, endElem,
                                       sendBufferConnectivity.data()));

  for (int i = 0; i < size; ++i) {
    if (i == MASTER_NODE) continue;
    /*--- In CGNS numbering starts form 1 and ranges are inclusive ---*/
    firstElem = endElem + 1;
    endElem += static_cast<cgsize_t>(distElem[i]);
    const auto recvSize = static_cast<int>((endElem - firstElem + 1) * nPointsElem);
    recvBufferConnectivity.resize(recvSize);

    const auto recvByte = static_cast<int>(recvBufferConnectivity.size() * sizeof(cgsize_t));
    SU2_MPI::Recv(recvBufferConnectivity.data(), recvByte, MPI_CHAR, i, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

    if (!recvBufferConnectivity.empty())
      CallCGNS(cg_elements_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsSection, firstElem, endElem,
                                         recvBufferConnectivity.data()));
  }
  cumulative += static_cast<cgsize_t>(nTotElem);
}

void CCGNSFileWriter::InitializeFields() {
  /*--- Create "Fields" node to store solution. ---*/
  if (rank == MASTER_NODE) CallCGNS(cg_sol_write(cgnsFileID, cgnsBase, cgnsZone, "Fields", Vertex, &cgnsFields));
}
#endif  // HAVE_CGNS
