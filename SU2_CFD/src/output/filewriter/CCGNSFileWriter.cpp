/*!
 * \file CCGNSFileWriter.cpp
 * \brief Filewriter class for CGNS format.
 * \author G. Baldan
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

CCGNSFileWriter::CCGNSFileWriter(string valFileName, CParallelDataSorter* valDataSorter, bool isSurf)
    : CFileWriter(std::move(valFileName), valDataSorter, fileExt) {
  isSurface = isSurf;
}

CCGNSFileWriter::~CCGNSFileWriter() {}

void CCGNSFileWriter::Write_Data() {
#ifdef HAVE_CGNS
  /*--- Open the CGNS file for writing.  ---*/
  initializeMeshFile();

  /*--- Write point coordinates. ---*/
  WriteCoordinate(1, "CoordinateX");
  WriteCoordinate(2, "CoordinateY");
  if (nDim == 3) WriteCoordinate(3, "CoordinateZ");

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

  const vector<string> fieldNames = dataSorter->GetFieldNames();
  for (unsigned long i = nDim; i < fieldNames.size(); ++i) {
    WriteField(i, fieldNames[i]);
  }

  /*--- Close the CGNS file. ---*/
  if (rank == MASTER_NODE) CallCGNS(cg_close(cgnsFileID));

#endif
}

#ifdef HAVE_CGNS
void CCGNSFileWriter::initializeMeshFile() {
  if (!dataSorter->GetConnectivitySorted()) {
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  int nCell;

  nLocalPoints = static_cast<int>(dataSorter->GetnPoints());
  nDim = dataSorter->GetnDim();
  GlobalElem = static_cast<cgsize_t>(dataSorter->GetnElemGlobal());
  GlobalPoint = static_cast<cgsize_t>(dataSorter->GetnPointsGlobal());
  cumulative = 0;

  /*--- If surface file cell dimension is decreased. ---*/
  nCell = nDim - isSurface;

  if (rank == MASTER_NODE) {
    /*--- Remove the previous file if present. ---*/
    remove(fileName.c_str());

    /*--- Create CGNS file and open in write mode. ---*/
    CallCGNS(cg_open(fileName.c_str(), CG_MODE_WRITE, &cgnsFileID));

    /*--- Create Base. ---*/
    CallCGNS(cg_base_write(cgnsFileID, "Base", nCell, nDim, &cgnsBase));

    /*--- Create Zone. ---*/
    vector<cgsize_t> zoneData(3);

    zoneData[0] = GlobalPoint;
    zoneData[1] = GlobalElem;
    zoneData[2] = 0;

    CallCGNS(cg_zone_write(cgnsFileID, cgnsBase, "Zone", zoneData.data(), Unstructured, &cgnsZone));
  }
}

void CCGNSFileWriter::WriteCoordinate(int CoordinateNumber, const string& CoordinateName) {
  cgsize_t nodeBegin;
  cgsize_t nodeEnd;
  vector<float> sendBufferFloat, recvBufferFloat;

  /*--- Create send buffer. ---*/
  sendBufferFloat.resize(nLocalPoints);

  unsigned short iCoord = CoordinateNumber - 1;
  for (int iPoint = 0; iPoint < nLocalPoints; iPoint++) {
    sendBufferFloat[iPoint] = static_cast<float>(dataSorter->GetData(iCoord, iPoint));
  }

#ifdef HAVE_MPI
  SU2_MPI::Request sendReq;
  if (rank != MASTER_NODE)
    SU2_MPI::Isend(sendBufferFloat.data(), nLocalPoints, MPI_FLOAT, MASTER_NODE, 0, SU2_MPI::GetComm(), &sendReq);
#endif

  /*--- Coordinate vector is written in blocks, one for each process. ---*/
  if (rank == MASTER_NODE) {
    nodeBegin = 1;
    nodeEnd = static_cast<cgsize_t>(nLocalPoints);
    CallCGNS(cg_coord_partial_write(cgnsFileID, cgnsBase, cgnsZone, RealSingle, CoordinateName.c_str(), &nodeBegin,
                                    &nodeEnd, sendBufferFloat.data(), &CoordinateNumber));
  }

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) {
    for (int i = 1; i < size; ++i) {
      nodeBegin = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i) + 1);
      nodeEnd = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i + 1));

      int recvSize = static_cast<int>(nodeEnd - nodeBegin + 1);
      recvBufferFloat.resize(recvSize);

      SU2_MPI::Recv(recvBufferFloat.data(), recvSize, MPI_FLOAT, i, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      if (recvSize <= 0) continue;
      CallCGNS(cg_coord_partial_write(cgnsFileID, cgnsBase, cgnsZone, RealSingle, CoordinateName.c_str(), &nodeBegin,
                                      &nodeEnd, recvBufferFloat.data(), &CoordinateNumber));
    }
  }

  if (rank != MASTER_NODE) SU2_MPI::Wait(&sendReq, MPI_STATUS_IGNORE);
#endif
}

void CCGNSFileWriter::WriteConnectivity(GEO_TYPE type, const string& SectionName) {
  cgsize_t firstElem, endElem;

  unsigned long iElem, iPoint;
  unsigned long nLocalElem, nTotElem;
  unsigned short nPointsElem;
  int cgnsSection;

  ElementType_t elementType;

#ifdef HAVE_MPI
  SU2_MPI::Request sendReq;
#endif

  nLocalElem = dataSorter->GetnElem(type);
  nTotElem = dataSorter->GetnElemGlobal(type);

  nPointsElem = nPointsOfElementType(type);

  elementType = GetCGNSType(type);

  if (nTotElem > 0) {
    firstElem = cumulative + 1;
    endElem = cumulative + static_cast<cgsize_t>(nTotElem);

    /*--- Create a new CGNS node to store connectivity. ---*/
    if (rank == MASTER_NODE)
      CallCGNS(cg_section_partial_write(cgnsFileID, cgnsBase, cgnsZone, SectionName.c_str(), elementType, firstElem,
                                        endElem, 0, &cgnsSection));

      /*--- Retrieve element distribution among processes. ---*/
#ifdef HAVE_MPI
    vector<unsigned long> distElem;
    distElem.resize(size);

    SU2_MPI::Allgather(&nLocalElem, 1, MPI_UNSIGNED_LONG, distElem.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

    firstElem = cumulative + 1;
    endElem = cumulative + static_cast<cgsize_t>(distElem[rank]);

    vector<cgsize_t> sendBufferConnectivity, recvBufferConnectivity;
#else
    firstElem = cumulative + 1;
    endElem = cumulative + static_cast<cgsize_t>(nLocalElem);
    vector<cgsize_t> sendBufferConnectivity;
#endif

    /*--- Connectivity is stored in send buffer. ---*/
    sendBufferConnectivity.resize(nLocalElem * nPointsElem);

    for (iElem = 0; iElem < nLocalElem; iElem++) {
      for (iPoint = 0; iPoint < nPointsElem; iPoint++) {
        sendBufferConnectivity[iPoint + nPointsElem * iElem] =
            static_cast<cgsize_t>(dataSorter->GetElem_Connectivity(type, iElem, iPoint));
      }
    }

#ifdef HAVE_MPI
    auto bufferSize = static_cast<int>(nLocalElem * nPointsElem * sizeof(cgsize_t));
    if (rank != MASTER_NODE)
      SU2_MPI::Isend(sendBufferConnectivity.data(), bufferSize, MPI_CHAR, MASTER_NODE, rank, SU2_MPI::GetComm(),
                     &sendReq);
#endif

    /*--- Connectivity vector is written in blocks, one for each process. ---*/
    if (rank == MASTER_NODE) {
      if (nLocalElem > 0)
        CallCGNS(cg_elements_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsSection, firstElem, endElem,
                                           sendBufferConnectivity.data()));
    }
#ifdef HAVE_MPI
    if (rank == MASTER_NODE) {
      for (int i = 1; i < size; ++i) {
        firstElem += static_cast<cgsize_t>(distElem[i - 1]);
        endElem += static_cast<cgsize_t>(distElem[i]);
        int recvSize = static_cast<int>(endElem - firstElem + 1);
        recvBufferConnectivity.resize(recvSize * nPointsElem);

        auto recvByte = static_cast<int>(recvSize * nPointsElem * sizeof(cgsize_t));
        SU2_MPI::Recv(recvBufferConnectivity.data(), recvByte, MPI_CHAR, i, i, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

        if (!recvBufferConnectivity.empty())
          CallCGNS(cg_elements_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsSection, firstElem, endElem,
                                             recvBufferConnectivity.data()));
      }
      cumulative += static_cast<cgsize_t>(nTotElem);
    }

    if (rank != MASTER_NODE) SU2_MPI::Wait(&sendReq, MPI_STATUS_IGNORE);
#endif
  }
}

void CCGNSFileWriter::InitializeFields() {
  /*--- Create "Fields" node to store solution. ---*/
  if (rank == MASTER_NODE) CallCGNS(cg_sol_write(cgnsFileID, cgnsBase, cgnsZone, "Fields", Vertex, &cgnsFields));
}

void CCGNSFileWriter::WriteField(int iField, const string& FieldName) {
  cgsize_t nodeBegin;
  cgsize_t nodeEnd;
  vector<float> sendBufferFloat, recvBufferFloat;

  /*--- Solution is stored in send buffer. ---*/
  sendBufferFloat.resize(nLocalPoints);

  for (int iPoint = 0; iPoint < nLocalPoints; iPoint++) {
    sendBufferFloat[iPoint] = static_cast<float>(dataSorter->GetData(iField, iPoint));
  }

  int fieldNumber;

#ifdef HAVE_MPI
  SU2_MPI::Request sendReq;
  if (rank != MASTER_NODE)
    SU2_MPI::Isend(sendBufferFloat.data(), nLocalPoints, MPI_FLOAT, MASTER_NODE, 0, SU2_MPI::GetComm(), &sendReq);
#endif

  /*--- Field vector is written in blocks, one for each process. ---*/
  if (rank == MASTER_NODE) {
    nodeBegin = 1;
    nodeEnd = static_cast<cgsize_t>(nLocalPoints);

    CallCGNS(cg_field_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, RealSingle, FieldName.c_str(),
                                    &nodeBegin, &nodeEnd, sendBufferFloat.data(), &fieldNumber));
  }

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) {
    for (int i = 1; i < size; ++i) {
      nodeBegin = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i) + 1);
      nodeEnd = static_cast<cgsize_t>(dataSorter->GetnPointCumulative(i + 1));

      int recvSize = static_cast<int>(nodeEnd - nodeBegin + 1);

      recvBufferFloat.resize(recvSize);

      SU2_MPI::Recv(recvBufferFloat.data(), recvSize, MPI_FLOAT, i, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      if (recvSize <= 0) continue;
      CallCGNS(cg_field_partial_write(cgnsFileID, cgnsBase, cgnsZone, cgnsFields, RealSingle, FieldName.c_str(),
                                      &nodeBegin, &nodeEnd, recvBufferFloat.data(), &fieldNumber));
    }
  }

  if (rank != MASTER_NODE) SU2_MPI::Wait(&sendReq, MPI_STATUS_IGNORE);
#endif
}
#endif  // HAVE_CGNS
