/*!
 * \file CSU2BinaryMeshReaderBase.cpp
 * \brief Helper class for the reading of a native SU2 binary grid file.
 * \author T. Economon, E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CSU2BinaryMeshReaderBase.hpp"

#include <set>

CSU2BinaryMeshReaderBase::CSU2BinaryMeshReaderBase(CConfig* val_config, unsigned short val_iZone,
                                                   unsigned short val_nZone)
    : CSU2MeshReaderBase(val_config, val_iZone, val_nZone) {}

CSU2BinaryMeshReaderBase::~CSU2BinaryMeshReaderBase(void) = default;

void CSU2BinaryMeshReaderBase::ReadConnectivityType() {
  /*--- Initialize the byte swapping to false and
        read the size of the connectivity type. ---*/
  swap_bytes = false;
  ReadBinaryData<int>(&size_conn_type, 1);

  /*--- Check if byte swapping must be applied. ---*/
  if ((size_conn_type != 4) && (size_conn_type != 8)) {
    SwapBytes((char*)&size_conn_type, sizeof(int), 1);
    swap_bytes = true;
  }

  /*--- The size of the connectivity type must be either 4 or 8. ---*/
  if ((size_conn_type != 4) && (size_conn_type != 8))
    SU2_MPI::Error(string("The file ") + meshFilename + string(" is not a valid SU2 binary file"), CURRENT_FUNCTION);
}

void CSU2BinaryMeshReaderBase::ReadMetadata(CConfig* config) {
  const bool harmonic_balance = config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE;
  const bool multizone_file = config->GetMultizone_Mesh();

  /*--- Open the grid file and check if it went OK. ---*/
  mesh_file = fopen(meshFilename.c_str(), "rb");
  if (!mesh_file)
    SU2_MPI::Error(
        string("Error opening SU2 binary grid file ") + meshFilename + string(". Check that the file exists"),
        CURRENT_FUNCTION);

  /*--- Read the size of the connectivity type and check if byte swapping
        must be applied. ---*/
  ReadConnectivityType();

  /*--- Check for a harmonic balance simulation. If So, the data of the
        first zone can be read. Otherwise jump to the location of the
        current zone. ---*/
  if (harmonic_balance) {
    if (rank == MASTER_NODE) cout << "Reading time instance " << config->GetiInst() + 1 << "." << endl;
    fseek(mesh_file, 2 * sizeof(int), SEEK_SET);
  } else {
    FastForwardToMyZone();
    if (nZones > 1 && multizone_file) {
      if (rank == MASTER_NODE) cout << "Reading zone " << myZone << " from native SU2 binary mesh." << endl;
    }
  }

  /*--- Read the meta data from the current position. ---*/
  ReadMetadataZone();

  /*--- Close the grid file again. ---*/
  fclose(mesh_file);
}

void CSU2BinaryMeshReaderBase::ReadPointCoordinates() {
  /* No support yet for actuator disks */
  if (actuator_disk) SU2_MPI::Error("No support for actuator disks yet", CURRENT_FUNCTION);

  /* Jump over the number of points, because it is already known, and
     determine the position in the file where the point section ends. */
  fseek(mesh_file, size_conn_type, SEEK_CUR);
  auto pos_end_point = ftell(mesh_file) + numberOfGlobalPoints * (dimension * sizeof(double) + size_conn_type);

  /* Define a linear partitioner for the points. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints, 0);

  /* Jump to the position in the file where the points are stored
     that this rank must read. */
  const auto firstIndex = pointPartitioner.GetFirstIndexOnRank(rank);
  fseek(mesh_file, firstIndex * (dimension * sizeof(double) + size_conn_type), SEEK_CUR);

  /* Determine the number of local points and prepare the local
     data structure to store the point coordinates. */
  numberOfLocalPoints = pointPartitioner.GetSizeOnRank(rank);
  localPointCoordinates.resize(dimension);
  for (int k = 0; k < dimension; k++) localPointCoordinates[k].resize(numberOfLocalPoints);

  /*--- Read the point coordinates into our data structure. ---*/
  for (unsigned long i = 0; i < numberOfLocalPoints; ++i) {
    double Coords[3];
    ReadBinaryData<double>(Coords, dimension);
    fseek(mesh_file, size_conn_type, SEEK_CUR);

    for (unsigned short iDim = 0; iDim < dimension; iDim++) {
      localPointCoordinates[iDim][i] = Coords[iDim];
    }
  }

  /* Jump to the end of the coordinate section. */
  fseek(mesh_file, pos_end_point, SEEK_SET);
}

void CSU2BinaryMeshReaderBase::ReadVolumeElementConnectivity() {
  /* Jump over the zone ID, number of dimensions and number of elements,
     because this information is already known. */
  fseek(mesh_file, size_conn_type + 2 * sizeof(int), SEEK_CUR);

  /* Get a linear partitioner of the elements. */
  CLinearPartitioner elemPartitioner(numberOfGlobalElements, 0);

  /* Determine the position at the end of the offset array, where
     the total size of the connectivity is stored. */
  const auto first_index = elemPartitioner.GetFirstIndexOnRank(rank);
  const auto pos_size_global_conn = ftell(mesh_file) + numberOfGlobalElements * size_conn_type;

  /* Jump to position in the file where the offset data is stored for
     the element ragne this rank will read. Allocate the memory for
     this offset array. */
  fseek(mesh_file, first_index * size_conn_type, SEEK_CUR);
  numberOfLocalElements = elemPartitioner.GetSizeOnRank(rank);
  vector<uint64_t> offset(numberOfLocalElements + 1);

  /* Read the offset array from the file. It will be stored as uint64_t,
     but it may be stored differently in the file. */
  if (size_conn_type == 4) {
    vector<uint32_t> tmp(offset.size());
    ReadBinaryData<uint32_t>(tmp.data(), tmp.size());
    for (size_t i = 0; i < tmp.size(); ++i) offset[i] = static_cast<uint64_t>(tmp[i]);
  } else {
    ReadBinaryData<uint64_t>(offset.data(), offset.size());
  }

  /* Jump to the location where the total size of the connectivity array
     is stored and read the size. Determine the location of the end of
     the connectivity array. */
  fseek(mesh_file, pos_size_global_conn, SEEK_SET);
  const auto size__global_conn = ReadBinaryNEntities();
  const auto pos_end_conn = ftell(mesh_file) + size__global_conn * size_conn_type;

  /* Read the connectivity data of the elements this rank should read.
     Store the data in uint64_t. */
  const auto size_conn = offset.back() - offset[0];
  vector<uint64_t> conn_buff(size_conn);
  fseek(mesh_file, offset[0] * size_conn_type, SEEK_CUR);

  if (size_conn_type == 4) {
    vector<uint32_t> tmp(conn_buff.size());
    ReadBinaryData<uint32_t>(tmp.data(), tmp.size());
    for (size_t i = 0; i < tmp.size(); ++i) conn_buff[i] = static_cast<uint64_t>(tmp[i]);
  } else {
    ReadBinaryData<uint64_t>(conn_buff.data(), conn_buff.size());
  }

  /* Jump to the end of the connectivity data. */
  fseek(mesh_file, pos_end_conn, SEEK_SET);

#ifdef HAVE_MPI

  /* Update the offset, such that it corresponds to the data in the
     local connectivity buffer. */
  for (size_t i = 1; i < offset.size(); ++i) offset[i] -= offset[0];
  offset[0] = 0;

  /* Get a linear partitioner of the points. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints, 0);

  /*--- Determine the ranks on which the elements must actually be stored.
        Note that an element can be stored on multiple ranks, as the points
        must be surrounded by all its elements. ---*/
  std::vector<int> ranks_elements;
  ranks_elements.reserve(numberOfLocalElements);
  std::vector<uint64_t> number_of_ranks_elements(numberOfLocalElements + 1);
  number_of_ranks_elements[0] = 0;

  for (unsigned long i = 0; i < numberOfLocalElements; ++i) {
    /* Determine the ranks where this element must be stored by looping
       over its nodes. */
    set<int> ranks_this_elem;
    for (uint64_t j = offset[i] + 1; j < (offset[i + 1] - 1); ++j) {
      auto rank_node = pointPartitioner.GetRankContainingIndex(conn_buff[j]);
      ranks_this_elem.insert(static_cast<int>(rank_node));
    }

    /* Store the data. */
    number_of_ranks_elements[i + 1] = number_of_ranks_elements[i] + ranks_this_elem.size();
    for (auto rank_elem : ranks_this_elem) ranks_elements.push_back(rank_elem);
  }

  /* Create the send buffers. Both the size of each connectivity
     information and the connectivity information itself is stored.*/
  std::vector<std::vector<unsigned long>> send_buf;
  send_buf.resize(size);

  for (unsigned long i = 0; i < numberOfLocalElements; ++i) {
    for (uint64_t j = number_of_ranks_elements[i]; j < number_of_ranks_elements[i + 1]; ++j) {
      const int ii = ranks_elements[j];

      auto size_this_conn = offset[i + 1] - offset[i];
      send_buf[ii].push_back(size_this_conn);
      for (uint64_t k = offset[i]; k < offset[i + 1]; ++k) send_buf[ii].push_back(conn_buff[k]);
    }
  }

  /* Determine the number of ranks from which this rank will receive data.
     Allow for self communication. */
  int nRankRecv;
  vector<int> sendToRank(size, 0), sizeSend(size, 1);
  for (auto rank_elem : ranks_elements) sendToRank[rank_elem] = 1;
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeSend.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /* Explicitly delete the memory that is not needed anymore. */
  vector<uint64_t>().swap(offset);
  vector<uint64_t>().swap(conn_buff);
  vector<int>().swap(sizeSend);

  /* Determine the number of ranks to which this rank will send data
     and allocate the memory for the send requests. */
  int nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i]) ++nRankSend;
  }

  vector<SU2_MPI::Request> sendReqs(nRankSend);

  /* Send the data using non-blocking sends. */
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i]) {
      SU2_MPI::Isend(send_buf[i].data(), send_buf[i].size(), MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(),
                     &sendReqs[nRankSend]);
      ++nRankSend;
    }
  }

  /* Define the receive buffers and receive the messages. */
  std::vector<std::vector<unsigned long>> recv_buf;
  recv_buf.resize(size);

  for (int i = 0; i < nRankRecv; ++i) {
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int rankRecv = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);
    recv_buf[rankRecv].resize(sizeMess);
    SU2_MPI::Recv(recv_buf[rankRecv].data(), sizeMess, MPI_UNSIGNED_LONG, rankRecv, rank, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends and release the memory of the send buffers. */
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i]) vector<unsigned long>().swap(send_buf[i]);
  }

  /* Synchronize the MPI ranks, because wild cards have been used. */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

  /*--- Store the information in the receive buffers in the offset and conn_buff
        vectors, such that it is consistent with the information without MPI.
        Release the memory of the receive buffers afterwards. ---*/
  offset.push_back(0);
  for (int i = 0; i < size; ++i) {
    if (recv_buf[i].size() > 0) {
      size_t ind = 0;
      while (ind < recv_buf[i].size()) {
        auto n_items = recv_buf[i][ind++];
        offset.push_back(offset.back() + n_items);
        for (unsigned long j = 0; j < n_items; ++j, ++ind) conn_buff.push_back(recv_buf[i][ind]);
      }
      vector<unsigned long>().swap(recv_buf[i]);
    }
  }

#endif

  /*--- Extract the connectivity data from conn_buf and store the data
        in the appropriate member variables. ---*/
  numberOfLocalElements = offset.size() - 1;
  array<unsigned long, N_POINTS_HEXAHEDRON> connectivity{};

  for (unsigned long i = 0; i < numberOfLocalElements; ++i) {
    auto ind = offset[i];
    auto size_this_elem = static_cast<unsigned short>(offset[i + 1] - offset[i]);
    auto VTK_Type = conn_buff[ind++];
    const auto nPointsElem = nPointsOfElementType(static_cast<unsigned short>(VTK_Type));
    if (size_this_elem < (nPointsElem + 2)) SU2_MPI::Error("Not enough items in volume connectivity", CURRENT_FUNCTION);

    for (unsigned short j = 0; j < nPointsElem; ++j, ++ind) connectivity[j] = conn_buff[ind];
    auto GlobalIndex = conn_buff[ind];

    localVolumeElementConnectivity.push_back(GlobalIndex);
    localVolumeElementConnectivity.push_back(VTK_Type);
    /// TODO: Use a compressed format.
    for (unsigned short j = 0; j < N_POINTS_HEXAHEDRON; ++j) {
      localVolumeElementConnectivity.push_back(connectivity[j]);
    }
  }
}

void CSU2BinaryMeshReaderBase::ReadSurfaceElementConnectivity() {
  /* The number of surface markers is already known, so jump over it. */
  fseek(mesh_file, sizeof(int), SEEK_CUR);

  /* Allocate the memory for the first index of the connectivty of the
     surface elements and the marker names. Note that all ranks store
     the entire surface connectivity. */
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);

  array<unsigned long, N_POINTS_HEXAHEDRON> connectivity{};

  /* Loop over the number of markers. */
  for (unsigned long iMarker = 0; iMarker < numberOfMarkers; ++iMarker) {
    /* Read the name of the surface marker. */
    char charStr[SU2_STRING_SIZE];
    ReadBinaryData<char>(charStr, SU2_STRING_SIZE);
    markerNames[iMarker] = string(charStr);

    /*--- Throw an error if we find deprecated references to SEND_RECEIVE
           boundaries in the mesh. ---*/
    if (markerNames[iMarker] == "SEND_RECEIVE")
      SU2_MPI::Error(
          "Mesh file contains deprecated SEND_RECEIVE marker!\n"
          "Please remove any SEND_RECEIVE markers from the SU2 binary mesh.",
          CURRENT_FUNCTION);

    /* Read the number of elements for this boundary marker. */
    const auto nElem_Bound = ReadBinaryNEntities();

    /*--- Read the offset array from the file. It will be stored as
          uint64_t, but it may be stored differently in the file. ---*/
    vector<uint64_t> offset(nElem_Bound + 1);
    if (size_conn_type == 4) {
      vector<uint32_t> tmp(offset.size());
      ReadBinaryData<uint32_t>(tmp.data(), tmp.size());
      for (size_t i = 0; i < tmp.size(); ++i) offset[i] = static_cast<uint64_t>(tmp[i]);
    } else {
      ReadBinaryData<uint64_t>(offset.data(), offset.size());
    }

    /*--- Read the connectivity and store it in a buffer.
          Always use uint64_t for this internally. ---*/
    vector<uint64_t> conn_buff(offset.back());
    if (size_conn_type == 4) {
      vector<uint32_t> tmp(conn_buff.size());
      ReadBinaryData<uint32_t>(tmp.data(), tmp.size());
      for (size_t i = 0; i < tmp.size(); ++i) conn_buff[i] = static_cast<uint64_t>(tmp[i]);
    } else {
      ReadBinaryData<uint64_t>(conn_buff.data(), conn_buff.size());
    }

    /*--- Loop over the surface elements to store the connectivity
          in the reequired data structures. ---*/
    for (unsigned long i = 0; i < nElem_Bound; ++i) {
      auto ind = offset[i];
      auto size_this_elem = static_cast<unsigned short>(offset[i + 1] - offset[i]);
      auto VTK_Type = conn_buff[ind++];
      const auto nPointsElem = nPointsOfElementType(static_cast<unsigned short>(VTK_Type));
      if (size_this_elem < (nPointsElem + 1))
        SU2_MPI::Error("Not enough items in surface connectivity", CURRENT_FUNCTION);

      if (dimension == 3 && VTK_Type == LINE) {
        SU2_MPI::Error(
            "Line boundary conditions are not possible for 3D calculations.\n"
            "Please check the SU2 binary file.",
            CURRENT_FUNCTION);
      }

      for (unsigned short j = 0; j < nPointsElem; ++j, ++ind) connectivity[j] = conn_buff[ind];

      surfaceElementConnectivity[iMarker].push_back(0);
      surfaceElementConnectivity[iMarker].push_back(VTK_Type);
      for (unsigned short j = 0; j < N_POINTS_HEXAHEDRON; ++j) {
        surfaceElementConnectivity[iMarker].push_back(connectivity[j]);
      }
    }
  }
}

void CSU2BinaryMeshReaderBase::FastForwardToMyZone() {
  /*--- Jump to the position where the data starts for the first zone. ---*/
  fseek(mesh_file, 2 * sizeof(int), SEEK_SET);

  /*--- Loop over the lower numbered zones and read their meta data. ---*/
  for (int zone = 0; zone < myZone; ++zone) ReadMetadataZone();
}

uint64_t CSU2BinaryMeshReaderBase::ReadBinaryNEntities() {
  /*--- Define the return value as an uint64_t. ---*/
  uint64_t nEntities;

  /*--- Read the actual data, depending on the connectivity type. ---*/
  if (size_conn_type == 4) {
    uint32_t dummy;
    ReadBinaryData<uint32_t>(&dummy, 1);
    nEntities = static_cast<uint64_t>(dummy);
  } else {
    ReadBinaryData<uint64_t>(&nEntities, 1);
  }

  return nEntities;
}

void CSU2BinaryMeshReaderBase::ReadMetadataZone() {
  /*--- Skip the zone ID and read the number of dimensions. ---*/
  int nDim;
  fseek(mesh_file, sizeof(int), SEEK_CUR);
  ReadBinaryData<int>(&nDim, 1);
  dimension = static_cast<unsigned short>(nDim);

  /*--- Read the number of elements. ---*/
  const auto nElem = ReadBinaryNEntities();
  numberOfGlobalElements = static_cast<unsigned long>(nElem);

  /*--- Jump to the end of the offset section, read the size of
        the connectivity and jump over it. ---*/
  fseek(mesh_file, nElem * size_conn_type, SEEK_CUR);
  const auto size_conn = ReadBinaryNEntities();
  fseek(mesh_file, size_conn * size_conn_type, SEEK_CUR);

  /*--- Read the number of points and jump over the coordinate section. ---*/
  const auto nPoints = ReadBinaryNEntities();
  numberOfGlobalPoints = static_cast<unsigned long>(nPoints);
  fseek(mesh_file, nPoints * (nDim * sizeof(double) + size_conn_type), SEEK_CUR);

  /*--- Read the number of markers and loop over them. ---*/
  int nMark;
  ReadBinaryData<int>(&nMark, 1);
  numberOfMarkers = static_cast<unsigned long>(nMark);

  for (int mark = 0; mark < nMark; ++mark) {
    /*--- Jump over the name of the marker and read
          the number of surface elements. ---*/
    fseek(mesh_file, SU2_STRING_SIZE * sizeof(char), SEEK_CUR);
    const auto nElemMark = ReadBinaryNEntities();

    /*--- Jump to the end of the offset section of this marker, read
          the size of the connectivity and jump over it. ---*/
    fseek(mesh_file, nElemMark * size_conn_type, SEEK_CUR);
    const auto size_conn_mark = ReadBinaryNEntities();
    fseek(mesh_file, size_conn_mark * size_conn_type, SEEK_CUR);
  }
}
