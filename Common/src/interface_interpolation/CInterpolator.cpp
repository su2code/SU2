/*!
 * \file CInterpolator.cpp
 * \brief Definition of the base class for interface interpolation.
 * \author H. Kline
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

#include "../../include/interface_interpolation/CInterpolator.hpp"

#include <set>

#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"

CInterpolator::CInterpolator(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                             unsigned int jZone)
    : rank(SU2_MPI::GetRank()),
      size(SU2_MPI::GetSize()),
      donorZone(iZone),
      targetZone(jZone),
      Geometry(geometry_container),
      donor_geometry(geometry_container[iZone][INST_0][MESH_0]),
      target_geometry(geometry_container[jZone][INST_0][MESH_0]) {}

bool CInterpolator::CheckInterfaceBoundary(int markDonor, int markTarget) {
  /*--- Determine whether the boundary is not on the rank because of
   *    the partition or because it is not part of the zone. ---*/
  int donorCheck = -1, targetCheck = -1;
  SU2_MPI::Allreduce(&markDonor, &donorCheck, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&markTarget, &targetCheck, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
  return (donorCheck != -1) && (targetCheck != -1);
}

bool CInterpolator::CheckZonesInterface(const CConfig* donor, const CConfig* target) {
  /*--- Loop over all interface markers to find if the 2 zones share any interface boundary. ---*/
  for (auto iInter = 0; iInter < (donor->GetMarker_n_ZoneInterface() / 2); iInter++) {
    if (CheckInterfaceBoundary(donor->FindInterfaceMarker(iInter), target->FindInterfaceMarker(iInter))) return true;
  }
  return false;
}

void CInterpolator::Determine_ArraySize(int markDonor, int markTarget, unsigned long nVertexDonor,
                                        unsigned short nDim) {
  /*--- Count donor vertices. ---*/
  auto nLocalVertex_Donor = 0ul;
  for (auto iVertex = 0ul; iVertex < nVertexDonor; iVertex++) {
    auto iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    nLocalVertex_Donor += donor_geometry->nodes->GetDomain(iPointDonor);
  }

  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;

  /*--- Send Interface vertex information --*/
  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1,
                     MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
}

void CInterpolator::Collect_VertexInfo(int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {
  Buffer_Send_Coord.resize(MaxLocalVertex_Donor, nDim);
  Buffer_Send_GlobalPoint.resize(MaxLocalVertex_Donor) = -1;

  /*--- Copy coordinates and point to the auxiliar vector --*/
  unsigned long iLocalVertexDonor = 0;

  for (unsigned long iVertex = 0; iVertex < nVertexDonor; iVertex++) {
    auto iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    if (donor_geometry->nodes->GetDomain(iPointDonor)) {
      Buffer_Send_GlobalPoint[iLocalVertexDonor] = donor_geometry->nodes->GetGlobalIndex(iPointDonor);
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord(iLocalVertexDonor, iDim) = donor_geometry->nodes->GetCoord(iPointDonor, iDim);
      iLocalVertexDonor++;
    }
  }

  SU2_MPI::Allgather(Buffer_Send_Coord.data(), MaxLocalVertex_Donor * nDim, MPI_DOUBLE, Buffer_Receive_Coord.data(),
                     MaxLocalVertex_Donor * nDim, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_GlobalPoint.data(), MaxLocalVertex_Donor, MPI_UNSIGNED_LONG,
                     Buffer_Receive_GlobalPoint.data(), MaxLocalVertex_Donor, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
}

unsigned long CInterpolator::Collect_ElementInfo(int markDonor, unsigned short nDim, bool compress,
                                                 vector<unsigned long>& allNumElem, vector<unsigned short>& numNodes,
                                                 su2matrix<long>& idxNodes) const {
  const auto maxElemNodes = (nDim == 2u) ? 2u : 4u;  // line and quad respectively

  unsigned long nElemDonor = 0;
  if (markDonor != -1) nElemDonor = donor_geometry->GetnElem_Bound(markDonor);

  allNumElem.resize(size);
  SU2_MPI::Allgather(&nElemDonor, 1, MPI_UNSIGNED_LONG, allNumElem.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  auto nMaxElemDonor = *max_element(allNumElem.begin(), allNumElem.end());

  vector<unsigned short> bufferSendNum(nMaxElemDonor);
  su2matrix<long> bufferSendIdx(nMaxElemDonor, maxElemNodes);

  numNodes.resize(nMaxElemDonor * size);
  idxNodes.resize(nMaxElemDonor * size, maxElemNodes);

  for (auto iElem = 0ul; iElem < nElemDonor; ++iElem) {
    const auto nNode = donor_geometry->bound[markDonor][iElem]->GetnNodes();
    bufferSendNum[iElem] = nNode;
    assert(nNode <= maxElemNodes && "Donor element has too many nodes.");

    for (auto iNode = 0u; iNode < nNode; ++iNode) {
      auto iPoint = donor_geometry->bound[markDonor][iElem]->GetNode(iNode);
      auto iPointGlobal = donor_geometry->nodes->GetGlobalIndex(iPoint);
      bufferSendIdx(iElem, iNode) = iPointGlobal;
    }
  }

  SU2_MPI::Allgather(bufferSendNum.data(), bufferSendNum.size(), MPI_UNSIGNED_SHORT, numNodes.data(),
                     bufferSendNum.size(), MPI_UNSIGNED_SHORT, SU2_MPI::GetComm());
  SU2_MPI::Allgather(bufferSendIdx.data(), bufferSendIdx.size(), MPI_LONG, idxNodes.data(), bufferSendIdx.size(),
                     MPI_LONG, SU2_MPI::GetComm());

  if (!compress) return accumulate(allNumElem.begin(), allNumElem.end(), 0ul);

  /*--- Compress the information (overlapping copy do not use memcpy). ---*/

  unsigned long dstIdx = 0;
  for (int iProcessor = 0; iProcessor < size; ++iProcessor) {
    auto srcOffset = iProcessor * nMaxElemDonor;
    for (auto idx = 0u; idx < allNumElem[iProcessor]; ++idx) {
      numNodes[dstIdx] = numNodes[srcOffset + idx];
      for (auto iNode = 0u; iNode < maxElemNodes; ++iNode) idxNodes(dstIdx, iNode) = idxNodes(srcOffset + idx, iNode);
      ++dstIdx;
    }
  }

  return dstIdx;
}

void CInterpolator::ReconstructBoundary(unsigned long val_zone, int val_marker) {
  const CGeometry* geom = Geometry[val_zone][INST_0][MESH_0];
  const auto nDim = geom->GetnDim();

  /*--- If this zone has no parts of the marker, it will not participate
   * in the first part of this function (collection). ---*/
  unsigned long nVertex, nElems;
  if (val_marker != -1) {
    nVertex = geom->GetnVertex(val_marker);
    nElems = geom->GetnElem_Bound(val_marker);
  } else {
    nVertex = 0;
    nElems = 0;
  }

  /*--- Get the number of domain vertices on the marker, and a mapping
   * (iVertex) -> (iLocalVertex, non-domain points being ignored). ---*/
  su2vector<unsigned long> iVertex_to_iLocalVertex(nVertex);
  unsigned long nLocalVertex = 0;
  for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
    const auto iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      iVertex_to_iLocalVertex[iVertex] = nLocalVertex;
      nLocalVertex++;
    } else {
      iVertex_to_iLocalVertex[iVertex] = numeric_limits<unsigned long>::max();
    }
  }

  // coordinates of all domain vertices on the marker
  su2activematrix Buffer_Send_Coord(nLocalVertex, nDim);
  // global point IDs of all domain vertices on the marker
  su2vector<unsigned long> Buffer_Send_GlobalPoint(nVertex);

  /*--- Assign to each domain vertex on the marker, identified by local point ID,
   * a set of surface-neighbor vertices on the marker, identified by global point ID. ---*/
  unordered_map<unsigned long, set<unsigned long> > neighbors;

  /*--- Define or initialize them. ---*/
  for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
    const auto iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      const auto iLocalVertex = iVertex_to_iLocalVertex[iVertex];
      Buffer_Send_GlobalPoint[iLocalVertex] = geom->nodes->GetGlobalIndex(iPoint);
      for (unsigned long iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord(iLocalVertex, iDim) = geom->nodes->GetCoord(iPoint, iDim);
      neighbors.insert(pair<unsigned long, set<unsigned long> >(iPoint, set<unsigned long>()));
    }
  }

  /*--- Define the neighbors map. ---*/
  for (unsigned long iElem = 0; iElem < nElems; iElem++) {
    const CPrimalGrid* elem = geom->bound[val_marker][iElem];
    for (unsigned short iNode = 0; iNode < elem->GetnNodes(); iNode++) {
      const auto iPoint = elem->GetNode(iNode);
      if (geom->nodes->GetDomain(iPoint)) {
        set<unsigned long>& neighb = neighbors.at(iPoint);
        for (unsigned short iNeighbor = 0; iNeighbor < elem->GetnNeighbor_Nodes(iNode); iNeighbor++) {
          unsigned long jPoint = elem->GetNode(elem->GetNeighbor_Nodes(iNode, iNeighbor));
          unsigned long jPoint_global = geom->nodes->GetGlobalIndex(jPoint);
          neighb.insert(jPoint_global);
        }
      }
    }
  }

  // numbers of surface-neighbors of all domain vertices on the marker
  su2vector<unsigned long> Buffer_Send_nLinkedNodes(nLocalVertex);
  // cumsum of Buffer_Send_nLinkedNodes
  su2vector<unsigned long> Buffer_Send_StartLinkedNodes(nLocalVertex);
  unsigned long nLocalLinkedNodes = 0;
  for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
    const auto iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      const auto iLocalVertex = iVertex_to_iLocalVertex[iVertex];
      Buffer_Send_nLinkedNodes[iLocalVertex] = neighbors.at(iPoint).size();
      Buffer_Send_StartLinkedNodes[iLocalVertex] = nLocalLinkedNodes;
      nLocalLinkedNodes += Buffer_Send_nLinkedNodes[iLocalVertex];
    }
  }
  // global point IDs of surface-neighbors of all domain vertices on the marker
  su2vector<unsigned long> Buffer_Send_LinkedNodes(nLocalLinkedNodes);
  unsigned long index = 0;
  for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
    const auto iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      for (unsigned long jPoint_global : neighbors.at(iPoint)) {
        Buffer_Send_LinkedNodes[index] = jPoint_global;
        index++;
      }
    }
  }

  /*--- Now these arrays of all processes must be joined to a single/global arrays. For this,
   * the entries of StartLinkedNodes must be shifted.
   * Furthermore, the global point IDs in LinkedNodes are replaced by global vertex IDs.
   * For this, the master process collects the data from all processes, joins them and broadcasts them again. ---*/

  /*--- Allocate global arrays. ---*/
  unsigned long nGlobalLinkedNodes;
  SU2_MPI::Allreduce(&nLocalVertex, &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nLocalLinkedNodes, &nGlobalLinkedNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  Buffer_Receive_Coord.resize(nGlobalVertex, nDim);
  Buffer_Receive_GlobalPoint.resize(nGlobalVertex);
  Buffer_Receive_Proc.resize(nGlobalVertex);

  Buffer_Receive_nLinkedNodes.resize(nGlobalVertex);
  Buffer_Receive_LinkedNodes.resize(nGlobalLinkedNodes);
  Buffer_Receive_StartLinkedNodes.resize(nGlobalVertex);

  /*--- Master process gathers data from all ranks ---*/
#ifdef HAVE_MPI
  int nProcessor = size, iRank;
  unsigned long received_nLocalLinkedNodes, received_nLocalVertex;

  if (rank == MASTER_NODE) {
    /*--- "Receive" from master process, i.e. copy. ---*/
    for (unsigned long iVertex = 0; iVertex < nLocalVertex; iVertex++)
      for (unsigned long iDim = 0; iDim < nDim; iDim++)
        Buffer_Receive_Coord(iVertex, iDim) = Buffer_Send_Coord(iVertex, iDim);

    for (unsigned long iVertex = 0; iVertex < nLocalVertex; iVertex++) {
      Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Proc[iVertex] = MASTER_NODE;
      Buffer_Receive_nLinkedNodes[iVertex] = Buffer_Send_nLinkedNodes[iVertex];
      Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
    }

    for (unsigned long iVertex = 0; iVertex < nLocalLinkedNodes; iVertex++)
      Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];

    unsigned long received_nLocalVertex_sum = nLocalVertex;
    unsigned long received_nLocalLinkedNodes_sum = nLocalLinkedNodes;

    /*--- Receive from other processes to the appropriate position in the buffers and shift StartLinkedNodes indices.
     * ---*/
    for (iRank = 1; iRank < nProcessor; iRank++) {
      SU2_MPI::Recv(&received_nLocalLinkedNodes, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_LinkedNodes[received_nLocalLinkedNodes_sum], received_nLocalLinkedNodes,
                    MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      SU2_MPI::Recv(&received_nLocalVertex, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(Buffer_Receive_Coord[received_nLocalVertex_sum], nDim * received_nLocalVertex, MPI_DOUBLE, iRank, 1,
                    SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      SU2_MPI::Recv(&Buffer_Receive_GlobalPoint[received_nLocalVertex_sum], received_nLocalVertex, MPI_UNSIGNED_LONG,
                    iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_nLinkedNodes[received_nLocalVertex_sum], received_nLocalVertex, MPI_UNSIGNED_LONG,
                    iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_StartLinkedNodes[received_nLocalVertex_sum], received_nLocalVertex,
                    MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      for (unsigned long iVertex = 0; iVertex < received_nLocalVertex; iVertex++) {
        Buffer_Receive_Proc[received_nLocalVertex_sum + iVertex] = iRank;
        Buffer_Receive_StartLinkedNodes[received_nLocalVertex_sum + iVertex] += received_nLocalLinkedNodes_sum;
      }

      received_nLocalVertex_sum += received_nLocalVertex;
      received_nLocalLinkedNodes_sum += received_nLocalLinkedNodes;
    }
  } else {
    SU2_MPI::Send(&nLocalLinkedNodes, 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_LinkedNodes.data(), nLocalLinkedNodes, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(&nLocalVertex, 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_Coord.data(), Buffer_Send_Coord.size(), MPI_DOUBLE, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(Buffer_Send_GlobalPoint.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_nLinkedNodes.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_StartLinkedNodes.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
  }
#else
  for (unsigned long iVertex = 0; iVertex < nLocalVertex; iVertex++)
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Buffer_Receive_Coord(iVertex, iDim) = Buffer_Send_Coord(iVertex, iDim);

  for (unsigned long iVertex = 0; iVertex < nGlobalVertex; iVertex++) {
    Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
    Buffer_Receive_Proc[iVertex] = MASTER_NODE;
    Buffer_Receive_nLinkedNodes[iVertex] = Buffer_Send_nLinkedNodes[iVertex];
    Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
  }

  for (unsigned long iVertex = 0; iVertex < nGlobalLinkedNodes; iVertex++)
    Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];
#endif

  /*--- Master process replaces global point indices in Buffer_Receive_LinkedNodes by their indices in
   * Buffer_Receive_GlobalPoint, Buffer_Receive_nLinkedNodes etc. ---*/
  if (rank == MASTER_NODE) {
    for (unsigned long iVertex = 0; iVertex < nGlobalVertex; iVertex++) {
      unsigned long* uptr = &Buffer_Receive_LinkedNodes[Buffer_Receive_StartLinkedNodes[iVertex]];

      for (unsigned long jLinkedNode = 0; jLinkedNode < Buffer_Receive_nLinkedNodes[iVertex]; jLinkedNode++) {
        const auto jPoint = uptr[jLinkedNode];
        bool found = false;  // Global point index has been found
        for (unsigned long kVertex = 0; kVertex < nGlobalVertex; kVertex++) {
          if (Buffer_Receive_GlobalPoint[kVertex] == jPoint) {
            uptr[jLinkedNode] = kVertex;
            found = true;
            break;
          }
        }

        if (!found) {  // remove from list
          for (unsigned long kVertex = jLinkedNode; kVertex < Buffer_Receive_nLinkedNodes[iVertex] - 1; kVertex++) {
            uptr[kVertex] = uptr[kVertex + 1];
          }
          Buffer_Receive_nLinkedNodes[iVertex]--;
          jLinkedNode--;
        }
      }
    }
  }

  SU2_MPI::Bcast(Buffer_Receive_GlobalPoint.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_Coord.data(), Buffer_Receive_Coord.size(), MPI_DOUBLE, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_Proc.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());

  SU2_MPI::Bcast(Buffer_Receive_nLinkedNodes.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_StartLinkedNodes.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_LinkedNodes.data(), nGlobalLinkedNodes, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
}
