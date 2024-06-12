/*!
 * \file CInterpolator.cpp
 * \brief Definition of the base class for interface interpolation.
 * \author H. Kline
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

void CInterpolator::ReconstructBoundary_Extended(const CConfig* const* config, unsigned long val_zone, int val_marker) {
  const CGeometry* geom = Geometry[val_zone][INST_0][MESH_0];

  unsigned long iVertex, jVertex, kVertex;

  unsigned long count, iTmp, *uptr, dPoint, EdgeIndex, jEdge, nEdges, nNodes, nVertex, iDim, nDim, iPoint;

  unsigned long nGlobalLinkedNodes, nLocalVertex, nLocalLinkedNodes;

  /*---modification made here---*/
  unsigned long nGlobalSurroundNodes, nLocalSurroundNodes;
  unsigned long nGlobalLinkedElems, nLocalLinkedElems;
  unsigned long iElem, nElems, nSurroundNodes, mElem, mElemNodes, jNode, jPoint, kNode;

  nDim = geom->GetnDim();

  if (val_marker != -1)
    nVertex = geom->GetnVertex(val_marker);
  else
    nVertex = 0;

  /*--- Get the number of domain vertices on the marker. ---*/
  su2vector<unsigned long> iVertex_to_iLocalVertex(nVertex);
  nLocalVertex = 0;
  for (unsigned long iVertex = 0; iVertex < nVertex; iVertex++) {
    const auto iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      nLocalVertex++;
    }
  }

  su2activematrix Buffer_Send_Coord(nLocalVertex, nDim);
  su2vector<unsigned long> Buffer_Send_GlobalPoint(nVertex);
  su2vector<unsigned long> Buffer_Send_nLinkedNodes(nVertex);
  su2vector<unsigned long> Buffer_Send_StartLinkedNodes(nVertex);
  unsigned long** Aux_Send_Map = new unsigned long*[nVertex];
  /*---modification made here---*/
  short* Buffer_Send_PerBound = new short[nVertex];
  su2vector<unsigned long> Buffer_Send_nSurroundNodes(nVertex);
  su2vector<unsigned long> Buffer_Send_StartSurroundNodes(nVertex);
  unsigned long** Aux_Send_Map_SurroundNodes = new unsigned long*[nVertex];

  su2vector<unsigned long> Buffer_Send_nLinkedElems(nVertex);
  su2vector<unsigned long> Buffer_Send_StartLinkedElems(nVertex);
  unsigned long** Aux_Send_Map_Elem = new unsigned long*[nVertex];
  // unsigned long **Aux_Send_Map_ElemNodes        = new unsigned long*[ nVertex ];

#ifdef HAVE_MPI
  int nProcessor = size, iRank;
  unsigned long iTmp2, tmp_index, tmp_index_2;
  /*---modification made here---*/
  unsigned long iTmp3, tmp_index_3, iTmp4, tmp_index_4;
#endif

  /*--- Copy coordinates and point to the auxiliar vector ---*/

  nGlobalVertex = 0;
  nLocalVertex = 0;
  nLocalLinkedNodes = 0;
  /*---modification made here---*/
  nLocalLinkedElems = 0;
  nLocalSurroundNodes = 0;

  unsigned short PeriodicBoundary, nPeriodic;
  nPeriodic = config[val_zone]->GetnMarker_Periodic();
  /*---modification made above---*/

  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    Buffer_Send_nLinkedNodes[iVertex] = 0;
    Aux_Send_Map[iVertex] = nullptr;
    /*---modification made here---*/
    Buffer_Send_PerBound[iVertex] = -1;
    Buffer_Send_nLinkedElems[iVertex] = 0;
    Aux_Send_Map_Elem[iVertex] = nullptr;
    Buffer_Send_nSurroundNodes[iVertex] = 0;
    Aux_Send_Map_SurroundNodes[iVertex] = nullptr;
    // Aux_Send_Map_ElemNodes[iVertex]   = nullptr;
    /*---modification made above---*/

    iPoint = geom->vertex[val_marker][iVertex]->GetNode();

    if (geom->nodes->GetDomain(iPoint)) {
      Buffer_Send_GlobalPoint[nLocalVertex] = geom->nodes->GetGlobalIndex(iPoint);

      for (iDim = 0; iDim < nDim; iDim++) Buffer_Send_Coord(nLocalVertex, iDim) = geom->nodes->GetCoord(iPoint, iDim);

      // compute the number of linked neighbours nNodes for iPoint
      nNodes = 0;
      nEdges = geom->nodes->GetnPoint(iPoint);

      for (jEdge = 0; jEdge < nEdges; jEdge++) {
        EdgeIndex = geom->nodes->GetEdge(iPoint, jEdge);

        if (iPoint == geom->edges->GetNode(EdgeIndex, 0))
          dPoint = geom->edges->GetNode(EdgeIndex, 1);
        else
          dPoint = geom->edges->GetNode(EdgeIndex, 0);
        // find the dPoint on the marker
        if (geom->nodes->GetVertex(dPoint, val_marker) != -1) nNodes++;
      }

      Buffer_Send_StartLinkedNodes[nLocalVertex] = nLocalLinkedNodes;
      Buffer_Send_nLinkedNodes[nLocalVertex] = nNodes;

      nLocalLinkedNodes += nNodes;

      // the same loop as before, but this time assign global index into auxiliary matrix Aux_Send_Map
      Aux_Send_Map[nLocalVertex] = new unsigned long[nNodes];
      nNodes = 0;

      for (jEdge = 0; jEdge < nEdges; jEdge++) {
        EdgeIndex = geom->nodes->GetEdge(iPoint, jEdge);

        if (iPoint == geom->edges->GetNode(EdgeIndex, 0))
          dPoint = geom->edges->GetNode(EdgeIndex, 1);
        else
          dPoint = geom->edges->GetNode(EdgeIndex, 0);

        if (geom->nodes->GetVertex(dPoint, val_marker) != -1) {
          Aux_Send_Map[nLocalVertex][nNodes] = geom->nodes->GetGlobalIndex(dPoint);
          nNodes++;
        }
      }

      /*---modification made here---*/
      // determine if this point belongs to periodic boundary, useful for turbomachinery mode
      // loop through all markers to find periodic boundary
      for (int jMarker = 0; jMarker < config[val_zone]->GetnMarker_All(); jMarker++) {
        if (config[val_zone]->GetMarker_All_KindBC(jMarker) != PERIODIC_BOUNDARY) continue;
        // determine if iVertex/iPoint is on this periodic boundary
        PeriodicBoundary = config[val_zone]->GetMarker_All_PerBound(jMarker);
        jVertex = geom->nodes->GetVertex(iPoint, jMarker);
        // cout<<"PeriodicBoundary: "<<PeriodicBoundary<<endl;
        //  we have two candidates for periodic boundary
        //  the index of PerBound starts from 1
        if ((jVertex != -1) && (PeriodicBoundary == (val_zone + 1))) {
          Buffer_Send_PerBound[nLocalVertex] = PeriodicBoundary;
        }
        if ((jVertex != -1) && (PeriodicBoundary == (val_zone + 1 + nPeriodic / 2))) {
          Buffer_Send_PerBound[nLocalVertex] = PeriodicBoundary;
        }
      }

      // get the number of surrounding nodes(share the same element) for iPoint
      nSurroundNodes = nNodes;
      nElems = geom->nodes->GetnElem(iPoint);
      // get the number of surrounding nodes for iPoint, redundant points exist.
      for (iElem = 0; iElem < nElems; iElem++) {
        mElem = geom->nodes->GetElem(iPoint, iElem);
        mElemNodes = geom->elem[mElem]->GetnNodes();
        for (jNode = 0; jNode < mElemNodes; jNode++) {
          jPoint = geom->elem[mElem]->GetNode(jNode);
          // cout <<"jPoint: "<<jPoint<<endl;
          if (geom->nodes->GetVertex(jPoint, val_marker) != -1 && jPoint != iPoint) {
            nSurroundNodes++;
            // cout<<"jPoint on marker: "<<jPoint<<endl;
            //  loop through all linked points and exclude them
            for (kNode = 0; kNode < nNodes; kNode++) {
              // cout<<"Neighbours of iVertex: "<<Aux_Send_Map[nLocalVertex][kNode]<<endl;
              if (geom->nodes->GetGlobalIndex(jPoint) == Aux_Send_Map[nLocalVertex][kNode]) {
                // cout<<"jPoint on edge: "<<jPoint<<endl;
                nSurroundNodes--;
                break;
              }
            }
          }
        }
      }

      Buffer_Send_StartSurroundNodes[nLocalVertex] = nLocalSurroundNodes;
      Buffer_Send_nSurroundNodes[nLocalVertex] = nSurroundNodes;
      nLocalSurroundNodes += nSurroundNodes;
      // the same loop as before, but this time assign global index into auxiliary matrix Aux_Send_Map
      // assign one more position for list, cause we first assign its value then we decide if it exists already.
      Aux_Send_Map_SurroundNodes[nLocalVertex] = new unsigned long[nSurroundNodes + 1];
      // cout <<"##########reconstruct 1.1"<<endl;
      //  assign linked points to list first
      for (kNode = 0; kNode < nNodes; kNode++) {
        Aux_Send_Map_SurroundNodes[nLocalVertex][kNode] = Aux_Send_Map[nLocalVertex][kNode];
      }

      nSurroundNodes = nNodes;

      // get the number of linked Elements nElems for iPoint
      nElems = geom->nodes->GetnElem(iPoint);
      // get the number of surrounding nodes for iPoint, redundant points exist.
      for (iElem = 0; iElem < nElems; iElem++) {
        mElem = geom->nodes->GetElem(iPoint, iElem);
        mElemNodes = geom->elem[mElem]->GetnNodes();
        for (jNode = 0; jNode < mElemNodes; jNode++) {
          jPoint = geom->elem[mElem]->GetNode(jNode);
          if (geom->nodes->GetVertex(jPoint, val_marker) != -1 && jPoint != iPoint) {
            // assign this point first, next step determine whether keep it or not.
            Aux_Send_Map_SurroundNodes[nLocalVertex][nSurroundNodes] = geom->nodes->GetGlobalIndex(jPoint);
            nSurroundNodes++;
            // loop through all linked points and exclude them
            for (kNode = 0; kNode < nNodes; kNode++) {
              if (geom->nodes->GetGlobalIndex(jPoint) == Aux_Send_Map[nLocalVertex][kNode]) {
                nSurroundNodes--;
                break;
              }
            }
          }
        }
      }
      /*---modification made above---*/

      /*---modification made here---*/

      // get the number of linked Elements nElems for iPoint
      nElems = geom->nodes->GetnElem(iPoint);

      Buffer_Send_StartLinkedElems[nLocalVertex] = nLocalLinkedElems;
      Buffer_Send_nLinkedElems[nLocalVertex] = nElems;
      nLocalLinkedElems += nElems;

      // loop local elem indices iElem, get global index and assign global index into auxiliary matrix Aux_Send_Map
      Aux_Send_Map_Elem[nLocalVertex] = new unsigned long[nElems];

      for (iElem = 0; iElem < nElems; iElem++) {
        unsigned long iElem_local = geom->nodes->GetElem(iPoint, iElem);
        // Aux_Send_Map_Elem[nLocalVertex][iElem] = geom->nodes->GetElem(iPoint, iElem);
        Aux_Send_Map_Elem[nLocalVertex][iElem] = geom->elem[iElem_local]->GetGlobalIndex();
        ;
      }
      /*---modification made above---*/

      nLocalVertex++;
    }
  }

  su2vector<unsigned long> Buffer_Send_LinkedNodes(nLocalLinkedNodes);

  nLocalLinkedNodes = 0;

  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    for (jEdge = 0; jEdge < Buffer_Send_nLinkedNodes[iVertex]; jEdge++) {
      Buffer_Send_LinkedNodes[nLocalLinkedNodes] = Aux_Send_Map[iVertex][jEdge];
      nLocalLinkedNodes++;
    }
  }

  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    if (Aux_Send_Map[iVertex] != nullptr) delete[] Aux_Send_Map[iVertex];
  }
  delete[] Aux_Send_Map;
  Aux_Send_Map = nullptr;

  /*---modification made here---*/
  su2vector<unsigned long> Buffer_Send_SurroundNodes(nLocalSurroundNodes);
  nLocalSurroundNodes = 0;

  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    for (kNode = 0; kNode < Buffer_Send_nSurroundNodes[iVertex]; kNode++) {
      Buffer_Send_SurroundNodes[nLocalSurroundNodes] = Aux_Send_Map_SurroundNodes[iVertex][kNode];
      nLocalSurroundNodes++;
    }
  }

  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    if (Aux_Send_Map_SurroundNodes[iVertex] != nullptr) delete[] Aux_Send_Map_SurroundNodes[iVertex];
  }
  delete[] Aux_Send_Map_SurroundNodes;
  Aux_Send_Map_SurroundNodes = nullptr;

  su2vector<unsigned long> Buffer_Send_LinkedElems(nLocalLinkedElems);
  nLocalLinkedElems = 0;

  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    for (iElem = 0; iElem < Buffer_Send_nLinkedElems[iVertex]; iElem++) {
      Buffer_Send_LinkedElems[nLocalLinkedElems] = Aux_Send_Map_Elem[iVertex][iElem];
      nLocalLinkedElems++;
    }
  }

  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    if (Aux_Send_Map_Elem[iVertex] != nullptr) delete[] Aux_Send_Map_Elem[iVertex];
  }
  delete[] Aux_Send_Map_Elem;
  Aux_Send_Map_Elem = nullptr;
  /*---modification made above---*/

  /*--- Reconstruct  boundary by gathering data from all ranks ---*/

  SU2_MPI::Allreduce(&nLocalVertex, &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nLocalLinkedNodes, &nGlobalLinkedNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  /*---modification made here---*/
  SU2_MPI::Allreduce(&nLocalLinkedElems, &nGlobalLinkedElems, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nLocalSurroundNodes, &nGlobalSurroundNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  Buffer_Receive_Coord.resize(nGlobalVertex, nDim);
  Buffer_Receive_GlobalPoint.resize(nGlobalVertex);
  Buffer_Receive_Proc.resize(nGlobalVertex);

  Buffer_Receive_nLinkedNodes.resize(nGlobalVertex);
  Buffer_Receive_LinkedNodes.resize(nGlobalLinkedNodes);
  Buffer_Receive_StartLinkedNodes.resize(nGlobalVertex);
  /*---modification made here---*/
  Buffer_Receive_PerBound = new short[nGlobalVertex];
  Buffer_Receive_nLinkedElems.resize(nGlobalVertex);
  Buffer_Receive_LinkedElems.resize(nGlobalLinkedElems);
  Buffer_Receive_StartLinkedElems.resize(nGlobalVertex);
  Buffer_Receive_nSurroundNodes.resize(nGlobalVertex);
  Buffer_Receive_SurroundNodes.resize(nGlobalSurroundNodes);
  Buffer_Receive_StartSurroundNodes.resize(nGlobalVertex);
  /*---modification made above---*/

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) {
    for (unsigned long iVertex = 0; iVertex < nLocalVertex; iVertex++)
      for (unsigned long iDim = 0; iDim < nDim; iDim++)
        Buffer_Receive_Coord(iVertex, iDim) = Buffer_Send_Coord(iVertex, iDim);

    for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
      Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Proc[iVertex] = MASTER_NODE;
      Buffer_Receive_nLinkedNodes[iVertex] = Buffer_Send_nLinkedNodes[iVertex];
      Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
      /*---modification made here---*/
      Buffer_Receive_PerBound[iVertex] = Buffer_Send_PerBound[iVertex];
      Buffer_Receive_nLinkedElems[iVertex] = Buffer_Send_nLinkedElems[iVertex];
      Buffer_Receive_StartLinkedElems[iVertex] = Buffer_Send_StartLinkedElems[iVertex];
      Buffer_Receive_nSurroundNodes[iVertex] = Buffer_Send_nSurroundNodes[iVertex];
      Buffer_Receive_StartSurroundNodes[iVertex] = Buffer_Send_StartSurroundNodes[iVertex];
      /*---modification made above---*/
    }

    for (iVertex = 0; iVertex < nLocalLinkedNodes; iVertex++)
      Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];

    /*---modification made here---*/
    for (iVertex = 0; iVertex < nLocalLinkedElems; iVertex++)
      Buffer_Receive_LinkedElems[iVertex] = Buffer_Send_LinkedElems[iVertex];

    for (iVertex = 0; iVertex < nLocalSurroundNodes; iVertex++)
      Buffer_Receive_SurroundNodes[iVertex] = Buffer_Send_SurroundNodes[iVertex];
    /*---modification made above---*/

    tmp_index = nLocalVertex;
    tmp_index_2 = nLocalLinkedNodes;
    tmp_index_3 = nLocalLinkedElems;
    tmp_index_4 = nLocalSurroundNodes;

    for (iRank = 1; iRank < nProcessor; iRank++) {
      SU2_MPI::Recv(&iTmp2, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_LinkedNodes[tmp_index_2], iTmp2, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);

      SU2_MPI::Recv(&iTmp, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(Buffer_Receive_Coord[tmp_index], nDim * iTmp, MPI_DOUBLE, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);

      SU2_MPI::Recv(&Buffer_Receive_GlobalPoint[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_nLinkedNodes[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_StartLinkedNodes[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);

      /*---modification made here---*/
      SU2_MPI::Recv(&Buffer_Receive_PerBound[tmp_index], iTmp, MPI_SHORT, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);

      SU2_MPI::Recv(&iTmp4, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_SurroundNodes[tmp_index_4], iTmp4, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&iTmp3, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_LinkedElems[tmp_index_3], iTmp3, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_nLinkedElems[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_StartLinkedElems[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_nSurroundNodes[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(),
                    MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_StartSurroundNodes[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1,
                    SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      /*---modification made above---*/

      for (iVertex = 0; iVertex < iTmp; iVertex++) {
        Buffer_Receive_Proc[tmp_index + iVertex] = iRank;
        // set offset for StartLinkedNodes/ELems
        Buffer_Receive_StartLinkedNodes[tmp_index + iVertex] += tmp_index_2;
        /*---modification made here---*/
        Buffer_Receive_StartLinkedElems[tmp_index + iVertex] += tmp_index_3;
        Buffer_Receive_StartSurroundNodes[tmp_index + iVertex] += tmp_index_4;
        /*---modification made above---*/
      }

      tmp_index += iTmp;
      tmp_index_2 += iTmp2;
      /*---modification made here---*/
      tmp_index_3 += iTmp3;
      tmp_index_4 += iTmp4;
      /*---modification made above---*/
    }
  } else {
    SU2_MPI::Send(&nLocalLinkedNodes, 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_LinkedNodes.data(), nLocalLinkedNodes, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(&nLocalVertex, 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_Coord.data(), Buffer_Send_Coord.size(), MPI_DOUBLE, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(Buffer_Send_GlobalPoint.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_nLinkedNodes.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_StartLinkedNodes.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());

    /*---modification made here---*/
    SU2_MPI::Send(Buffer_Send_PerBound, nLocalVertex, MPI_SHORT, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(&nLocalSurroundNodes, 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_SurroundNodes.data(), nLocalSurroundNodes, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(&nLocalLinkedElems, 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_LinkedElems.data(), nLocalLinkedElems, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_nLinkedElems.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_StartLinkedElems.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_nSurroundNodes.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_StartSurroundNodes.data(), nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    /*---modification made above---*/
  }
#else
  for (unsigned long iVertex = 0; iVertex < nLocalVertex; iVertex++)
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Buffer_Receive_Coord(iVertex, iDim) = Buffer_Send_Coord(iVertex, iDim);

  for (iVertex = 0; iVertex < nGlobalVertex; iVertex++) {
    Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
    Buffer_Receive_Proc[iVertex] = MASTER_NODE;
    Buffer_Receive_nLinkedNodes[iVertex] = Buffer_Send_nLinkedNodes[iVertex];
    Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];

    /*---modification made here---*/
    Buffer_Receive_PerBound[iVertex] = Buffer_Send_PerBound[iVertex];

    Buffer_Receive_nLinkedElems[iVertex] = Buffer_Send_nLinkedElems[iVertex];
    Buffer_Receive_StartLinkedElems[iVertex] = Buffer_Send_StartLinkedElems[iVertex];
    Buffer_Receive_nSurroundNodes[iVertex] = Buffer_Send_nSurroundNodes[iVertex];
    Buffer_Receive_StartSurroundNodes[iVertex] = Buffer_Send_StartSurroundNodes[iVertex];
    /*---modification made above---*/
  }

  for (iVertex = 0; iVertex < nGlobalLinkedNodes; iVertex++)
    Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];
  /*---modification made here---*/
  for (iVertex = 0; iVertex < nGlobalLinkedElems; iVertex++)
    Buffer_Receive_LinkedElems[iVertex] = Buffer_Send_LinkedElems[iVertex];
  for (iVertex = 0; iVertex < nGlobalSurroundNodes; iVertex++)
    Buffer_Receive_SurroundNodes[iVertex] = Buffer_Send_SurroundNodes[iVertex];
    /*---modification made above---*/
#endif

  if (rank == MASTER_NODE) {
    for (iVertex = 0; iVertex < nGlobalVertex; iVertex++) {
      count = 0;
      uptr = &Buffer_Receive_LinkedNodes[Buffer_Receive_StartLinkedNodes[iVertex]];
      // option for linkedNodes list
      // replace the global index with the reconstructed list index, i.e. nGlobaVertex index
      for (jVertex = 0; jVertex < Buffer_Receive_nLinkedNodes[iVertex]; jVertex++) {
        iTmp = uptr[jVertex];
        for (kVertex = 0; kVertex < nGlobalVertex; kVertex++) {
          if (Buffer_Receive_GlobalPoint[kVertex] == long(iTmp)) {
            uptr[jVertex] = kVertex;
            count++;
            break;
          }
        }
        // if we do not find this neighbour of iVertex in reconstructed list nGlobaVertex
        // remove it from the linked list, but the start list is not changed
        if (count != (jVertex + 1)) {
          for (kVertex = jVertex; kVertex < Buffer_Receive_nLinkedNodes[iVertex] - 1; kVertex++) {
            uptr[kVertex] = uptr[kVertex + 1];
          }
          Buffer_Receive_nLinkedNodes[iVertex]--;
          jVertex--;
        }
      }
    }

    // option for surroundNodes list
    // replace the global index with the reconstructed list index, i.e. nGlobaVertex index
    for (iVertex = 0; iVertex < nGlobalVertex; iVertex++) {
      count = 0;
      uptr = &Buffer_Receive_SurroundNodes[Buffer_Receive_StartSurroundNodes[iVertex]];
      // make direct connection using new index nGlobalVertex

      for (jVertex = 0; jVertex < Buffer_Receive_nSurroundNodes[iVertex]; jVertex++) {
        iTmp = uptr[jVertex];
        for (kVertex = 0; kVertex < nGlobalVertex; kVertex++) {
          if (Buffer_Receive_GlobalPoint[kVertex] == long(iTmp)) {
            uptr[jVertex] = kVertex;
            count++;
            break;
          }
        }

        // if we do not find this neighbour of iVertex in reconstructed list nGlobaVertex
        // remove it from the surrounding list, but the start list is not changed

        if (count != (jVertex + 1)) {
          for (kVertex = jVertex; kVertex < Buffer_Receive_nSurroundNodes[iVertex] - 1; kVertex++) {
            uptr[kVertex] = uptr[kVertex + 1];
          }
          Buffer_Receive_nSurroundNodes[iVertex]--;
          jVertex--;
        }
      }
      // auto    iPoint = geom->vertex[val_marker][iVertex]->GetNode();
      // auto dPoint = geom->nodes->GetGlobalIndex(iPoint);
      // if(iVertex == 0)
      // cout <<"iVertex: "<<iVertex<<", GlobalIdx: "<<dPoint<<", nSurroundNodes:
      // "<<Buffer_Receive_nSurroundNodes[iVertex]<<endl;
    }
  }
  SU2_MPI::Bcast(Buffer_Receive_GlobalPoint.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_Coord.data(), Buffer_Receive_Coord.size(), MPI_DOUBLE, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_Proc.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());

  SU2_MPI::Bcast(Buffer_Receive_nLinkedNodes.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_StartLinkedNodes.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_LinkedNodes.data(), nGlobalLinkedNodes, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());

  /*---modification made here---*/
  SU2_MPI::Bcast(Buffer_Receive_PerBound, nGlobalVertex, MPI_SHORT, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_nLinkedElems.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_StartLinkedElems.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_LinkedElems.data(), nGlobalLinkedElems, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_nSurroundNodes.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_StartSurroundNodes.data(), nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_SurroundNodes.data(), nGlobalSurroundNodes, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  /*---modification made above---*/
}

void CInterpolator::Collect_BundaryPitchInfo(const CConfig* const* config, int val_marker_donor,
                                             int val_marker_target) {
  unsigned short iSpan;
  unsigned long TimeIter;
  TimeIter = config[donorZone]->GetTimeIter();
  if (TimeIter > 0) {
    /* this way for sure get nSpanDonor for each rank */
    nSpanDonor = config[donorZone]->GetnSpanWiseSections();
    nSpanTarget = config[targetZone]->GetnSpanWiseSections();
  }

  Pitch = new su2double[nSpanDonor];
  InitMaxAng_donor = new su2double[nSpanDonor];
  InitMinAng_donor = new su2double[nSpanDonor];
  MaxAng_donor = new su2double[nSpanDonor];
  MinAng_donor = new su2double[nSpanDonor];
  InitMaxAng_target = new su2double[nSpanTarget];
  InitMinAng_target = new su2double[nSpanTarget];
  MaxAng_target = new su2double[nSpanTarget];
  MinAng_target = new su2double[nSpanTarget];
  SpanValuesDonor = new su2double[nSpanDonor];
  SpanValuesTarget = new su2double[nSpanTarget];
  SpanValuesDonor_shadow = new su2double[nSpanDonor];
  SpanValuesTarget_shadow = new su2double[nSpanTarget];

  for (iSpan = 0; iSpan < nSpanDonor; iSpan++) {
    Pitch[iSpan] = -1.0;
    InitMaxAng_donor[iSpan] = -1.0;
    InitMinAng_donor[iSpan] = -1.0;
    SpanValuesDonor_shadow[iSpan] = -1.0;
  }
  for (iSpan = 0; iSpan < nSpanTarget; iSpan++) {
    InitMaxAng_target[iSpan] = -1.0;
    InitMinAng_target[iSpan] = -1.0;
    SpanValuesTarget_shadow[iSpan] = -1.0;
  }

  /*--- get target side max/min angle locally on specific rank ---*/
  /* on the rank that has markDonor>=0 */
  if (TimeIter > 0) {
    if (val_marker_donor != -1) {
      SpanValuesDonor =
          donor_geometry->GetSpanWiseValue(config[donorZone]->GetMarker_All_TurbomachineryFlag(val_marker_donor));
      for (iSpan = 0; iSpan < nSpanDonor; iSpan++) {
        SpanValuesDonor_shadow[iSpan] = SpanValuesDonor[iSpan];
        InitMaxAng_donor[iSpan] = donor_geometry->GetMaxAngularCoord(val_marker_donor, iSpan);
        InitMinAng_donor[iSpan] = donor_geometry->GetMinAngularCoord(val_marker_donor, iSpan);
        Pitch[iSpan] = InitMaxAng_donor[iSpan] - InitMinAng_donor[iSpan];
      }
    }
    /*--- get donor side max/min angle locally on specific rank ---*/
    /* on the rank that has markTarget>=0 */
    if (val_marker_target != -1) {
      SpanValuesTarget =
          target_geometry->GetSpanWiseValue(config[targetZone]->GetMarker_All_TurbomachineryFlag(val_marker_target));
      for (iSpan = 0; iSpan < nSpanTarget; iSpan++) {
        SpanValuesTarget_shadow[iSpan] = SpanValuesTarget[iSpan];
        InitMaxAng_target[iSpan] = target_geometry->GetMaxAngularCoord(val_marker_target, iSpan);
        InitMinAng_target[iSpan] = target_geometry->GetMinAngularCoord(val_marker_target, iSpan);

        /*---currently we assume each zone has same pitch angle ---*/
        // Pitch = MaxAng_target - MinAng_target;
      }
    }
  }
/*--- collect the max/min/pitch angle info and redistribute to each rank ---*/
#ifdef HAVE_MPI

  int iSize;
  bool DonorFound = false, TargetFound = false;
  su2double *BuffPitch = NULL, *BuffInitMaxAng_donor = NULL, *BuffInitMinAng_donor = NULL,
            *BuffInitMaxAng_target = NULL, *BuffInitMinAng_target = NULL, *BuffSpanValuesDonor = NULL,
            *BuffSpanValuesTarget = NULL;
  int nSpanSize_donor, nSpanSize_target;

  nSpanSize_donor = size * nSpanDonor;
  nSpanSize_target = size * nSpanTarget;
  BuffPitch = new su2double[nSpanSize_donor];
  BuffInitMaxAng_donor = new su2double[nSpanSize_donor];
  BuffInitMinAng_donor = new su2double[nSpanSize_donor];
  BuffInitMaxAng_target = new su2double[nSpanSize_target];
  BuffInitMinAng_target = new su2double[nSpanSize_target];
  BuffSpanValuesDonor = new su2double[nSpanSize_donor];
  BuffSpanValuesTarget = new su2double[nSpanSize_target];

  /*--- set default value ---*/
  for (iSpan = 0; iSpan < nSpanSize_donor; iSpan++) {
    BuffPitch[iSpan] = -1.0;
    BuffInitMaxAng_donor[iSpan] = -1.0;
    BuffInitMinAng_donor[iSpan] = -1.0;
    BuffSpanValuesDonor[iSpan] = -1.0;
  }
  for (iSpan = 0; iSpan < nSpanSize_target; iSpan++) {
    BuffInitMaxAng_target[iSpan] = -1.0;
    BuffInitMinAng_target[iSpan] = -1.0;
    BuffSpanValuesTarget[iSpan] = -1.0;
  }

  SU2_MPI::Allgather(Pitch, nSpanDonor, MPI_DOUBLE, BuffPitch, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(InitMaxAng_donor, nSpanDonor, MPI_DOUBLE, BuffInitMaxAng_donor, nSpanDonor, MPI_DOUBLE,
                     MPI_COMM_WORLD);
  SU2_MPI::Allgather(InitMinAng_donor, nSpanDonor, MPI_DOUBLE, BuffInitMinAng_donor, nSpanDonor, MPI_DOUBLE,
                     MPI_COMM_WORLD);
  SU2_MPI::Allgather(InitMaxAng_target, nSpanTarget, MPI_DOUBLE, BuffInitMaxAng_target, nSpanTarget, MPI_DOUBLE,
                     MPI_COMM_WORLD);
  SU2_MPI::Allgather(InitMinAng_target, nSpanTarget, MPI_DOUBLE, BuffInitMinAng_target, nSpanTarget, MPI_DOUBLE,
                     MPI_COMM_WORLD);
  SU2_MPI::Allgather(SpanValuesDonor_shadow, nSpanDonor, MPI_DOUBLE, BuffSpanValuesDonor, nSpanDonor, MPI_DOUBLE,
                     MPI_COMM_WORLD);
  SU2_MPI::Allgather(SpanValuesTarget_shadow, nSpanTarget, MPI_DOUBLE, BuffSpanValuesTarget, nSpanTarget, MPI_DOUBLE,
                     MPI_COMM_WORLD);

  for (iSpan = 0; iSpan < nSpanDonor; iSpan++) {
    Pitch[iSpan] = -1.0;
    InitMaxAng_donor[iSpan] = -1.0;
    InitMinAng_donor[iSpan] = -1.0;
    SpanValuesDonor_shadow[iSpan] = -1.0;
  }
  for (iSpan = 0; iSpan < nSpanTarget; iSpan++) {
    InitMaxAng_target[iSpan] = -1.0;
    InitMinAng_target[iSpan] = -1.0;
    SpanValuesTarget_shadow[iSpan] = -1.0;
  }

  for (iSize = 0; iSize < size; iSize++) {
    /*--- more than one rank can have correct boundary mark and each of them contain whole
    data set we need. if one rank is found, download data to every rank and then just break the loop ---*/
    DonorFound = false, TargetFound = false;
    if (BuffPitch[nSpanDonor * iSize] > 0.0) {
      for (iSpan = 0; iSpan < nSpanDonor; iSpan++) {
        Pitch[iSpan] = BuffPitch[nSpanDonor * iSize + iSpan];
        InitMaxAng_donor[iSpan] = BuffInitMaxAng_donor[nSpanDonor * iSize + iSpan];
        InitMinAng_donor[iSpan] = BuffInitMinAng_donor[nSpanDonor * iSize + iSpan];
        // InitMaxAng_target[iSpan]			= BuffInitMaxAng_target[nSpanDonor*iSize + iSpan];
        // InitMinAng_target[iSpan] 			= BuffInitMinAng_target[nSpanDonor*iSize + iSpan];
        SpanValuesDonor_shadow[iSpan] = BuffSpanValuesDonor[nSpanDonor * iSize + iSpan];
        // SpanValuesTarget_shadow[iSpan]			= BuffSpanValuesTarget[nSpanDonor*iSize + iSpan];
      }
      DonorFound = true;
    }

    if (BuffInitMaxAng_target[nSpanTarget * iSize] > 0.0) {
      for (iSpan = 0; iSpan < nSpanTarget; iSpan++) {
        InitMaxAng_target[iSpan] = BuffInitMaxAng_target[nSpanTarget * iSize + iSpan];
        InitMinAng_target[iSpan] = BuffInitMinAng_target[nSpanTarget * iSize + iSpan];
        SpanValuesTarget_shadow[iSpan] = BuffSpanValuesTarget[nSpanTarget * iSize + iSpan];
      }
      TargetFound = true;
    }

    if (DonorFound && TargetFound) break;
  }
  delete[] BuffPitch;
  delete[] BuffInitMaxAng_donor;
  delete[] BuffInitMinAng_donor;
  delete[] BuffInitMaxAng_target;
  delete[] BuffInitMinAng_target;
  delete[] BuffSpanValuesDonor;
  delete[] BuffSpanValuesTarget;

#endif
}
