/*!
 * \file CInterpolator.cpp
 * \brief Definition of the base class for interface interpolation.
 * \author H. Kline
 * \version 7.1.1 "Blackbird"
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

#include "../../include/interface_interpolation/CInterpolator.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"


CInterpolator::CInterpolator(CGeometry ****geometry_container, const CConfig* const* config,
                             unsigned int iZone, unsigned int jZone) :
  rank(SU2_MPI::GetRank()),
  size(SU2_MPI::GetSize()),
  donorZone(iZone),
  targetZone(jZone),
  Geometry(geometry_container),
  donor_geometry(geometry_container[iZone][INST_0][MESH_0]),
  target_geometry(geometry_container[jZone][INST_0][MESH_0]) {
}

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
  for (auto iInter = 0; iInter < (donor->GetMarker_n_ZoneInterface()/2); iInter++) {
    if (CheckInterfaceBoundary(donor->FindInterfaceMarker(iInter), target->FindInterfaceMarker(iInter)))
      return true;
  }
  return false;
}

void CInterpolator::Determine_ArraySize(int markDonor, int markTarget,
                    unsigned long nVertexDonor, unsigned short nDim) {

  /*--- Count donor vertices. ---*/
  auto nLocalVertex_Donor = 0ul;
  for (auto iVertex = 0ul; iVertex < nVertexDonor; iVertex++) {
    auto iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    nLocalVertex_Donor += donor_geometry->nodes->GetDomain(iPointDonor);
  }

  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;

  /*--- Send Interface vertex information --*/
  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG,
                     Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
}

void CInterpolator::Collect_VertexInfo(int markDonor, int markTarget,
                    unsigned long nVertexDonor, unsigned short nDim) {

  unsigned long iVertex;
  unsigned short iDim;

  for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) Buffer_Send_GlobalPoint[iVertex] = -1;

  for (iVertex = 0; iVertex < MaxLocalVertex_Donor*nDim; iVertex++) Buffer_Send_Coord[iVertex] = 0.0;

  /*--- Copy coordinates and point to the auxiliar vector --*/
  auto iLocalVertexDonor = 0ul;

  for (iVertex = 0; iVertex < nVertexDonor; iVertex++) {
    auto iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    if (donor_geometry->nodes->GetDomain(iPointDonor)) {
      Buffer_Send_GlobalPoint[iLocalVertexDonor] = donor_geometry->nodes->GetGlobalIndex(iPointDonor);
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iLocalVertexDonor*nDim+iDim] = donor_geometry->nodes->GetCoord(iPointDonor, iDim);
      iLocalVertexDonor++;
    }
  }
  auto nBuffer_Coord = MaxLocalVertex_Donor*nDim;
  auto nBuffer_Point = MaxLocalVertex_Donor;

  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE,
                     Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_LONG,
                     Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_LONG, SU2_MPI::GetComm());
}

unsigned long CInterpolator::Collect_ElementInfo(int markDonor, unsigned short nDim, bool compress,
                             vector<unsigned long>& allNumElem, vector<unsigned short>& numNodes,
                             su2matrix<long>& idxNodes) const {

  const auto maxElemNodes = (nDim == 2u)? 2u : 4u; // line and quad respectively

  unsigned long nElemDonor = 0;
  if (markDonor != -1) nElemDonor = donor_geometry->GetnElem_Bound(markDonor);

  allNumElem.resize(size);
  SU2_MPI::Allgather(&nElemDonor, 1, MPI_UNSIGNED_LONG, allNumElem.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  auto nMaxElemDonor = *max_element(allNumElem.begin(), allNumElem.end());

  vector<unsigned short> bufferSendNum(nMaxElemDonor);
  su2matrix<long> bufferSendIdx(nMaxElemDonor, maxElemNodes);

  numNodes.resize(nMaxElemDonor*size);
  idxNodes.resize(nMaxElemDonor*size, maxElemNodes);

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

  SU2_MPI::Allgather(bufferSendNum.data(), bufferSendNum.size(), MPI_UNSIGNED_SHORT,
                     numNodes.data(),      bufferSendNum.size(), MPI_UNSIGNED_SHORT, SU2_MPI::GetComm());
  SU2_MPI::Allgather(bufferSendIdx.data(), bufferSendIdx.size(), MPI_LONG,
                     idxNodes.data(),      bufferSendIdx.size(), MPI_LONG, SU2_MPI::GetComm());

  if (!compress)
    return accumulate(allNumElem.begin(), allNumElem.end(), 0ul);

  /*--- Compress the information (overlapping copy do not use memcpy). ---*/

  unsigned long dstIdx = 0;
  for (int iProcessor = 0; iProcessor < size; ++iProcessor) {
    auto srcOffset = iProcessor * nMaxElemDonor;
    for (auto idx = 0u; idx < allNumElem[iProcessor]; ++idx) {
      numNodes[dstIdx] = numNodes[srcOffset+idx];
      for (auto iNode = 0u; iNode < maxElemNodes; ++iNode)
        idxNodes(dstIdx, iNode) = idxNodes(srcOffset+idx, iNode);
      ++dstIdx;
    }
  }

  return dstIdx;
}

void CInterpolator::ReconstructBoundary(unsigned long val_zone, int val_marker){

  CGeometry *geom = Geometry[val_zone][INST_0][MESH_0];

  unsigned long iVertex, kVertex;

  unsigned long *uptr, nVertex, nElems, iDim, nDim, iPoint;

  unsigned long nGlobalLinkedNodes, nLocalVertex, nLocalLinkedNodes;

  nDim = geom->GetnDim();

  /*--- If this zone has no parts of the marker, it will not participate
   * in the first part of this function (collection). ---*/
  if( val_marker != -1 ){
    nVertex  = geom->GetnVertex(  val_marker  );
    nElems   = geom->GetnElem_Bound(val_marker);
  } else {
    nVertex  = 0;
    nElems   = 0;
  }

  /*--- Get the number of domain vertices on the marker, and a mapping
  * (iVertex) -> (iLocalVertex, non-domain points being ignored). ---*/
  unsigned long* iVertex_to_iLocalVertex = new unsigned long[ nVertex ];
  nLocalVertex = 0;
  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)){
      iVertex_to_iLocalVertex[iVertex] = nLocalVertex;
      nLocalVertex++;
    }  else {
      iVertex_to_iLocalVertex[iVertex] = numeric_limits<unsigned long>::max();
    }
  }

  // coordinates of all domain vertices on the marker
  su2double *Buffer_Send_Coord           = new su2double     [ nLocalVertex * nDim ];
  // global point IDs of all domain vertices on the marker
  long *Buffer_Send_GlobalPoint = new long [ nVertex ];

  // Assign to each domain vertex on the marker, identified by local point ID,
  // a set of surface-neighbor vertices on the marker, identified by global point ID.
  map<unsigned long, forward_list<unsigned long>*> neighbors;

  /*--- Define or initialize them. ---*/
  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      unsigned long iLocalVertex = iVertex_to_iLocalVertex[iVertex];
      Buffer_Send_GlobalPoint[iLocalVertex] = (long) geom->nodes->GetGlobalIndex(iPoint);
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iLocalVertex*nDim+iDim] = geom->nodes->GetCoord(iPoint, iDim);
      neighbors.insert(pair<unsigned long, forward_list<unsigned long>*>(iPoint, new forward_list<unsigned long>));
    }
  }

  /*--- Define the neighbors map. ---*/
  for(unsigned long iElem=0; iElem < nElems; iElem++){
    CPrimalGrid* elem = geom->bound[val_marker][iElem];
    for(unsigned short iNode=0; iNode<elem->GetnNodes(); iNode++){
      iPoint = elem->GetNode(iNode);
      if (geom->nodes->GetDomain(iPoint)) {
        forward_list<unsigned long>* neighb = neighbors.at(iPoint);
        for(unsigned short iNeighbor=0; iNeighbor<elem->GetnNeighbor_Nodes(iNode); iNeighbor++){
          unsigned long jPoint = elem->GetNode( elem->GetNeighbor_Nodes(iNode,iNeighbor) );
          unsigned long jPoint_global = geom->nodes->GetGlobalIndex(jPoint);
          if( std::find(std::begin(*neighb), std::end(*neighb), jPoint_global) == std::end(*neighb) ){
            neighb->emplace_front( jPoint_global );
          }
        }
      }
    }
  }

  // numbers of surface-neighbors of all domain vertices on the marker
  unsigned long *Buffer_Send_nLinkedNodes       = new unsigned long [ nLocalVertex ];
  // cumsum of Buffer_Send_nLinkedNodes
  unsigned long *Buffer_Send_StartLinkedNodes   = new unsigned long [ nLocalVertex ];
  nLocalLinkedNodes = 0;
  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      unsigned long iLocalVertex = iVertex_to_iLocalVertex[iVertex];
      Buffer_Send_nLinkedNodes[iLocalVertex] = std::count_if(
        std::begin(*neighbors[iPoint]), std::end(*neighbors[iPoint]),
        [](unsigned long i){return true;} );
      Buffer_Send_StartLinkedNodes[iLocalVertex] = nLocalLinkedNodes;
      nLocalLinkedNodes += Buffer_Send_nLinkedNodes[iLocalVertex];
    }
  }
  // global point IDs of surface-neighbors of all domain vertices on the marker
  unsigned long *Buffer_Send_LinkedNodes        = new unsigned long [ nLocalLinkedNodes ];
  unsigned long index = 0;
  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    if (geom->nodes->GetDomain(iPoint)) {
      for(unsigned long jPoint_global : *(neighbors[iPoint])){
        Buffer_Send_LinkedNodes[index] = jPoint_global;
        index++;
      }
      delete neighbors[iPoint];
      neighbors[iPoint]=nullptr;
    }
  }

  delete[] iVertex_to_iLocalVertex;


  /*--- Now these arrays of all processes must be joined to a single/global arrays. For this,
   * the entries of StartLinkedNodes must be shifted.
   * Furthermore, the global point IDs in LinkedNodes are replaced by global vertex IDs.
   * For this, the master process collects the data from all processes, joins them and broadcasts them again. ---*/

  /*--- Allocate global arrays. ---*/
  SU2_MPI::Allreduce(     &nLocalVertex,      &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nLocalLinkedNodes, &nGlobalLinkedNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  Buffer_Receive_Coord       = new su2double    [ nGlobalVertex * nDim ];
  Buffer_Receive_GlobalPoint = new long[ nGlobalVertex ];
  Buffer_Receive_Proc        = new unsigned long[ nGlobalVertex ];

  Buffer_Receive_nLinkedNodes     = new unsigned long[ nGlobalVertex ];
  Buffer_Receive_LinkedNodes      = new unsigned long[ nGlobalLinkedNodes   ];
  Buffer_Receive_StartLinkedNodes = new unsigned long[ nGlobalVertex ];

  /*--- Master process gathers data from all ranks ---*/
#ifdef HAVE_MPI
  int nProcessor = size, iRank;
  unsigned long received_nLocalLinkedNodes, received_nLocalVertex;

  if (rank == MASTER_NODE){
    /*--- "Receive" from master process, i.e. copy. ---*/
    for (iVertex = 0; iVertex < nDim*nLocalVertex; iVertex++)
      Buffer_Receive_Coord[iVertex]  = Buffer_Send_Coord[iVertex];

    for (iVertex = 0; iVertex < nLocalVertex; iVertex++){
      Buffer_Receive_GlobalPoint[iVertex]      = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Proc[iVertex]             = MASTER_NODE;
      Buffer_Receive_nLinkedNodes[iVertex]     = Buffer_Send_nLinkedNodes[iVertex];
      Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
    }

    for (iVertex = 0; iVertex < nLocalLinkedNodes; iVertex++)
      Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];

    unsigned long received_nLocalVertex_sum   = nLocalVertex;
    unsigned long received_nLocalLinkedNodes_sum = nLocalLinkedNodes;

    /*--- Receive from other processes to the appropriate position in the buffers and shift StartLinkedNodes indices. ---*/
    for(iRank = 1; iRank < nProcessor; iRank++){

      SU2_MPI::Recv(&received_nLocalLinkedNodes, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_LinkedNodes[received_nLocalLinkedNodes_sum], received_nLocalLinkedNodes, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      SU2_MPI::Recv(&received_nLocalVertex, 1, MPI_UNSIGNED_LONG, iRank, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_Coord[received_nLocalVertex_sum*nDim], nDim*received_nLocalVertex,        MPI_DOUBLE, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      SU2_MPI::Recv(     &Buffer_Receive_GlobalPoint[received_nLocalVertex_sum], received_nLocalVertex, MPI_LONG, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(    &Buffer_Receive_nLinkedNodes[received_nLocalVertex_sum], received_nLocalVertex, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_StartLinkedNodes[received_nLocalVertex_sum], received_nLocalVertex, MPI_UNSIGNED_LONG, iRank, 1, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      for (iVertex = 0; iVertex < received_nLocalVertex; iVertex++){
        Buffer_Receive_Proc[ received_nLocalVertex_sum + iVertex ] = iRank;
        Buffer_Receive_StartLinkedNodes[ received_nLocalVertex_sum + iVertex ] += received_nLocalLinkedNodes_sum;
      }

      received_nLocalVertex_sum   += received_nLocalVertex;
      received_nLocalLinkedNodes_sum += received_nLocalLinkedNodes;
    }
  }
  else{
    SU2_MPI::Send(     &nLocalLinkedNodes,                 1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_LinkedNodes, nLocalLinkedNodes, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(    &nLocalVertex,                   1, MPI_UNSIGNED_LONG, 0, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_Coord, nDim * nLocalVertex,        MPI_DOUBLE, 0, 1, SU2_MPI::GetComm());

    SU2_MPI::Send(     Buffer_Send_GlobalPoint, nLocalVertex, MPI_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(    Buffer_Send_nLinkedNodes, nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(Buffer_Send_StartLinkedNodes, nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, SU2_MPI::GetComm());
  }
#else
  for (iVertex = 0; iVertex < nDim * nGlobalVertex; iVertex++)
    Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];

  for (iVertex = 0; iVertex < nGlobalVertex; iVertex++){
    Buffer_Receive_GlobalPoint[iVertex]      = Buffer_Send_GlobalPoint[iVertex];
    Buffer_Receive_Proc[iVertex]             = MASTER_NODE;
    Buffer_Receive_nLinkedNodes[iVertex]     = Buffer_Send_nLinkedNodes[iVertex];
    Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
  }

  for (iVertex = 0; iVertex < nGlobalLinkedNodes; iVertex++)
    Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];
#endif

  /*--- Master process replaced global point indices in Buffer_Receive_LinkedNodes by their indices in
   * Buffer_Receive_GlobalPoint, Buffer_Receive_nLinkedNodes etc. ---*/
  if (rank == MASTER_NODE){
    for (iVertex = 0; iVertex < nGlobalVertex; iVertex++){
      uptr = &Buffer_Receive_LinkedNodes[ Buffer_Receive_StartLinkedNodes[iVertex] ];

      for (unsigned long jLinkedNode = 0; jLinkedNode < Buffer_Receive_nLinkedNodes[iVertex]; jLinkedNode++){
        unsigned long jPoint = uptr[ jLinkedNode ];
        bool found = false; // Global point index has been found
        for (kVertex = 0; kVertex < nGlobalVertex; kVertex++){
          if( Buffer_Receive_GlobalPoint[kVertex] == long(jPoint) ){
            uptr[ jLinkedNode ] = kVertex;
            found = true;
            break;
          }
        }

        if( !found ){ // remove from list
          for (kVertex = jLinkedNode; kVertex < Buffer_Receive_nLinkedNodes[iVertex]-1; kVertex++){
            uptr[ kVertex ] = uptr[ kVertex + 1];
          }
          Buffer_Receive_nLinkedNodes[iVertex]--;
          jLinkedNode--;
        }
      }
    }
  }

  SU2_MPI::Bcast(Buffer_Receive_GlobalPoint, nGlobalVertex, MPI_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_Coord, nGlobalVertex*nDim, MPI_DOUBLE, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_Proc, nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());

  SU2_MPI::Bcast(Buffer_Receive_nLinkedNodes, nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_StartLinkedNodes, nGlobalVertex, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());
  SU2_MPI::Bcast(Buffer_Receive_LinkedNodes, nGlobalLinkedNodes, MPI_UNSIGNED_LONG, 0, SU2_MPI::GetComm());

  delete [] Buffer_Send_Coord;            Buffer_Send_Coord            = nullptr;
  delete [] Buffer_Send_GlobalPoint;      Buffer_Send_GlobalPoint      = nullptr;
  delete [] Buffer_Send_LinkedNodes;      Buffer_Send_LinkedNodes      = nullptr;
  delete [] Buffer_Send_nLinkedNodes;     Buffer_Send_nLinkedNodes     = nullptr;
  delete [] Buffer_Send_StartLinkedNodes; Buffer_Send_StartLinkedNodes = nullptr;

}
