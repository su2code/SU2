/*!
 * \file CMultiGridGeometry.cpp
 * \brief Implementation of the multigrid geometry class.
 * \author F. Palacios, T. Economon
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

#include "../../include/geometry/CMultiGridGeometry.hpp"
#include "../../include/geometry/CMultiGridQueue.hpp"
#include "../../include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CMultiGridGeometry::CMultiGridGeometry(CGeometry* fine_grid, CConfig* config, unsigned short iMesh) : CGeometry() {
  nDim = fine_grid->GetnDim();  // Write the number of dimensions of the coarse grid.

  /*--- Maximum agglomeration size in 2D is 4 nodes, in 3D is 8 nodes. ---*/
  const short int maxAgglomSize = (nDim == 2) ? 4 : 8;

  /*--- Inherit boundary properties from fine grid ---*/
  boundIsStraight = fine_grid->boundIsStraight;

  /*--- Agglomeration Scheme II (Nishikawa, Diskin, Thomas)
        Create a queue system to do the agglomeration
   1st) More than two markers ---> Vertices (never agglomerate)
   2nd) Two markers ---> Edges (agglomerate if same BC, never agglomerate if different BC)
   3rd) One marker ---> Surface (always agglomerate)
   4th) No marker ---> Internal Volume (always agglomerate) ---*/

  // Note that for MPI, we introduce interfaces and we can choose to have agglomeration over
  // the interface or not. Nishikawa chooses not to agglomerate over interfaces.

  /*--- Set a marker to indicate indirect agglomeration, for quads and hexs,
   i.e. consider up to neighbors of neighbors.
   For other levels this information is propagated down during their construction. ---*/

  if (iMesh == MESH_1) {
    for (auto iPoint = 0ul; iPoint < fine_grid->GetnPoint(); iPoint++)
      fine_grid->nodes->SetAgglomerate_Indirect(iPoint, false);

    for (auto iElem = 0ul; iElem < fine_grid->GetnElem(); iElem++) {
      if ((fine_grid->elem[iElem]->GetVTK_Type() == HEXAHEDRON) ||
          (fine_grid->elem[iElem]->GetVTK_Type() == QUADRILATERAL)) {
        for (auto iNode = 0u; iNode < fine_grid->elem[iElem]->GetnNodes(); iNode++) {
          const auto iPoint = fine_grid->elem[iElem]->GetNode(iNode);
          fine_grid->nodes->SetAgglomerate_Indirect(iPoint, true);
        }
      }
    }
  }

  /*--- Create the coarse grid structure using as baseline the fine grid ---*/

  CMultiGridQueue MGQueue_InnerCV(fine_grid->GetnPoint());
  vector<unsigned long> Suitable_Indirect_Neighbors;

  nodes = new CPoint(fine_grid->GetnPoint(), nDim, iMesh, config);

  unsigned long Index_CoarseCV = 0;

  /*--- Statistics for Euler wall agglomeration ---*/
  map<unsigned short, unsigned long> euler_wall_agglomerated, euler_wall_rejected_curvature,
      euler_wall_rejected_straight;
  for (unsigned short iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      euler_wall_agglomerated[iMarker] = 0;
      euler_wall_rejected_curvature[iMarker] = 0;
      euler_wall_rejected_straight[iMarker] = 0;
    }
  }

  /*--- STEP 0: agglomerate the stretched layer above viscous walls along implicit lines, wall CV
   *    included. This runs before the general boundary agglomeration below so that the wall control
   *    volume and the layers stacked on top of it share one footprint; letting the general scheme
   *    claim the wall first would fix a footprint chosen without any knowledge of the lines, and the
   *    stack above it could then only be misaligned with its own base. Everything it claims is
   *    already marked agglomerated, so the boundary and interior passes below simply skip it. ---*/
  if (config->GetMGOptions().MG_Implicit_Lines) {
    AgglomerateImplicitLines(Index_CoarseCV, fine_grid, config, MGQueue_InnerCV);
  }

  /*--- STEP 1: The first step is the boundary agglomeration. ---*/
  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    /*--- Skip periodic boundaries: do not agglomerate on periodic markers. ---*/
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) continue;

    /*--- Skip SEND_RECEIVE markers. Carrying one does not put a point on a boundary, it only
     *    records that the point is mirrored on another rank. A point whose only markers are
     *    SEND_RECEIVE is an interior point, and is left to the domain pass (STEP 2) which is
     *    where a serial run would agglomerate it too. Points that do sit on a physical boundary
     *    are still reached here through their physical marker. ---*/
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;

    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();

      /*--- If the element has not been previously agglomerated and it
      belongs to this physical domain, and it meets the geometrical
      criteria, the agglomeration is studied. ---*/
      vector<short> marker_seed;

      if ((!fine_grid->nodes->GetAgglomerate(iPoint)) && (fine_grid->nodes->GetDomain(iPoint)) &&
          (GeometricalCheck(iPoint, fine_grid, config))) {
        unsigned short nChildren = 1;

        /*--- We set an index for the parent control volume, this
         also marks it as agglomerated. ---*/

        fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);

        /*--- We add the seed point (child) to the parent control volume ---*/

        nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);
        bool agglomerate_seed = false;
        auto counter = 0;
        unsigned short copy_marker[3] = {};
        marker_seed.push_back(iMarker);

        /*--- For a particular point in the fine grid we save all the physical markers that are in
         that point. SEND_RECEIVE markers are deliberately not counted: including them would make
         an ordinary wall point look like a ridge, and a wall/symmetry ridge look like a corner
         (which the counter > 2 rule below then refuses to agglomerate at all), so a point would be
         classified differently depending only on where the partition happens to cut. ---*/

        for (auto jMarker = 0u; jMarker < fine_grid->GetnMarker(); jMarker++) {
          if (config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) continue;
          if (fine_grid->nodes->GetVertex(iPoint, jMarker) != -1) {
            /*--- Count every physical marker, the counter > 2 test needs the true count, but only
             store the first few, which is all the matching rules ever look at. ---*/
            if (counter < 3) copy_marker[counter] = jMarker;
            counter++;

            if (jMarker != iMarker) {
              marker_seed.push_back(jMarker);
            }
          }
        }

        /*--- To agglomerate a vertex it must have only one physical bc.
         This can be improved. If there is only one marker, it is a good
         candidate for agglomeration ---*/

        /*--- 1 BC, so either an edge in 2D or the interior of a plane in 3D ---*/
        /*--- Valley -> Valley : conditionally allowed when both points are on the same marker. ---*/
        /*--- ! Note that in the case of MPI SEND_RECEIVE markers, we might need other conditions ---*/
        if (counter == 1) {
          // The seed/parent is one valley, so we set this part to true
          // if the child is on the same valley, we set it to true as well.
          agglomerate_seed = true;

          /*--- Euler walls: check curvature-based agglomeration criterion ---*/
          if (config->GetMarker_All_KindBC(marker_seed[0]) == EULER_WALL) {
            /*--- Allow agglomeration if marker is straight OR local curvature is small ---*/
            if (!boundIsStraight[marker_seed[0]]) {
              /*--- Compute local curvature at this point ---*/
              su2double local_curvature = ComputeLocalCurvature(fine_grid, iPoint, marker_seed[0]);
              // limit to 45 degrees
              if (local_curvature >= 45.0) {
                agglomerate_seed = false;  // High curvature: do not agglomerate
                euler_wall_rejected_curvature[marker_seed[0]]++;
              } else {
                euler_wall_agglomerated[marker_seed[0]]++;
              }
            } else {
              /*--- Straight wall: agglomerate ---*/
              euler_wall_agglomerated[marker_seed[0]]++;
            }
          }
        }

        /*--- Two physical markers meet here. ---*/
        if (counter == 2) {
          /*--- In 2D that is a genuine corner in the geometry, which is never agglomerated. A wall
                point merely split by a partition interface no longer reaches this branch: it counts
                one physical marker and is handled above as the valley point it is. ---*/
          if (nDim == 2) agglomerate_seed = false;
          /*--- In 3D, this is a ridge point (an edge feature where two surface markers meet).
                Always allow it to attempt agglomeration here; SetBoundAgglomeration() enforces
                the actual ridge-ridge rule downstream: it may only pair with a neighboring ridge
                point that carries the identical physical marker pair. A mismatched marker pair
                usually indicates a genuine sharp corner in the geometry and is correctly left
                un-merged (falls through to the singleton leftover loop). ---*/
          if (nDim == 3) agglomerate_seed = true;

          /*--- Euler walls: check curvature-based agglomeration criterion for both markers ---*/
          // only in 3d because in 2d it's a corner
          bool euler_wall_rejected_here = false;
          for (unsigned short i = 0; i < 2; i++) {
            if ((nDim == 3) && (config->GetMarker_All_KindBC(copy_marker[i]) == EULER_WALL)) {
              if (!boundIsStraight[copy_marker[i]]) {
                /*--- Compute local curvature at this point ---*/
                su2double local_curvature = ComputeLocalCurvature(fine_grid, iPoint, copy_marker[i]);
                // limit to 45 degrees
                if (local_curvature >= 45.0) {
                  agglomerate_seed = false;  // High curvature: do not agglomerate
                  euler_wall_rejected_curvature[copy_marker[i]]++;
                  euler_wall_rejected_here = true;
                }
              }
              /*--- Track agglomeration if not rejected ---*/
              if (agglomerate_seed && !euler_wall_rejected_here) {
                euler_wall_agglomerated[copy_marker[i]]++;
              }
            }
          }
        }

        /*--- If there are more than 2 markers, the aglomeration will be discarded ---*/

        if (counter > 2) agglomerate_seed = false;

        /*--- If the seed (parent) can be agglomerated, we try to agglomerate connected childs to the parent ---*/
        /*--- Note that in 2D we allow a maximum of 4 nodes to be agglomerated ---*/

        if (agglomerate_seed) {
          /*--- Now we do a sweep over all the nodes that surround the seed point ---*/

          for (auto CVPoint : fine_grid->nodes->GetPoints(iPoint)) {
            /*--- The new point can be agglomerated ---*/

            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config)) {
              /*--- We set the value of the parent ---*/

              fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

              /*--- We set the value of the child ---*/

              nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
              nChildren++;
              /*--- In 2D, we agglomerate exactly 2 nodes if the nodes are on the line edge. ---*/
              if ((nDim == 2) && (counter == 1)) break;
              /*--- In 3D, we agglomerate exactly 2 nodes if the nodes are on the surface edge. ---*/
              if ((nDim == 3) && (counter == 2)) break;
              /*--- Apply maxAgglomSize limit for 3D internal boundary face nodes (counter==1 in 3D). ---*/
              if (nChildren == maxAgglomSize) break;
            }
          }

          /*--- Only take into account indirect neighbors for 3D faces, not 2D. ---*/
          if (nDim == 3) {
            Suitable_Indirect_Neighbors.clear();

            if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint))
              SetSuitableNeighbors(Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);

            /*--- Now we do a sweep over all the indirect nodes that can be added ---*/

            for (auto CVPoint : Suitable_Indirect_Neighbors) {
              /*--- The new point can be agglomerated ---*/

              if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config)) {
                /*--- We set the value of the parent ---*/

                fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

                /*--- We set the indirect agglomeration information of the corse point
                based on its children in the fine grid. ---*/

                if (fine_grid->nodes->GetAgglomerate_Indirect(CVPoint))
                  nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);

                /*--- We set the value of the child ---*/

                nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
                nChildren++;
                /*--- Apply maxAgglomSize limit for 3D internal boundary face nodes. ---*/
                if (nChildren == maxAgglomSize) break;
              }
            }
          }
        }

        /*--- Update the number of children of the coarse control volume. ---*/

        nodes->SetnChildren_CV(Index_CoarseCV, nChildren);
        Index_CoarseCV++;
      }
    }
  }

  /*--- Do not agglomerate any leftover node with more than one physical boundary condition,
   i.e. make one coarse CV with a single child. ---*/

  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    /*--- As in STEP 1, a SEND_RECEIVE marker does not make a point a boundary point. Turning the
     leftovers of those markers into single-child coarse CVs here would strand every interior point
     along a partition interface before the domain pass below ever gets to see it. ---*/
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;

    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();

      if ((!fine_grid->nodes->GetAgglomerate(iPoint)) && (fine_grid->nodes->GetDomain(iPoint))) {
        fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);
        nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);
        nodes->SetnChildren_CV(Index_CoarseCV, 1);
        Index_CoarseCV++;
      }
    }
  }

  /*--- Update the queue with the results from the boundary agglomeration ---*/

  for (auto iPoint = 0ul; iPoint < fine_grid->GetnPoint(); iPoint++) {
    if (fine_grid->nodes->GetAgglomerate(iPoint)) {
      MGQueue_InnerCV.RemoveCV(iPoint);

    } else {
      /*--- Count the number of agglomerated neighbors, and modify the queue,
       Points with more agglomerated neighbors are processed first. ---*/

      short priority = 0;
      for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) {
        priority += fine_grid->nodes->GetAgglomerate(jPoint);
      }
      MGQueue_InnerCV.MoveCV(iPoint, priority);
    }
  }

  /*--- STEP 2: Agglomerate the domain points. ---*/

  auto iteration = 0ul;
  while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {
    const auto iPoint = MGQueue_InnerCV.NextCV();
    iteration++;

    /*--- If the element has not been previously agglomerated, belongs to the physical domain,
     and satisfies several geometrical criteria then the seed CV is accepted for agglomeration. ---*/

    if ((!fine_grid->nodes->GetAgglomerate(iPoint)) && (fine_grid->nodes->GetDomain(iPoint)) &&
        (GeometricalCheck(iPoint, fine_grid, config))) {
      unsigned short nChildren = 1;

      /*--- We set an index for the parent control volume ---*/

      fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);

      /*--- We add the seed point (child) to the parent control volume ---*/

      nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);

      /*--- Update the queue with the seed point (remove the seed and
       increase the priority of its neighbors) ---*/

      MGQueue_InnerCV.Update(iPoint, fine_grid);

      /*--- Now we do a sweep over all the nodes that surround the seed point ---*/

      for (auto CVPoint : fine_grid->nodes->GetPoints(iPoint)) {
        /*--- Determine if the CVPoint can be agglomerated ---*/

        if ((!fine_grid->nodes->GetAgglomerate(CVPoint)) && (fine_grid->nodes->GetDomain(CVPoint)) &&
            (GeometricalCheck(CVPoint, fine_grid, config))) {
          /*--- We set the value of the parent ---*/

          fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

          /*--- We set the value of the child ---*/

          nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
          nChildren++;

          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of its neighbors) ---*/

          MGQueue_InnerCV.Update(CVPoint, fine_grid);
        }
        if (nChildren == maxAgglomSize) break;
      }

      /*--- Identify the indirect neighbors ---*/

      Suitable_Indirect_Neighbors.clear();
      if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint))
        SetSuitableNeighbors(Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);

      /*--- Now we do a sweep over all the indirect nodes that can be added ---*/

      for (auto CVPoint : Suitable_Indirect_Neighbors) {
        // if we have reached the maximum, get out.
        if (nChildren == maxAgglomSize) break;
        /*--- The new point can be agglomerated ---*/

        if ((!fine_grid->nodes->GetAgglomerate(CVPoint)) && (fine_grid->nodes->GetDomain(CVPoint))) {
          /*--- We set the value of the parent ---*/

          fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

          /*--- We set the indirect agglomeration information ---*/

          if (fine_grid->nodes->GetAgglomerate_Indirect(CVPoint)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);

          /*--- We set the value of the child ---*/

          nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
          nChildren++;

          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/

          MGQueue_InnerCV.Update(CVPoint, fine_grid);
        }
      }

      /*--- Update the number of control of childrens ---*/

      nodes->SetnChildren_CV(Index_CoarseCV, nChildren);
      Index_CoarseCV++;
    } else {
      /*--- The seed point can not be agglomerated because of size, domain, streching, etc.
       move the point to the lowest priority ---*/

      MGQueue_InnerCV.MoveCV(iPoint, -1);
    }
  }

  /*--- Convert any point that was not agglomerated into a coarse point. ---*/

  for (auto iPoint = 0ul; iPoint < fine_grid->GetnPoint(); iPoint++) {
    if ((!fine_grid->nodes->GetAgglomerate(iPoint)) && (fine_grid->nodes->GetDomain(iPoint))) {
      fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);
      if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);
      nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);
      nodes->SetnChildren_CV(Index_CoarseCV, 1);
      Index_CoarseCV++;
    }
  }

  nPointDomain = Index_CoarseCV;
  nPoint = nPointDomain;

  /*--- Check that there are no hanging nodes. Detect isolated points
   (only 1 neighbor), and merge their children CV's with the neighbor. ---*/

  SetPoint_Connectivity(fine_grid);

  /*--- The connectivity just built only knows about coarse CVs of this rank: the halo CVs do not
   exist until the MPI relay below runs, and the relay cannot run earlier because it broadcasts the
   parent indices that the merge here is still free to change. So a CV touching a partition boundary
   may look isolated while actually having neighbors on the other rank, and merging it would be
   wrong. Mark those CVs and leave them alone; a genuinely isolated CV in the interior is unaffected.
   Note this deliberately keeps the conservative outcome the sentinel used to produce by accident,
   but only for the CVs that really do border another rank rather than for every CV near one. ---*/

  vector<bool> touchesPartition(nPointDomain, false);
  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++) {
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
      const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      for (auto iFinePoint_Neighbor : fine_grid->nodes->GetPoints(iFinePoint)) {
        if (fine_grid->nodes->GetParent_CV(iFinePoint_Neighbor) == std::numeric_limits<unsigned long>::max()) {
          touchesPartition[iCoarsePoint] = true;
          break;
        }
      }
      if (touchesPartition[iCoarsePoint]) break;
    }
  }

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++) {
    if ((nodes->GetnPoint(iCoarsePoint) == 1) && !touchesPartition[iCoarsePoint]) {
      /*--- Find the neighbor of the isolated point. This neighbor is the right control volume ---*/

      const auto iCoarsePoint_Complete = nodes->GetPoint(iCoarsePoint, 0);

      /*--- Check if merging would exceed the maximum agglomeration size ---*/
      auto nChildren_Target = nodes->GetnChildren_CV(iCoarsePoint_Complete);
      auto nChildren_Isolated = nodes->GetnChildren_CV(iCoarsePoint);
      auto nChildren_Total = nChildren_Target + nChildren_Isolated;

      /*--- If the total would exceed maxAgglomSize, try to redistribute children to neighbors ---*/
      if (nChildren_Total > maxAgglomSize) {
        /*--- Find neighbors of the target coarse point that have room ---*/
        unsigned short nChildrenToRedistribute = nChildren_Total - maxAgglomSize;

        for (auto jCoarsePoint : nodes->GetPoints(iCoarsePoint_Complete)) {
          if (nChildrenToRedistribute == 0) break;

          auto nChildren_Neighbor = nodes->GetnChildren_CV(jCoarsePoint);
          if (nChildren_Neighbor < maxAgglomSize) {
            unsigned short nCanTransfer =
                min(nChildrenToRedistribute, static_cast<unsigned short>(maxAgglomSize - nChildren_Neighbor));

            /*--- Transfer children from target to neighbor ---*/
            for (unsigned short iTransfer = 0; iTransfer < nCanTransfer; iTransfer++) {
              /*--- Take from the end of the target's children list ---*/
              auto nChildren_Current = nodes->GetnChildren_CV(iCoarsePoint_Complete);
              if (nChildren_Current > 0) {
                auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint_Complete, nChildren_Current - 1);

                /*--- Add to neighbor ---*/
                auto nChildren_Neighbor_Current = nodes->GetnChildren_CV(jCoarsePoint);
                nodes->SetChildren_CV(jCoarsePoint, nChildren_Neighbor_Current, iFinePoint);
                nodes->SetnChildren_CV(jCoarsePoint, nChildren_Neighbor_Current + 1);

                /*--- Update parent ---*/
                fine_grid->nodes->SetParent_CV(iFinePoint, jCoarsePoint);

                /*--- Remove from target (by reducing count) ---*/
                nodes->SetnChildren_CV(iCoarsePoint_Complete, nChildren_Current - 1);

                nChildrenToRedistribute--;
              }
            }
          }
        }

        /*--- Update the target's child count after redistribution ---*/
        nChildren_Target = nodes->GetnChildren_CV(iCoarsePoint_Complete);
      }

      /*--- Add the isolated point's children to the target control volume ---*/
      auto nChildren = nChildren_Target;
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
        const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        nodes->SetChildren_CV(iCoarsePoint_Complete, nChildren, iFinePoint);
        nChildren++;
        fine_grid->nodes->SetParent_CV(iFinePoint, iCoarsePoint_Complete);
      }

      /*--- Update the number of children control volumes ---*/

      nodes->SetnChildren_CV(iCoarsePoint_Complete, nChildren);
      nodes->SetnChildren_CV(iCoarsePoint, 0);
    }
  }

  /*--- Reset the neighbor information. ---*/

  nodes->ResetPoints();

#ifdef HAVE_MPI
  /*--- Reset halo point parents before MPI agglomeration.
   When creating level N from level N-1, the fine grid (level N-1)
   already has Parent_CV set from when it was created from level N-2.
   Those parent indices point to level N, but when creating level N+1, they would be
   incorrectly interpreted as level N+1 indices. ---*/

  for (auto iPoint = fine_grid->GetnPointDomain(); iPoint < fine_grid->GetnPoint(); iPoint++) {
    fine_grid->nodes->SetParent_CV(iPoint, std::numeric_limits<unsigned long>::max());
  }

  /*--- Dealing with MPI parallelization, the objective is that the received nodes must be agglomerated
   in the same way as the donor (send) nodes. Send the node agglomeration information of the donor
   (parent and children). The agglomerated halos of this rank are set according to the rank where
   they are domain points. ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      const auto MarkerS = iMarker;      // sending marker
      const auto MarkerR = iMarker + 1;  // receiving marker

      const auto send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
      const auto receive_from = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;

      const auto nVertexS = fine_grid->nVertex[MarkerS];
      const auto nVertexR = fine_grid->nVertex[MarkerR];

      /*--- Allocate Receive and Send buffers  ---*/

      vector<unsigned long> Buffer_Receive_Children(nVertexR);
      vector<unsigned long> Buffer_Send_Children(nVertexS);

      vector<unsigned long> Buffer_Receive_Parent(nVertexR);
      vector<unsigned long> Buffer_Send_Parent(nVertexS);

      /*--- Copy the information that should be sent, child and parent indices. ---*/

      for (auto iVertex = 0ul; iVertex < nVertexS; iVertex++) {
        const auto iPoint = fine_grid->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Children[iVertex] = iPoint;
        Buffer_Send_Parent[iVertex] = fine_grid->nodes->GetParent_CV(iPoint);
      }

      /*--- Send/Receive information. ---*/

      SU2_MPI::Sendrecv(Buffer_Send_Children.data(), nVertexS, MPI_UNSIGNED_LONG, send_to, 0,
                        Buffer_Receive_Children.data(), nVertexR, MPI_UNSIGNED_LONG, receive_from, 0,
                        SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Sendrecv(Buffer_Send_Parent.data(), nVertexS, MPI_UNSIGNED_LONG, send_to, 1,
                        Buffer_Receive_Parent.data(), nVertexR, MPI_UNSIGNED_LONG, receive_from, 1, SU2_MPI::GetComm(),
                        MPI_STATUS_IGNORE);

      /*--- Create a list of the parent nodes without duplicates. ---*/

      auto Aux_Parent = Buffer_Receive_Parent;

      sort(Aux_Parent.begin(), Aux_Parent.end());
      auto it1 = unique(Aux_Parent.begin(), Aux_Parent.end());
      Aux_Parent.resize(it1 - Aux_Parent.begin());

      /*--- Create the local and remote vector for the parents and children CVs. ---*/

      const auto& Parent_Remote = Buffer_Receive_Parent;
      vector<unsigned long> Parent_Local(nVertexR);
      vector<unsigned long> Children_Local(nVertexR);

      /*--- First pass: Determine which parents will actually be used (have non-skipped children).
       This prevents creating orphaned halo CVs that have coordinates (0,0,0). ---*/
      vector<bool> parent_used(Aux_Parent.size(), false);
      vector<unsigned long> parent_local_index(Aux_Parent.size(), std::numeric_limits<unsigned long>::max());

      for (auto iVertex = 0ul; iVertex < nVertexR; iVertex++) {
        const auto iPoint_Fine = fine_grid->vertex[MarkerR][iVertex]->GetNode();
        auto existing_parent = fine_grid->nodes->GetParent_CV(iPoint_Fine);

        /*--- Skip if already agglomerated (first-wins policy) ---*/
        if (existing_parent != std::numeric_limits<unsigned long>::max()) continue;

        /*--- Skip if received parent is invalid (sending rank didn't agglomerate this point) ---*/
        if (Parent_Remote[iVertex] == std::numeric_limits<unsigned long>::max()) continue;

        /*--- Find which parent this vertex maps to ---*/
        for (auto jVertex = 0ul; jVertex < Aux_Parent.size(); jVertex++) {
          if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
            parent_used[jVertex] = true;
            break;
          }
        }
      }

      /*--- Assign local indices only to used parents ---*/
      unsigned long nUsedParents = 0;
      for (auto jVertex = 0ul; jVertex < Aux_Parent.size(); jVertex++) {
        if (parent_used[jVertex]) {
          parent_local_index[jVertex] = Index_CoarseCV + nUsedParents;
          nUsedParents++;
        }
      }

      /*--- Now map each received vertex to its local parent ---*/
      for (auto iVertex = 0ul; iVertex < nVertexR; iVertex++) {
        Parent_Local[iVertex] = std::numeric_limits<unsigned long>::max();
        for (auto jVertex = 0ul; jVertex < Aux_Parent.size(); jVertex++) {
          if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
            Parent_Local[iVertex] = parent_local_index[jVertex];
            break;
          }
        }

        Children_Local[iVertex] = fine_grid->vertex[MarkerR][iVertex]->GetNode();
      }

      /*--- Only increment by the number of parents that will actually be used ---*/
      Index_CoarseCV += nUsedParents;

      /*--- Create the final structure ---*/
      for (auto iVertex = 0ul; iVertex < nVertexR; iVertex++) {
        const auto iPoint_Fine = Children_Local[iVertex];

        /*--- Skip if this halo point was already agglomerated (first-wins policy) ---*/
        auto existing_parent = fine_grid->nodes->GetParent_CV(iPoint_Fine);
        if (existing_parent != std::numeric_limits<unsigned long>::max()) continue;

        /*--- Skip if parent mapping is invalid (sender didn't agglomerate) ---*/
        const auto iPoint_Coarse = Parent_Local[iVertex];
        if (iPoint_Coarse == std::numeric_limits<unsigned long>::max()) continue;

        /*--- Append to existing children, don't overwrite ---*/
        auto existing_children_count = nodes->GetnChildren_CV(iPoint_Coarse);

        fine_grid->nodes->SetParent_CV(iPoint_Fine, iPoint_Coarse);
        nodes->SetChildren_CV(iPoint_Coarse, existing_children_count, iPoint_Fine);
        nodes->SetnChildren_CV(iPoint_Coarse, existing_children_count + 1);
        nodes->SetDomain(iPoint_Coarse, false);
      }
    }
  }

#endif  // HAVE_MPI

  /*--- Update the number of points after the MPI agglomeration ---*/

  nPoint = Index_CoarseCV;

  /*--- Console output with the summary of the agglomeration ---*/
  unsigned long nPointFine = fine_grid->GetnPointDomain();
  unsigned long Global_nPointCoarse, Global_nPointFine, Min_nPointCoarse;

  SU2_MPI::Allreduce(&nPointDomain, &Global_nPointCoarse, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nPointFine, &Global_nPointFine, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nPointDomain, &Min_nPointCoarse, 1, MPI_UNSIGNED_LONG, MPI_MIN, SU2_MPI::GetComm());

  SetGlobal_nPointDomain(Global_nPointCoarse);

  if (iMesh != MESH_0) {
    /*--- Initialize coarse-level CFL from config. MG_CFL_SCALING will
          apply per-level reductions during the multigrid cycle. ---*/
    config->SetCFL(iMesh, config->GetCFL(MESH_0));
  }

  const su2double ratio = su2double(Global_nPointFine) / su2double(Global_nPointCoarse);

  /*--- Stop coarsening once the smallest per-rank partition falls below the minimum,
        not just the summed total, since each rank runs its own MG hierarchy locally
        and a partition that agglomerates down to too few CVs degenerates the operator
        on that rank even if other ranks still have plenty of points. ---*/
  if (Min_nPointCoarse < config->GetMGOptions().MG_Min_MeshSize) {
    if (rank == MASTER_NODE)
      cout << "MG level " << iMesh << " has only " << Min_nPointCoarse
           << " CVs on the smallest partition (< MG_MIN_MESHSIZE=" << config->GetMGOptions().MG_Min_MeshSize
           << "). Reducing MG levels to " << iMesh - 1 << "." << endl;
    config->SetMGLevels(iMesh - 1);
  } else if (rank == MASTER_NODE) {
    PrintingToolbox::CTablePrinter MGTable(&std::cout);
    MGTable.AddColumn("MG Level", 10);
    MGTable.AddColumn("CVs", 10);
    MGTable.AddColumn("Aggl. Rate", 10);
    MGTable.AddColumn("CFL", 10);
    MGTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);

    if (iMesh == MESH_1) {
      MGTable.PrintHeader();
      MGTable << iMesh - 1 << Global_nPointFine << "1/1.00" << config->GetCFL(iMesh - 1);
    }
    stringstream ss;
    ss << "1/" << std::setprecision(3) << ratio;
    MGTable << iMesh << Global_nPointCoarse << ss.str() << config->GetCFL(iMesh);
    if (iMesh == config->GetnMGLevels()) {
      MGTable.PrintFooter();
    }
  }

  /*--- Output Euler wall agglomeration statistics ---*/
  if (rank == MASTER_NODE) {
    /*--- Gather global statistics for Euler walls ---*/
    bool has_euler_walls = false;
    for (unsigned short iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
        has_euler_walls = true;
        break;
      }
    }

    if (has_euler_walls) {
      cout << endl;
      cout << "Euler Wall Agglomeration Statistics (45° curvature threshold):" << endl;
      cout << "----------------------------------------------------------------" << endl;

      for (unsigned short iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
          string marker_name = config->GetMarker_All_TagBound(iMarker);
          unsigned long agglomerated = euler_wall_agglomerated[iMarker];
          unsigned long rejected = euler_wall_rejected_curvature[iMarker];
          unsigned long total = agglomerated + rejected;

          if (total > 0) {
            su2double accept_rate = 100.0 * su2double(agglomerated) / su2double(total);
            cout << "  Marker: " << marker_name << endl;
            cout << "    Seeds agglomerated:       " << agglomerated << " (" << std::setprecision(1) << std::fixed
                 << accept_rate << "%)" << endl;
            cout << "    Seeds rejected (>45° curv): " << rejected << " (" << std::setprecision(1) << std::fixed
                 << (100.0 - accept_rate) << "%)" << endl;
          }
        }
      }
      cout << "----------------------------------------------------------------" << endl;
    }
  }

  edgeColorGroupSize = config->GetEdgeColoringGroupSize();
}

bool CMultiGridGeometry::GeometricalCheck(unsigned long iPoint, const CGeometry* fine_grid,
                                          const CConfig* config) const {
  su2double max_dimension = 1.2;

  /*--- Evaluate the total size of the element ---*/

  bool Volume = true;
  su2double ratio = pow(fine_grid->nodes->GetVolume(iPoint), 1.0 / su2double(nDim)) * max_dimension;
  su2double limit = pow(config->GetDomainVolume(), 1.0 / su2double(nDim));
  if (ratio > limit) {
    Volume = false;
    cout << "Volume limit reached!" << endl;
  }

  return (Volume);
}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, vector<short> marker_seed,
                                               const CGeometry* fine_grid, const CConfig* config) const {
  bool agglomerate_CV = false;

  /*--- Basic condition, the point has not been previously agglomerated, it belongs to the domain,
   and has passed some basic geometrical checks. ---*/

  if ((!fine_grid->nodes->GetAgglomerate(CVPoint)) && (fine_grid->nodes->GetDomain(CVPoint)) &&
      (GeometricalCheck(CVPoint, fine_grid, config))) {
    /*--- If the point belongs to a boundary, its type must be compatible with the seed marker. ---*/

    int counter = 0;
    unsigned short copy_marker[3] = {};

    if (fine_grid->nodes->GetBoundary(CVPoint)) {
      /*--- Identify the physical markers of the vertex that we want to agglomerate. SEND_RECEIVE
       markers are skipped for the same reason as on the seed side: they say nothing about the
       boundary condition the candidate carries, only that it is mirrored on another rank. A
       candidate whose markers are all SEND_RECEIVE therefore ends up with counter == 0 and is
       rejected below, which is the answer a serial run gives for the interior point it really is. ---*/

      for (auto jMarker = 0u; jMarker < fine_grid->GetnMarker(); jMarker++) {
        if (config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) continue;
        if (fine_grid->nodes->GetVertex(CVPoint, jMarker) != -1) {
          if (counter < 3) copy_marker[counter] = jMarker;
          counter++;
        }
      }

      /*--- The basic condition is that the agglomerated vertex must have the same physical marker
       as the seed. ---*/

      /*--- Only one marker in the vertex that is going to be agglomerated ---*/

      /*--- Valley -> Valley: only if of the same type---*/
      if (counter == 1) {
        /*--- We agglomerate if there is only one marker and it is the same marker as the seed marker ---*/
        // So this is the case when in 2D we are on an edge, and in 3D we are in the interior of a surface.
        // note that this should be the same marker id, not just the same marker type.
        /*--- Because neither side counts SEND_RECEIVE any more, this test is symmetric: a wall point
         on a partition interface is now just as absorbable by an ordinary wall seed as the other way
         around. Previously only an interface seed could absorb an interface candidate, so whenever an
         ordinary seed reached the shared neighbours first the interface point was stranded. ---*/
        if ((marker_seed.size() == 1) && (copy_marker[0] == marker_seed[0])) agglomerate_CV = true;
      }

      /*--- If there are two markers in the vertex that is going to be aglomerated ---*/
      if (counter == 2) {
        /*--- In 2D this is a corner and we do not agglomerate ---*/
        if (nDim == 2) {
          agglomerate_CV = false;
        }
        /*--- Both markers have to be the same. ---*/
        else if (marker_seed.size() == 2) {
          if (((copy_marker[0] == marker_seed[0]) && (copy_marker[1] == marker_seed[1])) ||
              ((copy_marker[0] == marker_seed[1]) && (copy_marker[1] == marker_seed[0]))) {
            agglomerate_CV = true;
          }
        }
      }
    }
    /*--- If the element belongs to the domain, it is never agglomerated with a boundary node. ---*/
    else {
      agglomerate_CV = false;
    }
  }

  return agglomerate_CV;
}

/*--- ---*/

void CMultiGridGeometry::SetSuitableNeighbors(vector<unsigned long>& Suitable_Indirect_Neighbors, unsigned long iPoint,
                                              unsigned long Index_CoarseCV, const CGeometry* fine_grid) const {
  /*--- Create a list with the first neighbors, including the seed. ---*/

  vector<unsigned long> First_Neighbor_Points;
  First_Neighbor_Points.push_back(iPoint);
  for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) First_Neighbor_Points.push_back(jPoint);

  /*--- Create a list with the second neighbors, without first, and seed neighbors. ---*/

  vector<unsigned long> Second_Neighbor_Points, Second_Origin_Points, Suitable_Second_Neighbors;

  for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) {
    for (auto kPoint : fine_grid->nodes->GetPoints(jPoint)) {
      /*--- Check that the second neighbor does not belong to the first neighbors or the seed. ---*/

      auto end = First_Neighbor_Points.end();
      if (find(First_Neighbor_Points.begin(), end, kPoint) == end) {
        Second_Neighbor_Points.push_back(kPoint);  // neighbor of a neighbor, not connected to original ipoint
        Second_Origin_Points.push_back(jPoint);    // the neighbor that is connected to ipoint
      }
    }
  }

  /*--- Identify those second neighbors that are repeated (candidates to be added).
   For a mesh of quads this produces a 9-point stencil from the 5-point of direct
   neighbors, and for hexs it produces a 27-point stencil. ---*/

  for (auto iNeighbor = 0ul; iNeighbor < Second_Neighbor_Points.size(); iNeighbor++) {
    for (auto jNeighbor = iNeighbor + 1; jNeighbor < Second_Neighbor_Points.size(); jNeighbor++) {
      /*--- Repeated second neighbor with different origin ---*/

      if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) &&
          (Second_Origin_Points[iNeighbor] != Second_Origin_Points[jNeighbor])) {
        Suitable_Indirect_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);

        /*--- Create a list of suitable second neighbors, that we will use
         to compute the third neighbors. --*/

        Suitable_Second_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);
      }
    }
  }

  /*--- Remove duplicates ---*/

  sort(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  auto it1 = unique(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  Suitable_Second_Neighbors.resize(it1 - Suitable_Second_Neighbors.begin());
}

void CMultiGridGeometry::SetPoint_Connectivity(const CGeometry* fine_grid) {
  /*--- Temporary, CPoint (nodes) then compresses this structure. ---*/
  vector<vector<unsigned long>> points(nPoint);

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPoint; iCoarsePoint++) {
    /*--- For each child CV (of the fine grid), ---*/
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
      const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      /*--- loop over the parent CVs (coarse grid) of its (fine) neighbors. ---*/
      for (auto iFinePoint_Neighbor : fine_grid->nodes->GetPoints(iFinePoint)) {
        const auto iParent = fine_grid->nodes->GetParent_CV(iFinePoint_Neighbor);
        /*--- Skip neighbors whose parent is not known yet. The first call to this function happens
         during construction, before the MPI relay has assigned parents to the fine grid's halo
         points, so those still hold the sentinel. Letting it through would add one fake neighbor to
         every coarse CV along a partition boundary, which both corrupts nNeighbor and hides the CV
         from the isolated-CV repair. The driver calls this again once the relay has run. ---*/
        if (iParent == std::numeric_limits<unsigned long>::max()) continue;
        /*--- If it is not the target coarse point, it is a coarse neighbor. ---*/
        if (iParent != iCoarsePoint) {
          /*--- Avoid duplicates. ---*/
          auto End = points[iCoarsePoint].end();
          if (find(points[iCoarsePoint].begin(), End, iParent) == End) points[iCoarsePoint].push_back(iParent);
        }
      }
    }

    /*--- See CPhysicalGeometry::SetPoint_Connectivity for why we sort. ---*/
    sort(points[iCoarsePoint].begin(), points[iCoarsePoint].end());

    /*--- Set the number of neighbors variable, this is
     important for JST and multigrid in parallel ---*/
    nodes->SetnNeighbor(iCoarsePoint, points[iCoarsePoint].size());
  }

  nodes->SetPoints(points);
}

void CMultiGridGeometry::SetVertex(const CGeometry* fine_grid, const CConfig* config) {
  nMarker = fine_grid->GetnMarker();
  unsigned short nMarker_Max = config->GetnMarker_Max();

  /*--- If any children node belong to the boundary then the entire control
   volume will belong to the boundary ---*/
  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPoint; iCoarsePoint++)
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
      auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      if (fine_grid->nodes->GetBoundary(iFinePoint)) {
        nodes->SetBoundary(iCoarsePoint, nMarker);
        break;
      }
    }

  vertex = new CVertex**[nMarker];
  nVertex = new unsigned long[nMarker];

  Tag_to_Marker = new string[nMarker_Max];
  for (auto iMarker_Tag = 0u; iMarker_Tag < nMarker_Max; iMarker_Tag++)
    Tag_to_Marker[iMarker_Tag] = fine_grid->GetMarker_Tag(iMarker_Tag);

  /*--- Compute the number of vertices to do the dimensionalization ---*/
  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPoint; iCoarsePoint++) {
    if (nodes->GetBoundary(iCoarsePoint)) {
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
        auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
          if ((fine_grid->nodes->GetVertex(iFinePoint, iMarker) != -1) &&
              (nodes->GetVertex(iCoarsePoint, iMarker) == -1)) {
            auto iVertex = nVertex[iMarker];
            nodes->SetVertex(iCoarsePoint, iVertex, iMarker);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex*[fine_grid->GetnVertex(iMarker) + 1];
    nVertex[iMarker] = 0;
  }

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPoint; iCoarsePoint++)
    if (nodes->GetBoundary(iCoarsePoint))
      for (auto iMarker = 0u; iMarker < nMarker; iMarker++) nodes->SetVertex(iCoarsePoint, -1, iMarker);

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPoint; iCoarsePoint++) {
    if (nodes->GetBoundary(iCoarsePoint)) {
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
        auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
          if ((fine_grid->nodes->GetVertex(iFinePoint, iMarker) != -1) &&
              (nodes->GetVertex(iCoarsePoint, iMarker) == -1)) {
            auto iVertex = nVertex[iMarker];
            vertex[iMarker][iVertex] = new CVertex(iCoarsePoint, nDim);
            nodes->SetVertex(iCoarsePoint, iVertex, iMarker);

            /*--- Set the transformation to apply ---*/
            unsigned long ChildVertex = fine_grid->nodes->GetVertex(iFinePoint, iMarker);
            unsigned short RotationKind = fine_grid->vertex[iMarker][ChildVertex]->GetRotation_Type();
            vertex[iMarker][iVertex]->SetRotation_Type(RotationKind);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }
}

void CMultiGridGeometry::MatchActuator_Disk(const CConfig* config) {
  int iProcessor = size;

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
        (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
      for (auto iVertex = 0u; iVertex < nVertex[iMarker]; iVertex++) {
        auto iPoint = vertex[iMarker][iVertex]->GetNode();
        if (nodes->GetDomain(iPoint)) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, nodes->GetGlobalIndex(iPoint), iVertex, iMarker, iProcessor);
        }
      }
    }
  }
}

void CMultiGridGeometry::SetControlVolume(const CGeometry* fine_grid, unsigned short action) {
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Compute the area of the coarse volume ---*/
    for (auto iCoarsePoint = 0u; iCoarsePoint < nPoint; iCoarsePoint++) {
      nodes->SetVolume(iCoarsePoint, 0.0);
      su2double Coarse_Volume = 0.0;
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
        auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        Coarse_Volume += fine_grid->nodes->GetVolume(iFinePoint);
      }
      nodes->SetVolume(iCoarsePoint, max(Coarse_Volume, EPS));
    }

    /*--- Update or not the values of faces at the edge ---*/
    if (action != ALLOCATE) {
      edges->SetZeroValues();
    }

    for (auto iCoarsePoint = 0u; iCoarsePoint < nPoint; iCoarsePoint++)
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
        auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);

        for (auto iFinePoint_Neighbor : fine_grid->nodes->GetPoints(iFinePoint)) {
          auto iParent = fine_grid->nodes->GetParent_CV(iFinePoint_Neighbor);
          if ((iParent != iCoarsePoint) && (iParent < iCoarsePoint)) {
            auto FineEdge = fine_grid->FindEdge(iFinePoint, iFinePoint_Neighbor);

            bool change_face_orientation = false;
            if (iFinePoint < iFinePoint_Neighbor) change_face_orientation = true;

            auto CoarseEdge = FindEdge(iParent, iCoarsePoint);

            const auto Normal = fine_grid->edges->GetNormal(FineEdge);

            if (change_face_orientation) {
              edges->SubNormal(CoarseEdge, Normal);
            } else {
              edges->AddNormal(CoarseEdge, Normal);
            }
          }
        }
      }

    /*--- Check if there is a normal with null area ---*/

    for (auto iEdge = 0u; iEdge < nEdge; iEdge++) {
      const auto NormalFace = edges->GetNormal(iEdge);
      su2double Area = GeometryToolbox::Norm(nDim, NormalFace);
      if (Area == 0.0) {
        su2double DefaultNormal[3] = {EPS * EPS};
        edges->SetNormal(iEdge, DefaultNormal);
      }
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CMultiGridGeometry::SetBoundControlVolume(const CGeometry* fine_grid, const CConfig* config,
                                               unsigned short action) {
  su2double Normal[MAXNDIM] = {0.0}, *NormalFace = nullptr;

  if (action != ALLOCATE) {
    SU2_OMP_FOR_DYN(1)
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++)
      for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) vertex[iMarker][iVertex]->SetZeroValues();
    END_SU2_OMP_FOR
  }

  SU2_OMP_FOR_DYN(1)
  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
      auto iCoarsePoint = vertex[iMarker][iVertex]->GetNode();
      for (auto iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
        auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        if (fine_grid->nodes->GetVertex(iFinePoint, iMarker) != -1) {
          auto FineVertex = fine_grid->nodes->GetVertex(iFinePoint, iMarker);
          fine_grid->vertex[iMarker][FineVertex]->GetNormal(Normal);
          vertex[iMarker][iVertex]->AddNormal(Normal);
        }
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- Check if there is a normal with null area ---*/
  SU2_OMP_FOR_DYN(1)
  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      su2double Area = GeometryToolbox::Norm(nDim, NormalFace);
      if (Area == 0.0)
        for (auto iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS * EPS;
    }
  }
  END_SU2_OMP_FOR

  SU2_OMP_SAFE_GLOBAL_ACCESS(ComputeModifiedSymmetryNormals(config);)
}

void CMultiGridGeometry::SetCoord(const CGeometry* fine_grid) {
  SU2_OMP_FOR_STAT(roundUpDiv(nPoint, omp_get_max_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < nPoint; Point_Coarse++) {
    auto Area_Parent = nodes->GetVolume(Point_Coarse);
    su2double Coordinates[3] = {0.0};
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
      auto Area_Children = fine_grid->nodes->GetVolume(Point_Fine);
      auto Coordinates_Fine = fine_grid->nodes->GetCoord(Point_Fine);
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Coordinates[iDim] += Coordinates_Fine[iDim] * Area_Children / Area_Parent;
    }
    nodes->SetCoord(Point_Coarse, Coordinates);
  }
  END_SU2_OMP_FOR
}

void CMultiGridGeometry::SetMultiGridWallHeatFlux(const CGeometry* fine_grid, unsigned short val_marker) {
  struct {
    const CGeometry* fine_grid;
    unsigned short marker;
    su2double* target;

    su2double Get(unsigned long iVertex) const { return fine_grid->GetCustomBoundaryHeatFlux(marker, iVertex); }
    void Set(unsigned long iVertex, const su2double& val) const { target[iVertex] = val; }

  } wall_heat_flux;

  wall_heat_flux.fine_grid = fine_grid;
  wall_heat_flux.marker = val_marker;
  wall_heat_flux.target = CustomBoundaryHeatFlux[val_marker];

  SetMultiGridMarkerQuantity(fine_grid, val_marker, wall_heat_flux);
}

void CMultiGridGeometry::SetMultiGridWallTemperature(const CGeometry* fine_grid, unsigned short val_marker) {
  struct {
    const CGeometry* fine_grid;
    unsigned short marker;
    su2double* target;

    su2double Get(unsigned long iVertex) const { return fine_grid->GetCustomBoundaryTemperature(marker, iVertex); }
    void Set(unsigned long iVertex, const su2double& val) const { target[iVertex] = val; }

  } wall_temperature;

  wall_temperature.fine_grid = fine_grid;
  wall_temperature.marker = val_marker;
  wall_temperature.target = CustomBoundaryTemperature[val_marker];

  SetMultiGridMarkerQuantity(fine_grid, val_marker, wall_temperature);
}

void CMultiGridGeometry::SetRestricted_GridVelocity(const CGeometry* fine_grid) {
  /*--- Loop over all coarse mesh points. ---*/
  SU2_OMP_FOR_STAT(roundUpDiv(nPoint, omp_get_max_threads()))
  for (unsigned long Point_Coarse = 0; Point_Coarse < nPoint; Point_Coarse++) {
    su2double Area_Parent = nodes->GetVolume(Point_Coarse);

    /*--- Initialize coarse grid velocity to zero. ---*/
    su2double Grid_Vel[3] = {0.0, 0.0, 0.0};

    /*--- Loop over all of the children for this coarse CV and compute
     a grid velocity based on the values in the child CVs (fine mesh). ---*/
    for (unsigned short iChild = 0; iChild < nodes->GetnChildren_CV(Point_Coarse); iChild++) {
      unsigned long Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChild);
      su2double Area_Child = fine_grid->nodes->GetVolume(Point_Fine);
      const su2double* Grid_Vel_Fine = fine_grid->nodes->GetGridVel(Point_Fine);
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Grid_Vel[iDim] += Grid_Vel_Fine[iDim] * Area_Child / Area_Parent;
    }

    /*--- Set the grid velocity for this coarse node. ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) nodes->SetGridVel(Point_Coarse, iDim, Grid_Vel[iDim]);
  }
  END_SU2_OMP_FOR
}

void CMultiGridGeometry::FindNormal_Neighbor(const CConfig* config) {
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) {
      for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
        auto iPoint = vertex[iMarker][iVertex]->GetNode();

        /*--- If the node belong to the domain ---*/
        if (nodes->GetDomain(iPoint)) {
          /*--- Compute closest normal neighbor ---*/
          unsigned long Point_Normal = 0;
          su2double* Normal = vertex[iMarker][iVertex]->GetNormal();
          su2double cos_max = -1.0;
          for (auto jPoint : nodes->GetPoints(iPoint)) {
            su2double scalar_prod = 0.0;
            su2double norm_vect = 0.0;
            su2double norm_Normal = 0.0;
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              su2double diff_coord = nodes->GetCoord(jPoint, iDim) - nodes->GetCoord(iPoint, iDim);
              scalar_prod += diff_coord * Normal[iDim];
              norm_vect += diff_coord * diff_coord;
              norm_Normal += Normal[iDim] * Normal[iDim];
            }
            norm_vect = sqrt(norm_vect);
            norm_Normal = sqrt(norm_Normal);
            su2double cos_alpha = scalar_prod / (norm_vect * norm_Normal);

            /*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
            if (cos_alpha >= cos_max) {
              Point_Normal = jPoint;
              cos_max = cos_alpha;
            }
          }
          vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
        }
      }
    }
  }
}

su2double CMultiGridGeometry::ComputeLocalCurvature(const CGeometry* fine_grid, unsigned long iPoint,
                                                    unsigned short iMarker) const {
  /*--- Compute local curvature (maximum angle between adjacent face normals) at a boundary vertex.
        This is used to determine if agglomeration is safe based on a curvature threshold. ---*/

  /*--- Get the vertex index for this point on this marker ---*/
  long iVertex = fine_grid->nodes->GetVertex(iPoint, iMarker);
  if (iVertex < 0) return 0.0;  // Point not on this marker

  /*--- Get the normal at this vertex ---*/
  su2double Normal_i[MAXNDIM] = {0.0};
  fine_grid->vertex[iMarker][iVertex]->GetNormal(Normal_i);
  su2double Area_i = GeometryToolbox::Norm(int(nDim), Normal_i);

  if (Area_i < EPS) return 0.0;  // Skip degenerate vertices

  /*--- Normalize the normal ---*/
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Normal_i[iDim] /= Area_i;
  }

  /*--- Find maximum angle with neighboring vertices on the same marker ---*/
  su2double max_angle = 0.0;

  /*--- Loop over edges connected to this point ---*/
  for (unsigned short iEdge = 0; iEdge < fine_grid->nodes->GetnPoint(iPoint); iEdge++) {
    unsigned long jPoint = fine_grid->nodes->GetPoint(iPoint, iEdge);

    /*--- Check if neighbor is also on this marker ---*/
    long jVertex = fine_grid->nodes->GetVertex(jPoint, iMarker);
    if (jVertex < 0) continue;  // Not on this marker

    /*--- Get normal at neighbor vertex ---*/
    su2double Normal_j[MAXNDIM] = {0.0};
    fine_grid->vertex[iMarker][jVertex]->GetNormal(Normal_j);
    su2double Area_j = GeometryToolbox::Norm(int(nDim), Normal_j);

    if (Area_j < EPS) continue;  // Skip degenerate neighbor

    /*--- Normalize the neighbor normal ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Normal_j[iDim] /= Area_j;
    }

    /*--- Compute dot product: cos(angle) = n_i · n_j ---*/
    su2double dot_product = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      dot_product += Normal_i[iDim] * Normal_j[iDim];
    }

    /*--- Clamp to [-1, 1] to avoid numerical issues with acos ---*/
    dot_product = max(-1.0, min(1.0, dot_product));

    /*--- Compute angle in degrees ---*/
    su2double angle_rad = acos(dot_product);
    su2double angle_deg = angle_rad * 180.0 / PI_NUMBER;

    /*--- Track maximum angle ---*/
    max_angle = max(max_angle, angle_deg);
  }

  return max_angle;
}

void CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseCV, const CGeometry* fine_grid,
                                                  const CConfig* config, CMultiGridQueue& MGQueue_InnerCV) {
  /*--- Parameters ---*/
  const su2double ANGLE_THRESHOLD_DEG = 20.0; /*!< Stop a line if the direction deviates more than this. */
  const unsigned long MAX_LINE_LENGTH = config->GetMGOptions().MG_Implicit_Lines_MaxLength;
  const su2double cos_threshold = cos(ANGLE_THRESHOLD_DEG * PI_NUMBER / 180.0);
  const bool ISOTROPIC = config->GetMGOptions().MG_Implicit_Lines_Isotropic;

  const unsigned long nPointFine = fine_grid->GetnPoint();
  const unsigned long starting_Index_CoarseCV = Index_CoarseCV;
  const bool DEBUG_OUTPUT = (rank == MASTER_NODE);

  /*--- How many parallel lines one coarse CV may span tangential to the wall. In 2D a wall "face" is
   *    an edge with 2 end nodes, in 3D a quadrilateral with 4 corner nodes, which is the number of
   *    lines that must be bundled to coarsen by 2 in every wall-tangential direction. ---*/
  unsigned long max_group = config->GetMGOptions().MG_Implicit_Lines_Max_Group;
  if (max_group == 0) max_group = (nDim == 2) ? 2 : 4;

  /*==================================================================================================
   *  PHASE A - build the implicit lines.
   *
   *  Lines are grown one step at a time across ALL lines simultaneously rather than one line to
   *  completion at a time, and a node is claimed globally the moment any line takes it. Growing them
   *  one-at-a-time lets an early line run the full depth of the layer and consume nodes that a later,
   *  neighbouring line needed, so that later line terminates after a step or two; the lines then have
   *  wildly different lengths and cannot be bundled into columns of uniform depth. Advancing in
   *  lockstep makes all lines compete for each layer on equal terms, which on an extruded prismatic
   *  layer reproduces the mesh's own structure: every line reaches the same depth.
   *================================================================================================*/
  vector<vector<unsigned long>> lines;      /*!< lines[i] = [wall_node, interior_1, interior_2, ...] */
  vector<su2double> dir;                    /*!< Current marching direction of each line, nDim per line. */
  vector<char> claimed(nPointFine, 0);      /*!< Node already belongs to some line. */

  /*--- Nodes that sit on a boundary with a boundary condition on it. A line must not grow into one,
   *    because those nodes belong to the boundary agglomeration and a stack that absorbed one would
   *    straddle two boundaries. CPoint's Boundary flag cannot answer this on its own: it is set for
   *    every marker a node belongs to, SEND_RECEIVE included, so on a partitioned mesh it is also
   *    true for the ordinary interior nodes of the send fringe. Testing it directly would stop every
   *    line that reaches the fringe one layer short of the partition, leaving the top of those
   *    columns to isotropic agglomeration purely because of where the mesh was cut. ---*/
  vector<char> onPhysicalBoundary(nPointFine, 0);
  for (auto iPoint = 0ul; iPoint < nPointFine; ++iPoint) {
    if (!fine_grid->nodes->GetBoundary(iPoint)) continue;
    for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
      if (fine_grid->nodes->GetVertex(iPoint, iMarker) != -1) {
        onPhysicalBoundary[iPoint] = 1;
        break;
      }
    }
  }

  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    /*--- Seed only from viscous (no-slip) walls: those are the boundaries with a stretched layer
     *    above them. Seeding from farfield/inlet/outlet/symmetry would claim layer nodes before the
     *    wall lines could reach them. ---*/
    const auto bc = config->GetMarker_All_KindBC(iMarker);
    if (bc != HEAT_FLUX && bc != ISOTHERMAL && bc != CHT_WALL_INTERFACE && bc != SMOLUCHOWSKI_MAXWELL) continue;

    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if (!fine_grid->nodes->GetDomain(iPoint)) continue;
      if (fine_grid->nodes->GetAgglomerate(iPoint)) continue;
      if (claimed[iPoint]) continue; /*--- A node on two wall markers must seed only one line. ---*/

      const long ChildVertex = fine_grid->nodes->GetVertex(iPoint, iMarker);
      if (ChildVertex == -1) continue;
      su2double Normal[MAXNDIM] = {0.0};
      fine_grid->vertex[iMarker][ChildVertex]->GetNormal(Normal);
      const su2double nrm = GeometryToolbox::Norm(nDim, Normal);
      if (nrm <= 0.0) continue;

      lines.push_back({iPoint});
      claimed[iPoint] = 1;
      for (unsigned short d = 0; d < nDim; ++d) dir.push_back(Normal[d] / nrm);
    }
  }

  if (lines.empty()) return;

  vector<char> growing(lines.size(), 1);
  for (bool any_grew = true; any_grew;) {
    any_grew = false;
    for (unsigned long li = 0; li < lines.size(); ++li) {
      if (!growing[li]) continue;
      if (lines[li].size() >= MAX_LINE_LENGTH) {
        growing[li] = 0;
        continue;
      }
      const auto current = lines[li].back();
      su2double best_dot = -2.0;
      unsigned long best_neighbor = ULONG_MAX;

      for (auto jPoint : fine_grid->nodes->GetPoints(current)) {
        /*--- Halo nodes stay out: their parent is dictated by the rank that owns them and arrives
         *    through the MPI relay, so a line claiming one here would fight that assignment. A line
         *    therefore still ends at the partition itself, but now only there, instead of one layer
         *    earlier at the fringe of owned nodes. ---*/
        if (!fine_grid->nodes->GetDomain(jPoint)) continue;
        if (onPhysicalBoundary[jPoint]) continue;
        if (fine_grid->nodes->GetAgglomerate(jPoint)) continue;
        if (claimed[jPoint]) continue;

        su2double vec[MAXNDIM] = {0.0};
        GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jPoint), fine_grid->nodes->GetCoord(current), vec);
        const su2double len = GeometryToolbox::Norm(nDim, vec);
        if (len <= 0.0) continue;
        for (unsigned short d = 0; d < nDim; ++d) vec[d] /= len;

        const su2double dot = GeometryToolbox::DotProduct(nDim, vec, &dir[li * nDim]);
        if (dot > best_dot) {
          best_dot = dot;
          best_neighbor = jPoint;
        }
      }

      if (best_neighbor == ULONG_MAX || best_dot < cos_threshold) {
        growing[li] = 0;
        continue;
      }

      su2double step[MAXNDIM] = {0.0};
      GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(best_neighbor), fine_grid->nodes->GetCoord(current),
                                step);
      const su2double slen = GeometryToolbox::Norm(nDim, step);
      if (slen <= 0.0) {
        growing[li] = 0;
        continue;
      }
      for (unsigned short d = 0; d < nDim; ++d) dir[li * nDim + d] = step[d] / slen;

      lines[li].push_back(best_neighbor);
      claimed[best_neighbor] = 1;
      any_grew = true;
    }
  }

  /*--- A line needs at least one interior node to contribute anything. Drop the rest and release
   *    their wall seed so ordinary boundary agglomeration can treat it normally. ---*/
  {
    vector<vector<unsigned long>> kept;
    kept.reserve(lines.size());
    for (auto& L : lines) {
      if (L.size() >= 2)
        kept.push_back(std::move(L));
      else
        claimed[L[0]] = 0;
    }
    lines = std::move(kept);
  }
  if (lines.empty()) return;

  /*==================================================================================================
   *  PHASE B - partition the lines into bundles.
   *
   *  Every line must end up in exactly one bundle, and a bundle must be a compact patch on the wall:
   *  in 3D the four lines rising from the corners of one wall quadrilateral, in 2D the two lines from
   *  the ends of one wall edge. Selecting, for each line independently, some set of neighbours to
   *  merge with does not do this - the relation is not symmetric, so line 1 claiming {2,3,4} does not
   *  stop line 2 from claiming {1,3,5}, and the bundles overlap and fight over nodes.
   *
   *  Building the partition by repeated pairwise matching avoids that by construction. One matching
   *  round pairs adjacent lines into 2-bundles (the wall edge); a second round pairs adjacent
   *  2-bundles into 4-bundles (the wall quadrilateral). Each round is a matching, so membership is
   *  mutually exclusive at every stage and the result is a true partition. It also needs nothing but
   *  point-to-point connectivity, so it works identically on every multigrid level - boundary face
   *  connectivity does not exist on agglomerated grids, so a literal "same quadrilateral" test would
   *  only ever work for the first coarsening.
   *
   *  Lines may only be bundled when their wall nodes carry the same set of physical markers, so that
   *  a bundle never straddles a boundary-condition change (the same rule ordinary agglomeration uses:
   *  ridges merge only with ridges, valleys only with valleys).
   *================================================================================================*/
  const unsigned long nLines = lines.size();

  /*--- Marker signature of each line's wall node. ---*/
  vector<vector<unsigned short>> sig(nLines);
  for (unsigned long li = 0; li < nLines; ++li) {
    for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
      if (fine_grid->nodes->GetVertex(lines[li][0], iMarker) != -1) sig[li].push_back(iMarker);
    }
    sort(sig[li].begin(), sig[li].end());
  }

  /*--- Line adjacency, inherited from the wall nodes' mesh connectivity. ---*/
  unordered_map<unsigned long, unsigned long> lineOfWallNode;
  lineOfWallNode.reserve(nLines);
  for (unsigned long li = 0; li < nLines; ++li) lineOfWallNode[lines[li][0]] = li;

  vector<vector<unsigned long>> adj(nLines);
  for (unsigned long li = 0; li < nLines; ++li) {
    for (auto jPoint : fine_grid->nodes->GetPoints(lines[li][0])) {
      const auto it = lineOfWallNode.find(jPoint);
      if (it != lineOfWallNode.end() && it->second != li) adj[li].push_back(it->second);
    }
  }

  /*--- Round 1: match adjacent lines into pairs. ---*/
  vector<vector<unsigned long>> bundles;
  vector<long> bundleOf(nLines, -1);
  bundles.reserve(nLines);
  for (unsigned long li = 0; li < nLines; ++li) {
    if (bundleOf[li] >= 0) continue;
    const long b = static_cast<long>(bundles.size());
    bundles.push_back({li});
    bundleOf[li] = b;
    if (max_group < 2) continue;
    for (auto lj : adj[li]) {
      if (bundleOf[lj] >= 0 || sig[lj] != sig[li]) continue;
      bundles[b].push_back(lj);
      bundleOf[lj] = b;
      break;
    }
  }

  /*--- Further rounds: merge adjacent bundles while they still fit. In 3D this turns pairs into
   *    quadrilaterals; in 2D max_group is 2 so it only absorbs leftover singletons. ---*/
  for (bool changed = true; changed;) {
    changed = false;

    vector<vector<long>> badj(bundles.size());
    for (unsigned long li = 0; li < nLines; ++li)
      for (auto lj : adj[li])
        if (bundleOf[li] != bundleOf[lj]) badj[bundleOf[li]].push_back(bundleOf[lj]);
    for (auto& v : badj) {
      sort(v.begin(), v.end());
      v.erase(unique(v.begin(), v.end()), v.end());
    }

    vector<char> consumed(bundles.size(), 0);
    vector<vector<unsigned long>> merged;
    merged.reserve(bundles.size());
    for (unsigned long b = 0; b < bundles.size(); ++b) {
      if (consumed[b]) continue;
      consumed[b] = 1;
      auto group = bundles[b];
      for (auto h : badj[b]) {
        if (consumed[h]) continue;
        if (group.size() + bundles[h].size() > max_group) continue;
        if (sig[bundles[h].front()] != sig[group.front()]) continue;
        consumed[h] = 1;
        group.insert(group.end(), bundles[h].begin(), bundles[h].end());
        changed = true;
        break;
      }
      merged.push_back(std::move(group));
    }

    bundles = std::move(merged);
    for (unsigned long b = 0; b < bundles.size(); ++b)
      for (auto li : bundles[b]) bundleOf[li] = static_cast<long>(b);
  }

  /*==================================================================================================
   *  PHASE C - extrude each bundle into a stack of coarse control volumes.
   *
   *  The bundle's wall nodes become one coarse CV, and each successive layer of the bundle's lines
   *  becomes the next, so the coarse grid inherits the layer structure of the fine grid and a line
   *  relaxation remains meaningful on it. Because Phase A made the lines node-disjoint and Phase B
   *  made the bundles a partition, no two bundles can ever contend for the same node, so a stack is
   *  never interrupted part way up.
   *
   *  The multigrid queue is deliberately not updated here: the sync loop that follows the boundary
   *  agglomeration removes every point already marked agglomerated, so removing them a second time
   *  from this function would be an error.
   *================================================================================================*/
  const unsigned long nBlock = ISOTROPIC ? 2 : 1;

  map<unsigned long, unsigned long> bundle_size_histogram;
  unsigned long nStacks = 0, nTruncated = 0;

  for (const auto& members : bundles) {
    bundle_size_histogram[members.size()]++;

    /*--- The wall CV. Claiming it here, before ordinary boundary agglomeration runs, is what keeps
     *    the whole stack aligned: the layer above a wall CV has exactly the same footprint. ---*/
    bool valid = true;
    for (auto li : members)
      if (!GeometricalCheck(lines[li][0], fine_grid, config)) valid = false;
    if (!valid) continue;

    for (unsigned long c = 0; c < members.size(); ++c) {
      const auto p = lines[members[c]][0];
      fine_grid->nodes->SetParent_CV(p, Index_CoarseCV);
      nodes->SetChildren_CV(Index_CoarseCV, c, p);
    }
    nodes->SetnChildren_CV(Index_CoarseCV, static_cast<unsigned short>(members.size()));
    Index_CoarseCV++;
    nStacks++;

    /*--- The interior layers, in lockstep, to the end of the shortest line in the bundle. ---*/
    unsigned long shortest = ULONG_MAX;
    for (auto li : members) shortest = std::min<unsigned long>(shortest, lines[li].size());

    unsigned long placed = 1;
    for (unsigned long first = 1; first + nBlock <= shortest; first += nBlock) {
      vector<unsigned long> group;
      group.reserve(members.size() * nBlock);
      for (auto li : members)
        for (unsigned long b = 0; b < nBlock; ++b) group.push_back(lines[li][first + b]);

      valid = true;
      for (auto p : group)
        if (!GeometricalCheck(p, fine_grid, config)) valid = false;
      if (!valid) break;

      for (unsigned long c = 0; c < group.size(); ++c) {
        fine_grid->nodes->SetParent_CV(group[c], Index_CoarseCV);
        nodes->SetChildren_CV(Index_CoarseCV, c, group[c]);
      }
      nodes->SetnChildren_CV(Index_CoarseCV, static_cast<unsigned short>(group.size()));
      Index_CoarseCV++;
      placed = first + nBlock;
    }

    /*--- Nodes above the shortest line, or above a block that did not divide evenly, are left to
     *    ordinary domain agglomeration. ---*/
    for (auto li : members) nTruncated += lines[li].size() - placed;
  }

  if (DEBUG_OUTPUT) {
    unsigned long nLineNodes = 0;
    for (const auto& L : lines) nLineNodes += L.size();
    cout << "  Implicit lines: " << nLines << " lines, " << nStacks << " stacks, bundle sizes ";
    for (const auto& h : bundle_size_histogram) cout << h.first << "x" << h.second << " ";
    cout << "\n  Coarse CVs from lines: " << (Index_CoarseCV - starting_Index_CoarseCV) << " covering "
         << (nLineNodes - nTruncated) << "/" << nLineNodes << " line nodes";
    if (nTruncated > 0) cout << " (" << nTruncated << " left to domain agglomeration)";
    cout << endl;
  }
}
