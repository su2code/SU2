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
    AgglomerateImplicitLines(Index_CoarseCV, fine_grid, config);
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

  /*--- STEP 2: Agglomerate the domain points.
   *
   *  A seed is grown one node at a time, and at every step the node that joins is the one most tied
   *  to what the CV already holds, i.e. the candidate sharing the most edges with its current members.
   *  Taking the seed's neighbours in whatever order the connectivity happens to list them, as this
   *  used to, ignores the shape of the result: on a hex mesh the seed has one neighbour per axis, so
   *  the first few are picked from different axes and the CV grows into a star. The leftovers around
   *  it then have to be swept up by later seeds, which is where the ragged pieces come from - a CV
   *  holding three nodes in one mesh plane and a single node in the next reads, in a cross-section
   *  through the second plane, as an isolated node sitting in the corner of an L.
   *
   *  Counting shared edges instead closes those shapes off by construction. Once a seed has taken two
   *  nodes along different axes, the node diagonally between them touches two members while every
   *  other option still touches one, so it wins and completes the square; the same argument then
   *  repeats one axis up and completes the cube. On a structured hex mesh the result is an exact
   *  2x2x2 block, which is what the implicit-line stacks already produce and what the isotropic
   *  region should match. Ties are settled by distance to the centroid of the current members, which
   *  keeps the growth compact on meshes with no preferred axis to reason about.
   *
   *  The candidate set grows with the CV rather than being fixed to the seed's own neighbours: the
   *  far corner of a cube is not adjacent to the seed, it only becomes reachable once the nodes
   *  between them have joined. That also makes SetSuitableNeighbors unnecessary here, since the nodes
   *  it used to supply are reached through ordinary edges as the frontier advances. ---*/

  /*--- Scratch shared by all seeds. The markers are cleared per CV, touching only what was used. ---*/
  vector<char> inCV(fine_grid->GetnPoint(), 0);
  vector<char> isCandidate(fine_grid->GetnPoint(), 0);
  vector<unsigned long> members, candidates;
  members.reserve(maxAgglomSize);

  auto iteration = 0ul;
  while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {
    const auto iPoint = MGQueue_InnerCV.NextCV();
    iteration++;

    /*--- If the element has not been previously agglomerated, belongs to the physical domain,
     and satisfies several geometrical criteria then the seed CV is accepted for agglomeration. ---*/

    if ((!fine_grid->nodes->GetAgglomerate(iPoint)) && (fine_grid->nodes->GetDomain(iPoint)) &&
        (GeometricalCheck(iPoint, fine_grid, config))) {
      members.clear();
      candidates.clear();
      unsigned short nChildren = 0;

      /*--- Take a node into the CV and let the frontier grow with it. ---*/
      auto addMember = [&](unsigned long CVPoint) {
        fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);
        nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
        nChildren++;

        if (fine_grid->nodes->GetAgglomerate_Indirect(CVPoint)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);

        /*--- Remove it from the queue and raise the priority of its neighbours. ---*/
        MGQueue_InnerCV.Update(CVPoint, fine_grid);

        members.push_back(CVPoint);
        inCV[CVPoint] = 1;

        for (auto jPoint : fine_grid->nodes->GetPoints(CVPoint)) {
          if (inCV[jPoint] || isCandidate[jPoint]) continue;
          if (fine_grid->nodes->GetAgglomerate(jPoint) || !fine_grid->nodes->GetDomain(jPoint)) continue;
          if (!GeometricalCheck(jPoint, fine_grid, config)) continue;
          isCandidate[jPoint] = 1;
          candidates.push_back(jPoint);
        }
      };

      addMember(iPoint);

      while (nChildren < maxAgglomSize) {
        /*--- Centroid of what the CV holds so far, used only to break ties. ---*/
        su2double centroid[MAXNDIM] = {0.0};
        for (auto jPoint : members) {
          const auto* coord = fine_grid->nodes->GetCoord(jPoint);
          for (auto iDim = 0u; iDim < nDim; iDim++) centroid[iDim] += coord[iDim] / su2double(members.size());
        }

        unsigned long best = std::numeric_limits<unsigned long>::max();
        unsigned short best_shared = 0;
        su2double best_dist = std::numeric_limits<su2double>::max();

        for (auto CVPoint : candidates) {
          if (inCV[CVPoint]) continue;

          unsigned short shared = 0;
          for (auto jPoint : fine_grid->nodes->GetPoints(CVPoint)) shared += inCV[jPoint];

          const su2double dist =
              GeometryToolbox::SquaredDistance(nDim, fine_grid->nodes->GetCoord(CVPoint), centroid);

          if ((shared > best_shared) || ((shared == best_shared) && (dist < best_dist))) {
            best = CVPoint;
            best_shared = shared;
            best_dist = dist;
          }
        }

        if (best == std::numeric_limits<unsigned long>::max()) break;
        addMember(best);
      }

      for (auto jPoint : members) inCV[jPoint] = 0;
      for (auto jPoint : candidates) isCandidate[jPoint] = 0;

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

  /*--- The pass above only rescues a coarse CV that has a single coarse neighbor, i.e. one enclosed
   entirely within another. It does nothing for a CV that was left with a single fine-grid child but
   sits between several coarse neighbors, which happens routinely with the implicit-line stacks above:
   a bundle whose members reach the top of the boundary layer at slightly different heights (graded
   spacing, local curvature) narrows one line at a time (see AgglomerateImplicitLines, Phase C), and
   the fine node just past where a line dropped out is left to fend for itself. It usually has no
   unclaimed neighbor left to agglomerate with by the time the ordinary domain pass (STEP 2) reaches
   it, so it becomes a coarse CV of its own, surrounded on multiple sides by CVs it cannot join under
   the single-neighbor rule above. Left alone, this is exactly the "hole" that turns what should be a
   compact block into an L shape with an isolated singleton sitting in the missing corner. Merge such
   a CV into whichever neighbor currently has the fewest children (and still has room under
   maxAgglomSize): the smallest neighbor is the one most likely to be the under-filled block the
   singleton belongs with, e.g. the 3-node L that this merge completes into a 4-node square, and
   folding into it keeps the coarse grid from accumulating disproportionately large CVs.

   Unlike the pass above this one deliberately does not skip CVs that border another rank. That guard
   exists because "has exactly one coarse neighbour" is measured on connectivity this rank cannot yet
   see in full, so it misfires on a CV whose other neighbours merely live across the partition. The
   trigger here is a child count, which is complete locally: a CV owning one fine point owns one fine
   point no matter how the mesh was cut. The merge is local too - every candidate neighbour comes from
   the connectivity built above and is therefore an owned CV, every child is an owned point, and the
   MPI relay below broadcasts the result afterwards, so the ranks stay consistent. Requiring at least
   two local neighbours below also rules out the one case that could not be repaired locally, a CV
   whose neighbours are all halo: a halo CV mirrors another rank's decision and must not be added to
   here. Keeping the guard cost 455 of 473 unrepaired singletons on a four-rank 3D bump. ---*/

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++) {
    if (nodes->GetnChildren_CV(iCoarsePoint) != 1) continue;
    if (nodes->GetnPoint(iCoarsePoint) <= 1) continue; /*--- Already handled above, or truly islanded. ---*/

    /*--- The smallest neighbour is the one most likely to be an under-filled block, and taking the
     smallest also means a neighbour with room is always preferred over a full one: anything below
     maxAgglomSize necessarily has fewer children than anything at it. When every neighbour is full
     the smallest is taken anyway, one child over the limit. maxAgglomSize is a quality target rather
     than a capacity - Children_CV grows on demand - and a CV of nine is a far smaller blemish on the
     coarse grid than the CV of one it removes. ---*/
    unsigned long best_neighbor = std::numeric_limits<unsigned long>::max();
    unsigned short best_nChildren = 0;
    for (auto jCoarsePoint : nodes->GetPoints(iCoarsePoint)) {
      const auto nChildren_Neighbor = nodes->GetnChildren_CV(jCoarsePoint);
      /*--- Skip neighbors already emptied by an earlier merge in this same pass. ---*/
      if (nChildren_Neighbor == 0) continue;
      if ((best_neighbor == std::numeric_limits<unsigned long>::max()) || (nChildren_Neighbor < best_nChildren)) {
        best_nChildren = nChildren_Neighbor;
        best_neighbor = jCoarsePoint;
      }
    }
    if (best_neighbor == std::numeric_limits<unsigned long>::max()) continue; /*--- Every neighbor was emptied. ---*/

    const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, 0);
    nodes->SetChildren_CV(best_neighbor, best_nChildren, iFinePoint);
    nodes->SetnChildren_CV(best_neighbor, best_nChildren + 1);
    fine_grid->nodes->SetParent_CV(iFinePoint, best_neighbor);
    nodes->SetnChildren_CV(iCoarsePoint, 0);
  }

  /*--- Compact the coarse numbering. Both repair passes empty a control volume by moving its children
   elsewhere, but nPointDomain was fixed before them, so an emptied index survives as a control volume
   that owns no fine points at all. It ends up with EPS volume, (0,0,0) coordinates and, because
   neighbours are derived through children, no neighbours either. Nothing is then adjacent to it on the
   next level, so it cannot be agglomerated there and becomes a one-child CV again, and so on down the
   hierarchy: the repair would otherwise trade a bad CV on this level for a degenerate one on every
   level above. Squeezing the empty slots out keeps the CV count honest and stops that propagation.

   Renumbering is safe here because only three things refer to coarse indices at this point - the
   children lists, the indirect-agglomeration flags, and the fine grid's parent indices - and the MPI
   relay below has not run yet, so no other rank has seen these numbers. Only parents of *owned* fine
   points are remapped; halo points are assigned by the relay afterwards. Every surviving parent keeps
   a valid new index because a CV that still has a child is never removed. ---*/

  {
    constexpr auto NO_INDEX = std::numeric_limits<unsigned long>::max();
    vector<unsigned long> newIndex(nPointDomain, NO_INDEX);
    unsigned long nKept = 0;
    for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++)
      if (nodes->GetnChildren_CV(iCoarsePoint) > 0) newIndex[iCoarsePoint] = nKept++;

    if (nKept < nPointDomain) {
      /*--- Move the survivors down. newIndex[i] <= i, so a forward sweep never lands on a CV that has
       not been moved yet. ---*/
      for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++) {
        const auto iNew = newIndex[iCoarsePoint];
        if ((iNew == NO_INDEX) || (iNew == iCoarsePoint)) continue;
        const auto nChildren = nodes->GetnChildren_CV(iCoarsePoint);
        for (auto iChildren = 0u; iChildren < nChildren; iChildren++)
          nodes->SetChildren_CV(iNew, iChildren, nodes->GetChildren_CV(iCoarsePoint, iChildren));
        nodes->SetnChildren_CV(iNew, nChildren);
        nodes->SetAgglomerate_Indirect(iNew, nodes->GetAgglomerate_Indirect(iCoarsePoint));
      }
      for (auto iCoarsePoint = nKept; iCoarsePoint < nPointDomain; iCoarsePoint++)
        nodes->SetnChildren_CV(iCoarsePoint, 0);

      for (auto iFinePoint = 0ul; iFinePoint < fine_grid->GetnPointDomain(); iFinePoint++) {
        const auto iParent = fine_grid->nodes->GetParent_CV(iFinePoint);
        if (iParent != NO_INDEX) fine_grid->nodes->SetParent_CV(iFinePoint, newIndex[iParent]);
      }

      nPointDomain = nKept;
      nPoint = nKept;
      Index_CoarseCV = nKept;
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

CMultiGridGeometry::CNodeStiffness CMultiGridGeometry::ComputeNodeStiffness(const CGeometry* fine_grid) const {
  /*--- Strength of the coupling across the dual face between a node and one of its neighbours. For a cell
   *    of streamwise size dx and wall-normal size dy this is 1/dy across the wall-normal face and 1/dx
   *    across the tangential one, so the ratio of the largest weight at a node to the smallest is the
   *    local cell aspect ratio. That makes the aspect ratio available from the dual grid alone, which
   *    SetControlVolume builds on every multigrid level, whereas CGeometry::Aspect_Ratio exists only on
   *    MESH_0. The same quantity decides where LINELET preconditioner lines stop, in GetLineletInfo.
   *
   *    Measuring the whole grid once keeps this out of the line-growth loop, which previously rescanned
   *    every neighbour of a node to find its weakest edge on every step of every line. ---*/
  const auto nPointFine = fine_grid->GetnPoint();

  CNodeStiffness stiff;
  stiff.wMin.assign(nPointFine, 0.0);
  stiff.wMax.assign(nPointFine, 0.0);
  stiff.jStiffest.assign(nPointFine, std::numeric_limits<unsigned long>::max());

  for (auto iPoint = 0ul; iPoint < nPointFine; ++iPoint) {
    su2double wmin = std::numeric_limits<su2double>::max(), wmax = 0.0;
    auto jStiffest = std::numeric_limits<unsigned long>::max();

    for (auto iNeigh = 0u; iNeigh < fine_grid->nodes->GetnPoint(iPoint); ++iNeigh) {
      const auto jPoint = fine_grid->nodes->GetPoint(iPoint, iNeigh);
      const auto iEdge = fine_grid->nodes->GetEdge(iPoint, iNeigh);
      const su2double area = GeometryToolbox::Norm(nDim, fine_grid->edges->GetNormal(iEdge));
      const su2double w =
          0.5 * area * (1.0 / fine_grid->nodes->GetVolume(iPoint) + 1.0 / fine_grid->nodes->GetVolume(jPoint));
      if (w > wmax) {
        wmax = w;
        jStiffest = jPoint;
      }
      wmin = std::min(wmin, w);
    }

    /*--- A node with no neighbours keeps the zeroed defaults, so AspectRatio reads 1. ---*/
    if (jStiffest != std::numeric_limits<unsigned long>::max()) {
      stiff.wMin[iPoint] = wmin;
      stiff.wMax[iPoint] = wmax;
      stiff.jStiffest[iPoint] = jStiffest;
    }
  }
  return stiff;
}

vector<vector<unsigned long>> CMultiGridGeometry::BuildImplicitLines(const CGeometry* fine_grid, const CConfig* config,
                                                                    const CNodeStiffness& stiff) const {
  /*--- Stop a line where its direction deviates by more than this from the previous step. ---*/
  constexpr su2double ANGLE_THRESHOLD_DEG = 20.0;
  /*--- Fraction of a marker's nodes that must sit in a layer before the whole marker may seed. ---*/
  constexpr su2double QUALIFIED_FRACTION = 0.5;

  const su2double cos_threshold = cos(ANGLE_THRESHOLD_DEG * PI_NUMBER / 180.0);
  const unsigned long MAX_LINE_LENGTH = config->GetMGOptions().MG_Implicit_Lines_MaxLength;
  const su2double MIN_AR = config->GetMGOptions().MG_Implicit_Lines_Min_AR;
  const bool USE_AR = (MIN_AR > 1.0);

  const auto nPointFine = fine_grid->GetnPoint();
  const auto nMarkerFine = fine_grid->GetnMarker();
  constexpr auto NO_POINT = std::numeric_limits<unsigned long>::max();

  vector<vector<unsigned long>> lines;
  vector<std::array<su2double, MAXNDIM>> dir; /*!< Current marching direction of each line. */
  vector<char> claimed(nPointFine, 0);

  auto isWall = [&](unsigned short bc) {
    return (bc == HEAT_FLUX) || (bc == ISOTHERMAL) || (bc == CHT_WALL_INTERFACE) || (bc == SMOLUCHOWSKI_MAXWELL);
  };

  /*--- Nodes on a boundary that carries a boundary condition. A line must not grow into one: those
   *    nodes belong to the boundary agglomeration and a stack absorbing one would straddle two
   *    boundaries. CPoint's Boundary flag cannot answer this, as it is also set by SEND_RECEIVE, so on
   *    a partitioned mesh it is true for ordinary interior nodes of the send fringe and every line
   *    would stop one layer short of the partition. Walking the markers' own vertex lists is both
   *    exact and cheaper than testing every point against every marker. ---*/
  vector<char> onPhysicalBoundary(nPointFine, 0);
  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++)
      onPhysicalBoundary[fine_grid->vertex[iMarker][iVertex]->GetNode()] = 1;
  }

  /*--- Unit normal of a boundary at a vertex, false if the marker does not reach iPoint. ---*/
  auto vertexNormal = [&](unsigned long iPoint, unsigned short iMarker, su2double* unitNormal) {
    const long ChildVertex = fine_grid->nodes->GetVertex(iPoint, iMarker);
    if (ChildVertex == -1) return false;
    fine_grid->vertex[iMarker][ChildVertex]->GetNormal(unitNormal);
    const su2double nrm = GeometryToolbox::Norm(nDim, unitNormal);
    if (nrm <= 0.0) return false;
    for (unsigned short d = 0; d < nDim; ++d) unitNormal[d] /= nrm;
    return true;
  };

  /*--- True if the mesh at iPoint is stretched along the boundary normal, i.e. this boundary has a
   *    layer growing off it the way a viscous wall does. On the side planes of a bump the mesh is just
   *    as stretched, but in a direction running ALONG the plane, and those nodes belong to the wall. ---*/
  auto hasLayerNormalTo = [&](unsigned long iPoint, const su2double* unitNormal) {
    const auto jStiffest = stiff.jStiffest[iPoint];
    if (jStiffest == NO_POINT) return false;
    if (stiff.AspectRatio(iPoint) < MIN_AR) return false;

    su2double vec[MAXNDIM] = {0.0};
    GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jStiffest), fine_grid->nodes->GetCoord(iPoint), vec);
    const su2double len = GeometryToolbox::Norm(nDim, vec);
    if (len <= 0.0) return false;
    for (unsigned short d = 0; d < nDim; ++d) vec[d] /= len;
    return fabs(GeometryToolbox::DotProduct(nDim, vec, unitNormal)) >= cos_threshold;
  };

  /*--- Seed a line at every eligible node of a marker. ---*/
  auto seedMarker = [&](unsigned short iMarker, bool requireLayer) {
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if (!fine_grid->nodes->GetDomain(iPoint)) continue;
      if (fine_grid->nodes->GetAgglomerate(iPoint)) continue;
      if (claimed[iPoint]) continue; /*--- A node on two markers must seed only one line. ---*/

      su2double Normal[MAXNDIM] = {0.0};
      if (!vertexNormal(iPoint, iMarker, Normal)) continue;
      if (requireLayer && !hasLayerNormalTo(iPoint, Normal)) continue;

      lines.push_back({iPoint});
      claimed[iPoint] = 1;
      std::array<su2double, MAXNDIM> d0{};
      for (unsigned short d = 0; d < nDim; ++d) d0[d] = Normal[d];
      dir.push_back(d0);
    }
  };

  /*--- Viscous walls always carry a stretched layer, so they seed unconditionally. Running them first
   *    also settles nodes where a wall meets another boundary: the wall claims them. ---*/
  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++)
    if (isWall(config->GetMarker_All_KindBC(iMarker))) seedMarker(iMarker, false);

  /*--- Boundaries other than viscous walls that nevertheless carry a stretched layer normal to
   *    themselves, such as a symmetry plane laid in the same surface as a wall - the floor either side
   *    of a bump, or ahead of a flat plate's leading edge. Seeding only from walls leaves the mesh
   *    above those to isotropic agglomeration, so the coarse grid changes character across a line the
   *    fine grid does not single out.
   *
   *    The verdict is taken per marker, not per node: seeding isolated qualifying nodes on a marker
   *    that mostly does not qualify scatters one-line bundles, i.e. coarse CVs that do not coarsen
   *    tangentially at all. The two populations are far apart in practice - boundaries with a layer
   *    normal to them qualify at 100%, side planes, inlets and far fields at 13% and below - so any
   *    threshold near a half separates them. ---*/
  if (USE_AR) {
    /*--- Counted per configuration-file marker, not per local marker: ranks agree on neither the
     *    number nor the order of local markers, because each partition appends its own SEND_RECEIVE
     *    markers, so the same index means a different boundary elsewhere. ---*/
    const auto nMarkerCfg = config->GetnMarker_CfgFile();
    vector<unsigned long> nValid(nMarkerCfg, 0), nQualified(nMarkerCfg, 0);

    auto canSeed = [&](unsigned short iMarker) {
      const auto bc = config->GetMarker_All_KindBC(iMarker);
      /*--- Periodic boundaries are left out: the two halves are the same physical location under a
       *    transform and have their own matching, which a line running into one would disturb. ---*/
      return (bc != SEND_RECEIVE) && (bc != PERIODIC_BOUNDARY) && !isWall(bc);
    };

    vector<unsigned short> cfgOfMarker(nMarkerFine, 0);
    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++)
      if (canSeed(iMarker))
        cfgOfMarker[iMarker] = config->GetMarker_CfgFile_TagBound(config->GetMarker_All_TagBound(iMarker));

    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
      if (!canSeed(iMarker)) continue;
      for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
        const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
        if (!fine_grid->nodes->GetDomain(iPoint)) continue;
        su2double Normal[MAXNDIM] = {0.0};
        if (!vertexNormal(iPoint, iMarker, Normal)) continue;
        nValid[cfgOfMarker[iMarker]]++;
        if (hasLayerNormalTo(iPoint, Normal)) nQualified[cfgOfMarker[iMarker]]++;
      }
    }

    /*--- A marker is generally split over several ranks, so the verdict must be taken on all of it. ---*/
    if (nMarkerCfg > 0) {
      vector<unsigned long> tmp(nMarkerCfg);
      SU2_MPI::Allreduce(nValid.data(), tmp.data(), nMarkerCfg, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      nValid.swap(tmp);
      SU2_MPI::Allreduce(nQualified.data(), tmp.data(), nMarkerCfg, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      nQualified.swap(tmp);
    }

    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
      if (!canSeed(iMarker)) continue;
      const auto iCfg = cfgOfMarker[iMarker];
      if (nValid[iCfg] == 0) continue;
      if (su2double(nQualified[iCfg]) < QUALIFIED_FRACTION * su2double(nValid[iCfg])) continue;
      seedMarker(iMarker, true);
    }
  }

  /*--- Grow every line one step per sweep rather than one line to completion at a time, claiming a
   *    node globally the moment any line takes it. Growing them one at a time lets an early line run
   *    the full depth of the layer and consume nodes a neighbouring line needed, so that line stops
   *    after a step or two and the lengths become too uneven to bundle into columns of uniform depth.
   *    In lockstep all lines compete for each layer on equal terms, which on an extruded prismatic
   *    layer reproduces the mesh's own structure. ---*/
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
      su2double best_dot = -2.0, best_dir[MAXNDIM] = {0.0};
      auto best = NO_POINT;
      unsigned short best_neigh = 0;

      for (auto iNeigh = 0u; iNeigh < fine_grid->nodes->GetnPoint(current); ++iNeigh) {
        const auto jPoint = fine_grid->nodes->GetPoint(current, iNeigh);
        /*--- Halo nodes stay out: their parent is dictated by the rank that owns them and arrives
         *    through the MPI relay, so a line claiming one would fight that assignment. ---*/
        if (!fine_grid->nodes->GetDomain(jPoint)) continue;
        if (onPhysicalBoundary[jPoint]) continue;
        if (fine_grid->nodes->GetAgglomerate(jPoint)) continue;
        if (claimed[jPoint]) continue;

        su2double vec[MAXNDIM] = {0.0};
        GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jPoint), fine_grid->nodes->GetCoord(current), vec);
        const su2double len = GeometryToolbox::Norm(nDim, vec);
        if (len <= 0.0) continue;
        for (unsigned short d = 0; d < nDim; ++d) vec[d] /= len;

        const su2double dot = GeometryToolbox::DotProduct(nDim, vec, dir[li].data());
        if (dot > best_dot) {
          best_dot = dot;
          best = jPoint;
          best_neigh = iNeigh;
          for (unsigned short d = 0; d < nDim; ++d) best_dir[d] = vec[d];
        }
      }

      if ((best == NO_POINT) || (best_dot < cos_threshold)) {
        growing[li] = 0;
        continue;
      }

      /*--- End the line where the mesh stops being stretched along it. Without this the only limits
       *    are the direction cone and MAX_LINE_LENGTH, so a line leaves the boundary layer and keeps
       *    stacking coarse CVs along a direction the fine grid does not single out. Taking the weight
       *    of this one edge over the weakest edge at the node keeps the measure directional: a mesh
       *    graded in the streamwise direction reads as stretched to an undirected min/max test even
       *    far from any wall. ---*/
      if (USE_AR && (stiff.wMin[current] > 0.0)) {
        const auto jPoint = fine_grid->nodes->GetPoint(current, best_neigh);
        const auto iEdge = fine_grid->nodes->GetEdge(current, best_neigh);
        const su2double area = GeometryToolbox::Norm(nDim, fine_grid->edges->GetNormal(iEdge));
        const su2double w =
            0.5 * area * (1.0 / fine_grid->nodes->GetVolume(current) + 1.0 / fine_grid->nodes->GetVolume(jPoint));
        if (w / stiff.wMin[current] < MIN_AR) {
          growing[li] = 0;
          continue;
        }
      }

      for (unsigned short d = 0; d < nDim; ++d) dir[li][d] = best_dir[d];
      lines[li].push_back(best);
      claimed[best] = 1;
      any_grew = true;
    }
  }

  /*--- A line needs at least one interior node to contribute anything. ---*/
  vector<vector<unsigned long>> kept;
  kept.reserve(lines.size());
  for (auto& L : lines)
    if (L.size() >= 2) kept.push_back(std::move(L));

  return kept;
}

vector<vector<unsigned long>> CMultiGridGeometry::BundleImplicitLines(const vector<vector<unsigned long>>& lines,
                                                                     const CGeometry* fine_grid,
                                                                     const CConfig* config,
                                                                     vector<vector<unsigned long>>& adj) const {
  /*--- Every line must end up in exactly one bundle, and a bundle must be a compact patch on the wall:
   *    in 3D the four lines rising from the corners of one wall quadrilateral, in 2D the two lines from
   *    the ends of one wall edge. Choosing, for each line independently, a set of neighbours to merge
   *    with does not do this - the relation is not symmetric, so line 1 claiming {2,3} does not stop
   *    line 2 claiming {1,4}, and the bundles overlap and fight over nodes.
   *
   *    Repeated pairwise matching avoids that by construction. One round pairs adjacent lines into the
   *    wall edge, a second pairs adjacent pairs into the wall quadrilateral. Each round is a matching,
   *    so membership stays mutually exclusive and the result is a true partition. It needs nothing but
   *    point-to-point connectivity, so it works identically on every multigrid level - boundary face
   *    connectivity does not survive agglomeration, so a literal "same quadrilateral" test would only
   *    ever work for the first coarsening. Two rounds reach 4, which is the 3D group size, so the
   *    number of rounds follows from max_group instead of iterating to a fixed point. ---*/
  const auto nLines = lines.size();
  unsigned long max_group = config->GetMGOptions().MG_Implicit_Lines_Max_Group;
  if (max_group == 0) max_group = (nDim == 2) ? 2 : 4;

  /*--- Marker signature of each line's boundary node, as a bitmask over the physical markers. Lines
   *    may only be bundled when these match, so a bundle never straddles a change of boundary
   *    condition - the rule ordinary agglomeration uses for ridges and valleys. ---*/
  const auto nMarkerFine = fine_grid->GetnMarker();
  vector<short> physBit(nMarkerFine, -1);
  unsigned nPhys = 0;
  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) physBit[iMarker] = static_cast<short>(nPhys++);

  const unsigned nWords = std::max(1u, (nPhys + 63u) / 64u);
  vector<uint64_t> sig(nLines * nWords, 0);
  for (unsigned long li = 0; li < nLines; ++li)
    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
      if (physBit[iMarker] < 0) continue;
      if (fine_grid->nodes->GetVertex(lines[li][0], iMarker) == -1) continue;
      const auto b = static_cast<unsigned>(physBit[iMarker]);
      sig[li * nWords + b / 64] |= (uint64_t(1) << (b % 64));
    }
  auto sameSig = [&](unsigned long la, unsigned long lb) {
    for (unsigned w = 0; w < nWords; ++w)
      if (sig[la * nWords + w] != sig[lb * nWords + w]) return false;
    return true;
  };

  /*--- Line adjacency, inherited from the boundary nodes' mesh connectivity. ---*/
  vector<long> lineOfWallNode(fine_grid->GetnPoint(), -1);
  for (unsigned long li = 0; li < nLines; ++li) lineOfWallNode[lines[li][0]] = static_cast<long>(li);

  adj.assign(nLines, {});
  for (unsigned long li = 0; li < nLines; ++li)
    for (auto jPoint : fine_grid->nodes->GetPoints(lines[li][0])) {
      const auto lj = lineOfWallNode[jPoint];
      if ((lj >= 0) && (static_cast<unsigned long>(lj) != li)) adj[li].push_back(static_cast<unsigned long>(lj));
    }

  vector<vector<unsigned long>> groups(nLines);
  vector<unsigned long> groupOf(nLines);
  for (unsigned long li = 0; li < nLines; ++li) {
    groups[li] = {li};
    groupOf[li] = li;
  }

  const unsigned nRounds = (max_group <= 1) ? 0 : ((max_group <= 2) ? 1 : 2);
  vector<std::pair<unsigned long, unsigned long>> shared;

  for (unsigned round = 0; round < nRounds; ++round) {
    vector<char> consumed(groups.size(), 0);
    vector<vector<unsigned long>> merged;
    merged.reserve(groups.size());

    for (unsigned long g = 0; g < groups.size(); ++g) {
      if (consumed[g]) continue;
      consumed[g] = 1;
      auto group = groups[g];

      /*--- Count how many line-to-line adjacencies this group shares with each candidate. Merging the
       *    candidate that shares the most keeps the patch square: a pair lying alongside this one
       *    touches it along its whole length and shares two adjacencies, whereas a pair continuing in
       *    the same direction touches at one end and shares one. Taking the first candidate that fits
       *    instead, as this did, produces a 1x4 strip of lines about as often as the mesh offers one,
       *    and those extrude into coarse CVs elongated in one wall-tangential direction. ---*/
      shared.clear();
      for (auto li : group)
        for (auto lj : adj[li]) {
          const auto h = groupOf[lj];
          if ((h == g) || consumed[h]) continue;
          if (group.size() + groups[h].size() > max_group) continue;
          if (!sameSig(groups[h].front(), group.front())) continue;

          bool seen = false;
          for (auto& s : shared)
            if (s.first == h) {
              s.second++;
              seen = true;
              break;
            }
          if (!seen) shared.emplace_back(h, 1);
        }

      auto bestH = std::numeric_limits<unsigned long>::max();
      unsigned long bestShared = 0;
      for (const auto& s : shared)
        if (s.second > bestShared) {
          bestShared = s.second;
          bestH = s.first;
        }

      if (bestH != std::numeric_limits<unsigned long>::max()) {
        consumed[bestH] = 1;
        group.insert(group.end(), groups[bestH].begin(), groups[bestH].end());
      }
      merged.push_back(std::move(group));
    }

    groups = std::move(merged);
    for (unsigned long g = 0; g < groups.size(); ++g)
      for (auto li : groups[g]) groupOf[li] = g;
  }

  return groups;
}

void CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseCV, const CGeometry* fine_grid,
                                                  const CConfig* config) {
  const auto starting_Index_CoarseCV = Index_CoarseCV;

  const auto stiff = ComputeNodeStiffness(fine_grid);

  /*--- PHASE A and B. Neither may be skipped on a rank with no lines: PHASE A takes a collective to
   *    agree on which markers carry a layer, and the summary at the end takes another. ---*/
  const auto lines = BuildImplicitLines(fine_grid, config, stiff);

  vector<vector<unsigned long>> adj, bundles;
  if (!lines.empty()) bundles = BundleImplicitLines(lines, fine_grid, config, adj);

  /*==================================================================================================
   *  PHASE C - extrude each bundle into a stack of coarse control volumes.
   *
   *  The bundle's boundary nodes become one coarse CV and each successive layer becomes the next, so
   *  the coarse grid inherits the layer structure of the fine grid and a line relaxation stays
   *  meaningful on it. PHASE A made the lines node-disjoint and PHASE B made the bundles a partition,
   *  so no two bundles can contend for a node and a stack is never interrupted part way up.
   *
   *  How many fine layers go into one coarse CV decides what is coarsened. One layer coarsens only
   *  tangentially and leaves the wall-normal line intact, which is what a line-implicit smoother needs
   *  where the cells are stretched; two layers coarsen in every direction, which is what the far field
   *  wants. MG_IMPLICIT_LINES_ISO_AR switches between them by the local aspect ratio as the stack
   *  rises, so the same stack can start semi-coarsened at the wall and finish isotropic - and because
   *  each semi-coarsening halves the stretching, the switch happens by itself at whatever level the
   *  mesh stops being anisotropic, rather than being tied to a multigrid level. With it at 0 the old
   *  behaviour stands and MG_IMPLICIT_LINES_ISOTROPIC decides for the whole stack.
   *
   *  The multigrid queue is deliberately not touched here: the sync loop after the boundary
   *  agglomeration removes every point already marked agglomerated, so doing it here too would be an
   *  error.
   *================================================================================================*/
  const bool ISOTROPIC = config->GetMGOptions().MG_Implicit_Lines_Isotropic;
  const su2double ISO_AR = config->GetMGOptions().MG_Implicit_Lines_Iso_AR;
  const bool HYBRID = (ISO_AR > 0.0);

  /*--- Are the still-growing lines one connected patch on the wall? A set that is no longer connected
   *    would put two separated columns into a single CV. Scans PHASE B's adjacency, over at most
   *    max_group entries, so the quadratic form is cheap. ---*/
  auto isConnected = [&adj](const vector<unsigned long>& members, const vector<unsigned long>& act) {
    if (act.size() <= 1) return true;
    vector<char> seen(act.size(), 0);
    vector<unsigned long> stack{0};
    seen[0] = 1;
    unsigned long nSeen = 1;
    while (!stack.empty()) {
      const auto cur = stack.back();
      stack.pop_back();
      const auto& neighbors = adj[members[act[cur]]];
      for (unsigned long k = 0; k < act.size(); ++k) {
        if (seen[k]) continue;
        if (find(neighbors.begin(), neighbors.end(), members[act[k]]) != neighbors.end()) {
          seen[k] = 1;
          nSeen++;
          stack.push_back(k);
        }
      }
    }
    return nSeen == act.size();
  };

  unsigned long nStacks = 0, nTruncated = 0, nSemiCV = 0, nFullCV = 0;
  unsigned long histogram[9] = {0};
  vector<unsigned long> placed, active, group;

  for (const auto& members : bundles) {
    histogram[std::min<size_t>(members.size(), 8)]++;

    /*--- The boundary CV. Claiming it before ordinary boundary agglomeration runs is what keeps the
     *    stack aligned: the layer above has exactly the same footprint. ---*/
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

    /*--- The interior layers, in lockstep. A line that runs out stops contributing and the others
     *    carry on, so the stack narrows as it rises instead of being cut to its shortest member.
     *    Dropping below two lines would extrude a column one line wide, thinner than anything the
     *    domain pass would build there, so the stack ends instead; a bundle that only ever had one
     *    line is exempt, being one line wide by construction. ---*/
    const unsigned long minActive = std::min<unsigned long>(2, members.size());
    placed.assign(members.size(), 1);

    auto collectActive = [&](unsigned long first, unsigned long blk) {
      active.clear();
      for (unsigned long m = 0; m < members.size(); ++m)
        if (first + blk <= lines[members[m]].size()) active.push_back(m);
    };

    for (unsigned long first = 1;;) {
      unsigned long nBlock = ISOTROPIC ? 2 : 1;
      if (HYBRID) {
        /*--- Stay semi-coarsened while any line of this bundle is still in stretched mesh at this
         *    height; the wall-normal direction is shared by the whole stack. ---*/
        su2double arHere = 0.0;
        for (unsigned long m = 0; m < members.size(); ++m)
          if (first < lines[members[m]].size())
            arHere = std::max(arHere, stiff.AspectRatio(lines[members[m]][first]));
        nBlock = (arHere > ISO_AR) ? 1 : 2;
      }

      collectActive(first, nBlock);
      /*--- One layer short of a full block at the top, take it as a single layer rather than drop it. ---*/
      if ((nBlock == 2) && (active.size() < minActive)) {
        nBlock = 1;
        collectActive(first, nBlock);
      }

      if (active.size() < minActive) break;
      if (!isConnected(members, active)) break;

      group.clear();
      group.reserve(active.size() * nBlock);
      for (auto m : active)
        for (unsigned long b = 0; b < nBlock; ++b) group.push_back(lines[members[m]][first + b]);

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
      ((nBlock == 1) ? nSemiCV : nFullCV)++;

      for (auto m : active) placed[m] = first + nBlock;
      first += nBlock;
    }

    for (unsigned long m = 0; m < members.size(); ++m) nTruncated += lines[members[m]].size() - placed[m];
  }

  /*--- Summary over all ranks. Reporting rank 0's own lines, as this used to, makes a partitioned run
   *    look like a fraction of the mesh it is not, and hides how much of the layer the partitioning
   *    cost: a line stops at the partition, so the count of nodes left to ordinary agglomeration is
   *    the number to watch when adding ranks. Every rank must reach these collectives. ---*/
  unsigned long nLineNodes = 0;
  for (const auto& L : lines) nLineNodes += L.size();

  unsigned long local[6] = {lines.size(), nStacks, nLineNodes, nTruncated, nSemiCV, nFullCV};
  unsigned long total[6] = {0};
  SU2_MPI::Allreduce(local, total, 6, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  unsigned long localCV = Index_CoarseCV - starting_Index_CoarseCV, totalCV = 0;
  SU2_MPI::Allreduce(&localCV, &totalCV, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  unsigned long histTotal[9] = {0};
  SU2_MPI::Allreduce(histogram, histTotal, 9, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  if (rank == MASTER_NODE) {
    cout << "  Implicit lines: " << total[0] << " lines, " << total[1] << " stacks, bundle sizes ";
    for (unsigned s = 1; s <= 8; ++s)
      if (histTotal[s] > 0) cout << s << "x" << histTotal[s] << " ";
    cout << "\n  Coarse CVs from lines: " << totalCV << " covering " << (total[2] - total[3]) << "/" << total[2]
         << " line nodes";
    if (total[3] > 0) cout << " (" << total[3] << " left to domain agglomeration)";
    if (total[4] + total[5] > 0)
      cout << "\n  Stack layers: " << total[4] << " semi-coarsened, " << total[5] << " isotropic";
    cout << endl;
  }
}
