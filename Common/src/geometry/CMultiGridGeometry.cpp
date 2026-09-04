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

  /*--- STEP 0: pave the domain with advancing fronts rising from the boundaries, wall CV included.
   *    This runs before the general boundary agglomeration below so that the wall control volume and
   *    the layers stacked on top of it share one footprint; letting the general scheme claim the wall
   *    first would fix a footprint chosen without any knowledge of the fronts, and the stack above it
   *    could then only be misaligned with its own base. Everything it claims is
   *    already marked agglomerated, so the boundary and interior passes below simply skip it. ---*/
  const auto starting_idx_lines_DBG = Index_CoarseCV;
  if (config->GetMGOptions().MG_Implicit_Lines) {
    AgglomerateImplicitLines(Index_CoarseCV, fine_grid, config);
  }
  const auto idx_after_lines_DBG = Index_CoarseCV;

  /*--- Points carrying a physical boundary condition. SEND_RECEIVE is not one: it only records that
   *    the point is mirrored on another rank. Used below to tell a genuine interior point, which a
   *    boundary CV may absorb, from a point on another boundary, which it may not. ---*/
  vector<char> onPhysBoundary(fine_grid->GetnPoint(), 0);
  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++)
      onPhysBoundary[fine_grid->vertex[iMarker][iVertex]->GetNode()] = 1;
  }
  /*--- Nodes where two different boundary conditions meet. The rule has to hold for every phase or
   *    the same node is treated one way by the paving and another here. In 2D the corner test below
   *    already refused them; in 3D a ridge of such nodes carries one identical marker PAIR all along
   *    it and would otherwise pair up with itself quite happily. ---*/
  const auto mixedBC = FindMixedBoundaryNodes(fine_grid, config);

  vector<unsigned long> bmembers;
  vector<unsigned long> nThicken_DBG(fine_grid->GetnMarker(), 0), nFlat_DBG(fine_grid->GetnMarker(), 0);

  /*--- Whether a boundary coarse CV may be grown into the interior, see MG_BOUNDARY_THICKEN_AR. The
   *    measurement is only needed when it can be, so a run that leaves the option off pays nothing. ---*/
  const su2double thickenAR = config->GetMGOptions().MG_Boundary_Thicken_AR;
  const bool THICKEN = (thickenAR > 0.0);
  const CNodeStiffness boundStiff = THICKEN ? ComputeNodeStiffness(fine_grid) : CNodeStiffness();

  /*--- True where the mesh carries a stretched layer running normal to this boundary, the situation a
   *    flat boundary CV exists to preserve. Both halves of the test matter, and testing the aspect
   *    ratio alone is what made a fixed threshold useless on a real mesh: the ratio is undirected, so
   *    a boundary lying in mesh that is merely graded ALONG itself - a symmetry plane with streamwise
   *    stretching, say - reads as strongly stretched and never gets thickened, even though nothing
   *    normal to it would be lost. Requiring the stiffest direction at the node to line up with the
   *    boundary normal separates the two: only a boundary with cells stacked against it qualifies. ---*/
  constexpr su2double THICKEN_ANGLE_DEG = 20.0;
  const su2double thickenCos = cos(THICKEN_ANGLE_DEG * PI_NUMBER / 180.0);

  auto boundaryHasLayer = [&](unsigned long iPoint, unsigned short iMarker) {
    const auto jStiffest = boundStiff.jStiffest[iPoint];
    if (jStiffest == std::numeric_limits<unsigned long>::max()) return false;
    if (boundStiff.AspectRatio(iPoint) < thickenAR) return false;

    const long iVertex = fine_grid->nodes->GetVertex(iPoint, iMarker);
    if (iVertex == -1) return false;
    su2double normal[MAXNDIM] = {0.0};
    fine_grid->vertex[iMarker][iVertex]->GetNormal(normal);
    const su2double nrm = GeometryToolbox::Norm(nDim, normal);
    if (nrm <= 0.0) return false;

    su2double vec[MAXNDIM] = {0.0};
    GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jStiffest), fine_grid->nodes->GetCoord(iPoint), vec);
    const su2double len = GeometryToolbox::Norm(nDim, vec);
    if (len <= 0.0) return false;

    su2double dot = 0.0;
    for (unsigned short d = 0; d < nDim; ++d) dot += (vec[d] / len) * (normal[d] / nrm);
    return fabs(dot) >= thickenCos;
  };

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
        bmembers.clear();
        bmembers.push_back(iPoint);
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

        /*--- ...and so is a node where two markers of DIFFERENT type meet, whatever the count and
         *    whatever the dimension. A coarse CV holding such a node would average two conditions the
         *    fine grid applies separately. ---*/
        if (mixedBC[iPoint]) agglomerate_seed = false;

        /*--- If the seed (parent) can be agglomerated, we try to agglomerate connected childs to the parent ---*/
        /*--- Note that in 2D we allow a maximum of 4 nodes to be agglomerated ---*/

        if (agglomerate_seed) {
          /*--- Now we do a sweep over all the nodes that surround the seed point ---*/

          for (auto CVPoint : fine_grid->nodes->GetPoints(iPoint)) {
            /*--- The new point can be agglomerated ---*/

            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config, mixedBC)) {
              /*--- We set the value of the parent ---*/

              fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

              /*--- We set the value of the child ---*/

              nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
              nChildren++;
              bmembers.push_back(CVPoint);
              /*--- In 2D, we agglomerate exactly 2 nodes if the nodes are on the line edge. ---*/
              if ((nDim == 2) && (counter == 1)) break;
              /*--- In 3D, we agglomerate exactly 2 nodes if the nodes are on the surface edge. ---*/
              if ((nDim == 3) && (counter == 2)) break;
              /*--- Apply maxAgglomSize limit for 3D internal boundary face nodes (counter==1 in 3D). ---*/
              if (nChildren >= maxAgglomSize) break;
            }
          }

          /*--- Only take into account indirect neighbors for 3D faces, not 2D. The size test has to be
           *    made on entry as well: the sweep above leaves with the CV exactly full whenever it hit
           *    the limit, and without a guard here the first indirect candidate pushed it to nine
           *    children, after which the equality test below could never match again and the CV grew
           *    without a bound at all. ---*/
          if ((nDim == 3) && (nChildren < maxAgglomSize)) {
            Suitable_Indirect_Neighbors.clear();

            if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint))
              SetSuitableNeighbors(Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);

            /*--- Now we do a sweep over all the indirect nodes that can be added ---*/

            for (auto CVPoint : Suitable_Indirect_Neighbors) {
              /*--- The new point can be agglomerated ---*/

              if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config, mixedBC)) {
                /*--- We set the value of the parent ---*/

                fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

                /*--- We set the indirect agglomeration information of the corse point
                based on its children in the fine grid. ---*/

                if (fine_grid->nodes->GetAgglomerate_Indirect(CVPoint))
                  nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);

                /*--- We set the value of the child ---*/

                nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
                nChildren++;
                bmembers.push_back(CVPoint);
                /*--- Apply maxAgglomSize limit for 3D internal boundary face nodes. ---*/
                if (nChildren >= maxAgglomSize) break;
              }
            }
          }

          /*--- Thicken the surface patch into the interior. Everything above only ever considers
           *    candidates that lie on the boundary themselves, because SetBoundAgglomeration refuses
           *    an interior point outright, so a boundary coarse CV comes out as a film one fine cell
           *    thick: 2x2 nodes at best on a surface, 2 on a ridge, never the 2x2x2 block the domain
           *    pass builds everywhere else. That is a coarse CV of 4 nodes where 8 were available,
           *    and since every boundary of the mesh is covered in them they make up a large share of
           *    the coarse grid by count while holding very few nodes each.
           *
           *    Growing into the interior fixes that without touching which surface nodes belong
           *    together: the footprint on the boundary is already decided above, this only adds the
           *    layer underneath it. Candidates are restricted to genuinely interior points, so the CV
           *    still cannot straddle two boundary conditions, and the node that joins is the one
           *    sharing the most faces with what the CV already holds - the same rule STEP 2 uses,
           *    which is what makes it close into blocks rather than grow into stars. ---*/
          const bool thickenThis = THICKEN && !boundaryHasLayer(iPoint, iMarker);
          (thickenThis ? nThicken_DBG : nFlat_DBG)[iMarker]++;

          while (thickenThis && (nChildren < maxAgglomSize)) {
            auto best = std::numeric_limits<unsigned long>::max();
            unsigned short best_shared = 0;

            for (auto mPoint : bmembers) {
              for (auto jPoint : fine_grid->nodes->GetPoints(mPoint)) {
                if (fine_grid->nodes->GetAgglomerate(jPoint)) continue;
                if (!fine_grid->nodes->GetDomain(jPoint)) continue;
                if (onPhysBoundary[jPoint]) continue;
                if (!GeometricalCheck(jPoint, fine_grid, config)) continue;

                unsigned short shared = 0;
                for (auto kPoint : fine_grid->nodes->GetPoints(jPoint))
                  shared += (find(bmembers.begin(), bmembers.end(), kPoint) != bmembers.end());

                if (shared > best_shared) {
                  best_shared = shared;
                  best = jPoint;
                }
              }
            }

            if (best == std::numeric_limits<unsigned long>::max()) break;

            fine_grid->nodes->SetParent_CV(best, Index_CoarseCV);
            if (fine_grid->nodes->GetAgglomerate_Indirect(best)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);
            nodes->SetChildren_CV(Index_CoarseCV, nChildren, best);
            nChildren++;
            bmembers.push_back(best);
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
   *  region should match.
   *
   *  Ties are settled by distance to the centroid, but that distance has to be measured in cells and
   *  not in metres. Until the CV holds an L there is nothing for the shared count to prefer - on a hex
   *  graph the node diagonal to two members is not adjacent to the seed, so every candidate shares
   *  exactly one face and the distance decides alone. In a stretched cell the neighbour across the
   *  thin direction is nearer than any neighbour along the layer by whatever the aspect ratio happens
   *  to be, so the CV steps that way, and then finds the next step in the SAME direction nearer still.
   *  It walks the boundary layer end to end: measured on a real mesh, 63% of the CVs built here came
   *  out as eight nodes in a straight line. That is the wrong shape, and worse, the wrong direction to
   *  coarsen in, since it merges exactly the wall-normal cells the implicit lines exist to keep apart.
   *
   *  The yardstick is the seed's own incident edges: a candidate's offset is divided by the length of
   *  the edge pointing most nearly the same way. A step across the layer and a step along it then both
   *  come to about one, whatever the stretching, and the shared count takes over from there. Choosing
   *  the edge by direction rather than by which Cartesian axis it is closest to is what keeps this
   *  usable on a curved boundary, where the wall-normal direction is not any one axis and a per-axis
   *  spacing would be measuring the wrong thing over most of the surface.
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

  /*--- A local frame at the current seed, up to nDim of its incident edges chosen to be as mutually
   *    orthogonal as possible, each with its own length. It is the yardstick described above: an
   *    offset is resolved onto these directions and each component divided by that direction's own
   *    spacing. Being built from the edges themselves it turns with the mesh, so it still measures
   *    the wall-normal direction correctly where a boundary curves away from any Cartesian axis, and
   *    on an axis-aligned mesh it reduces to dividing x, y and z by their own spacings. ---*/
  vector<std::array<su2double, MAXNDIM>> frameDir, edgeDir;
  vector<su2double> frameLen, edgeLen;
  vector<char> edgeUsed;

  const auto idx_after_bound_DBG = Index_CoarseCV;
  unsigned long nRejectedSeed_DBG = 0;

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

      edgeDir.clear();
      edgeLen.clear();
      su2double seedShortest = std::numeric_limits<su2double>::max();
      for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) {
        su2double e[MAXNDIM] = {0.0};
        GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jPoint), fine_grid->nodes->GetCoord(iPoint), e);
        const su2double len = GeometryToolbox::Norm(nDim, e);
        if (len <= 0.0) continue;
        std::array<su2double, MAXNDIM> u{};
        for (unsigned short d = 0; d < nDim; ++d) u[d] = e[d] / len;
        edgeDir.push_back(u);
        edgeLen.push_back(len);
        seedShortest = std::min(seedShortest, len);
      }
      /*--- A seed with no usable edge cannot be measured against anything, leave the scale at one. ---*/
      if (seedShortest == std::numeric_limits<su2double>::max()) seedShortest = 1.0;

      /*--- The shortest edge goes in first, so that the direction the mesh is stretched along is the
       *    one measured against its own small spacing; the rest are taken in order of how orthogonal
       *    they are to what is already in the frame. ---*/
      frameDir.clear();
      frameLen.clear();
      edgeUsed.assign(edgeDir.size(), 0);
      while (frameDir.size() < nDim) {
        long pick = -1;
        su2double bestScore = -1.0;
        for (size_t i = 0; i < edgeDir.size(); ++i) {
          if (edgeUsed[i]) continue;
          su2double score;
          if (frameDir.empty()) {
            score = 1.0 / edgeLen[i];
          } else {
            score = 2.0;
            for (const auto& f : frameDir) {
              su2double a = 0.0;
              for (unsigned short d = 0; d < nDim; ++d) a += f[d] * edgeDir[i][d];
              score = std::min(score, 1.0 - fabs(a));
            }
          }
          if (score > bestScore) {
            bestScore = score;
            pick = static_cast<long>(i);
          }
        }
        if (pick < 0) break;
        edgeUsed[pick] = 1;
        frameDir.push_back(edgeDir[pick]);
        frameLen.push_back(edgeLen[pick]);
      }
      /*--- Without a full frame the offset cannot be resolved, fall back to one isotropic scale. ---*/
      const bool haveFrame = (frameDir.size() == nDim);

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

          /*--- Distance to the centroid in cells rather than in metres: the offset is scaled by the
           *    seed edge pointing most nearly along it. ---*/
          su2double off[MAXNDIM] = {0.0};
          for (unsigned short d = 0; d < nDim; ++d) off[d] = fine_grid->nodes->GetCoord(CVPoint)[d] - centroid[d];
          const su2double offLen = GeometryToolbox::Norm(nDim, off);

          su2double dist = 0.0;
          if (haveFrame) {
            for (size_t k = 0; k < frameDir.size(); ++k) {
              su2double p = 0.0;
              for (unsigned short d = 0; d < nDim; ++d) p += off[d] * frameDir[k][d];
              const su2double q = p / frameLen[k];
              dist += q * q;
            }
          } else {
            const su2double r = offLen / seedShortest;
            dist = r * r;
          }

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

      nRejectedSeed_DBG++;
      MGQueue_InnerCV.MoveCV(iPoint, -1);
    }
  }

  const auto idx_after_domain_DBG = Index_CoarseCV;

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

  /*--- TEMPORARY DIAGNOSTIC: size distribution of the coarse CVs, by the phase that created them.
   *
   *    The child count on its own does not say what a CV looks like: eight children can be the 2x2x2
   *    block that is wanted, or a 4x2x1 slab, or a 1x8 strip, and in a planar slice through the mesh
   *    those are indistinguishable from a CV that really is small - a cube shows only its four nodes
   *    that lie in the slice, a slab edge-on shows two. Counting the fine edges whose two ends fall in
   *    the same CV separates them without reference to any coordinate direction: a 2x2x2 block has 12
   *    internal edges, a flat 2x2 has 4, a 1x4 strip has 3, a pair has 1. So "8 children, 12 edges" is
   *    a cube and "8 children, 10 edges" is not. ---*/
  vector<unsigned short> intEdges_DBG(Index_CoarseCV, 0);
  for (auto iPoint = 0ul; iPoint < fine_grid->GetnPoint(); iPoint++) {
    if (!fine_grid->nodes->GetDomain(iPoint)) continue;
    const auto pi = fine_grid->nodes->GetParent_CV(iPoint);
    if (pi >= Index_CoarseCV) continue;
    for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) {
      if (jPoint <= iPoint) continue;
      if (!fine_grid->nodes->GetDomain(jPoint)) continue;
      if (fine_grid->nodes->GetParent_CV(jPoint) == pi) intEdges_DBG[pi]++;
    }
  }

  /*--- Every rank owns a slice of the mesh, so every rank has its own share of these counts. Summing
   *    them and printing once is the only version that means anything: rank 0's slice is not the
   *    grid, and printing from all ranks interleaves the lines into each other - at eight ranks the
   *    output was literally unreadable, with one rank's histogram spliced through the middle of
   *    another's. Values are packed into one buffer so the whole report costs two collectives. ---*/
  {
    constexpr unsigned NPHASE = 5, NPV = 30; /*--- Phases, and values accumulated per phase. ---*/
    /*--- Layout within a phase: [0..9] children histogram, [10] CVs, [11] nodes, [12] cubes,
     *    [13] slabs, [14] squares, [15] strips, [16..29] internal-edge histogram of the 8s. ---*/
    vector<unsigned long> acc(NPHASE * NPV + 5, 0);

    auto histOf = [&](unsigned long lo, unsigned long hi, unsigned phase) {
      auto* a = &acc[phase * NPV];
      for (auto c = lo; c < hi; ++c) {
        const auto n = nodes->GetnChildren_CV(c);
        const auto e = intEdges_DBG[c];
        if (n == 8) {
          a[(e >= 12) ? 12 : 13]++;
          a[16 + std::min<unsigned short>(e, 13)]++;
        }
        if (n == 4) a[(e >= 4) ? 14 : 15]++;
        a[std::min<unsigned short>(n, 9)]++;
        a[10]++;
        a[11] += n;
      }
    };
    histOf(0, starting_idx_lines_DBG, 0);
    histOf(starting_idx_lines_DBG, idx_after_lines_DBG, 1);
    histOf(idx_after_lines_DBG, idx_after_bound_DBG, 2);
    histOf(idx_after_bound_DBG, idx_after_domain_DBG, 3);
    histOf(idx_after_domain_DBG, Index_CoarseCV, 4);

    for (auto iPoint = 0ul; iPoint < fine_grid->GetnPoint(); iPoint++) {
      if (!fine_grid->nodes->GetDomain(iPoint)) continue;
      acc[NPHASE * NPV + 0] += fine_grid->nodes->GetnPoint(iPoint);
      acc[NPHASE * NPV + 1]++;
    }
    acc[NPHASE * NPV + 2] = nRejectedSeed_DBG;
    acc[NPHASE * NPV + 3] = iteration;
    acc[NPHASE * NPV + 4] = fine_grid->GetnPointDomain();

    vector<unsigned long> tot(acc.size(), 0);
    SU2_MPI::Allreduce(acc.data(), tot.data(), static_cast<int>(acc.size()), MPI_UNSIGNED_LONG, MPI_SUM,
                       SU2_MPI::GetComm());

    /*--- Thickening is counted per configuration-file marker, not per local marker: ranks agree on
     *    neither the number nor the order of local markers, because each appends its own
     *    SEND_RECEIVE markers, so the same index means a different boundary elsewhere. ---*/
    const auto nMarkerCfg = config->GetnMarker_CfgFile();
    vector<unsigned long> mk(2 * nMarkerCfg, 0), mkTot(2 * nMarkerCfg, 0);
    for (auto m = 0u; m < fine_grid->GetnMarker(); ++m) {
      if (config->GetMarker_All_KindBC(m) == SEND_RECEIVE) continue;
      const auto c = config->GetMarker_CfgFile_TagBound(config->GetMarker_All_TagBound(m));
      mk[2 * c] += nThicken_DBG[m];
      mk[2 * c + 1] += nFlat_DBG[m];
    }
    if (nMarkerCfg > 0)
      SU2_MPI::Allreduce(mk.data(), mkTot.data(), static_cast<int>(mk.size()), MPI_UNSIGNED_LONG, MPI_SUM,
                         SU2_MPI::GetComm());

    if (rank == MASTER_NODE) {
      const char* phaseName[NPHASE] = {"pre-lines      ", "implicit lines ", "boundary STEP1 ", "domain STEP2   ",
                                       "leftover single"};
      const auto degN = std::max(tot[NPHASE * NPV + 1], 1ul);
      cout << "  CV size distribution by phase (maxAgglomSize=" << maxAgglomSize
           << ", mean fine-graph degree=" << (su2double(tot[NPHASE * NPV + 0]) / su2double(degN))
           << ", a 2x2x2 block has 12 internal edges):" << endl;

      for (unsigned ph = 0; ph < NPHASE; ++ph) {
        const auto* a = &tot[ph * NPV];
        if (a[10] == 0) continue;
        cout << "    " << phaseName[ph] << ": " << a[10] << " CVs, " << a[11] << " nodes, avg "
             << (su2double(a[11]) / su2double(a[10])) << "  sizes";
        for (unsigned sz = 1; sz <= 9; ++sz)
          if (a[sz] > 0) cout << " " << sz << ":" << a[sz];
        if (a[8] > 0) {
          cout << "  [of the 8s: " << a[12] << " are 2x2x2 cubes, " << a[13] << " are slabs/strips; internal edges";
          for (unsigned e = 7; e <= 13; ++e)
            if (a[16 + e] > 0) cout << " " << (e == 13 ? ">=13" : to_string(e)) << ":" << a[16 + e];
          cout << "]";
        }
        if (a[4] > 0) cout << "  [of the 4s: " << a[14] << " square, " << a[15] << " strip]";
        cout << endl;
      }
      cout << "    STEP2 rejected seeds: " << tot[NPHASE * NPV + 2] << ", iterations used: " << tot[NPHASE * NPV + 3]
           << "/" << tot[NPHASE * NPV + 4] << endl;
      for (auto c = 0u; c < nMarkerCfg; ++c) {
        if (mkTot[2 * c] + mkTot[2 * c + 1] == 0) continue;
        cout << "    marker " << config->GetMarker_CfgFile_TagBound(c) << ": thickened " << mkTot[2 * c]
             << ", kept flat " << mkTot[2 * c + 1] << endl;
      }
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

  /*--- Coarse CVs that must be left exactly as the agglomeration made them: those holding a node
   *    where two different boundary conditions meet. Both repair passes below exist to get rid of
   *    one-child control volumes, and a deliberately isolated junction IS a one-child control volume,
   *    so without this they simply undo the isolation. They have to be protected as a TARGET as well
   *    as a source: pass two merges a singleton into its smallest neighbour, and a one-child junction
   *    CV is by construction the smallest neighbour there is. ---*/
  vector<bool> mustStayAlone(nPointDomain, false);
  /*--- ...and which coarse CVs hold a boundary node at all. A boundary node is never agglomerated
   *    with an interior one, which the agglomeration itself already respects - the first advance of
   *    every front is a single layer, so the boundary layer is its own CV. The repair passes below
   *    would otherwise undo it from the other end: they merge a one-child CV into a neighbour, and a
   *    boundary CV's neighbours include the interior CV sitting on top of it. On the flat plate they
   *    happened to pick a boundary neighbour every time, which is luck rather than a rule. ---*/
  vector<bool> cvOnBoundary(nPointDomain, false);
  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++)
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
      const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      if (mixedBC[iFinePoint]) mustStayAlone[iCoarsePoint] = true;
      if (onPhysBoundary[iFinePoint]) cvOnBoundary[iCoarsePoint] = true;
    }

  /*--- Which physical boundaries each coarse CV sits on, one bit per marker. cvOnBoundary above only
   *    records THAT a CV touches a boundary, and the repair passes below compare nothing else, so a
   *    one-child CV on boundary A is free to be merged into a neighbour that lies on boundary B. That
   *    hands the target CV a marker none of its own children carried: SetVertex gives a coarse CV
   *    every marker of every child, so the merged CV becomes a vertex of A while its centroid sits
   *    wherever the B stack put it, and the boundary condition for A is then applied to a control
   *    volume that is mostly not on A. Comparing the marker SETS - a strictly stronger test than the
   *    boolean, since an interior CV has an empty set - keeps a merge inside one boundary. ---*/
  vector<unsigned long long> cvMarkerMask(nPointDomain, 0);
  {
    vector<int> bitOfMarker(fine_grid->GetnMarker(), -1);
    int nBits = 0;
    for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
      if (nBits < 64) bitOfMarker[iMarker] = nBits;
      nBits++;
    }
    if (nBits <= 64) {
      for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++)
        for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
          const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
          for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
            if (bitOfMarker[iMarker] < 0) continue;
            if (fine_grid->nodes->GetVertex(iFinePoint, iMarker) >= 0)
              cvMarkerMask[iCoarsePoint] |= 1ULL << bitOfMarker[iMarker];
          }
        }
    } else {
      /*--- More markers than bits: fall back to the boolean test. ---*/
      for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++)
        cvMarkerMask[iCoarsePoint] = cvOnBoundary[iCoarsePoint] ? 1ULL : 0ULL;
    }
  }

  /*--- A boundary CV built by the paving is the BASE of a stack, and the paving's whole guarantee is
   *    that every layer above it was given that same footprint. The repair passes must not take the
   *    footprint away: merging a stack base into the neighbouring stack leaves the column above it
   *    headless and the merged base wider than either column, so base and first layer no longer line
   *    up. On the next level that misalignment is fatal rather than cosmetic - the stiffest neighbour
   *    of the widened base is then a LATERAL one, SeedFrontNodes' hasLayerNormalTo rejects it, the CV
   *    does not seed a front, and it is swallowed mid-stack by an interior front instead. The coarse
   *    CV that results carries the boundary marker with its body off the boundary, and the boundary
   *    condition is applied to it. A one-node patch marching as a one-wide stack is a deliberate
   *    choice in AgglomerateImplicitLines, not damage for these passes to repair; what they are for
   *    is the INTERIOR singleton left where a line narrows, and that is untouched by this. ---*/
  auto isStackBase = [&](unsigned long iCoarsePoint) {
    return cvOnBoundary[iCoarsePoint] && (iCoarsePoint >= starting_idx_lines_DBG) &&
           (iCoarsePoint < idx_after_lines_DBG);
  };

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
    if (mustStayAlone[iCoarsePoint]) continue;
    if ((nodes->GetnPoint(iCoarsePoint) == 1) && !touchesPartition[iCoarsePoint]) {
      /*--- Find the neighbor of the isolated point. This neighbor is the right control volume ---*/

      const auto iCoarsePoint_Complete = nodes->GetPoint(iCoarsePoint, 0);
      if (mustStayAlone[iCoarsePoint_Complete]) continue;
      if (isStackBase(iCoarsePoint) || isStackBase(iCoarsePoint_Complete)) continue;
      if (cvMarkerMask[iCoarsePoint] != cvMarkerMask[iCoarsePoint_Complete]) continue;

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
          if (mustStayAlone[jCoarsePoint]) continue;
          if (isStackBase(jCoarsePoint)) continue;
          if (cvMarkerMask[jCoarsePoint] != cvMarkerMask[iCoarsePoint_Complete]) continue;

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
    if (mustStayAlone[iCoarsePoint]) continue;
    if (isStackBase(iCoarsePoint)) continue;
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
      if (mustStayAlone[jCoarsePoint]) continue;
      if (isStackBase(jCoarsePoint)) continue;
      if (cvMarkerMask[jCoarsePoint] != cvMarkerMask[iCoarsePoint]) continue;
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

  /*--- Optional dump of the finished agglomeration: one line per owned fine point, "x y z parentCV".
   *    It has to come from the END of the constructor, after both repair passes and the renumbering,
   *    because those still move fine points between coarse CVs - a dump taken before them shows what
   *    the agglomeration intended rather than what the solver will actually use. ---*/
  if (getenv("DUMPAGGLOM") != nullptr) {
    ofstream fdump(string("agglom_level") + to_string(iMesh) + "_r" + to_string(rank) + ".dat");
    for (auto iPoint = 0ul; iPoint < fine_grid->GetnPoint(); iPoint++) {
      if (!fine_grid->nodes->GetDomain(iPoint)) continue;
      const auto* c = fine_grid->nodes->GetCoord(iPoint);
      fdump << c[0] << " " << c[1] << " " << (nDim == 3 ? c[2] : 0.0) << " " << fine_grid->nodes->GetParent_CV(iPoint)
            << "\n";
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

vector<char> CMultiGridGeometry::FindMixedBoundaryNodes(const CGeometry* fine_grid, const CConfig* config) const {
  vector<char> mixed(fine_grid->GetnPoint(), 0);

  /*--- The first physical condition seen at each node, -1 until one is. A second one of a different
   *    KIND is what makes the node mixed. A second marker of the same kind - two wall patches meeting
   *    - is not, and neither is SEND_RECEIVE, which records where the partition runs and says nothing
   *    about the boundary condition. ---*/
  vector<short> firstBC(fine_grid->GetnPoint(), -1);

  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    const auto bc = static_cast<short>(config->GetMarker_All_KindBC(iMarker));
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if (firstBC[iPoint] < 0)
        firstBC[iPoint] = bc;
      else if (firstBC[iPoint] != bc)
        mixed[iPoint] = 1;
    }
  }
  return mixed;
}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, vector<short> marker_seed,
                                               const CGeometry* fine_grid, const CConfig* config,
                                               const vector<char>& mixedBC) const {
  /*--- A node where two boundary conditions of different type meet is never merged with anything. ---*/
  if (mixedBC[CVPoint]) return false;

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

namespace {

/*--- Unit normal of a boundary at a vertex, false if the marker does not reach iPoint. Boundary
 *    normals point INTO the domain, so this doubles as a starting direction to march in. ---*/
bool VertexUnitNormal(const CGeometry* grid, unsigned short nDim, unsigned long iPoint, unsigned short iMarker,
                      su2double* unitNormal) {
  const long iVertex = grid->nodes->GetVertex(iPoint, iMarker);
  if (iVertex == -1) return false;
  grid->vertex[iMarker][iVertex]->GetNormal(unitNormal);
  const su2double nrm = GeometryToolbox::Norm(nDim, unitNormal);
  if (nrm <= 0.0) return false;
  for (unsigned short d = 0; d < nDim; ++d) unitNormal[d] /= nrm;
  return true;
}

}  // namespace

CMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeometry* fine_grid, const CConfig* config,
                                                                   const CNodeStiffness& stiff) const {
  /*--- Fraction of a marker's nodes that must sit in a layer before the whole marker may seed. ---*/
  constexpr su2double QUALIFIED_FRACTION = 0.5;
  constexpr su2double ANGLE_THRESHOLD_DEG = 30.0;
  const su2double cos_threshold = cos(ANGLE_THRESHOLD_DEG * PI_NUMBER / 180.0);

  const su2double MIN_AR = config->GetMGOptions().MG_Implicit_Lines_Min_AR;
  const bool USE_AR = (MIN_AR > 1.0);

  const auto nMarkerFine = fine_grid->GetnMarker();
  constexpr auto NO_POINT = std::numeric_limits<unsigned long>::max();

  CFrontSeeds seeds;
  vector<char> taken(fine_grid->GetnPoint(), 0);

  auto isWall = [&](unsigned short bc) {
    return (bc == HEAT_FLUX) || (bc == ISOTHERMAL) || (bc == CHT_WALL_INTERFACE) || (bc == SMOLUCHOWSKI_MAXWELL);
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

  auto seedMarker = [&](unsigned short iMarker, bool requireLayer) {
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if (!fine_grid->nodes->GetDomain(iPoint)) continue;
      if (fine_grid->nodes->GetAgglomerate(iPoint)) continue;
      if (taken[iPoint]) continue; /*--- A node on two markers must seed only one front. ---*/

      su2double Normal[MAXNDIM] = {0.0};
      if (!VertexUnitNormal(fine_grid, nDim, iPoint, iMarker, Normal)) continue;
      if (requireLayer && !hasLayerNormalTo(iPoint, Normal)) continue;

      std::array<su2double, MAXNDIM> n0{};
      for (unsigned short d = 0; d < nDim; ++d) n0[d] = Normal[d];
      seeds.node.push_back(iPoint);
      seeds.normal.push_back(n0);
      taken[iPoint] = 1;
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
   *    that mostly does not qualify scatters one-node patches, i.e. coarse CVs that do not coarsen
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
       *    transform and have their own matching, which a front running into one would disturb. ---*/
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
        if (!VertexUnitNormal(fine_grid, nDim, iPoint, iMarker, Normal)) continue;
        nValid[cfgOfMarker[iMarker]]++;
        if (hasLayerNormalTo(iPoint, Normal)) nQualified[cfgOfMarker[iMarker]]++;
      }
    }

    /*--- A marker is generally split over several ranks, so the verdict must be taken on all of it.
     *    Every rank reaches these collectives, including one that owns no boundary at all. ---*/
    if (nMarkerCfg > 0) {
      vector<unsigned long> tmp(nMarkerCfg);
      SU2_MPI::Allreduce(nValid.data(), tmp.data(), nMarkerCfg, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      nValid.swap(tmp);
      SU2_MPI::Allreduce(nQualified.data(), tmp.data(), nMarkerCfg, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      nQualified.swap(tmp);
    }

    if ((rank == MASTER_NODE) && (getenv("DUMPAGGLOM") != nullptr))
      for (auto iCfg = 0u; iCfg < nMarkerCfg; ++iCfg)
        if (nValid[iCfg] > 0)
          cout << "  Seed vote " << config->GetMarker_CfgFile_TagBound(iCfg) << ": " << nQualified[iCfg] << "/"
               << nValid[iCfg] << " = " << (su2double(nQualified[iCfg]) / su2double(nValid[iCfg])) << endl;

    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
      if (!canSeed(iMarker)) continue;
      const auto iCfg = cfgOfMarker[iMarker];
      if (nValid[iCfg] == 0) continue;
      if (su2double(nQualified[iCfg]) < QUALIFIED_FRACTION * su2double(nValid[iCfg])) continue;
      seedMarker(iMarker, true);
    }
  }

  return seeds;
}

vector<vector<unsigned long>> CMultiGridGeometry::BuildFrontPatches(const CFrontSeeds& seeds,
                                                                    const CGeometry* fine_grid, const CConfig* config,
                                                                    const vector<char>& mixedBC) const {
  /*--- Every seed must end up in exactly one patch, and a patch must be compact: in 3D the four nodes
   *    of one boundary quadrilateral, in 2D the two ends of one boundary edge. Choosing, for each seed
   *    independently, a set of neighbours to merge with does not do this - the relation is not
   *    symmetric, so seed 1 claiming {2,3} does not stop seed 2 claiming {1,4}, and the patches
   *    overlap and fight over nodes.
   *
   *    Repeated pairwise matching avoids that by construction. One round pairs adjacent seeds into the
   *    boundary edge, a second pairs adjacent pairs into the boundary quadrilateral. Each round is a
   *    matching, so membership stays mutually exclusive and the result is a true partition. It needs
   *    nothing but point-to-point connectivity, so it works identically on every multigrid level -
   *    boundary face connectivity does not survive agglomeration, so a literal "same quadrilateral"
   *    test would only ever work for the first coarsening. Two rounds reach 4, which is the 3D patch
   *    size, so the number of rounds follows from max_group instead of iterating to a fixed point.
   *
   *    This patch is the ONLY thing that decides the footprint of the stack above it. The front that
   *    rises from it keeps exactly these nodes' successors, layer after layer, so getting the patch
   *    square is what makes the coarse CVs square all the way up. ---*/
  const auto nSeeds = seeds.node.size();
  unsigned long max_group = config->GetMGOptions().MG_Implicit_Lines_Max_Group;
  if (max_group == 0) max_group = (nDim == 2) ? 2 : 4;

  /*--- Marker signature of each seed, as a bitmask over the physical markers. Seeds may only be
   *    matched when these agree, so a patch never straddles a change of boundary condition - the rule
   *    ordinary agglomeration uses for ridges and valleys. ---*/
  const auto nMarkerFine = fine_grid->GetnMarker();
  vector<short> physBit(nMarkerFine, -1);
  unsigned nPhys = 0;
  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) physBit[iMarker] = static_cast<short>(nPhys++);

  const unsigned nWords = std::max(1u, (nPhys + 63u) / 64u);
  vector<uint64_t> sig(nSeeds * nWords, 0);
  for (unsigned long si = 0; si < nSeeds; ++si)
    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
      if (physBit[iMarker] < 0) continue;
      if (fine_grid->nodes->GetVertex(seeds.node[si], iMarker) == -1) continue;
      const auto b = static_cast<unsigned>(physBit[iMarker]);
      sig[si * nWords + b / 64] |= (uint64_t(1) << (b % 64));
    }
  auto sameSig = [&](unsigned long sa, unsigned long sb) {
    for (unsigned w = 0; w < nWords; ++w)
      if (sig[sa * nWords + w] != sig[sb * nWords + w]) return false;
    return true;
  };

  /*--- Seed-to-seed adjacency, inherited from the boundary nodes' mesh connectivity. ---*/
  vector<long> seedOfNode(fine_grid->GetnPoint(), -1);
  for (unsigned long si = 0; si < nSeeds; ++si) seedOfNode[seeds.node[si]] = static_cast<long>(si);

  vector<vector<unsigned long>> adj(nSeeds);
  for (unsigned long si = 0; si < nSeeds; ++si)
    for (auto jPoint : fine_grid->nodes->GetPoints(seeds.node[si])) {
      const auto sj = seedOfNode[jPoint];
      if ((sj >= 0) && (static_cast<unsigned long>(sj) != si)) adj[si].push_back(static_cast<unsigned long>(sj));
    }

  /*--- Global point index of each seed, used everywhere below as the deterministic sort key. Local
   *    indices depend on the partitioning, so ordering by them would make the coarse grid depend on
   *    the rank count. ---*/
  vector<unsigned long> sgkey(nSeeds);
  for (unsigned long si = 0; si < nSeeds; ++si) sgkey[si] = fine_grid->nodes->GetGlobalIndex(seeds.node[si]);

  vector<vector<unsigned long>> groups;
  vector<unsigned long> groupOf(nSeeds);
  groups.reserve(nSeeds);

  /*--- Every seed starts as its own group; the rounds below merge them. ---*/
  for (unsigned long si = 0; si < nSeeds; ++si) {
    groupOf[si] = si;
    groups.push_back({si});
  }

  const unsigned nRounds = (max_group <= 1) ? 0 : ((max_group <= 2) ? 1 : 2);

  /*--- One admissible merge of two groups, weighted by how many seed-to-seed adjacencies they share.
   *    A group lying ALONGSIDE this one touches it along its whole length and shares two, whereas one
   *    continuing in the same direction touches at an end and shares one. So weight 2 is the square
   *    and weight 1 is the strip, and the strip extrudes into a coarse CV elongated in one
   *    boundary-tangential direction. ---*/
  struct CMerge {
    unsigned long g, h;       /*!< \brief The two groups, g < h. */
    unsigned long weight;     /*!< \brief Shared adjacencies: 2 makes a square, 1 makes a strip. */
    unsigned long keyG, keyH; /*!< \brief Their global-index keys, the deterministic tie-break. */
  };

  vector<CMerge> merges;
  vector<std::pair<unsigned long, unsigned long>> shared;
  vector<unsigned long> gkey;
  vector<char> consumed;

  for (unsigned round = 0; round < nRounds; ++round) {
    const auto nGroups = groups.size();

    /*--- Sort key of each group: the smallest global point index it holds. ---*/
    gkey.assign(nGroups, std::numeric_limits<unsigned long>::max());
    for (unsigned long g = 0; g < nGroups; ++g)
      for (auto si : groups[g]) gkey[g] = std::min(gkey[g], sgkey[si]);

    /*--- Every merge this round could make. Counted once per unordered pair: adjacency is symmetric,
     *    so the count from g's side equals the count from h's, and taking only h > g avoids both. ---*/
    merges.clear();
    for (unsigned long g = 0; g < nGroups; ++g) {
      /*--- A node where two different boundary conditions meet stays a patch of its own, so the front
       *    rising from it is one node wide. Both sides of a merge are tested: skipping only the mixed
       *    group would still let an ordinary group reach out and take it. ---*/
      if (mixedBC[seeds.node[groups[g].front()]]) continue;
      shared.clear();
      for (auto si : groups[g])
        for (auto sj : adj[si]) {
          const auto h = groupOf[sj];
          if (h <= g) continue;
          if (mixedBC[seeds.node[groups[h].front()]]) continue;
          if (groups[g].size() + groups[h].size() > max_group) continue;
          if (!sameSig(groups[h].front(), groups[g].front())) continue;

          bool seen = false;
          for (auto& t : shared)
            if (t.first == h) {
              t.second++;
              seen = true;
              break;
            }
          if (!seen) shared.emplace_back(h, 1);
        }
      for (const auto& t : shared) merges.push_back({g, t.first, t.second, gkey[g], gkey[t.first]});
    }

    /*--- Best first, over ALL groups at once. Sweeping the groups in index order instead and letting
     *    each take its own best partner is what produced the strips: a group with no square partner
     *    left would take a weight-1 merge and consume a group that a later one needed for its square,
     *    and the failures cascade - on a structured 3D wall that came to 19% of the footprints.
     *    Ordering the merges globally means every square in the mesh is made before the first strip is
     *    even considered, so a strip only forms where the surface genuinely offers nothing better. ---*/
    std::sort(merges.begin(), merges.end(), [](const CMerge& a, const CMerge& b) {
      if (a.weight != b.weight) return a.weight > b.weight;
      if (a.keyG != b.keyG) return a.keyG < b.keyG;
      return a.keyH < b.keyH;
    });

    consumed.assign(nGroups, 0);
    vector<vector<unsigned long>> merged;
    merged.reserve(nGroups);

    for (const auto& m : merges) {
      if (consumed[m.g] || consumed[m.h]) continue;
      consumed[m.g] = consumed[m.h] = 1;
      auto group = groups[m.g];
      group.insert(group.end(), groups[m.h].begin(), groups[m.h].end());
      merged.push_back(std::move(group));
    }
    /*--- Whatever found no partner passes through unchanged. ---*/
    for (unsigned long g = 0; g < nGroups; ++g)
      if (!consumed[g]) merged.push_back(std::move(groups[g]));

    groups = std::move(merged);
    for (unsigned long g = 0; g < groups.size(); ++g)
      for (auto si : groups[g]) groupOf[si] = g;
  }

  return groups;
}

void CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseCV, const CGeometry* fine_grid,
                                                  const CConfig* config) {
  /*==================================================================================================
   *  Paving by advancing fronts.
   *
   *  The boundary is agglomerated first into patches (PHASE 1), and every patch then rises into the
   *  domain as a FRONT that keeps its footprint: at each round every front node picks a successor,
   *  and the front advances only if all of them succeeded and none was lost to another front. So the
   *  layers of one stack are congruent by construction, and a coarse CV can never contain a node that
   *  belongs over a neighbouring patch.
   *
   *  This is the part the earlier line-based version could not guarantee, and the reason it is gone.
   *  There, each wall node grew its own line first and the lines were grouped into bundles only
   *  afterwards, so nothing tied a line to the footprint it would later be asked to share: the
   *  marching direction was re-set to the last step taken every step, which let a line random-walk
   *  tangentially one legal 30-degree step at a time into a neighbouring column, and the connectivity
   *  test that should have caught it was applied to the lines' WALL roots, which stay adjacent no
   *  matter how far apart their tops drift. Ragged line lengths then made it worse, because a bundle
   *  was allowed to carry on with whichever subset of its lines was still long enough, so the stack
   *  changed footprint as it rose.
   *
   *  Three rules replace all of that:
   *
   *   - The front is the primitive. Nothing marches except a whole patch, so there is no such thing as
   *     an individual line to drift.
   *   - All-or-nothing layers. A front that cannot fill an entire layer retires and leaves the rest to
   *     ordinary agglomeration, instead of continuing narrower.
   *   - Contention resolved globally, once per layer, from bids collected before any is granted. Two
   *     fronts reaching for the same node is exactly the event "the fronts have met", and it is caught
   *     in the layer where it happens rather than fifteen layers later.
   *
   *  A front is stopped by TWO things and nothing else: reaching a boundary, or being unable to lay a
   *  layer topologically identical to the one it is standing on (see layerIsIsomorphic below). There
   *  is no limit on how far it may turn and no threshold on how stretched the mesh has to be. Those
   *  limits used to exist and they were the wrong instrument: an aspect-ratio cut in particular stops
   *  each front on a contour of the LOCAL cell shape, which on a flat plate is a contour of the
   *  streamwise spacing, so fronts died at different heights and the paved region ended in a
   *  staircase - and, being driven by dx rather than by the boundary layer, a staircase running the
   *  wrong way. Without it the same case paves to the far boundary at a uniform depth and the domain
   *  pass has nothing left to do.
   *
   *  Each coarse CV is the front's footprint taken TWO layers deep, so the paving coarsens by the same
   *  factor in every direction: a 2x2 wall patch and the two layers above it make the 2x2x2 block, and
   *  the next CV of the stack starts from the footprint the front already has. Taking one layer per CV
   *  instead leaves the wall-normal direction uncoarsened entirely, which is what a line-implicit
   *  smoother wants and not what this is for. The single-layer CV survives in one place only: the top
   *  of a stack that retires with one layer buffered, where the alternative is dropping it back to
   *  ordinary agglomeration.
   *
   *  The multigrid queue is deliberately not touched here: the sync loop after the boundary
   *  agglomeration removes every point already marked agglomerated, so doing it here too would be an
   *  error.
   *================================================================================================*/
  const auto starting_Index_CoarseCV = Index_CoarseCV;
  const auto nPointFine = fine_grid->GetnPoint();
  const auto nMarkerFine = fine_grid->GetnMarker();
  constexpr auto NO_POINT = std::numeric_limits<unsigned long>::max();
  const short int maxAgglomSize = (nDim == 2) ? 4 : 8;

  /*--- How nearly parallel a step must be to a boundary's normal to count as running INTO that
   *    boundary rather than along it. This is not a limit on where a front may go - it is only how
   *    "the front has reached a boundary" is recognised. ---*/
  constexpr su2double BOUNDARY_ALIGN_DEG = 30.0;
  const su2double cos_boundary = cos(BOUNDARY_ALIGN_DEG * PI_NUMBER / 180.0);
  /*--- Weight of the new step direction when the front's marching direction is updated. The direction
   *    only ever RANKS candidates, it never rejects one, so this is a preference and not a limit. ---*/
  constexpr su2double DIR_BLEND = 0.5;

  /*--- Safety cap on stack depth, 0 for none. A front is meant to run until it reaches a boundary or
   *    the mesh stops offering a clean extrusion, so this is off by default. ---*/
  const unsigned long MAX_LINE_LENGTH = config->GetMGOptions().MG_Implicit_Lines_MaxLength;
  unsigned long max_group = config->GetMGOptions().MG_Implicit_Lines_Max_Group;
  if (max_group == 0) max_group = (nDim == 2) ? 2 : 4;

  const auto stiff = ComputeNodeStiffness(fine_grid);

  /*--- PHASE 1. SeedFrontNodes must be reached by every rank, including one that owns no boundary:
   *    it takes a collective to agree on which markers carry a layer. ---*/
  const auto seeds = SeedFrontNodes(fine_grid, config, stiff);
  const auto mixedBC = FindMixedBoundaryNodes(fine_grid, config);
  const auto patches = BuildFrontPatches(seeds, fine_grid, config, mixedBC);

  /*--- Nodes on a boundary that carries a boundary condition. A front must not grow into one: those
   *    nodes belong to the boundary agglomeration and a stack absorbing one would straddle two
   *    boundaries. CPoint's Boundary flag cannot answer this, as it is also set by SEND_RECEIVE, so on
   *    a partitioned mesh it is true for ordinary interior nodes of the send fringe and every front
   *    would stop one layer short of the partition. Walking the markers' own vertex lists is both
   *    exact and cheaper than testing every point against every marker. ---*/
  vector<char> onPhysicalBoundary(nPointFine, 0);
  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++)
      onPhysicalBoundary[fine_grid->vertex[iMarker][iVertex]->GetNode()] = 1;
  }

  /*--- A node on a physical boundary only ENDS a front if the front is stepping INTO that boundary,
   *    i.e. the step direction is roughly parallel to the boundary's own normal. A node that merely
   *    runs ALONG a boundary - a column on a spanwise symmetry plane, say - has that boundary's normal
   *    roughly PERPENDICULAR to the step, and is a legitimate interior node of the stack, not its end:
   *    onPhysicalBoundary alone cannot tell these apart, since it only records marker membership, not
   *    which direction the marker's surface runs in. Without this a front that happens to sit on a
   *    tangential boundary the whole way up dies at its very first step. ---*/
  auto entersBoundary = [&](unsigned long jPoint, const su2double* stepDir) {
    for (unsigned short iMarker = 0; iMarker < nMarkerFine; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
      su2double n[MAXNDIM] = {0.0};
      if (!VertexUnitNormal(fine_grid, nDim, jPoint, iMarker, n)) continue;
      if (fabs(GeometryToolbox::DotProduct(nDim, n, stepDir)) >= cos_boundary) return true;
    }
    return false;
  };

  /*==================================================================================================
   *  PHASE 2 - advance every front, one layer per round.
   *================================================================================================*/

  /*--- One front node's chosen successor, before contention over it is resolved. ---*/
  struct CStep {
    unsigned long node;     /*!< Candidate successor. */
    unsigned long from;     /*!< The front node that proposed it. */
    unsigned long key;      /*!< Global index of "from", the final deterministic tie-break. */
    su2double score;        /*!< Alignment of the step with the front's marching direction. */
    su2double dist;         /*!< Length of the step. */
    su2double dir[MAXNDIM]; /*!< Unit step direction, reused to update the front's direction. */
  };

  /*--- Fronts are not a fixed set: one handed over from a neighbouring rank is appended while the
   *    rounds are running, so every per-front array grows and the loops below are bounded by
   *    front.size() rather than by the number of patches. ---*/
  vector<vector<unsigned long>> front, pending, handTo;
  vector<std::array<su2double, MAXNDIM>> dirNow;
  vector<char> alive;
  vector<unsigned long> depth, nBlock, pendingLayers, tag;
  vector<vector<CStep>> prop;

  /*--- A name for a front that means the same thing on every rank, so a stack handed across a
   *    partition can be recognised on the far side and so two ranks reaching for the same node can be
   *    separated the same way by both. The smallest global point index of the patch it grew from is
   *    unique, since a seed belongs to exactly one patch; the +1 leaves 0 free to mean "nothing". ---*/
  auto addFront = [&](const vector<unsigned long>& layer, const std::array<su2double, MAXNDIM>& dir,
                      unsigned long frontTag, unsigned long block) {
    front.push_back(layer);
    pending.push_back(layer);
    handTo.emplace_back();
    dirNow.push_back(dir);
    alive.push_back(1);
    depth.push_back(0);
    nBlock.push_back(block);
    pendingLayers.push_back(1);
    tag.push_back(frontTag);
    prop.emplace_back();
    return front.size() - 1;
  };

  vector<char> claimed(nPointFine, 0);
  /*--- Confirmed owner of a claimed node, -1 while free. Only ever written when a layer is accepted,
   *    so a bid that is still being contested never appears here. ---*/
  vector<int> frontOf(nPointFine, -1);

  vector<char> failed;
  vector<unsigned short> failReason;
  vector<unsigned long> stopCounts(N_STOP_REASONS, 0);

  /*--- Set when a front hands only PART of its footprint over and goes on marching here with what is
   *    left of it, so the retirement pass at the end of the round knows not to kill it. ---*/
  vector<char> keepLocal;
  /*--- The name the handed-over piece travels under. After a split this is NOT the name of the front
   *    it came from: the two pieces are separate stacks from here on, and giving them one name would
   *    let the far side group a piece of this stack with a piece of another one. ---*/
  vector<unsigned long> handTag;

  /*--- A name for a set of nodes that both ranks sharing them would compute identically. The nodes of
   *    a footprint are claimed by one front and by no other, so the smallest global index in it is a
   *    unique name for that front; the +1 leaves 0 free to mean "nothing". ---*/
  auto tagOfSet = [&](const vector<unsigned long>& set) {
    unsigned long t = std::numeric_limits<unsigned long>::max();
    for (auto p : set) t = std::min(t, fine_grid->nodes->GetGlobalIndex(p));
    return t + 1;
  };

  /*--- The bid table. Only the index is kept per mesh point, and the bids themselves live in a
   *    compact vector holding one entry per candidate actually bid on this round - a few per front,
   *    against one entry per point in the mesh. Storing a whole CStep per point instead costs about
   *    sixty bytes times nPoint, which on the meshes this code is meant for is hundreds of megabytes
   *    of table that is empty almost everywhere. ---*/
  constexpr unsigned NOBID = std::numeric_limits<unsigned>::max();
  vector<unsigned> bidIdx(nPointFine, NOBID);
  vector<CStep> bids;
  vector<unsigned long> bidOwner;

  /*--- Scratch for the layer under construction, hoisted so a front does not allocate per layer. ---*/
  vector<unsigned long> newLayer;

  unsigned long nStacks = 0, nSemiCV = 0, nFullCV = 0, nLayers = 0, nCovered = 0;
  unsigned long nHandedOut = 0, nHandedIn = 0;
  unsigned long nSplit = 0, nSplitLocal = 0, nSplitHanded = 0, nSplitDropped = 0, nNoHalo = 0;

  /*--- One footprint node arriving from a neighbouring rank, to be regrouped by tag. ---*/
  struct CInherited {
    unsigned long tag;      /*!< \brief The front it belongs to, the same name on both ranks. */
    unsigned long node;     /*!< \brief Local index of the node, owned by this rank. */
    su2double dir[MAXNDIM]; /*!< \brief The marching direction the stack arrives with. */
  };
  vector<CInherited> inherited;

  /*--- Where each halo node sits in this rank's SEND_RECEIVE receive lists, so a handover can be
   *    packed against the right vertex without searching. A halo node appears in exactly one. ---*/
  vector<int> haloMarker(nPointFine, -1);
  vector<unsigned long> haloVertex(nPointFine, 0);
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (!((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)))
      continue;
    const auto MarkerR = iMarker + 1;
    for (auto iVertex = 0ul; iVertex < fine_grid->nVertex[MarkerR]; iVertex++) {
      const auto p = fine_grid->vertex[MarkerR][iVertex]->GetNode();
      haloMarker[p] = static_cast<int>(MarkerR);
      haloVertex[p] = iVertex;
    }
  }
  unsigned long histogram[9] = {0};

  auto markFail = [&](unsigned long f, unsigned short why) {
    if (!failed[f]) {
      failed[f] = 1;
      failReason[f] = why;
    }
  };

  /*--- How many fine layers the next coarse CV of this front holds: always two, so the stack coarsens
   *    by the same factor along the marching direction as the footprint does across it. The only
   *    exception is a footprint already so wide that a second layer would exceed the agglomeration
   *    size limit, which can only happen if MG_IMPLICIT_LINES_MAX_GROUP was raised past the
   *    dimension's default. ---*/
  auto blockFor = [&](const vector<unsigned long>& layer) -> unsigned long {
    return (layer.size() * 2 > static_cast<size_t>(maxAgglomSize)) ? 1 : 2;
  };

  /*--- Turn everything buffered for this front into one coarse control volume. ---*/
  auto emit = [&](unsigned long f) {
    if (pending[f].empty()) return;
    for (unsigned long c = 0; c < pending[f].size(); ++c) {
      const auto p = pending[f][c];
      fine_grid->nodes->SetParent_CV(p, Index_CoarseCV);
      nodes->SetChildren_CV(Index_CoarseCV, c, p);
      if (fine_grid->nodes->GetAgglomerate_Indirect(p)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);
    }
    nodes->SetnChildren_CV(Index_CoarseCV, static_cast<unsigned short>(pending[f].size()));
    Index_CoarseCV++;
    nCovered += pending[f].size();
    ((pendingLayers[f] == 1) ? nSemiCV : nFullCV)++;

    pending[f].clear();
    pendingLayers[f] = 0;
    nBlock[f] = blockFor(front[f]);
  };

  auto isAdjacent = [&](unsigned long a, unsigned long b) {
    const auto& pts = fine_grid->nodes->GetPoints(a);
    return std::find(pts.begin(), pts.end(), b) != pts.end();
  };

  /*--- Is a footprint one connected patch? A set that falls into pieces is not the extrusion of
   *    anything, so a footprint arriving from a neighbour and a piece left behind by a split both
   *    have to pass this before they are allowed to carry a stack. ---*/
  auto isConnectedLayer = [&](const vector<unsigned long>& layer) {
    if (layer.size() < 2) return true;
    vector<char> seen(layer.size(), 0);
    vector<size_t> stk{0};
    seen[0] = 1;
    size_t nSeen = 1;
    while (!stk.empty()) {
      const auto cur = stk.back();
      stk.pop_back();
      for (size_t k = 0; k < layer.size(); ++k) {
        if (seen[k] || !isAdjacent(layer[cur], layer[k])) continue;
        seen[k] = 1;
        nSeen++;
        stk.push_back(k);
      }
    }
    return nSeen == layer.size();
  };

  /*--- The one test that decides whether a front may advance: is the layer it is about to lay down
   *    TOPOLOGICALLY IDENTICAL to the layer it is standing on?
   *
   *    "old" and "new" are index-aligned, so phi maps old[k] to new[k], and phi is the extrusion the
   *    front is proposing. It is a valid layer exactly when phi is an isomorphism of the two induced
   *    subgraphs AND matches them up one for one:
   *
   *      - same number of cells, which index alignment already gives;
   *      - every new cell is adjacent to exactly one old cell, and that one is its own preimage;
   *      - every old cell is adjacent to exactly one new cell, and that one is its own image;
   *      - the same edges: old[k]-old[l] is an edge if and only if new[k]-new[l] is, which makes the
   *        edge counts equal and carries connectivity across from the old layer for free.
   *
   *    Nothing else stops a front. There is no cone on how far it may turn and no threshold on how
   *    stretched the mesh has to be: it runs until it reaches a boundary or until the mesh stops
   *    offering a clean extrusion, and this is what "stops offering" means. ---*/
  auto layerIsIsomorphic = [&](const vector<unsigned long>& oldL, const vector<unsigned long>& newL) {
    const auto n = oldL.size();
    if (newL.size() != n) return false;

    for (size_t k = 0; k < n; ++k) {
      unsigned nOld = 0, nNew = 0;
      for (size_t l = 0; l < n; ++l) {
        nOld += isAdjacent(newL[k], oldL[l]);
        nNew += isAdjacent(oldL[k], newL[l]);
      }
      /*--- Exactly one partner each way, and it has to be the one phi names. ---*/
      if ((nOld != 1) || (nNew != 1)) return false;
      if (!isAdjacent(oldL[k], newL[k])) return false;
    }

    for (size_t k = 0; k < n; ++k)
      for (size_t l = k + 1; l < n; ++l)
        if (isAdjacent(oldL[k], oldL[l]) != isAdjacent(newL[k], newL[l])) return false;

    return true;
  };

  /*--- The boundary layer of every front. Claiming it before ordinary boundary agglomeration runs is
   *    what keeps the stack aligned: every layer above has exactly the same footprint. ---*/
  for (const auto& patch : patches) {
    /*--- A one-node patch is allowed to seed a front and marches as a stack one node wide. It looks
     *    like a poor footprint, but the alternative - leaving it to ordinary agglomeration - is far
     *    worse: a seed that cannot pair is one whose whole COLUMN then goes unpaved, from the wall to
     *    wherever the front would have stopped, and those columns land in exactly the places that
     *    cannot pair for a reason. On the flat plate they were the two nodes at the leading edge,
     *    where the wall meets the symmetry plane and the marker signatures differ, plus the two
     *    domain corners: five full-height stripes cut through the paved region, the worst of them
     *    right at the leading edge. ---*/
    bool valid = !patch.empty();
    for (auto si : patch) {
      const auto p = seeds.node[si];
      if (!GeometricalCheck(p, fine_grid, config) || fine_grid->nodes->GetAgglomerate(p)) valid = false;
    }
    if (!valid) continue;

    vector<unsigned long> layer0;
    std::array<su2double, MAXNDIM> n0{};
    unsigned long frontTag = std::numeric_limits<unsigned long>::max();
    for (auto si : patch) {
      layer0.push_back(seeds.node[si]);
      frontTag = std::min(frontTag, fine_grid->nodes->GetGlobalIndex(seeds.node[si]));
      for (unsigned short d = 0; d < nDim; ++d) n0[d] += seeds.normal[si][d];
    }
    const su2double nrm = GeometryToolbox::Norm(nDim, n0.data());
    /*--- A patch whose members' normals cancel has no direction to march in. ---*/
    if (nrm <= 0.0) continue;
    for (unsigned short d = 0; d < nDim; ++d) n0[d] /= nrm;

    /*--- The boundary layer becomes a coarse CV on its own: a boundary node is never merged with an
     *    interior one, so the first advance of every front is a single layer. Only the first - emit()
     *    then asks blockFor again for what is by then an ordinary interior layer, and the rest of the
     *    stack rises two nodes at a time. This is also what isolates a junction node without needing
     *    a rule of its own: such a node is a patch of one, so its first CV holds it and nothing else. ---*/
    const auto f = addFront(layer0, n0, frontTag + 1, 1);
    for (auto p : layer0) {
      claimed[p] = 1;
      frontOf[p] = static_cast<int>(f);
    }
    histogram[std::min<size_t>(layer0.size(), 8)]++;
    nStacks++;
    nLayers++;
    emit(f);
  }

  for (unsigned long layer = 1;; ++layer) {
    /*--- Whether ANY rank still has a live front, not just this one. Every rank has to run the same
     *    number of rounds because each round ends in a handover exchange that they all take part in:
     *    a rank whose own fronts are long finished may still be about to receive a stack from a
     *    neighbour, and a rank that dropped out of the loop early would hang the ones that did not. ---*/
    int aliveLocal = 0;
    for (unsigned long f = 0; f < front.size(); ++f) aliveLocal |= alive[f];
    int aliveGlobal = 0;
    SU2_MPI::Allreduce(&aliveLocal, &aliveGlobal, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
    if (aliveGlobal == 0) break;

    failed.assign(front.size(), 0);
    failReason.assign(front.size(), 0);
    keepLocal.assign(front.size(), 0);
    handTag.assign(front.size(), 0);
    for (const auto& b : bids) bidIdx[b.node] = NOBID;
    bids.clear();
    bidOwner.clear();

    /*--- (a) Every alive front proposes a successor for each of its nodes. A front that cannot fill a
     *    whole layer proposes NOTHING: it is retiring this round anyway, and letting its partial bids
     *    stand would let a dying front displace a healthy one out of nodes it can still use. ---*/
    for (unsigned long f = 0; f < front.size(); ++f) {
      if (!alive[f]) continue;
      prop[f].clear();
      handTo[f].clear();

      if ((MAX_LINE_LENGTH > 0) && (depth[f] + 1 >= MAX_LINE_LENGTH)) {
        markFail(f, STOP_MAX_LENGTH);
        continue;
      }

      for (auto n : front[f]) {
        auto best = NO_POINT;
        su2double best_dot = -2.0, best_len = 0.0, best_dir[MAXNDIM] = {0.0};
        /*--- The best step onto a node this rank does NOT own, kept separately. It cannot be claimed
         *    here, but it is where the stack would go next, so it is what gets handed over. ---*/
        auto bestHalo = NO_POINT;
        su2double bestHalo_dot = -2.0;
        bool sawCollision = false, sawBoundary = false, sawAgglom = false, sawPartition = false;
        bool sawGeom = false;

        for (auto iNeigh = 0u; iNeigh < fine_grid->nodes->GetnPoint(n); ++iNeigh) {
          const auto jPoint = fine_grid->nodes->GetPoint(n, iNeigh);

          su2double vec[MAXNDIM] = {0.0};
          GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jPoint), fine_grid->nodes->GetCoord(n), vec);
          const su2double len = GeometryToolbox::Norm(nDim, vec);
          if (len <= 0.0) continue;
          for (unsigned short d = 0; d < nDim; ++d) vec[d] /= len;

          /*--- The marching direction RANKS the candidates and nothing more: whichever free neighbour
           *    lies most nearly ahead is the one proposed. There is no cone, so a front is never
           *    stopped for turning - only for running out of mesh to extrude into. ---*/
          const su2double dot = GeometryToolbox::DotProduct(nDim, vec, dirNow[f].data());

          /*--- Halo nodes stay out: their parent is dictated by the rank that owns them and arrives
           *    through the MPI relay, so a front claiming one would fight that assignment. This is
           *    what a front hits when it reaches a partition interface, and it needs a reason of its
           *    own, rather than being skipped silently and leaving some other candidate to explain a
           *    stop that was really the partitioning. ---*/
          if (!fine_grid->nodes->GetDomain(jPoint)) {
            sawPartition = true;
            /*--- Held as a handover candidate, subject to the same admissibility the owner would
             *    apply anyway; whether it is still free is the owner's to decide. ---*/
            if (dot > bestHalo_dot && !(onPhysicalBoundary[jPoint] && entersBoundary(jPoint, vec)) &&
                GeometricalCheck(jPoint, fine_grid, config)) {
              bestHalo_dot = dot;
              bestHalo = jPoint;
            }
            continue;
          }
          if (fine_grid->nodes->GetAgglomerate(jPoint)) {
            /*--- A node another front has already turned into a coarse CV also reads as agglomerated,
             *    so ask frontOf first: that is the fronts meeting, not an earlier phase. ---*/
            if ((frontOf[jPoint] >= 0) && (frontOf[jPoint] != static_cast<int>(f)))
              sawCollision = true;
            else if (frontOf[jPoint] < 0)
              sawAgglom = true;
            continue;
          }
          if (claimed[jPoint]) {
            if (frontOf[jPoint] != static_cast<int>(f)) sawCollision = true;
            continue;
          }

          if (onPhysicalBoundary[jPoint] && entersBoundary(jPoint, vec)) {
            sawBoundary = true;
            continue;
          }
          if (!GeometricalCheck(jPoint, fine_grid, config)) {
            sawGeom = true;
            continue;
          }

          if (dot > best_dot) {
            best_dot = dot;
            best = jPoint;
            best_len = len;
            for (unsigned short d = 0; d < nDim; ++d) best_dir[d] = vec[d];
          }
        }

        if ((best == NO_POINT) && (bestHalo != NO_POINT)) {
          /*--- Nowhere left on this rank, but the stack does continue - just on someone else's side
           *    of the interface. Remember where, and let the classification below decide whether the
           *    whole layer goes over. ---*/
          handTo[f].push_back(bestHalo);
          continue;
        }

        if (best == NO_POINT) {
          /*--- Priority order picks the cleanest explanation first: reaching a physical boundary is a
           *    correct, expected stop and takes priority even if some other, non-viable candidate also
           *    happened to be claimed. Only report a collision when no boundary was involved. ---*/
          if (sawBoundary)
            markFail(f, STOP_PHYS_BOUNDARY);
          else if (sawPartition) {
            nNoHalo++;
            markFail(f, STOP_PARTITION);
          } else if (sawCollision)
            markFail(f, STOP_COLLISION);
          else if (sawAgglom)
            markFail(f, STOP_AGGLOMERATED);
          else if (sawGeom)
            markFail(f, STOP_GEOMETRY);
          else
            markFail(f, STOP_NO_NEIGHBOR);
          break;
        }

        CStep s{};
        s.node = best;
        s.from = n;
        s.key = fine_grid->nodes->GetGlobalIndex(n);
        s.score = best_dot;
        s.dist = best_len;
        for (unsigned short d = 0; d < nDim; ++d) s.dir[d] = best_dir[d];
        prop[f].push_back(s);
      }

      /*--- A footprint that reaches an interface can be cut by it. If the WHOLE footprint crosses,
       *    the stack is handed over intact and this front is finished. If only part of it crosses -
       *    the interface running ALONG the stack rather than across it - the footprint is SPLIT: the
       *    piece whose successors this rank owns marches on here as a narrower stack, and the rest
       *    goes to the rank owning the nodes it was reaching for. Both pieces are renamed, because
       *    they are separate stacks from here on.
       *
       *    Retiring on a cut, which is what this did before, was the expensive part of partitioning:
       *    an interface that cuts a stack at layer k cuts it at every layer above k as well, so one
       *    straddle did not cost one layer, it cost the whole remaining column - about a hundred
       *    nodes each on the flat plate. ---*/
      if (failed[f]) {
        prop[f].clear();
        handTo[f].clear();
      } else if (!handTo[f].empty()) {
        /*--- prop[f] is built in the order of front[f], so this is the piece that stays, in the same
         *    order, and phi still runs index for index between it and the layer it proposes. ---*/
        vector<unsigned long> narrow;
        for (const auto& s : prop[f]) narrow.push_back(s.from);

        /*--- A cut can leave the local piece in two disconnected halves - a square footprint cut
         *    diagonally does exactly that - and that is not a layer. Drop it and hand over the rest;
         *    the stack still survives on the far side instead of ending here. ---*/
        if (!narrow.empty() && !isConnectedLayer(narrow)) {
          nSplitDropped += narrow.size();
          narrow.clear();
        }

        handTag[f] = tagOfSet(handTo[f]);

        if (narrow.empty()) {
          prop[f].clear();
        } else {
          nSplit++;
          nSplitLocal += narrow.size();
          nSplitHanded += handTo[f].size();
          /*--- Close the coarse CV that is open on the WIDE footprint before narrowing, so that no CV
           *    ever ends up holding two layers of different shape. ---*/
          front[f] = narrow;
          emit(f);
          nBlock[f] = blockFor(front[f]);
          tag[f] = tagOfSet(front[f]);
          keepLocal[f] = 1;
        }
      }
    }

    /*--- (b) Contention resolved from bids that were all collected before any was granted, so the
     *    outcome is a pure function of the proposals and does not depend on the order the fronts are
     *    visited in. That is what makes the coarse grid reproducible. ---*/
    auto better = [](const CStep& a, const CStep& b) {
      if (a.score != b.score) return a.score > b.score;
      if (a.dist != b.dist) return a.dist < b.dist;
      return a.key < b.key;
    };

    for (unsigned long f = 0; f < front.size(); ++f) {
      if (!alive[f] || failed[f]) continue;
      for (const auto& s : prop[f]) {
        /*--- A front that has just lost a bid is retiring, so it must not go on to place the rest and
         *    displace a front that is still healthy. What it placed BEFORE losing does stay in the
         *    table, and can still cost another front a candidate it would otherwise have won: the
         *    only way to avoid that entirely is to re-run the contention to a fixed point after every
         *    retirement. The residual is conservative - it retires a front near a seam one layer
         *    early, never merges anything it should not - and the seam goes to ordinary agglomeration
         *    either way, so it is not worth an inner iteration. ---*/
        if (failed[f]) break;

        if (bidIdx[s.node] == NOBID) {
          bidIdx[s.node] = static_cast<unsigned>(bids.size());
          bids.push_back(s);
          bidOwner.push_back(f);
          continue;
        }

        const auto k = bidIdx[s.node];
        const auto g = bidOwner[k];
        /*--- Two nodes of the SAME front reaching for one successor is a pinch: the layer would come
         *    out narrower than the front, which all-or-nothing does not allow. ---*/
        const unsigned short why = (g == f) ? STOP_PINCH : STOP_COLLISION;

        if (better(s, bids[k])) {
          markFail(g, why);
          bids[k] = s;
          bidOwner[k] = f;
        } else {
          markFail(f, why);
        }

        /*--- A head-on meeting stops BOTH fronts. Letting the winner carry on through the seam would
         *    push its stack into territory the other front had every right to, and the asymmetry
         *    shows up in the coarse grid as one stack overshooting the other. A glancing contact
         *    (directions not opposed) is not a meeting and only costs the loser. ---*/
        if ((g != f) && (GeometryToolbox::DotProduct(nDim, dirNow[f].data(), dirNow[g].data()) < 0.0)) {
          markFail(f, STOP_COLLISION);
          markFail(g, STOP_COLLISION);
        }
      }
    }

    /*--- (c) All-or-nothing acceptance: a front takes the whole layer or none of it and retires. A
     *    bid only becomes a claim here, so a retiring front never has to give anything back and the
     *    nodes it was reaching for stay available to ordinary agglomeration. ---*/
    for (unsigned long f = 0; f < front.size(); ++f) {
      if (!alive[f]) continue;

      newLayer.clear();
      if (!failed[f]) {
        /*--- Built in the order prop[f] was, which is the order of front[f], so newLayer[k] is the
         *    successor proposed by front[f][k] and the two vectors carry phi between them. ---*/
        for (const auto& s : prop[f]) {
          const auto k = bidIdx[s.node];
          if ((k != NOBID) && (bidOwner[k] == f)) newLayer.push_back(s.node);
        }
        /*--- Every bid of a front that was not marked failed must have been granted. ---*/
        if (newLayer.size() != front[f].size())
          markFail(f, STOP_COLLISION);
        else if (!layerIsIsomorphic(front[f], newLayer))
          markFail(f, STOP_TOPOLOGY);
      }

      if (failed[f]) {
        /*--- Nothing to give back: a bid only becomes a claim on acceptance below. ---*/
        alive[f] = 0;
        stopCounts[failReason[f]]++;
        /*--- One layer short of a full block at the top: take what is buffered as its own coarse CV
         *    rather than dropping it back to ordinary agglomeration. ---*/
        emit(f);
        continue;
      }

      /*--- Accept. The direction is blended rather than replaced, and it is the FRONT's direction,
       *    updated once from the mean of the steps its nodes just took, not one direction per node
       *    free to wander off on its own. ---*/
      su2double mean[MAXNDIM] = {0.0};
      for (const auto& s : prop[f])
        for (unsigned short d = 0; d < nDim; ++d) mean[d] += s.dir[d];
      const su2double meanNrm = GeometryToolbox::Norm(nDim, mean);
      if (meanNrm > 0.0) {
        su2double blended[MAXNDIM] = {0.0};
        for (unsigned short d = 0; d < nDim; ++d)
          blended[d] = (1.0 - DIR_BLEND) * dirNow[f][d] + DIR_BLEND * mean[d] / meanNrm;
        const su2double bNrm = GeometryToolbox::Norm(nDim, blended);
        if (bNrm > 0.0)
          for (unsigned short d = 0; d < nDim; ++d) dirNow[f][d] = blended[d] / bNrm;
      }

      for (auto p : newLayer) {
        claimed[p] = 1;
        frontOf[p] = static_cast<int>(f);
      }
      front[f] = std::move(newLayer);
      depth[f]++;
      nLayers++;

      pending[f].insert(pending[f].end(), front[f].begin(), front[f].end());
      pendingLayers[f]++;
      if (pendingLayers[f] >= nBlock[f]) emit(f);
    }

    /*==============================================================================================
     *  (d) Hand stacks across partition interfaces.
     *
     *  A front that has run into the halo cannot go on here: those nodes belong to another rank, and
     *  their parent is that rank's to assign. Without this the stack simply ended at the interface,
     *  and on four ranks that was 70 of 71 fronts - the paved fraction of the flat plate fell from
     *  99% to 48% purely because of where the partition happened to cut.
     *
     *  Instead the footprint is sent to the owner, which picks the stack up and carries on. What
     *  crosses is not a coarse CV - a CV belongs wholly to one rank - but the FOOTPRINT, so the two
     *  halves of the stack stay the same shape and the coarse grid reads the same across the seam.
     *  The handover travels the reverse of the usual halo direction: a node this rank sees as halo is
     *  one the neighbour owns, so it is packed against the RECEIVE marker and sent to the rank this
     *  marker normally receives from.
     *============================================================================================*/
    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
      if (!((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)))
        continue;

      const auto MarkerS = iMarker, MarkerR = iMarker + 1;
      const auto send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
      const auto receive_from = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
      const auto nVertexS = fine_grid->nVertex[MarkerS];
      const auto nVertexR = fine_grid->nVertex[MarkerR];

      /*--- Packed against the halo vertices, i.e. what this rank wants the neighbour to continue. ---*/
      vector<unsigned long> tagOut(nVertexR, 0), tagIn(nVertexS, 0);
      vector<su2double> dirOut(nVertexR * nDim, 0.0), dirIn(nVertexS * nDim, 0.0);

      for (unsigned long f = 0; f < front.size(); ++f) {
        if (handTo[f].empty()) continue;
        for (auto p : handTo[f]) {
          if (haloMarker[p] != static_cast<int>(MarkerR)) continue;
          const auto v = haloVertex[p];
          /*--- Two fronts of this rank reaching for the same node: the lower tag takes it, which is
           *    a decision both ranks would reach the same way. ---*/
          if ((tagOut[v] != 0) && (tagOut[v] <= handTag[f])) continue;
          tagOut[v] = handTag[f];
          for (unsigned short d = 0; d < nDim; ++d) dirOut[v * nDim + d] = dirNow[f][d];
        }
      }

      SU2_MPI::Sendrecv(tagOut.data(), nVertexR, MPI_UNSIGNED_LONG, receive_from, 2, tagIn.data(), nVertexS,
                        MPI_UNSIGNED_LONG, send_to, 2, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
      SU2_MPI::Sendrecv(dirOut.data(), nVertexR * nDim, MPI_DOUBLE, receive_from, 3, dirIn.data(), nVertexS * nDim,
                        MPI_DOUBLE, send_to, 3, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      for (auto iVertex = 0ul; iVertex < nVertexS; iVertex++) {
        if (tagIn[iVertex] == 0) continue;
        CInherited h;
        h.tag = tagIn[iVertex];
        h.node = fine_grid->vertex[MarkerS][iVertex]->GetNode();
        for (unsigned short d = 0; d < nDim; ++d) h.dir[d] = dirIn[iVertex * nDim + d];
        inherited.push_back(h);
      }
    }

    /*--- A front that handed its WHOLE footprint over is finished here; the neighbour owns the rest
     *    of the stack. One that handed over only a piece keeps marching on what was left of it. ---*/
    for (unsigned long f = 0; f < front.size(); ++f) {
      if (handTo[f].empty()) continue;
      handTo[f].clear();
      nHandedOut++;
      if (keepLocal[f]) continue;
      alive[f] = 0;
      stopCounts[STOP_PARTITION]++;
      emit(f);
    }

    /*--- Adopt what the neighbours sent. Tags are processed in ascending order so that two ranks
     *    handing stacks onto overlapping nodes are separated the same way whatever order the messages
     *    happened to arrive in. A footprint whose nodes are not all still free is dropped: the stack
     *    simply ends, exactly as it would have before. ---*/
    std::sort(inherited.begin(), inherited.end(),
              [](const CInherited& a, const CInherited& b) { return a.tag < b.tag; });

    for (size_t i = 0; i < inherited.size();) {
      size_t j = i;
      while ((j < inherited.size()) && (inherited[j].tag == inherited[i].tag)) ++j;

      vector<unsigned long> layer0;
      bool ok = true;
      for (size_t k = i; k < j; ++k) {
        const auto p = inherited[k].node;
        if (claimed[p] || fine_grid->nodes->GetAgglomerate(p) || !GeometricalCheck(p, fine_grid, config)) ok = false;
        layer0.push_back(p);
      }
      /*--- The footprint has to arrive whole and connected, the same test any other layer passes. ---*/
      if (ok && !isConnectedLayer(layer0)) ok = false;

      if (ok) {
        std::array<su2double, MAXNDIM> d0{};
        for (unsigned short d = 0; d < nDim; ++d) d0[d] = inherited[i].dir[d];
        /*--- An inherited layer is an interior one, so it is NOT subject to the single-layer rule the
         *    boundary layer gets: it opens an ordinary two-deep coarse CV and waits for its partner. ---*/
        const auto nf = addFront(layer0, d0, inherited[i].tag, blockFor(layer0));
        for (auto p : layer0) {
          claimed[p] = 1;
          frontOf[p] = static_cast<int>(nf);
        }
        nLayers++;
        nHandedIn++;
        if (pendingLayers[nf] >= nBlock[nf]) emit(nf);
      }
      i = j;
    }
    inherited.clear();
  }

  /*--- Nothing should be left buffered, but a front retired outside the loop would strand its nodes
   *    with a parent index that was never assigned. ---*/
  for (unsigned long f = 0; f < front.size(); ++f) emit(f);

  /*--- How far each front actually got. This is the number to watch when the paved region does not
   *    look like a front: fronts that all reach the same height leave a flat interface with ordinary
   *    agglomeration, and a spread here is that interface coming out as a staircase instead. ---*/
  unsigned long dmin = std::numeric_limits<unsigned long>::max(), dmax = 0;
  for (unsigned long f = 0; f < front.size(); ++f) {
    if (front[f].empty()) continue;
    dmin = std::min(dmin, depth[f]);
    dmax = std::max(dmax, depth[f]);
  }
  /*--- A rank with no fronts of its own must not drag the reported minimum to zero: leaving dmin at
   *    its sentinel keeps it out of the MPI_MIN below, so the range describes the fronts that exist
   *    rather than the ranks that have none. ---*/
  if (getenv("DUMPAGGLOM") != nullptr)
    cout << "  Paving rank " << rank << ": " << front.size() << " fronts, depth " << dmin << ".." << dmax
         << ", local nPointDomain " << fine_grid->GetnPointDomain() << endl;

  /*--- Summary over all ranks. Reporting rank 0's own fronts makes a partitioned run look like a
   *    fraction of the mesh it is not, and hides how much of the layer the partitioning cost: a front
   *    stops at the partition, so the number of nodes left to ordinary agglomeration is the number to
   *    watch when adding ranks. Every rank must reach these collectives. ---*/
  unsigned long nSeedNodes = seeds.node.size();
  unsigned long local[13] = {nSeedNodes, nStacks, nLayers,     nCovered,     nSemiCV, nFullCV,      nHandedOut,
                             nHandedIn,  nSplit,  nSplitLocal, nSplitHanded, nNoHalo, nSplitDropped};
  unsigned long total[13] = {0};
  SU2_MPI::Allreduce(local, total, 13, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  unsigned long localCV = Index_CoarseCV - starting_Index_CoarseCV, totalCV = 0;
  SU2_MPI::Allreduce(&localCV, &totalCV, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  unsigned long histTotal[9] = {0};
  SU2_MPI::Allreduce(histogram, histTotal, 9, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  vector<unsigned long> stopTotal(N_STOP_REASONS, 0);
  SU2_MPI::Allreduce(stopCounts.data(), stopTotal.data(), N_STOP_REASONS, MPI_UNSIGNED_LONG, MPI_SUM,
                     SU2_MPI::GetComm());

  unsigned long depthMin = 0, depthMax = 0;
  SU2_MPI::Allreduce(&dmin, &depthMin, 1, MPI_UNSIGNED_LONG, MPI_MIN, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&dmax, &depthMax, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  if (depthMin == std::numeric_limits<unsigned long>::max()) depthMin = 0; /*--- No fronts anywhere. ---*/

  if (rank == MASTER_NODE) {
    cout << "  Paving fronts: " << total[1] << " fronts from " << total[0] << " seed nodes, patch sizes ";
    for (unsigned s = 1; s <= 8; ++s)
      if (histTotal[s] > 0) cout << s << "x" << histTotal[s] << " ";
    cout << "\n  Coarse CVs from fronts: " << totalCV << " covering " << total[3] << " nodes in " << total[2]
         << " layers, front depth " << depthMin << " to " << depthMax;
    if (total[4] + total[5] > 0)
      cout << "\n  Coarse CVs by depth: " << total[5] << " two layers deep, " << total[4]
           << " one layer (top of a stack)";
    if (total[6] + total[7] > 0)
      cout << "\n  Stacks handed across partitions: " << total[6] << " sent, " << total[7] << " picked up";
    if (total[8] > 0)
      cout << "\n  Footprints split at partitions: " << total[8] << " cut by an interface (" << total[9]
           << " nodes marching on here, " << total[10] << " handed across)";
    if (total[11] + total[12] > 0)
      cout << "\n  Stacks lost at partitions: " << total[11] << " blocked with nowhere to hand to, " << total[12]
           << " nodes in split pieces that came apart";
    /*--- Why each front stopped. Reaching a boundary is the one correct stop; everything else is the
     *    mesh failing to offer a layer topologically identical to the current one, broken down by how
     *    it failed. COLLISION and PINCH are two fronts, or two nodes of one front, reaching for the
     *    same node; TOPOLOGY is a layer that was claimable but not an extrusion. ---*/
    cout << "\n  Front advance stopped due to: physical-boundary " << stopTotal[STOP_PHYS_BOUNDARY] << ", partition "
         << stopTotal[STOP_PARTITION] << ", front-collision " << stopTotal[STOP_COLLISION] << ", pinch "
         << stopTotal[STOP_PINCH] << ", already-agglomerated " << stopTotal[STOP_AGGLOMERATED] << ", dead-end "
         << stopTotal[STOP_NO_NEIGHBOR] << ", topology " << stopTotal[STOP_TOPOLOGY] << ", geometry "
         << stopTotal[STOP_GEOMETRY] << ", max-length " << stopTotal[STOP_MAX_LENGTH];
    cout << endl;
  }
}
