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

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace {

/*--- Euler wall nodes are not agglomerated where the surface turns by more than this, in
 *    degrees. ---*/
constexpr passivedouble EULER_WALL_MAX_CURVATURE = 45.0;

/*--- Equivalence-class id per entity for the set of physical markers it lies on: entities with the
 *    same set share an id, 0 means no marker. Only set equality is ever asked, so the sets are
 *    interned and compared as ids, which puts no limit on the number of markers. The (entity,
 *    marker) pairs may arrive in any order and are consumed. ---*/
vector<unsigned long> MarkerSetClasses(unsigned long nEntity, vector<std::pair<unsigned long, unsigned short>>& pairs) {
  vector<unsigned long> classOfEntity(nEntity, 0);

  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

  map<vector<unsigned short>, unsigned long> classOf;
  unsigned long nClass = 0;
  for (size_t i = 0; i < pairs.size();) {
    size_t j = i;
    vector<unsigned short> markerSet;
    while ((j < pairs.size()) && (pairs[j].first == pairs[i].first)) markerSet.push_back(pairs[j++].second);

    const auto res = classOf.emplace(std::move(markerSet), nClass + 1);
    if (res.second) nClass++;
    classOfEntity[pairs[i].first] = res.first->second;
    i = j;
  }
  return classOfEntity;
}

}  // namespace

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

  /*--- Points carrying a physical boundary condition. This does not include SEND_RECEIVE. ---*/
  vector<char> onPhysBoundary(fine_grid->GetnPoint(), 0);
  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++)
      onPhysBoundary[fine_grid->vertex[iMarker][iVertex]->GetNode()] = 1;
  }
  /*--- Nodes where two different boundary conditions meet.  ---*/
  const auto mixedBC = FindMixedBoundaryNodes(fine_grid, config);

  /*--- STEP 0: pave the domain with advancing fronts rising from the boundaries. The coarse CVs it
   *    creates occupy [firstLineCV, endLineCV). ---*/
  const auto firstLineCV = Index_CoarseCV;
  vector<unsigned long> neverGrewCV;
  if (config->GetMGOptions().MG_Implicit_Lines) {
    pavingReport = PaveAdvancingFronts(Index_CoarseCV, fine_grid, config, iMesh, mixedBC, onPhysBoundary, neverGrewCV);
  }
  const auto endLineCV = Index_CoarseCV;

  /*--- STEP 1: The first step is the boundary agglomeration. ---*/
  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    /*--- Skip periodic boundaries: do not agglomerate on periodic markers. ---*/
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) continue;

    /*--- Skip SEND_RECEIVE markers, those points are left to the domain pass. ---*/
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
         that point.  ---*/

        for (auto jMarker = 0u; jMarker < fine_grid->GetnMarker(); jMarker++) {
          if (config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) continue;
          if (fine_grid->nodes->GetVertex(iPoint, jMarker) != -1) {
            /*--- Count every physical marker, but store only the first few. ---*/
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
              if (local_curvature >= EULER_WALL_MAX_CURVATURE) {
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
          /*--- In 2D that is a corner in the geometry, which is never agglomerated. ---*/
          if (nDim == 2) agglomerate_seed = false;
          /*--- In 3D, this is a ridge point (an edge where two surface markers meet). ---*/
          if (nDim == 3) agglomerate_seed = true;

          /*--- Euler walls: check curvature-based agglomeration criterion for both markers ---*/
          // only in 3d because in 2d it's a corner
          bool euler_wall_rejected_here = false;
          for (unsigned short i = 0; i < 2; i++) {
            if ((nDim == 3) && (config->GetMarker_All_KindBC(copy_marker[i]) == EULER_WALL)) {
              if (!boundIsStraight[copy_marker[i]]) {
                /*--- Compute local curvature at this point ---*/
                su2double local_curvature = ComputeLocalCurvature(fine_grid, iPoint, copy_marker[i]);
                if (local_curvature >= EULER_WALL_MAX_CURVATURE) {
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

        /*--- ...and so is a node where two markers of different type meet, in any dimension. ---*/
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
              /*--- In 2D, we agglomerate exactly 2 nodes if the nodes are on the line edge. ---*/
              if ((nDim == 2) && (counter == 1)) break;
              /*--- In 3D, we agglomerate exactly 2 nodes if the nodes are on the surface edge. ---*/
              if ((nDim == 3) && (counter == 2)) break;
              /*--- Apply maxAgglomSize limit for 3D internal boundary face nodes (counter==1 in 3D). ---*/
              if (nChildren >= maxAgglomSize) break;
            }
          }

          /*--- Indirect neighbors only for 3D faces.  ---*/
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
                /*--- Apply maxAgglomSize limit for 3D internal boundary face nodes. ---*/
                if (nChildren >= maxAgglomSize) break;
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
    /*--- As in STEP 1, a SEND_RECEIVE marker does not make a point a boundary point. ---*/
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

  /*--- STEP 2: Agglomerate the domain points. A seed grows one node at a time, taking the candidate
   *    that shares the most edges with the current members. ---*/
  vector<char> inCV(fine_grid->GetnPoint(), 0);
  vector<char> isCandidate(fine_grid->GetnPoint(), 0);
  vector<unsigned long> members, candidates;
  members.reserve(maxAgglomSize);

  /*--- A local frame at the seed: up to nDim of its incident edges, as mutually orthogonal as
   *    possible, each with its own length. ---*/
  vector<std::array<su2double, MAXNDIM>> frameDir, edgeDir;
  vector<su2double> frameLen, edgeLen;
  vector<char> edgeUsed;

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

      /*--- The shortest edge goes in first, the rest in order of how orthogonal they are to what is
       *    already in the frame. ---*/
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

          /*--- Distance to the centroid in cells: the offset is scaled by the seed edge pointing
           *    most nearly along it. ---*/
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

  /*--- The connectivity just built only knows about coarse CVs of this rank, so a CV touching a
   partition boundary may look isolated. Mark those CVs and leave them alone. ---*/

  /*--- Coarse CVs to leave exactly as the agglomeration made them: those holding a node where two
   *    different boundary conditions meet. ---*/
  vector<bool> mustStayAlone(nPointDomain, false);
  /*--- ...and which coarse CVs hold a boundary node at all. ---*/
  vector<bool> cvOnBoundary(nPointDomain, false);
  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++)
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren++) {
      const auto iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      if (mixedBC[iFinePoint]) mustStayAlone[iCoarsePoint] = true;
      if (onPhysBoundary[iFinePoint]) cvOnBoundary[iCoarsePoint] = true;
    }

  /*--- Which physical boundaries each coarse CV sits on, as a marker-set class. Comparing classes
   *    keeps a merge inside one boundary. Halo points hold the parent sentinel and drop out. ---*/
  vector<std::pair<unsigned long, unsigned short>> cvMarker;
  for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    for (auto iVertex = 0ul; iVertex < fine_grid->nVertex[iMarker]; iVertex++) {
      const auto iParent = fine_grid->nodes->GetParent_CV(fine_grid->vertex[iMarker][iVertex]->GetNode());
      if (iParent < nPointDomain) cvMarker.push_back({iParent, static_cast<unsigned short>(iMarker)});
    }
  }
  const auto cvMarkerClass = MarkerSetClasses(nPointDomain, cvMarker);

  /*--- A CV whose front never advanced past its seed is not a stack: protecting it anyway can
   *    leave it many orders of magnitude smaller than its neighbours (its footprint is a single
   *    fine cell, most often on a thin near-wall layer), which is unstable once its correction is
   *    carried up through further coarsening. Let it fall through to ordinary repair instead. ---*/
  vector<bool> neverGrew(nPointDomain, false);
  for (auto iCV : neverGrewCV)
    if (iCV < nPointDomain) neverGrew[iCV] = true;

  /*--- A boundary CV built by the paving is the base of a stack and keeps its footprint, so the
   *    repair passes below leave it alone -- unless it never grew into one, see above. ---*/
  auto isStackBase = [&](unsigned long iCoarsePoint) {
    return cvOnBoundary[iCoarsePoint] && (iCoarsePoint >= firstLineCV) && (iCoarsePoint < endLineCV) &&
           !neverGrew[iCoarsePoint];
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
      if (cvMarkerClass[iCoarsePoint] != cvMarkerClass[iCoarsePoint_Complete]) continue;

      /*--- Check if merging would exceed the maximum agglomeration size ---*/
      auto nChildren_Target = nodes->GetnChildren_CV(iCoarsePoint_Complete);
      auto nChildren_Isolated = nodes->GetnChildren_CV(iCoarsePoint);
      auto nChildren_Total = nChildren_Target + nChildren_Isolated;

      /*--- If the total would exceed maxAgglomSize, try to redistribute children to neighbors. The
       merge below runs whether or not the quota is met, so the limit can still be exceeded. ---*/
      if (nChildren_Total > maxAgglomSize) {
        /*--- Find neighbors of the target coarse point that have room ---*/
        unsigned short nChildrenToRedistribute = nChildren_Total - maxAgglomSize;

        for (auto jCoarsePoint : nodes->GetPoints(iCoarsePoint_Complete)) {
          if (nChildrenToRedistribute == 0) break;
          /*--- The isolated CV is a neighbour of the target and hands anything it takes straight
           *    back below, spending the quota without lowering the count. ---*/
          if (jCoarsePoint == iCoarsePoint) continue;
          if (mustStayAlone[jCoarsePoint]) continue;
          if (isStackBase(jCoarsePoint)) continue;
          if (cvMarkerClass[jCoarsePoint] != cvMarkerClass[iCoarsePoint_Complete]) continue;

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

  /*--- Merge a coarse CV that still holds a single fine child into whichever coarse neighbor has the
   fewest children. Both the merged CV and the target are owned by this rank. ---*/

  for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++) {
    if (nodes->GetnChildren_CV(iCoarsePoint) != 1) continue;
    if (mustStayAlone[iCoarsePoint]) continue;
    if (isStackBase(iCoarsePoint)) continue;
    if (nodes->GetnPoint(iCoarsePoint) <= 1) continue; /*--- Already handled above, or truly islanded. ---*/

    /*--- Pick the neighbour with the fewest children. The target takes the child even when it is
     already at maxAgglomSize, so the limit is exceeded by one rather than a single-child CV
     surviving. ---*/
    unsigned long best_neighbor = std::numeric_limits<unsigned long>::max();
    unsigned short best_nChildren = 0;
    for (auto jCoarsePoint : nodes->GetPoints(iCoarsePoint)) {
      if (mustStayAlone[jCoarsePoint]) continue;
      if (isStackBase(jCoarsePoint)) continue;
      if (cvMarkerClass[jCoarsePoint] != cvMarkerClass[iCoarsePoint]) continue;
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

  /*--- Compact the coarse numbering, squeezing out the indices the repair passes emptied. The
   children lists, indirect-agglomeration flags and owned parent indices are remapped. ---*/

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
        nodes->SetChildren_CV(iNew, nodes->GetChildren_CV(iCoarsePoint));
        nodes->SetnChildren_CV(iNew, nodes->GetnChildren_CV(iCoarsePoint));
        nodes->SetAgglomerate_Indirect(iNew, nodes->GetAgglomerate_Indirect(iCoarsePoint));
      }
      for (auto iCoarsePoint = nKept; iCoarsePoint < nPointDomain; iCoarsePoint++)
        nodes->SetnChildren_CV(iCoarsePoint, 0);

      for (auto iFinePoint = 0ul; iFinePoint < fine_grid->GetnPointDomain(); iFinePoint++) {
        const auto iParent = fine_grid->nodes->GetParent_CV(iFinePoint);
        if (iParent == NO_INDEX) continue;
        /*--- A fine point may only reference a CV the compaction kept. ---*/
        if ((iParent >= nPointDomain) || (newIndex[iParent] == NO_INDEX))
          SU2_MPI::Error("Multigrid compaction: a fine point still references an emptied coarse CV.",
                         CURRENT_FUNCTION);
        fine_grid->nodes->SetParent_CV(iFinePoint, newIndex[iParent]);
      }

      nPointDomain = nKept;
      nPoint = nKept;
      Index_CoarseCV = nKept;
    }
  }

  /*--- Reset the neighbor information. ---*/

  nodes->ResetPoints();

#ifdef HAVE_MPI
  /*--- Reset halo point parents before MPI agglomeration, the fine grid still carries the parent
   indices it was given when it was itself built. ---*/

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

      /*--- First pass: determine which parents will actually be used, i.e. have non-skipped
       children. ---*/
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

  /*--- Stop coarsening once the smallest per-rank partition falls below the minimum, not just the
        summed total. ---*/
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

vector<char> CMultiGridGeometry::FindMixedBoundaryNodes(const CGeometry* fine_grid, const CConfig* config) const {
  vector<char> mixed(fine_grid->GetnPoint(), 0);

  /*--- The first physical condition seen at each node, -1 until one is. A second one of a different
   *    kind is what makes the node mixed. ---*/
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
      /*--- Identify the physical markers of the vertex that we want to agglomerate. A candidate whose
       markers are all SEND_RECEIVE ends up with counter == 0 and is rejected below. ---*/

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
        /*--- Skip neighbors whose parent is not known yet, halo points still hold the sentinel until
         the MPI relay has run. ---*/
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
  /*--- Local curvature is the maximum angle between adjacent face normals at a boundary vertex. ---*/

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

namespace {

/*--- Unit normal of a boundary at a vertex, false if the marker does not reach iPoint. Boundary
 *    normals point into the domain. ---*/
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

/*--- Are a and b neighbours in the fine grid? ---*/
bool IsAdjacent(const CGeometry* grid, unsigned long a, unsigned long b) {
  const auto& pts = grid->nodes->GetPoints(a);
  return std::find(pts.begin(), pts.end(), b) != pts.end();
}

/*--- Is a footprint one connected patch? ---*/
bool IsConnectedLayer(const CGeometry* grid, const vector<unsigned long>& layer) {
  if (layer.size() < 2) return true;
  vector<char> seen(layer.size(), 0);
  vector<size_t> stk{0};
  seen[0] = 1;
  size_t nSeen = 1;
  while (!stk.empty()) {
    const auto cur = stk.back();
    stk.pop_back();
    for (size_t k = 0; k < layer.size(); ++k) {
      if (seen[k] || !IsAdjacent(grid, layer[cur], layer[k])) continue;
      seen[k] = 1;
      nSeen++;
      stk.push_back(k);
    }
  }
  return nSeen == layer.size();
}

/*--- Is the new layer topologically identical to the old? The layers are index-aligned, so the map
 *    old[k] to new[k] has to be an isomorphism of the induced subgraphs. ---*/
bool LayerIsIsomorphic(const CGeometry* grid, const vector<unsigned long>& oldL, const vector<unsigned long>& newL) {
  const auto n = oldL.size();
  if (newL.size() != n) return false;

  for (size_t k = 0; k < n; ++k) {
    unsigned nOld = 0, nNew = 0;
    for (size_t l = 0; l < n; ++l) {
      nOld += IsAdjacent(grid, newL[k], oldL[l]);
      nNew += IsAdjacent(grid, oldL[k], newL[l]);
    }
    /*--- Exactly one partner each way, and it has to be the one phi names. ---*/
    if ((nOld != 1) || (nNew != 1)) return false;
    if (!IsAdjacent(grid, oldL[k], newL[k])) return false;
  }

  for (size_t k = 0; k < n; ++k)
    for (size_t l = k + 1; l < n; ++l)
      if (IsAdjacent(grid, oldL[k], oldL[l]) != IsAdjacent(grid, newL[k], newL[l])) return false;

  return true;
}

/*--- Does the step run into a boundary at jPoint, i.e. roughly along that boundary's normal? ---*/
bool EntersBoundary(const CGeometry* grid, const CConfig* config, unsigned short nDim, unsigned long jPoint,
                    const su2double* stepDir, su2double cosBoundary) {
  for (unsigned short iMarker = 0; iMarker < grid->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;
    su2double n[3] = {0.0}; /*--- nDim is at most 3. ---*/
    if (!VertexUnitNormal(grid, nDim, jPoint, iMarker, n)) continue;
    if (fabs(GeometryToolbox::DotProduct(nDim, n, stepDir)) >= cosBoundary) return true;
  }
  return false;
}

/*--- Fine layers the next coarse CV of this front holds: two, or one if that would exceed the
 *    agglomeration size limit. ---*/
unsigned long BlockFor(short int maxAgglomSize, const vector<unsigned long>& layer) {
  return (layer.size() * 2 > static_cast<size_t>(maxAgglomSize)) ? 1 : 2;
}

/*--- Rank-independent name for a set of nodes: the smallest global index in it, +1 so that 0 means
 *    "nothing". ---*/
unsigned long TagOfSet(const CGeometry* grid, const vector<unsigned long>& set) {
  unsigned long t = std::numeric_limits<unsigned long>::max();
  for (auto p : set) t = std::min(t, grid->nodes->GetGlobalIndex(p));
  return t + 1;
}

}  // namespace

CMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeometry* fine_grid,
                                                                   const CConfig* config) const {
  /*--- Fraction of a marker's nodes that must sit in a layer before the whole marker may seed. ---*/
  constexpr passivedouble QUALIFIED_FRACTION = 0.5;
  constexpr passivedouble ANGLE_THRESHOLD_DEG = 30.0;
  /*--- Smallest local cell aspect ratio for which a node still counts as part of a stretched
   *    layer. ---*/
  constexpr passivedouble MIN_AR = 2.0;
  const su2double cos_threshold = cos(ANGLE_THRESHOLD_DEG * PI_NUMBER / 180.0);

  const auto nMarkerFine = fine_grid->GetnMarker();
  constexpr auto NO_POINT = std::numeric_limits<unsigned long>::max();

  CFrontSeeds seeds;
  vector<char> taken(fine_grid->GetnPoint(), 0);

  auto isWall = [&](unsigned short bc) {
    return (bc == HEAT_FLUX) || (bc == ISOTHERMAL) || (bc == CHT_WALL_INTERFACE) || (bc == SMOLUCHOWSKI_MAXWELL);
  };


  /*--- True if the mesh at iPoint is stretched along the boundary normal, i.e. this boundary has a
   *    layer growing off it the way a viscous wall does. The coupling across the dual face between
   *    iPoint and a neighbour is measured here, so the ratio of largest to smallest weight is the
   *    local aspect ratio. Only boundary nodes are ever asked, at most twice each. ---*/
  auto hasLayerNormalTo = [&](unsigned long iPoint, const su2double* unitNormal) {
    su2double wMin = std::numeric_limits<su2double>::max(), wMax = 0.0;
    auto jStiffest = NO_POINT;
    for (auto iNeigh = 0u; iNeigh < fine_grid->nodes->GetnPoint(iPoint); ++iNeigh) {
      const auto jPoint = fine_grid->nodes->GetPoint(iPoint, iNeigh);
      const auto iEdge = fine_grid->nodes->GetEdge(iPoint, iNeigh);
      const su2double area = GeometryToolbox::Norm(nDim, fine_grid->edges->GetNormal(iEdge));
      const su2double w =
          0.5 * area * (1.0 / fine_grid->nodes->GetVolume(iPoint) + 1.0 / fine_grid->nodes->GetVolume(jPoint));
      if (w > wMax) {
        wMax = w;
        jStiffest = jPoint;
      }
      wMin = std::min(wMin, w);
    }

    /*--- A node with no neighbours has no aspect ratio to measure. ---*/
    if (jStiffest == NO_POINT) return false;
    if (((wMin > 0.0) ? wMax / wMin : su2double(1.0)) < MIN_AR) return false;

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

      /*--- A front claims its seed before the boundary pass runs, so the Euler wall curvature
       *    limit is applied here too, on the same terms. ---*/
      if ((config->GetMarker_All_KindBC(iMarker) == EULER_WALL) && !fine_grid->boundIsStraight[iMarker] &&
          (ComputeLocalCurvature(fine_grid, iPoint, iMarker) >= EULER_WALL_MAX_CURVATURE)) {
        seeds.nRefusedCurvature++;
        continue;
      }

      std::array<su2double, MAXNDIM> n0{};
      for (unsigned short d = 0; d < nDim; ++d) n0[d] = Normal[d];
      seeds.node.push_back(iPoint);
      seeds.normal.push_back(n0);
      taken[iPoint] = 1;
    }
  };

  /*--- Viscous walls always carry a stretched layer, so they seed unconditionally and first, which
   *    gives them the nodes where a wall meets another boundary. ---*/
  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++)
    if (isWall(config->GetMarker_All_KindBC(iMarker))) seedMarker(iMarker, false);

  /*--- Non-wall boundaries that still carry a layer normal to themselves, counted per
   *    configuration-file marker because local marker indices differ between ranks. ---*/
  const auto nMarkerCfg = config->GetnMarker_CfgFile();
  vector<unsigned long> nValid(nMarkerCfg, 0), nQualified(nMarkerCfg, 0);

  auto canSeed = [&](unsigned short iMarker) {
    const auto bc = config->GetMarker_All_KindBC(iMarker);
    /*--- Periodic boundaries are left out, they have their own matching. ---*/
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

  /*--- A marker is generally split over several ranks, so the verdict is taken on all of it. Every
   *    rank reaches these collectives. ---*/
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


  return seeds;
}

vector<vector<unsigned long>> CMultiGridGeometry::BuildFrontPatches(const CFrontSeeds& seeds,
                                                                    const CGeometry* fine_grid, const CConfig* config,
                                                                    const vector<char>& mixedBC) const {
  /*--- Repeated pairwise matching groups the seeds into connected patches of 1 to max_group seeds,
   *    smaller wherever no partner was found. A patch is any connected shape the boundary gives:
   *    a square, a strip, a triangle or a single node. What follows keys off the patch size and
   *    the layer isomorphism check, never off an assumed shape. ---*/
  const auto nSeeds = seeds.node.size();
  const unsigned long max_group = (nDim == 2) ? 2 : 4;

  /*--- Marker-set class of each seed. Seeds may only be matched when these agree. ---*/
  const auto nMarkerFine = fine_grid->GetnMarker();
  vector<std::pair<unsigned long, unsigned short>> seedMarker;
  for (unsigned long si = 0; si < nSeeds; ++si)
    for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++)
      if ((config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
          (fine_grid->nodes->GetVertex(seeds.node[si], iMarker) != -1))
        seedMarker.push_back({si, static_cast<unsigned short>(iMarker)});
  const auto sig = MarkerSetClasses(nSeeds, seedMarker);

  /*--- Seed-to-seed adjacency, inherited from the boundary nodes' mesh connectivity. ---*/
  vector<long> seedOfNode(fine_grid->GetnPoint(), -1);
  for (unsigned long si = 0; si < nSeeds; ++si) seedOfNode[seeds.node[si]] = static_cast<long>(si);

  vector<vector<unsigned long>> adj(nSeeds);
  for (unsigned long si = 0; si < nSeeds; ++si)
    for (auto jPoint : fine_grid->nodes->GetPoints(seeds.node[si])) {
      const auto sj = seedOfNode[jPoint];
      if ((sj >= 0) && (static_cast<unsigned long>(sj) != si)) adj[si].push_back(static_cast<unsigned long>(sj));
    }

  /*--- Global point index of each seed, used below as the partitioning-independent sort key. ---*/
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

  const unsigned nRounds = (max_group <= 2) ? 1 : 2;

  /*--- One admissible merge of two groups, weighted by how many seed-to-seed adjacencies they
   *    share: 2 for a group lying alongside, 1 for one continuing in the same direction. ---*/
  struct CMerge {
    unsigned long g, h;       /*!< \brief The two groups, g < h. */
    unsigned long weight;     /*!< \brief Shared adjacencies: 2 makes a square, 1 makes a strip. */
    unsigned long size;       /*!< \brief Seeds the merged group would hold. */
    unsigned long keyG, keyH; /*!< \brief Their global-index keys, the deterministic tie-break. */
  };

  vector<CMerge> merges;
  vector<unsigned long> touched, nShared(nSeeds, 0);
  vector<unsigned long> gkey;
  vector<char> consumed;

  for (unsigned round = 0; round < nRounds; ++round) {
    const auto nGroups = groups.size();

    /*--- Sort key of each group: the smallest global point index it holds. ---*/
    gkey.assign(nGroups, std::numeric_limits<unsigned long>::max());
    for (unsigned long g = 0; g < nGroups; ++g)
      for (auto si : groups[g]) gkey[g] = std::min(gkey[g], sgkey[si]);

    /*--- Every merge this round could make, counted once per unordered pair by taking only
     *    h > g. ---*/
    merges.clear();
    for (unsigned long g = 0; g < nGroups; ++g) {
      /*--- A node where two different boundary conditions meet stays a patch of its own, so both
       *    sides of a merge are tested for it. ---*/
      if (mixedBC[seeds.node[groups[g].front()]]) continue;
      touched.clear();
      for (auto si : groups[g])
        for (auto sj : adj[si]) {
          const auto h = groupOf[sj];
          if (h <= g) continue;
          if (mixedBC[seeds.node[groups[h].front()]]) continue;
          if (groups[g].size() + groups[h].size() > max_group) continue;
          if (sig[groups[h].front()] != sig[groups[g].front()]) continue;
          if (nShared[h]++ == 0) touched.push_back(h);
        }
      for (auto h : touched) {
        merges.push_back({g, h, nShared[h], static_cast<unsigned long>(groups[g].size() + groups[h].size()),
                          gkey[g], gkey[h]});
        nShared[h] = 0;
      }
    }

    /*--- Best merges first over all groups at once, so every square is considered before the first
     *    strip. At equal weight the bigger patch wins, to fill max_group before settling for less. ---*/
    std::sort(merges.begin(), merges.end(), [](const CMerge& a, const CMerge& b) {
      if (a.weight != b.weight) return a.weight > b.weight;
      if (a.size != b.size) return a.size > b.size;
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

string CMultiGridGeometry::PaveAdvancingFronts(unsigned long& Index_CoarseCV, const CGeometry* fine_grid,
                                                    const CConfig* config, unsigned short iMesh,
                                                    const vector<char>& mixedBC,
                                                    const vector<char>& onPhysBoundary,
                                                    vector<unsigned long>& neverGrewCV) {
  /*--- Paving by advancing fronts. Each boundary patch rises into the domain keeping its footprint,
   *    stopping at a boundary or where the next layer is not isomorphic to the current one. ---*/
  const auto starting_Index_CoarseCV = Index_CoarseCV;
  const auto nPointFine = fine_grid->GetnPoint();
  constexpr auto NO_POINT = std::numeric_limits<unsigned long>::max();
  const short int maxAgglomSize = (nDim == 2) ? 4 : 8;

  /*--- How nearly parallel a step must be to a boundary's normal to count as running into that
   *    boundary rather than along it. ---*/
  constexpr passivedouble BOUNDARY_ALIGN_DEG = 30.0;
  const su2double cos_boundary = cos(BOUNDARY_ALIGN_DEG * PI_NUMBER / 180.0);
  /*--- Weight of the new step direction when the front's marching direction is updated. ---*/
  constexpr passivedouble DIR_BLEND = 0.5;

  /*--- PHASE 1. SeedFrontNodes is collective and must be reached by every rank. ---*/
  const auto seeds = SeedFrontNodes(fine_grid, config);
  const auto patches = BuildFrontPatches(seeds, fine_grid, config, mixedBC);

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

  /*--- One advancing front. Fronts handed over from a neighbouring rank are appended while the
   *    rounds are running, so the loops below are bounded by fronts.size(). ---*/
  struct CFront {
    vector<unsigned long> nodes;          /*!< \brief Current footprint. */
    vector<unsigned long> pending;        /*!< \brief Nodes buffered for the coarse CV being built. */
    vector<unsigned long> handTo;         /*!< \brief Halo nodes the stack should continue onto. */
    vector<CStep> prop;                   /*!< \brief This round's proposed successors. */
    std::array<su2double, MAXNDIM> dir{}; /*!< \brief Marching direction. */
    unsigned long tag = 0;                /*!< \brief Rank-independent name, see TagOfSet. */
    unsigned long handTag = 0;            /*!< \brief Name the handed-over piece travels under, which
                                           *    after a split differs from tag. */
    unsigned long depth = 0;              /*!< \brief Layers laid. */
    unsigned long nBlock = 0;             /*!< \brief Fine layers the next coarse CV holds. */
    unsigned long pendingLayers = 0;      /*!< \brief Layers currently buffered. */
    char alive = 1;
    char failed = 0;
    char keepLocal = 0; /*!< \brief Handed over only part of its footprint, so it marches on here. */
    unsigned long seedCV = std::numeric_limits<unsigned long>::max(); /*!< \brief Coarse CV index of
        this front's first emitted layer, recorded so a front that never advances past it can be
        told apart, after the fact, from one that grew into a genuine stack. */
  };
  vector<CFront> fronts;

  auto addFront = [&](const vector<unsigned long>& layer, const std::array<su2double, MAXNDIM>& dir,
                      unsigned long frontTag, unsigned long block) {
    CFront F;
    F.nodes = F.pending = layer;
    F.dir = dir;
    F.tag = frontTag;
    F.nBlock = block;
    F.pendingLayers = 1;
    fronts.push_back(std::move(F));
    return fronts.size() - 1;
  };

  vector<char> claimed(nPointFine, 0);

  /*--- Best free neighbour of n to step onto, ranked by alignment with dir. Local and halo
   *    candidates are ranked separately, only the halo one can be handed over. ---*/
  struct CCandidate {
    unsigned long node = std::numeric_limits<unsigned long>::max();
    unsigned long halo = std::numeric_limits<unsigned long>::max();
    su2double dot = -2.0, len = 0.0, dir[MAXNDIM] = {0.0};
  };

  auto bestSuccessor = [&](unsigned long n, const su2double* marchDir) {
    CCandidate c;
    su2double haloDot = -2.0;

    for (auto jPoint : fine_grid->nodes->GetPoints(n)) {
      su2double vec[MAXNDIM] = {0.0};
      GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jPoint), fine_grid->nodes->GetCoord(n), vec);
      const su2double len = GeometryToolbox::Norm(nDim, vec);
      if (len <= 0.0) continue;
      for (unsigned short d = 0; d < nDim; ++d) vec[d] /= len;

      const su2double dot = GeometryToolbox::DotProduct(nDim, vec, marchDir);
      const bool admissible =
          !(onPhysBoundary[jPoint] && EntersBoundary(fine_grid, config, nDim, jPoint, vec, cos_boundary)) &&
          GeometricalCheck(jPoint, fine_grid, config);

      /*--- Halo parents are assigned by the owning rank through the MPI relay, so a halo node is
       *    only checked for admissibility here, never claimed. ---*/
      if (!fine_grid->nodes->GetDomain(jPoint)) {
        if (dot > haloDot && admissible) {
          haloDot = dot;
          c.halo = jPoint;
        }
        continue;
      }
      if (fine_grid->nodes->GetAgglomerate(jPoint) || claimed[jPoint] || !admissible) continue;

      if (dot > c.dot) {
        c.dot = dot;
        c.node = jPoint;
        c.len = len;
        for (unsigned short d = 0; d < nDim; ++d) c.dir[d] = vec[d];
      }
    }
    return c;
  };

  /*--- Bid table: an index per mesh point, the bids themselves in a compact vector. ---*/
  constexpr unsigned NOBID = std::numeric_limits<unsigned>::max();
  vector<unsigned> bidIdx(nPointFine, NOBID);
  vector<CStep> bids;
  vector<unsigned long> bidOwner;

  /*--- Scratch for the layer under construction. ---*/
  vector<unsigned long> newLayer;

  /*--- Summed over all ranks for the one-line report at the end. ---*/
  enum {
    P_STACKS,       /*!< \brief Fronts started. */
    P_LAYERS,       /*!< \brief Layers laid by all fronts. */
    P_COVERED,      /*!< \brief Fine nodes given a coarse parent by the paving. */
    P_RETIRED,      /*!< \brief Fronts that died on a failed layer rather than running out. */
    P_SEED_CURV,    /*!< \brief Euler wall seeds refused by the curvature limit. */
    P_HAND_SENT,    /*!< \brief Footprints offered across a partition interface. */
    P_HAND_SPLIT,   /*!< \brief ...of those, ones whose nodes are owned by more than one rank. */
    P_HAND_CONTEST, /*!< \brief Halo nodes two fronts reached for in the same round. */
    P_HAND_TAKEN,   /*!< \brief Inherited footprints adopted. */
    P_HAND_DROPPED, /*!< \brief Inherited footprints refused as not free or not connected. */
    P_PATCH1,       /*!< \brief Seed patches by width, 1 to 4. ---*/
    P_PATCH2,
    P_PATCH3,
    P_PATCH4,
    P_COUNT
  };
  unsigned long ct[P_COUNT] = {0};

  /*--- Patch widths, the footprint every front starts from. ---*/
  for (const auto& patch : patches) {
    if (patch.empty()) continue;
    ct[P_PATCH1 + std::min<size_t>(patch.size(), 4) - 1]++;
  }
  ct[P_SEED_CURV] = seeds.nRefusedCurvature;

  /*--- One footprint node arriving from a neighbouring rank, to be regrouped by tag. ---*/
  struct CInherited {
    unsigned long tag;      /*!< \brief The front it belongs to, the same name on both ranks. */
    unsigned long node;     /*!< \brief Local index of the node, owned by this rank. */
    su2double dir[MAXNDIM]; /*!< \brief The marching direction the stack arrives with. */
  };
  vector<CInherited> inherited;

  /*--- Where each halo node sits in this rank's SEND_RECEIVE receive lists, so a handover can be
   *    packed against the right vertex. ---*/
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

  auto markFail = [&](unsigned long f) { fronts[f].failed = 1; };

  /*--- Turn everything buffered for this front into one coarse control volume. ---*/
  auto emit = [&](unsigned long f) {
    if (fronts[f].pending.empty()) return;
    if (fronts[f].seedCV == std::numeric_limits<unsigned long>::max()) fronts[f].seedCV = Index_CoarseCV;
    nodes->SetChildren_CV(Index_CoarseCV, fronts[f].pending);
    for (auto p : fronts[f].pending) {
      fine_grid->nodes->SetParent_CV(p, Index_CoarseCV);
      if (fine_grid->nodes->GetAgglomerate_Indirect(p)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);
    }
    nodes->SetnChildren_CV(Index_CoarseCV, static_cast<unsigned short>(fronts[f].pending.size()));
    Index_CoarseCV++;
    ct[P_COVERED] += fronts[f].pending.size();

    fronts[f].pending.clear();
    fronts[f].pendingLayers = 0;
    fronts[f].nBlock = BlockFor(maxAgglomSize, fronts[f].nodes);
  };

  /*--- The boundary layer of every front, claimed before ordinary boundary agglomeration runs so
   *    that every layer above keeps the same footprint. ---*/
  for (const auto& patch : patches) {
    /*--- A one-node patch may seed a front and marches as a stack one node wide. ---*/
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

    /*--- The boundary layer is a coarse CV on its own, so only the first advance is a single
     *    layer. ---*/
    const auto f = addFront(layer0, n0, frontTag + 1, 1);
    for (auto p : layer0) {
      claimed[p] = 1;
    }
    ct[P_STACKS]++;
    ct[P_LAYERS]++;
    emit(f);
  }

  /*--- The handover carries a marching direction, but only to make discrete choices, so it is
   *    exchanged as a passive type and never taped. ---*/
  using CPassiveMPI = SelectMPIWrapper<passivedouble>::W;

  /*--- One SEND_RECEIVE marker pair per neighbour, with where its vertices sit in the flat
   *    exchange buffers. ---*/
  struct CHandoverPair {
    unsigned short markerS, markerR;
    int send_to, receive_from;
    unsigned long nVertexS, nVertexR, offS, offR;
  };
  vector<CHandoverPair> handPairs;
  unsigned long nSendTotal = 0, nRecvTotal = 0;
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (!((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)))
      continue;
    CHandoverPair hp;
    hp.markerS = iMarker;
    hp.markerR = iMarker + 1;
    hp.send_to = config->GetMarker_All_SendRecv(hp.markerS) - 1;
    hp.receive_from = abs(config->GetMarker_All_SendRecv(hp.markerR)) - 1;
    hp.nVertexS = fine_grid->nVertex[hp.markerS];
    hp.nVertexR = fine_grid->nVertex[hp.markerR];
    hp.offS = nSendTotal;
    hp.offR = nRecvTotal;
    nSendTotal += hp.nVertexS;
    nRecvTotal += hp.nVertexR;
    handPairs.push_back(hp);
  }

  /*--- Reused every round: the handovers bucketed by the receive marker they cross, and the
   *    exchange buffers packed against the halo vertices of every pair at once. ---*/
  vector<vector<std::pair<unsigned long, unsigned long>>> handByMarker(config->GetnMarker_All());
  vector<unsigned long> tagOut(nRecvTotal), tagIn(nSendTotal);
  vector<passivedouble> dirOut(nRecvTotal * nDim), dirIn(nSendTotal * nDim);
  vector<CPassiveMPI::Request> handReq(4 * handPairs.size());
  vector<CPassiveMPI::Status> handStat(4 * handPairs.size());

  for (unsigned long layer = 1;; ++layer) {
    /*--- Every rank runs the same number of rounds, each ends in a collective handover
     *    exchange. ---*/
    int aliveLocal = 0;
    for (unsigned long f = 0; f < fronts.size(); ++f) aliveLocal |= fronts[f].alive;
    int aliveGlobal = 0;
    SU2_MPI::Allreduce(&aliveLocal, &aliveGlobal, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
    if (aliveGlobal == 0) break;

    for (auto& F : fronts) F.failed = F.keepLocal = F.handTag = 0;
    for (const auto& b : bids) bidIdx[b.node] = NOBID;
    bids.clear();
    bidOwner.clear();

    /*--- (a) Every alive front proposes a successor for each of its nodes. A front that cannot fill
     *    a whole layer proposes nothing and retires this round. ---*/
    for (unsigned long f = 0; f < fronts.size(); ++f) {
      if (!fronts[f].alive) continue;
      fronts[f].prop.clear();
      fronts[f].handTo.clear();

      for (auto n : fronts[f].nodes) {
        const auto c = bestSuccessor(n, fronts[f].dir.data());

        /*--- Nothing free here but the stack continues across the interface. ---*/
        if ((c.node == NO_POINT) && (c.halo != NO_POINT)) {
          fronts[f].handTo.push_back(c.halo);
          continue;
        }
        /*--- No successor at all: a boundary, a partition, another front, or unusable mesh. ---*/
        if (c.node == NO_POINT) {
          markFail(f);
          break;
        }

        CStep s{c.node, n, fine_grid->nodes->GetGlobalIndex(n), c.dot, c.len, {}};
        for (unsigned short d = 0; d < nDim; ++d) s.dir[d] = c.dir[d];
        fronts[f].prop.push_back(s);
      }

      /*--- An interface can cut a footprint. All of it crossing hands the stack over intact, part of
       *    it crossing splits the footprint into a piece that stays and a piece handed across. ---*/
      if (fronts[f].failed) {
        fronts[f].prop.clear();
        fronts[f].handTo.clear();
      } else if (!fronts[f].handTo.empty()) {
        /*--- fronts[f].prop is built in the order of fronts[f].nodes, so this is the piece that
         *    stays, in the same order as the layer it proposes. ---*/
        vector<unsigned long> narrow;
        for (const auto& s : fronts[f].prop) narrow.push_back(s.from);

        /*--- A cut can leave the local piece disconnected, which is not a layer. Drop it and hand
         *    over the rest. ---*/
        if (!narrow.empty() && !IsConnectedLayer(fine_grid, narrow)) {
          narrow.clear();
        }

        fronts[f].handTag = TagOfSet(fine_grid, fronts[f].handTo);

        if (narrow.empty()) {
          fronts[f].prop.clear();
        } else {
          /*--- Close the coarse CV open on the wide footprint before narrowing, so that no CV holds
           *    two layers of different shape. ---*/
          fronts[f].nodes = narrow;
          emit(f);
          fronts[f].nBlock = BlockFor(maxAgglomSize, fronts[f].nodes);
          fronts[f].tag = TagOfSet(fine_grid, fronts[f].nodes);
          fronts[f].keepLocal = 1;
        }
      }
    }

    /*--- (b) Place every bid; nothing is granted until (c) reads the settled table. ---*/
    auto better = [](const CStep& a, const CStep& b) {
      if (a.score != b.score) return a.score > b.score;
      if (a.dist != b.dist) return a.dist < b.dist;
      return a.key < b.key;
    };

    for (unsigned long f = 0; f < fronts.size(); ++f) {
      if (!fronts[f].alive || fronts[f].failed) continue;
      for (const auto& s : fronts[f].prop) {
        /*--- A front that has lost a bid retires and places no more of its layer. ---*/
        if (fronts[f].failed) break;

        if (bidIdx[s.node] == NOBID) {
          bidIdx[s.node] = static_cast<unsigned>(bids.size());
          bids.push_back(s);
          bidOwner.push_back(f);
          continue;
        }

        const auto k = bidIdx[s.node];
        const auto g = bidOwner[k];
        /*--- Two nodes of the same front reaching for one successor is a pinch, the layer would come
         *    out narrower than the front. ---*/
        if (better(s, bids[k])) {
          markFail(g);
          bids[k] = s;
          bidOwner[k] = f;
        } else {
          markFail(f);
        }

        /*--- A head-on meeting stops both fronts, a glancing contact only costs the loser. ---*/
        if ((g != f) && (GeometryToolbox::DotProduct(nDim, fronts[f].dir.data(), fronts[g].dir.data()) < 0.0)) {
          markFail(f);
          markFail(g);
        }
      }
    }

    /*--- (c) All-or-nothing acceptance: a front takes the whole layer or none of it and retires. A
     *    bid only becomes a claim here. ---*/
    for (unsigned long f = 0; f < fronts.size(); ++f) {
      if (!fronts[f].alive) continue;

      newLayer.clear();
      if (!fronts[f].failed) {
        /*--- Built in proposal order, so newLayer[k] is the successor of nodes[k]. ---*/
        for (const auto& s : fronts[f].prop) {
          const auto k = bidIdx[s.node];
          if ((k != NOBID) && (bidOwner[k] == f)) newLayer.push_back(s.node);
        }
        /*--- Every bid of a front that was not marked failed must have been granted. ---*/
        if ((newLayer.size() != fronts[f].nodes.size()) || !LayerIsIsomorphic(fine_grid, fronts[f].nodes, newLayer))
          markFail(f);
      }

      if (fronts[f].failed) {
        fronts[f].alive = 0;
        ct[P_RETIRED]++;
        emit(f);
        continue;
      }

      /*--- Turn the marching direction towards the mean of the steps just taken. ---*/
      su2double mean[MAXNDIM] = {0.0};
      for (const auto& s : fronts[f].prop)
        for (unsigned short d = 0; d < nDim; ++d) mean[d] += s.dir[d];
      const su2double meanNrm = GeometryToolbox::Norm(nDim, mean);
      if (meanNrm > 0.0) {
        su2double blended[MAXNDIM] = {0.0};
        for (unsigned short d = 0; d < nDim; ++d)
          blended[d] = (1.0 - DIR_BLEND) * fronts[f].dir[d] + DIR_BLEND * mean[d] / meanNrm;
        const su2double bNrm = GeometryToolbox::Norm(nDim, blended);
        if (bNrm > 0.0)
          for (unsigned short d = 0; d < nDim; ++d) fronts[f].dir[d] = blended[d] / bNrm;
      }

      for (auto p : newLayer) {
        claimed[p] = 1;
      }
      fronts[f].nodes = std::move(newLayer);
      fronts[f].depth++;
      ct[P_LAYERS]++;

      fronts[f].pending.insert(fronts[f].pending.end(), fronts[f].nodes.begin(), fronts[f].nodes.end());
      fronts[f].pendingLayers++;
      if (fronts[f].pendingLayers >= fronts[f].nBlock) emit(f);
    }

    /*--- (d) Hand stacks across partition interfaces. The footprint is sent to the owning rank,
     *    packed against the receive marker and sent to the rank that marker receives from. ---*/
    for (auto& bucket : handByMarker) bucket.clear();
    for (unsigned long f = 0; f < fronts.size(); ++f) {
      int firstMarker = -1;
      bool split = false;
      for (auto p : fronts[f].handTo) {
        if (haloMarker[p] < 0) continue;
        if (firstMarker < 0) firstMarker = haloMarker[p];
        else if (haloMarker[p] != firstMarker) split = true;
        handByMarker[haloMarker[p]].push_back({f, p});
      }
      /*--- A footprint owned by more than one rank cannot travel whole. ---*/
      if (firstMarker >= 0) {
        ct[P_HAND_SENT]++;
        ct[P_HAND_SPLIT] += split;
      }
    }

    /*--- Pack every pair, then exchange them all at once so a round costs one wait, not one
     *    round trip per neighbour. Tag and direction go separately, there is no byte type to
     *    send a struct with. ---*/
    std::fill(tagOut.begin(), tagOut.end(), 0ul);
    std::fill(dirOut.begin(), dirOut.end(), passivedouble(0.0));

    for (const auto& hp : handPairs)
      for (const auto& hand : handByMarker[hp.markerR]) {
        const auto& F = fronts[hand.first];
        const auto v = hp.offR + haloVertex[hand.second];
        /*--- Two fronts reaching for one node: the lower tag takes it. ---*/
        if (tagOut[v] != 0) ct[P_HAND_CONTEST]++;
        if ((tagOut[v] != 0) && (tagOut[v] <= F.handTag)) continue;
        tagOut[v] = F.handTag;
        for (unsigned short d = 0; d < nDim; ++d) dirOut[v * nDim + d] = SU2_TYPE::GetValue(F.dir[d]);
      }

    unsigned long nReq = 0;
    for (const auto& hp : handPairs) {
      CPassiveMPI::Irecv(&tagIn[hp.offS], hp.nVertexS, MPI_UNSIGNED_LONG, hp.send_to, 2, CPassiveMPI::GetComm(),
                         &handReq[nReq++]);
      CPassiveMPI::Irecv(&dirIn[hp.offS * nDim], hp.nVertexS * nDim, MPI_DOUBLE, hp.send_to, 3,
                         CPassiveMPI::GetComm(), &handReq[nReq++]);
      CPassiveMPI::Isend(&tagOut[hp.offR], hp.nVertexR, MPI_UNSIGNED_LONG, hp.receive_from, 2,
                         CPassiveMPI::GetComm(), &handReq[nReq++]);
      CPassiveMPI::Isend(&dirOut[hp.offR * nDim], hp.nVertexR * nDim, MPI_DOUBLE, hp.receive_from, 3,
                         CPassiveMPI::GetComm(), &handReq[nReq++]);
    }
    if (nReq > 0) CPassiveMPI::Waitall(static_cast<int>(nReq), handReq.data(), handStat.data());

    for (const auto& hp : handPairs)
      for (auto iVertex = 0ul; iVertex < hp.nVertexS; iVertex++) {
        const auto v = hp.offS + iVertex;
        if (tagIn[v] == 0) continue;
        inherited.push_back({tagIn[v], fine_grid->vertex[hp.markerS][iVertex]->GetNode(), {}});
        for (unsigned short d = 0; d < nDim; ++d) inherited.back().dir[d] = dirIn[v * nDim + d];
      }

    /*--- A front that handed its whole footprint over is finished here, one that handed over only a
     *    piece keeps marching on what was left. ---*/
    for (unsigned long f = 0; f < fronts.size(); ++f) {
      if (fronts[f].handTo.empty()) continue;
      fronts[f].handTo.clear();
      if (fronts[f].keepLocal) continue;
      fronts[f].alive = 0;
      emit(f);
    }

    /*--- Adopt what the neighbours sent, tags ascending so arrival order cannot change the outcome.
     *    A footprint whose nodes are not all free is dropped. ---*/
    std::sort(inherited.begin(), inherited.end(),
              [](const CInherited& a, const CInherited& b) { return a.tag < b.tag; });

    for (size_t i = 0; i < inherited.size();) {
      size_t j = i;
      while ((j < inherited.size()) && (inherited[j].tag == inherited[i].tag)) ++j;

      vector<unsigned long> layer0;
      bool ok = true;
      for (size_t k = i; k < j; ++k) {
        const auto p = inherited[k].node;
        /*--- One node can be handed over by two neighbours under the same tag. ---*/
        if (std::find(layer0.begin(), layer0.end(), p) != layer0.end()) continue;
        if (claimed[p] || fine_grid->nodes->GetAgglomerate(p) || !GeometricalCheck(p, fine_grid, config)) ok = false;
        layer0.push_back(p);
      }
      /*--- Whatever arrived must be free and form one connected layer. ---*/
      if (ok && !IsConnectedLayer(fine_grid, layer0)) ok = false;

      if (ok) {
        std::array<su2double, MAXNDIM> d0{};
        for (unsigned short d = 0; d < nDim; ++d) d0[d] = inherited[i].dir[d];
        /*--- An inherited layer is an interior one, so it opens an ordinary two-deep coarse CV. ---*/
        const auto nf = addFront(layer0, d0, inherited[i].tag, BlockFor(maxAgglomSize, layer0));
        for (auto p : layer0) {
          claimed[p] = 1;
        }
        ct[P_LAYERS]++;
        ct[P_HAND_TAKEN]++;
        if (fronts[nf].pendingLayers >= fronts[nf].nBlock) emit(nf);
      } else {
        ct[P_HAND_DROPPED]++;
      }
      i = j;
    }
    inherited.clear();
  }

  /*--- Emit whatever is still buffered, so no node is left without a parent index. This also
   *    catches an adopted front (block size 2) whose very first round after adoption fails: it
   *    dies at depth 0 without ever having called emit() on its own, so its seed CV is only
   *    assigned here. ---*/
  for (unsigned long f = 0; f < fronts.size(); ++f) emit(f);

  /*--- A front that never advanced past its seed did not build a stack at all: it is
   *    indistinguishable from ordinary boundary agglomeration and gains nothing from the
   *    protection below, which exists to keep a genuine line intact. ---*/
  for (const auto& F : fronts)
    if ((F.depth == 0) && (F.seedCV != std::numeric_limits<unsigned long>::max())) neverGrewCV.push_back(F.seedCV);

  /*--- A rank with no fronts leaves dmin at its sentinel, keeping it out of the MPI_MIN. ---*/
  unsigned long dmin = std::numeric_limits<unsigned long>::max(), dmax = 0;
  for (unsigned long f = 0; f < fronts.size(); ++f) {
    if (fronts[f].nodes.empty()) continue;
    dmin = std::min(dmin, fronts[f].depth);
    dmax = std::max(dmax, fronts[f].depth);
  }

  /*--- Every rank must reach these. ---*/
  unsigned long tot[P_COUNT] = {0}, pairTot[2] = {0}, depthMin = 0, depthMax = 0;
  unsigned long pair[2] = {Index_CoarseCV - starting_Index_CoarseCV, seeds.node.size()};
  SU2_MPI::Allreduce(ct, tot, P_COUNT, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(pair, pairTot, 2, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&dmin, &depthMin, 1, MPI_UNSIGNED_LONG, MPI_MIN, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&dmax, &depthMax, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  if (depthMin == std::numeric_limits<unsigned long>::max()) depthMin = 0;

  if (rank != MASTER_NODE) return {};

  stringstream out;
  out << "  MG level " << iMesh << " paving: " << tot[P_STACKS] << " fronts from " << pairTot[1] << " seeds, "
      << pairTot[0] << " CVs covering " << tot[P_COVERED] << " nodes in " << tot[P_LAYERS] << " layers, depth "
      << depthMin << " to " << depthMax << "\n";
  out << "    patches 1/2/3/4 wide: " << tot[P_PATCH1] << "/" << tot[P_PATCH2] << "/" << tot[P_PATCH3] << "/"
      << tot[P_PATCH4] << ", " << tot[P_SEED_CURV] << " seeds refused on curvature, " << tot[P_RETIRED]
      << " fronts retired early\n";
  if (tot[P_HAND_SENT] > 0)
    out << "    handovers: " << tot[P_HAND_SENT] << " offered (" << tot[P_HAND_SPLIT] << " split across ranks, "
        << tot[P_HAND_CONTEST] << " nodes contested), " << tot[P_HAND_TAKEN] << " adopted, " << tot[P_HAND_DROPPED]
        << " dropped\n";
  return out.str();
}
