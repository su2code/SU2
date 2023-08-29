/*!
 * \file CSlidingMesh.cpp
 * \brief Implementation of sliding mesh interpolation.
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

#include "../../include/interface_interpolation/CSlidingMesh.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

CSlidingMesh::CSlidingMesh(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                           unsigned int jZone)
    : CInterpolator(geometry_container, config, iZone, jZone) {
  SetTransferCoeff(config);
}

void CSlidingMesh::SetTransferCoeff(const CConfig* const* config) {
  /* 0 - Variable declaration */

  /* --- General variables --- */

  bool check;

  unsigned short iDim;

  unsigned long ii, jj, *uptr;
  unsigned long vPoint;
  unsigned long iEdgeVisited, nEdgeVisited, iNodeVisited;
  unsigned long nAlreadyVisited, nToVisit, StartVisited;

  unsigned long *alreadyVisitedDonor, *ToVisit, *tmpVect;
  unsigned long *storeProc, *tmp_storeProc;

  su2double dTMP;
  su2double *Coeff_Vect, *tmp_Coeff_Vect;

  /* --- Geometrical variables --- */

  su2double *Coord_i, *Coord_j, dist, mindist, *Normal;
  su2double Area, Area_old, tmp_Area;
  su2double LineIntersectionLength, *Direction, length;

  /* --- Markers Variables --- */

  unsigned short iMarkerInt, nMarkerInt;

  unsigned long iVertex, nVertexTarget;

  int markDonor, markTarget;

  /* --- Target variables --- */

  unsigned long target_iPoint, jVertexTarget;
  unsigned long nEdges_target, nNode_target;

  su2vector<unsigned long> Target_nLinkedNodes;
  su2vector<unsigned long> Target_StartLinkedNodes;
  unsigned long* target_segment;
  su2vector<unsigned long> Target_LinkedNodes;
  su2vector<unsigned long> Target_GlobalPoint, Donor_GlobalPoint;

  su2double *target_iMidEdge_point, *target_jMidEdge_point, **target_element;
  su2activematrix TargetPoint_Coord;

  /* --- Donor variables --- */

  unsigned long donor_StartIndex, donor_forward_point, donor_backward_point, donor_iPoint, donor_OldiPoint;
  unsigned long nEdges_donor, nNode_donor, nGlobalVertex_Donor;

  unsigned long nDonorPoints, iDonor;
  unsigned long *Donor_Vect, *tmp_Donor_Vect;
  su2vector<unsigned long> Donor_nLinkedNodes;
  su2vector<unsigned long> Donor_StartLinkedNodes;
  su2vector<unsigned long> Donor_LinkedNodes;
  su2vector<unsigned long> Donor_Proc;

  su2double *donor_iMidEdge_point, *donor_jMidEdge_point;
  su2double** donor_element;
  su2activematrix DonorPoint_Coord;

  targetVertices.resize(config[targetZone]->GetnMarker_All());

  /* 1 - Variable pre-processing */

  const unsigned short nDim = donor_geometry->GetnDim();

  /*--- Setting up auxiliary vectors ---*/

  Donor_Vect = nullptr;
  Coeff_Vect = nullptr;
  storeProc = nullptr;

  tmp_Donor_Vect = nullptr;
  tmp_Coeff_Vect = nullptr;
  tmp_storeProc = nullptr;

  Normal = new su2double[nDim];
  Direction = new su2double[nDim];

  /* 2 - Find boundary tag between touching grids */

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt = (int)(config[donorZone]->GetMarker_n_ZoneInterface()) / 2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt = 0; iMarkerInt < nMarkerInt; iMarkerInt++) {
    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor = config[donorZone]->FindInterfaceMarker(iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = config[targetZone]->FindInterfaceMarker(iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    nVertexTarget = 0;
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*
    3 -Reconstruct the boundaries from parallel partitioning
    */

    /*--- Target boundary ---*/
    ReconstructBoundary(targetZone, markTarget);

    nGlobalVertex_Target = nGlobalVertex;

    TargetPoint_Coord = Buffer_Receive_Coord;
    Target_GlobalPoint = Buffer_Receive_GlobalPoint;
    Target_nLinkedNodes = Buffer_Receive_nLinkedNodes;
    Target_StartLinkedNodes = Buffer_Receive_StartLinkedNodes;
    Target_LinkedNodes = Buffer_Receive_LinkedNodes;

    /*--- Donor boundary ---*/
    ReconstructBoundary(donorZone, markDonor);

    nGlobalVertex_Donor = nGlobalVertex;

    DonorPoint_Coord = Buffer_Receive_Coord;
    Donor_GlobalPoint = Buffer_Receive_GlobalPoint;
    Donor_nLinkedNodes = Buffer_Receive_nLinkedNodes;
    Donor_StartLinkedNodes = Buffer_Receive_StartLinkedNodes;
    Donor_LinkedNodes = Buffer_Receive_LinkedNodes;
    Donor_Proc = Buffer_Receive_Proc;

    /*--- Starts building the supermesh layer (2D or 3D) ---*/
    /* - For each target node, it first finds the closest donor point
     * - Then it creates the supermesh in the close proximity of the target point:
     * - Starting from the closest donor node, it expands the supermesh by including
     * donor elements neighboring the initial one, until the overall target area is fully covered.
     */
    if (nVertexTarget) targetVertices[markTarget].resize(nVertexTarget);

    if (nDim == 2) {
      target_iMidEdge_point = new su2double[nDim];
      target_jMidEdge_point = new su2double[nDim];

      donor_iMidEdge_point = new su2double[nDim];
      donor_jMidEdge_point = new su2double[nDim];

      /*--- Starts with supermesh reconstruction ---*/

      target_segment = new unsigned long[2];

      for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {
        nDonorPoints = 0;

        /*--- Stores coordinates of the target node ---*/

        target_iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();

        if (target_geometry->nodes->GetDomain(target_iPoint)) {
          Coord_i = target_geometry->nodes->GetCoord(target_iPoint);

          /*--- Brute force to find the closest donor_node ---*/

          mindist = 1E6;
          donor_StartIndex = 0;

          for (donor_iPoint = 0; donor_iPoint < nGlobalVertex_Donor; donor_iPoint++) {
            Coord_j = DonorPoint_Coord[donor_iPoint];

            dist = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

            if (dist < mindist) {
              mindist = dist;
              donor_StartIndex = donor_iPoint;
            }

            if (dist == 0.0) {
              donor_StartIndex = donor_iPoint;
              break;
            }
          }

          donor_iPoint = donor_StartIndex;
          donor_OldiPoint = donor_iPoint;

          /*--- Contruct information regarding the target cell ---*/

          auto dPoint = target_geometry->nodes->GetGlobalIndex(target_iPoint);
          for (jVertexTarget = 0; jVertexTarget < nGlobalVertex_Target; jVertexTarget++)
            if (dPoint == Target_GlobalPoint[jVertexTarget]) break;

          if (Target_nLinkedNodes[jVertexTarget] == 1) {
            target_segment[0] = Target_LinkedNodes[Target_StartLinkedNodes[jVertexTarget]];
            target_segment[1] = jVertexTarget;
          } else {
            target_segment[0] = Target_LinkedNodes[Target_StartLinkedNodes[jVertexTarget]];
            target_segment[1] = Target_LinkedNodes[Target_StartLinkedNodes[jVertexTarget] + 1];
          }

          dTMP = 0;
          for (iDim = 0; iDim < nDim; iDim++) {
            target_iMidEdge_point[iDim] =
                (TargetPoint_Coord(target_segment[0], iDim) + target_geometry->nodes->GetCoord(target_iPoint, iDim)) /
                2.;
            target_jMidEdge_point[iDim] =
                (TargetPoint_Coord(target_segment[1], iDim) + target_geometry->nodes->GetCoord(target_iPoint, iDim)) /
                2.;

            Direction[iDim] = target_jMidEdge_point[iDim] - target_iMidEdge_point[iDim];
            dTMP += Direction[iDim] * Direction[iDim];
          }

          dTMP = sqrt(dTMP);
          for (iDim = 0; iDim < nDim; iDim++) Direction[iDim] /= dTMP;

          length = GeometryToolbox::Distance(nDim, target_iMidEdge_point, target_jMidEdge_point);

          check = false;

          /*--- Proceeds along the forward direction (depending on which connected boundary node is found first) ---*/

          while (!check) {
            /*--- Proceeds until the value of the intersection area is null ---*/

            if (Donor_nLinkedNodes[donor_iPoint] == 1) {
              donor_forward_point = Donor_LinkedNodes[Donor_StartLinkedNodes[donor_iPoint]];
              donor_backward_point = donor_iPoint;
            } else {
              uptr = &Donor_LinkedNodes[Donor_StartLinkedNodes[donor_iPoint]];

              if (donor_OldiPoint != uptr[0]) {
                donor_forward_point = uptr[0];
                donor_backward_point = uptr[1];
              } else {
                donor_forward_point = uptr[1];
                donor_backward_point = uptr[0];
              }
            }

            if (donor_iPoint >= nGlobalVertex_Donor) {
              check = true;
              continue;
            }

            for (iDim = 0; iDim < nDim; iDim++) {
              donor_iMidEdge_point[iDim] =
                  (DonorPoint_Coord(donor_forward_point, iDim) + DonorPoint_Coord(donor_iPoint, iDim)) / 2.;
              donor_jMidEdge_point[iDim] =
                  (DonorPoint_Coord(donor_backward_point, iDim) + DonorPoint_Coord(donor_iPoint, iDim)) / 2.;
            }

            LineIntersectionLength =
                ComputeLineIntersectionLength(nDim, target_iMidEdge_point, target_jMidEdge_point, donor_iMidEdge_point,
                                              donor_jMidEdge_point, Direction);

            if (LineIntersectionLength == 0.0) {
              check = true;
              continue;
            }

            /*--- In case the element intersects the target cell, update the auxiliary communication data structure
             * ---*/

            tmp_Coeff_Vect = new su2double[nDonorPoints + 1];
            tmp_Donor_Vect = new unsigned long[nDonorPoints + 1];
            tmp_storeProc = new unsigned long[nDonorPoints + 1];

            for (iDonor = 0; iDonor < nDonorPoints; iDonor++) {
              tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
              tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
              tmp_storeProc[iDonor] = storeProc[iDonor];
            }

            tmp_Donor_Vect[nDonorPoints] = donor_iPoint;
            tmp_Coeff_Vect[nDonorPoints] = LineIntersectionLength / length;
            tmp_storeProc[nDonorPoints] = Donor_Proc[donor_iPoint];

            delete[] Donor_Vect;
            delete[] Coeff_Vect;
            delete[] storeProc;

            Donor_Vect = tmp_Donor_Vect;
            Coeff_Vect = tmp_Coeff_Vect;
            storeProc = tmp_storeProc;

            donor_OldiPoint = donor_iPoint;
            donor_iPoint = donor_forward_point;

            nDonorPoints++;
          }

          if (Donor_nLinkedNodes[donor_StartIndex] == 2) {
            check = false;

            uptr = &Donor_LinkedNodes[Donor_StartLinkedNodes[donor_StartIndex]];

            donor_iPoint = uptr[1];
            donor_OldiPoint = donor_StartIndex;
          } else
            check = true;

          /*--- Proceeds along the backward direction (depending on which connected boundary node is found first) ---*/

          while (!check) {
            /*--- Proceeds until the value of the intersection length is null ---*/
            if (Donor_nLinkedNodes[donor_iPoint] == 1) {
              donor_forward_point = donor_OldiPoint;
              donor_backward_point = donor_iPoint;
            } else {
              uptr = &Donor_LinkedNodes[Donor_StartLinkedNodes[donor_iPoint]];

              if (donor_OldiPoint != uptr[0]) {
                donor_forward_point = uptr[0];
                donor_backward_point = uptr[1];
              } else {
                donor_forward_point = uptr[1];
                donor_backward_point = uptr[0];
              }
            }

            if (donor_iPoint >= nGlobalVertex_Donor) {
              check = true;
              continue;
            }

            for (iDim = 0; iDim < nDim; iDim++) {
              donor_iMidEdge_point[iDim] =
                  (DonorPoint_Coord(donor_forward_point, iDim) + DonorPoint_Coord(donor_iPoint, iDim)) / 2.;
              donor_jMidEdge_point[iDim] =
                  (DonorPoint_Coord(donor_backward_point, iDim) + DonorPoint_Coord(donor_iPoint, iDim)) / 2.;
            }

            LineIntersectionLength =
                ComputeLineIntersectionLength(nDim, target_iMidEdge_point, target_jMidEdge_point, donor_iMidEdge_point,
                                              donor_jMidEdge_point, Direction);

            if (LineIntersectionLength == 0.0) {
              check = true;
              continue;
            }

            /*--- In case the element intersects the target cell, update the auxiliary communication data structure
             * ---*/

            tmp_Coeff_Vect = new su2double[nDonorPoints + 1];
            tmp_Donor_Vect = new unsigned long[nDonorPoints + 1];
            tmp_storeProc = new unsigned long[nDonorPoints + 1];

            for (iDonor = 0; iDonor < nDonorPoints; iDonor++) {
              tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
              tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
              tmp_storeProc[iDonor] = storeProc[iDonor];
            }

            tmp_Coeff_Vect[nDonorPoints] = LineIntersectionLength / length;
            tmp_Donor_Vect[nDonorPoints] = donor_iPoint;
            tmp_storeProc[nDonorPoints] = Donor_Proc[donor_iPoint];

            delete[] Donor_Vect;
            delete[] Coeff_Vect;
            delete[] storeProc;

            Donor_Vect = tmp_Donor_Vect;
            Coeff_Vect = tmp_Coeff_Vect;
            storeProc = tmp_storeProc;

            donor_OldiPoint = donor_iPoint;
            donor_iPoint = donor_forward_point;

            nDonorPoints++;
          }

          /*--- Set the communication data structure and copy data from the auxiliary vectors ---*/

          targetVertices[markTarget][iVertex].resize(nDonorPoints);

          for (iDonor = 0; iDonor < nDonorPoints; iDonor++) {
            targetVertices[markTarget][iVertex].coefficient[iDonor] = Coeff_Vect[iDonor];
            targetVertices[markTarget][iVertex].globalPoint[iDonor] = Donor_GlobalPoint[Donor_Vect[iDonor]];
            targetVertices[markTarget][iVertex].processor[iDonor] = storeProc[iDonor];
          }
        }
      }

      delete[] target_segment;

      delete[] target_iMidEdge_point;
      delete[] target_jMidEdge_point;

      delete[] donor_iMidEdge_point;
      delete[] donor_jMidEdge_point;
    } else {
      /* --- 3D geometry, creates a superficial super-mesh --- */

      for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {
        nDonorPoints = 0;

        /*--- Stores coordinates of the target node ---*/

        target_iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();

        if (!target_geometry->nodes->GetDomain(target_iPoint)) continue;

        Coord_i = target_geometry->nodes->GetCoord(target_iPoint);

        target_geometry->vertex[markTarget][iVertex]->GetNormal(Normal);

        /*--- The value of Area computed here includes also portion of boundary belonging to different marker ---*/
        Area = GeometryToolbox::Norm(nDim, Normal);

        for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] /= Area;

        for (iDim = 0; iDim < nDim; iDim++) Coord_i[iDim] = target_geometry->nodes->GetCoord(target_iPoint, iDim);

        auto dPoint = target_geometry->nodes->GetGlobalIndex(target_iPoint);
        for (target_iPoint = 0; target_iPoint < nGlobalVertex_Target; target_iPoint++) {
          if (dPoint == Target_GlobalPoint[target_iPoint]) break;
        }

        /*--- Build local surface dual mesh for target element ---*/

        nEdges_target = Target_nLinkedNodes[target_iPoint];

        nNode_target = 2 * (nEdges_target + 1);

        target_element = new su2double*[nNode_target];
        for (ii = 0; ii < nNode_target; ii++) target_element[ii] = new su2double[nDim];

        nNode_target = Build_3D_surface_element(Target_LinkedNodes, Target_StartLinkedNodes, Target_nLinkedNodes,
                                                TargetPoint_Coord, target_iPoint, target_element);

        /*--- Brute force to find the closest donor_node ---*/

        mindist = 1E6;
        donor_StartIndex = 0;

        for (donor_iPoint = 0; donor_iPoint < nGlobalVertex_Donor; donor_iPoint++) {
          Coord_j = DonorPoint_Coord[donor_iPoint];

          dist = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

          if (dist < mindist) {
            mindist = dist;
            donor_StartIndex = donor_iPoint;
          }

          if (dist == 0.0) {
            donor_StartIndex = donor_iPoint;
            break;
          }
        }

        donor_iPoint = donor_StartIndex;

        nEdges_donor = Donor_nLinkedNodes[donor_iPoint];

        donor_element = new su2double*[2 * nEdges_donor + 2];
        for (ii = 0; ii < 2 * nEdges_donor + 2; ii++) donor_element[ii] = new su2double[nDim];

        nNode_donor = Build_3D_surface_element(Donor_LinkedNodes, Donor_StartLinkedNodes, Donor_nLinkedNodes,
                                               DonorPoint_Coord, donor_iPoint, donor_element);

        Area = 0;
        for (ii = 1; ii < nNode_target - 1; ii++) {
          for (jj = 1; jj < nNode_donor - 1; jj++) {
            Area += Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii + 1],
                                                  donor_element[0], donor_element[jj], donor_element[jj + 1], Normal);
          }
        }

        for (ii = 0; ii < 2 * nEdges_donor + 2; ii++) delete[] donor_element[ii];
        delete[] donor_element;

        nDonorPoints = 1;

        /*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/

        Coeff_Vect = new su2double[nDonorPoints];
        Donor_Vect = new unsigned long[nDonorPoints];
        storeProc = new unsigned long[nDonorPoints];

        Coeff_Vect[0] = Area;
        Donor_Vect[0] = donor_iPoint;
        storeProc[0] = Donor_Proc[donor_iPoint];

        alreadyVisitedDonor = new unsigned long[1];

        alreadyVisitedDonor[0] = donor_iPoint;
        nAlreadyVisited = 1;
        StartVisited = 0;

        Area_old = -1;

        while (Area > Area_old) {
          /*
           * - Starting from the closest donor_point, it expands the supermesh by a countour search pattern.
           * - The closest donor element becomes the core, at each iteration a new layer of elements around the core is
           * taken into account
           */

          Area_old = Area;

          ToVisit = nullptr;
          nToVisit = 0;

          for (iNodeVisited = StartVisited; iNodeVisited < nAlreadyVisited; iNodeVisited++) {
            vPoint = alreadyVisitedDonor[iNodeVisited];

            nEdgeVisited = Donor_nLinkedNodes[vPoint];

            for (iEdgeVisited = 0; iEdgeVisited < nEdgeVisited; iEdgeVisited++) {
              donor_iPoint = Donor_LinkedNodes[Donor_StartLinkedNodes[vPoint] + iEdgeVisited];

              /*--- Check if the node to visit is already listed in the data structure to avoid double visits ---*/

              check = false;

              for (jj = 0; jj < nAlreadyVisited; jj++) {
                if (donor_iPoint == alreadyVisitedDonor[jj]) {
                  check = true;
                  break;
                }
              }

              if (check == 0 && ToVisit != nullptr) {
                for (jj = 0; jj < nToVisit; jj++)
                  if (donor_iPoint == ToVisit[jj]) {
                    check = true;
                    break;
                  }
              }

              if (check == 0) {
                /*--- If the node was not already visited, visit it and list it into data structure ---*/

                tmpVect = new unsigned long[nToVisit + 1];

                for (jj = 0; jj < nToVisit; jj++) tmpVect[jj] = ToVisit[jj];
                tmpVect[nToVisit] = donor_iPoint;

                delete[] ToVisit;

                ToVisit = tmpVect;
                tmpVect = nullptr;

                nToVisit++;

                /*--- Find the value of the intersection area between the current donor element and the target element
                 * --- */

                nEdges_donor = Donor_nLinkedNodes[donor_iPoint];

                donor_element = new su2double*[2 * nEdges_donor + 2];
                for (ii = 0; ii < 2 * nEdges_donor + 2; ii++) donor_element[ii] = new su2double[nDim];

                nNode_donor = Build_3D_surface_element(Donor_LinkedNodes, Donor_StartLinkedNodes, Donor_nLinkedNodes,
                                                       DonorPoint_Coord, donor_iPoint, donor_element);

                tmp_Area = 0;
                for (ii = 1; ii < nNode_target - 1; ii++)
                  for (jj = 1; jj < nNode_donor - 1; jj++)
                    tmp_Area += Compute_Triangle_Intersection(target_element[0], target_element[ii],
                                                              target_element[ii + 1], donor_element[0],
                                                              donor_element[jj], donor_element[jj + 1], Normal);

                for (ii = 0; ii < 2 * nEdges_donor + 2; ii++) delete[] donor_element[ii];
                delete[] donor_element;

                /*--- In case the element intersect the target cell update the auxiliary communication data structure
                 * ---*/

                tmp_Coeff_Vect = new su2double[nDonorPoints + 1];
                tmp_Donor_Vect = new unsigned long[nDonorPoints + 1];
                tmp_storeProc = new unsigned long[nDonorPoints + 1];

                for (iDonor = 0; iDonor < nDonorPoints; iDonor++) {
                  tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
                  tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
                  tmp_storeProc[iDonor] = storeProc[iDonor];
                }

                tmp_Coeff_Vect[nDonorPoints] = tmp_Area;
                tmp_Donor_Vect[nDonorPoints] = donor_iPoint;
                tmp_storeProc[nDonorPoints] = Donor_Proc[donor_iPoint];

                delete[] Donor_Vect;
                delete[] Coeff_Vect;
                delete[] storeProc;

                Donor_Vect = tmp_Donor_Vect;
                Coeff_Vect = tmp_Coeff_Vect;
                storeProc = tmp_storeProc;

                tmp_Coeff_Vect = nullptr;
                tmp_Donor_Vect = nullptr;
                tmp_storeProc = nullptr;

                nDonorPoints++;

                Area += tmp_Area;
              }
            }
          }

          /*--- Update auxiliary data structure ---*/

          StartVisited = nAlreadyVisited;

          tmpVect = new unsigned long[nAlreadyVisited + nToVisit];

          for (jj = 0; jj < nAlreadyVisited; jj++) tmpVect[jj] = alreadyVisitedDonor[jj];

          for (jj = 0; jj < nToVisit; jj++) tmpVect[nAlreadyVisited + jj] = ToVisit[jj];

          delete[] alreadyVisitedDonor;

          alreadyVisitedDonor = tmpVect;

          nAlreadyVisited += nToVisit;

          delete[] ToVisit;
        }

        delete[] alreadyVisitedDonor;

        /*--- Set the communication data structure and copy data from the auxiliary vectors ---*/

        targetVertices[markTarget][iVertex].resize(nDonorPoints);

        for (iDonor = 0; iDonor < nDonorPoints; iDonor++) {
          targetVertices[markTarget][iVertex].coefficient[iDonor] = Coeff_Vect[iDonor] / Area;
          targetVertices[markTarget][iVertex].globalPoint[iDonor] = Donor_GlobalPoint[Donor_Vect[iDonor]];
          targetVertices[markTarget][iVertex].processor[iDonor] = storeProc[iDonor];
        }

        for (ii = 0; ii < 2 * nEdges_target + 2; ii++) delete[] target_element[ii];
        delete[] target_element;

        delete[] Donor_Vect;
        Donor_Vect = nullptr;
        delete[] Coeff_Vect;
        Coeff_Vect = nullptr;
        delete[] storeProc;
        storeProc = nullptr;
      }
    }
  }

  delete[] Normal;
  delete[] Direction;

  delete[] Donor_Vect;
  delete[] Coeff_Vect;
  delete[] storeProc;
}

int CSlidingMesh::Build_3D_surface_element(const su2vector<unsigned long>& map,
                                           const su2vector<unsigned long>& startIndex,
                                           const su2vector<unsigned long>& nNeighbor, su2activematrix const& coord,
                                           unsigned long centralNode, su2double** element) {
  /*--- Given a node "centralNode", this routines reconstruct the vertex centered
   *    surface element around the node and store it into "element" ---*/

  constexpr unsigned short nDim = 3;

  const unsigned long* OuterNodes;

  /* --- Store central node as element first point --- */

  for (unsigned short iDim = 0; iDim < nDim; iDim++) element[0][iDim] = coord(centralNode, iDim);

  unsigned long nOuterNodes = nNeighbor[centralNode];

  OuterNodes = &map[startIndex[centralNode]];

  // For each neighbor n of centralNode, store <=2 neighbors of centralNode that are neighbors of n.
  su2matrix<int> OuterNodesNeighbour(nOuterNodes, 2);
  OuterNodesNeighbour = -1;
  // Typically there are exactly 2 such neighbors, and all the neighbors of centralNode can be
  // arranged into a closed chain. However at 1D boundaries of 2D markers, the neighbors might
  // be an open chain, with the two ends having only 1 such neighbor.
  // StartNode is a node where we can start to iterate through the chain. I.e. it's any node in
  // case of a closed chain, or one of the two ends in case of an open chain.
  int StartNode = 0;
  for (unsigned long iNode = 0; iNode < nOuterNodes; iNode++) {
    int count = 0;  // number of neighboring outer nodes already found
    const unsigned long iPoint = OuterNodes[iNode];
    const unsigned long* ptr = &map[startIndex[iPoint]];

    for (unsigned long jNode = 0; jNode < nNeighbor[iPoint]; jNode++) {
      const unsigned long jPoint = ptr[jNode];
      for (unsigned long kNode = 0; kNode < nOuterNodes; kNode++) {
        if (jPoint == OuterNodes[kNode] && jPoint != centralNode) {
          OuterNodesNeighbour(iNode, count) = static_cast<int>(kNode);
          count++;
          break;
        }
      }
    }
    if (count == 1) StartNode = static_cast<int>(iNode);
  }

  /* --- Build element, starts from one outer node and loops along the external edges until the element is reconstructed
   * --- */

  int CurrentNode = StartNode;
  int NextNode = OuterNodesNeighbour(CurrentNode, 0);
  unsigned long iElementNode = 1;

  while (NextNode !=
         -1) {  // We finished iterating through the chain if it is an open chain and we reached the other end.

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      element[iElementNode][iDim] = (element[0][iDim] + coord(OuterNodes[CurrentNode], iDim)) / 2.;

    iElementNode++;

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      element[iElementNode][iDim] =
          (element[0][iDim] + coord[OuterNodes[CurrentNode]][iDim] + coord(OuterNodes[NextNode], iDim)) / 3.;
    iElementNode++;

    // "Place the next domino piece in the correct orientation."
    if (OuterNodesNeighbour(NextNode, 0) == CurrentNode) {
      CurrentNode = NextNode;
      NextNode = OuterNodesNeighbour(NextNode, 1);
    } else {
      CurrentNode = NextNode;
      NextNode = OuterNodesNeighbour(NextNode, 0);
    }

    // We finished iterating through the chain if it is closed and we reached the beginning again.
    if (CurrentNode == StartNode) break;
  }

  if (CurrentNode ==
      StartNode) {  // This is a closed element, so add again element 1 to the end of the structure, useful later
    for (unsigned short iDim = 0; iDim < nDim; iDim++) element[iElementNode][iDim] = element[1][iDim];
    iElementNode++;
  } else {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      element[iElementNode][iDim] = (element[0][iDim] + coord(OuterNodes[CurrentNode], iDim)) / 2.;
    iElementNode++;
  }

  return static_cast<int>(iElementNode);
}

su2double CSlidingMesh::ComputeLineIntersectionLength(unsigned short nDim, const su2double* A1, const su2double* A2,
                                                      const su2double* B1, const su2double* B2,
                                                      const su2double* Direction) {
  /*--- Given 2 segments, each defined by 2 points, it projects them along a given direction
   *    and it computes the length of the segment resulting from their intersection ---*/
  /*--- The algorithm works for both 2D and 3D problems ---*/

  unsigned short iDim;

  su2double dotA2, dotB1, dotB2;

  dotA2 = 0;
  for (iDim = 0; iDim < nDim; iDim++) dotA2 += (A2[iDim] - A1[iDim]) * Direction[iDim];

  if (dotA2 >= 0) {
    dotB1 = 0;
    dotB2 = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      dotB1 += (B1[iDim] - A1[iDim]) * Direction[iDim];
      dotB2 += (B2[iDim] - A1[iDim]) * Direction[iDim];
    }
  } else {
    dotA2 *= -1;

    dotB1 = 0;
    dotB2 = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      dotB1 -= (B1[iDim] - A1[iDim]) * Direction[iDim];
      dotB2 -= (B2[iDim] - A1[iDim]) * Direction[iDim];
    }
  }

  if (dotB1 >= 0 && dotB1 <= dotA2) {
    if (dotB2 < 0) return fabs(dotB1);
    if (dotB2 > dotA2) return fabs(dotA2 - dotB1);

    return fabs(dotB1 - dotB2);
  }

  if (dotB2 >= 0 && dotB2 <= dotA2) {
    if (dotB1 < 0) return fabs(dotB2);
    if (dotB1 > dotA2) return fabs(dotA2 - dotB2);
  }

  if ((dotB1 <= 0 && dotA2 <= dotB2) || (dotB2 <= 0 && dotA2 <= dotB1)) return fabs(dotA2);

  return 0.0;
}

su2double CSlidingMesh::Compute_Triangle_Intersection(const su2double* A1, const su2double* A2, const su2double* A3,
                                                      const su2double* B1, const su2double* B2, const su2double* B3,
                                                      const su2double* Direction) {
  /* --- This routine is ONLY for 3D grids --- */
  /* --- Projects triangle points onto a plane, specified by its normal "Direction", and calls the
   * ComputeIntersectionArea routine --- */

  unsigned short iDim;
  constexpr unsigned short nDim = 3;

  su2double I[3], J[3], K[3];
  su2double a1[3], a2[3], a3[3];
  su2double b1[3], b2[3], b3[3];
  su2double m1, m2;

  /* --- Reference frame is determined by: x = A1A2 y = x ^ ( -Direction ) --- */

  for (iDim = 0; iDim < 3; iDim++) {
    a1[iDim] = 0;
    a2[iDim] = 0;
    a3[iDim] = 0;

    b1[iDim] = 0;
    b2[iDim] = 0;
    b3[iDim] = 0;
  }

  m1 = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    K[iDim] = Direction[iDim];

    m1 += K[iDim] * K[iDim];
  }

  for (iDim = 0; iDim < nDim; iDim++) K[iDim] /= sqrt(m1);

  m2 = 0;
  for (iDim = 0; iDim < nDim; iDim++) m2 += (A2[iDim] - A1[iDim]) * K[iDim];

  m1 = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    I[iDim] = (A2[iDim] - A1[iDim]) - m2 * K[iDim];
    m1 += I[iDim] * I[iDim];
  }

  for (iDim = 0; iDim < nDim; iDim++) I[iDim] /= sqrt(m1);

  // Cross product to find Y
  J[0] = K[1] * I[2] - K[2] * I[1];
  J[1] = -(K[0] * I[2] - K[2] * I[0]);
  J[2] = K[0] * I[1] - K[1] * I[0];

  /* --- Project all points on the plane specified by Direction and change their reference frame taking A1 as origin ---
   */

  for (iDim = 0; iDim < nDim; iDim++) {
    a2[0] += (A2[iDim] - A1[iDim]) * I[iDim];
    a2[1] += (A2[iDim] - A1[iDim]) * J[iDim];
    a2[2] += (A2[iDim] - A1[iDim]) * K[iDim];

    a3[0] += (A3[iDim] - A1[iDim]) * I[iDim];
    a3[1] += (A3[iDim] - A1[iDim]) * J[iDim];
    a3[2] += (A3[iDim] - A1[iDim]) * K[iDim];

    b1[0] += (B1[iDim] - A1[iDim]) * I[iDim];
    b1[1] += (B1[iDim] - A1[iDim]) * J[iDim];
    b1[2] += (B1[iDim] - A1[iDim]) * K[iDim];

    b2[0] += (B2[iDim] - A1[iDim]) * I[iDim];
    b2[1] += (B2[iDim] - A1[iDim]) * J[iDim];
    b2[2] += (B2[iDim] - A1[iDim]) * K[iDim];

    b3[0] += (B3[iDim] - A1[iDim]) * I[iDim];
    b3[1] += (B3[iDim] - A1[iDim]) * J[iDim];
    b3[2] += (B3[iDim] - A1[iDim]) * K[iDim];
  }

  /*--- Compute intersection area ---*/

  return ComputeIntersectionArea(a1, a2, a3, b1, b2, b3);
}

su2double CSlidingMesh::ComputeIntersectionArea(const su2double* P1, const su2double* P2, const su2double* P3,
                                                const su2double* Q1, const su2double* Q2, const su2double* Q3) {
  /* --- This routines computes the area of the polygonal element generated by the superimposition of 2 planar triangle
   * --- */
  /* --- The 2 triangle must lie on the same plane --- */

  unsigned short iDim, nPoints = 0, i, j, k;
  unsigned short min_theta_index;

  su2double points[16][2], IntersectionPoint[2], theta[6];
  su2double TriangleP[4][2], TriangleQ[4][2];
  su2double Area, det, dot1, dot2, dtmp, min_theta;

  constexpr unsigned short nDim = 2;

  for (iDim = 0; iDim < nDim; iDim++) {
    TriangleP[0][iDim] = 0;
    TriangleP[1][iDim] = P2[iDim] - P1[iDim];
    TriangleP[2][iDim] = P3[iDim] - P1[iDim];
    TriangleP[3][iDim] = 0;

    TriangleQ[0][iDim] = Q1[iDim] - P1[iDim];
    TriangleQ[1][iDim] = Q2[iDim] - P1[iDim];
    TriangleQ[2][iDim] = Q3[iDim] - P1[iDim];
    TriangleQ[3][iDim] = Q1[iDim] - P1[iDim];
  }

  for (j = 0; j < 3; j++) {
    if (CheckPointInsideTriangle(TriangleP[j], TriangleQ[0], TriangleQ[1], TriangleQ[2])) {
      // Then P1 is also inside triangle Q, so store it
      for (iDim = 0; iDim < nDim; iDim++) points[nPoints][iDim] = TriangleP[j][iDim];

      nPoints++;
    }
  }

  for (j = 0; j < 3; j++) {
    if (CheckPointInsideTriangle(TriangleQ[j], TriangleP[0], TriangleP[1], TriangleP[2])) {
      // Then Q1 is also inside triangle P, so store it
      for (iDim = 0; iDim < nDim; iDim++) points[nPoints][iDim] = TriangleQ[j][iDim];

      nPoints++;
    }
  }

  // Compute all edge intersections

  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      det = (TriangleP[j][0] - TriangleP[j + 1][0]) * (TriangleQ[i][1] - TriangleQ[i + 1][1]) -
            (TriangleP[j][1] - TriangleP[j + 1][1]) * (TriangleQ[i][0] - TriangleQ[i + 1][0]);

      if (det != 0.0) {
        ComputeLineIntersectionPoint(TriangleP[j], TriangleP[j + 1], TriangleQ[i], TriangleQ[i + 1], IntersectionPoint);

        dot1 = 0;
        dot2 = 0;
        for (iDim = 0; iDim < nDim; iDim++) {
          dot1 += (TriangleP[j][iDim] - IntersectionPoint[iDim]) * (TriangleP[j + 1][iDim] - IntersectionPoint[iDim]);
          dot2 += (TriangleQ[i][iDim] - IntersectionPoint[iDim]) * (TriangleQ[i + 1][iDim] - IntersectionPoint[iDim]);
        }

        if (dot1 <= 0 && dot2 <= 0) {  // It found one intersection

          // Store temporarily the intersection point

          for (iDim = 0; iDim < nDim; iDim++) points[nPoints][iDim] = IntersectionPoint[iDim];

          nPoints++;
        }
      }
    }
  }

  // Remove double points, if any

  for (i = 0; i < nPoints; i++) {
    for (j = i + 1; j < nPoints; j++) {
      if (points[j][0] == points[i][0] && points[j][1] == points[i][1]) {
        for (k = j; k < nPoints - 1; k++) {
          points[k][0] = points[k + 1][0];
          points[k][1] = points[k + 1][1];
        }
        nPoints--;
        j--;
      }
    }
  }

  // Re-order nodes

  for (i = 1; i < nPoints; i++) {  // Change again reference frame
    for (iDim = 0; iDim < nDim; iDim++) points[i][iDim] -= points[0][iDim];

    // Compute polar azimuth for each node but the first
    theta[i] = atan2(points[i][1], points[i][0]);
  }

  for (iDim = 0; iDim < nDim; iDim++) points[0][iDim] = 0;

  for (i = 1; i < nPoints; i++) {
    min_theta = theta[i];
    min_theta_index = 0;

    for (j = i + 1; j < nPoints; j++) {
      if (theta[j] < min_theta) {
        min_theta = theta[j];
        min_theta_index = j;
      }
    }

    if (min_theta_index != 0) {
      dtmp = theta[i];
      theta[i] = theta[min_theta_index];
      theta[min_theta_index] = dtmp;

      dtmp = points[i][0];
      points[i][0] = points[min_theta_index][0];
      points[min_theta_index][0] = dtmp;

      dtmp = points[i][1];
      points[i][1] = points[min_theta_index][1];
      points[min_theta_index][1] = dtmp;
    }
  }

  // compute area using cross product rule, points position are referred to the
  // 2-dimensional, local, reference frame centered in points[0]

  Area = 0;

  if (nPoints > 2) {
    for (i = 1; i < nPoints - 1; i++) {
      // Ax*By
      Area += (points[i][0] - points[0][0]) * (points[i + 1][1] - points[0][1]);

      // Ay*Bx
      Area -= (points[i][1] - points[0][1]) * (points[i + 1][0] - points[0][0]);
    }
  }

  return fabs(Area) / 2;
}

void CSlidingMesh::ComputeLineIntersectionPoint(const su2double* A1, const su2double* A2, const su2double* B1,
                                                const su2double* B2, su2double* IntersectionPoint) {
  /* --- Uses determinant rule to compute the intersection point between 2 straight segments --- */
  /* This works only for lines on a 2D plane, A1, A2 and B1, B2 are respectively the head and the tail points of each
   * segment, since they're on a 2D plane they are defined by a 2-elements array containing their coordinates */

  su2double det;

  det = (A1[0] - A2[0]) * (B1[1] - B2[1]) - (A1[1] - A2[1]) * (B1[0] - B2[0]);

  if (det != 0.0) {  // else there is no intersection point
    IntersectionPoint[0] =
        ((A1[0] * A2[1] - A1[1] * A2[0]) * (B1[0] - B2[0]) - (B1[0] * B2[1] - B1[1] * B2[0]) * (A1[0] - A2[0])) / det;
    IntersectionPoint[1] =
        ((A1[0] * A2[1] - A1[1] * A2[0]) * (B1[1] - B2[1]) - (B1[0] * B2[1] - B1[1] * B2[0]) * (A1[1] - A2[1])) / det;
  }
}

bool CSlidingMesh::CheckPointInsideTriangle(const su2double* Point, const su2double* T1, const su2double* T2,
                                            const su2double* T3) {
  /* --- Check whether a point "Point" lies inside or outside a triangle defined by 3 points "T1", "T2", "T3" --- */
  /* For each edge it checks on which side the point lies:
   * - Computes the unit vector pointing at the internal side of the edge
   * - Comutes the vector that connects the point to a point along the edge
   * - If the dot product is positive it means that the point is on the internal side of the edge
   * - If the check is positive for all the 3 edges, then the point lies within the triangle
   */

  unsigned short iDim, check = 0;

  su2double vect1[2], vect2[2], r[2];
  su2double dot;

  constexpr unsigned short nDim = 2;

  /* --- Check first edge --- */

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vect1[iDim] = T3[iDim] - T1[iDim];  // vec 1 is aligned to the edge
    vect2[iDim] = T2[iDim] - T1[iDim];  // vect 2 is the vector connecting one edge point to the third triangle vertex

    r[iDim] = Point[iDim] - T1[iDim];  // Connects point to vertex T1

    dot += vect2[iDim] * vect2[iDim];
  }
  dot = sqrt(dot);

  for (iDim = 0; iDim < nDim; iDim++) vect2[iDim] /= dot;

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) dot += vect1[iDim] * vect2[iDim];

  for (iDim = 0; iDim < nDim; iDim++)
    vect1[iDim] = T3[iDim] - (T1[iDim] + dot * vect2[iDim]);  // Computes the inward unit vector

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++)  // Checs that the point lies on the internal plane
    dot += vect1[iDim] * r[iDim];

  if (dot >= 0) check++;

  /* --- Check second edge --- */

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vect1[iDim] = T1[iDim] - T2[iDim];
    vect2[iDim] = T3[iDim] - T2[iDim];

    r[iDim] = Point[iDim] - T2[iDim];

    dot += vect2[iDim] * vect2[iDim];
  }
  dot = sqrt(dot);

  for (iDim = 0; iDim < nDim; iDim++) vect2[iDim] /= dot;

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) dot += vect1[iDim] * vect2[iDim];

  for (iDim = 0; iDim < nDim; iDim++) vect1[iDim] = T1[iDim] - (T2[iDim] + dot * vect2[iDim]);

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) dot += vect1[iDim] * r[iDim];

  if (dot >= 0) check++;

  /* --- Check third edge --- */

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vect1[iDim] = T2[iDim] - T3[iDim];
    vect2[iDim] = T1[iDim] - T3[iDim];

    r[iDim] = Point[iDim] - T3[iDim];

    dot += vect2[iDim] * vect2[iDim];
  }
  dot = sqrt(dot);

  for (iDim = 0; iDim < nDim; iDim++) vect2[iDim] /= dot;

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) dot += vect1[iDim] * vect2[iDim];

  for (iDim = 0; iDim < nDim; iDim++) vect1[iDim] = T2[iDim] - (T3[iDim] + dot * vect2[iDim]);

  dot = 0;
  for (iDim = 0; iDim < nDim; iDim++) dot += vect1[iDim] * r[iDim];

  if (dot >= 0) check++;

  return (check == 3);
}
