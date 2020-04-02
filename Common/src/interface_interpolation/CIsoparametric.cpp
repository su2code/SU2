/*!
 * \file CIsoparametric.cpp
 * \brief Implementation isoparametric interpolation (using FE shape functions).
 * \author H. Kline, P. Gomes
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/interface_interpolation/CIsoparametric.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/geometry/elements/CElement.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"
#include <unordered_map>

using namespace GeometryToolbox;


CIsoparametric::CIsoparametric(CGeometry ****geometry_container, const CConfig* const* config, unsigned int iZone,
                               unsigned int jZone) : CInterpolator(geometry_container, config, iZone, jZone) {
  Set_TransferCoeff(config);
}

void CIsoparametric::PrintStatistics(void) const {
  if (rank != MASTER_NODE) return;
  cout << " Maximum distance to closest donor element: " << MaxDistance << ".\n"
       << " Interpolation mitigated on " << ErrorCounter << " (" << ErrorRate << "%) target vertices." << endl;
}

void CIsoparametric::Set_TransferCoeff(const CConfig* const* config) {

  /*--- Angle between target and donor below which we trigger fallback measures. ---*/
  const su2double thetaMin = 80.0/180.0*PI_NUMBER;

  const int nProcessor = size;
  const auto nMarkerInt = config[donorZone]->GetMarker_n_ZoneInterface()/2;
  const auto nDim = donor_geometry->GetnDim();

  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];

  /*--- Init stats. ---*/
  MaxDistance = 0.0; ErrorCounter = 0;
  unsigned long nGlobalVertexTarget = 0;

  /*--- Cycle over nMarkersInt interface to determine communication pattern. ---*/

  for (unsigned short iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

    /* High level procedure:
     * - Loop through vertices of the target grid;
     * - Find nearest element;
     * - Compute and set the transfer coefficients.
     */

    /*--- On the donor side: find the tag of the boundary sharing the interface. ---*/
    const auto markDonor = Find_InterfaceMarker(config[donorZone], iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface. ---*/
    const auto markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step. ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0, nVertexTarget = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex(markDonor);
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor. ---*/
    Determine_ArraySize(markDonor, markTarget, nVertexDonor, nDim);

    const auto nGlobalVertexDonor = accumulate(Buffer_Receive_nVertex_Donor,
                                    Buffer_Receive_nVertex_Donor+nProcessor, 0ul);

    Buffer_Send_Coord = new su2double [ MaxLocalVertex_Donor * nDim ];
    Buffer_Send_GlobalPoint = new long [ MaxLocalVertex_Donor ];
    Buffer_Receive_Coord = new su2double [ nProcessor * MaxLocalVertex_Donor * nDim ];
    Buffer_Receive_GlobalPoint = new long [ nProcessor * MaxLocalVertex_Donor ];

    /*--- Collect coordinates and global point indices. ---*/
    Collect_VertexInfo(markDonor, markTarget, nVertexDonor, nDim);

    /*--- Compress the vertex information, and build a map of global point to "compressed
     *    index" to then reconstruct the donor elements in local index space. ---*/

    su2activematrix donorCoord(nGlobalVertexDonor, nDim);
    vector<long> donorPoint(nGlobalVertexDonor);
    vector<int> donorProc(nGlobalVertexDonor);
    unordered_map<long, unsigned long> globalToLocalMap;

    auto iCount = 0ul;
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor) {
      auto offset = iProcessor * MaxLocalVertex_Donor;
      for (auto iVertex = 0ul; iVertex < Buffer_Receive_nVertex_Donor[iProcessor]; ++iVertex) {
        for (int iDim = 0; iDim < nDim; ++iDim)
          donorCoord(iCount,iDim) = Buffer_Receive_Coord[(offset+iVertex)*nDim + iDim];
        donorPoint[iCount] = Buffer_Receive_GlobalPoint[offset+iVertex];
        donorProc[iCount] = iProcessor;
        assert((globalToLocalMap.count(donorPoint[iCount]) == 0) && "Duplicate donor point found.");
        globalToLocalMap[donorPoint[iCount]] = iCount;
        ++iCount;
      }
    }
    assert((iCount == nGlobalVertexDonor) && "Global donor point count mismatch.");

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_GlobalPoint;
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_GlobalPoint;

    /*--- Collect donor element (face) information. ---*/

    vector<unsigned long> allNumElem;
    vector<unsigned short> elemNumNodes;
    su2matrix<long> elemIdxNodes;

    auto nGlobalElemDonor = Collect_ElementInfo(markDonor, nDim, true,
                              allNumElem, elemNumNodes, elemIdxNodes);

    su2activematrix elemCentroid(nGlobalElemDonor, nDim);

    SU2_OMP_PARALLEL
    {
    /*--- Compute element centroids to then find the closest one to a given target point. ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(nGlobalElemDonor,omp_get_max_threads()))
    for (auto iElem = 0u; iElem < nGlobalElemDonor; ++iElem) {

      const auto nNode = elemNumNodes[iElem];

      for (auto iDim = 0u; iDim < nDim; ++iDim)
        elemCentroid(iElem, iDim) = 0.0;

      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        /*--- Map the node to "local" coordinates. ---*/
        assert(globalToLocalMap.count(elemIdxNodes(iElem,iNode)) &&
               "Unknown donor point referenced by donor element.");
        const auto iVertex = globalToLocalMap.at(elemIdxNodes(iElem,iNode));
        elemIdxNodes(iElem,iNode) = iVertex;

        for (auto iDim = 0u; iDim < nDim; ++iDim)
          elemCentroid(iElem, iDim) += donorCoord(iVertex, iDim) / nNode;
      }
    }

    /*--- Compute transfer coefficients for each target point. ---*/
    su2double maxDist = 0.0;
    unsigned long errorCount = 0, totalCount = 0;

    SU2_OMP_FOR_DYN(roundUpDiv(nVertexTarget,2*omp_get_max_threads()))
    for (auto iVertex = 0u; iVertex < nVertexTarget; ++iVertex) {

      auto target_vertex = target_geometry->vertex[markTarget][iVertex];
      const auto iPoint = target_vertex->GetNode();

      if (!target_geometry->node[iPoint]->GetDomain()) continue;
      totalCount += 1;

      /*--- Coordinates and normal of the target point. ---*/
      const su2double* coord_i = target_geometry->node[iPoint]->GetCoord();
      const su2double* normal_i = target_vertex->GetNormal();

      /*--- Find closest element (the naive way). ---*/
      su2double minDist = 1e20;
      auto iElemDonor = 0u;
      for (auto iElem = 0u; iElem < nGlobalElemDonor; ++iElem) {
        su2double dist = SquaredDistance(nDim, coord_i, elemCentroid[iElem]);
        if (dist < minDist) {
          minDist = dist;
          iElemDonor = iElem;
        }
      }

      /*--- Fetch donor element info. ---*/
      int procList[4] = {0};
      long nodeList[4] = {0};
      su2double coords[4][3] = {{0.0}};

      const auto nNode = elemNumNodes[iElemDonor];

      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        const auto iVertex = elemIdxNodes(iElemDonor, iNode);
        procList[iNode] = donorProc[iVertex];
        nodeList[iNode] = donorPoint[iVertex];
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          coords[iNode][iDim] = donorCoord(iVertex,iDim);
      }

      const su2double* coord_j = elemCentroid[iElemDonor];
      su2double normal_j[3] = {0.0};

      switch (nNode) {
        case 2: LineNormal(coords, normal_j); break;
        case 3: TriangleNormal(coords, normal_j); break;
        case 4: QuadrilateralNormal(coords, normal_j); break;
      }

      /*--- Project the target onto the donor plane, as grids may not be exactly conformal.
       *    For quadrilaterals this is approximate as the element may be warped. ---*/
      su2double proj = DotProduct(nDim, normal_i, normal_j);
      su2double theta = acos(fabs(proj) / (Norm(nDim, normal_i)*Norm(nDim, normal_j)));

      su2double projCoord[3] = {0.0};

      if (theta >= thetaMin) {
        su2double dist_ij[3] = {0.0};
        Distance(nDim, coord_j, coord_i, dist_ij);
        proj = DotProduct(nDim, dist_ij, normal_j) / proj;
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          projCoord[iDim] = coord_i[iDim] + proj*dist_ij[iDim];

        maxDist = max(maxDist, proj*Norm(nDim, dist_ij));
      }
      else {
        /*--- Target and donor are too out of alignement, as fallback
         *    use the element centroid as the projected coordinate. ---*/
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          projCoord[iDim] = coord_j[iDim];

        errorCount += 1;
        maxDist = max(maxDist, sqrt(minDist));
      }

      /*--- Compute and set interpolation coefficients. ---*/

      su2double isoparams[4] = {0.0};
      switch (nNode) {
        case 2: errorCount += LineIsoparameters(coords, projCoord, isoparams); break;
        case 3: errorCount += TriangleIsoparameters(coords, projCoord, isoparams); break;
        case 4: errorCount += QuadrilateralIsoparameters(coords, projCoord, isoparams); break;
      }

      target_vertex->Allocate_DonorInfo(nNode);

      for (auto iDonor = 0u; iDonor < nNode; ++iDonor) {
        target_vertex->SetDonorCoeff(iDonor, isoparams[iDonor]);
        target_vertex->SetInterpDonorPoint(iDonor, nodeList[iDonor]);
        target_vertex->SetInterpDonorProcessor(iDonor, procList[iDonor]);
      }

    }
    SU2_OMP_CRITICAL
    {
      MaxDistance = max(MaxDistance,maxDist);
      ErrorCounter += errorCount;
      nGlobalVertexTarget += totalCount;
    }
    } // end SU2_OMP_PARALLEL

  } // end nMarkerInt loop

  /*--- Final reduction of statistics. ---*/
  su2double tmp = MaxDistance;
  unsigned long tmp1 = ErrorCounter, tmp2 = nGlobalVertexTarget;
  SU2_MPI::Allreduce(&tmp, &MaxDistance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&tmp1, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&tmp2, &nGlobalVertexTarget, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  ErrorRate = 100*su2double(ErrorCounter) / nGlobalVertexTarget;

}

int CIsoparametric::LineIsoparameters(const su2double X[][3], const su2double *xj, su2double *isoparams) {

  su2double l01 = Distance(2, X[0], X[1]);
  su2double l0j = Distance(2, X[0], xj);
  su2double lj1 = Distance(2, xj, X[1]);

  /*--- Detect out of bounds point. ---*/

  const int outOfBounds = (l0j+lj1) > l01;

  if (outOfBounds) {
    l0j = (l0j > lj1)? l01 : 0.0; // which ever is closest becomes the donor
    lj1 = l01 - l0j;
  }

  isoparams[0] = lj1 / l01;
  isoparams[1] = 1.0 - isoparams[0];

  return outOfBounds;
}

int CIsoparametric::TriangleIsoparameters(const su2double X[][3], const su2double *xj, su2double *isoparams) {

  /*--- The isoparameters are the solution to the determined system X^T * isoparams = xj.
   *    This is consistent with the shape functions of the linear triangular element. ---*/

  su2double A[3][3] = {{0.0}}; // = X^T

  for (int i = 0; i < 3; ++i) {
    isoparams[i] = xj[i]; // use isoparams as rhs
    for (int j = 0; j < 3; ++j)
      A[i][j] = X[j][i];
  }

  /*--- Solve normal system by in-place Gaussian elimination without pivoting. ---*/

  /*--- Transform system in Upper Matrix. ---*/
  for (int i = 1; i < 3; ++i) {
    for (int j = 0; j < i; ++j) {
      su2double w = A[i][j] / A[j][j];
      for (int k = j; k < 3; ++k)
        A[i][k] -= w * A[j][k];
      isoparams[i] -= w * isoparams[j];
    }
  }

  /*--- Backwards substitution. ---*/
  for (int i = 2; i >= 0; --i) {
    for (int j = i+1; j < 3; ++j)
      isoparams[i] -= A[i][j] * isoparams[j];
    isoparams[i] /= A[i][i];
  }

  /*--- Detect out of bounds point. ---*/
  const int outOfBounds = (isoparams[0] < 0.0) || (isoparams[1] < 0.0) || (isoparams[2] < 0.0);

  /*--- Simple mitigation. ---*/
  if (outOfBounds)
    isoparams[0] = isoparams[1] = isoparams[2] = 1.0/3.0;

  return outOfBounds;
}

int CIsoparametric::QuadrilateralIsoparameters(const su2double X[][3], const su2double *xj, su2double *isoparams) {

  /*--- The isoparameters are the shape functions (Ni) evaluated at xj, for that we need
   *    the corresponding Xi and Eta, which are obtained by solving the overdetermined
   *    nonlinear system xj - X^T * Ni(Xi,Eta) = 0 via the Gauss-Newton method. ---*/

  constexpr int NITER = 10;
  const su2double tol = 1e-12;

  su2double Xi = 0.0, Eta = 0.0, eps;

  /*--- Finding Xi and Eta is a "third order" effect that we do not
   *    differentiate (also because we need to iterate). ---*/
  const bool tapeActive = AD::TapeActive();
  AD::StopRecording();

  for (int iter = 0; iter < NITER; ++iter) {

    /*--- Evaluate the residual. ---*/
    su2double r[3] = {xj[0], xj[1], xj[2]};
    su2double Ni[4] = {0.0};
    CQUAD4::ShapeFunctions(Xi, Eta, Ni);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 4; ++j)
        r[i] -= X[j][i] * Ni[j];

    /*--- Evaluate the residual Jacobian. ---*/
    su2double dNi[4][2] = {{0.0}};
    CQUAD4::ShapeFunctionJacobian(Xi, Eta, dNi);

    su2double jac[3][2];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 4; ++k)
          jac[i][j] -= X[k][i] * dNi[k][j];

    /*--- Compute the correction (normal equations and Cramer's rule). ---*/
    su2double A[2][2] = {{0.0}}, b[2] = {0.0};
    for (int i = 0; i < 2; ++i) {
      for (int j = i; j < 2; ++j)
        for (int k = 0; k < 3; ++k)
          A[i][j] += jac[k][i] * jac[k][j];
      A[1][0] = A[0][1];

      for (int k = 0; k < 3; ++k)
        b[i] += jac[k][i] * r[k];
    }

    su2double detA = 1.0 / (A[0][0]*A[1][1] - A[0][1]*A[1][0]);
    su2double dXi = (b[0]*A[1][1] - b[1]*A[0][1]) * detA;
    su2double dEta = (A[0][0]*b[1] - A[1][0]*b[0]) * detA;
    Xi -= dXi;
    Eta -= dEta;

    eps = fabs(dXi)+fabs(dEta);
    if (eps < tol) break;
  }

  if (tapeActive) AD::StartRecording();

  int outOfBounds = 0;

  if (eps > 0.01) {
    /*--- Iteration did not converge. ---*/
    Xi = Eta = 0.0;
    outOfBounds = 1;
  }
  else {
    /*--- Check bounds. ---*/
    outOfBounds = (fabs(Xi) > 1.0) || (fabs(Eta) > 1.0);

    /*--- Mitigate by clamping coordinates. ---*/
    if (outOfBounds) {
      Xi = max(-1.0, min(Xi, 1.0));
      Eta = max(-1.0, min(Eta, 1.0));
    }
  }

  /*--- Evaluate isoparameters. ---*/
  CQUAD4::ShapeFunctions(Xi, Eta, isoparams);

  return outOfBounds;
}
