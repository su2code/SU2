/*!
 * \file CIsoparametric.cpp
 * \brief Implementation isoparametric interpolation (using FE shape functions).
 * \author P. Gomes
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

#include "../../include/interface_interpolation/CIsoparametric.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/geometry/elements/CElement.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"
#include <unordered_map>

using namespace GeometryToolbox;

CIsoparametric::CIsoparametric(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                               unsigned int jZone)
    : CInterpolator(geometry_container, config, iZone, jZone) {
  SetTransferCoeff(config);
}

void CIsoparametric::PrintStatistics() const {
  if (rank != MASTER_NODE) return;
  cout << "  Maximum distance to closest donor element: " << MaxDistance << ".\n"
       << "  Interpolation clipped for " << ErrorCounter << " (" << ErrorRate << "%) target vertices." << endl;
}

void CIsoparametric::SetTransferCoeff(const CConfig* const* config) {
  const su2double matchingVertexTol = 1e-12;  // 1um^2

  const int nProcessor = size;
  const auto nMarkerInt = config[donorZone]->GetMarker_n_ZoneInterface() / 2;
  const auto nDim = donor_geometry->GetnDim();

  Buffer_Receive_nVertex_Donor = new unsigned long[nProcessor];

  /*--- Make space for donor info. ---*/

  targetVertices.resize(config[targetZone]->GetnMarker_All());

  /*--- Init stats. ---*/
  MaxDistance = 0.0;
  ErrorCounter = 0;
  unsigned long nGlobalVertexTarget = 0;

  /*--- Cycle over nMarkersInt interface to determine communication pattern. ---*/

  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; iMarkerInt++) {
    /* High level procedure:
     * - Loop through vertices of the target grid;
     * - Find nearest element;
     * - Compute and set the transfer coefficients.
     */

    /*--- On the donor side: find the tag of the boundary sharing the interface. ---*/
    const auto markDonor = config[donorZone]->FindInterfaceMarker(iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface. ---*/
    const auto markTarget = config[targetZone]->FindInterfaceMarker(iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step. ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0, nVertexTarget = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex(markDonor);
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor. ---*/
    Determine_ArraySize(markDonor, markTarget, nVertexDonor, nDim);

    const auto nGlobalVertexDonor =
        accumulate(Buffer_Receive_nVertex_Donor, Buffer_Receive_nVertex_Donor + nProcessor, 0ul);

    Buffer_Send_Coord.resize(MaxLocalVertex_Donor, nDim);
    Buffer_Send_GlobalPoint.resize(MaxLocalVertex_Donor);
    Buffer_Receive_Coord.resize(nProcessor * MaxLocalVertex_Donor, nDim);
    Buffer_Receive_GlobalPoint.resize(nProcessor * MaxLocalVertex_Donor);

    /*--- Collect coordinates and global point indices. ---*/
    Collect_VertexInfo(markDonor, markTarget, nVertexDonor, nDim);
    if (nVertexTarget) targetVertices[markTarget].resize(nVertexTarget);

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
        for (int iDim = 0; iDim < nDim; ++iDim) donorCoord(iCount, iDim) = Buffer_Receive_Coord(offset + iVertex, iDim);
        donorPoint[iCount] = Buffer_Receive_GlobalPoint[offset + iVertex];
        donorProc[iCount] = iProcessor;
        assert((globalToLocalMap.count(donorPoint[iCount]) == 0) && "Duplicate donor point found.");
        globalToLocalMap[donorPoint[iCount]] = iCount;
        ++iCount;
      }
    }
    assert((iCount == nGlobalVertexDonor) && "Global donor point count mismatch.");

    /*--- Collect donor element (face) information. ---*/

    vector<unsigned long> allNumElem;
    vector<unsigned short> elemNumNodes;
    su2matrix<long> elemIdxNodes;

    const auto nGlobalElemDonor = Collect_ElementInfo(markDonor, nDim, true, allNumElem, elemNumNodes, elemIdxNodes);

    /*--- Map the node to "local" indices and create a list of connected elements for each vertex. ---*/

    vector<vector<unsigned> > vertexElements(nGlobalVertexDonor);

    for (auto iElem = 0u; iElem < nGlobalElemDonor; ++iElem) {
      const auto nNode = elemNumNodes[iElem];

      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        assert(globalToLocalMap.count(elemIdxNodes(iElem, iNode)) &&
               "Unknown donor point referenced by donor element.");
        const auto iVertex = globalToLocalMap.at(elemIdxNodes(iElem, iNode));
        elemIdxNodes(iElem, iNode) = iVertex;

        vertexElements[iVertex].push_back(iElem);
      }
    }

    /*--- Compute transfer coefficients for each target point. ---*/
    SU2_OMP_PARALLEL {
      su2double maxDist = 0.0;
      unsigned long errorCount = 0, totalCount = 0;

      SU2_OMP_FOR_DYN(roundUpDiv(nVertexTarget, 2 * omp_get_max_threads()))
      for (auto iVertexTarget = 0u; iVertexTarget < nVertexTarget; ++iVertexTarget) {
        auto& target_vertex = targetVertices[markTarget][iVertexTarget];
        const auto iPoint = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

        if (!target_geometry->nodes->GetDomain(iPoint)) continue;
        totalCount += 1;

        /*--- Coordinates of the target point. ---*/
        const su2double* coord_i = target_geometry->nodes->GetCoord(iPoint);

        /*--- Find the closest donor vertex. ---*/
        su2double minDist = 1e9;
        unsigned iClosestVertex = 0;
        for (auto iVertexDonor = 0u; iVertexDonor < nGlobalVertexDonor; ++iVertexDonor) {
          su2double d = SquaredDistance(nDim, coord_i, donorCoord[iVertexDonor]);
          if (d < minDist) {
            minDist = d;
            iClosestVertex = iVertexDonor;
          }
        }

        if (minDist < matchingVertexTol) {
          /*--- Perfect match. ---*/
          target_vertex.resize(1);
          target_vertex.coefficient[0] = 1.0;
          target_vertex.globalPoint[0] = donorPoint[iClosestVertex];
          target_vertex.processor[0] = donorProc[iClosestVertex];
          continue;
        }

        /*--- Evaluate interpolation for the elements connected to the closest vertex. ---*/
        DonorInfo donor;
        donor.error = 2;
        donor.distance = 1e9;
        for (auto iElem : vertexElements[iClosestVertex]) {
          /*--- Fetch element info. ---*/
          DonorInfo candidate;
          candidate.iElem = iElem;
          const auto nNode = elemNumNodes[iElem];
          su2double coords[4][3] = {{0.0}};

          for (auto iNode = 0u; iNode < nNode; ++iNode) {
            const auto iVertex = elemIdxNodes(iElem, iNode);
            for (auto iDim = 0u; iDim < nDim; ++iDim) coords[iNode][iDim] = donorCoord(iVertex, iDim);
          }

          /*--- Compute the interpolation coefficients. ---*/
          switch (nNode) {
            case 2:
              candidate.error = LineIsoparameters(coords, coord_i, candidate.isoparams);
              break;
            case 3:
              candidate.error = TriangleIsoparameters(coords, coord_i, candidate.isoparams);
              break;
            case 4:
              candidate.error = QuadrilateralIsoparameters(coords, coord_i, candidate.isoparams);
              break;
          }

          /*--- Evaluate distance from target to final mapped point. ---*/
          su2double finalCoord[3] = {0.0};
          for (auto iDim = 0u; iDim < nDim; ++iDim)
            for (auto iNode = 0u; iNode < nNode; ++iNode)
              finalCoord[iDim] += coords[iNode][iDim] * candidate.isoparams[iNode];

          candidate.distance = Distance(nDim, coord_i, finalCoord);

          /*--- Detect a very bad candidate (NaN). ---*/
          if (candidate.distance != candidate.distance) continue;

          /*--- Check if the candidate is an improvement, update donor if so. ---*/
          if (candidate < donor) donor = candidate;
        }

        if (donor.error > 1) SU2_MPI::Error("Isoparametric interpolation failed, NaN detected.", CURRENT_FUNCTION);

        errorCount += donor.error;
        maxDist = max(maxDist, donor.distance);

        const auto nNode = elemNumNodes[donor.iElem];

        target_vertex.resize(nNode);

        for (auto iNode = 0u; iNode < nNode; ++iNode) {
          const auto iVertex = elemIdxNodes(donor.iElem, iNode);
          target_vertex.coefficient[iNode] = donor.isoparams[iNode];
          target_vertex.globalPoint[iNode] = donorPoint[iVertex];
          target_vertex.processor[iNode] = donorProc[iVertex];
        }
      }
      END_SU2_OMP_FOR
      SU2_OMP_CRITICAL {
        MaxDistance = max(MaxDistance, maxDist);
        ErrorCounter += errorCount;
        nGlobalVertexTarget += totalCount;
      }
      END_SU2_OMP_CRITICAL
    }
    END_SU2_OMP_PARALLEL

  }  // end nMarkerInt loop

  /*--- Final reduction of statistics. ---*/
  su2double tmp = MaxDistance;
  unsigned long tmp1 = ErrorCounter, tmp2 = nGlobalVertexTarget;
  SU2_MPI::Allreduce(&tmp, &MaxDistance, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&tmp1, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&tmp2, &nGlobalVertexTarget, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  ErrorRate = 100 * su2double(ErrorCounter) / nGlobalVertexTarget;
}

int CIsoparametric::LineIsoparameters(const su2double X[][3], const su2double* xj, su2double* isoparams) {
  /*--- Project the target point onto the line. ---*/

  su2double normal[2] = {0.0};
  LineNormal(X, normal);
  su2double xprj[2] = {0.0};
  PointPlaneProjection<su2double, 2>(xj, X[0], normal, xprj);

  su2double l01 = Distance(2, X[0], X[1]);
  su2double l0j = Distance(2, X[0], xprj);
  su2double lj1 = Distance(2, xprj, X[1]);

  /*--- Detect out of bounds point. ---*/

  const int outOfBounds = (l0j + lj1) > (2 * l01);

  isoparams[0] = max(-0.5, min(lj1 / l01, 1.5));
  isoparams[1] = 1.0 - isoparams[0];

  return outOfBounds;
}

int CIsoparametric::TriangleIsoparameters(const su2double X[][3], const su2double* xj, su2double* isoparams) {
  /*--- The isoparameters are the solution to the determined system X^T * isoparams = xj.
   *    For which we solve the normal equations to avoid divisions by zero.
   *    This is consistent with the shape functions of the linear triangular element. ---*/

  const su2double extrapTol = -0.5;

  /*--- Project the target point onto the triangle. ---*/

  su2double normal[3] = {0.0}, xproj[3] = {0.0};
  TriangleNormal(X, normal);
  PointPlaneProjection<su2double, 3>(xj, X[0], normal, xproj);

  su2double A[3][3] = {{0.0}};  // = X*X^T
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) A[i][j] += X[i][k] * X[j][k];

    isoparams[i] = 0.0;  // use isoparams as rhs
    for (int k = 0; k < 3; ++k) isoparams[i] += X[i][k] * xproj[k];
  }

  /*--- Solve system by in-place Gaussian elimination without pivoting. ---*/

  /*--- Transform system in Upper Matrix. ---*/
  for (int i = 1; i < 3; ++i) {
    for (int j = 0; j < i; ++j) {
      su2double w = A[i][j] / A[j][j];
      for (int k = j; k < 3; ++k) A[i][k] -= w * A[j][k];
      isoparams[i] -= w * isoparams[j];
    }
  }
  /*--- Backwards substitution. ---*/
  for (int i = 2; i >= 0; --i) {
    for (int j = i + 1; j < 3; ++j) isoparams[i] -= A[i][j] * isoparams[j];
    isoparams[i] /= A[i][i];
  }

  /*--- Detect out of bounds point. ---*/
  const int outOfBounds = (isoparams[0] < extrapTol) || (isoparams[1] < extrapTol) || (isoparams[2] < extrapTol);

  /*--- Mitigation. ---*/
  if (outOfBounds) {
    /*--- Clamp. ---*/
    su2double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
      isoparams[i] = max(isoparams[i], extrapTol);
      sum += isoparams[i];
    }
    /*--- Enforce unit sum. ---*/
    for (int i = 0; i < 3; ++i) isoparams[i] /= sum;
  }

  return outOfBounds;
}

int CIsoparametric::QuadrilateralIsoparameters(const su2double X[][3], const su2double* xj, su2double* isoparams) {
  /*--- The isoparameters are the shape functions (Ni) evaluated at xj, for that we need
   *    the corresponding Xi and Eta, which are obtained by solving the overdetermined
   *    nonlinear system r = xj - X^T * Ni(Xi,Eta) = 0 via the modified Marquardt method.
   *    It is not necessary to project the point as minimizing ||r|| is equivalent. ---*/

  const su2double extrapTol = 3.0;

  constexpr int NITER = 20;
  const su2double tol = 1e-10, lambda = 0.05;

  su2double Xi = 0.0, Eta = 0.0, eps;

  /*--- Finding Xi and Eta is a "third order" effect that we do not
   *    differentiate (also because we need to iterate). ---*/

  const bool wasActive = AD::BeginPassive();

  for (int iter = 0; iter < NITER; ++iter) {
    /*--- Evaluate the residual. ---*/
    su2double r[3] = {xj[0], xj[1], xj[2]};
    su2double Ni[4] = {0.0};
    CQUAD4::ShapeFunctions(Xi, Eta, Ni);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 4; ++j) r[i] -= X[j][i] * Ni[j];

    /*--- Evaluate the residual Jacobian. ---*/
    su2double dNi[4][2] = {{0.0}};
    CQUAD4::ShapeFunctionJacobian(Xi, Eta, dNi);

    su2double jac[3][2] = {{0.0}};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 4; ++k) jac[i][j] -= X[k][i] * dNi[k][j];

    /*--- Compute the correction (normal equations and Cramer's rule). ---*/
    su2double A[2][2] = {{0.0}}, b[2] = {0.0};
    for (int i = 0; i < 2; ++i) {
      for (int j = i; j < 2; ++j)
        for (int k = 0; k < 3; ++k) A[i][j] += jac[k][i] * jac[k][j];

      A[i][i] *= (1.0 + lambda);

      for (int k = 0; k < 3; ++k) b[i] += jac[k][i] * r[k];
    }
    A[1][0] = A[0][1];

    su2double detA = 1.0 / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
    su2double dXi = (b[0] * A[1][1] - b[1] * A[0][1]) * detA;
    su2double dEta = (A[0][0] * b[1] - A[1][0] * b[0]) * detA;
    Xi -= dXi;
    Eta -= dEta;

    eps = fabs(dXi) + fabs(dEta);
    if (eps < tol) break;
  }

  AD::EndPassive(wasActive);

  int outOfBounds = 0;

  if (eps > 0.01) {
    /*--- Iteration diverged, hard fallback. ---*/
    Xi = Eta = 0.0;
    outOfBounds = 1;
  } else {
    /*--- Check bounds. ---*/
    outOfBounds = (fabs(Xi) > extrapTol) || (fabs(Eta) > extrapTol);

    /*--- Mitigate by clamping coordinates. ---*/
    if (outOfBounds) {
      Xi = max(-extrapTol, min(Xi, extrapTol));
      Eta = max(-extrapTol, min(Eta, extrapTol));
    }
  }

  /*--- Evaluate isoparameters at final Xi and Eta. ---*/
  CQUAD4::ShapeFunctions(Xi, Eta, isoparams);

  return outOfBounds;
}
