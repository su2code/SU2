/*!
 * \file computeGradientsLeastSquares.hpp
 * \brief Generic implementation of Least-Squares gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
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

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "correctGradientsSymmetry.hpp"

namespace detail {

/*!
 * \brief Flattened index of entry (iDim,jDim), iDim <= jDim, of an upper triangular
 *        matrix stored row-wise (the layout of the cached LSQ metric terms).
 */
FORCEINLINE constexpr size_t lsqCacheIdx(size_t nDim, size_t iDim, size_t jDim) {
  return iDim * nDim - (iDim * (iDim - 1)) / 2 + (jDim - iDim);
}

/*!
 * \brief Prepare Smatrix for 2D.
 * \ingroup FvmAlgos
 */
FORCEINLINE void computeSmatrix(su2double r11, su2double r12, su2double r13,
                                su2double r22, su2double r23, su2double r33,
                                su2double detR2, su2double Smatrix[][2]) {
  Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
  Smatrix[0][1] = -r11*r12/detR2;
  Smatrix[1][1] = r11*r11/detR2;
}

/*!
 * \brief Prepare Smatrix for 3D.
 * \ingroup FvmAlgos
 */
FORCEINLINE void computeSmatrix(su2double r11, su2double r12, su2double r13,
                                su2double r22, su2double r23, su2double r33,
                                su2double detR2, su2double Smatrix[][3]) {
  su2double z11 = r22*r33;
  su2double z12 =-r12*r33;
  su2double z13 = r12*r23-r13*r22;
  su2double z22 = r11*r33;
  su2double z23 =-r11*r23;
  su2double z33 = r11*r22;

  Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
  Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
  Smatrix[0][2] = (z13*z33)/detR2;
  Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
  Smatrix[1][2] = (z23*z33)/detR2;
  Smatrix[2][2] = (z33*z33)/detR2;
}

/*!
 * \brief Factorize the accumulated normal matrix A of the least-squares problem
 *        (Cholesky) and form S = inv(A), the entries r* are the unique entries of A.
 *        A (nearly) singular matrix results in S = 0, i.e. a zero gradient.
 * \ingroup FvmAlgos
 */
template<size_t nDim>
FORCEINLINE void invertNormalMatrix(su2double r11, su2double r12, su2double r13,
                                    su2double r22, su2double r23_a, su2double r23_b,
                                    su2double r33, su2double Smatrix[][nDim])
{
  const auto eps = pow(std::numeric_limits<passivedouble>::epsilon(),2);

  r11 = sqrt(max(r11, eps));
  r12 /= r11;
  r22 = sqrt(max(r22 - r12*r12, eps));

  su2double r23 = 0.0;
  if (nDim == 3) {
    r13 /= r11;
    r23 = r23_a/r22 - r23_b*r12/(r11*r22);
    r33 = sqrt(max(r33 - r23*r23 - r13*r13, eps));
  }
  else {
    r13 = 0.0;
    r33 = 1.0;
  }

  /*--- Compute determinant ---*/

  const su2double detR2 = pow(r11*r22*r33, 2);

  /*--- S matrix := inv(R)*traspose(inv(R)), detect singular matrix ---*/

  if (detR2 > eps) {
    computeSmatrix(r11, r12, r13, r22, r23, r33, detR2, Smatrix);
  }
}

/*!
 * \brief Solve the least-squares problem for one point.
 * \ingroup FvmAlgos
 * \note See detail::computeGradientsLeastSquares for the
 *       purpose of template "nDim" and "periodic".
 */
template<size_t nDim, bool periodic, class GradientType, class RMatrixType>
FORCEINLINE void solveLeastSquares(size_t iPoint,
                                   size_t varBegin,
                                   size_t varEnd,
                                   const RMatrixType& Rmatrix,
                                   GradientType& gradient)
{
  /*--- Entries of the normal matrix A. ---*/

  if (periodic) {
    AD::StartPreacc();
    AD::SetPreaccIn(Rmatrix(iPoint,0,0));
    AD::SetPreaccIn(Rmatrix(iPoint,0,1));
    AD::SetPreaccIn(Rmatrix(iPoint,1,1));
  }

  const su2double r11 = Rmatrix(iPoint,0,0);
  const su2double r12 = Rmatrix(iPoint,0,1);
  const su2double r22 = Rmatrix(iPoint,1,1);
  su2double r13 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;

  if (nDim == 3) {
    if (periodic) {
      AD::SetPreaccIn(Rmatrix(iPoint,0,2));
      AD::SetPreaccIn(Rmatrix(iPoint,1,2));
      AD::SetPreaccIn(Rmatrix(iPoint,2,1));
      AD::SetPreaccIn(Rmatrix(iPoint,2,2));
    }

    r13 = Rmatrix(iPoint,0,2);
    r33 = Rmatrix(iPoint,2,2);
    r23_a = Rmatrix(iPoint,1,2);
    r23_b = Rmatrix(iPoint,2,1);
  }

  su2double Smatrix[nDim][nDim] = {{0.0}};

  invertNormalMatrix<nDim>(r11, r12, r13, r22, r23_a, r23_b, r33, Smatrix);

  if (periodic) {
    /*--- Stop preacc here as gradient is in/out. ---*/
    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = iDim; jDim < nDim; ++jDim)
        AD::SetPreaccOut(Smatrix[iDim][jDim]);
    AD::EndPreacc();
  }

  /*--- Computation of the gradient: S*c ---*/

  for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
  {
    su2double Cvector[nDim] = {0.0};

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = 0; jDim < nDim; ++jDim)
        Cvector[iDim] += Smatrix[min(iDim,jDim)][max(iDim,jDim)] * gradient(iPoint, iVar, jDim);

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      gradient(iPoint, iVar, iDim) = Cvector[iDim];
  }

  if (!periodic) {
    /*--- Stop preacc here instead as gradient is only out. ---*/
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        AD::SetPreaccOut(gradient(iPoint, iVar, iDim));
    AD::EndPreacc();
  }
}

/*!
 * \brief Rotation matrix of a periodic marker (rotation about x, then y, then z axis).
 * \ingroup FvmAlgos
 */
FORCEINLINE void periodicRotationMatrix(const su2double* angles, su2double rotMatrix[][2]) {
  GeometryToolbox::RotationMatrix(angles[2], rotMatrix);
}
FORCEINLINE void periodicRotationMatrix(const su2double* angles, su2double rotMatrix[][3]) {
  GeometryToolbox::RotationMatrix(angles[0], angles[1], angles[2], rotMatrix);
}

/*!
 * \brief Add the contributions of the neighbors across the periodic boundaries to the
 *        (not yet factorized) least-squares normal matrices of the periodic points.
 * \ingroup FvmAlgos
 * \note This is the grid-dependent part of the least-squares cases of CSolver::Initiate/
 *       CompletePeriodicComms (which, with cached metrics, exchange only the right-hand
 *       sides). It needs only the geometry, owner of the periodic communication framework,
 *       and the config (rotation/translation of the markers), hence the cached metrics can
 *       be completed during the geometry preprocessing. Like the solver comms, the exchange
 *       is done once per pair of periodic markers and the points on the current pair are
 *       updated, the sends always cover all pairs. Must be called by all threads if inside
 *       an OpenMP parallel region.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem.
 * \param[in] weighted - Use inverse-distance weights.
 * \param[in,out] normalMatrix - Unique entries of the normal matrix of each point (upper
 *                triangle row-wise, see lsqCacheIdx), the periodic terms are added to it.
 */
template<size_t nDim>
void addPeriodicLSQMetricTerms(CGeometry& geometry, const CConfig& config, bool weighted,
                               su2activematrix& normalMatrix)
{
  constexpr unsigned short nEntries = nDim*(nDim+1)/2;
  const auto nPeriodic = config.GetnMarker_Periodic();

  /*--- Make sure the communication buffers are large enough. ---*/

  geometry.AllocatePeriodicComms(nEntries);

  /*--- Status is global so all threads can see the result of Waitany. ---*/

  static SU2_MPI::Status status;

  for (unsigned short iPair = 1; iPair <= nPeriodic/2; ++iPair) {

    /*--- Load the periodic terms of the owned periodic points into the send buffers. ---*/

    if (geometry.nPeriodicSend > 0) {

      geometry.PostPeriodicRecvs(&geometry, &config, COMM_TYPE::DOUBLE, nEntries);

      for (int iMessage = 0; iMessage < geometry.nPeriodicSend; ++iMessage) {

        const auto msg_offset = geometry.nPoint_PeriodicSend[iMessage];
        const auto nSend = geometry.nPoint_PeriodicSend[iMessage+1] - msg_offset;

        SU2_OMP_FOR_STAT(32)
        for (int iSend = 0; iSend < nSend; ++iSend) {

          const auto iPoint = geometry.Local_Point_PeriodicSend[msg_offset + iSend];
          const auto iMarker = geometry.Local_Marker_PeriodicSend[msg_offset + iSend];

          /*--- Transformation of the marker, the point and its neighbors are moved to
           *    their location on the donor marker (rotation about the center followed by
           *    translation). ---*/

          const auto& markerTag = config.GetMarker_All_TagBound(iMarker);
          const su2double* center = config.GetPeriodicRotCenter(markerTag);
          const su2double* angles = config.GetPeriodicRotAngles(markerTag);
          const su2double* trans = config.GetPeriodicTranslation(markerTag);

          su2double rotMatrix[nDim][nDim] = {{0.0}};
          periodicRotationMatrix(angles, rotMatrix);

          su2double translation[nDim] = {0.0};
          for (size_t iDim = 0; iDim < nDim; ++iDim) translation[iDim] = center[iDim] + trans[iDim];

          auto transform = [&](const su2double* coord, su2double* rotCoord) {
            su2double distance[nDim] = {0.0};
            GeometryToolbox::Distance(nDim, coord, center, distance);
            GeometryToolbox::Rotate(rotMatrix, translation, distance, rotCoord);
          };

          su2double rotCoord_i[nDim] = {0.0};
          transform(geometry.nodes->GetCoord(iPoint), rotCoord_i);

          /*--- Accumulate the normal matrix terms of the transformed neighbors. ---*/

          su2double A[nEntries] = {0.0};

          for (auto jPoint : geometry.nodes->GetPoints(iPoint)) {

            /*--- Avoid periodic boundary points so that we do not
             *    duplicate edges on both sides of the periodic BC. ---*/

            if (geometry.nodes->GetPeriodicBoundary(jPoint)) continue;

            su2double rotCoord_j[nDim] = {0.0};
            transform(geometry.nodes->GetCoord(jPoint), rotCoord_j);

            su2double dist_ij[nDim] = {0.0};
            GeometryToolbox::Distance(nDim, rotCoord_j, rotCoord_i, dist_ij);

            su2double weight = 1.0;
            if (weighted) {
              weight = GeometryToolbox::SquaredNorm(nDim, dist_ij);
              if (weight == 0.0) continue;
            }

            for (size_t iDim = 0; iDim < nDim; ++iDim)
              for (size_t jDim = iDim; jDim < nDim; ++jDim)
                A[lsqCacheIdx(nDim, iDim, jDim)] += dist_ij[iDim]*dist_ij[jDim]/weight;
          }

          const auto buf_offset = (msg_offset + iSend)*nEntries;

          for (unsigned short k = 0; k < nEntries; ++k)
            geometry.bufD_PeriodicSend[buf_offset + k] = A[k];
        }
        END_SU2_OMP_FOR

        geometry.PostPeriodicSends(&geometry, &config, COMM_TYPE::DOUBLE, nEntries, iMessage);
      }
    }

    /*--- Accumulate the received terms in the normal matrices of the points on the
     *    current pair of periodic markers. ---*/

    if (geometry.nPeriodicRecv > 0) {

      for (int iMessage = 0; iMessage < geometry.nPeriodicRecv; ++iMessage) {

        /*--- Receive the messages dynamically based on the order they arrive. ---*/

        int source = 0;
#ifdef HAVE_MPI
        int ind = 0;
        SU2_OMP_SAFE_GLOBAL_ACCESS(SU2_MPI::Waitany(geometry.nPeriodicRecv, geometry.req_PeriodicRecv, &ind, &status);)
        source = status.MPI_SOURCE;
#else
        source = SU2_MPI::GetRank();
        SU2_OMP_BARRIER
#endif
        const auto jRecv = geometry.PeriodicRecv2Neighbor.at(source);
        const auto msg_offset = geometry.nPoint_PeriodicRecv[jRecv];
        const auto nRecv = geometry.nPoint_PeriodicRecv[jRecv+1] - msg_offset;

        SU2_OMP_FOR_STAT(32)
        for (int iRecv = 0; iRecv < nRecv; ++iRecv) {

          const auto iPoint = geometry.Local_Point_PeriodicRecv[msg_offset + iRecv];
          const auto iPeriodic = geometry.Local_Marker_PeriodicRecv[msg_offset + iRecv];

          if ((iPeriodic != iPair) && (iPeriodic != iPair + nPeriodic/2)) continue;

          const auto buf_offset = (msg_offset + iRecv)*nEntries;

          for (unsigned short k = 0; k < nEntries; ++k)
            normalMatrix(iPoint, k) += geometry.bufD_PeriodicRecv[buf_offset + k];
        }
        END_SU2_OMP_FOR
      }

      /*--- Verify that all non-blocking sends have finished. ---*/

#ifdef HAVE_MPI
      SU2_OMP_SAFE_GLOBAL_ACCESS(SU2_MPI::Waitall(geometry.nPeriodicSend, geometry.req_PeriodicSend, MPI_STATUS_IGNORE);)
#endif
    }
  }
}

/*!
 * \brief Assemble, factorize, and store the least-squares gradient metric terms
 *        (S = inv(A), upper triangle row-wise) of a grid in the geometry cache.
 * \ingroup FvmAlgos
 * \note The metric terms depend only on the node coordinates and the weighting. They are
 *       computed during the geometry preprocessing (see CDriver::InitializeGeometry) and,
 *       on moving/deforming grids, recomputed on the first gradient evaluation after the
 *       dual grid update invalidates them (CGeometry::SetControlVolume). On grids with
 *       periodic boundaries the normal matrix also accumulates the contributions of the
 *       periodic neighbors (see addPeriodicLSQMetricTerms). The function is safe to call
 *       from inside or outside an OpenMP parallel region (by all threads), and returns
 *       immediately if the cache is already valid.
 * \param[in] geometry - Geometric grid properties, owner of the cache.
 * \param[in] config - Configuration of the problem.
 * \param[in] weighted - Use inverse-distance weights.
 */
template<size_t nDim>
void computeLSQMetrics(CGeometry& geometry, const CConfig& config, bool weighted)
{
  if (geometry.LSQMetricCacheIsValid(weighted)) return;

  const size_t nPoint = geometry.GetnPoint();
  const size_t nPointDomain = geometry.GetnPointDomain();
  constexpr size_t nEntries = nDim*(nDim+1)/2;

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const size_t chunkSize = computeStaticChunkSize(nPointDomain,
                           omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  auto& metricCache = geometry.GetLSQMetricCache(weighted);

  /*--- The periodic comms may deliver contributions to halo points, hence the size nPoint. ---*/

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    if (metricCache.size() == 0) metricCache.resize(nPoint, nEntries);
  } END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Accumulate the unique entries of the normal matrix A (upper triangle, row-wise)
   *    in the cache, the factorization is done in place after the periodic exchange. ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    const auto coord_i = geometry.nodes->GetCoord(iPoint);

    su2double A[nEntries] = {0.0};

    for (auto jPoint : geometry.nodes->GetPoints(iPoint)) {
      su2double dist_ij[nDim] = {0.0};
      GeometryToolbox::Distance(nDim, geometry.nodes->GetCoord(jPoint), coord_i, dist_ij);

      su2double weight = 1.0;
      if (weighted) {
        const su2double dist2 = GeometryToolbox::SquaredNorm(nDim, dist_ij);
        if (dist2 <= 0.0) continue;
        weight = 1.0 / dist2;
      }

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        for (size_t jDim = iDim; jDim < nDim; ++jDim)
          A[lsqCacheIdx(nDim, iDim, jDim)] += dist_ij[iDim]*dist_ij[jDim]*weight;
    }

    for (size_t k = 0; k < nEntries; ++k) metricCache(iPoint, k) = A[k];
  }
  END_SU2_OMP_FOR

  /*--- Add the contributions of the neighbors across the periodic boundaries. ---*/

  if (config.GetnMarker_Periodic() > 0) {
    addPeriodicLSQMetricTerms<nDim>(geometry, config, weighted, metricCache);
  }

  /*--- Factorize and invert the normal matrix in place, S = inv(A). ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    const su2double r11 = metricCache(iPoint, lsqCacheIdx(nDim, 0, 0));
    const su2double r12 = metricCache(iPoint, lsqCacheIdx(nDim, 0, 1));
    const su2double r22 = metricCache(iPoint, lsqCacheIdx(nDim, 1, 1));
    su2double r13 = 0.0, r23 = 0.0, r33 = 0.0;

    if (nDim == 3) {
      r13 = metricCache(iPoint, lsqCacheIdx(nDim, 0, nDim-1));
      r23 = metricCache(iPoint, lsqCacheIdx(nDim, 1, nDim-1));
      r33 = metricCache(iPoint, lsqCacheIdx(nDim, nDim-1, nDim-1));
    }

    su2double Smatrix[nDim][nDim] = {{0.0}};

    /*--- A is symmetric, its (2,1) entry is r13. ---*/

    invertNormalMatrix<nDim>(r11, r12, r13, r22, r23, r13, r33, Smatrix);

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = iDim; jDim < nDim; ++jDim)
        metricCache(iPoint, lsqCacheIdx(nDim, iDim, jDim)) = Smatrix[iDim][jDim];
  }
  END_SU2_OMP_FOR

  /*--- Declare the cache valid and make sure the edge coloring is available before the
   *    first cached evaluation, building it here avoids a race on its lazy construction. ---*/

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    geometry.GetEdgeColoring();
    geometry.SetLSQMetricCacheValid(weighted);
  } END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

/*!
 * \brief Fast least-squares gradient evaluation reusing the cached metric terms.
 * \ingroup FvmAlgos
 * \note Requires valid metric terms for this weighting in the geometry cache (see
 *       computeLSQMetrics, shared by all solvers, one slot per weighting). Only the
 *       right-hand side b = sum_k w*dist*(u_k - u_i) is accumulated (in an edge loop,
 *       i.e. each edge is visited once since its contribution is identical for both
 *       end points), followed by the product S*b per point. On grids with periodic
 *       boundaries the right-hand sides are completed with the contributions of the
 *       periodic neighbors, the metric terms of the cache already include them so the
 *       periodic communications must exchange only the right-hand sides (which is what
 *       they do when CConfig::GetLSQMetricCaching() is true).
 */
template<size_t nDim, class FieldType, class GradientType>
void computeGradientsLeastSquaresCached(CSolver* solver,
                                        MPI_QUANTITIES kindMpiComm,
                                        PERIODIC_QUANTITIES kindPeriodicComm,
                                        CGeometry& geometry,
                                        const CConfig& config,
                                        bool weighted,
                                        const FieldType& field,
                                        const size_t varBegin,
                                        const size_t varEnd,
                                        const int idxVel,
                                        GradientType& gradient)
{
  const auto& metricCache = geometry.GetLSQMetricCache(weighted);
  const size_t nPoint = geometry.GetnPoint();
  const size_t nPointDomain = geometry.GetnPointDomain();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const size_t chunkSize = computeStaticChunkSize(nPointDomain,
                           omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- Clear the right-hand-side accumulators, including halo points, which
   *    receive edge contributions (discarded when halos are communicated). ---*/

  SU2_OMP_FOR_STAT(2048)
  for (size_t iPoint = 0; iPoint < nPoint; ++iPoint)
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = 0.0;
  END_SU2_OMP_FOR

  /*--- Accumulate the RHS in a loop over the edges: the contribution of edge {i,j} is
   *    w*dist_ij*(u_j - u_i) for BOTH end points. A race-free edge coloring is required
   *    with multiple threads, the "natural" coloring (single color, used with the
   *    reducer strategy) forces the fallback to a thread-safe loop over nodes. ---*/

  const auto& coloring = geometry.GetEdgeColoring();

  const bool safeColoring = (omp_get_max_threads() == 1) || (coloring.getOuterSize() > 1);

  if (safeColoring) {
    const size_t groupSize = geometry.GetEdgeColorGroupSize();

    for (auto iColor = 0ul; iColor < coloring.getOuterSize(); ++iColor) {
      const auto* edgeIndices = coloring.innerIdx(iColor);
      const auto nEdgesColor = coloring.getNumNonZeros(iColor);

      SU2_OMP_FOR_DYN(nextMultiple(size_t(32), groupSize))
      for (auto k = 0ul; k < nEdgesColor; ++k) {
        const auto iEdge = edgeIndices[k];
        const auto iPoint = geometry.edges->GetNode(iEdge, 0);
        const auto jPoint = geometry.edges->GetNode(iEdge, 1);

        su2double dist_ij[nDim] = {0.0};
        GeometryToolbox::Distance(nDim, geometry.nodes->GetCoord(jPoint),
                                  geometry.nodes->GetCoord(iPoint), dist_ij);

        su2double weight = 1.0;
        if (weighted) {
          const su2double dist2 = GeometryToolbox::SquaredNorm(nDim, dist_ij);
          if (dist2 <= 0.0) continue;
          weight = 1.0 / dist2;
        }

        for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
          const su2double delta_ij = weight * (field(jPoint,iVar) - field(iPoint,iVar));

          for (size_t iDim = 0; iDim < nDim; ++iDim) {
            const su2double contrib = dist_ij[iDim] * delta_ij;
            gradient(iPoint, iVar, iDim) += contrib;
            gradient(jPoint, iVar, iDim) += contrib;
          }
        }
      }
      END_SU2_OMP_FOR
    }
  }
  else {
    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      const auto coord_i = geometry.nodes->GetCoord(iPoint);

      for (auto jPoint : geometry.nodes->GetPoints(iPoint)) {
        su2double dist_ij[nDim] = {0.0};
        GeometryToolbox::Distance(nDim, geometry.nodes->GetCoord(jPoint), coord_i, dist_ij);

        su2double weight = 1.0;
        if (weighted) {
          const su2double dist2 = GeometryToolbox::SquaredNorm(nDim, dist_ij);
          if (dist2 <= 0.0) continue;
          weight = 1.0 / dist2;
        }

        for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
          const su2double delta_ij = weight * (field(jPoint,iVar) - field(iPoint,iVar));

          for (size_t iDim = 0; iDim < nDim; ++iDim)
            gradient(iPoint, iVar, iDim) += dist_ij[iDim] * delta_ij;
        }
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Complete the RHS with the contributions from across the periodic boundaries. ---*/

  if ((solver != nullptr) && (config.GetnMarker_Periodic() > 0)) {
    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic) {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    }
  }

  /*--- Multiply the RHS by the cached S matrix. ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
      su2double Cvector[nDim] = {0.0};

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        for (size_t jDim = 0; jDim < nDim; ++jDim)
          Cvector[iDim] += metricCache(iPoint, lsqCacheIdx(nDim, min(iDim,jDim), max(iDim,jDim))) *
                           gradient(iPoint, iVar, jDim);

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = Cvector[iDim];
    }
  }
  END_SU2_OMP_FOR

  /*--- Compute the corrections for symmetry planes and Euler walls. ---*/

  correctGradientsSymmetry<nDim>(geometry, config, varBegin, varEnd, idxVel, gradient);

  /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

  if (solver != nullptr) {
    solver->InitiateComms(&geometry, &config, kindMpiComm);
    solver->CompleteComms(&geometry, &config, kindMpiComm);
  }
}

/*!
 * \brief Compute the gradient of a field using inverse-distance-weighted or
 *        unweighted Least-Squares approximation.
 * \ingroup FvmAlgos
 * \note See notes from computeGradientsGreenGauss.hpp.
 * \param[in] solver - Optional, solver associated with the field (used only for MPI).
 * \param[in] kindMpiComm - Type of MPI communication required.
 * \param[in] kindPeriodicComm - Type of periodic communication required.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] weighted - Use inverse-distance weights.
 * \param[in] config - Configuration of the problem, used to identify types of boundaries.
 * \param[in] field - Generic object implementing operator (iPoint, iVar).
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[in] idxVel - Index to velocity, -1 if no velocity is present in the solver.
 * \param[out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 * \param[out] Rmatrix - Generic object implementing operator (iPoint, iDim, iDim), not used
 *             (and not written) with the cached metric terms.
 * \param[in] useCaching - Use the cached metric terms (see computeLSQMetrics). On grids with
 *            periodic boundaries the format of the periodic communications depends on the
 *            caching, hence it is governed by CConfig::GetLSQMetricCaching() instead.
 */
template<size_t nDim, class FieldType, class GradientType, class RMatrixType>
void computeGradientsLeastSquares(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  const CConfig& config,
                                  bool weighted,
                                  const FieldType& field,
                                  const size_t varBegin,
                                  const size_t varEnd,
                                  const int idxVel,
                                  GradientType& gradient,
                                  RMatrixType& Rmatrix,
                                  bool useCaching)
{
  const bool periodic = (solver != nullptr) && (config.GetnMarker_Periodic() > 0);

  /*--- Use the cached metric terms, rebuilding them if the coordinates changed since the
   *    geometry preprocessing (or if this combination was not covered by it). On grids with
   *    periodic boundaries the cached metrics include the periodic contributions and the
   *    periodic comms of the solver exchange only the LSQ right-hand sides (as dictated by
   *    the config), the caching is thus only usable if those comms are performed. ---*/

  bool useCachedMetrics = useCaching;

  if (config.GetnMarker_Periodic() > 0) {
    useCachedMetrics = config.GetLSQMetricCaching() && periodic && (kindPeriodicComm != PERIODIC_NONE);
  }

  if (useCachedMetrics) {
    computeLSQMetrics<nDim>(geometry, config, weighted);
    computeGradientsLeastSquaresCached<nDim>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                             weighted, field, varBegin, varEnd, idxVel, gradient);
    return;
  }

  const size_t nPointDomain = geometry.GetnPointDomain();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  size_t chunkSize = computeStaticChunkSize(nPointDomain,
                     omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- First loop over non-halo points of the grid. ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
  {
    auto nodes = geometry.nodes;
    const auto coord_i = nodes->GetCoord(iPoint);

    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();
    AD::SetPreaccIn(coord_i, nDim);

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      AD::SetPreaccIn(field(iPoint,iVar));

    /*--- Clear gradient and Rmatrix. ---*/

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        gradient(iPoint, iVar, iDim) = 0.0;

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = 0; jDim < nDim; ++jDim)
        Rmatrix(iPoint, iDim, jDim) = 0.0;


    for (auto jPoint : nodes->GetPoints(iPoint))
    {
      const auto coord_j = geometry.nodes->GetCoord(jPoint);
      AD::SetPreaccIn(coord_j, nDim);


      /*--- Distance vector from iPoint to jPoint ---*/

      su2double dist_ij[nDim] = {0.0};
      GeometryToolbox::Distance(nDim, coord_j, coord_i, dist_ij);


      /*--- Compute inverse weight, default 1 (unweighted). ---*/

      su2double weight = 1.0;
      if(weighted) weight = GeometryToolbox::SquaredNorm(nDim, dist_ij);

      /*--- Summations for entries of upper triangular matrix R. ---*/

      if (weight > 0.0)
      {
        weight = 1.0 / weight;

        for (size_t iDim = 0; iDim < nDim; ++iDim)
          for (size_t jDim = iDim; jDim < nDim; ++jDim)
            Rmatrix(iPoint,iDim,jDim) += dist_ij[iDim]*dist_ij[jDim]*weight;

        if (nDim == 3)
          Rmatrix(iPoint,2,1) += dist_ij[0]*dist_ij[nDim-1]*weight;

        /*--- Entries of c:= transpose(A)*b ---*/

        for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        {
          AD::SetPreaccIn(field(jPoint,iVar));

          su2double delta_ij = weight * (field(jPoint,iVar) - field(iPoint,iVar));

          for (size_t iDim = 0; iDim < nDim; ++iDim)
            gradient(iPoint, iVar, iDim) += dist_ij[iDim] * delta_ij;
        }
      }
    }

    if (periodic)
    {
      /*--- A second loop is required after periodic comms, checkpoint the preacc. ---*/

      for (size_t iDim = 0; iDim < nDim; ++iDim)
        for (size_t jDim = 0; jDim < nDim; ++jDim)
          AD::SetPreaccOut(Rmatrix(iPoint, iDim, jDim));

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        for (size_t iDim = 0; iDim < nDim; ++iDim)
          AD::SetPreaccOut(gradient(iPoint, iVar, iDim));

      AD::EndPreacc();
    }
    else {
      /*--- Periodic comms are not needed, solve the LS problem for iPoint. ---*/

      solveLeastSquares<nDim, false>(iPoint, varBegin, varEnd, Rmatrix, gradient);
    }
  }
  END_SU2_OMP_FOR

  /*--- Correct the gradient values across any periodic boundaries. ---*/

  if (periodic)
  {
    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
    {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    }

    /*--- Second loop over points of the grid to compute final gradient. ---*/

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
      solveLeastSquares<nDim, true>(iPoint, varBegin, varEnd, Rmatrix, gradient);
    END_SU2_OMP_FOR
  }

  /* --- compute the corrections for symmetry planes and Euler walls. --- */

  correctGradientsSymmetry<nDim>(geometry, config, varBegin, varEnd, idxVel, gradient);

  /*--- If no solver was provided we do not communicate ---*/

  if (solver != nullptr)
  {
    /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

    solver->InitiateComms(&geometry, &config, kindMpiComm);
    solver->CompleteComms(&geometry, &config, kindMpiComm);
  }

}
} // end namespace

/*!
 * \brief Instantiations for 2D and 3D.
 * \ingroup FvmAlgos
 */
template<class FieldType, class GradientType, class RMatrixType>
void computeGradientsLeastSquares(CSolver* solver,
                                  MPI_QUANTITIES kindMpiComm,
                                  PERIODIC_QUANTITIES kindPeriodicComm,
                                  CGeometry& geometry,
                                  const CConfig& config,
                                  bool weighted,
                                  const FieldType& field,
                                  const size_t varBegin,
                                  const size_t varEnd,
                                  const int idxVel,
                                  GradientType& gradient,
                                  RMatrixType& Rmatrix,
                                  bool useCaching = false) {
  switch (geometry.GetnDim()) {
  case 2:
    detail::computeGradientsLeastSquares<2>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                            weighted, field, varBegin, varEnd, idxVel, gradient, Rmatrix, useCaching);
    break;
  case 3:
    detail::computeGradientsLeastSquares<3>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                            weighted, field, varBegin, varEnd, idxVel, gradient, Rmatrix, useCaching);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}

/*!
 * \brief Compute (if not already valid) the cached least-squares gradient metric terms of
 *        a grid for one type of weighting, see detail::computeLSQMetrics.
 * \ingroup FvmAlgos
 */
inline void computeLSQGradientMetrics(CGeometry& geometry, const CConfig& config, bool weighted) {
  switch (geometry.GetnDim()) {
  case 2:
    detail::computeLSQMetrics<2>(geometry, config, weighted);
    break;
  case 3:
    detail::computeLSQMetrics<3>(geometry, config, weighted);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}
