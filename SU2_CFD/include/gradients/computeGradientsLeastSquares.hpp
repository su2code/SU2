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

/*!
 * \brief Caching strategy for the least-squares metric terms (S := inv(A)).
 * \note The metric terms depend only on the grid coordinates and the type of weighting,
 *       they can be computed once (BUILD) and reused (APPLY) while the coordinates do
 *       not change, which reduces the work per evaluation to the right-hand-side sums
 *       and one small matrix-vector product per point. On moving/deforming grids the
 *       cache is invalidated when the dual grid is updated (CGeometry::SetControlVolume)
 *       and rebuilt on the next evaluation.
 */
enum class LSQ_METRIC_CACHE {
  NONE,   /*!< \brief Assemble and factorize the metric terms on the fly (no reuse). */
  BUILD,  /*!< \brief Assemble and factorize, then store S in the geometry cache for reuse. */
  APPLY,  /*!< \brief Reuse the S matrix stored in the geometry cache. */
};

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
                                   GradientType& gradient,
                                   su2activematrix* metricCache = nullptr)
{
  const auto eps = pow(std::numeric_limits<passivedouble>::epsilon(),2);

  /*--- Entries of upper triangular matrix R. ---*/

  if (periodic) {
    AD::StartPreacc();
    AD::SetPreaccIn(Rmatrix(iPoint,0,0));
    AD::SetPreaccIn(Rmatrix(iPoint,0,1));
    AD::SetPreaccIn(Rmatrix(iPoint,1,1));
  }

  su2double r11 = Rmatrix(iPoint,0,0);
  su2double r12 = Rmatrix(iPoint,0,1);
  su2double r22 = Rmatrix(iPoint,1,1);
  su2double r13 = 0.0, r23 = 0.0, r33 = 1.0;

  r11 = sqrt(max(r11, eps));
  r12 /= r11;
  r22 = sqrt(max(r22 - r12*r12, eps));

  if (nDim == 3) {
    if (periodic) {
      AD::SetPreaccIn(Rmatrix(iPoint,0,2));
      AD::SetPreaccIn(Rmatrix(iPoint,1,2));
      AD::SetPreaccIn(Rmatrix(iPoint,2,1));
      AD::SetPreaccIn(Rmatrix(iPoint,2,2));
    }

    r13 = Rmatrix(iPoint,0,2);
    r33 = Rmatrix(iPoint,2,2);
    const auto r23_a = Rmatrix(iPoint,1,2);
    const auto r23_b = Rmatrix(iPoint,2,1);

    r13 /= r11;
    r23 = r23_a/r22 - r23_b*r12/(r11*r22);
    r33 = sqrt(max(r33 - r23*r23 - r13*r13, eps));
  }

  /*--- Compute determinant ---*/

  const su2double detR2 = pow(r11*r22*r33, 2);

  /*--- S matrix := inv(R)*traspose(inv(R)) ---*/

  su2double Smatrix[nDim][nDim] = {{0.0}};

  /*--- Detect singular matrix ---*/

  if (detR2 > eps) {
    computeSmatrix(r11, r12, r13, r22, r23, r33, detR2, Smatrix);
  }

  /*--- Store the metric terms in the geometry cache for reuse (a singular point stores
   *    S = 0, i.e. its gradient will be zero, as in the on-the-fly path). Only used in
   *    primal mode, hence no AD handling. ---*/

  if (metricCache != nullptr) {
    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = iDim; jDim < nDim; ++jDim)
        (*metricCache)(iPoint, lsqCacheIdx(nDim, iDim, jDim)) = Smatrix[iDim][jDim];
  }

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
 * \brief Fast least-squares gradient evaluation reusing cached metric terms.
 * \ingroup FvmAlgos
 * \note Requires a previous call in BUILD mode with the same weighting, which stores
 *       S = inv(A) in the geometry cache (shared by all solvers, one slot per weighting).
 *       Only the right-hand side b = sum_k w*dist*(u_k - u_i) is accumulated (in an edge
 *       loop, i.e. each edge is visited once since its contribution is identical for both
 *       end points), followed by the product S*b per point. Not compatible with periodic
 *       boundaries.
 */
template<size_t nDim, class FieldType, class GradientType>
void computeGradientsLeastSquaresCached(CSolver* solver,
                                        MPI_QUANTITIES kindMpiComm,
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
 * \param[out] Rmatrix - Generic object implementing operator (iPoint, iDim, iDim).
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
                                  LSQ_METRIC_CACHE cacheMode)
{
  const bool periodic = (solver != nullptr) && (config.GetnMarker_Periodic() > 0);

  /*--- The metric caching does not support the mid-computation periodic accumulations. ---*/

  if (periodic) cacheMode = LSQ_METRIC_CACHE::NONE;

  if (cacheMode == LSQ_METRIC_CACHE::APPLY) {
    computeGradientsLeastSquaresCached<nDim>(solver, kindMpiComm, geometry, config, weighted,
                                             field, varBegin, varEnd, idxVel, gradient);
    return;
  }

  const size_t nPointDomain = geometry.GetnPointDomain();

  /*--- In BUILD mode, allocate the geometry cache for this weighting (all threads
   *    synchronize on the master doing the allocation). ---*/

  su2activematrix* metricCache = nullptr;

  if (cacheMode == LSQ_METRIC_CACHE::BUILD) {
    metricCache = &geometry.GetLSQMetricCache(weighted);
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      if (metricCache->size() == 0) metricCache->resize(nPointDomain, nDim*(nDim+1)/2);
    } END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

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

      solveLeastSquares<nDim, false>(iPoint, varBegin, varEnd, Rmatrix, gradient, metricCache);
    }
  }
  END_SU2_OMP_FOR

  /*--- Declare the cache valid and make sure the edge coloring is available before the
   *    first cached evaluation, building it here avoids a race on its lazy construction. ---*/

  if (metricCache != nullptr) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      geometry.GetEdgeColoring();
      geometry.SetLSQMetricCacheValid(weighted);
    } END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

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
                                  LSQ_METRIC_CACHE cacheMode = LSQ_METRIC_CACHE::NONE) {
  switch (geometry.GetnDim()) {
  case 2:
    detail::computeGradientsLeastSquares<2>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                            weighted, field, varBegin, varEnd, idxVel, gradient, Rmatrix, cacheMode);
    break;
  case 3:
    detail::computeGradientsLeastSquares<3>(solver, kindMpiComm, kindPeriodicComm, geometry, config,
                                            weighted, field, varBegin, varEnd, idxVel, gradient, Rmatrix, cacheMode);
    break;
  default:
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
    break;
  }
}
