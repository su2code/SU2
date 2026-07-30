/*!
 * \file computeGradientsGreenGauss.hpp
 * \brief Generic implementation of Green-Gauss gradient computation.
 * \note This allows the same implementation to be used for conservative
 *       and primitive variables of any solver.
 * \author P. Gomes
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

#pragma once

#include <vector>
#include <algorithm>

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "correctGradientsSymmetry.hpp"
#include "../limiters/fieldMinMax.hpp"

namespace detail {

/*!
 * \brief Accumulate the Green-Gauss contributions of the edges around one point.
 * \ingroup FvmAlgos
 * \note With Difference=false the usual (f_i+f_j)/2 face value is used, which requires the
 *       boundary faces of the control volume to be accounted for separately (see
 *       greenGaussBoundaryTerms). With Difference=true the face value is (f_j-f_i)/2, which
 *       is equivalent because the normals of a closed control volume sum to zero, and has
 *       the advantage of making the boundary contributions vanish identically.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] iPoint - Point for which to compute the gradient.
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[in] field - Generic object implementing operator (iPoint, iVar).
 * \param[out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 * \param[in,out] extrema - Optionally track the min/max of the field over the neighbors, which
 *                comes for free here since their values are needed for the gradient anyway
 *                (see the limited version of this method). Pass CNoMinMax to skip it.
 */
template <size_t nDim, bool Difference, class FieldType, class GradientType, class MinMaxType = CNoMinMax>
FORCEINLINE void greenGaussPointContributions(const CGeometry& geometry, size_t iPoint, size_t varBegin,
                                              size_t varEnd, const FieldType& field, GradientType& gradient,
                                              MinMaxType extrema = {}) {
  const auto nodes = geometry.nodes;

  AD::SetPreaccIn(nodes->GetVolume(iPoint));
  AD::SetPreaccIn(nodes->GetPeriodicVolume(iPoint));

  for (size_t iVar = varBegin; iVar < varEnd; ++iVar) AD::SetPreaccIn(field(iPoint, iVar));

  /*--- Clear the gradient. --*/

  for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
    for (size_t iDim = 0; iDim < nDim; ++iDim) gradient(iPoint, iVar, iDim) = 0.0;

  /*--- Handle averaging and division by volume in one constant. ---*/

  const su2double halfOnVol = 0.5 / (nodes->GetVolume(iPoint) + nodes->GetPeriodicVolume(iPoint));

  /*--- Add a contribution due to each neighbor. ---*/

  for (size_t iNeigh = 0; iNeigh < nodes->GetnPoint(iPoint); ++iNeigh) {
    const size_t iEdge = nodes->GetEdge(iPoint, iNeigh);
    const size_t jPoint = nodes->GetPoint(iPoint, iNeigh);

    /*--- Determine if edge points inwards or outwards of iPoint.
     *    If inwards we need to flip the area vector. ---*/

    const su2double dir = (iPoint < jPoint) ? 1.0 : -1.0;
    const su2double weight = dir * halfOnVol;

    const auto area = geometry.edges->GetNormal(iEdge);
    AD::SetPreaccIn(area, nDim);

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
      AD::SetPreaccIn(field(jPoint, iVar));
      const su2double faceValue =
          Difference ? (field(jPoint, iVar) - field(iPoint, iVar)) : (field(jPoint, iVar) + field(iPoint, iVar));
      const su2double flux = weight * faceValue;

      for (size_t iDim = 0; iDim < nDim; ++iDim) gradient(iPoint, iVar, iDim) += flux * area[iDim];

      extrema.update(iVar, field(jPoint, iVar));
    }
  }
}

/*!
 * \brief Add the contributions of the boundary faces of the control volumes, which close the
 *        Green-Gauss integration when the (f_i+f_j)/2 face value is used for interior faces.
 * \ingroup FvmAlgos
 * \note Markers that do not correspond to actual boundaries of the domain are skipped, the
 *       points on those markers see their neighbors across the marker as regular (halo) points.
 */
template <size_t nDim, class FieldType, class GradientType>
void greenGaussBoundaryTerms(CGeometry& geometry, const CConfig& config, size_t varBegin, size_t varEnd,
                             const FieldType& field, GradientType& gradient) {
  static constexpr size_t MAXNVAR = 20;
  su2double flux[MAXNVAR] = {0.0};

  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if ((config.GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config.GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
        (config.GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      /*--- Work is shared in inner loop as two markers
       *    may try to update the same point. ---*/

      SU2_OMP_FOR_STAT(32)
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
        auto nodes = geometry.nodes;

        /*--- Halo points do not need to be considered. ---*/

        if (!nodes->GetDomain(iPoint)) continue;

        su2double volume = nodes->GetVolume(iPoint) + nodes->GetPeriodicVolume(iPoint);
        const auto area = geometry.vertex[iMarker][iVertex]->GetNormal();

        for (size_t iVar = varBegin; iVar < varEnd; iVar++)
          flux[iVar] = field(iPoint,iVar) / volume;

        for (size_t iVar = varBegin; iVar < varEnd; iVar++) {
          for (size_t iDim = 0; iDim < nDim; iDim++) {
            gradient(iPoint, iVar, iDim) -= flux[iVar] * area[iDim];
          }
        } // loop over variables
      } // vertices
      END_SU2_OMP_FOR
    } //found right marker
  } // iMarkers
}

/*!
 * \brief Compute the gradient of a field using the Green-Gauss theorem.
 * \ingroup FvmAlgos
 * \note Template nDim to allow efficient unrolling of inner loops.
 * \note Gradients can be computed only for a contiguous range of variables, defined
 *       by [varBegin, varEnd[ (e.g. 0,1 computes the gradient of the 1st variable).
 *       This can be used, for example, to compute only velocity gradients.
 * \note The function uses an optional solver object to perform communications, if
 *       none (nullptr) is provided the function does not fail (the objective of
 *       this is to improve test-ability).
 * \param[in] solver - Optional, solver associated with the field (used only for MPI).
 * \param[in] kindMpiComm - Type of MPI communication required.
 * \param[in] kindPeriodicComm - Type of periodic communication required.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem, used to identify types of boundaries.
 * \param[in] field - Generic object implementing operator (iPoint, iVar).
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[in] idxVel - Index of velocity, or -1 if no velocity present.
 * \param[out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 */
template <size_t nDim, class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver, MPI_QUANTITIES kindMpiComm, PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry, const CConfig& config, const FieldType& field,
                                const size_t varBegin, const size_t varEnd, const int idxVel, GradientType& gradient) {
  const size_t nPointDomain = geometry.GetnPointDomain();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- For each (non-halo) volume integrate over its faces (edges). ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();

    greenGaussPointContributions<nDim, false>(geometry, iPoint, varBegin, varEnd, field, gradient);

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim) AD::SetPreaccOut(gradient(iPoint, iVar, iDim));

    AD::EndPreacc();
  }
  END_SU2_OMP_FOR

  /*--- Add edges of markers that contribute to the gradients ---*/

  greenGaussBoundaryTerms<nDim>(geometry, config, varBegin, varEnd, field, gradient);

  /*--- Compute the corrections for symmetry planes and Euler walls. ---*/

  correctGradientsSymmetry<nDim>(geometry, config, varBegin, varEnd, idxVel, gradient);

  /*--- If no solver was provided we do not communicate ---*/

  if (solver == nullptr) return;

  /*--- Account for periodic contributions. ---*/

  for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic() / 2; ++iPeriodic) {
    solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
  }

  /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);
}
}  // namespace detail



/*!
 * \brief Instantiations for 2D and 3D.
 * \ingroup FvmAlgos
 */
template <class FieldType, class GradientType>
void computeGradientsGreenGauss(CSolver* solver, MPI_QUANTITIES kindMpiComm, PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry, const CConfig& config, const FieldType& field,
                                const size_t varBegin, const size_t varEnd, const int idxVel, GradientType& gradient) {
  switch (geometry.GetnDim()) {
    case 2:
      detail::computeGradientsGreenGauss<2>(solver, kindMpiComm, kindPeriodicComm, geometry, config, field, varBegin,
                                            varEnd, idxVel, gradient);
      break;
    case 3:
      detail::computeGradientsGreenGauss<3>(solver, kindMpiComm, kindPeriodicComm, geometry, config, field, varBegin,
                                            varEnd, idxVel, gradient);
      break;
    default:
      SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);
      break;
  }
}
