/*!
 * \file computeGradientsGreenGaussLimited.hpp
 * \brief Green-Gauss gradients limited on the fly, in a single pass over the points.
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

#include "computeGradientsGreenGauss.hpp"
#include "../limiters/computeLimiters.hpp"

namespace detail {

/*!
 * \brief Compute Green-Gauss gradients and scale them by a point-based limiter, in one pass.
 * \ingroup FvmAlgos
 * \note Compared to computing the gradients and then the limiters, this avoids a second pass
 *       over the mesh (the neighbors of each point are visited twice, but they are in cache
 *       the second time), and it avoids storing the limiters altogether.
 * \note The gradient uses the (f_j-f_i)/2 face value, which is equivalent to (f_i+f_j)/2 for
 *       a closed control volume (the face normals sum to zero) but makes the boundary
 *       contributions vanish, and so the loop over boundary markers is not needed.
 * \note With U-MUSCL (nonzero MUSCL_KAPPA) the result is not identical to limiting during the
 *       reconstruction, because the reconstruction is affine and not linear in the gradient,
 *       its centered part ends up being limited as well.
 * \note Periodic points are handled, but only approximately, their limiter is based on the
 *       extrema of the neighbors on this side of the interface (the ones on the other side
 *       would require the limiter communications that this method exists to avoid). That
 *       limiter is too weak to be relied on, which is why CConfig::LimitedGradientPossible
 *       does not allow this method on periodic meshes.
 * \param[in] solver - Optional, solver associated with the field (used only for MPI).
 * \param[in] kindMpiComm - Type of MPI communication required.
 * \param[in] kindPeriodicComm - Type of periodic communication required.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem.
 * \param[in] umusclKappa - Blending parameter for U-MUSCL reconstruction.
 * \param[in] field - Generic object implementing operator (iPoint, iVar).
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[in] idxVel - Index of velocity, or -1 if no velocity present.
 * \param[out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 */
template <size_t nDim, LIMITER LimiterKind, class FieldType, class GradientType>
void computeGradientsGreenGaussLimited(CSolver* solver, MPI_QUANTITIES kindMpiComm,
                                       PERIODIC_QUANTITIES kindPeriodicComm, CGeometry& geometry,
                                       const CConfig& config, su2double umusclKappa, const FieldType& field,
                                       const size_t varBegin, const size_t varEnd, const int idxVel,
                                       GradientType& gradient) {
  constexpr size_t MAXNVAR = 32;

  if (varEnd > MAXNVAR) SU2_MPI::Error("Number of variables is too large, increase MAXNVAR.", CURRENT_FUNCTION);

  const size_t nPointDomain = geometry.GetnPointDomain();

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  CLimiterDetails<LimiterKind> limiterDetails;

  limiterDetails.preprocess(geometry, config, varBegin, varEnd, field);

  /*--- Points on periodic boundaries only own part of their control volume, the rest is held
   *    by their periodic partner, so their gradient is only complete after the periodic
   *    communications. They are limited in a second (small) pass, limiting the partial
   *    gradient would be meaningless. Points not on a periodic boundary have no partner
   *    volume, which makes them cheap to tell apart. ---*/

  const bool periodic = (solver != nullptr) && (config.GetnMarker_Periodic() > 0);

  auto isPeriodic = [&](size_t iPoint) {
    return periodic && (geometry.nodes->GetPeriodicVolume(iPoint) != 0.0);
  };

  /*--- Scale the gradient of a point by its limiter, which is not stored. The extrema of the
   *    field over the neighbors are taken as input because they are gathered by the gradient
   *    loop, only the projections require a second sweep over the neighbors, since they need
   *    the gradient to be complete. Pass a tracker for "extrema" to gather them here instead. ---*/

  auto limitPoint = [&](size_t iPoint, const su2double* minVals, const su2double* maxVals, auto extrema) {
    su2double projMin[MAXNVAR], projMax[MAXNVAR];

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar) projMin[iVar] = projMax[iVar] = 0.0;

    limiterExtrema<nDim>(geometry, iPoint, varBegin, varEnd, umusclKappa, field, gradient, extrema, projMin, projMax);

    const su2double geoFactor = limiterDetails.geometricFactor(iPoint, geometry);

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
      const su2double limMax = limiterDetails.limiterFunction(iVar, projMax[iVar], maxVals[iVar] - field(iPoint, iVar));
      const su2double limMin = limiterDetails.limiterFunction(iVar, projMin[iVar], minVals[iVar] - field(iPoint, iVar));
      const su2double limiter = geoFactor * min(limMax, limMin);

      for (size_t iDim = 0; iDim < nDim; ++iDim) {
        gradient(iPoint, iVar, iDim) *= limiter;
        AD::SetPreaccOut(gradient(iPoint, iVar, iDim));
      }
    }
  };

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();
    AD::SetPreaccIn(geometry.nodes->GetCoord(iPoint), nDim);

    su2double minVals[MAXNVAR], maxVals[MAXNVAR];

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar) minVals[iVar] = maxVals[iVar] = field(iPoint, iVar);

    /*--- Integrate over the faces of the control volume, the difference formulation makes the
     *    boundary faces drop out of the summation. The min/max of the field over the neighbors
     *    are determined by the same loop, their values are loaded for the gradient anyway. ---*/

    greenGaussPointContributions<nDim, true>(geometry, iPoint, varBegin, varEnd, field, gradient,
                                             CPointMinMax{minVals, maxVals});

    if (isPeriodic(iPoint)) {
      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        for (size_t iDim = 0; iDim < nDim; ++iDim) AD::SetPreaccOut(gradient(iPoint, iVar, iDim));
    } else {
      limitPoint(iPoint, minVals, maxVals, CNoMinMax{});
    }

    AD::EndPreacc();
  }
  END_SU2_OMP_FOR

  /*--- Account for periodic contributions, this completes the gradient of periodic points,
   *    which can then be limited. Their extrema are still those of the neighbors on this
   *    side of the boundary, no limiter communications are needed. ---*/

  if (periodic) {
    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic() / 2; ++iPeriodic) {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm);
    }

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      if (!isPeriodic(iPoint)) continue;

      if (omp_get_num_threads() == 1) AD::StartPreacc();
      AD::SetPreaccIn(geometry.nodes->GetCoord(iPoint), nDim);

      su2double minVals[MAXNVAR], maxVals[MAXNVAR];

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar) {
        AD::SetPreaccIn(field(iPoint, iVar));
        minVals[iVar] = maxVals[iVar] = field(iPoint, iVar);
        for (size_t iDim = 0; iDim < nDim; ++iDim) AD::SetPreaccIn(gradient(iPoint, iVar, iDim));
      }

      /*--- The extrema of these points were gathered before the communications, but the
       *    gradient they were meant for was incomplete, so they are gathered again here. ---*/

      limitPoint(iPoint, minVals, maxVals, CPointMinMax{minVals, maxVals});

      AD::EndPreacc();
    }
    END_SU2_OMP_FOR
  }

  /*--- Compute the corrections for symmetry planes and Euler walls. Unlike in the two-pass
   *    procedure these are applied to the limited gradients, but the property they enforce
   *    (no normal gradient of scalars, mirrored velocity gradients) is not affected by the
   *    scaling, only by the correction being applied last. ---*/

  correctGradientsSymmetry<nDim>(geometry, config, varBegin, varEnd, idxVel, gradient);

  /*--- If no solver was provided we do not communicate ---*/

  if (solver == nullptr) return;

  /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

  solver->InitiateComms(&geometry, &config, kindMpiComm);
  solver->CompleteComms(&geometry, &config, kindMpiComm);
}
}  // namespace detail

/*!
 * \brief Instantiations for 2D and 3D, and for each kind of limiter.
 * \ingroup FvmAlgos
 * \note LIMITER::NONE and edge-based limiters are not supported, callers are expected to use
 *       the plain Green-Gauss gradients in those cases (see CConfig::GetLimitedGradientRecon).
 */
template <class FieldType, class GradientType>
void computeGradientsGreenGaussLimited(LIMITER LimiterKind, CSolver* solver, MPI_QUANTITIES kindMpiComm,
                                       PERIODIC_QUANTITIES kindPeriodicComm, CGeometry& geometry,
                                       const CConfig& config, su2double umusclKappa, const FieldType& field,
                                       const size_t varBegin, const size_t varEnd, const int idxVel,
                                       GradientType& gradient) {
  if (geometry.GetnDim() != 2 && geometry.GetnDim() != 3)
    SU2_MPI::Error("Too many dimensions to compute gradients.", CURRENT_FUNCTION);

#define INSTANTIATE(KIND)                                                                                       \
  if (geometry.GetnDim() == 2) {                                                                                \
    detail::computeGradientsGreenGaussLimited<2, KIND>(solver, kindMpiComm, kindPeriodicComm, geometry, config,  \
                                                       umusclKappa, field, varBegin, varEnd, idxVel, gradient);  \
  } else {                                                                                                      \
    detail::computeGradientsGreenGaussLimited<3, KIND>(solver, kindMpiComm, kindPeriodicComm, geometry, config,  \
                                                       umusclKappa, field, varBegin, varEnd, idxVel, gradient);  \
  }
  switch (LimiterKind) {
    case LIMITER::BARTH_JESPERSEN: INSTANTIATE(LIMITER::BARTH_JESPERSEN); break;
    case LIMITER::VENKATAKRISHNAN: INSTANTIATE(LIMITER::VENKATAKRISHNAN); break;
    case LIMITER::NISHIKAWA_R3: INSTANTIATE(LIMITER::NISHIKAWA_R3); break;
    case LIMITER::NISHIKAWA_R4: INSTANTIATE(LIMITER::NISHIKAWA_R4); break;
    case LIMITER::NISHIKAWA_R5: INSTANTIATE(LIMITER::NISHIKAWA_R5); break;
    case LIMITER::VENKATAKRISHNAN_WANG: INSTANTIATE(LIMITER::VENKATAKRISHNAN_WANG); break;
    case LIMITER::WALL_DISTANCE: INSTANTIATE(LIMITER::WALL_DISTANCE); break;
    case LIMITER::SHARP_EDGES: INSTANTIATE(LIMITER::SHARP_EDGES); break;
    default: SU2_MPI::Error("Unknown limiter type.", CURRENT_FUNCTION); break;
  }
#undef INSTANTIATE
}
