/*!
 * \file computeLimiters_impl.hpp
 * \brief Generic computation of limiters.
 * \note Common methods are derived by defining small details
 *       via specialization of CLimiterDetails.
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


/*!
 * \brief Generic limiter computation for methods based on one limiter
 *        value per point (as opposed to one per edge) and per variable.
 * \ingroup FvmAlgos
 * \note This implementation can be used to derive most common methods
 *       by specializing the limiter functions (e.g. Venkatakrishnan)
 *       and the geometric modifications (e.g. sharp edges), this is done
 *       via specialization of "CLimiterDetails" (all its methods). Then,
 *       a call to this function with the specialized "LimiterKind" should
 *       be added to the body of "computeLimiters()".
 *       See also the notes in computeGradientsGreenGauss.hpp
 *
 * Arguments:
 * \param[in] solver - Optional, solver associated with the field (used only for MPI).
 * \param[in] kindMpiComm - Type of MPI communication required.
 * \param[in] kindPeriodicComm1 - Type of periodic comm. to determine min/max.
 * \param[in] kindPeriodicComm2 - Type of periodic comm. to adjust limiters.
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem.
 * \param[in] varBegin - First variable index for which to compute limiters.
 * \param[in] varEnd - End of computation range (nVar = end-begin).
 * \param[in] field - Variable field.
 * \param[in] gradient - Gradient of the field.
 * \param[out] fieldMin - Minimum field values over direct neighbors of each point.
 * \param[out] fieldMax - As above but maximum values.
 * \param[out] limiter - Reconstruction limiter for the field.
 *
 * Template parameters:
 * \param nDim - Number of dimensions.
 * \param LimiterKind - Used to instantiate the right details class.
 * \param FieldType - Generic object with operator (iPoint,iVar).
 * \param GradientType - Generic object with operator (iPoint,iVar,iDim).
 */
template<size_t nDim, LIMITER LimiterKind, class FieldType, class GradientType>
void computeLimiters_impl(CSolver* solver,
                          MPI_QUANTITIES kindMpiComm,
                          PERIODIC_QUANTITIES kindPeriodicComm1,
                          PERIODIC_QUANTITIES kindPeriodicComm2,
                          CGeometry& geometry,
                          const CConfig& config,
                          size_t varBegin,
                          size_t varEnd,
                          const FieldType& field,
                          const GradientType& gradient,
                          FieldType& fieldMin,
                          FieldType& fieldMax,
                          FieldType& limiter)
{
  constexpr size_t MAXNVAR = 32;

  if (varEnd > MAXNVAR)
    SU2_MPI::Error("Number of variables is too large, increase MAXNVAR.", CURRENT_FUNCTION);

  const size_t nPointDomain = geometry.GetnPointDomain();
  const size_t nPoint = geometry.GetnPoint();

  /*--- If we do not have periodicity we can use a
   *    more efficient access pattern to memory. ---*/

  const bool periodic = (solver != nullptr) &&
                        (kindPeriodicComm1 != PERIODIC_NONE) &&
                        (config.GetnMarker_Periodic() > 0);

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  const auto chunkSize = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

  /*--- If limiters are frozen do not record the computation ---*/
  bool wasActive = false;
  if (config.GetDiscrete_Adjoint() && config.GetFrozen_Limiter_Disc()) {
    wasActive = AD::BeginPassive();
  }

  CLimiterDetails<LimiterKind> limiterDetails;

  limiterDetails.preprocess(geometry, config, varBegin, varEnd, field);

  /*--- Initialize all min/max field values if we have
   *    periodic comms. otherwise do it inside main loop. ---*/

  if (periodic)
  {
    SU2_OMP_FOR_STAT(chunkSize)
    for (size_t iPoint = 0; iPoint < nPoint; ++iPoint)
      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
        fieldMax(iPoint,iVar) = fieldMin(iPoint,iVar) = field(iPoint,iVar);
    END_SU2_OMP_FOR

    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
    {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm1);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm1);
    }
  }

  /*--- Compute limiter for each point. ---*/

  SU2_OMP_FOR_DYN(chunkSize)
  for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
  {
    auto nodes = geometry.nodes;
    const auto coord_i = nodes->GetCoord(iPoint);

    /*--- Cannot preaccumulate if hybrid parallel due to shared reading. ---*/
    if (omp_get_num_threads() == 1) AD::StartPreacc();
    AD::SetPreaccIn(coord_i, nDim);

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
    {
      AD::SetPreaccIn(field(iPoint,iVar));

      if (periodic) {
        /*--- Started outside loop, so counts as input. ---*/
        AD::SetPreaccIn(fieldMax(iPoint,iVar));
        AD::SetPreaccIn(fieldMin(iPoint,iVar));
      }
      else {
        /*--- Initialize min/max now for iPoint if not periodic. ---*/
        fieldMax(iPoint,iVar) = field(iPoint,iVar);
        fieldMin(iPoint,iVar) = field(iPoint,iVar);
      }

      for(size_t iDim = 0; iDim < nDim; ++iDim)
        AD::SetPreaccIn(gradient(iPoint,iVar,iDim));
    }

    /*--- Initialize min/max projection out of iPoint. ---*/

    su2double projMax[MAXNVAR], projMin[MAXNVAR];

    for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      projMax[iVar] = projMin[iVar] = 0.0;

    /*--- Compute max/min projection and values over direct neighbors. ---*/

    for (auto jPoint : geometry.nodes->GetPoints(iPoint)) {

      const auto coord_j = geometry.nodes->GetCoord(jPoint);
      AD::SetPreaccIn(coord_j, nDim);

      /*--- Distance vector from iPoint to face (middle of the edge). ---*/

      su2double dist_ij[nDim] = {0.0};

      for(size_t iDim = 0; iDim < nDim; ++iDim)
        dist_ij[iDim] = 0.5 * (coord_j[iDim] - coord_i[iDim]);

      /*--- Project each variable, update min/max. ---*/

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        su2double proj = 0.0;

        for(size_t iDim = 0; iDim < nDim; ++iDim)
          proj += dist_ij[iDim] * gradient(iPoint,iVar,iDim);

        projMax[iVar] = max(projMax[iVar], proj);
        projMin[iVar] = min(projMin[iVar], proj);

        AD::SetPreaccIn(field(jPoint,iVar));

        fieldMax(iPoint,iVar) = max(fieldMax(iPoint,iVar), field(jPoint,iVar));
        fieldMin(iPoint,iVar) = min(fieldMin(iPoint,iVar), field(jPoint,iVar));
      }
    }

    /*--- Compute the geometric factor. ---*/

    su2double geoFactor = limiterDetails.geometricFactor(iPoint, geometry);

    /*--- Final limiter computation for each variable, get the min limiter
     *    out of the positive/negative projections and deltas. ---*/

    for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
    {
      su2double limMax = limiterDetails.limiterFunction(iVar, projMax[iVar],
                         fieldMax(iPoint,iVar) - field(iPoint,iVar));

      su2double limMin = limiterDetails.limiterFunction(iVar, projMin[iVar],
                         fieldMin(iPoint,iVar) - field(iPoint,iVar));

      limiter(iPoint,iVar) = geoFactor * min(limMax, limMin);

      AD::SetPreaccOut(limiter(iPoint,iVar));
    }

    AD::EndPreacc();
  }
  END_SU2_OMP_FOR

  /*--- Account for periodic effects, take the minimum limiter on each periodic pair. ---*/
  if (periodic)
  {
    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
    {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm2);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm2);
    }
  }

  /*--- Obtain the limiters at halo points from the MPI ranks that own them.
   *    If no solver was provided we do not communicate. ---*/
  if (solver != nullptr)
  {
    solver->InitiateComms(&geometry, &config, kindMpiComm);
    solver->CompleteComms(&geometry, &config, kindMpiComm);
  }

  AD::EndPassive(wasActive);

}
