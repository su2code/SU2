/*!
 * \file generalizedLimiter.hpp
 * \brief Generic and general computation of limiters.
 * \note Common methods are derived by defining only the details.
 * \author P. Gomes
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/omp_structure.hpp"


/*!
 * \brief A traits class for limiters, see notes for "generalizedLimiter()".
 */
template<ENUM_LIMITER LimiterKind>
struct LimiterDetails
{
  static_assert(false, "The default LimiterDetails should not be used.");

  /*!
   * \brief Compute any global value that may be needed by the other functions.
   * \note This function is called once by multiple threads.
   */
  template<class FieldType>
  inline void preprocess(CGeometry&, CConfig&, size_t varBegin,
                         size_t varEnd, const FieldType&) {}

  /*!
   * \brief Geometric modifier (e.g. increase limiting near sharp edges).
   * \note This function is called once per point inside an AD pre-
   *       -accumulation region, newly used variables should be registered.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry&) const {return 1.0;}

  /*!
   * \brief Smooth (usually) function of the minimum ratio (delta/projection)
   *        and of the maximum absolute delta, over direct neighbors of a point.
   * \note This function is called once per point and per variable
   *       (also inside an AD pre-accumulation region).
   */
  inline su2double smoothFunction(size_t iVar, su2double ratio,
                                  su2double delta) const {return 1.0;}
};


/*!
 * \brief General limiter computation for methods based on one limiter
 *        value per point (as opposed to one per edge) and per variable.
 * \note This implementation can be used to derive most common methods
 *       by specializing the limiter functions (e.g. Venkatakrishnan)
 *       and the geometric modifications (e.g. sharp edges), this is done
 *       via specialization of "LimiterDetails" (all its methods). Then,
 *       a function that wraps the call to this one with the specialized
 *       "LimiterKind" should be created. If that function is used by more
 *       than one solver do not inline the seemingly trivial definition,
 *       put it in its own .cpp instead so it gets compiled only once.
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
 * \param FieldType - Generic object with operator (iPoint,iVar)
 * \param GradientType - Generic object with operator (iPoint,iVar,iDim)
 * \param LimiterKind - Used to instantiate the right details class.
 */
template<class FieldType, class GradientType, ENUM_LIMITER LimiterKind>
void generalizedLimiter(CSolver* solver,
                        MPI_QUANTITIES kindMpiComm,
                        PERIODIC_QUANTITIES kindPeriodicComm1,
                        PERIODIC_QUANTITIES kindPeriodicComm2,
                        CGeometry& geometry,
                        CConfig& config,
                        size_t varBegin,
                        size_t varEnd,
                        const FieldType& field,
                        const GradientType& gradient,
                        FieldType& fieldMin,
                        FieldType& fieldMax,
                        FieldType& limiter)
{
  constexpr size_t MAXNDIM = 3;
  constexpr size_t MAXNVAR = 8;

  if (varEnd > MAXNVAR)
    SU2_MPI::Error("Number of variables is too large, increase MAXNVAR.", CURRENT_FUNCTION);

  LimiterDetails<LimiterKind> details;

  su2double eps = std::numeric_limits<passivedouble>::epsilon();

  size_t nPointDomain = geometry.GetnPointDomain();
  size_t nPoint = geometry.GetnPoint();
  size_t nDim = geometry.GetnDim();

  /*--- If we do not have periodicity we can use a
   *    more efficient access pattern to memory. ---*/

  bool periodic = (solver != nullptr) &&
                  (kindPeriodicComm1 != PERIODIC_NONE) &&
                  (config.GetnMarker_Periodic() > 0);

#ifdef HAVE_OMP
  constexpr size_t OMP_MAX_CHUNK = 512;

  size_t chunkSize = computeStaticChunkSize(nPointDomain,
                     omp_get_max_threads(), OMP_MAX_CHUNK);
#endif

#ifdef CODI_REVERSE_TYPE
  bool tapeActive = false;

  if (config.GetDiscrete_Adjoint() && config.GetFrozen_Limiter_Disc()) {
    /*--- If limiters are frozen do not record the computation ---*/
    tapeActive = AD::globalTape.isActive();
    AD::StopRecording();
  }
#endif

  /*--- Start OpenMP parallel section. ---*/

  SU2_OMP_PARALLEL
  {
    details.preprocess(geometry, config, varBegin, varEnd, field);
    SU2_OMP_BARRIER

    /*--- Initialize all min/max field values if we have
     *    periodic comms. otherwise do it inside main loop. ---*/

    if (periodic)
    {
      SU2_OMP_FOR_STAT(chunkSize)
      for (size_t iPoint = 0; iPoint < nPoint; ++iPoint)
        for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
          fieldMax(iPoint,iVar) = fieldMin(iPoint,iVar) = field(iPoint,iVar);

      SU2_OMP_MASTER
      {
        for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
        {
          solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm1);
          solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm1);
        }
      }
      SU2_OMP_BARRIER
    }

    /*--- Compute limiter for each point. ---*/

    SU2_OMP_FOR_DYN(chunkSize)
    for (size_t iPoint = 0; iPoint < nPointDomain; ++iPoint)
    {
      auto node = geometry.node[iPoint];
      const su2double* coord_i = node->GetCoord();

      AD::StartPreacc();
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

      /*--- Initialize min/max projection out of iPoint, this is done
       *    with a small value to prevent eventual division by 0. ---*/

      su2double projMax[MAXNVAR], projMin[MAXNVAR];

      for (size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        projMax[iVar] = eps;
        projMin[iVar] =-eps;
      }

      /*--- Compute max/min projection and values over direct neighbors. ---*/

      for(size_t iNeigh = 0; iNeigh < node->GetnPoint(); ++iNeigh)
      {
        size_t jPoint = node->GetPoint(iNeigh);

        const su2double* coord_j = geometry.node[jPoint]->GetCoord();
        AD::SetPreaccIn(coord_j, nDim);

        /*--- Distance vector from iPoint to face (middle of the edge). ---*/

        su2double dist_ij[MAXNDIM];

        for(size_t iDim = 0; iDim < nDim; ++iDim)
          dist_ij[iDim] = 0.5 * (coord_j[iDim] - coord_i[iDim]);

        /*--- Project each variable, update min/max. ---*/

        for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        {
          su2double proj = 0.0;

          for(size_t iDim = 0; iDim < nDim; ++iDim)
            proj += dist_ij[iDim] * gradient(iPoint,iVar,iDim);

          projMax[iVar] = std::max(projMax[iVar], proj);
          projMin[iVar] = std::min(projMin[iVar], proj);

          AD::SetPreaccIn(field(jPoint,iVar));

          fieldMax(iPoint,iVar) = std::max(fieldMax(iPoint,iVar), field(jPoint,iVar));
          fieldMin(iPoint,iVar) = std::min(fieldMin(iPoint,iVar), field(jPoint,iVar));
        }
      }

      /*--- Compute the geometric factor. ---*/

      su2double geoFactor = details.geometricFactor(iPoint, geometry);

      /*--- Final limiter computation for each variable. ---*/

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        su2double deltaMax = fieldMax(iPoint,iVar) - field(iPoint,iVar);
        su2double deltaMin = fieldMin(iPoint,iVar) - field(iPoint,iVar);
        su2double delta = std::max(std::max(deltaMax, -deltaMin), eps);

        su2double ratioMax = deltaMax / projMax[iVar];
        su2double ratioMin = deltaMin / projMin[iVar];
        su2double ratio = std::min(limMax, limMin);

        limiter(iPoint,iVar) = geoFactor * details.smoothFunction(iVar, ratio, delta);

        AD::SetPreaccOut(limiter(iPoint,iVar));
      }

      AD::EndPreacc();
    }

  } // end SU2_OMP_PARALLEL

  /*--- If no solver was provided we do not communicate. ---*/

  if (solver != nullptr)
  {
    /*--- Account for periodic effects, take the minimum limiter on each periodic pair.
     *    Recall that only the field min/max were communicated and not the projections, the
     *    resulting limiter is only correct if the smooth limiter function is monotonic. ---*/

    for (size_t iPeriodic = 1; iPeriodic <= config.GetnMarker_Periodic()/2; ++iPeriodic)
    {
      solver->InitiatePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm2);
      solver->CompletePeriodicComms(&geometry, &config, iPeriodic, kindPeriodicComm2);
    }

    /*--- Obtain the gradients at halo points from the MPI ranks that own them. ---*/

    solver->InitiateComms(&geometry, &config, kindMpiComm);
    solver->CompleteComms(&geometry, &config, kindMpiComm);
  }

#ifdef CODI_REVERSE_TYPE
  if (tapeActive) AD::StartRecording();
#endif
}
