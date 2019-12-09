/*!
 * \file computeLimitersVenkatWang.hpp
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

#include "generalizedLimiter.hpp"


/*!
 * \brief Venkatakrishnan-Wang specialization.
 */
template<>
struct LimiterDetails<VENKATAKRISHNAN_WANG>
{
  su2activematrix fieldMin, fieldMax;
  su2activevector eps2;

  /*!
   * \brief Store the solution range based eps^2 parameter.
   */
  template<class FieldType>
  inline void preprocess(CGeometry& geometry, CConfig& config, size_t varBegin,
                         size_t varEnd, const FieldType& field)
  {
    /*--- Determine the max and min global value for each variable.
     *    Each thread initially works on one row of fieldMin/Max (the
     *    rows are padded to a multiple of 8 to avoid false sharing),
     *    then the master thread performs a reduction over all threads
     *    and mpi ranks onto row 0, the final result. ---*/

    size_t nThread = omp_get_num_threads();
    size_t nCols = roundUpDiv(varEnd,8)*8;

    SU2_OMP_MASTER
    {
      su2double largeNum = 0.1*std::numeric_limits<passivedouble>::max();
      fieldMin.resize(nThreads, nCols) = largeNum;
      fieldMax.resize(nThreads, nCols) =-largeNum;
      eps2.resize(nCols);
    }
    SU2_OMP_BARRIER

    /*--- Per thread reduction. ---*/

    SU2_OMP_FOR_STAT(512)
    for(size_t iPoint = 0; iPoint < geometry.GetnPointDomain(); ++iPoint)
    {
      size_t iThread = omp_get_thread_num();

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        fieldMin(iThread, iVar) = std::min(fieldMin(iThread, iVar), field(iPoint, iVar));
        fieldMax(iThread, iVar) = std::max(fieldMax(iThread, iVar), field(iPoint, iVar));
      }
    }

    SU2_OMP_MASTER
    {
      /*--- Per rank reduction. ---*/

      for(size_t iThread = 1; iThread < nThread; ++iThread)
      {
        for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        {
          fieldMin(0,iVar) = std::min(fieldMin(0,iVar), fieldMin(iThread, iVar));
          fieldMax(0,iVar) = std::max(fieldMax(0,iVar), fieldMax(iThread, iVar));
        }
      }

      /*--- Global reduction. ---*/

      SU2_MPI::Allreduce(fieldMin[0], eps2.data(), nCols, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        fieldMin(0,iVar) = eps2(iVar);

      SU2_MPI::Allreduce(fieldMax[0], eps2.data(), nCols, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        fieldMax(0,iVar) = eps2(iVar);

      /*--- Compute eps^2 ---*/

      su2double K = config.GetVenkat_LimiterCoeff();

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        eps2(iVar) = std::pow(K*(fieldMax(0,iVar) - fieldMin(0,iVar)), 2);
    }
    SU2_OMP_BARRIER

  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double smoothFunction(size_t, su2double ratio, su2double delta) const
  {
    AD::SetPreaccIn(eps2[iVar]);

    su2double lim = ratio*ratio*(1.0+eps2[iVar]/(delta*delta)) + ratio
    return (lim + ratio) / (lim + 2.0);
  }
};


/*!
 * \brief Create a version of the generalized limiter for Venkatakrishnan.
 */
template<class FieldType, class GradientType>
void computeLimitersVenkatWang(CSolver* solver,
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
  generalizedLimiter<FieldType, GradientType, VENKATAKRISHNAN_WANG>(
    solver, kindMpiComm, kindPeriodicComm1, kindPeriodicComm2, geometry,
    varBegin, varEnd, field, gradient, fieldMin, fieldMax, limiter);
}
