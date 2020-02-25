/*!
 * \file CLimiterDetails.hpp
 * \brief A template class that allows defining limiters via
 *        specialization of particular details.
 * \author P. Gomes
 * \version 7.0.1 "Blackbird"
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


/*!
 * \brief A traits class for limiters, see notes for "computeLimiters_impl()".
 * \note There is no default implementation (the code will compile but not
 *       link) specialization is mandatory.
 */
template<ENUM_LIMITER LimiterKind>
struct CLimiterDetails
{
  /*!
   * \brief Compute any global value that may be needed by the other functions.
   * \note This function is called once by multiple threads.
   */
  template<class FieldType>
  inline void preprocess(CGeometry&, CConfig&, size_t varBegin,
                         size_t varEnd, const FieldType&);

  /*!
   * \brief Geometric modifier (e.g. increase limiting near sharp edges).
   * \note This function is called once per point inside an AD pre-
   *       -accumulation region, newly used variables should be registered.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry&) const;

  /*!
   * \brief Smooth (usually) function of the maximum/minimum (positive/negative)
   *        gradient projections onto the edges, and the deltas over direct neighbors.
   *        Both proj and delta may be 0.0, beware of divisions.
   * \note This function is called twice (min/max) per point per variable
   *       (also inside an AD pre-accumulation region).
   */
  inline su2double limiterFunction(size_t iVar, su2double proj, su2double delta) const;
};


/*!
 * \brief Common small functions used by limiters.
 */
namespace LimiterHelpers
{
  inline passivedouble epsilon() {return std::numeric_limits<passivedouble>::epsilon();}

  inline su2double venkatFunction(su2double proj, su2double delta, su2double eps2)
  {
    su2double y = delta*(delta+proj) + eps2;
    return (y + delta*proj) / (y + 2*proj*proj);
  }

  inline su2double raisedSine(su2double dist)
  {
    su2double factor = 0.5*(1.0+dist+sin(PI_NUMBER*dist)/PI_NUMBER);
    return max(0.0, min(factor, 1.0));
  }
}


/*!
 * \brief Barth-Jespersen specialization.
 */
template<>
struct CLimiterDetails<BARTH_JESPERSEN>
{
  su2double eps2;

  /*!
   * \brief Set a small epsilon to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(Ts&...) {eps2 = LimiterHelpers::epsilon();}

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Venkatakrishnan function with a numerical epsilon.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers::venkatFunction(proj, delta, eps2);
  }
};


/*!
 * \brief Venkatakrishnan specialization.
 */
template<>
struct CLimiterDetails<VENKATAKRISHNAN>
{
  su2double eps2;

  /*!
   * \brief Store the reference lenght based eps^2 parameter,
   *        limited to a small number to avoid divisions by 0.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, CConfig& config, Ts&...)
  {
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    su2double eps1 = fabs(L*K);
    eps2 = max(eps1*eps1*eps1, LimiterHelpers::epsilon());
  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers::venkatFunction(proj, delta, eps2);
  }
};


/*!
 * \brief Venkatakrishnan-Wang specialization.
 */
template<>
struct CLimiterDetails<VENKATAKRISHNAN_WANG>
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
      fieldMin.resize(nThread, nCols) = largeNum;
      fieldMax.resize(nThread, nCols) =-largeNum;
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
        fieldMin(iThread, iVar) = min(fieldMin(iThread, iVar), field(iPoint, iVar));
        fieldMax(iThread, iVar) = max(fieldMax(iThread, iVar), field(iPoint, iVar));
      }
    }

    SU2_OMP_MASTER
    {
      /*--- Per rank reduction. ---*/

      for(size_t iThread = 1; iThread < nThread; ++iThread)
      {
        for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        {
          fieldMin(0,iVar) = min(fieldMin(0,iVar), fieldMin(iThread, iVar));
          fieldMax(0,iVar) = max(fieldMax(0,iVar), fieldMax(iThread, iVar));
        }
      }

      /*--- Global reduction, (re)using eps2 as the recv buffer. ---*/

      SU2_MPI::Allreduce(fieldMin[0], eps2.data(), nCols, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        fieldMin(0,iVar) = eps2(iVar);

      SU2_MPI::Allreduce(fieldMax[0], eps2.data(), nCols, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
        fieldMax(0,iVar) = eps2(iVar);

      /*--- Compute eps^2 ---*/

      su2double K = config.GetVenkat_LimiterCoeff();

      for(size_t iVar = varBegin; iVar < varEnd; ++iVar)
      {
        su2double range = fieldMax(0,iVar) - fieldMin(0,iVar);
        eps2(iVar) = max(pow(K*range, 2), LimiterHelpers::epsilon());
      }
    }
    SU2_OMP_BARRIER

  }

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts&...) const {return 1.0;}

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t iVar, su2double proj, su2double delta) const
  {
    AD::SetPreaccIn(eps2(iVar));
    return LimiterHelpers::venkatFunction(proj, delta, eps2(iVar));
  }
};


/*!
 * \brief Venkatakrishnan with sharp edge modification.
 */
template<>
struct CLimiterDetails<SHARP_EDGES>
{
  su2double eps1, eps2, sharpCoeff;

  /*!
   * \brief Store the reference lenght based eps^2 parameter.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, CConfig& config, Ts&...)
  {
    sharpCoeff = config.GetAdjSharp_LimiterCoeff();
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps1 = fabs(L*K);
    eps2 = max(eps1*eps1*eps1, LimiterHelpers::epsilon());
  }

  /*!
   * \brief Full limiting (1st order) near sharp edges.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry& geometry) const
  {
    AD::SetPreaccIn(geometry.node[iPoint]->GetSharpEdge_Distance());
    su2double dist = geometry.node[iPoint]->GetSharpEdge_Distance()/(sharpCoeff*eps1)-1.0;
    return LimiterHelpers::raisedSine(dist);
  }

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers::venkatFunction(proj, delta, eps2);
  }
};


/*!
 * \brief Venkatakrishnan with wall distance modification.
 */
template<>
struct CLimiterDetails<WALL_DISTANCE>
{
  su2double eps1, eps2, sharpCoeff;

  /*!
   * \brief Store the reference lenght based eps^2 parameter.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, CConfig& config, Ts&...)
  {
    sharpCoeff = config.GetAdjSharp_LimiterCoeff();
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps1 = fabs(L*K);
    eps2 = max(eps1*eps1*eps1, LimiterHelpers::epsilon());
  }

  /*!
   * \brief Full limiting (1st order) near walls.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry& geometry) const
  {
    AD::SetPreaccIn(geometry.node[iPoint]->GetWall_Distance());
    su2double dist = geometry.node[iPoint]->GetWall_Distance()/(sharpCoeff*eps1)-1.0;
    return LimiterHelpers::raisedSine(dist);
  }

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double limiterFunction(size_t, su2double proj, su2double delta) const
  {
    return LimiterHelpers::venkatFunction(proj, delta, eps2);
  }
};
