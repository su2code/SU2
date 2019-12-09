/*!
 * \file CLimiterDetails.hpp
 * \brief A template class that allows defining limiters via
 *        specialization of particular details.
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


/*!
 * \brief A traits class for limiters, see notes for "computeLimiters_impl()".
 * \note There is no default implementation the code will compile but not
 *       link, specialization is mandatory.
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
   * \brief Smooth (usually) function of the minimum ratio (delta/projection)
   *        and of the maximum absolute delta, over direct neighbors of a point.
   * \note This function is called once per point and per variable
   *       (also inside an AD pre-accumulation region).
   */
  inline su2double smoothFunction(size_t iVar, su2double ratio, su2double delta) const;
};


/*!
 * \brief Barth-Jespersen specialization.
 */
template<>
struct CLimiterDetails<BARTH_JESPERSEN>
{
  /*!
   * \brief Nothing to preprocess.
   */
  template<class... Ts>
  inline void preprocess(Ts...) {}

  /*!
   * \brief No geometric modification for this kind of limiter.
   */
  template<class... Ts>
  inline su2double geometricFactor(Ts...) const {return 1.0;}

  /*!
   * \brief Smooth function of min(2,ratio).
   */
  inline su2double smoothFunction(size_t, su2double ratio, su2double) const
  {
    ratio = std::min(su2double(2.0), ratio);
    su2double lim = ratio*ratio + ratio;
    return (lim + ratio) / (lim + 2.0);
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
   * \brief Store the reference lenght based eps^2 parameter.
   */
  template<class... Ts>
  inline void preprocess(CGeometry&, CConfig& config, Ts...)
  {
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps2 = std::pow(L*K, 3);
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
    su2double lim = ratio*ratio*(1.0+eps2/(delta*delta)) + ratio;
    return (lim + ratio) / (lim + 2.0);
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
  inline su2double smoothFunction(size_t iVar, su2double ratio, su2double delta) const
  {
    AD::SetPreaccIn(eps2(iVar));

    su2double lim = ratio*ratio*(1.0+eps2(iVar)/(delta*delta)) + ratio;
    return (lim + ratio) / (lim + 2.0);
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
  inline void preprocess(CGeometry&, CConfig& config, Ts...)
  {
    sharpCoeff = config.GetAdjSharp_LimiterCoeff();
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps1 = L*K;
    eps2 = eps1*eps1*eps1;
  }

  /*!
   * \brief Full limiting (1st order) near sharp edges.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry& geometry) const
  {
    AD::SetPreaccIn(geometry.node[iPoint]->GetSharpEdge_Distance());

    su2double dist = (geometry.node[iPoint]->GetSharpEdge_Distance() - sharpCoeff*eps1);

    if(dist < -eps1) return 0.0;
    else if (dist > eps1) return 1.0;
    else return 0.5*(1.0+(dist/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*dist/eps1));
  }

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double smoothFunction(size_t, su2double ratio, su2double delta) const
  {
    su2double lim = ratio*ratio*(1.0+eps2/(delta*delta)) + ratio;
    return (lim + ratio) / (lim + 2.0);
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
  inline void preprocess(CGeometry&, CConfig& config, Ts...)
  {
    sharpCoeff = config.GetAdjSharp_LimiterCoeff();
    su2double L = config.GetRefElemLength();
    su2double K = config.GetVenkat_LimiterCoeff();
    eps1 = L*K;
    eps2 = eps1*eps1*eps1;
  }

  /*!
   * \brief Full limiting (1st order) near sharp edges.
   */
  inline su2double geometricFactor(size_t iPoint, CGeometry& geometry) const
  {
    AD::SetPreaccIn(geometry.node[iPoint]->GetWall_Distance());

    su2double dist = (geometry.node[iPoint]->GetWall_Distance() - sharpCoeff*eps1);

    if(dist < -eps1) return 0.0;
    else if (dist > eps1) return 1.0;
    else return 0.5*(1.0+(dist/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*dist/eps1));
  }

  /*!
   * \brief Smooth function that disables limiting in smooth regions.
   */
  inline su2double smoothFunction(size_t, su2double ratio, su2double delta) const
  {
    su2double lim = ratio*ratio*(1.0+eps2/(delta*delta)) + ratio;
    return (lim + ratio) / (lim + 2.0);
  }
};
