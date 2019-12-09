/*!
 * \file computeLimitersWallDistance.hpp
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
 * \brief Venkatakrishnan with wall distance modification.
 */
template<>
struct LimiterDetails<WALL_DISTANCE>
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
    su2double lim = ratio*ratio*(1.0+eps2/(delta*delta)) + ratio
    return (lim + ratio) / (lim + 2.0);
  }
};


/*!
 * \brief Create a version of the generalized limiter for wall distance.
 */
template<class FieldType, class GradientType>
void computeLimitersWallDistance(CSolver* solver,
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
  generalizedLimiter<FieldType, GradientType, WALL_DISTANCE>(solver,
    kindMpiComm, kindPeriodicComm1, kindPeriodicComm2, geometry,
    varBegin, varEnd, field, gradient, fieldMin, fieldMax, limiter);
}
