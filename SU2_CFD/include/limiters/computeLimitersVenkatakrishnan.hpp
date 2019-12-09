/*!
 * \file computeLimitersVenkatakrishnan.hpp
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
 * \brief Venkatakrishnan specialization.
 */
template<>
struct LimiterDetails<VENKATAKRISHNAN>
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
    su2double lim = ratio*ratio*(1.0+eps2/(delta*delta)) + ratio
    return (lim + ratio) / (lim + 2.0);
  }
};


/*!
 * \brief Create a version of the generalized limiter for Venkatakrishnan.
 */
template<class FieldType, class GradientType>
void computeLimitersVenkatakrishnan(CSolver* solver,
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
  generalizedLimiter<FieldType, GradientType, VENKATAKRISHNAN>(solver,
    kindMpiComm, kindPeriodicComm1, kindPeriodicComm2, geometry,
    varBegin, varEnd, field, gradient, fieldMin, fieldMax, limiter);
}
