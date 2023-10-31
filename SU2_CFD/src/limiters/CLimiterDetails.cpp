/*!
 * \file CLimiterDetails.cpp
 * \brief A class template that allows defining limiters via
 *        specialization of particular details.
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

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../include/limiters/CLimiterDetails.hpp"

/*--- Definition of the static members of the Venkatakrishnan-Wang
 * specialization of CLimiterDetails, need to be here due to ODR. ---*/
su2activevector CLimiterDetails<LIMITER::VENKATAKRISHNAN_WANG>::sharedMin;
su2activevector CLimiterDetails<LIMITER::VENKATAKRISHNAN_WANG>::sharedMax;
