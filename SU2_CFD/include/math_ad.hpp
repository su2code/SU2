/*!
 * \file math_ad.hpp
 * \brief Differentiated math routines
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START fabs_d

double fabs_d(double a, double ad, double *fabs) ;

//SU2_DIFF END fabs_d

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START pow_d

double pow_d(double a, double ad, double b, double *c) ;

//SU2_DIFF END pow_d
