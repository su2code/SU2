/*!
 * \file math_ad.cpp
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

#include "../include/math_ad.hpp"

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START fabs_d

double fabs_d(double a, double ad, double *fabs) {
    double b;
    double bd;
    double epsilon;
    double arg1;
    double arg1d;
    epsilon = 1e-10;
    arg1d = ad*a + a*ad;
    arg1 = a*a + epsilon*epsilon;
    bd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
    b = sqrt(arg1);
    *fabs = b;
    return bd;
}

//SU2_DIFF END fabs_d

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START pow_d

// Exact pow_d()
double pow_d(double a, double ad, double b, double *c) {
    double cd;

    c[0] = pow(a, b);
    cd = b*ad*pow(a, (b-1.0));
    return cd;
}

//SU2_DIFF END pow_d
