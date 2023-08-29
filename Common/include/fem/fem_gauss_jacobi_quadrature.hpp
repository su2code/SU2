/*!
 * \file fem_gauss_jacobi_quadrature.hpp
 * \brief Headers of the functions to compute the integration points of the
          Gauss Jacobi quadrature rules.
          The functions are in the <i>fem_gauss_jacobi_quadrature.cpp</i> file.
          All the functions in this class are based on the program JACOBI_RULE
          of John Burkardt.
 * \author E. van der Weide
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

/*
From: "Burkardt, John" <jburkardt@fsu.edu>
Subject: Re: Gauss-Jacobi Quadrature Code
Date: July 6, 2016 at 4:14:17 AM PDT
To: "Thomas D. Economon" <economon@stanford.edu>
Cc: Edwin van der Weide <E.T.A.vanderWeide@ctw.utwente.nl>

Dear Thomas Economon,

Thank you for your note about jacobi_rule.cpp; I am glad you find it of use.

You are welcome to use this code in any way you choose.

The only reason I put a Gnu LGPL on it was because, from a number of inquiries,
it became clear to me that people were hesitant to use a piece of software if
there was no license information on it at all.  At the time, LGPL seemed the least
restrictive.  I didn't pay attention to the vexing fact that now there are several
versions of LGPL, because I really don't care.

If this statement is not enough for you, let me know what would be better.
I have no objection to you taking a copy of jacobi_rule.cpp and altering the
license to agree with your preferred version of LGPL, for instance.

John Burkardt
________________________________________
From: Thomas D. Economon <economon@stanford.edu>
Sent: Tuesday, July 5, 2016 11:04:34 PM
To: Burkardt, John
Cc: Edwin van der Weide
Subject: Gauss-Jacobi Quadrature Code

Dear John,

I am the lead developer for the open-source SU2 suite, which is primarily for
computational fluid dynamics, but also for solving other PDE systems
(http://su2.stanford.edu, https://github.com/su2code/SU2). Over the past few months,
we have been working on a new high-order discontinuous Galerkin fluid solver within
the codebase, led by Prof. Edwin van der Weide at the University of Twente (cc’d).

We found your attached code to be very helpful in resolving some issues with integration
for our pyramid elements, and we would really like to reuse the implementation if possible.
First, we would like to check if this is ok with you personally, and ask how we can
properly attribute your work in the code? Second, we are curious just what version of the
GNU LGPL license you are using to make sure that we don’t have any licensing issues.
For SU2, we are using GNU LGPL v2.1.

Thanks for the time and take care,
Tom
*/

#pragma once

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

#include "../CConfig.hpp"

using namespace std;

/*!
 * \class CGaussJacobiQuadrature
 * \brief Class used to determine the quadrature points of the Gauss Jacobi
          integration rules.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CGaussJacobiQuadrature {
 public:
  /*!
   * \brief Function, which serves as the API to compute the integration points
            and weights.
   * \param[in]     alpha     Parameter in the weighting function (b-x)^alpha*(x-a)^beta
                              in the Gauss Jacobi rule.
   * \param[in]     beta      Parameter in the weighting function (b-x)^alpha*(x-a)^beta
                              in the Gauss Jacobi rule.
   * \param[in]     a         Lower bound of the integration interval, usually -1.0.
   * \param[in]     b         Upper bound of the integration interval, usually  1.0.
   * \param[in,out] GJPoints  Location of the Gauss-Jacobi integration points.
   * \param[in,out] GJWeights Weights of the Gauss-Jacobi integration points.
   */
  void GetQuadraturePoints(const passivedouble alpha, const passivedouble beta, const passivedouble a,
                           const passivedouble b, vector<passivedouble>& GJPoints, vector<passivedouble>& GJWeights);

 private:
  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void cdgqf(int nt, int kind, passivedouble alpha, passivedouble beta, passivedouble t[], passivedouble wts[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void cgqf(int nt, int kind, passivedouble alpha, passivedouble beta, passivedouble a, passivedouble b,
            passivedouble t[], passivedouble wts[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  passivedouble class_matrix(int kind, int m, passivedouble alpha, passivedouble beta, passivedouble aj[],
                             passivedouble bj[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void imtqlx(int n, passivedouble d[], passivedouble e[], passivedouble z[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void parchk(int kind, int m, passivedouble alpha, passivedouble beta);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  passivedouble r8_epsilon();

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  passivedouble r8_sign(passivedouble x);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void scqf(int nt, const passivedouble t[], const int mlt[], const passivedouble wts[], int nwts, int ndx[],
            passivedouble swts[], passivedouble st[], int kind, passivedouble alpha, passivedouble beta,
            passivedouble a, passivedouble b);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void sgqf(int nt, const passivedouble aj[], passivedouble bj[], passivedouble zemu, passivedouble t[],
            passivedouble wts[]);
};
