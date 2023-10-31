/*!
 * \file fem_gauss_jacobi_quadrature.cpp
 * \brief Functions to compute the points and weights for the Gauss-Jacobi
          quadrature rules. All the functions in this file are based on the
          program JACOBI_RULE of John Burkardt.
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

#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

void CGaussJacobiQuadrature::GetQuadraturePoints(const passivedouble alpha, const passivedouble beta,
                                                 const passivedouble a, const passivedouble b,
                                                 vector<passivedouble>& GJPoints, vector<passivedouble>& GJWeights) {
  /*--- Determine the number of integration points. Check if the number makes sense. ---*/
  auto nIntPoints = (unsigned int)GJPoints.size();
  if (nIntPoints < 1 || nIntPoints > 100)
    SU2_MPI::Error("Invalid number of Gauss Jacobi integration points", CURRENT_FUNCTION);

  /*--- Call the function cgqf to do the actual work. ---*/
  int kind = 4;
  cgqf(nIntPoints, kind, alpha, beta, a, b, GJPoints.data(), GJWeights.data());
}

//****************************************************************************80

void CGaussJacobiQuadrature::cdgqf(int nt, int kind, passivedouble alpha, passivedouble beta, passivedouble t[],
                                   passivedouble wts[])

//****************************************************************************80
//
//  Purpose:
//
//    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with a classical weight function with default values for A and B,
//    and only simple knots.
//
//    There are no moments checks and no printing is done.
//
//    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, passivedouble ALPHA, the value of Alpha, if needed.
//
//    Input, passivedouble BETA, the value of Beta, if needed.
//
//    Output, passivedouble T[NT], the knots.
//
//    Output, passivedouble WTS[NT], the weights.
//
{
  passivedouble* aj;
  passivedouble* bj;
  passivedouble zemu;

  parchk(kind, 2 * nt, alpha, beta);
  //
  //  Get the Jacobi matrix and zero-th moment.
  //
  aj = new passivedouble[nt];
  bj = new passivedouble[nt];

  zemu = class_matrix(kind, nt, alpha, beta, aj, bj);
  //
  //  Compute the knots and weights.
  //
  sgqf(nt, aj, bj, zemu, t, wts);

  delete[] aj;
  delete[] bj;
}
//****************************************************************************80

void CGaussJacobiQuadrature::cgqf(int nt, int kind, passivedouble alpha, passivedouble beta, passivedouble a,
                                  passivedouble b, passivedouble t[], passivedouble wts[])

//****************************************************************************80
//
//  Purpose:
//
//    CGQF computes knots and weights of a Gauss quadrature formula.
//
//  Discussion:
//
//    The user may specify the interval (A,B).
//
//    Only simple knots are produced.
//
//    Use routine EIQFS to evaluate this quadrature formula.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, passivedouble ALPHA, the value of Alpha, if needed.
//
//    Input, passivedouble BETA, the value of Beta, if needed.
//
//    Input, passivedouble A, B, the interval endpoints, or
//    other parameters.
//
//    Output, passivedouble T[NT], the knots.
//
//    Output, passivedouble WTS[NT], the weights.
//
{
  int i;
  int* mlt;
  int* ndx;
  //
  //  Compute the Gauss quadrature formula for default values of A and B.
  //
  cdgqf(nt, kind, alpha, beta, t, wts);
  //
  //  Prepare to scale the quadrature formula to other weight function with
  //  valid A and B.
  //
  mlt = new int[nt];
  for (i = 0; i < nt; i++) {
    mlt[i] = 1;
  }
  ndx = new int[nt];
  for (i = 0; i < nt; i++) {
    ndx[i] = i + 1;
  }
  scqf(nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);

  delete[] mlt;
  delete[] ndx;
}
//****************************************************************************80

passivedouble CGaussJacobiQuadrature::class_matrix(int kind, int m, passivedouble alpha, passivedouble beta,
                                                   passivedouble aj[], passivedouble bj[])

//****************************************************************************80
//
//  Purpose:
//
//    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
//
//  Discussion:
//
//    This routine computes the diagonal AJ and sub-diagonal BJ
//    elements of the order M tridiagonal symmetric Jacobi matrix
//    associated with the polynomials orthogonal with respect to
//    the weight function specified by KIND.
//
//    For weight functions 1-7, M elements are defined in BJ even
//    though only M-1 are needed.  For weight function 8, BJ(M) is
//    set to zero.
//
//    The zero-th moment of the weight function is returned in ZEMU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the Jacobi matrix.
//
//    Input, passivedouble ALPHA, the value of Alpha, if needed.
//
//    Input, passivedouble BETA, the value of Beta, if needed.
//
//    Output, passivedouble AJ[M], BJ[M], the diagonal and subdiagonal
//    of the Jacobi matrix.
//
//    Output, passivedouble CLASS_MATRIX, the zero-th moment.
//
{
  passivedouble a2b2;
  passivedouble ab;
  passivedouble aba;
  passivedouble abi;
  passivedouble abj;
  passivedouble abti;
  passivedouble apone;
  int i;
  passivedouble pi = 3.14159265358979323846264338327950;
  passivedouble temp;
  passivedouble temp2;
  passivedouble zemu;

  temp = r8_epsilon();

  parchk(kind, 2 * m - 1, alpha, beta);

  temp2 = 0.5;

  if (500.0 * temp < fabs(pow(tgamma(temp2), 2) - pi)) {
    cout << "\n";
    cout << "CLASS_MATRIX - Fatal error!\n";
    cout << "  Gamma function does not match machine parameters.\n";
    exit(1);
  }

  if (kind == 1) {
    ab = 0.0;

    zemu = 2.0 / (ab + 1.0);

    for (i = 0; i < m; i++) {
      aj[i] = 0.0;
    }

    for (i = 1; i <= m; i++) {
      abi = i + ab * (i % 2);
      abj = 2 * i + ab;
      bj[i - 1] = sqrt(abi * abi / (abj * abj - 1.0));
    }
  } else if (kind == 2) {
    zemu = pi;

    for (i = 0; i < m; i++) {
      aj[i] = 0.0;
    }

    bj[0] = sqrt(0.5);
    for (i = 1; i < m; i++) {
      bj[i] = 0.5;
    }
  } else if (kind == 3) {
    ab = alpha * 2.0;
    zemu = pow(2.0, ab + 1.0) * pow(tgamma(alpha + 1.0), 2) / tgamma(ab + 2.0);

    for (i = 0; i < m; i++) {
      aj[i] = 0.0;
    }

    bj[0] = sqrt(1.0 / (2.0 * alpha + 3.0));
    for (i = 2; i <= m; i++) {
      bj[i - 1] = sqrt(i * (i + ab) / (4.0 * pow(i + alpha, 2) - 1.0));
    }
  } else if (kind == 4) {
    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = pow(2.0, ab + 1.0) * tgamma(alpha + 1.0) * tgamma(beta + 1.0) / tgamma(abi);
    aj[0] = (beta - alpha) / abi;
    bj[0] = sqrt(4.0 * (1.0 + alpha) * (1.0 + beta) / ((abi + 1.0) * abi * abi));
    a2b2 = beta * beta - alpha * alpha;

    for (i = 2; i <= m; i++) {
      abi = 2.0 * i + ab;
      aj[i - 1] = a2b2 / ((abi - 2.0) * abi);
      abi = abi * abi;
      bj[i - 1] = sqrt(4.0 * i * (i + alpha) * (i + beta) * (i + ab) / ((abi - 1.0) * abi));
    }
  } else if (kind == 5) {
    zemu = tgamma(alpha + 1.0);

    for (i = 1; i <= m; i++) {
      aj[i - 1] = 2.0 * i - 1.0 + alpha;
      bj[i - 1] = sqrt(i * (i + alpha));
    }
  } else if (kind == 6) {
    zemu = tgamma((alpha + 1.0) / 2.0);

    for (i = 0; i < m; i++) {
      aj[i] = 0.0;
    }

    for (i = 1; i <= m; i++) {
      bj[i - 1] = sqrt((i + alpha * (i % 2)) / 2.0);
    }
  } else if (kind == 7) {
    ab = alpha;
    zemu = 2.0 / (ab + 1.0);

    for (i = 0; i < m; i++) {
      aj[i] = 0.0;
    }

    for (i = 1; i <= m; i++) {
      abi = i + ab * (i % 2);
      abj = 2 * i + ab;
      bj[i - 1] = sqrt(abi * abi / (abj * abj - 1.0));
    }
  } else  // if ( kind == 8 )
  {
    ab = alpha + beta;
    zemu = tgamma(alpha + 1.0) * tgamma(-(ab + 1.0)) / tgamma(-beta);
    apone = alpha + 1.0;
    aba = ab * apone;
    aj[0] = -apone / (ab + 2.0);
    bj[0] = -aj[0] * (beta + 1.0) / (ab + 2.0) / (ab + 3.0);
    for (i = 2; i <= m; i++) {
      abti = ab + 2.0 * i;
      aj[i - 1] = aba + 2.0 * (ab + i) * (i - 1);
      aj[i - 1] = -aj[i - 1] / abti / (abti - 2.0);
    }

    for (i = 2; i <= m - 1; i++) {
      abti = ab + 2.0 * i;
      bj[i - 1] = i * (alpha + i) / (abti - 1.0) * (beta + i) / (abti * abti) * (ab + i) / (abti + 1.0);
    }
    bj[m - 1] = 0.0;
    for (i = 0; i < m; i++) {
      bj[i] = sqrt(bj[i]);
    }
  }

  return zemu;
}
//****************************************************************************80

void CGaussJacobiQuadrature::imtqlx(int n, passivedouble d[], passivedouble e[], passivedouble z[])

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine.
//
//    It has been modified to produce the product Q' * Z, where Z is an input
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.
//    The changes consist (essentialy) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, passivedouble D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, passivedouble E(N), the subdiagonal entries of the
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, passivedouble Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  passivedouble b;
  passivedouble c;
  passivedouble f;
  passivedouble g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m = 0;  // To avoid a compiler warning.
  int mml;
  passivedouble p;
  passivedouble prec;
  passivedouble r;
  passivedouble s;

  prec = r8_epsilon();

  if (n == 1) {
    return;
  }

  e[n - 1] = 0.0;

  for (l = 1; l <= n; l++) {
    j = 0;
    for (;;) {
      for (m = l; m <= n; m++) {
        if (m == n) {
          break;
        }

        if (fabs(e[m - 1]) <= prec * (fabs(d[m - 1]) + fabs(d[m]))) {
          break;
        }
      }
      p = d[l - 1];
      if (m == l) {
        break;
      }
      if (itn <= j) {
        cout << "\n";
        cout << "IMTQLX - Fatal error!\n";
        cout << "  Iteration limit exceeded\n";
        exit(1);
      }
      j = j + 1;
      g = (d[l] - p) / (2.0 * e[l - 1]);
      r = sqrt(g * g + 1.0);
      g = d[m - 1] - p + e[l - 1] / (g + fabs(r) * r8_sign(g));
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for (ii = 1; ii <= mml; ii++) {
        i = m - ii;
        f = s * e[i - 1];
        b = c * e[i - 1];

        if (fabs(g) <= fabs(f)) {
          c = g / f;
          r = sqrt(c * c + 1.0);
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        } else {
          s = f / g;
          r = sqrt(s * s + 1.0);
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = (d[i - 1] - g) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i - 1] + c * f;
        z[i - 1] = c * z[i - 1] - s * f;
      }
      d[l - 1] = d[l - 1] - p;
      e[l - 1] = g;
      e[m - 1] = 0.0;
    }
  }
  //
  //  Sorting.
  //
  for (ii = 2; ii <= m; ii++) {
    i = ii - 1;
    k = i;
    p = d[i - 1];

    for (j = ii; j <= n; j++) {
      if (d[j - 1] < p) {
        k = j;
        p = d[j - 1];
      }
    }

    if (k != i) {
      d[k - 1] = d[i - 1];
      d[i - 1] = p;
      p = z[i - 1];
      z[i - 1] = z[k - 1];
      z[k - 1] = p;
    }
  }
}
//****************************************************************************80

void CGaussJacobiQuadrature::parchk(int kind, int m, passivedouble alpha, passivedouble beta)

//****************************************************************************80
//
//  Purpose:
//
//    PARCHK checks parameters ALPHA and BETA for classical weight functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the highest moment to
//    be calculated.  This value is only needed when KIND = 8.
//
//    Input, passivedouble ALPHA, BETA, the parameters, if required
//    by the value of KIND.
//
{
  passivedouble tmp;

  if (kind <= 0) {
    cout << "\n";
    cout << "PARCHK - Fatal error!\n";
    cout << "  KIND <= 0.\n";
    exit(1);
  }
  //
  //  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
  //
  if (3 <= kind && alpha <= -1.0) {
    cout << "\n";
    cout << "PARCHK - Fatal error!\n";
    cout << "  3 <= KIND and ALPHA <= -1.\n";
    exit(1);
  }
  //
  //  Check BETA for Jacobi.
  //
  if (kind == 4 && beta <= -1.0) {
    cout << "\n";
    cout << "PARCHK - Fatal error!\n";
    cout << "  KIND == 4 and BETA <= -1.0.\n";
    exit(1);
  }
  //
  //  Check ALPHA and BETA for rational.
  //
  if (kind == 8) {
    tmp = alpha + beta + m + 1.0;
    if (0.0 <= tmp || tmp <= beta) {
      cout << "\n";
      cout << "PARCHK - Fatal error!\n";
      cout << "  KIND == 8 but condition on ALPHA and BETA fails.\n";
      exit(1);
    }
  }
}
//****************************************************************************80

passivedouble CGaussJacobiQuadrature::r8_epsilon()

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, passivedouble R8_EPSILON, the R8 round-off unit.
//
{
  const passivedouble value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

passivedouble CGaussJacobiQuadrature::r8_sign(passivedouble x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, passivedouble X, the number whose sign is desired.
//
//    Output, passivedouble R8_SIGN, the sign of X.
//
{
  passivedouble value;

  if (x < 0.0) {
    value = -1.0;
  } else {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

void CGaussJacobiQuadrature::scqf(int nt, const passivedouble t[], const int mlt[], const passivedouble wts[], int nwts,
                                  int ndx[], passivedouble swts[], passivedouble st[], int kind, passivedouble alpha,
                                  passivedouble beta, passivedouble a, passivedouble b)

//****************************************************************************80
//
//  Purpose:
//
//    SCQF scales a quadrature formula to a nonstandard interval.
//
//  Discussion:
//
//    The arrays WTS and SWTS may coincide.
//
//    The arrays T and ST may coincide.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, passivedouble T[NT], the original knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, passivedouble WTS[NWTS], the weights.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.
//    For more details see the comments in CAWIQ.
//
//    Output, passivedouble SWTS[NWTS], the scaled weights.
//
//    Output, passivedouble ST[NT], the scaled knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, passivedouble ALPHA, the value of Alpha, if needed.
//
//    Input, passivedouble BETA, the value of Beta, if needed.
//
//    Input, passivedouble A, B, the interval endpoints.
//
{
  passivedouble al;
  passivedouble be;
  int i;
  int k;
  int l;
  passivedouble p;
  passivedouble shft;
  passivedouble slp;
  passivedouble temp;
  passivedouble tmp;

  temp = r8_epsilon();

  parchk(kind, 1, alpha, beta);

  if (kind == 1) {
    al = 0.0;
    be = 0.0;
    if (fabs(b - a) <= temp) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit(1);
    }
    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;
  } else if (kind == 2) {
    al = -0.5;
    be = -0.5;
    if (fabs(b - a) <= temp) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit(1);
    }
    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;
  } else if (kind == 3) {
    al = alpha;
    be = alpha;
    if (fabs(b - a) <= temp) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit(1);
    }
    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;
  } else if (kind == 4) {
    al = alpha;
    be = beta;

    if (fabs(b - a) <= temp) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit(1);
    }
    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;
  } else if (kind == 5) {
    if (b <= 0.0) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  B <= 0\n";
      exit(1);
    }
    shft = a;
    slp = 1.0 / b;
    al = alpha;
    be = 0.0;
  } else if (kind == 6) {
    if (b <= 0.0) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  B <= 0.\n";
      exit(1);
    }
    shft = a;
    slp = 1.0 / sqrt(b);
    al = alpha;
    be = 0.0;
  } else if (kind == 7) {
    al = alpha;
    be = 0.0;
    if (fabs(b - a) <= temp) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit(1);
    }
    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;
  } else if (kind == 8) {
    if (a + b <= 0.0) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  A + B <= 0.\n";
      exit(1);
    }
    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;
  } else  // if ( kind == 9 )
  {
    al = 0.5;
    be = 0.5;
    if (fabs(b - a) <= temp) {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit(1);
    }
    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;
  }

  p = pow(slp, al + be + 1.0);

  for (k = 0; k < nt; k++) {
    st[k] = shft + slp * t[k];
    l = abs(ndx[k]);

    if (l != 0) {
      tmp = p;
      for (i = l - 1; i <= l - 1 + mlt[k] - 1; i++) {
        swts[i] = wts[i] * tmp;
        tmp = tmp * slp;
      }
    }
  }
}
//****************************************************************************80

void CGaussJacobiQuadrature::sgqf(int nt, const passivedouble aj[], passivedouble bj[], passivedouble zemu,
                                  passivedouble t[], passivedouble wts[])

//****************************************************************************80
//
//  Purpose:
//
//    SGQF computes knots and weights of a Gauss Quadrature formula.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with simple knots from the Jacobi matrix and the zero-th
//    moment of the weight function, using the Golub-Welsch technique.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, passivedouble AJ[NT], the diagonal of the Jacobi matrix.
//
//    Input/output, passivedouble BJ[NT], the subdiagonal of the Jacobi
//    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
//
//    Input, passivedouble ZEMU, the zero-th moment of the weight function.
//
//    Output, passivedouble T[NT], the knots.
//
//    Output, passivedouble WTS[NT], the weights.
//
{
  int i;
  //
  //  Exit if the zero-th moment is not positive.
  //
  if (zemu <= 0.0) {
    cout << "\n";
    cout << "SGQF - Fatal error!\n";
    cout << "  ZEMU <= 0.\n";
    exit(1);
  }
  //
  //  Set up vectors for IMTQLX.
  //
  for (i = 0; i < nt; i++) {
    t[i] = aj[i];
  }
  wts[0] = sqrt(zemu);
  for (i = 1; i < nt; i++) {
    wts[i] = 0.0;
  }
  //
  //  Diagonalize the Jacobi matrix.
  //
  imtqlx(nt, t, bj, wts);

  for (i = 0; i < nt; i++) {
    wts[i] = wts[i] * wts[i];
  }
}
