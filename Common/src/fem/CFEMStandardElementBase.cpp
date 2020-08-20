/*!
 * \file CFEMStandardElementBase.cpp
 * \brief Functions for the base class CFEMStandardElementBase.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fem/CFEMStandardElementBase.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardElementBase.                     */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElementBase::GetNDOFsStatic(unsigned short VTK_Type,
                                                       unsigned short nPoly) {
  unsigned short nDOFsEdge = nPoly + 1;
  unsigned short nDOFs;

  switch(VTK_Type) {

    case LINE:
      nDOFs = nDOFsEdge;
      break;

    case TRIANGLE:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case QUADRILATERAL:
      nDOFs = nDOFsEdge*nDOFsEdge;
      break;

    case TETRAHEDRON:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(nDOFsEdge+2)/6;
      break;

    case HEXAHEDRON:
      nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
      break;

    case PRISM:
      nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case PYRAMID:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
      break;

    default:
      nDOFs = 0;  // Indicates an invalid element.
      break;
  }

  return nDOFs;
}

void CFEMStandardElementBase::IntegrationPointsLine(const unsigned short nPoints,
                                                    vector<su2double>    &rLine,
                                                    vector<su2double>    &wLine) {

  /*--- Allocate the memory for rLine and wLine. ---*/
  rLine.resize(nPoints);
  wLine.resize(nPoints);

  /*--- The distribution of points is symmetric. Hence only half the number of
        integration points need to be computed. This means that for an odd number
        the point half way is identically 0. ---*/
  const unsigned short nn = nPoints/2;

  if(2*nn < nPoints) rLine[nn] = 0.0;

  /*--- Constants used in the initial guess of the roots in the loop below. ---*/
  const su2double t1 = 1.0 - (nPoints-1)/(8.0*nPoints*nPoints*nPoints);
  const su2double t2 = PI_NUMBER/(4.0*nPoints + 2.0);

  /*--- The remaing points must be computed. These are the roots of P_n(x),
        P_n is the classic Legendre polynomial of order n. ---*/
  unsigned short ii = nPoints -1;
  for(unsigned i=0; i<nn; ++i, --ii) {

    /*--- Initial guess of this root. ---*/
    su2double x = t1*cos(t2*(4*i+3));

    /*--- Determine the Legendre Polynomials P_n and P_{n-1} and the value f = P_n. ---*/
    su2double Pn   = Legendre(nPoints,   x);
    su2double Pnm1 = Legendre(nPoints-1, x);
    su2double f    = Pn;

    /*--- Solve the root using Halley's method.
          Loop until machine precision has been reached. ---*/
    for(;;)
    {
      /*--- Determine the value of the first and second derivative of f. ---*/
      const su2double df  = nPoints*(Pnm1 - x*Pn)/(1.0-x*x);
      const su2double d2f = (2.0*x*df - nPoints*(nPoints+1)*Pn)/(1.0-x*x);

      /*--- Compute the new value of the root. ---*/
      x = x - 2.0*f*df/(2.0*df*df - f*d2f);

      /*--- Determine the new value of the Legendre polynomials and
            compute the new value of f. Store the old value. ---*/
      const su2double fOld = f;
      Pn   = Legendre(nPoints,   x);
      Pnm1 = Legendre(nPoints-1, x);
      f    = Pn;

      /*--- Convergence criterion. ---*/
      if(fabs(fOld) <= fabs(f)) break;
    }

    /*--- Store the value as well as the symmetric equivalent. ---*/
    rLine[ii] =  x;
    rLine[i]  = -x;
  }

  /*--- Compute the integration weights of the points.
        Make sure the sum is exactly 2. ---*/
  su2double f = 0.0;
  for(unsigned short i=0; i<nPoints; ++i)
  {
    const su2double Pnm1 = Legendre(nPoints-1, rLine[i]);
    wLine[i] = 2.0*(1.0-rLine[i]*rLine[i])/(nPoints*nPoints*Pnm1*Pnm1);
    f       += wLine[i];
  }

  if(fabs(f-2.0) > 1.e-6)
    SU2_MPI::Error(string("Something wrong in computing the weights"),
                   CURRENT_FUNCTION);

  f = 2.0/f;
  for(unsigned short i=0; i<nPoints; ++i) wLine[i] *= f;
}

void CFEMStandardElementBase::Location1DGridDOFsEquidistant(vector<su2double> &r) {

  /*--- Allocate the memory and set the location of the DOFs using
        equidistant spacing. ---*/
  r.resize(nPoly+1);
  const passivedouble dh = 2.0/nPoly;

  for(unsigned short i=0; i<=nPoly; ++i)
    r[i] = -1.0 + i*dh;
}

void CFEMStandardElementBase::Location1DGridDOFsLGL(vector<su2double> &r) {

  /*--- Allocate the memory. ---*/
  const unsigned short nPoints = nPoly+1;
  r.resize(nPoints);

  /*--- The distribution of points is symmetric. Hence only half the number
        of points must be computed. The first and last are at the end of
        the interval and for and add number the mid point is zero. ---*/
  const unsigned short nn = nPoints/2;

  r[0]     = -1.0;
  r[nPoly] =  1.0;

  if(2*nn < nPoints) r[nn] = 0.0;

  /*--- Constants used in the initial guess of the roots in the loop below. ---*/
  const su2double t1 = 1.0 - 3.0*(nPoly-1)/(8.0*nPoly*nPoly*nPoly);
  const su2double t2 = PI_NUMBER/(4.0*nPoly + 1.0);

  /*--- The remaing points must be computed. These are the roots of P'_{n-1}(x),
        P_n is the classic Legendre polynomial of order n. Loop over roots to
        be computed. ---*/
  unsigned short ii = nPoly-1;
  for(unsigned short i=1; i<nn; ++i, --ii) {

    /*--- Initial guess of this root. ---*/
    su2double x = t1*cos(t2*(4*i+1));

    /*--- Determine the Legendre Polynomials P_{n-2} and P_{n-1} and
          the value f = P'_{n-1}(x). ---*/
    su2double Pnm1 = Legendre(nPoly,   x);
    su2double Pnm2 = Legendre(nPoly-1, x);
    su2double f    = nPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

    /*--- Solve the root using Halley's method.
          Loop until machine precision has been reached. ---*/
    for(;;) {

      /*--- Determine the value of the first and second derivative of f. ---*/
      const su2double df  = (2.0*x*f - nPoints*nPoly*Pnm1)/(1.0-x*x);
      const su2double d2f = (2.0*x*df - (nPoints*nPoly-2)*f)/(1.0-x*x);

      /*--- Compute the new value of the root. ---*/
      x = x - 2.0*f*df/(2.0*df*df - f*d2f);

      /*--- Determine the new value of the Legendre polynomials and
            compute the new value of f. Store the old value. ---*/
      const su2double fOld = f;

      Pnm1 = Legendre(nPoly,   x);
      Pnm2 = Legendre(nPoly-1, x);
      f    = nPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

      /*--- Convergence criterion. ---*/
      if(fabs(fOld) <= fabs(f)) break; 
    }

    /*--- Store the value as well as the symmetric equivalent. ---*/
    r[ii] =  x;
    r[i]  = -x;
  }
}

su2double CFEMStandardElementBase::Legendre(unsigned short n,
                                            su2double      x) {

  /*--- Initialization of the polynomials Pnm1 and Pn. ---*/
  su2double Pnm1 = 1.0;
  su2double Pn   = x;

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else {

    /*--- Recursive definition of Pn. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      const su2double tmp = Pnm1;
      Pnm1 = Pn;

      Pn = ((2*i-1)*x*Pn - (i-1)*tmp)/i;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}

su2double CFEMStandardElementBase::NormJacobi(unsigned short n,
                                              unsigned short alpha,
                                              unsigned short beta,
                                              su2double      x) {
  /*--- Some abbreviations. ---*/
  su2double ap1   = alpha + 1;
  su2double bp1   = beta  + 1;
  su2double apb   = alpha + beta;
  su2double apbp1 = apb + 1;
  su2double apbp2 = apb + 2;
  su2double apbp3 = apb + 3;
  su2double b2ma2 = beta*beta - alpha*alpha;

  /*--- Determine the term, which involves the gamma function for P0. As the
        arguments are integers, this term can be computed easily, because
        Gamma(n+1) = n!. ---*/
  su2double Gamap1 = 1.0, Gambp1 = 1.0, Gamapbp2 = 1.0;
  for(unsigned short i=2; i<=alpha; ++i)          Gamap1   *= i;
  for(unsigned short i=2; i<=beta; ++i)           Gambp1   *= i;
  for(unsigned short i=2; i<=(alpha+beta+1); ++i) Gamapbp2 *= i;

  /*--- Initialize the normalized polynomials. ---*/
  su2double Pnm1 = sqrt(pow(0.5,apbp1)*Gamapbp2/(Gamap1*Gambp1));
  su2double Pn   = 0.5*Pnm1*(apbp2*x + alpha - beta)*sqrt(apbp3/(ap1*bp1));

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else
  {
    /*--- The value of the normalized Jacobi polynomial must be obtained via recursion. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      /*--- Compute the coefficients a for i and i-1 and the coefficient bi. ---*/
      unsigned short j = i-1;
      su2double   tmp  = 2*j + apb;
      su2double   aim1 = 2.0*sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))
                       / tmp;

      su2double bi = b2ma2/(tmp*(tmp+2.0));

      tmp          = 2*i + apb;
      su2double ai = 2.0*sqrt(i*(i+apb)*(i+alpha)*(i+beta)/((tmp-1.0)*(tmp+1.0)))
                   / tmp;

      /*--- Compute the new value of Pn and make sure to store Pnm1 correctly. ---*/
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((x-bi)*Pn - aim1*tmp)/ai;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}
