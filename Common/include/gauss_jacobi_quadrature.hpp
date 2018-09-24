/*!
 * \file gauss_jacobi_quadrature.hpp
 * \brief Headers of the functions to compute the integration points of the
          Gauss Jacobi quadrature rules.
          The functions are in the <i>gauss_jacobi_quadrature.cpp</i> file.
          All the functions in this class are based on the program JACOBI_RULE
          of John Burkardt.
 * \author E. van der Weide
 * \version 6.1.0 "Falcon"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#pragma once

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

#include "config_structure.hpp"

using namespace std;

/*!
 * \class CGaussJacobiQuadrature
 * \brief Class used to determine the quadrature points of the Gauss Jacobi
          integration rules.
 * \author E. van der Weide
 * \version 6.1.0 "Falcon"
 */
class CGaussJacobiQuadrature {
public:

  /*!
   * \brief Constructor of the class, nothing to be done.
   */
  CGaussJacobiQuadrature();

  /*!
   * \brief Destructor of the class, nothing to be done.
   */
  ~CGaussJacobiQuadrature();

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
  void GetQuadraturePoints(const passivedouble   alpha,     const passivedouble   beta,
                           const passivedouble   a,         const passivedouble   b,
                           vector<passivedouble> &GJPoints, vector<passivedouble> &GJWeights);
private:
  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void cdgqf(int nt, int kind, passivedouble alpha, passivedouble beta, passivedouble t[],
             passivedouble wts[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void cgqf(int nt, int kind, passivedouble alpha, passivedouble beta, passivedouble a,
            passivedouble b, passivedouble t[], passivedouble wts[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  passivedouble class_matrix(int kind, int m, passivedouble alpha, passivedouble beta,
                             passivedouble aj[], passivedouble bj[]);

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
  void scqf(int nt, passivedouble t[], int mlt[], passivedouble wts[], int nwts, int ndx[],
            passivedouble swts[], passivedouble st[], int kind, passivedouble alpha,
            passivedouble beta, passivedouble a, passivedouble b);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void sgqf(int nt, passivedouble aj[], passivedouble bj[], passivedouble zemu, passivedouble t[],
            passivedouble wts[]);
};

#include "gauss_jacobi_quadrature.inl"
