/*!
 * \file fem_standard_element.hpp
 * \brief Headers of the main functions for the FEM standard elements.
 *        The functions are in the <i>fem_standard_element.cpp</i> file.
 * \author E. van der Weide
 * \version 4.1.2 "Cardinal"
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
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include <iostream>
#include <vector>
#include <cstdlib>

#include "config_structure.hpp"

using namespace std;

/*!
 * \class FEMStandardElementBaseClass
 * \brief Base class for a FEM standard element.
 * \author E. van der Weide
 * \version 4.1.2 "Cardinal"
 */
class FEMStandardElementBaseClass {
protected:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short orderExact;   /*!< \brief Polynomial order that must be integrated exactly by the integration rule. */
  unsigned short nIntegration; /*!< \brief Number of points used in the numerical integration. */

  bool constJacobian;          /*!< \brief Whether or not the element has a constant Jacobian. */

  vector<su2double> rIntegration; /*!< \brief r-location of the integration points for this standard element. */
  vector<su2double> sIntegration; /*!< \brief s-location of the integration points for this standard element, if needed. */
  vector<su2double> tIntegration; /*!< \brief t-location of the integration points for this standard element, if needed. */
  vector<su2double> wIntegration; /*!< \brief The weights of the integration points for this standard element. */

public:
  /*!
  * \brief Standard Constructor. Nothing to be done.
  */
  FEMStandardElementBaseClass();

  /*!
  * \brief Destructor. Nothing to be done, because the vectors are deleted automatically.
  */
  virtual ~FEMStandardElementBaseClass();

protected:
  /*!
  * \brief Alternative constructor.
  * \param[in] val_VTK_Type   - Type of the element using the VTK convention.
  * \param[in] val_nPoly      - Polynomial degree of the element.
  * \param[in] val_constJac   - Whether or not the Jacobians are constant.
  * \param[in] config         - Object, which contains the input parameters.
  * \param[in] val_orderExact - The order of the polynomials that must be integrated
                                exactly by the integration rule. If 0, it will
                                be determined from the polynomial degree and the
                                parameters in config.
  */
  FEMStandardElementBaseClass(unsigned short val_VTK_Type,
                              unsigned short val_nPoly,
                              bool           val_constJac,
                              CConfig        *config,
                              unsigned short val_orderExact);

public:
  /*!
  * \brief Static function, which makes available the number of DOFs for an element
           corresponding to the arguments.
  * \param[in] VTK_Type         - Type of the element using the VTK convention.
  * \param[in] nPoly            - Polynomial degree of the element.
  * \param[in] typeErrorMessage - Default argument used to write a good error message.
  * \return The number of DOFs
  */
  static unsigned short GetNDOFsStatic(unsigned short VTK_Type,
                                       unsigned short nPoly,
                                       unsigned long  typeErrorMessage = 0);

  /*!
  * \brief Function, which makes available the number of integration points for this standard element.
  * \return  The number of integration points of this standard element.
  */
  unsigned short GetNIntegration(void);

  /*!
  * \brief Static function, which makes available the number of integration points for an element
           corresponding to the arguments.
  * \param[in] VTK_Type   - Type of the element using the VTK convention.
  * \param[in] orderExact - Polynomial degree that must be integrated exactly.
  * \param[in] config     - Object, which contains the input parameters.
  * \return  The number of integration points.
  */
  static unsigned short GetNIntegrationStatic(unsigned short VTK_Type,
                                              unsigned short nPoly,
                                              CConfig        *config);
  /*!
  * \brief Function, which makes available the polynomial order that must be integrated exactly.
  * \return  The polynomial order that must be integrated exactly.
  */
  unsigned short GetOrderExact(void);

protected:
  /*!
  * \brief Function, which checks if the sum of the given derivatives of the
           Lagrangian interpolation functions is 0 in the points.
  * \param[in] nPoints         - Number of points to be checked.
  * \param[in] nDOFs           - Number of DOFs of the element.
  * \param[in] dLagBasisPoints - Values of the derivatives of the Lagrangian
                                 interpolation functions in the given points.
  */
  void CheckSumDerivativesLagrangianBasisFunctions(const unsigned short    nPoints,
                                                   const unsigned short    nDOFs,
                                                   const vector<su2double> &dLagBasisPoints);

  /*!
  * \brief Function, which checks if the sum of the given Lagrangian interpolation
           functions is 1 in the points.
  * \param[in]     nPoints        - Number of points to be checked.
  * \param[in]     nDOFs          - Number of DOFs of the element.
  * \param[in,out] lagBasisPoints - Values of the Lagrangian interpolation
                                    functions in the given points.
  */
  void CheckSumLagrangianBasisFunctions(const unsigned short nPoints,
                                        const unsigned short nDOFs,
                                        vector<su2double>    &lagBasisPoints);

  /*!
  * \brief Function, which copies the data of the given object into the current object.
  * \param[in] other - Object, whose data is copied.
  */
  void Copy(const FEMStandardElementBaseClass &other);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a line.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the line element.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the line element.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesLine(const unsigned short    nPoly,
                                                 const vector<su2double> &rPoints,
                                                 unsigned short          &nDOFs,
                                                 vector<su2double>       &rDOFs,
                                                 vector<su2double>       &lagBasisPoints,
                                                 vector<su2double>       &drLagBasisPoints);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a triangle.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  sPoints          - s-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the triangle.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the triangle.
  * \param[out] sDOFs            - s-coordinates of the DOFs of the triangle.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesTriangle(const unsigned short    nPoly,
                                                     const vector<su2double> &rPoints,
                                                     const vector<su2double> &sPoints,
                                                     unsigned short          &nDOFs,
                                                     vector<su2double>       &rDOFs,
                                                     vector<su2double>       &sDOFs,
                                                     vector<su2double>       &lagBasisPoints,
                                                     vector<su2double>       &drLagBasisPoints,
                                                     vector<su2double>       &dsLagBasisPoints);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a quadrilateral.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  sPoints          - s-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the quadrilateral.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the quadrilateral.
  * \param[out] sDOFs            - s-coordinates of the DOFs of the quadrilateral.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesQuadrilateral(const unsigned short    nPoly,
                                                          const vector<su2double> &rPoints,
                                                          const vector<su2double> &sPoints,
                                                          unsigned short          &nDOFs,
                                                          vector<su2double>       &rDOFs,
                                                          vector<su2double>       &sDOFs,
                                                          vector<su2double>       &lagBasisPoints,
                                                          vector<su2double>       &drLagBasisPoints,
                                                          vector<su2double>       &dsLagBasisPoints);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a tetrahedron.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  sPoints          - s-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  tPoints          - t-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the quadrilateral.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the tetrahedron.
  * \param[out] sDOFs            - s-coordinates of the DOFs of the tetrahedron.
  * \param[out] tDOFs            - t-coordinates of the DOFs of the tetrahedron.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesTetrahedron(const unsigned short    nPoly,
                                                        const vector<su2double> &rPoints,
                                                        const vector<su2double> &sPoints,
                                                        const vector<su2double> &tPoints,
                                                        unsigned short          &nDOFs,
                                                        vector<su2double>       &rDOFs,
                                                        vector<su2double>       &sDOFs,
                                                        vector<su2double>       &tDOFs,
                                                        vector<su2double>       &lagBasisPoints,
                                                        vector<su2double>       &drLagBasisPoints,
                                                        vector<su2double>       &dsLagBasisPoints,
                                                        vector<su2double>       &dtLagBasisPoints);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a pyramid.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  sPoints          - s-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  tPoints          - t-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the quadrilateral.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the pyramid.
  * \param[out] sDOFs            - s-coordinates of the DOFs of the pyramid.
  * \param[out] tDOFs            - t-coordinates of the DOFs of the pyramid.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesPyramid(const unsigned short    nPoly,
                                                    const vector<su2double> &rPoints,
                                                    const vector<su2double> &sPoints,
                                                    const vector<su2double> &tPoints,
                                                    unsigned short          &nDOFs,
                                                    vector<su2double>       &rDOFs,
                                                    vector<su2double>       &sDOFs,
                                                    vector<su2double>       &tDOFs,
                                                    vector<su2double>       &lagBasisPoints,
                                                    vector<su2double>       &drLagBasisPoints,
                                                    vector<su2double>       &dsLagBasisPoints,
                                                    vector<su2double>       &dtLagBasisPoints);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a prism.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  sPoints          - s-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  tPoints          - t-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the quadrilateral.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the prism.
  * \param[out] sDOFs            - s-coordinates of the DOFs of the prism.
  * \param[out] tDOFs            - t-coordinates of the DOFs of the prism.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesPrism(const unsigned short    nPoly,
                                                  const vector<su2double> &rPoints,
                                                  const vector<su2double> &sPoints,
                                                  const vector<su2double> &tPoints,
                                                  unsigned short          &nDOFs,
                                                  vector<su2double>       &rDOFs,
                                                  vector<su2double>       &sDOFs,
                                                  vector<su2double>       &tDOFs,
                                                  vector<su2double>       &lagBasisPoints,
                                                  vector<su2double>       &drLagBasisPoints,
                                                  vector<su2double>       &dsLagBasisPoints,
                                                  vector<su2double>       &dtLagBasisPoints);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a hexahedron.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  sPoints          - s-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[in]  tPoints          - t-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the quadrilateral.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the hexahedron.
  * \param[out] sDOFs            - s-coordinates of the DOFs of the hexahedron.
  * \param[out] tDOFs            - t-coordinates of the DOFs of the hexahedron.
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesHexahedron(const unsigned short    nPoly,
                                                       const vector<su2double> &rPoints,
                                                       const vector<su2double> &sPoints,
                                                       const vector<su2double> &tPoints,
                                                       unsigned short          &nDOFs,
                                                       vector<su2double>       &rDOFs,
                                                       vector<su2double>       &sDOFs,
                                                       vector<su2double>       &tDOFs,
                                                       vector<su2double>       &lagBasisPoints,
                                                       vector<su2double>       &drLagBasisPoints,
                                                       vector<su2double>       &dsLagBasisPoints,
                                                       vector<su2double>       &dtLagBasisPoints);
private:
  /*!
  * \brief Function, which determines the 1D Gauss Legendre integration points and weights.
  * \param[in,out] GLPoints  - The location of the Gauss-Legendre integration points.
  * \param[in,out] GLWeights - The weights of the Gauss-Legendre integration points.
  */
  void GaussLegendrePoints1D(vector<su2double> &GLPoints,
                             vector<su2double> &GLWeights);

  /*!
  * \brief Function, which computes the value of the gradient of the Jacobi polynomial for the given x-coordinate.
  * \param[in] n     - Order of the Jacobi polynomial.
  * \param[in] alpha - Alpha coefficient of the Jacobi polynomial.
  * \param[in] beta  - Beta coefficient of the Jacobi polynomial.
  * \param[in] x     - Coordinate (-1 <= x <= 1) for which the gradient of the Jacobi polynomial must be evaluated.
  * \return            The value of the gradient of the normalized Jacobi polynomial f order n for the given value of x.
  */
  su2double GradNormJacobi(unsigned short n,
                           unsigned short alpha,
                           unsigned short beta,
                           su2double      x);

  /*!
  * \brief Function, which computes the gradient of the Vandermonde matrix for a standard 1D edge.
  * \param[in]  nDOFs - Number of DOFs, which in 1D is the polynomial degree + 1.
  * \param[in]  r     - Parametric coordinates for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient of the Vandermonde matrix in all r-locations.
  */
  void GradVandermonde1D(unsigned short          nDOFs,
                         const vector<su2double> &r,
                         vector<su2double>       &VDr);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard triangle.
  * \param[in]  nPoly - Polynomial degree of the triangle.
  * \param[in]  nDOFs - Number of DOFs of the triangle.
  * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r- and s-locations.
  * \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r- and s-locations.
  */
  void GradVandermonde2D_Triangle(unsigned short          nPoly,
                                  unsigned short          nDOFs,
                                  const vector<su2double> &r,
                                  const vector<su2double> &s,
                                  vector<su2double>       &VDr,
                                  vector<su2double>       &VDs);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard quadrilateral.
  * \param[in]  nPoly - Polynomial degree of the quadrilateral.
  * \param[in]  nDOFs - Number of DOFs of the quadrilateral.
  * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r- and s-locations.
  * \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r- and s-locations.
  */
  void GradVandermonde2D_Quadrilateral(unsigned short          nPoly,
                                       unsigned short          nDOFs,
                                       const vector<su2double> &r,
                                       const vector<su2double> &s,
                                       vector<su2double>       &VDr,
                                       vector<su2double>       &VDs);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard tetrahedron.
  * \param[in]  nPoly - Polynomial degree of the tetrahedron.
  * \param[in]  nDOFs - Number of DOFs of the tetrahedron.
  * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt   - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Tetrahedron(unsigned short          nPoly,
                                     unsigned short          nDOFs,
                                     const vector<su2double> &r,
                                     const vector<su2double> &s,
                                     const vector<su2double> &t,
                                     vector<su2double>       &VDr,
                                     vector<su2double>       &VDs,
                                     vector<su2double>       &VDt);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard pyramid.
  * \param[in]  nPoly - Polynomial degree of the pyramid.
  * \param[in]  nDOFs - Number of DOFs of the pyramid.
  * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt   - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Pyramid(unsigned short          nPoly,
                                 unsigned short          nDOFs,
                                 const vector<su2double> &r,
                                 const vector<su2double> &s,
                                 const vector<su2double> &t,
                                 vector<su2double>       &VDr,
                                 vector<su2double>       &VDs,
                                 vector<su2double>       &VDt);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard prism.
  * \param[in]  nPoly - Polynomial degree of the prism.
  * \param[in]  nDOFs - Number of DOFs of the prism.
  * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt   - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Prism(unsigned short          nPoly,
                               unsigned short          nDOFs,
                               const vector<su2double> &r,
                               const vector<su2double> &s,
                               const vector<su2double> &t,
                               vector<su2double>       &VDr,
                               vector<su2double>       &VDs,
                               vector<su2double>       &VDt);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard hexahedron.
  * \param[in]  nPoly - Polynomial degree of the hexahedron.
  * \param[in]  nDOFs - Number of DOFs of the hexahedron.
  * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt   - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Hexahedron(unsigned short          nPoly,
                                    unsigned short          nDOFs,
                                    const vector<su2double> &r,
                                    const vector<su2double> &s,
                                    const vector<su2double> &t,
                                    vector<su2double>       &VDr,
                                    vector<su2double>       &VDs,
                                    vector<su2double>       &VDt);

  /*!
  * \brief Function, which determines the integration points for a line
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsLine(void);

  /*!
  * \brief Function, which determines the integration points for a triangle
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsTriangle(void);

  /*!
  * \brief Function, which determines the integration points for a quadrilateral
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsQuadrilateral(void);

  /*!
  * \brief Function, which determines the integration points for a tetrahedron
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsTetrahedron(void);

  /*!
  * \brief Function, which determines the integration points for a pyramid
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsPyramid(void);

  /*!
  * \brief Function, which determines the integration points for a prism
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsPrism(void);

  /*!
  * \brief Function, which determines the integration points for a hexahedron
           such that polynomials of orderExact are integrated exactly.
  */
  void IntegrationPointsHexahedron(void);

  /*!
  * \brief Function, which computes the inverse of the given square matrix.
  * \param[in]     n - Number of rows/columns of the square matrix A.
  * \param[in,out] A - On input the square matrix to be inverted. On output the inverse.
  */
  void InverseMatrix(unsigned short    n,
                     vector<su2double> &A);

  /*!
  * \brief Function, which computes the value of the Legendre polynomials
           Pn and Pnm1 for given x (-1 <= x <= 1) and n.
  * \param[in]  x    - x-coordinate for which the Legendre polynomial must be computed.
  * \param[in]  n    - order of the Legendre polynomial.
  * \param[out] Pnm1 - Legendre polynomials of order n-1, which must be computed.
  * \param[out] Pn   - Legendre polynomials of order n, which must be computed.
  */
  void Legendre(su2double      x,
                unsigned short n,
                su2double      &Pnm1,
                su2double      &Pn);

  /*!
  * \brief Function, which carries out a matrix matrix multiplication to obtain
           data in points and stores the transpose of the result.
  * \param[in]  nDOFs   - Dimension of the matrices. This typically corresponds
                          to the number of DOFs that are considered.
  * \param[in]  nPoints - Dimension of the matrices. This typically corresponds
                          to the number of integration points used.
  * \param[in]  A       - First matrix in the matrix matrix product A*B, dimension nPoints X nDOFs.
  * \param[in]  B       - Second matrix in the matrix matrix product A*B, dimension nDOFs X nDOFs.
  * \param[out] C       - Result of A*B. The transpose of the result is stored, dimension nDOFs X nPoints.
  */
  void MatMulTranspose(unsigned short nDOFs,
                       unsigned short nPoints,
                       vector<su2double> &A,
                       vector<su2double> &B,
                       vector<su2double> &C);

  /*!
  * \brief Function, which computes the value of the Jacobi polynomial for the given x-coordinate.
  * \param[in] n     - Order of the Jacobi polynomial.
  * \param[in] alpha - Alpha coefficient of the Jacobi polynomial.
  * \param[in] beta  - Beta coefficient of the Jacobi polynomial.
  * \param[in] x     - Coordinate (-1 <= x <= 1) for which the Jacobi polynomial must be evaluated.
  * \return            The value of the normalized Jacobi polynomial f order n for the given value of x.
  */
  su2double NormJacobi(unsigned short n,
                       unsigned short alpha,
                       unsigned short beta,
                       su2double      x);

  /*!
  * \brief Function, which computes the Vandermonde matrix for a standard 1D edge.
  * \param[in]  nDOFs - Number of DOFs, which in 1D is the polynomial degree + 1.
  * \param[in]  r     - Parametric coordinates for which the Vandermonde matrix must be computed.
  * \param[out] V     - Matrix to store the Vandermonde matrix in all r-locations.
  */
  void Vandermonde1D(unsigned short          nDOFs,
                     const vector<su2double> &r,
                     vector<su2double>       &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard triangle.
   * \param[in]  nPoly - Polynomial degree of the triangle.
   * \param[in]  nDOFs - Number of DOFs of the triangle.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r- and s-locations.
  */
  void Vandermonde2D_Triangle(unsigned short          nPoly,
                              unsigned short          nDOFs,
                              const vector<su2double> &r,
                              const vector<su2double> &s,
                              vector<su2double>       &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard quadrilateral.
   * \param[in]  nPoly - Polynomial degree of the quadrilateral.
   * \param[in]  nDOFs - Number of DOFs of the quadrilateral.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r- and s-locations.
  */
  void Vandermonde2D_Quadrilateral(unsigned short          nPoly,
                                   unsigned short          nDOFs,
                                   const vector<su2double> &r,
                                   const vector<su2double> &s,
                                   vector<su2double>       &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  nPoly - Polynomial degree of the tetrahedron.
   * \param[in]  nDOFs - Number of DOFs of the tetrahedron.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Tetrahedron(unsigned short          nPoly,
                                 unsigned short          nDOFs,
                                 const vector<su2double> &r,
                                 const vector<su2double> &s,
                                 const vector<su2double> &t,
                                 vector<su2double>       &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard pyramid.
   * \param[in]  nPoly - Polynomial degree of the pyramid.
   * \param[in]  nDOFs - Number of DOFs of the pyramid.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Pyramid(unsigned short          nPoly,
                             unsigned short          nDOFs,
                             const vector<su2double> &r,
                             const vector<su2double> &s,
                             const vector<su2double> &t,
                             vector<su2double>       &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard prism.
   * \param[in]  nPoly - Polynomial degree of the prism.
   * \param[in]  nDOFs - Number of DOFs of the prism.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Prism(unsigned short          nPoly,
                           unsigned short          nDOFs,
                           const vector<su2double> &r,
                           const vector<su2double> &s,
                           const vector<su2double> &t,
                           vector<su2double>       &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard hexahedron.
   * \param[in]  nPoly - Polynomial degree of the hexahedron.
   * \param[in]  nDOFs - Number of DOFs of the hexahedron.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Hexahedron(unsigned short          nPoly,
                                unsigned short          nDOFs,
                                const vector<su2double> &r,
                                const vector<su2double> &s,
                                const vector<su2double> &t,
                                vector<su2double>       &V);
};

/*!
 * \class FEMStandardElementClass
 * \brief Class to define a FEM standard element.
 * \author E. van der Weide
 * \version 4.1.2 "Cardinal"
 */
class FEMStandardElementClass : public FEMStandardElementBaseClass {
private:
  unsigned short nPoly;        /*!< \brief Polynomial degree of the element. */
  unsigned short nDOFs;        /*!< \brief Number of DOFs of the element. */

  vector<su2double> rDOFs;     /*!< \brief r-location of the DOFs for this standard element. */
  vector<su2double> sDOFs;     /*!< \brief s-location of the DOFs for this standard element, if needed. */
  vector<su2double> tDOFs;     /*!< \brief t-location of the DOFs for this standard element, if needed. */

  vector<su2double> lagBasisIntegration; /*!< \brief Lagrangian basis functions in the integration points. */

  vector<su2double> drLagBasisIntegration; /*!< \brief r-derivatives of the Lagrangian basis functions in the integration points. */
  vector<su2double> dsLagBasisIntegration; /*!< \brief s-derivatives of the Lagrangian basis functions in the integration points. */
  vector<su2double> dtLagBasisIntegration; /*!< \brief t-derivatives of the Lagrangian basis functions in the integration points. */

  vector<unsigned short> connFace0; /*!< \brief Local connectivity of face 0 of the element. The numbering of the DOFs is
                                                such that the element is to the left of the face. */
  vector<unsigned short> connFace1; /*!< \brief Local connectivity of face 1 of the element. The numbering of the DOFs is
                                                such that the element is to the left of the face. */
  vector<unsigned short> connFace2; /*!< \brief Local connectivity of face 2 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */
  vector<unsigned short> connFace3; /*!< \brief Local connectivity of face 3 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */
  vector<unsigned short> connFace4; /*!< \brief Local connectivity of face 4 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */
  vector<unsigned short> connFace5; /*!< \brief Local connectivity of face 5 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */

  vector<unsigned short> subConn1ForPlotting; /*!< \brief Local subconnectivity of element type 1 of the high order element.
                                                          Used for plotting. */
  vector<unsigned short> subConn2ForPlotting; /*!< \brief Local subconnectivity of element type 2 of the high order element.
                                                          Used for plotting. */
public:
  /*!
  * \brief Standard Constructor. Nothing to be done.
  */
  FEMStandardElementClass();

  /*!
  * \brief Destructor. Nothing to be done, because the vectors are deleted automatically.
  */
  ~FEMStandardElementClass();

  /*!
  * \brief Alternative constructor.
  * \param[in] val_VTK_Type   - Type of the element using the VTK convention.
  * \param[in] val_nPoly      - Polynomial degree of the element.
  * \param[in] val_constJac   - Whether or not the Jacobians are constant.
  * \param[in] config         - Object, which contains the input parameters.
  * \param[in] val_orderExact - Default argument. If specified, it contains the
                                order of the polynomials that must be integrated
                                exactly by the integration rule.
  */
  FEMStandardElementClass(unsigned short val_VTK_Type,
                          unsigned short val_nPoly,
                          bool           val_constJac,
                          CConfig        *config,
                          unsigned short val_orderExact = 0);
  /*!
  * \brief Copy constructor.
  * \param[in] other - Object, whose data must be copied.
  */
  FEMStandardElementClass(const FEMStandardElementClass &other);

  /*!
  * \brief Assignment operator.
  * \param[in] other - Object, to which this object must be assigned.
  * \return The current object, after the member variables were assigned the correct value.
  */
  FEMStandardElementClass& operator=(const FEMStandardElementClass &other);

  /*!
  * \brief Function, which makes available the r-derivatives of the basis functions in the integration points.
  * \return  The pointer to data, which stores the r-derivatives of the basis functions.
  */
  su2double *GetDrBasisFunctionsIntegration(void);

  /*!
  * \brief Function, which makes available the s-derivatives of the basis functions in the integration points.
  * \return  The pointer to data, which stores the s-derivatives of the basis functions.
  */
  su2double *GetDsBasisFunctionsIntegration(void);

  /*!
  * \brief Function, which makes available the t-derivatives of the basis functions in the integration points.
  * \return  The pointer to data, which stores the t-derivatives of the basis functions.
  */
  su2double *GetDtBasisFunctionsIntegration(void);

  /*!
  * \brief Function, which makes available the connectivity of face 0.
  * \return  The pointer to data, which stores the connectivity of face 0.
  */
  unsigned short *GetConnFace0(void);

  /*!
  * \brief Function, which makes available the connectivity of face 1.
  * \return  The pointer to data, which stores the connectivity of face 1.
  */
  unsigned short *GetConnFace1(void);

  /*!
  * \brief Function, which makes available the connectivity of face 2.
  * \return  The pointer to data, which stores the connectivity of face 2.
  */
  unsigned short *GetConnFace2(void);

  /*!
  * \brief Function, which makes available the connectivity of face 3.
  * \return  The pointer to data, which stores the connectivity of face 3.
  */
  unsigned short *GetConnFace3(void);

  /*!
  * \brief Function, which makes available the connectivity of face 4.
  * \return  The pointer to data, which stores the connectivity of face 4.
  */
  unsigned short *GetConnFace4(void);

  /*!
  * \brief Function, which makes available the connectivity of face 5.
  * \return  The pointer to data, which stores the connectivity of face 5.
  */
  unsigned short *GetConnFace5(void);

  /*!
  * \brief Function, which makes available the number of DOFs for this standard element.
  * \return  The number of DOFs of this standard element.
  */
  unsigned short GetNDOFs(void);

  /*!
  * \brief Function, which checks if the function arguments correspond to this standard element.
  * \param[in] val_VTK_Type - Type of the element using the VTK convention.
  * \param[in] val_nPoly    - Polynomial degree of the element.
  * \param[in] val_constJac - Whether or not the Jacobians are constant.
  * \return Whether or not the function arguments correspond to this standard element.
  */
  bool SameStandardElement(unsigned short val_VTK_Type,
                           unsigned short val_nPoly,
                           bool           val_constJac);

private:
  /*!
  * \brief Function, which changes the given quadrilateral connectivity, such that the direction coincides
           with the direction corresponding to corner vertices vert0, vert1, vert2, vert3.
  * \param[in,out] connQuad - Connectivity of the quadrilateral, that must be adapted.
  * \param[in]     vert0    - Corner vertex 0 of the desired sequence.
  * \param[in]     vert1    - Corner vertex 1 of the desired sequence.
  * \param[in]     vert2    - Corner vertex 2 of the desired sequence.
  * \param[in]     vert3    - Corner vertex 3 of the desired sequence.
  */
  void ChangeDirectionQuadConn(std::vector<unsigned short> &connQuad,
                               unsigned short              vert0,
                               unsigned short              vert1,
                               unsigned short              vert2,
                               unsigned short              vert3);

  /*!
  * \brief Function, which changes the given triangular connectivity, such that the direction coincides
           with the direction corresponding to corner vertices vert0, vert1, vert2.
  * \param[in,out] connTriangle - Connectivity of the triangle, that must be adapted.
  * \param[in]     vert0        - Corner vertex 0 of the desired sequence.
  * \param[in]     vert1        - Corner vertex 1 of the desired sequence.
  * \param[in]     vert2        - Corner vertex 2 of the desired sequence.
  */
  void ChangeDirectionTriangleConn(std::vector<unsigned short> &connTriangle,
                                   unsigned short              vert0,
                                   unsigned short              vert1,
                                   unsigned short              vert2);

  /*!
  * \brief Function, which copies the data of the given object into the current object.
  * \param[in] other - Object, whose data is copied.
  */
  void Copy(const FEMStandardElementClass &other);

  /*!
  * \brief Function, which creates all the data for a line element.
  */
  void DataStandardLine(void);

  /*!
  * \brief Function, which creates all the data for a triangular element.
  */
  void DataStandardTriangle(void);

  /*!
  * \brief Function, which creates all the data for a quadrilateral element.
  */
  void DataStandardQuadrilateral(void);

  /*!
  * \brief Function, which creates all the data for a tetrahedral element.
  */
  void DataStandardTetrahedron(void);

  /*!
  * \brief Function, which creates all the data for a pyramid element.
  */
  void DataStandardPyramid(void);

  /*!
  * \brief Function, which creates all the data for a prism element.
  */
  void DataStandardPrism(void);

  /*!
  * \brief Function, which creates all the data for a hexahedral element.
  */
  void DataStandardHexahedron(void);

  /*!
  * \brief Function, which determines the connectivity of the linear subtetrahedra for a high
           order tetrahedron.
  */
  void SubConnTetrahedron(void);

  /*!
  * \brief Function, which determines the connectivity of the linear subpyramids and
           subtetrahedra for a high order pyramid.
  */
  void SubConnPyramid(void);

  /*!
  * \brief Function, which determines the connectivity of the linear subprisms for a high
           order prism.
  */
  void SubConnPrism(void);

  /*!
  * \brief Function, which determines the connectivity of the linear subhexahedra for a high
           order hexahedron.
  */
  void SubConnHexahedron(void);
};

/*!
 * \class FEMStandardInternalFaceClass
 * \brief Class to define a FEM standard internal face.
 * \author E. van der Weide
 * \version 4.1.2 "Cardinal"
 */
class FEMStandardInternalFaceClass : public FEMStandardElementBaseClass {
private:
  unsigned short nDOFsFaceSide0;    /*!< \brief Number of DOFs on side 0 of the face. */
  unsigned short nDOFsFaceSide1;    /*!< \brief Number of DOFs on side 1 of the face. */

  unsigned short nPolyElemSide0;    /*!< \brief Polynomial degree of the element on side 0 of the face. */
  unsigned short nPolyElemSide1;    /*!< \brief Polynomial degree of the element on side 1 of the face. */
  unsigned short nDOFsElemSide0;    /*!< \brief Number of DOFs of the element on side 0 of the face. */
  unsigned short nDOFsElemSide1;    /*!< \brief Number of DOFs of the element on side 1 of the face. */
  unsigned short VTK_TypeElemSide0; /*!< \brief Type of the element on side 0 of the face using the VTK convention. */
  unsigned short VTK_TypeElemSide1; /*!< \brief Type of the element on side 1 of the face using the VTK convention. */

  vector<su2double> rDOFsFaceSide0;   /*!< \brief r-location of the DOFs on side 0 of the face. */
  vector<su2double> rDOFsFaceSide1;   /*!< \brief r-location of the DOFs on side 1 of the face. */
  vector<su2double> sDOFsFaceSide0;   /*!< \brief s-location of the DOFs on side 0 of the face, if needed. */
  vector<su2double> sDOFsFaceSide1;   /*!< \brief s-location of the DOFs on side 1 of the face, if needed. */

  vector<su2double> lagBasisIntegrationSide0; /*!< \brief Lagrangian basis functions in the integration points
                                                          of side0 of the face. */
  vector<su2double> lagBasisIntegrationSide1; /*!< \brief Lagrangian basis functions in the integration points
                                                          of side1 of the face. */

  vector<su2double> drLagBasisIntegrationSide0; /*!< \brief r-derivatives of the Lagrangian basis functions in the
                                                            integration points of the element on side 0 of the face. */
  vector<su2double> drLagBasisIntegrationSide1; /*!< \brief r-derivatives of the Lagrangian basis functions in the
                                                            integration points of the element on side 1 of the face. */
  vector<su2double> dsLagBasisIntegrationSide0; /*!< \brief s-derivatives of the Lagrangian basis functions in the
                                                            integration points of the element on side 0 of the face. */
  vector<su2double> dsLagBasisIntegrationSide1; /*!< \brief s-derivatives of the Lagrangian basis functions in the
                                                            integration points of the element on side 1 of the face. */
  vector<su2double> dtLagBasisIntegrationSide0; /*!< \brief t-derivatives of the Lagrangian basis functions in the
                                                            integration points of the element on side 0 of the face. */
  vector<su2double> dtLagBasisIntegrationSide1; /*!< \brief t-derivatives of the Lagrangian basis functions in the
                                                            integration points of the element on side 1 of the face. */
public:
  /*!
  * \brief Standard Constructor. Nothing to be done.
  */
  FEMStandardInternalFaceClass();

  /*!
  * \brief Destructor. Nothing to be done, because the vectors are deleted automatically.
  */
  ~FEMStandardInternalFaceClass();

  /*!
  * \brief Alternative constructor.
  * \param[in] val_VTK_TypeFace    - The type of the face using the VTK convention.
  * \param[in] val_VTK_TypeSide0   - Type of the element adjacent to side 0 of the face
                                     using the VTK convention.
  * \param[in] val_nPolySide0      - Polynomial degree of the element adjacent to side 0
                                     of the face.
  * \param[in] val_VTK_TypeSide1   - Type of the element adjacent to side 1 of the face
                                     using the VTK convention.
  * \param[in] val_nPolySide1      - Polynomial degree of the element adjacent to side 1
                                     of the face.
  * \param[in] val_constJac        - Whether or not the Jacobians are constant.
  * \param[in] config              - Object, which contains the input parameters.
  * \param[in] val_orderExact      - Default argument. If specified, it contains the
                                     order of the polynomials that must be integrated
                                     exactly by the integration rule.
  */
  FEMStandardInternalFaceClass(unsigned short val_VTK_TypeFace,
                               unsigned short val_VTK_TypeSide0,
                               unsigned short val_nPolySide0,
                               unsigned short val_VTK_TypeSide1,
                               unsigned short val_nPolySide1,
                               bool           val_constJac,
                               CConfig        *config,
                               unsigned short val_orderExact = 0);
  /*!
  * \brief Copy constructor.
  * \param[in] other - Object, whose data must be copied.
  */
  FEMStandardInternalFaceClass(const FEMStandardInternalFaceClass &other);

  /*!
  * \brief Assignment operator.
  * \param[in] other - Object, to which this object must be assigned.
  * \return The current object, after the member variables were assigned the correct value.
  */
  FEMStandardInternalFaceClass& operator=(const FEMStandardInternalFaceClass &other);

  /*!
  * \brief Function, which checks if the function arguments correspond to this standard face.
  * \param[in] val_VTK_TypeFace    - The type of the face using the VTK convention.
  * \param[in] val_constJac        - Whether or not the Jacobians are constant.
  * \param[in] val_VTK_TypeSide0   - Type of the element adjacent to side 0 of the face
                                     using the VTK convention.
  * \param[in] val_nPolySide0      - Polynomial degree of the element adjacent to side 0
                                     of the face.
  * \param[in] val_VTK_TypeSide1   - Type of the element adjacent to side 1 of the face
                                     using the VTK convention.
  * \param[in] val_nPolySide1      - Polynomial degree of the element adjacent to side 1
                                     of the face.
  */
  bool SameStandardMatchingFace(unsigned short val_VTK_TypeFace,
                                bool           val_constJac,
                                unsigned short val_VTK_TypeSide0,
                                unsigned short val_nPolySide0,
                                unsigned short val_VTK_TypeSide1,
                                unsigned short val_nPolySide1);
private:
  /*!
  * \brief Function, which copies the data of the given object into the current object.
  * \param[in] other - Object, whose data is copied.
  */
  void Copy(const FEMStandardInternalFaceClass &other);

  /*!
  * \brief Function, which computes the values of the derivatives of the basis functions
           of the adjacent elements in the integration points of the face.
  * \param[in]  VTK_TypeElem          - Type of the element adjacent to the face using the VTK convention.
  * \param[in]  nPolyElem             - Polynomial degree of the element adjacent to the face.
  * \param[out] nDOFsElem             - Number of DOFs of the element adjacent to the face.
  * \param[out] drLagBasisIntegration - r-derivatives of the basis functions in the integration points
                                        of the face.
  * \param[out] dsLagBasisIntegration - s-derivatives of the basis functions in the integration points
                                        of the face.
  * \param[out] dtLagBasisIntegration - t-derivatives of the basis functions in the integration points
                                        of the face.
  */
  void DerivativesBasisFunctionsAdjacentElement(unsigned short    VTK_TypeElem,
                                                unsigned short    nPolyElem,
                                                unsigned short    &nDOFsElem,
                                                vector<su2double> &drLagBasisIntegration,
                                                vector<su2double> &dsLagBasisIntegration,
                                                vector<su2double> &dtLagBasisIntegration);
};

#include "fem_standard_element.inl"
