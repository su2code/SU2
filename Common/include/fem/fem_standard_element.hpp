/*!
 * \file fem_standard_element.hpp
 * \brief Headers of the main functions for the FEM standard elements.
 *        The functions are in the <i>fem_standard_element.cpp</i> file.
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

#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>

#include "../CConfig.hpp"

using namespace std;

/*!
 * \class CFEMStandardElementBase
 * \brief Base class for a FEM standard element.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CFEMStandardElementBase {
 protected:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short orderExact;   /*!< \brief Polynomial order that must be integrated exactly by the integration rule. */
  unsigned short nIntegration; /*!< \brief Number of points used in the numerical integration. */

  bool constJacobian; /*!< \brief Whether or not the element has a constant Jacobian. */

  vector<su2double> rIntegration; /*!< \brief r-location of the integration points for this standard element. */
  vector<su2double>
      sIntegration; /*!< \brief s-location of the integration points for this standard element, if needed. */
  vector<su2double>
      tIntegration; /*!< \brief t-location of the integration points for this standard element, if needed. */
  vector<su2double> wIntegration; /*!< \brief The weights of the integration points for this standard element. */

 public:
  /*!
   * \brief Constructor. Nothing to be done.
   */
  CFEMStandardElementBase() = default;

  /*!
   * \brief Destructor. Nothing to be done, because the vectors are deleted automatically.
   */
  virtual ~CFEMStandardElementBase() = default;

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
  CFEMStandardElementBase(unsigned short val_VTK_Type, unsigned short val_nPoly, bool val_constJac, CConfig* config,
                          unsigned short val_orderExact);

 public:
  /*!
   * \brief Function, which makes available the type of the element.
   * \return  The type of the element using the VTK convention.
   */
  inline unsigned short GetVTK_Type(void) const { return VTK_Type; }

  /*!
   * \brief Function, which makes available the weights in the integration points.
   * \return  The const pointer to data, which stores the weights in the integration points.
   */
  inline const su2double* GetWeightsIntegration(void) const { return wIntegration.data(); }

  /*!
  * \brief Static function, which makes available the number of DOFs for an element
           corresponding to the arguments.
  * \param[in] VTK_Type         - Type of the element using the VTK convention.
  * \param[in] nPoly            - Polynomial degree of the element.
  * \param[in] typeErrorMessage - Default argument used to write a good error message.
  * \return The number of DOFs
  */
  static unsigned short GetNDOFsStatic(unsigned short VTK_Type, unsigned short nPoly,
                                       unsigned long typeErrorMessage = 0);

  /*!
   * \brief Function, which makes available the number of integration points for this standard element.
   * \return  The number of integration points of this standard element.
   */
  inline unsigned short GetNIntegration(void) const { return nIntegration; }

  /*!
   * \brief Function, which makes available the polynomial order that must be integrated exactly.
   * \return  The polynomial order that must be integrated exactly.
   */
  inline unsigned short GetOrderExact(void) const { return orderExact; }

  /*!
   * \brief Static function, which computes the inverse of the given square matrix.
   * \param[in]     n - Number of rows/columns of the square matrix A.
   * \param[in,out] A - On input the square matrix to be inverted. On output the inverse.
   */
  static void InverseMatrix(unsigned short n, vector<su2double>& A);

  /*!
   * \brief Function, which computes the gradient of the Vandermonde matrix for a standard 1D edge.
   * \param[in]  nDOFs - Number of DOFs, which in 1D is the polynomial degree + 1.
   * \param[in]  r     - Parametric coordinates for which the gradient of the Vandermonde matrix must be computed.
   * \param[out] VDr   - Matrix to store the gradient of the Vandermonde matrix in all r-locations.
   */
  void GradVandermonde1D(unsigned short nDOFs, const vector<su2double>& r, vector<su2double>& VDr);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard 1D edge.
   * \param[in]  nDOFs - Number of DOFs, which in 1D is the polynomial degree + 1.
   * \param[in]  r     - Parametric coordinates for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-locations.
   */
  void Vandermonde1D(unsigned short nDOFs, const vector<su2double>& r, vector<su2double>& V);

 protected:
  /*!
  * \brief Function, which checks if the sum of the given derivatives of the
           Lagrangian interpolation functions is 0 in the points.
  * \param[in] nPoints         - Number of points to be checked.
  * \param[in] nDOFs           - Number of DOFs of the element.
  * \param[in] dLagBasisPoints - Values of the derivatives of the Lagrangian
                                 interpolation functions in the given points.
  */
  void CheckSumDerivativesLagrangianBasisFunctions(const unsigned short nPoints, const unsigned short nDOFs,
                                                   const vector<su2double>& dLagBasisPoints);

  /*!
  * \brief Function, which checks if the sum of the given Lagrangian interpolation
           functions is 1 in the points.
  * \param[in]     nPoints        - Number of points to be checked.
  * \param[in]     nDOFs          - Number of DOFs of the element.
  * \param[in,out] lagBasisPoints - Values of the Lagrangian interpolation
                                    functions in the given points.
  */
  void CheckSumLagrangianBasisFunctions(const unsigned short nPoints, const unsigned short nDOFs,
                                        vector<su2double>& lagBasisPoints);

  /*!
   * \brief Function, which copies the data of the given object into the current object.
   * \param[in] other - Object, whose data is copied.
   */
  void Copy(const CFEMStandardElementBase& other);

  /*!
  * \brief Function, which computes the values of the derivatives of the basis functions
           of the adjacent elements in the integration points of the face.
  * \param[in]  VTK_TypeElem          - Type of the element adjacent to the face using the VTK convention.
  * \param[in]  nPolyElem             - Polynomial degree of the element adjacent to the face.
  * \param[in]  swapFaceInElement     - Whether or not the connectivity of the face must be swapped
                                        to get the desired numbering compared to the numbering used in the
                                        corresponding face of the element.
  * \param[out] nDOFsElem             - Number of DOFs of the element adjacent to the face.
  * \param[out] drLagBasisIntegration - r-derivatives of the basis functions in the integration points
                                        of the face.
  * \param[out] dsLagBasisIntegration - s-derivatives of the basis functions in the integration points
                                        of the face.
  * \param[out] dtLagBasisIntegration - t-derivatives of the basis functions in the integration points
                                        of the face.
  */
  void DerivativesBasisFunctionsAdjacentElement(unsigned short VTK_TypeElem, unsigned short nPolyElem,
                                                const bool swapFaceInElement, unsigned short& nDOFsElem,
                                                vector<su2double>& drLagBasisIntegration,
                                                vector<su2double>& dsLagBasisIntegration,
                                                vector<su2double>& dtLagBasisIntegration);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard triangle.
   * \param[in]  nPoly - Polynomial degree of the triangle.
   * \param[in]  nDOFs - Number of DOFs of the triangle.
   * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be
   * computed. \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix
   * must be computed. \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all
   * r- and s-locations. \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in
   * all r- and s-locations.
   */
  void GradVandermonde2D_Triangle(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                  const vector<su2double>& s, vector<su2double>& VDr, vector<su2double>& VDs);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard quadrilateral.
   * \param[in]  nPoly - Polynomial degree of the quadrilateral.
   * \param[in]  nDOFs - Number of DOFs of the quadrilateral.
   * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be
   * computed. \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix
   * must be computed. \param[out] VDr   - Matrix to store the gradient in r-direction of the Vandermonde matrix in all
   * r- and s-locations. \param[out] VDs   - Matrix to store the gradient in s-direction of the Vandermonde matrix in
   * all r- and s-locations.
   */
  void GradVandermonde2D_Quadrilateral(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                       const vector<su2double>& s, vector<su2double>& VDr, vector<su2double>& VDs);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  nPoly - Polynomial degree of the tetrahedron.
   * \param[in]  nDOFs - Number of DOFs of the tetrahedron.
   * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be
   * computed. \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix
   * must be computed. \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the
   * Vandermonde matrix must be computed. \param[out] VDr   - Matrix to store the gradient in r-direction of the
   * Vandermonde matrix in all r-, s- and t-locations. \param[out] VDs   - Matrix to store the gradient in s-direction
   * of the Vandermonde matrix in all r-, s- and t-locations. \param[out] VDt   - Matrix to store the gradient in
   * t-direction of the Vandermonde matrix in all r-, s- and t-locations.
   */
  void GradVandermonde3D_Tetrahedron(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                     const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& VDr,
                                     vector<su2double>& VDs, vector<su2double>& VDt);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard pyramid.
   * \param[in]  nPoly - Polynomial degree of the pyramid.
   * \param[in]  nDOFs - Number of DOFs of the pyramid.
   * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be
   * computed. \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix
   * must be computed. \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the
   * Vandermonde matrix must be computed. \param[out] VDr   - Matrix to store the gradient in r-direction of the
   * Vandermonde matrix in all r-, s- and t-locations. \param[out] VDs   - Matrix to store the gradient in s-direction
   * of the Vandermonde matrix in all r-, s- and t-locations. \param[out] VDt   - Matrix to store the gradient in
   * t-direction of the Vandermonde matrix in all r-, s- and t-locations.
   */
  void GradVandermonde3D_Pyramid(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                 const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& VDr,
                                 vector<su2double>& VDs, vector<su2double>& VDt);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard prism.
   * \param[in]  nPoly - Polynomial degree of the prism.
   * \param[in]  nDOFs - Number of DOFs of the prism.
   * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be
   * computed. \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix
   * must be computed. \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the
   * Vandermonde matrix must be computed. \param[out] VDr   - Matrix to store the gradient in r-direction of the
   * Vandermonde matrix in all r-, s- and t-locations. \param[out] VDs   - Matrix to store the gradient in s-direction
   * of the Vandermonde matrix in all r-, s- and t-locations. \param[out] VDt   - Matrix to store the gradient in
   * t-direction of the Vandermonde matrix in all r-, s- and t-locations.
   */
  void GradVandermonde3D_Prism(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                               const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& VDr,
                               vector<su2double>& VDs, vector<su2double>& VDt);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard hexahedron.
   * \param[in]  nPoly - Polynomial degree of the hexahedron.
   * \param[in]  nDOFs - Number of DOFs of the hexahedron.
   * \param[in]  r     - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be
   * computed. \param[in]  s     - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix
   * must be computed. \param[in]  t     - Parametric coordinate in t-direction for which the gradient of the
   * Vandermonde matrix must be computed. \param[out] VDr   - Matrix to store the gradient in r-direction of the
   * Vandermonde matrix in all r-, s- and t-locations. \param[out] VDs   - Matrix to store the gradient in s-direction
   * of the Vandermonde matrix in all r-, s- and t-locations. \param[out] VDt   - Matrix to store the gradient in
   * t-direction of the Vandermonde matrix in all r-, s- and t-locations.
   */
  void GradVandermonde3D_Hexahedron(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                    const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& VDr,
                                    vector<su2double>& VDs, vector<su2double>& VDt);

  /*!
  * \brief Function, which determines the values of the Lagrangian interpolation
           functions and its derivatives in the given set of points for a line.
  * \param[in]  nPoly            - Polynomial degree of the interpolation functions.
  * \param[in]  rPoints          - r-coordinates of the points where the
                                   interpolation functions must be evaluated.
  * \param[out] nDOFs            - Number of DOFs of the line element.
  * \param[out] rDOFs            - r-coordinates of the DOFs of the line element.
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesLine(const unsigned short nPoly, const vector<su2double>& rPoints,
                                                 unsigned short& nDOFs, vector<su2double>& rDOFs,
                                                 vector<su2double>& matVandermondeInv,
                                                 vector<su2double>& lagBasisPoints,
                                                 vector<su2double>& drLagBasisPoints);

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
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesTriangle(
      const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
      unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs, vector<su2double>& matVandermondeInv,
      vector<su2double>& lagBasisPoints, vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints);

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
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesQuadrilateral(
      const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
      unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs, vector<su2double>& matVandermondeInv,
      vector<su2double>& lagBasisPoints, vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints);

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
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesTetrahedron(
      const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
      const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
      vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
      vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints);

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
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesPyramid(
      const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
      const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
      vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
      vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints);

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
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesPrism(
      const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
      const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
      vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
      vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints);

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
  * \param[out] matVandermondeInv- Values of the inverse matrix of Vandermonde matrix
  * \param[out] lagBasisPoints   - Values of the Lagrangian interpolation
                                   functions in the given points.
  * \param[out] drLagBasisPoints - Values of the r-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dsLagBasisPoints - Values of the s-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  * \param[out] dtLagBasisPoints - Values of the t-derivatives of the Lagrangian
                                   interpolation functions in the given points.
  */
  void LagrangianBasisFunctionAndDerivativesHexahedron(
      const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
      const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
      vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
      vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints);

  /*!
  * \brief Function, which carries out a matrix matrix multiplication to obtain
           data in points and stores the result row major order.
  * \param[in]  nDOFs   - Dimension of the matrices. This typically corresponds
                          to the number of DOFs that are considered.
  * \param[in]  nPoints - Dimension of the matrices. This typically corresponds
                          to the number of integration points used.
  * \param[in]  A       - First matrix in the matrix matrix product A*B, dimension nPoints X nDOFs.
                          A is stored in column major order.
  * \param[in]  B       - Second matrix in the matrix matrix product A*B, dimension nDOFs X nDOFs.
                          B is stored in column major order.
  * \param[out] C       - Result of A*B, dimension nPoints X nDOFs. The result is stored in row
                          major order
  */
  void MatMulRowMajor(const unsigned short nDOFs, const unsigned short nPoints, const vector<su2double>& A,
                      const vector<su2double>& B, vector<su2double>& C);

  /*!
  * \brief Function, which computes the local connectivity of linear subelements of
           a line, which can be used for plotting.
  * \param[in]  nPoly    - Polynomial degree of the line.
  * \param[out] subConn  - The local subconnectivity of a line element.
  */
  void SubConnForPlottingLine(const unsigned short nPoly, vector<unsigned short>& subConn);

  /*!
  * \brief Function, which computes the local connectivity of linear subelements of
           a quadrilateral, which can be used for plotting.
  * \param[in]  nPoly    - Polynomial degree of the quadrilateral.
  * \param[out] subConn  - The local subconnectivity of a triangle element.
  */
  void SubConnForPlottingQuadrilateral(const unsigned short nPoly, vector<unsigned short>& subConn);

  /*!
  * \brief Function, which computes the local connectivity of linear subelements of
           a triangle, which can be used for plotting.
  * \param[in]  nPoly    - Polynomial degree of the triangle.
  * \param[out] subConn  - The local subconnectivity of a triangle element.
  */
  void SubConnForPlottingTriangle(const unsigned short nPoly, vector<unsigned short>& subConn);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard triangle.
   * \param[in]  nPoly - Polynomial degree of the triangle.
   * \param[in]  nDOFs - Number of DOFs of the triangle.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r- and s-locations.
   */
  void Vandermonde2D_Triangle(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                              const vector<su2double>& s, vector<su2double>& V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard quadrilateral.
   * \param[in]  nPoly - Polynomial degree of the quadrilateral.
   * \param[in]  nDOFs - Number of DOFs of the quadrilateral.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r- and s-locations.
   */
  void Vandermonde2D_Quadrilateral(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                   const vector<su2double>& s, vector<su2double>& V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  nPoly - Polynomial degree of the tetrahedron.
   * \param[in]  nDOFs - Number of DOFs of the tetrahedron.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
   */
  void Vandermonde3D_Tetrahedron(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                 const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard pyramid.
   * \param[in]  nPoly - Polynomial degree of the pyramid.
   * \param[in]  nDOFs - Number of DOFs of the pyramid.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
   */
  void Vandermonde3D_Pyramid(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                             const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard prism.
   * \param[in]  nPoly - Polynomial degree of the prism.
   * \param[in]  nDOFs - Number of DOFs of the prism.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
   */
  void Vandermonde3D_Prism(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                           const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard hexahedron.
   * \param[in]  nPoly - Polynomial degree of the hexahedron.
   * \param[in]  nDOFs - Number of DOFs of the hexahedron.
   * \param[in]  r     - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t     - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
   */
  void Vandermonde3D_Hexahedron(unsigned short nPoly, unsigned short nDOFs, const vector<su2double>& r,
                                const vector<su2double>& s, const vector<su2double>& t, vector<su2double>& V);

  /*!
  * \brief Function, which computes the constant in the penalty terms for a
           a viscous discretization.
  * \param[in] VTK_TypeElem - The element type, using the VTK convention,
                              adjacent to the face.
  * \param[in] nPolyElem    - The polynomial degree of the adjacent element.
  * \return                   The value of the viscous penalty parameter.
  */
  su2double ViscousPenaltyParameter(const unsigned short VTK_TypeElem, const unsigned short nPolyElem) const;

 private:
  /*!
   * \brief Function, which determines the 1D Gauss Legendre integration points and weights.
   * \param[in,out] GLPoints  - The location of the Gauss-Legendre integration points.
   * \param[in,out] GLWeights - The weights of the Gauss-Legendre integration points.
   */
  void GaussLegendrePoints1D(vector<su2double>& GLPoints, vector<su2double>& GLWeights);

  /*!
   * \brief Function, which computes the value of the gradient of the Jacobi polynomial for the given x-coordinate.
   * \param[in] n     - Order of the Jacobi polynomial.
   * \param[in] alpha - Alpha coefficient of the Jacobi polynomial.
   * \param[in] beta  - Beta coefficient of the Jacobi polynomial.
   * \param[in] x     - Coordinate (-1 <= x <= 1) for which the gradient of the Jacobi polynomial must be evaluated.
   * \return            The value of the gradient of the normalized Jacobi polynomial f order n for the given value of
   * x.
   */
  su2double GradNormJacobi(unsigned short n, unsigned short alpha, unsigned short beta, su2double x);

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
   * \brief Function, which computes the value of the Jacobi polynomial for the given x-coordinate.
   * \param[in] n     - Order of the Jacobi polynomial.
   * \param[in] alpha - Alpha coefficient of the Jacobi polynomial.
   * \param[in] beta  - Beta coefficient of the Jacobi polynomial.
   * \param[in] x     - Coordinate (-1 <= x <= 1) for which the Jacobi polynomial must be evaluated.
   * \return            The value of the normalized Jacobi polynomial f order n for the given value of x.
   */
  su2double NormJacobi(unsigned short n, unsigned short alpha, unsigned short beta, su2double x);
};

/*!
 * \class CFEMStandardElement
 * \brief Class to define a FEM standard element.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CFEMStandardElement : public CFEMStandardElementBase {
 private:
  unsigned short nPoly; /*!< \brief Polynomial degree of the element. */
  unsigned short nDOFs; /*!< \brief Number of DOFs of the element. */

  unsigned short VTK_Type1; /*!< \brief VTK type for elements of type 1 in subConn1ForPlotting. */
  unsigned short VTK_Type2; /*!< \brief VTK type for elements of type 2 in subConn2ForPlotting. */

  vector<su2double> rDOFs; /*!< \brief r-location of the DOFs for this standard element. */
  vector<su2double> sDOFs; /*!< \brief s-location of the DOFs for this standard element, if needed. */
  vector<su2double> tDOFs; /*!< \brief t-location of the DOFs for this standard element, if needed. */

  vector<su2double> lagBasisIntegration;      /*!< \brief Lagrangian basis functions in the integration points. */
  vector<su2double> lagBasisIntegrationTrans; /*!< \brief Transpose of lagBasisIntegration. It is stored such that
                                                          in the ADER-DG predictor step the residual is obtained
                                                          by one matrix multiplication. */
  vector<su2double> lagBasisSolDOFs; /*!< \brief Lagrangian basis functions in the solution DOFs. Only different
                                                 from 1 if the polynomial degree of the grid and solution differs. */

  vector<su2double>
      drLagBasisIntegration; /*!< \brief r-derivatives of the Lagrangian basis functions in the integration points. */
  vector<su2double>
      dsLagBasisIntegration; /*!< \brief s-derivatives of the Lagrangian basis functions in the integration points. */
  vector<su2double>
      dtLagBasisIntegration; /*!< \brief t-derivatives of the Lagrangian basis functions in the integration points. */

  vector<su2double> matVandermondeInv; /*!< \brief Inverse matrix of Vandermonde matrix in the DOFs for this standard
                                          element. This data is needed for the computation of shock sensing. */

  vector<su2double>
      matBasisIntegration; /*!< \brief Matrix of lagBasisIntegration, drLagBasisIntegration, dsLagBasisIntegration
                                       and dtLagBasisIntegration combined for efficiency when using BLAS routines. */
  vector<su2double> matDerBasisIntTrans; /*!< \brief Matrix of the transpose of the derivative part of
                                            matBasisIntegration. It is stored such that the volume residual can be
                                            computed in one matrix multiplication. */
  vector<su2double> matDerBasisSolDOFs;  /*!< \brief Matrix of the derivatives of the Lagrangian basis functions in the
                                            solution  DOFs. Needed to compute the metric terms in the solution DOFs. */
  vector<su2double> matDerBasisOwnDOFs;  /*!< \brief Matrix of the derivatives of the Lagrangian basis functions in the
                                            owned  DOFs. This differs from matDerBasisSolDOFs when the grid DOFs and the
                                                     solution DOFs do not coincide. This data is needed for the
                                            computation  of the derivatives of the metric terms. */
  vector<su2double> mat2ndDerBasisInt;   /*!< \brief Matrix which contains all possible second derivatives of the basis
                                            functions   in the integration points. As such second derivatives can be
                                            computed   using one call to the BLAS routines. */

  vector<unsigned short> connFace0; /*!< \brief Local connectivity of face 0 of the element. The numbering of the DOFs
                                       is such that the element is to the left of the face. */
  vector<unsigned short> connFace1; /*!< \brief Local connectivity of face 1 of the element. The numbering of the DOFs
                                       is such that the element is to the left of the face. */
  vector<unsigned short> connFace2; /*!< \brief Local connectivity of face 2 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */
  vector<unsigned short> connFace3; /*!< \brief Local connectivity of face 3 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */
  vector<unsigned short> connFace4; /*!< \brief Local connectivity of face 4 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */
  vector<unsigned short> connFace5; /*!< \brief Local connectivity of face 5 of the element, if present. The numbering
                                                of the DOFs is such that the element is to the left of the face. */

  vector<unsigned short> subConn1ForPlotting; /*!< \brief Local subconnectivity of element type 1 of the high order
                                                 element. Used for plotting. */
  vector<unsigned short> subConn2ForPlotting; /*!< \brief Local subconnectivity of element type 2 of the high order
                                                 element. Used for plotting. */
 public:
  /*!
  * \brief Alternative constructor.
  * \param[in] val_VTK_Type   - Type of the element using the VTK convention.
  * \param[in] val_nPoly      - Polynomial degree of the element.
  * \param[in] val_constJac   - Whether or not the Jacobians are constant.
  * \param[in] config         - Object, which contains the input parameters.
  * \param[in] val_orderExact - Default argument. If specified, it contains the
                                order of the polynomials that must be integrated
                                exactly by the integration rule.
  * \param[in] rLocSolDOFs    - Default argument. If specified, it contains the
                                parametric r location of the solution DOFs.
  * \param[in] sLocSolDOFs    - Default argument. If specified, it contains the
                                parametric s location of the solution DOFs.
  * \param[in] tLocSolDOFs    - Default argument. If specified, it contains the
                                parametric t location of the solution DOFs.
  */
  CFEMStandardElement(unsigned short val_VTK_Type, unsigned short val_nPoly, bool val_constJac, CConfig* config,
                      unsigned short val_orderExact = 0, const vector<su2double>* rLocSolDOFs = nullptr,
                      const vector<su2double>* sLocSolDOFs = nullptr, const vector<su2double>* tLocSolDOFs = nullptr);
  /*!
   * \brief Copy constructor.
   * \param[in] other - Object, whose data must be copied.
   */
  CFEMStandardElement(const CFEMStandardElement& other) : CFEMStandardElementBase(other) { Copy(other); }

  /*!
   * \brief Assignment operator.
   * \param[in] other - Object, to which this object must be assigned.
   * \return The current object, after the member variables were assigned the correct value.
   */
  CFEMStandardElement& operator=(const CFEMStandardElement& other) {
    Copy(other);
    return (*this);
  }

  /*!
  * \brief Function, which computes the Lagrangian basis functions for the
           given parametric coordinates.
  * \param[in]  parCoor   - Parametric coordinates for which the basis functions
                            and derivatives must be computed.
  * \param[out] lagBasis  - The values of the Lagrangian basis functions in parCoor.
  */
  void BasisFunctionsInPoint(const su2double* parCoor, vector<su2double>& lagBasis);
  /*!
  * \brief Function, which computes the Lagrangian basis functions and its
           derivatives for the given parametric coordinates.
  * \param[in]  parCoor   - Parametric coordinates for which the basis functions
                            and derivatives must be computed.
  * \param[out] lagBasis  - The values of the Lagrangian basis functions in parCoor.
  * \param[out] dLagBasis - The values of the derivatives of the basis functions
                            in parCoor.
  */
  void BasisFunctionsAndDerivativesInPoint(const su2double* parCoor, vector<su2double>& lagBasis,
                                           vector<vector<su2double> >& dLagBasis);

  /*!
   * \brief Function, which makes available the values of the basis functions in the integration points.
   * \return The pointer to data, which stores the basis functions in the integration points.
   */
  inline su2double* GetBasisFunctionsIntegration(void) { return lagBasisIntegration.data(); }

  /*!
   * \brief Function, which makes available the transpose of the basis functions in the integration points.
   * \return The pointer to data, which stores the transpose matrix of the basis functions.
   */
  inline const su2double* GetBasisFunctionsIntegrationTrans(void) const { return lagBasisIntegrationTrans.data(); }

  /*!
   * \brief Function, which makes available the values of the basis functions in the solution DOFs.
   * \return The pointer to data, which stores the basis functions in the solution DOFs.
   */
  inline const su2double* GetBasisFunctionsSolDOFs(void) const { return lagBasisSolDOFs.data(); }

  /*!
   * \brief Function, which makes available the r-derivatives of the basis functions in the integration points.
   * \return  The pointer to data, which stores the r-derivatives of the basis functions.
   */
  inline su2double* GetDrBasisFunctionsIntegration(void) { return drLagBasisIntegration.data(); }

  /*!
   * \brief Function, which makes available the s-derivatives of the basis functions in the integration points.
   * \return  The pointer to data, which stores the s-derivatives of the basis functions.
   */
  inline su2double* GetDsBasisFunctionsIntegration(void) { return dsLagBasisIntegration.data(); }

  /*!
   * \brief Function, which makes available the t-derivatives of the basis functions in the integration points.
   * \return  The pointer to data, which stores the t-derivatives of the basis functions.
   */
  inline su2double* GetDtBasisFunctionsIntegration(void) { return dtLagBasisIntegration.data(); }

  /*!
   * \brief Function, which makes available the matrix storage of the inverse of Vandermonde matrix of solution DOFs.
   * \return  The pointer to matVandermondeInv.
   */
  inline const su2double* GetMatVandermondeInv(void) const { return matVandermondeInv.data(); }

  /*!
   * \brief Function, which makes available the matrix storage of the basis functions in the integration points.
   * \return  The pointer to matBasisIntegration.
   */
  inline const su2double* GetMatBasisFunctionsIntegration(void) const { return matBasisIntegration.data(); }

  /*!
   * \brief Function, which makes available the transpose matrix of the derivative of the basis functions in the
   * integration points. \return  The pointer to matDerBasisIntTrans;
   */
  inline const su2double* GetDerMatBasisFunctionsIntTrans(void) const { return matDerBasisIntTrans.data(); }

  /*!
   * \brief Function, which makes available the matrix storage of the derivative of the basis functions in the own DOFs.
   * \return  The pointer to matDerBasisOwnDOFs.
   */
  inline const su2double* GetMatDerBasisFunctionsOwnDOFs(void) const { return matDerBasisOwnDOFs.data(); }

  /*!
   * \brief Function, which makes available the matrix storage of the derivative of the basis functions in the solution
   * DOFs. \return  The pointer to matDerBasisSolDOFs.
   */
  inline const su2double* GetMatDerBasisFunctionsSolDOFs(void) const { return matDerBasisSolDOFs.data(); }

  /*!
  * \brief Function, which makes available the matrix storage of the second derivativex of the basis functions
           in the integration points.
  * \return  The pointer to mat2ndDerBasisInt.
  */
  inline const su2double* GetMat2ndDerBasisFunctionsInt(void) const { return mat2ndDerBasisInt.data(); }

  /*!
   * \brief Function, which makes available the connectivity of face 0.
   * \return  The pointer to data, which stores the connectivity of face 0.
   */
  inline unsigned short* GetConnFace0(void) { return connFace0.data(); }

  /*!
   * \brief Function, which makes available the connectivity of face 1.
   * \return  The pointer to data, which stores the connectivity of face 1.
   */
  inline unsigned short* GetConnFace1(void) { return connFace1.data(); }

  /*!
   * \brief Function, which makes available the connectivity of face 2.
   * \return  The pointer to data, which stores the connectivity of face 2.
   */
  inline unsigned short* GetConnFace2(void) { return connFace2.data(); }

  /*!
   * \brief Function, which makes available the connectivity of face 3.
   * \return  The pointer to data, which stores the connectivity of face 3.
   */
  inline unsigned short* GetConnFace3(void) { return connFace3.data(); }

  /*!
   * \brief Function, which makes available the connectivity of face 4.
   * \return  The pointer to data, which stores the connectivity of face 4.
   */
  inline unsigned short* GetConnFace4(void) { return connFace4.data(); }

  /*!
   * \brief Function, which makes available the connectivity of face 5.
   * \return  The pointer to data, which stores the connectivity of face 5.
   */
  inline unsigned short* GetConnFace5(void) { return connFace5.data(); }

  /*!
   * \brief Function, which makes available the number of DOFs for this standard element.
   * \return  The number of DOFs of this standard element.
   */
  inline unsigned short GetNDOFs(void) const { return nDOFs; }

  /*!
   * \brief Function, which makes available the polynomial degree for this standard element.
   * \return  The polynomial degree of this standard element.
   */
  inline unsigned short GetNPoly(void) const { return nPoly; }

  /*!
   * \brief Function, which makes available the type of the element in subConn1ForPlotting.
   * \return  The type of the elements in subConn1ForPlotting using the VTK convention.
   */
  inline unsigned short GetVTK_Type1(void) const { return VTK_Type1; }

  /*!
   * \brief Function, which makes available the number of sub-elements of type 1 for plotting.
   * \return  The number of sub-elements of type 1 for plotting.
   */
  inline unsigned short GetNSubElemsType1(void) const {
    return subConn1ForPlotting.size() / GetNDOFsPerSubElem(GetVTK_Type1());
  }

  /*!
   * \brief Function, which makes available the the connectivity of the linear elements of type 1 as a const pointer.
   * \return  The pointer to the local connectivity of the linear elements of type 1.
   */
  inline const unsigned short* GetSubConnType1(void) const { return subConn1ForPlotting.data(); }

  /*!
   * \brief Function, which makes available the type of the element in subConn2ForPlotting.
   * \return  The type of the elements in subConn2ForPlotting using the VTK convention.
   */
  inline unsigned short GetVTK_Type2(void) const { return VTK_Type2; }

  /*!
   * \brief Function, which makes available the number of sub-elements of type 2 for plotting.
   * \return  The number of sub-elements of type 2 for plotting.
   */
  inline unsigned short GetNSubElemsType2(void) const {
    return subConn2ForPlotting.size() / GetNDOFsPerSubElem(GetVTK_Type2());
  }

  /*!
   * \brief Function, which makes available the the connectivity of the linear elements of type 2 as a const pointer.
   * \return  The pointer to the local connectivity of the linear elements of type 2.
   */
  inline const unsigned short* GetSubConnType2(void) const { return subConn2ForPlotting.data(); }

  /*!
   * \brief Function, which makes available the number of DOFs of a linear element, used for plotting.
   * \return  The number of DOFs of the linear elements.
   */
  unsigned short GetNDOFsPerSubElem(unsigned short val_VTK_Type) const;

  /*!
   * \brief Function, which makes available the r-location of the DOFs as a
            const pointer to the vector.
   * \return  The address of the vector, which stores the r-location of the DOFs.
   */
  inline const vector<su2double>* GetRDOFs(void) const { return &rDOFs; }

  /*!
   * \brief Function, which makes available the s-location of the DOFs as a
            const pointer to the vector.
   * \return  The address of the vector, which stores the s-location of the DOFs.
   */
  inline const vector<su2double>* GetSDOFs(void) const { return &sDOFs; }

  /*!
   * \brief Function, which makes available the t-location of the DOFs as a
            const pointer to the vector.
   * \return  The address of the vector, which stores the t-location of the DOFs.
   */
  inline const vector<su2double>* GetTDOFs(void) const { return &tDOFs; }

  /*!
   * \brief Function, which checks if the function arguments correspond to this standard element.
   * \param[in] val_VTK_Type - Type of the element using the VTK convention.
   * \param[in] val_nPoly    - Polynomial degree of the element.
   * \param[in] val_constJac - Whether or not the Jacobians are constant.
   * \return Whether or not the function arguments correspond to this standard element.
   */
  bool SameStandardElement(unsigned short val_VTK_Type, unsigned short val_nPoly, bool val_constJac);

  /*!
  * \brief Function, which estimates the amount of work for an element of this
           type. This information is used to determine a well balanced partition.
  * \param[in] config - Object, which contains the input parameters.
  */
  su2double WorkEstimateMetis(CConfig* config);

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
  void ChangeDirectionQuadConn(vector<unsigned short>& connQuad, unsigned short vert0, unsigned short vert1,
                               unsigned short vert2, unsigned short vert3) const;

  /*!
  * \brief Function, which changes the given triangular connectivity, such that the direction coincides
           with the direction corresponding to corner vertices vert0, vert1, vert2.
  * \param[in,out] connTriangle - Connectivity of the triangle, that must be adapted.
  * \param[in]     vert0        - Corner vertex 0 of the desired sequence.
  * \param[in]     vert1        - Corner vertex 1 of the desired sequence.
  * \param[in]     vert2        - Corner vertex 2 of the desired sequence.
  */
  void ChangeDirectionTriangleConn(vector<unsigned short>& connTriangle, unsigned short vert0, unsigned short vert1,
                                   unsigned short vert2) const;

  /*!
   * \brief Function, which copies the data of the given object into the current object.
   * \param[in] other - Object, whose data is copied.
   */
  void Copy(const CFEMStandardElement& other);

  /*!
  * \brief Function, which creates the basis functions and the matrix containing
           the derivatives of the basis functions in the given location of the
           parametric coordinates.
  * \param[in]  rLoc         - r-locations of the given points.
  * \param[in]  sLoc         - s-locations of the given points.
  * \param[in]  tLoc         - t-locations of the given points, if relevant.
  * \param[out] lagBasis     - Lagrangian basis functions in the given
                               parametric locations.
  * \param[out] matDerBasis  - Matrix to store the derivatives of the basis
                               functions.
  */
  void CreateBasisFunctionsAndMatrixDerivatives(const vector<su2double>& rLoc, const vector<su2double>& sLoc,
                                                const vector<su2double>& tLoc, vector<su2double>& matVandermondeInv,
                                                vector<su2double>& lagBasis, vector<su2double>& matDerBasis);
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
 * \class CFEMStandardInternalFace
 * \brief Class to define a FEM standard internal face.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CFEMStandardInternalFace : public CFEMStandardElementBase {
 private:
  unsigned short nDOFsFaceSide0; /*!< \brief Number of DOFs on side 0 of the face. */
  unsigned short nDOFsFaceSide1; /*!< \brief Number of DOFs on side 1 of the face. */

  unsigned short nPolyElemSide0;    /*!< \brief Polynomial degree of the element on side 0 of the face. */
  unsigned short nPolyElemSide1;    /*!< \brief Polynomial degree of the element on side 1 of the face. */
  unsigned short nDOFsElemSide0;    /*!< \brief Number of DOFs of the element on side 0 of the face. */
  unsigned short nDOFsElemSide1;    /*!< \brief Number of DOFs of the element on side 1 of the face. */
  unsigned short VTK_TypeElemSide0; /*!< \brief Type of the element on side 0 of the face using the VTK convention. */
  unsigned short VTK_TypeElemSide1; /*!< \brief Type of the element on side 1 of the face using the VTK convention. */

  bool swapFaceInElementSide0; /*!< \brief Whether or not the connectivity of the face must
                                           be swapped compared to the face of the corresponding
                                           standard element on side 0 of the face. */
  bool swapFaceInElementSide1; /*!< \brief Whether or not the connectivity of the face must
                                           be swapped compared to the face of the corresponding
                                           standard element on side 1 of the face. */

  su2double penaltyConstantFace; /*!< \brief The constant of the penalty parameter of the face,
                                             which is used in the viscous discretization. */

  vector<su2double> rDOFsFaceSide0; /*!< \brief r-location of the DOFs on side 0 of the face. */
  vector<su2double> rDOFsFaceSide1; /*!< \brief r-location of the DOFs on side 1 of the face. */
  vector<su2double> sDOFsFaceSide0; /*!< \brief s-location of the DOFs on side 0 of the face, if needed. */
  vector<su2double> sDOFsFaceSide1; /*!< \brief s-location of the DOFs on side 1 of the face, if needed. */

  vector<su2double> lagBasisFaceIntegrationSide0; /*!< \brief Lagrangian basis functions in the integration points
                                                              of side0 of the face. */
  vector<su2double> lagBasisFaceIntegrationSide1; /*!< \brief Lagrangian basis functions in the integration points
                                                              of side1 of the face. */

  vector<su2double> lagBasisFaceIntegrationTransposeSide0; /*!< \brief Transpose of lagBasisFaceIntegrationSide0. */
  vector<su2double> lagBasisFaceIntegrationTransposeSide1; /*!< \brief Transpose of lagBasisFaceIntegrationSide1. */

  vector<su2double> drLagBasisFaceIntegrationSide0; /*!< \brief r-derivatives of the face Lagrangian basis functions
                                                                of side 0 in the integration points. */
  vector<su2double> drLagBasisFaceIntegrationSide1; /*!< \brief r-derivatives of the face Lagrangian basis functions
                                                                of side 1 in the integration points. */
  vector<su2double> dsLagBasisFaceIntegrationSide0; /*!< \brief s-derivatives of the face Lagrangian basis functions
                                                                of side 0 in the integration points. */
  vector<su2double> dsLagBasisFaceIntegrationSide1; /*!< \brief s-derivatives of the face Lagrangian basis functions
                                                                of side 1 in the integration points. */

  vector<su2double> drLagBasisElemIntegrationSide0; /*!< \brief r-derivatives of the element Lagrangian basis functions
                                                                of side 0 in the integration points. */
  vector<su2double> drLagBasisElemIntegrationSide1; /*!< \brief r-derivatives of the element Lagrangian basis functions
                                                                of side 1 in the integration points. */
  vector<su2double> dsLagBasisElemIntegrationSide0; /*!< \brief s-derivatives of the element Lagrangian basis functions
                                                                of side 0 in the integration points. */
  vector<su2double> dsLagBasisElemIntegrationSide1; /*!< \brief s-derivatives of the element Lagrangian basis functions
                                                                of side 1 in the integration points. */
  vector<su2double> dtLagBasisElemIntegrationSide0; /*!< \brief t-derivatives of the element Lagrangian basis functions
                                                                of side 0 in the integration points. */
  vector<su2double> dtLagBasisElemIntegrationSide1; /*!< \brief t-derivatives of the element Lagrangian basis functions
                                                                of side 1 in the integration points. */

  vector<su2double>
      matDerBasisElemIntegrationSide0; /*!< \brief Matrix of drLagBasisElemIntegrationSide0,
                                          dsLagBasisElemIntegrationSide0 and dtLagBasisElemIntegrationSide0 combined for
                                          efficiency when using BLAS routines. */
  vector<su2double>
      matDerBasisElemIntegrationSide1; /*!< \brief Matrix of drLagBasisElemIntegrationSide1,
                                          dsLagBasisElemIntegrationSide1 and dtLagBasisElemIntegrationSide1 combined for
                                          efficiency when using BLAS routines. */

  vector<su2double> matDerBasisElemIntegrationTransposeSide0; /*!< \brief Transpose of matDerBasisElemIntegrationSide0,
                                                                 such that the residuals of the symmetrizing terms can
                                                                 be computed with a single matrix multiplication. */
  vector<su2double> matDerBasisElemIntegrationTransposeSide1; /*!< \brief Transpose of matDerBasisElemIntegrationSide1,
                                                                 such that the residuals of the symmetrizing terms can
                                                                 be computed with a single matrix multiplication. */
 public:
  /*!
  * \brief Alternative constructor.
  * \param[in] val_VTK_TypeFace            - The type of the face using the VTK convention.
  * \param[in] val_VTK_TypeSide0           - Type of the element adjacent to side 0 of the face
                                             using the VTK convention.
  * \param[in] val_nPolySide0              - Polynomial degree of the element adjacent to side 0
                                             of the face.
  * \param[in] val_VTK_TypeSide1           - Type of the element adjacent to side 1 of the face
                                             using the VTK convention.
  * \param[in] val_nPolySide1              - Polynomial degree of the element adjacent to side 1
                                             of the face.
  * \param[in] val_constJac                - Whether or not the Jacobians are constant.
  * \param[in] val_swapFaceInElementSide0  - Whether or not the connectivity of the face must be
                                             swapped compared to the face of the corresponding
                                             standard element on side 0 of the face.
  * \param[in] val_swapFaceInElementSide1  - Whether or not the connectivity of the face must be
                                             swapped compared to the face of the corresponding
                                             standard element on side 1 of the face.
  * \param[in] config                      - Object, which contains the input parameters.
  * \param[in] val_orderExact              - Default argument. If specified, it contains the
                                             order of the polynomials that must be integrated
                                             exactly by the integration rule.
  */
  CFEMStandardInternalFace(unsigned short val_VTK_TypeFace, unsigned short val_VTK_TypeSide0,
                           unsigned short val_nPolySide0, unsigned short val_VTK_TypeSide1,
                           unsigned short val_nPolySide1, bool val_constJac, bool val_swapFaceInElementSide0,
                           bool val_swapFaceInElementSide1, CConfig* config, unsigned short val_orderExact = 0);
  /*!
   * \brief Copy constructor.
   * \param[in] other - Object, whose data must be copied.
   */
  CFEMStandardInternalFace(const CFEMStandardInternalFace& other) : CFEMStandardElementBase(other) { Copy(other); }

  /*!
   * \brief Assignment operator.
   * \param[in] other - Object, to which this object must be assigned.
   * \return The current object, after the member variables were assigned the correct value.
   */
  CFEMStandardInternalFace& operator=(const CFEMStandardInternalFace& other) {
    Copy(other);
    return (*this);
  }

  /*!
  * \brief Function, which makes available the r-derivatives of the elements
           basis functions of side 0 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDrBasisElemIntegrationSide0(void) { return drLagBasisElemIntegrationSide0.data(); }

  /*!
  * \brief Function, which makes available the r-derivatives of the elements
           basis functions of side 1 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDrBasisElemIntegrationSide1(void) { return drLagBasisElemIntegrationSide1.data(); }

  /*!
  * \brief Function, which makes available the s-derivatives of the elements
           basis functions of side 0 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDsBasisElemIntegrationSide0(void) { return dsLagBasisElemIntegrationSide0.data(); }

  /*!
  * \brief Function, which makes available the s-derivatives of the elements
           basis functions of side 1 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDsBasisElemIntegrationSide1(void) { return dsLagBasisElemIntegrationSide1.data(); }

  /*!
  * \brief Function, which makes available the t-derivatives of the elements
           basis functions of side 0 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDtBasisElemIntegrationSide0(void) { return dtLagBasisElemIntegrationSide0.data(); }

  /*!
  * \brief Function, which makes available the t-derivatives of the elements
           basis functions of side 1 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDtBasisElemIntegrationSide1(void) { return dtLagBasisElemIntegrationSide1.data(); }

  /*!
  * \brief Function, which makes available the matrix with the derivatives of
           the element basis functions of side 0 in the integration points.
  * \return  The const pointer to data, which stores this information.
  */
  inline const su2double* GetMatDerBasisElemIntegrationSide0(void) const {
    return matDerBasisElemIntegrationSide0.data();
  }

  /*!
  * \brief Function, which makes available the matrix with the derivatives of
           the element basis functions of side 1 in the integration points.
  * \return  The const pointer to data, which stores this information.
  */
  inline const su2double* GetMatDerBasisElemIntegrationSide1(void) const {
    return matDerBasisElemIntegrationSide1.data();
  }

  /*!
  * \brief Function, which makes available the transpose of the matrix with
           the derivatives of the element basis functions of side 0 in the
           integration points.
  * \return  The const pointer to data, which stores this information.
  */
  inline const su2double* GetMatDerBasisElemIntegrationTransposeSide0(void) const {
    return matDerBasisElemIntegrationTransposeSide0.data();
  }

  /*!
  * \brief Function, which makes available the transpose of the matrix with
           the derivatives of the element basis functions of side 1 in the
           integration points.
  * \return  The const pointer to data, which stores this information.
  */
  inline const su2double* GetMatDerBasisElemIntegrationTransposeSide1(void) const {
    return matDerBasisElemIntegrationTransposeSide1.data();
  }

  /*!
  * \brief Function, which makes available the face basis functions of side 0
           in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetBasisFaceIntegrationSide0(void) const { return lagBasisFaceIntegrationSide0.data(); }

  /*!
  * \brief Function, which makes available the face basis functions of side 1
           in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetBasisFaceIntegrationSide1(void) const { return lagBasisFaceIntegrationSide1.data(); }

  /*!
  * \brief Function, which makes available transpose matrix of the face basis
           functions of side 0 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetBasisFaceIntegrationTransposeSide0(void) const {
    return lagBasisFaceIntegrationTransposeSide0.data();
  }

  /*!
  * \brief Function, which makes available transpose matrix of the face basis
           functions of side 1 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetBasisFaceIntegrationTransposeSide1(void) const {
    return lagBasisFaceIntegrationTransposeSide1.data();
  }

  /*!
  * \brief Function, which makes available the r-derivatives of the face basis
           functions of side 0 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDrBasisFaceIntegrationSide0(void) { return drLagBasisFaceIntegrationSide0.data(); }

  /*!
  * \brief Function, which makes available the r-derivatives of the face basis
           functions of side 1 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDrBasisFaceIntegrationSide1(void) { return drLagBasisFaceIntegrationSide1.data(); }

  /*!
  * \brief Function, which makes available the s-derivatives of the face basis
           functions of side 0 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDsBasisFaceIntegrationSide0(void) { return dsLagBasisFaceIntegrationSide0.data(); }

  /*!
  * \brief Function, which makes available the s-derivatives of the face basis
           functions of side 1 in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline su2double* GetDsBasisFaceIntegrationSide1(void) { return dsLagBasisFaceIntegrationSide1.data(); }

  /*!
  * \brief Function, which makes available the number of DOFs of the element
           on side 0 of the face.
  * \return  The number of DOFs of the element on side 0.
  */
  inline unsigned short GetNDOFsElemSide0(void) const { return nDOFsElemSide0; }

  /*!
  * \brief Function, which makes available the number of DOFs of the element
           on side 1 of the face.
  * \return  The number of DOFs on side 1.
  */
  inline unsigned short GetNDOFsElemSide1(void) const { return nDOFsElemSide1; }

  /*!
   * \brief Function, which makes available the number of DOFs on side 0 of the face.
   * \return  The number of DOFs on side 0.
   */
  inline unsigned short GetNDOFsFaceSide0(void) const { return nDOFsFaceSide0; }

  /*!
   * \brief Function, which makes available the number of DOFs on side 1 of the face.
   * \return  The number of DOFs on side 1.
   */
  inline unsigned short GetNDOFsFaceSide1(void) const { return nDOFsFaceSide1; }

  /*!
   * \brief Function, which makes available the penalty constant for this standard face.
   * \return  The penalty constant.
   */
  inline su2double GetPenaltyConstant(void) const { return penaltyConstantFace; }

  /*!
  * \brief Function, which checks if the function arguments correspond to this standard face.
  * \param[in] val_VTK_TypeFace           - The type of the face using the VTK convention.
  * \param[in] val_constJac               - Whether or not the Jacobians are constant.
  * \param[in] val_VTK_TypeSide0          - Type of the element adjacent to side 0
                                            of the face using the VTK convention.
  * \param[in] val_nPolySide0             - Polynomial degree of the element adjacent
                                            to side 0 of the face.
  * \param[in] val_VTK_TypeSide1          - Type of the element adjacent to side 1
                                            of the face using the VTK convention.
  * \param[in] val_nPolySide1             - Polynomial degree of the element adjacent
                                            to side 1 of the face.
  * \param[in] val_swapFaceInElementSide0 - Whether or not the connectivity of the face must
                                            be swapped w.r.t. the connectivity of face of the
                                            element on side 0.
  * \param[in] val_swapFaceInElementSide1 - Whether or not the connectivity of the face must
                                            be swapped w.r.t. the connectivity of face of the
                                            element on side 1.
  */
  bool SameStandardMatchingFace(unsigned short val_VTK_TypeFace, bool val_constJac, unsigned short val_VTK_TypeSide0,
                                unsigned short val_nPolySide0, unsigned short val_VTK_TypeSide1,
                                unsigned short val_nPolySide1, bool val_swapFaceInElementSide0,
                                bool val_swapFaceInElementSide1);

  /*!
  * \brief Function, which estimates the amount of work for an element of this
           type. This information is used to determine a well balanced partition.
  * \param[in] config - Object, which contains the input parameters.
  */
  su2double WorkEstimateMetis(CConfig* config);

 private:
  /*!
   * \brief Function, which copies the data of the given object into the current object.
   * \param[in] other - Object, whose data is copied.
   */
  void Copy(const CFEMStandardInternalFace& other);
};

/*!
 * \class CFEMStandardBoundaryFace
 * \brief Class to define a FEM standard boundary face.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CFEMStandardBoundaryFace : public CFEMStandardElementBase {
 private:
  unsigned short nDOFsFace; /*!< \brief Number of DOFs of the face. */

  unsigned short nPolyElem;    /*!< \brief Polynomial degree of the element adjacent to the face. */
  unsigned short nDOFsElem;    /*!< \brief Number of DOFs of the element adjacent element to the face. */
  unsigned short VTK_TypeElem; /*!< \brief Type of the element adjacent to the face using the VTK convention. */

  bool swapFaceInElement; /*!< \brief Whether or not the connectivity of the face must be swapped compared
                                      to the face of the corresponding standard element adjacent to the face. */

  su2double penaltyConstantFace; /*!< \brief The constant of the penalty parameter of the face,
                                             which is used in the viscous discretization. */

  vector<su2double> rDOFsFace; /*!< \brief r-location of the DOFs of the face. */
  vector<su2double> sDOFsFace; /*!< \brief s-location of the DOFs of the face, if needed. */

  vector<su2double> lagBasisFaceIntegration;          /*!< \brief Lagrangian basis functions in the integration
                                                                  points of the face. */
  vector<su2double> lagBasisFaceIntegrationTranspose; /*!< \brief Transpose of lagBasisFaceIntegration. */

  vector<su2double> drLagBasisFaceIntegration; /*!< \brief r-derivatives of the face Lagrangian basis functions
                                                           in the integration points. */
  vector<su2double> dsLagBasisFaceIntegration; /*!< \brief s-derivatives of the face Lagrangian basis functions
                                                           in the integration points. */

  vector<su2double> drLagBasisElemIntegration; /*!< \brief r-derivatives of the Lagrangian basis functions in the
                                                           integration points of the element adjacent to the face. */
  vector<su2double> dsLagBasisElemIntegration; /*!< \brief s-derivatives of the Lagrangian basis functions in the
                                                           integration points of the element adjacent to the face. */
  vector<su2double> dtLagBasisElemIntegration; /*!< \brief t-derivatives of the Lagrangian basis functions in the
                                                           integration points of the element adjacent to the face. */

  vector<su2double> matDerBasisElemIntegration; /*!< \brief Matrix of drLagBasisElemIntegration,
                                                   dsLagBasisElemIntegration and dtLagBasisElemIntegration combined for
                                                   efficiency when using BLAS routines. */

  vector<su2double> matDerBasisElemIntegrationTranspose; /*!< \brief Transpose of matDerBasisElemIntegration, such that
                                                                     the residuals of the symmetrizing terms can be
                                                            computed with a single matrix multiplication. */

  vector<unsigned short> subConnForPlotting; /*!< \brief Local subconnectivity of the high order element.
                                                         Used for plotting. */
 public:
  /*!
  * \brief Alternative constructor.
  * \param[in] val_VTK_TypeFace       - The type of the face using the VTK convention.
  * \param[in] val_VTK_TypeElem       - Type of the element adjacent to the face using the VTK convention.
  * \param[in] val_nPolyElem          - Polynomial degree of the element adjacent to the face.
  * \param[in] val_constJac           - Whether or not the Jacobians are constant.
  * \param[in] val_swapFaceInElement  - Whether or not the connectivity of the face must be swapped compared
                                        to the face of the corresponding standard element adjacent to the face.
  * \param[in] config                 - Object, which contains the input parameters.
  * \param[in] val_orderExact         - Default argument. If specified, it contains the order of the
                                        polynomials that must be integrated exactly by the integration rule.
  */
  CFEMStandardBoundaryFace(unsigned short val_VTK_TypeFace, unsigned short val_VTK_TypeElem,
                           unsigned short val_nPolyElem, bool val_constJac, bool val_swapFaceInElement, CConfig* config,
                           unsigned short val_orderExact = 0);

  /*!
   * \brief Copy constructor.
   * \param[in] other - Object, whose data must be copied.
   */
  CFEMStandardBoundaryFace(const CFEMStandardBoundaryFace& other) : CFEMStandardElementBase(other) { Copy(other); }

  /*!
   * \brief Assignment operator.
   * \param[in] other - Object, to which this object must be assigned.
   * \return The current object, after the member variables were assigned the correct value.
   */
  CFEMStandardBoundaryFace& operator=(const CFEMStandardBoundaryFace& other) {
    Copy(other);
    return (*this);
  }

  /*!
  * \brief Function, which makes available the r-derivatives of the element
           basis functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetDrBasisElemIntegration(void) const { return drLagBasisElemIntegration.data(); }

  /*!
  * \brief Function, which makes available the s-derivatives of the element
           basis functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetDsBasisElemIntegration(void) const { return dsLagBasisElemIntegration.data(); }

  /*!
  * \brief Function, which makes available the t-derivatives of the element
           basis functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetDtBasisElemIntegration(void) const { return dtLagBasisElemIntegration.data(); }

  /*!
  * \brief Function, which makes available the matrix with the derivatives of
           the element basis functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetMatDerBasisElemIntegration(void) const { return matDerBasisElemIntegration.data(); }

  /*!
  * \brief Function, which makes available the transpose of the matrix with
           the derivatives of the element basis functions in the integration
           points.
  * \return  The const pointer to data, which stores this information.
  */
  inline const su2double* GetMatDerBasisElemIntegrationTranspose(void) const {
    return matDerBasisElemIntegrationTranspose.data();
  }

  /*!
  * \brief Function, which makes available the face basis functions in the
           integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetBasisFaceIntegration(void) const { return lagBasisFaceIntegration.data(); }

  /*!
  * \brief Function, which makes available transpose matrix of the face basis
           functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetBasisFaceIntegrationTranspose(void) const {
    return lagBasisFaceIntegrationTranspose.data();
  }

  /*!
  * \brief Function, which makes available the r-derivatives of the face basis
           functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetDrBasisFaceIntegration(void) const { return drLagBasisFaceIntegration.data(); }

  /*!
  * \brief Function, which makes available the s-derivatives of the face basis
           functions in the integration points.
  * \return  The pointer to data, which stores this information.
  */
  inline const su2double* GetDsBasisFaceIntegration(void) const { return dsLagBasisFaceIntegration.data(); }

  /*!
  * \brief Function, which makes available the number of DOFs of the
           adjacent element.
  * \return  The number of DOFs of the element.
  */
  inline unsigned short GetNDOFsElem(void) const { return nDOFsElem; }

  /*!
   * \brief Function, which makes available the number of DOFs of the face.
   * \return  The number of DOFs of the face.
   */
  inline unsigned short GetNDOFsFace(void) const { return nDOFsFace; }

  /*!
  * \brief Function, which makes available the number of linear subfaces used
           for plotting, among others.
  * \return  The number of linear subfaces of the face.
  */
  inline unsigned short GetNSubFaces(void) const { return subConnForPlotting.size() / GetNDOFsPerSubFace(); }

  /*!
  * \brief Function, which makes available the number of DOFs of a linear subface, used
           for plotting, among others, plotting.
  * \return  The number of DOFs of a linear subfaces of the face.
  */
  unsigned short GetNDOFsPerSubFace(void) const;

  /*!
   * \brief Function, which makes available the penalty constant for this standard face.
   * \return  The penalty constant.
   */
  inline su2double GetPenaltyConstant(void) const { return penaltyConstantFace; }

  /*!
  * \brief Function, which makes available the the connectivity of the linear subfaces
           as a const pointer.
  * \return  The pointer to the local connectivity of the linear subfaces.
  */
  inline const unsigned short* GetSubFaceConn(void) const { return subConnForPlotting.data(); }

  /*!
  * \brief Function, which checks if the function arguments correspond to this standard face.
  * \param[in] val_VTK_TypeFace   - The type of the face using the VTK convention.
  * \param[in] val_constJac       - Whether or not the Jacobians are constant.
  * \param[in] val_VTK_TypeElem   - Type of the element adjacent to the face using the VTK convention.
  * \param[in] val_nPolyElem      - Polynomial degree of the element adjacent to the face.
  * \param[in] val_swapFaceInElem - Whether or not the connectivity of the face must be swapped w.r.t.
                                    the connectivity of face of the adjacent element.
  */
  bool SameStandardBoundaryFace(unsigned short val_VTK_TypeFace, bool val_constJac, unsigned short val_VTK_TypeElem,
                                unsigned short val_nPolyElem, bool val_swapFaceInElem);
  /*!
  * \brief Function, which estimates the amount of work for an element of this
           type. This information is used to determine a well balanced partition.
  * \param[in] config - Object, which contains the input parameters.
  */
  su2double WorkEstimateMetis(CConfig* config);

  /*!
  * \brief Function, which estimates the additional amount of work for an element
           of this type when a wall function treatment is used. This information
           is used to determine a well balanced partition.
  * \param[in] config    - Object, which contains the input parameters.
  * \param[in] nPointsWF - Number of points to discretize the wall model.
  */
  su2double WorkEstimateMetisWallFunctions(CConfig* config, const unsigned short nPointsWF);

 private:
  /*!
   * \brief Function, which copies the data of the given object into the current object.
   * \param[in] other - Object, whose data is copied.
   */
  void Copy(const CFEMStandardBoundaryFace& other);
};
