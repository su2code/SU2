/*!
 * \file fem_standard_element.hpp
 * \brief Headers of the main functions for the FEM standard elements.
 *        The functions are in the <i>fem_standard_element.cpp</i> file.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
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
 * \class FEMStandardElementClass
 * \brief Class to define an FEM standard element.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class FEMStandardElementClass {
private:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPoly;        /*!< \brief Polynomial degree of the element. */
  unsigned short nDOFs;        /*!< \brief Number of DOFs of the element. */
  unsigned short orderExact;   /*!< \brief Polynomial order that must be integrated exactly by the integration rule. */
  unsigned short nIntegration; /*!< \brief Number of points used in the numerical integration. */

  bool constJacobian;          /*!< \brief Whether or not the element has a constant Jacobian. */

  vector<su2double> rDOFs;     /*!< \brief r-location of the DOFs for this standard element. */
  vector<su2double> sDOFs;     /*!< \brief s-location of the DOFs for this standard element, if needed. */
  vector<su2double> tDOFs;     /*!< \brief t-location of the DOFs for this standard element, if needed. */

  vector<su2double> rIntegration; /*!< \brief r-location of the integration points for this standard element. */
  vector<su2double> sIntegration; /*!< \brief s-location of the integration points for this standard element, if needed. */
  vector<su2double> tIntegration; /*!< \brief t-location of the integration points for this standard element, if needed. */
  vector<su2double> wIntegration; /*!< \brief The weights of the integration points for this standard element. */

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
  * \param[in] val_VTK_Type - Type of the element using the VTK convention.
  * \param[in] val_nPoly    - Polynomial degree of the element.
  * \param[in] val_constJac - Whether or not the Jacobians are constant.
  * \param[in] config       - Object, which contains the input parameters.
  * \return Local (to the element) index of the nodes that compose the face.
  */
  FEMStandardElementClass(unsigned short val_VTK_Type,
                          unsigned short val_nPoly,
                          bool           val_constJac,
                          CConfig        *config); 
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
  * \brief Static function, which makes available the number of DOFs for an element
           corresponding to the arguments.
  * \param[in] VTK_Type         - Type of the element using the VTK convention.
  * \param[in] nPoly            - Polynomial degree of the element.
  * \param[in] typeErrorMessage - Default argument used to write a good error message.
  * \return  The number of DOFs
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
  * \brief Function, which checks if the functions arguments corresponds to this standard element.
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
  * \param[in]  r   - Parametric coordinates for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient of the Vandermonde matrix in all r-locations.
  */
  void GradVandermonde1D(vector<su2double> &r,
                         vector<su2double> &VDr);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard triangle.
  * \param[in]  r   - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s   - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r- and s-locations.
  * \param[out] VDs - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r- and s-locations.
  */
  void GradVandermonde2D_Triangle(vector<su2double> &r,
                                  vector<su2double> &s,
                                  vector<su2double> &VDr,
                                  vector<su2double> &VDs);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard quadrilateral.
  * \param[in]  r   - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s   - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r- and s-locations.
  * \param[out] VDs - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r- and s-locations.
  */
  void GradVandermonde2D_Quadrilateral(vector<su2double> &r,
                                       vector<su2double> &s,
                                       vector<su2double> &VDr,
                                       vector<su2double> &VDs);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard tetrahedron.
  * \param[in]  r   - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s   - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t   - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Tetrahedron(vector<su2double> &r,
                                     vector<su2double> &s,
                                     vector<su2double> &t,
                                     vector<su2double> &VDr,
                                     vector<su2double> &VDs,
                                     vector<su2double> &VDt);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard pyramid.
  * \param[in]  r   - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s   - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t   - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Pyramid(vector<su2double> &r,
                                 vector<su2double> &s,
                                 vector<su2double> &t,
                                 vector<su2double> &VDr,
                                 vector<su2double> &VDs,
                                 vector<su2double> &VDt);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard prism.
  * \param[in]  r   - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s   - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t   - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Prism(vector<su2double> &r,
                               vector<su2double> &s,
                               vector<su2double> &t,
                               vector<su2double> &VDr,
                               vector<su2double> &VDs,
                               vector<su2double> &VDt);

  /*!
  * \brief Function, which computes the gradients of the Vandermonde matrix for a standard hexahedron.
  * \param[in]  r   - Parametric coordinate in r-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  s   - Parametric coordinate in s-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[in]  t   - Parametric coordinate in t-direction for which the gradient of the Vandermonde matrix must be computed.
  * \param[out] VDr - Matrix to store the gradient in r-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDs - Matrix to store the gradient in s-direction of the Vandermonde matrix in all r-, s- and t-locations.
  * \param[out] VDt - Matrix to store the gradient in t-direction of the Vandermonde matrix in all r-, s- and t-locations.
  */
  void GradVandermonde3D_Hexahedron(vector<su2double> &r,
                                    vector<su2double> &s,
                                    vector<su2double> &t,
                                    vector<su2double> &VDr,
                                    vector<su2double> &VDs,
                                    vector<su2double> &VDt);

  /*!
  * \brief Function, which determines the integration points for a triangle such that polynomials of
           orderExact are integrated exactly.
  */
  void IntegrationPointsTriangle(void);

  /*!
  * \brief Function, which determines the integration points for a tetrahedron such that polynomials of
           orderExact are integrated exactly.
  */
  void IntegrationPointsTetrahedron(void);

  /*!
  * \brief Function, which determines the integration points for a pyramid such that polynomials of
           orderExact are integrated exactly.
  */
  void IntegrationPointsPyramid(void);

  /*!
  * \brief Function, which computes the inverse of the given square matrix.
  * \param[in]     n - Number of rows/columns of the square matrix A.
  * \param[in,out] A - On input the square matrix to be inverted. On output the inverse.
  */
  void InverseMatrix(unsigned short    n,
                     vector<su2double> &A);

  /*!
  * \brief Function, which computes the value of the Legendre polynomials Pn and Pnm1 for given x (-1 <= x <= 1) and n.
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
  * \brief Function, which carries out a matrix matrix multiplication to obtain data in the integration
           points and stores the transpose of the result.
  * \param[in]  A - First matrix in the matrix matrix product A*B, dimension nIntegration X nDOFs.
  * \param[in]  B - Second matrix in the matrix matrix product A*B, dimension nDOFs X nDOFs.
  * \param[out] C - Result of A*B. The transpose of the result is stored, dimension nDOFs X nIntegration.
  */
  void MatMulTranspose(vector<su2double> &A,
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

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard 1D edge.
   * \param[in]  r - Parametric coordinates for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r-locations.
  */
  void Vandermonde1D(vector<su2double> &r,
                     vector<su2double> &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard triangle.
   * \param[in]  r - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r- and s-locations.
  */
  void Vandermonde2D_Triangle(vector<su2double> &r,
                              vector<su2double> &s,
                              vector<su2double> &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard quadrilateral.
   * \param[in]  r - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r- and s-locations.
  */
  void Vandermonde2D_Quadrilateral(vector<su2double> &r,
                                   vector<su2double> &s,
                                   vector<su2double> &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  r - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Tetrahedron(vector<su2double> &r,
                                 vector<su2double> &s,
                                 vector<su2double> &t,
                                 vector<su2double> &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard pyramid.
   * \param[in]  r - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Pyramid(vector<su2double> &r,
                             vector<su2double> &s,
                             vector<su2double> &t,
                             vector<su2double> &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard prism.
   * \param[in]  r - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Prism(vector<su2double> &r,
                           vector<su2double> &s,
                           vector<su2double> &t,
                           vector<su2double> &V);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard hexahedron.
   * \param[in]  r - Parametric coordinates in r-direction for which the Vandermonde matrix must be computed.
   * \param[in]  s - Parametric coordinates in s-direction for which the Vandermonde matrix must be computed.
   * \param[in]  t - Parametric coordinates in t-direction for which the Vandermonde matrix must be computed.
   * \param[out] V - Matrix to store the Vandermonde matrix in all r-, s- and t-locations.
  */
  void Vandermonde3D_Hexahedron(vector<su2double> &r,
                                vector<su2double> &s,
                                vector<su2double> &t,
                                vector<su2double> &V);
};

#include "fem_standard_element.inl"
