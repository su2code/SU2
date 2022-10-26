/*!
 * \file CFEMStandardTriBase.hpp
 * \brief Base class for the FEM triangle standard element.
 *        The functions are in the <i>CFEMStandardTriBase.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#pragma once 

#include "CFEMStandardLineBase.hpp"

/*!
 * \class CFEMStandardTriBase
 * \brief Base class which defines methods and variables for
 *        the triangle standard element.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CFEMStandardTriBase: public virtual CFEMStandardLineBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CFEMStandardTriBase() = default;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardTriBase(const unsigned short val_nPoly,
                      const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardTriBase() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public static member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Static function, which makes available the number of integration
   *        points for a triangle.
   * \param[in] orderExact - Polynomial degree that must be integrated exactly.
   * \return The number of integration points for a triangle.
   */
  static unsigned short GetNIntTriangleStatic(unsigned short orderExact);

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which returns the number of faces of the volume element.
   * \return The number of faces of the volume element.
   */
  inline unsigned short GetNFaces(void) const override {return 3;}

  /*!
   * \brief Function, which returns the VTK type of the given face index.
   * \param[in] ind - Index of the face for which the VTK type must be returned.
   * \return The VTK type of the given face id of the element.
   */
  inline unsigned short GetVTK_Face(unsigned short ind) const override {return LINE;}

protected:

  vector<passivedouble> rTriangleInt;      /*!< \brief Parametric r-coordinates of the integration
                                                       points of the triangle. */
  vector<passivedouble> sTriangleInt;      /*!< \brief Parametric s-coordinates of the integration
                                                       points of the triangle. */
  vector<passivedouble> wTriangleInt;      /*!< \brief Weights of the integration points of the
                                                       triangle. */

  /*-----------------------------------------------------------------------------------*/
  /*---                         Protected member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which changes the given triangular connectivity, such that
   *        the direction coincides with the direction corresponding to corner
   *        vertices vert0, vert1, vert2.
   * \param[in,out] connTriangle - Connectivity of the triangle, that must be adapted.
   * \param[in]     vert0        - Corner vertex 0 of the desired sequence.
   * \param[in]     vert1        - Corner vertex 1 of the desired sequence.
   * \param[in]     vert2        - Corner vertex 2 of the desired sequence.
   */
  void ChangeDirectionTriangleConn(vector<unsigned short> &connTriangle,
                                   const unsigned short   vert0,
                                   const unsigned short   vert1,
                                   const unsigned short   vert2);

  /*!
   * \brief Function, which converts the 1D parametric coordinates of a face of
   *        the triangle to the 2D parametric coordinates of the actual triangle.
   * \param[in]  rLine       - 1D parametric coordinates of the face of the triangle.
   * \param[in]  faceID      - The corresponding faceID of the adjacent triangle.
   * \param[in]  orientation - Orientation of the line element relative to the triangle.
   * \param[out] rTriangle   - The parametric r-coordinates of the triangle
   *                           corresponding to rLine.
   * \param[out] sTriangle   - The parametric s-coordinates of the triangle
   *                           corresponding to rLine.
   */
  void ConvertCoor1DFaceTo2DTriangle(const vector<passivedouble> rLine,
                                     const unsigned short        faceID,
                                     const unsigned short        orientation,
                                     vector<passivedouble>       &rTriangle,
                                     vector<passivedouble>       &sTriangle);

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a triangle in the integration points for the given
   *        location of the DOFs.
   * \param[in]  mPoly  - Polynomial degree of the basis functions.
   * \param[in]  rDOFs  - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs  - Vector, which contains the parametric s-locations of the DOFs.
   * \param[in]  rInt   - Vector, which contains the parametric r-locations of the integration points.
   * \param[in]  sInt   - Vector, which contains the parametric s-locations of the integration points.
   * \param[out] derLag - Matrix, which contains the values of derivatives of all the
   *                      Lagrangian basis functions in all the integration points.
   */
  void DerLagBasisIntPointsTriangle(const unsigned short                   mPoly,
                                    const vector<passivedouble>            &rDOFs,
                                    const vector<passivedouble>            &sDOFs,
                                    const vector<passivedouble>            &rInt,
                                    const vector<passivedouble>            &sInt,
                                    vector<ColMajorMatrix<passivedouble> > &derLag);

  /*!
   * \brief Function, which computes the values of the 2nd derivatives of the Lagrangian
   *        basis functions of a triangle in the integration points for the given
   *        location of the DOFs.
   * \param[in]  mPoly  - Polynomial degree of the basis functions.
   * \param[in]  rDOFs  - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs  - Vector, which contains the parametric s-locations of the DOFs.
   * \param[in]  rInt   - Vector, which contains the parametric r-locations of the integration points.
   * \param[in]  sInt   - Vector, which contains the parametric s-locations of the integration points.
   * \param[out] hesLag - Matrix, which contains the values of 2nd derivatives of all the
   *                      Lagrangian basis functions in all the integration points.
   */
  void HesLagBasisIntPointsTriangle(const unsigned short                   mPoly,
                                    const vector<passivedouble>            &rDOFs,
                                    const vector<passivedouble>            &sDOFs,
                                    const vector<passivedouble>            &rInt,
                                    const vector<passivedouble>            &sInt,
                                    vector<ColMajorMatrix<passivedouble> > &hesLag);

  /*!
   * \brief Function, which computes the values of the Lagrangian basis functions of
   *        a triangle in the integration points for the given location of the DOFs.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  rDOFs - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs - Vector, which contains the parametric s-locations of the DOFs.
   * \param[in]  rInt  - Vector, which contains the parametric r-locations of the integration points.
   * \param[in]  sInt  - Vector, which contains the parametric s-locations of the integration points.
   * \param[out] lag   - Matrix, which contains the values of all the Lagrangian
   *                     basis functions in all the integration points.
   */
  void LagBasisIntPointsTriangle(const unsigned short          mPoly,
                                 const vector<passivedouble>   &rDOFs,
                                 const vector<passivedouble>   &sDOFs,
                                 const vector<passivedouble>   &rInt,
                                 const vector<passivedouble>   &sInt,
                                 ColMajorMatrix<passivedouble> &lag);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard triangle.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[out] VDr   - Matrix to store the derivative in r-direction of the Vandermonde matrix
   *                     in all r,s-locations.
   * \param[out] VDs   - Matrix to store the derivative in s-direction of the Vandermonde matrix
   *                     in all r,s-locations.
   */
  void GradVandermondeTriangle(const unsigned short          mPoly,
                               const vector<passivedouble>   &r,
                               const vector<passivedouble>   &s,
                               ColMajorMatrix<passivedouble> &VDr,
                               ColMajorMatrix<passivedouble> &VDs);

  /*!
   * \brief Function, which computes the Hessian of the Vandermonde matrix for a standard triangle.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[out] VDr2  - Matrix to store the 2nd derivative in r-direction of the Vandermonde matrix
   *                     in all r,s-locations.
   * \param[out] VDs2  - Matrix to store the 2nd derivative in s-direction of the Vandermonde matrix
   *                     in all r,s-locations.
   * \param[out] VDrs  - Matrix to store the cross derivative in r,s-direction of the Vandermonde matrix
   *                     in all r,s-locations.
   */
  void HesVandermondeTriangle(const unsigned short          mPoly,
                              const vector<passivedouble>   &r,
                              const vector<passivedouble>   &s,
                              ColMajorMatrix<passivedouble> &VDr2,
                              ColMajorMatrix<passivedouble> &VDs2,
                              ColMajorMatrix<passivedouble> &VDrs);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard triangle.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r,s-locations.
   */
  void VandermondeTriangle(const unsigned short          mPoly,
                           const vector<passivedouble>   &r,
                           const vector<passivedouble>   &s,
                           ColMajorMatrix<passivedouble> &V);

  /*!
   * \brief Function, which determines the integration points for a triangle
   *        such that polynomials of orderExact are integrated exactly.
   */
  void IntegrationPointsTriangle(void);

  /*!
   * \brief Function, which determines the location of the triangular grid DOFs for
   *        polynomial degree nPoly when an equidistant spacing is used.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[out] r     - Vector of the parametric r-coordinates of the DOFs.
   * \param[out] s     - Vector of the parametric s-coordinates of the DOFs.
   */
  void LocationTriangleGridDOFsEquidistant(const unsigned short  mPoly,
                                           vector<passivedouble> &r,
                                           vector<passivedouble> &s);

  /*!
   * \brief Function, which determines the location of the triangular grid DOFs for
   *        polynomial degree nPoly when the LGL distribution is used. The definition
   *        used is according to the book of Hesthaven and Warburton.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[out] r     - Vector of the parametric r-coordinates of the DOFs.
   * \param[out] s     - Vector of the parametric s-coordinates of the DOFs.
   */
  void LocationTriangleGridDOFsLGL(const unsigned short  mPoly,
                                   vector<passivedouble> &r,
                                   vector<passivedouble> &s);

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements.
   */
  void SubConnLinearElements(void);

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements. The splitting is done for one face
   *        w.r.t the volume. The output is stored in subConn1ForPlotting using node ID
   *        available in gridConnFaces.
   */
  void SubConnLinearElementsFace(int val_faceID_Elem);

    /*!
   * \brief Function, which creates the local grid connectivities of the faces
   *        of the volume element.
   */
  void LocalGridConnFaces(void);
  
private:
  /*-----------------------------------------------------------------------------------*/
  /*---                         Private member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes a scaled warp function for all DOFs of a standard triangle
   *        with a polynomial basis function of order nPoly. This info is used to construct
   *        the location of the DOFs inside a standard triangle.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[in]  rout  - Input coordinates for the warp factor.
   * \param[out] warp  - Warp factor to be computed.
   */
  void WarpFactor(const unsigned short        mPoly,
                  const vector<passivedouble> &rout,
                  vector<passivedouble>       &warp);
};
