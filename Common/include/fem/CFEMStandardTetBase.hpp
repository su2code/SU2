/*!
 * \file CFEMStandardTetBase.hpp
 * \brief Base class for the FEM tetrahedron standard element.
 *        The functions are in the <i>CFEMStandardTet.cpp</i> file.
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

#include "CFEMStandardTriBase.hpp"

/*!
 * \class CFEMStandardTetBase
 * \brief Base class which defines the variables and methods for the
 *        tetrahedron standard element.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CFEMStandardTetBase: public virtual CFEMStandardTriBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CFEMStandardTetBase() = default;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardTetBase(const unsigned short val_nPoly,
                      const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardTetBase() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public static member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Static function, which makes available the number of integration
   *        points for a tetrahedron.
   * \param[in] orderExact - Polynomial degree that must be integrated exactly.
   * \return The number of integration points for a tetrahedron.
   */
  static unsigned short GetNIntTetrahedronStatic(unsigned short orderExact);

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which returns the number of faces of the volume element.
   * \return The number of faces of the volume element.
   */
  inline unsigned short GetNFaces(void) const override {return 4;}

  /*!
   * \brief Function, which returns the VTK type of the given face index.
   * \param[in] ind - Index of the face for which the VTK type must be returned.
   * \return The VTK type of the given face id of the element.
   */
  inline unsigned short GetVTK_Face(unsigned short ind) const override {return TRIANGLE;}

protected:

  vector<passivedouble> rTetInt;      /*!< \brief Parametric r-coordinates of the integration
                                                  points of the tetrahedron. */
  vector<passivedouble> sTetInt;      /*!< \brief Parametric s-coordinates of the integration
                                                  points of the tetrahedron. */
  vector<passivedouble> tTetInt;      /*!< \brief Parametric t-coordinates of the integration
                                                  points of the tetrahedron. */
  vector<passivedouble> wTetInt;      /*!< \brief Weights of the integration points of the
                                                  tetrahedron. */

  /*-----------------------------------------------------------------------------------*/
  /*---                         Protected member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which converts the 2D parametric coordinates of a triangular face
   *        of the tetrahedron to the 3D parametric coordinates of the actual tetrahedron.
   * \param[in]  rF          - Parametric r-coordinates of the triangular face.
   * \param[in]  sF          - Parametric s-coordinates of the triangular face. 
   * \param[in]  faceID_Elem - The face ID of the element on which the face resides.
   * \param[in]  orientation - Orientation of the face w.r.t. the tetrahedron element.
   * \param[out] rTet        - Parametric r-coordinates of the face points in the pyramid.
   * \param[out] sTet        - Parametric s-coordinates of the face points in the pyramid.
   * \param[out] tTet        - Parametric t-coordinates of the face points in the pyramid.
   *                           for the face in the actual prism.
   */
  void ConvertCoor2DTriFaceTo3DTet(const vector<passivedouble> &rF,
                                   const vector<passivedouble> &sF,
                                   const unsigned short        faceID_Elem,
                                   const unsigned short        orientation,
                                   vector<passivedouble>       &rTet,
                                   vector<passivedouble>       &sTet,
                                   vector<passivedouble>       &tTet);

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a tetrahedron in the integration points for the given
   *        location of the DOFs.
   * \param[in]  mPoly  - Polynomial degree of the basis functions.
   * \param[in]  rDOFs  - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs  - Vector, which contains the parametric s-locations of the DOFs.
   * \param[in]  tDOFs  - Vector, which contains the parametric t-locations of the DOFs
   * \param[in]  rInt   - Vector, which contains the parametric r-locations of the integration points.
   * \param[in]  sInt   - Vector, which contains the parametric s-locations of the integration points.
   * \param[in]  tInt   - Vector, which contains the parametric t-locations of the integration points.
   * \param[out] derLag - Matrix, which contains the values of derivatives of all the
   *                      Lagrangian basis functions in all the integration points.
   */
  void DerLagBasisIntPointsTet(const unsigned short                   mPoly,
                               const vector<passivedouble>            &rDOFs,
                               const vector<passivedouble>            &sDOFs,
                               const vector<passivedouble>            &tDOFs,
                               const vector<passivedouble>            &rInt,
                               const vector<passivedouble>            &sInt,
                               const vector<passivedouble>            &tInt,
                               vector<ColMajorMatrix<passivedouble> > &derLag);

  /*!
   * \brief Function, which computes the values of the 2nd derivatives of the Lagrangian
   *        basis functions of a tetrahedron in the integration points for the given
   *        location of the DOFs.
   * \param[in]  mPoly  - Polynomial degree of the basis functions.
   * \param[in]  rDOFs  - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs  - Vector, which contains the parametric s-locations of the DOFs.
   * \param[in]  tDOFs  - Vector, which contains the parametric t-locations of the DOFs.
   * \param[in]  rInt   - Vector, which contains the parametric r-locations of the integration points.
   * \param[in]  sInt   - Vector, which contains the parametric s-locations of the integration points.
   * \param[in]  tInt   - Vector, which contains the parametric t-locations of the integration points.
   * \param[out] hesLag - Matrix, which contains the values of 2nd derivatives of all the
   *                      Lagrangian basis functions in all the integration points.
   */
  void HesLagBasisIntPointsTet(const unsigned short                   mPoly,
                               const vector<passivedouble>            &rDOFs,
                               const vector<passivedouble>            &sDOFs,
                               const vector<passivedouble>            &tDOFs,
                               const vector<passivedouble>            &rInt,
                               const vector<passivedouble>            &sInt,
                               const vector<passivedouble>            &tInt,
                               vector<ColMajorMatrix<passivedouble> > &hesLag); 

  /*!
   * \brief Function, which computes the values of the Lagrangian basis functions of
   *        a tetrahedron in the integration points for the given location of the DOFs.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  rDOFs - Vector, which contains the parametric r-locations of the DOFs.
   * \param[in]  sDOFs - Vector, which contains the parametric s-locations of the DOFs.
   * \param[in]  tDOFs - Vector, which contains the parametric t-locations of the DOFs.
   * \param[in]  rInt  - Vector, which contains the parametric r-locations of the integration points.
   * \param[in]  sInt  - Vector, which contains the parametric s-locations of the integration points.
   * \param[in]  tInt  - Vector, which contains the parametric t-locations of the integration points.
   * \param[out] lag   - Matrix, which contains the values of all the Lagrangian
   *                     basis functions in all the integration points.
   */
  void LagBasisIntPointsTet(const unsigned short          mPoly,
                            const vector<passivedouble>   &rDOFs,
                            const vector<passivedouble>   &sDOFs,
                            const vector<passivedouble>   &tDOFs,
                            const vector<passivedouble>   &rInt,
                            const vector<passivedouble>   &sInt,
                            const vector<passivedouble>   &tInt,
                            ColMajorMatrix<passivedouble> &lag);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t     - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] VDr   - Matrix to store the derivative in r-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDs   - Matrix to store the derivative in s-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDt   - Matrix to store the derivative in t-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   */
  void GradVandermondeTetrahedron(const unsigned short          mPoly,
                                  const vector<passivedouble>   &r,
                                  const vector<passivedouble>   &s,
                                  const vector<passivedouble>   &t,
                                  ColMajorMatrix<passivedouble> &VDr,
                                  ColMajorMatrix<passivedouble> &VDs,
                                  ColMajorMatrix<passivedouble> &VDt);

  /*!
   * \brief Function, which computes the Hessian of the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t     - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] VDr2  - Matrix to store the 2nd derivative in r-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDs2  - Matrix to store the 2nd derivative in s-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDt2  - Matrix to store the 2nd derivative in t-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDrs  - Matrix to store the cross derivative in rs-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDrt  - Matrix to store the cross derivative in rt-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   * \param[out] VDst  - Matrix to store the cross derivative in st-direction of the Vandermonde matrix
   *                     in all r,s,t-locations.
   */
  void HesVandermondeTetrahedron(const unsigned short          mPoly,
                                 const vector<passivedouble>   &r,
                                 const vector<passivedouble>   &s,
                                 const vector<passivedouble>   &t,
                                 ColMajorMatrix<passivedouble> &VDr2,
                                 ColMajorMatrix<passivedouble> &VDs2,
                                 ColMajorMatrix<passivedouble> &VDt2,
                                 ColMajorMatrix<passivedouble> &VDrs,
                                 ColMajorMatrix<passivedouble> &VDrt,
                                 ColMajorMatrix<passivedouble> &VDst);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard tetrahedron.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t     - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r,s,t-locations.
   */
  void VandermondeTetrahedron(const unsigned short          mPoly,
                              const vector<passivedouble>   &r,
                              const vector<passivedouble>   &s,
                              const vector<passivedouble>   &t,
                              ColMajorMatrix<passivedouble> &V);

  /*!
   * \brief Function, which determines the integration points for a tetrahedron
   *        such that polynomials of orderExact are integrated exactly.
   */
  void IntegrationPointsTetrahedron(void);

  /*!
   * \brief Function, which determines the location of the grid DOFs of a tetrahedron
   *        for polynomial degree nPoly when an equidistant spacing is used.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[out] rDOFs - Parametric r-coordinates of the DOFs.
   * \param[out] sDOFs - Parametric s-coordinates of the DOFs.
   * \param[out] tDOFs - Parametric t-coordinates of the DOFs.
   */
  void LocationTetGridDOFsEquidistant(const unsigned short  mPoly,
                                      vector<passivedouble> &rDOFs,
                                      vector<passivedouble> &sDOFs,
                                      vector<passivedouble> &tDOFs);

  /*!
   * \brief Function, which determines the location of the grid DOFs of a tetrahedron
   *        for polynomial degree nPoly when the LGL distribution is used.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[out] rDOFs - Parametric r-coordinates of the DOFs.
   * \param[out] sDOFs - Parametric s-coordinates of the DOFs.
   * \param[out] tDOFs - Parametric t-coordinates of the DOFs.
   */
  void LocationTetGridDOFsLGL(const unsigned short  mPoly,
                              vector<passivedouble> &rDOFs,
                              vector<passivedouble> &sDOFs,
                              vector<passivedouble> &tDOFs);

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements.
   */
  void SubConnLinearElements(void);

  /*!
   * \brief Function, which creates the local grid connectivities of the faces
   *        of the volume element.
   */
  void LocalGridConnFaces(void);

private:

  /*-----------------------------------------------------------------------------------*/
  /*---                          Private member functions.                          ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes the warping factors of a 3D triangular face,
            needed to compute the LGL distribution of a tetrahedron.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[in]  alpha - Blending factor for the given polynomial degree.
   * \param[in]  L1    - Barycentric coordinate of the DOFs (I think).
   * \param[in]  L2    - Barycentric coordinate of the DOFs (I think).
   * \param[in]  L3    - Barycentric coordinate of the DOFs (I think).
   * \param[out] dx    - Delta of the first parametric coordinate of the face.
   * \param[out] dy    - Delta of the second parametric coordinate of the face.
   */
  void EvalShift(const unsigned short        mPoly,
                 const passivedouble         alpha,
                 const vector<passivedouble> &L1,
                 const vector<passivedouble> &L2,
                 const vector<passivedouble> &L3,
                 vector<passivedouble>       &dx,
                 vector<passivedouble>       &dy);

  /*!
   * \brief Function, which computes the 1D edge warping.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[in]  xOut  - Coordinates to which the 1D edge warping must be applied.
   * \param[out] warp  - Edge warping value for all the nodal DOFs.
   */
  void EvalWarp(const unsigned short        mPoly,
                const vector<passivedouble> &xOut,
                vector<passivedouble>       &warp);
};
