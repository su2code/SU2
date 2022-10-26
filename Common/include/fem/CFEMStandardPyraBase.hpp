/*!
 * \file CFEMStandardPyraBase.hpp
 * \brief Base class for the FEM pyramid standard element.
 *        The functions are in the <i>CFEMStandardPyraBase.cpp</i> file.
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

#include "CFEMStandardQuadBase.hpp"
#include "CFEMStandardTriBase.hpp"

/*!
 * \class CFEMStandardPyraBase
 * \brief Base class which defines the variables and methods for the
 *        pyramid standard element.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CFEMStandardPyraBase: public virtual CFEMStandardQuadBase,
                            public virtual CFEMStandardTriBase  {
public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CFEMStandardPyraBase() = default;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardPyraBase(const unsigned short val_nPoly,
                       const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardPyraBase() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which returns the number of faces of the volume element.
   * \return The number of faces of the volume element.
   */
  inline unsigned short GetNFaces(void) const override {return 5;}

  /*!
   * \brief Function, which returns the VTK type of the given face index.
   * \param[in] ind - Index of the face for which the VTK type must be returned.
   * \return The VTK type of the given face id of the element.
   */
  inline unsigned short GetVTK_Face(unsigned short ind) const override {
    if(ind == 0) return QUADRILATERAL;
    else         return TRIANGLE;
  }

protected:

  vector<passivedouble> rLineIntGL;      /*!< \brief 1D parametric coordinates of the
                                                     Gauss-Legendre integration points. */
  vector<passivedouble> wLineIntGL;      /*!< \brief Weights of the 1D Gauss-Legendre
                                                     integration points. */

  vector<passivedouble> rLineIntGJ;      /*!< \brief 1D parametric coordinates of the
                                                     Gauss-Jacobi integration points. */
  vector<passivedouble> wLineIntGJ;      /*!< \brief Weights of the 1D Gauss-Jacobi
                                                     integration points. */

  /*-----------------------------------------------------------------------------------*/
  /*---                         Protected member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which converts the 2D parametric coordinates of a quadrilateral face
   *        of the pyramid to the 3D parametric coordinates of the actual pyramid.
   * \param[in]  rLine       - Parametric coordinates of the 1D reference line element.
   * \param[in]  faceID_Elem - The face ID of the element on which the face resides.
   * \param[in]  orientation - Orientation of the face w.r.t. the prism element.
   * \param[out] rPyra       - Parametric r-coordinates of the face points in the pyramid.
   * \param[out] sPyra       - Parametric s-coordinates of the face points in the pyramid.
   * \param[out] tPyra       - Parametric t-coordinates of the face points in the pyramid.
   */
  void ConvertCoor2DQuadFaceTo3DPyra(const vector<passivedouble> &rLine,
                                     const unsigned short        faceID_Elem,
                                     const unsigned short        orientation,
                                     vector<passivedouble>       &rPyra,
                                     vector<passivedouble>       &sPyra,
                                     vector<passivedouble>       &tPyra);

  /*!
   * \brief Function, which converts the 2D parametric coordinates of a triangular face
   *        of the pyramid to the 3D parametric coordinates of the actual pyramid.
   * \param[in]  rF          - Parametric r-coordinates of the triangular face.
   * \param[in]  sF          - Parametric s-coordinates of the triangular face. 
   * \param[in]  faceID_Elem - The face ID of the element on which the face resides.
   * \param[in]  orientation - Orientation of the face w.r.t. the pyramid element.
   * \param[out] rPyra       - Parametric r-coordinates of the face points in the pyramid.
   * \param[out] sPyra       - Parametric s-coordinates of the face points in the pyramid.
   * \param[out] tPyra       - Parametric t-coordinates of the face points in the pyramid.
   *                           for the face in the actual prism.
   */
  void ConvertCoor2DTriFaceTo3DPyra(const vector<passivedouble> &rF,
                                    const vector<passivedouble> &sF,
                                    const unsigned short        faceID_Elem,
                                    const unsigned short        orientation,
                                    vector<passivedouble>       &rPyra,
                                    vector<passivedouble>       &sPyra,
                                    vector<passivedouble>       &tPyra);

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a pyramid in the integration points for the given
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
  void DerLagBasisIntPointsPyra(const unsigned short                   mPoly,
                                const vector<passivedouble>            &rDOFs,
                                const vector<passivedouble>            &sDOFs,
                                const vector<passivedouble>            &tDOFs,
                                const vector<passivedouble>            &rInt,
                                const vector<passivedouble>            &sInt,
                                const vector<passivedouble>            &tInt,
                                vector<ColMajorMatrix<passivedouble> > &derLag);

  /*!
   * \brief Function, which computes the values of the 2nd derivatives of the Lagrangian
   *        basis functions of a pyramid in the integration points for the given
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
  void HesLagBasisIntPointsPyra(const unsigned short                   mPoly,
                                const vector<passivedouble>            &rDOFs,
                                const vector<passivedouble>            &sDOFs,
                                const vector<passivedouble>            &tDOFs,
                                const vector<passivedouble>            &rInt,
                                const vector<passivedouble>            &sInt,
                                const vector<passivedouble>            &tInt,
                                vector<ColMajorMatrix<passivedouble> > &hesLag);

  /*!
   * \brief Function, which computes the values of the Lagrangian basis functions of
   *        a pyramid in the integration points for the given location of the DOFs.
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
  void LagBasisIntPointsPyra(const unsigned short          mPoly,
                             const vector<passivedouble>   &rDOFs,
                             const vector<passivedouble>   &sDOFs,
                             const vector<passivedouble>   &tDOFs,
                             const vector<passivedouble>   &rInt,
                             const vector<passivedouble>   &sInt,
                             const vector<passivedouble>   &tInt,
                             ColMajorMatrix<passivedouble> &lag);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard pyramid.
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
  void GradVandermondePyramid(const unsigned short          mPoly,
                              const vector<passivedouble>   &r,
                              const vector<passivedouble>   &s,
                              const vector<passivedouble>   &t,
                              ColMajorMatrix<passivedouble> &VDr,
                              ColMajorMatrix<passivedouble> &VDs,
                              ColMajorMatrix<passivedouble> &VDt);

  /*!
   * \brief Function, which computes the Hessian of the Vandermonde matrix for a standard pyramid.
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
  void HesVandermondePyramid(const unsigned short          mPoly,
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
   * \brief Function, which computes the Vandermonde matrix for a standard pyramid.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t     - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r,s,t-locations.
   */
  void VandermondePyramid(const unsigned short          mPoly,
                          const vector<passivedouble>   &r,
                          const vector<passivedouble>   &s,
                          const vector<passivedouble>   &t,
                          ColMajorMatrix<passivedouble> &V);

  /*!
   * \brief Function, which computes the parametric coordinates of all
   *        the integration points.
   * \param[out] rInt - Parametric r-coordinates of the integration points.
   * \param[out] sInt - Parametric s-coordinates of the integration points.
   * \param[out] tInt - Parametric t-coordinates of the integration points.
   */
  void LocationAllIntegrationPoints(vector<passivedouble> &rInt,
                                    vector<passivedouble> &sInt,
                                    vector<passivedouble> &tInt);
  /*!
   * \brief Function, which determines the location of the grid DOFs of a pyramid for
   *        polynomial degree nPoly when an equidistant spacing is used.
   * \param[out] rDOFs - Parametric r-coordinates of the nodal DOFs.
   * \param[out] sDOFs - Parametric s-coordinates of the nodal DOFs.
   * \param[out] tDOFs - Parametric t-coordinates of the nodal DOFs.
   */
  void LocationPyramidGridDOFsEquidistant(const unsigned short  mPoly,
                                          vector<passivedouble> &rDOFs,
                                          vector<passivedouble> &sDOFs,
                                          vector<passivedouble> &tDOFs);

  /*!
   * \brief Function, which determines the location of the grid DOFs of a pyramid for
   *        polynomial degree nPoly when the LGL distribution is used.
   * \param[out] rDOFs - Parametric r-coordinates of the nodal DOFs.
   * \param[out] sDOFs - Parametric s-coordinates of the nodal DOFs.
   * \param[out] tDOFs - Parametric t-coordinates of the nodal DOFs.
   */
  void LocationPyramidGridDOFsLGL(const unsigned short  mPoly,
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
};
