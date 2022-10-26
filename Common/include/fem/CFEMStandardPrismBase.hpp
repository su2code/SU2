/*!
 * \file CFEMStandardPrismBase.hpp
 * \brief Base class for the FEM prism standard element.
 *        The functions are in the <i>CFEMStandardPrismBase.cpp</i> file.
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
 * \class CFEMStandardPrismBase
 * \brief Base class which defines the variables and methods for the
 *        prism standard element.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CFEMStandardPrismBase: public virtual CFEMStandardQuadBase,
                             public virtual CFEMStandardTriBase  {
public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CFEMStandardPrismBase() = default;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardPrismBase(const unsigned short val_nPoly,
                        const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardPrismBase() = default;

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
    if(ind < 2) return TRIANGLE;
    else        return QUADRILATERAL;
  }

protected:

  unsigned short nDOFsTriangle;  /*!< \brief Number of DOFs of the triangular parts. */
  unsigned short nIntTriangle;   /*!< \brief Number of integration points of the triangular parts. */

  /*-----------------------------------------------------------------------------------*/
  /*---                         Protected member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which converts the 2D parametric coordinates of a quadrilateral face
   *        of the prism to the 3D parametric coordinates of the actual prism.
   * \param[in]  rLine       - Parametric coordinates of the 1D reference line element.
   * \param[in]  faceID_Elem - The face ID of the element on which the face resides.
   * \param[in]  orientation - Orientation of the face w.r.t. the prism element.
   * \param[out] rPrism      - Parametric r-coordinates of the face points in the prism.
   * \param[out] sPrism      - Parametric s-coordinates of the face points in the prism.
   * \param[out] tPrism      - Parametric t-coordinates of the face points in the prism.
   */
  void ConvertCoor2DQuadFaceTo3DPrism(const vector<passivedouble> &rLine,
                                      const unsigned short        faceID_Elem,
                                      const unsigned short        orientation,
                                      vector<passivedouble>       &rPrism,
                                      vector<passivedouble>       &sPrism,
                                      vector<passivedouble>       &tPrism);

  /*!
   * \brief Function, which converts the 2D parametric coordinates of a triangular face
   *        of the prism to the 3D parametric coordinates of the actual prism.
   * \param[in]  rF             - Parametric r-coordinates of the triangular face.
   * \param[in]  sF             - Parametric s-coordinates of the triangular face. 
   * \param[in]  faceID_Elem    - The face ID of the element on which the face resides.
   * \param[in]  orientation    - Orientation of the face w.r.t. the prism element.
   * \param[out] rTrianglePrism - Parametric r-coordinates of the face points in the prism.
   * \param[out] sTrianglePrism - Parametric s-coordinates of the face points in the prism.
   * \param[out] rLinePrism     - Parametric coordinate in the structured direction
   *                              for the face in the actual prism.
   */
  void ConvertCoor2DTriFaceTo3DPrism(const vector<passivedouble> &rF,
                                     const vector<passivedouble> &sF,
                                     const unsigned short        faceID_Elem,
                                     const unsigned short        orientation,
                                     vector<passivedouble>       &rTrianglePrism,
                                     vector<passivedouble>       &sTrianglePrism,
                                     vector<passivedouble>       &rLinePrism);

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a prism in the integration points for the given
   *        location of the DOFs.
   * \param[in]  mPoly         - Polynomial degree of the basis functions.
   * \param[in]  rTriangleDOFs - Vector, which contains the parametric r-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  sTriangleDOFs - Vector, which contains the parametric s-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  rLineDOFs     - Vector, which contains the parametric r-locations of the DOFs
   *                             of the line normal to the triangles.
   * \param[in]  rTriangleInts - Vector, which contains the parametric r-locations of the
   *                             integration points of the base triangle.
   * \param[in]  sTriangleInts - Vector, which contains the parametric s-locations of the
   *                             integration points of the base triangle.
   * \param[in]  rLineInts     - Vector, which contains the parametric r-locations of the
   *                             integration points of the line normal to the triangles.
   * \param[out] derLag        - Matrix, which contains the values of derivatives of all the
   *                             Lagrangian basis functions in all the integration points.
   */
  void DerLagBasisIntPointsPrism(const unsigned short                   mPoly,
                                 const vector<passivedouble>            &rTriangleDOFs,
                                 const vector<passivedouble>            &sTriangleDOFs,
                                 const vector<passivedouble>            &rLineDOFs,
                                 const vector<passivedouble>            &rTriangleInts,
                                 const vector<passivedouble>            &sTrianglInts,
                                 const vector<passivedouble>            &rLineInts,
                                 vector<ColMajorMatrix<passivedouble> > &derLag);

  /*!
   * \brief Function, which computes the values of the 2nd derivatives of the Lagrangian
   *        basis functions of a prism in the integration points for the given
   *        location of the DOFs.
   * \param[in]  mPoly         - Polynomial degree of the basis functions.
   * \param[in]  rTriangleDOFs - Vector, which contains the parametric r-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  sTriangleDOFs - Vector, which contains the parametric s-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  rLineDOFs     - Vector, which contains the parametric r-locations of the DOFs
   *                             of the line normal to the triangles.
   * \param[in]  rTriangleInts - Vector, which contains the parametric r-locations of the
   *                             integration points of the base triangle.
   * \param[in]  sTriangleInts - Vector, which contains the parametric s-locations of the
   *                             integration points of the base triangle.
   * \param[in]  rLineInts     - Vector, which contains the parametric r-locations of the
   *                             integration points of the line normal to the triangles.
   * \param[out] hesLag        - Matrix, which contains the values of 2nd derivatives of all the
   *                             Lagrangian basis functions in all the integration points.
   */
  void HesLagBasisIntPointsPrism(const unsigned short                   mPoly,
                                 const vector<passivedouble>            &rTriangleDOFs,
                                 const vector<passivedouble>            &sTriangleDOFs,
                                 const vector<passivedouble>            &rLineDOFs,
                                 const vector<passivedouble>            &rTriangleInts,
                                 const vector<passivedouble>            &sTrianglInts,
                                 const vector<passivedouble>            &rLineInts,
                                 vector<ColMajorMatrix<passivedouble> > &hesLag);

  /*!
   * \brief Function, which computes the values of the Lagrangian basis functions of
   *        a prism in the integration points for the given location of the DOFs.
   * \param[in]  mPoly         - Polynomial degree of the basis functions.
   * \param[in]  rTriangleDOFs - Vector, which contains the parametric r-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  sTriangleDOFs - Vector, which contains the parametric s-locations of the DOFs
   *                             of the base triangle.
   * \param[in]  rLineDOFs     - Vector, which contains the parametric r-locations of the DOFs
   *                             of the line normal to the triangles.
   * \param[in]  rTriangleInts - Vector, which contains the parametric r-locations of the
   *                             integration points of the base triangle.
   * \param[in]  sTriangleInts - Vector, which contains the parametric s-locations of the
   *                             integration points of the base triangle.
   * \param[in]  rLineInts     - Vector, which contains the parametric r-locations of the
   *                             integration points of the line normal to the triangles.
   * \param[out] lag           - Matrix, which contains the values of all the Lagrangian
   *                             basis functions in all the integration points.
   */
  void LagBasisIntPointsPrism(const unsigned short          mPoly,
                              const vector<passivedouble>   &rTriangleDOFs,
                              const vector<passivedouble>   &sTriangleDOFs,
                              const vector<passivedouble>   &rLineDOFs,
                              const vector<passivedouble>   &rTriangleInts,
                              const vector<passivedouble>   &sTrianglInts,
                              const vector<passivedouble>   &rLineInts,
                              ColMajorMatrix<passivedouble> &lag);

  /*!
   * \brief Function, which computes the gradients of the Vandermonde matrix for a standard prism.
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
  void GradVandermondePrism(const unsigned short          mPoly,
                            const vector<passivedouble>   &r,
                            const vector<passivedouble>   &s,
                            const vector<passivedouble>   &t,
                            ColMajorMatrix<passivedouble> &VDr,
                            ColMajorMatrix<passivedouble> &VDs,
                            ColMajorMatrix<passivedouble> &VDt);

  /*!
   * \brief Function, which computes the Hessian of the Vandermonde matrix for a standard prism.
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
  void HesVandermondePrism(const unsigned short          mPoly,
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
   * \brief Function, which computes the Vandermonde matrix for a standard prism.
   * \param[in]  mPoly - Polynomial degree of the basis functions.
   * \param[in]  r     - Parametric r-coordinates for which the Vandermonde matrix must be computed.
   * \param[in]  s     - Parametric s-coordinates for which the Vandermonde matrix must be computed
   * \param[in]  t     - Parametric t-coordinates for which the Vandermonde matrix must be computed
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r,s,t-locations.
   */
  void VandermondePrism(const unsigned short          mPoly,
                        const vector<passivedouble>   &r,
                        const vector<passivedouble>   &s,
                        const vector<passivedouble>   &t,
                        ColMajorMatrix<passivedouble> &V);

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements.
   */
  void SubConnLinearElements(void);

  /*!
   * \brief Function, which computes the parametric coordinates of all the relevant
   *        points of a prism from the triangle and line data.
   * \param[in]  rTriangle - Vector, which contains the parametric r-locations of the base triangle.
   * \param[in]  sTriangle - Vector, which contains the parametric s-locations of the base triangle.
   * \param[in]  rLine     - Vector, which contains the parametric r-locations of the normal line.
   * \param[out] rPrism    - Parametric r-coordinates of the entire prism.
   * \param[out] sPrism    - Parametric s-coordinates of the entire prism.
   * \param[out] tPrism    - Parametric t-coordinates of the entire prism.
   */
  void LocationAllPointsPrism(const vector<passivedouble> &rTriangle,
                              const vector<passivedouble> &sTriangle,
                              const vector<passivedouble> &rLine,
                              vector<passivedouble>       &rPrism,
                              vector<passivedouble>       &sPrism,
                              vector<passivedouble>       &tPrism);

  /*!
   * \brief Function, which creates the local grid connectivities of the faces
   *        of the volume element.
   */
  void LocalGridConnFaces(void);
};
