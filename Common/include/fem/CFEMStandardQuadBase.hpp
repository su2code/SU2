/*!
 * \file CFEMStandardQuadBase.hpp
 * \brief Base class for the FEM quadrilateral standard element.
 *        The functions are in the <i>CFEMStandardQuadBase.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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
#include "../tensor_products/TensorProductVolumeIntPoints2D.hpp"

/*!
 * \class CFEMStandardQuadBase
 * \brief Base class which defines the variables and methods for the
 *        quadrilateral standard element.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardQuadBase: public virtual CFEMStandardLineBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CFEMStandardQuadBase() = default;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardQuadBase(const unsigned short val_nPoly,
                       const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardQuadBase() = default;

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
  inline unsigned short GetVTK_Face(unsigned short ind) const override {return LINE;}

protected:

  unsigned short nDOFs1D;   /*!< \brief Number of DOFs in one space direction. */
  unsigned short nInt1D;    /*!< \brief Number of integration points in one space direction. */

  /*-----------------------------------------------------------------------------------*/
  /*---                         Protected member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which changes the given quadrilateral connectivity, such that
   *        the direction coincides with the direction corresponding to corner
   *        vertices vert0, vert1, vert2, vert3.
   * \param[in,out] connQuad - Connectivity of the quadrilateral, that must be adapted.
   * \param[in]     vert0    - Corner vertex 0 of the desired sequence.
   * \param[in]     vert1    - Corner vertex 1 of the desired sequence.
   * \param[in]     vert2    - Corner vertex 2 of the desired sequence.
   * \param[in]     vert3    - Corner vertex 3 of the desired sequence.
   */
  void ChangeDirectionQuadConn(vector<unsigned short> &connQuad,
                               const unsigned short   vert0,
                               const unsigned short   vert1,
                               const unsigned short   vert2,
                               const unsigned short   vert3);

  /*!
   * \brief Function, which creates the components of the tensor to compute the data
   *        and derivatives on points on a line adjacent to a quad for the given
   *        face ID and orientation of the quadrilateral.
   * \param[in]  mPointsLine  - Number of points on the line.
   * \param[in]  mDOFs1D      - Number of DOFs in 1D of the quadrilateral.
   * \param[in]  faceID_Elem  - The face ID of the line in the adjacent quadrilateral.
   * \param[in]  orientation  - Orientation of the line w.r.t. the quadrilateral.
   * \param[in]  b1DPoints    - 1D basis functions of a line evaluated in the points.
   * \param[in]  derB1DPoints - Derivatives of bPoints1D.
   * \param[in]  b1DM1        - 1D basis functions of a line evaluated at r == -1.
   * \param[in]  b1DP1        - 1D basis functions of a line evaluated at r ==  1.
   * \param[in]  derB1DM1     - Derivatives of b1DM1.
   * \param[in]  derB1DP1     - Derivatives of b1DP1.
   * \param[out] tensorSol    - The two 1D components of the tensor to compute
   *                            the solution on the face of the quad.
   * \param[out] tensorDSolDr - The two 1D components of the tensor to compute the derivative
   *                            in r-direction of the solution on the face of the quad.
   * \param[out] tensorDSolDs - The two 1D components of the tensor to compute the derivative
   *                            in s-direction of the solution on the face of the quad.
   */
  void CreateTensorContributionsLineAdjQuad(const unsigned short                   mPointsLine,
                                            const unsigned short                   mDOFs1D,
                                            const unsigned short                   faceID_Elem,
                                            const unsigned short                   orientation,
                                            const ColMajorMatrix<passivedouble>    &b1DPoints,
                                            const ColMajorMatrix<passivedouble>    &derB1DPoints,
                                            const ColMajorMatrix<passivedouble>    &b1DM1,
                                            const ColMajorMatrix<passivedouble>    &b1DP1,
                                            const ColMajorMatrix<passivedouble>    &derB1DM1,
                                            const ColMajorMatrix<passivedouble>    &derB1DP1,
                                            vector<ColMajorMatrix<passivedouble> > &tensorSol,
                                            vector<ColMajorMatrix<passivedouble> > &tensorDSolDr,
                                            vector<ColMajorMatrix<passivedouble> > &tensorDSolDs);

  /*!
   * \brief Function, which sets the function pointer to carry out the tensor
   *        product to obtain the data in volume points.
   * \param[in]  K         - First dimensions of the input tensor.
   * \param[in]  M         - First dimensions of the tensor to be created.
   * \param[out] TPVolData - Function pointer to be set.
   */
  void SetFunctionPointerVolumeDataQuad(const unsigned short K,
                                        const unsigned short M,
                                        TPI2D                &TPVolData);

  /*!
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements.
   */
  void SubConnLinearElements(void);

  /*!
   * \brief Function, which carries out the tensor product C = A*B to obtain the
   *        volume data in the set of points determined by the arguments.
   * \param[in]  TPVolData - Function pointer to carry out the actual multiplication.
   * \param[in]  N         - Number of values to be computed in the integration points.
   * \param[in]  K         - First dimensions of the input tensor. Only used for timings.
   * \param[in]  M         - First dimensions of the tensor to be created. Only used for timings.
   * \param[in]  Ai        - 1D matrix for the i-component of the A tensor.
   * \param[in]  Aj        - 1D matrix for the j-component of the A tensor.
   * \param[in]  B         - B tensor stored as a matrix.
   * \param[out] C         - C tensor stored as a matrix.
   * \param[out] config    - Object used for the timing of the tensor product call.
   */
  void TensorProductVolumeDataQuad(TPI2D                               &TPVolData,
                                   const int                           N,
                                   const int                           K,
                                   const int                           M,
                                   const ColMajorMatrix<passivedouble> &Ai,
                                   const ColMajorMatrix<passivedouble> &Aj,
                                   const ColMajorMatrix<su2double>     &B,
                                   ColMajorMatrix<su2double>           &C,
                                   const CConfig                       *config);
};
