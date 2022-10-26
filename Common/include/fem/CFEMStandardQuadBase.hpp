/*!
 * \file CFEMStandardQuadBase.hpp
 * \brief Base class for the FEM quadrilateral standard element.
 *        The functions are in the <i>CFEMStandardQuadBase.cpp</i> file.
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
#include "../tensor_products/TensorProductVolumeIntPoints2D.hpp"

/*!
 * \class CFEMStandardQuadBase
 * \brief Base class which defines the variables and methods for the
 *        quadrilateral standard element.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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
   * \brief Function, which converts the 1D parametric coordinates of a face of
   *        the quad to the tangential and normal components of the actual quad.
   * \param[in]  rLine       - 1D parametric coordinates of a standard line element.
   * \param[in]  faceID      - The corresponding faceID of the adjacent quad.
   * \param[in]  orientation - Orientation of the line element relative to the quad.
   * \param[out] rNormal     - The parametric coordinate in the direction normal to
   *                           the face. This value is either -1 or 1.
   * \param[out] rTangential - The parametric coordinates in tangential to the face.
   */
  void ConvertCoor1DFaceTo2DQuad(const vector<passivedouble> rLine,
                                 const unsigned short        faceID,
                                 const unsigned short        orientation,
                                 vector<passivedouble>       &rNormal,
                                 vector<passivedouble>       &rTangential);

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
   * \brief Function, which creates the connectivity of the linear sub-elements when the
   *        high order element is split in such elements. The splitting is done for one face
   *        w.r.t the volume. The output is stored in subConn1ForPlotting using node ID
   *        available in gridConnFaces.
   */
  void SubConnLinearElementsFace(int val_faceID_Elem);

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
   * \param[in]  initZero  - Whether or not to initialize C to zero.
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
                                   const bool                          initZero,
                                   const CConfig                       *config);

  /*!
   * \brief Function, which creates the local grid connectivities of the faces
   *        of the volume element.
   */
  void LocalGridConnFaces(void);
};
