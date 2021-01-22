/*!
 * \file CFEMStandardTriAdjacentTetGrid.hpp
 * \brief Class for the FEM standard surface triangle element
 *        adjacent to a tetrahedron for the grid.
 *        The functions are in the <i>CFEMStandardTriAdjacentTetGrid.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

#include "CFEMStandardTetBase.hpp"
#include "CGemmStandard.hpp"

/*!
 * \class CFEMStandardTriAdjacentTetGrid.
 * \brief Class which defines the variables and methods for the triangle 
 *        standard surface element adjacent to a tetrahedron for the grid.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
 */
class CFEMStandardTriAdjacentTetGrid final: public CFEMStandardTetBase {

public:

  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardTriAdjacentTetGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial degree that must be integrated exactly.
   * \param[in] val_faceID_Elem - This is the face ID of the adjacent volume element
   *                              to which this surface element corresponds.
   * \param[in] val_orientation - Orientation of this surface element relative to
   *                              the adjacent volume element.
   * \param[in] val_useLGL      - Whether or not the LGL point distribution is used
   *                              for the DOFs.
   * \param[in] val_gemm        - Pointer to the gemm type used to carry out
   *                              the gemm functionality for this standard face.
   */
  CFEMStandardTriAdjacentTetGrid(const unsigned short val_nPoly,
                                 const unsigned short val_orderExact,
                                 const unsigned short val_faceID_Elem,
                                 const unsigned short val_orientation,
                                 const bool           val_useLGL,
                                 CGemmBase           *val_gemm);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardTriAdjacentTetGrid() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public member functions.                                ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes the coordinates in the integration points.
   * \param[in]  notUsed    - Argument present to be consistent with the base class
   *                          function, which is overwritten.
   * \param[in]  matCoorDOF - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matCoorInt - Matrix that contains the coordinates of the integration
   *                          points.
   */
  void CoorIntPoints(const bool                notUsed,
                     ColMajorMatrix<su2double> &matCoorDOF,
                     ColMajorMatrix<su2double> &matCoorInt) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        integration points.
   * \param[in]  notUsed       - Argument present to be consistent with the base class
   *                             function, which is overwritten.
   * \param[in]  matCoorDOF    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoorInt - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorIntPoints(const bool                         notUsed,
                                ColMajorMatrix<su2double>          &matCoorDOF,
                                vector<ColMajorMatrix<su2double> > &matDerCoorInt) override;
private:

  unsigned short faceID_Elem;   /*!< \brief Face ID of the adjacent tetrahedron, which corresponds to this face. */
  unsigned short orientation;   /*!< \brief Orientation of this face relative to the adjacent tetrahedron. */

  CGemmStandard *gemmDOFs2Int = nullptr; /*!< \brief Pointer to the gemm type used to to compute the data in the
                                                     integration points of the face from the volume DOFs. */

  ColMajorMatrix<passivedouble> lagBasisInt;             /*!< \brief The values of the Lagrangian basis functions
                                                                     in the integration points of the face. */
  vector<ColMajorMatrix<passivedouble> > derLagBasisInt; /*!< \brief The values of the derivatives of the Lagrangian
                                                                     basis functions in the integration points. It is a vector,
                                                                     because there are derivatives in three directions. */

  /*-----------------------------------------------------------------------------------*/
  /*---                     Private member functions.                               ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which convertes the gradients w.r.t. the parametric volume coordinates
   *        to the gradients w.r.t. the parametric surface coordinates of the face.
   * \param[in]  matDerVol  - Vector of matrices that contains the derivatives w.r.t.
   *                          the parametric volume coordinates.
   * \param[out] matDerFace - Vector of matrices that contains the derivatives w.r.t.
   *                          the parametric surface coordinates.
   */
  void ConvertVolumeToSurfaceGradients(vector<ColMajorMatrix<su2double> > &matDerVol,
                                       vector<ColMajorMatrix<su2double> > &matDerFace) override;
};
