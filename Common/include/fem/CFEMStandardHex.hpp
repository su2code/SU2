/*!
 * \file CFEMStandardHex.hpp
 * \brief Base class for the FEM hexahedron standard element.
 *        The functions are in the <i>CFEMStandardHex.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
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

#include "CFEMStandardElementBase.hpp"
#include "../tensor_products/TensorProductVolumeIntPoints3D.hpp"

/*!
 * \class CFEMStandardHex
 * \brief Base class which defines the variables and methods for the
 *        hexahedron standard element.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardHex: public CFEMStandardElementBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardHex() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardHex(const unsigned short val_nPoly,
                  const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardHex() = default;

protected:

  unsigned short nDOFs1D;   /*!< \brief Number of DOFs in one space direction. */
  unsigned short nInt1D;    /*!< \brief Number of integration points in one space direction. */

  vector<passivedouble> rLineDOFsEqui; /*!< \brief 1D parametric coordinates of the grid
                                                   DOFs when equidistant spacing is used. */
  vector<passivedouble> rLineDOFsLGL;  /*!< \brief 1D parametric coordinates of the grid
                                                   DOFs when the LGL distribution is used. */
  vector<passivedouble> rLineInt;      /*!< \brief 1D parametric coordinates of the
                                                   integration points. */
  vector<passivedouble> wLineInt;      /*!< \brief Weights of the 1D integration points. */

  /*!
   * \brief Function, which serves as an interface to carry out the tensor product C = A*B
   *        to obtain the data in the volume integration points.
   * \param[in]  N      - Number of values to be computed in the integration points.
   * \param[in]  Ai     - 1D matrix for the i-component of the A tensor.
   * \param[in]  Aj     - 1D matrix for the j-component of the A tensor.
   * \param[in]  Ak     - 1D matrix for the k-component of the A tensor.
   * \param[in]  B      - B tensor stored as a matrix.
   * \param[out] C      - C tensor stored as a matrix.
   * \param[out] config - Object used for the timing of the tensor product call.
   */
  void TensorProductIntegrationPoints(const int                           N,
                                      const ColMajorMatrix<passivedouble> &Ai,
                                      const ColMajorMatrix<passivedouble> &Aj,
                                      const ColMajorMatrix<passivedouble> &Ak,
                                      const ColMajorMatrix<su2double>     &B,
                                      ColMajorMatrix<su2double>           &C,
                                      const CConfig                       *config);

private:
  TPI3D TensorProductDataVolIntPoints; /*!< \brief Function pointer to carry out the tensor product
                                                   to compute the data in the volume integration points. */
};
