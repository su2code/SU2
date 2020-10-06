/*!
 * \file CFEMStandardTri.hpp
 * \brief Base class for the FEM triangle standard element.
 *        The functions are in the <i>CFEMStandardTri.cpp</i> file.
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

/*!
 * \class CFEMStandardTri
 * \brief Base class which defines the variables and methods for the
 *        triangle standard element.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardTri: public CFEMStandardElementBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardTri() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardTri(const unsigned short val_nPoly,
                  const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardTri() = default;

  /*!
   * \brief Function, which returns the number of faces of the volume element.
   * \return The number of faces of the volume element.
   */
  unsigned short GetNFaces(void) const override {return 3;}

  /*!
   * \brief Function, which returns the VTK type of the given face index.
   * \param[in] ind - Index of the face for which the VTK type must be returned.
   * \return The VTK type of the given face id of the element.
   */
  unsigned short GetVTK_Face(unsigned short ind) const override {return LINE;}

protected:

  vector<passivedouble> rTriangleDOFsEqui; /*!< \brief Parametric r-coordinates of the triangle grid
                                                       DOFs when equidistant spacing is used. */
  vector<passivedouble> sTriangleDOFsEqui; /*!< \brief Parametric s-coordinates of the triangle grid
                                                       DOFs when equidistant spacing is used. */
  vector<passivedouble> rTriangleDOFsLGL;  /*!< \brief Parametric r-coordinates of the triangle grid
                                                       DOFs when the LGL distribution is used. */
  vector<passivedouble> sTriangleDOFsLGL;  /*!< \brief Parametric s-coordinates of the triangle grid
                                                       DOFs when the LGL distribution is used. */

  vector<passivedouble> rTriangleInt;      /*!< \brief Parametric r-coordinates of the integration
                                                       points of the triangle. */
  vector<passivedouble> sTriangleInt;      /*!< \brief Parametric s-coordinates of the integration
                                                       points of the triangle. */
  vector<passivedouble> wTriangleInt;      /*!< \brief Weights of the integration points of the
                                                       triangle. */
};
