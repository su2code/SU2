/*!
 * \file CFEMStandardPyra.hpp
 * \brief Base class for the FEM pyramid standard element.
 *        The functions are in the <i>CFEMStandardPyra.cpp</i> file.
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
 * \class CFEMStandardPyra
 * \brief Base class which defines the variables and methods for the
 *        pyramid standard element.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardPyra: public CFEMStandardElementBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardPyra() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardPyra(const unsigned short val_nPoly,
                   const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardPyra() = default;

protected:

  unsigned short nDOFs1D;   /*!< \brief Number of DOFs along the line seqment of the base quad. */
  unsigned short nInt1D;    /*!< \brief Number of integration points in one space direction. */

  vector<passivedouble> rPyraDOFsEqui; /*!< \brief Parametric r-coordinates of the pyramid grid
                                                   DOFs when equidistant spacing is used. */
  vector<passivedouble> sPyraDOFsEqui; /*!< \brief Parametric s-coordinates of the pyramid grid
                                                   DOFs when equidistant spacing is used. */
  vector<passivedouble> tPyraDOFsEqui; /*!< \brief Parametric t-coordinates of the pyramid grid
                                                   DOFs when equidistant spacing is used. */

  vector<passivedouble> rPyraDOFsLGL;  /*!< \brief Parametric r-coordinates of the pyramid grid
                                                   DOFs when the LGL distribution is used. */
  vector<passivedouble> sPyraDOFsLGL;  /*!< \brief Parametric s-coordinates of the pyramid grid
                                                   DOFs when the LGL distribution is used. */
  vector<passivedouble> tPyraDOFsLGL;  /*!< \brief Parametric t-coordinates of the pyramid grid
                                                   DOFs when the LGL distribution is used. */

  vector<passivedouble> rLineIntGL;      /*!< \brief 1D parametric coordinates of the
                                                     Gauss-Legendre integration points. */
  vector<passivedouble> wLineIntGL;      /*!< \brief Weights of the 1D Gauss-Legendre
                                                     integration points. */

  vector<passivedouble> rLineIntGJ;      /*!< \brief 1D parametric coordinates of the
                                                     Gauss-Jacobi integration points. */
  vector<passivedouble> wLineIntGJ;      /*!< \brief Weights of the 1D Gauss-Jacobi
                                                     integration points. */
private:

  /*!
   * \brief Function, which determines the location of the grid DOFs of a pyramid for
   *        polynomial degree nPoly when an equidistant spacing is used.
   */
  void LocationPyramidGridDOFsEquidistant();

  /*!
   * \brief Function, which determines the location of the grid DOFs of a pyramid for
   *        polynomial degree nPoly when the LGL distribution is used.
   */
  void LocationPyramidGridDOFsLGL();
};
