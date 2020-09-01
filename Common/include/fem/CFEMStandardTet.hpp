/*!
 * \file CFEMStandardTet.hpp
 * \brief Base class for the FEM tetrahedron standard element.
 *        The functions are in the <i>CFEMStandardTet.cpp</i> file.
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
 * \class CFEMStandardTet
 * \brief Base class which defines the variables and methods for the
 *        tetrahedron standard element.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardTet: public CFEMStandardElementBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardTet() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardTet(const unsigned short val_nPoly,
                            const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardTet() = default;

protected:

  vector<passivedouble> rTetDOFsEqui; /*!< \brief Parametric r-coordinates of the tetrahedral grid
                                                  DOFs when equidistant spacing is used. */
  vector<passivedouble> sTetDOFsEqui; /*!< \brief Parametric s-coordinates of the tetrahedral grid
                                                  DOFs when equidistant spacing is used. */
  vector<passivedouble> tTetDOFsEqui; /*!< \brief Parametric t-coordinates of the tetrahedral grid
                                                  DOFs when equidistant spacing is used. */

  vector<passivedouble> rTetDOFsLGL;  /*!< \brief Parametric r-coordinates of the tetrahedral grid
                                                  DOFs when the LGL distribution is used. */
  vector<passivedouble> sTetDOFsLGL;  /*!< \brief Parametric s-coordinates of the tetrahedral grid
                                                  DOFs when the LGL distribution is used. */
  vector<passivedouble> tTetDOFsLGL;  /*!< \brief Parametric t-coordinates of the tetrahedral grid
                                                  DOFs when the LGL distribution is used. */

  vector<passivedouble> rTetInt;      /*!< \brief Parametric r-coordinates of the integration
                                                  points of the tetrahedron. */
  vector<passivedouble> sTetInt;      /*!< \brief Parametric s-coordinates of the integration
                                                  points of the tetrahedron. */
  vector<passivedouble> tTetInt;      /*!< \brief Parametric t-coordinates of the integration
                                                  points of the tetrahedron. */
  vector<passivedouble> wTetInt;      /*!< \brief Weights of the integration points of the
                                                  tetrahedron. */

private:

  /*!
   * \brief Function, which computes the warping factors of a 3D triangular face,
            needed to compute the LGL distribution of a tetrahedron.
   * \param[in]  alpha - Blending factor for the given polynomial degree.
   * \param[in]  L1    - Barycentric coordinate of the DOFs (I think).
   * \param[in]  L2    - Barycentric coordinate of the DOFs (I think).
   * \param[in]  L3    - Barycentric coordinate of the DOFs (I think).
   * \param[out] dx    - Delta of the first parametric coordinate of the face.
   * \param[out] dy    - Delta of the second parametric coordinate of the face.
   */
  void EvalShift(const passivedouble         alpha,
                 const vector<passivedouble> &L1,
                 const vector<passivedouble> &L2,
                 const vector<passivedouble> &L3,
                 vector<passivedouble>       &dx,
                 vector<passivedouble>       &dy);

  /*!
   * \brief Function, which computes the 1D edge warping.
   * \param[in]  xOut - Coordinates to which the 1D edge warping must be applied.
   * \param[out] warp - Edge warping value for all the nodal DOFs.
   */
  void EvalWarp(const vector<passivedouble> &xOut,
                vector<passivedouble>       &warp);

  /*!
   * \brief Function, which determines the location of the grid DOFs of a tetrahedron
   *        for polynomial degree nPoly when an equidistant spacing is used.
   * \param[out] r - Parametric r-coordinates of the DOFs.
   * \param[out] s - Parametric s-coordinates of the DOFs.
   * \param[out] t - Parametric t-coordinates of the DOFs.
   */
  void LocationTetGridDOFsEquidistant(vector<passivedouble> &r,
                                      vector<passivedouble> &s,
                                      vector<passivedouble> &t);

  /*!
   * \brief Function, which determines the location of the grid DOFs of a tetrahedron
   *        for polynomial degree nPoly when the LGL distribution is used.
   */
  void LocationTetGridDOFsLGL();
};
