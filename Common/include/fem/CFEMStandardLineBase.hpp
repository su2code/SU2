/*!
 * \file CFEMStandardLineBase.hpp
 * \brief Base class for the FEM line standard element.
 *        The functions are in the <i>CFEMStandardLineBase.cpp</i> file.
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

#include "CFEMStandardElementBase.hpp"

/*!
 * \class CFEMStandardLineBase
 * \brief Base class which defines variables and methods for the line standard element.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
class CFEMStandardLineBase: public CFEMStandardElementBase {

protected:
  vector<passivedouble> rLineInt;      /*!< \brief 1D parametric coordinates of the
                                                   integration points. */
  vector<passivedouble> wLineInt;      /*!< \brief Weights of the 1D integration points. */

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CFEMStandardLineBase() = default;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardLineBase(const unsigned short val_nPoly,
                       const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardLineBase() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                           Public member functions.                          ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes the values of the derivatives of the Lagrangian
   *        basis functions of a line in the given integration points for the
   *        given location of the DOFs.
   * \param[in]  rDOFs      - Vector, which contains the parametric locations of the DOFs.
   * \param[in]  rInt       - Vector, which contains the parametric locations of the
   *                          integration points.
   * \param[in]  usePadding - Whether or not to use padding.
   * \param[out] derLag     - Matrix, which contains the values of derivatives of all the
   *                          Lagrangian basis functions in all the integration points.
   */
  void DerLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                const vector<passivedouble>   &rInt,
                                const bool                    usePadding,
                                ColMajorMatrix<passivedouble> &derLag);

  /*!
   * \brief Function, which computes the values of the 2nd derivatives, Hessian,
   *        of the Lagrangian basis functions of a line in the given integration
   *        points for the given location of the DOFs.
   * \param[in]  rDOFs      - Vector, which contains the parametric locations of the DOFs.
   * \param[in]  rInt       - Vector, which contains the parametric locations of the
   *                          integration points.
   * \param[in]  usePadding - Whether or not to use padding.
   * \param[out] hesLag     - Matrix, which contains the values of 2nd derivatives of all
   *                          the Lagrangian basis functions in all the integration points.
   */
  void HesLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                const vector<passivedouble>   &rInt,
                                const bool                    usePadding,
                                ColMajorMatrix<passivedouble> &hesLag);

  /*!
   * \brief Static function, which computes the values of the Lagrangian basis functions
   *        of a line in the given integration points for the given location
   *        of the DOFs.
   * \param[in]  rDOFs      - Vector, which contains the parametric locations of the DOFs.
   * \param[in]  rInt       - Vector, which contains the parametric locations of the
   *                          integration points.
   * \param[in]  usePadding - Whether or not to use padding.
   * \param[out] lag        - Matrix, which contains the values of all the Lagrangian
   *                          basis functions in all the integration points.
   */
  void LagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                             const vector<passivedouble>   &rInt,
                             const bool                    usePadding,
                             ColMajorMatrix<passivedouble> &lag);

protected:

  /*-----------------------------------------------------------------------------------*/
  /*---                         Protected member functions.                         ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which determines the location of the 1D grid DOFs for polynomial
   *        degree mPoly when an equidistant spacing is used.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[out] r     - Vector of the parametric r-coordinates of the DOFs.
   */
  void Location1DGridDOFsEquidistant(const unsigned short  mPoly,
                                     vector<passivedouble> &r);

  /*!
   * \brief Function, which determines the location of the 1D grid DOFs for polynomial
   *        degree mPoly when the LGL distribution is used.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[out] r     - Vector of the parametric r-coordinates of the DOFs.
   */
  void Location1DGridDOFsLGL(const unsigned short  mPoly,
                             vector<passivedouble> &r);

    /*!
   * \brief Function, which computes the gradient of the Vandermonde matrix for a standard 1D edge.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[in]  r     - Parametric coordinates for which the gradient of the Vandermonde
   *                     matrix must be computed.
   * \param[out] VDr   - Matrix to store the gradient of the Vandermonde matrix in all r-locations.
   */
  void GradVandermonde1D(const unsigned short          mPoly,
                         const vector<passivedouble>   &r,
                         ColMajorMatrix<passivedouble> &VDr);

  /*!
   * \brief Function, which computes the Hessian (2nd derivative) of the Vandermonde
   *        matrix for a standard 1D edge.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[in]  r     - Parametric coordinates for which the Hessian of the Vandermonde
   *                     matrix must be computed.
   * \param[out] VD2r  - Matrix to store the Hessian of the Vandermonde matrix in all r-locations.
   */
  void HesVandermonde1D(const unsigned short          mPoly,
                        const vector<passivedouble>   &r,
                        ColMajorMatrix<passivedouble> &VD2r);

  /*!
   * \brief Function, which computes the Vandermonde matrix for a standard 1D edge.
   * \param[in]  mPoly - Polynomial order of the element.
   * \param[in]  r     - Parametric coordinates for which the Vandermonde matrix must be computed.
   * \param[out] V     - Matrix to store the Vandermonde matrix in all r-locations.
   */
  void Vandermonde1D(const unsigned short          mPoly,
                     const vector<passivedouble>   &r,
                     ColMajorMatrix<passivedouble> &V);

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

private:
  /*-----------------------------------------------------------------------------------*/
  /*---                          Private member functions.                          ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes the value of the unscaled Legendre polynomial for the given x-coordinate.
   * \param[in] n     - Order of the Jacobi polynomial.
   * \param[in] x     - Coordinate (-1 <= x <= 1) for which the Legendre polynomial must be evaluated
   * \return            The value of the unscaled Legendre polynomial of order n for the given value of x
   */
  passivedouble Legendre(unsigned short n,
                         passivedouble  x);
};
