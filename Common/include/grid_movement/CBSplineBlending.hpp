/*!
 * \file CBSplineBlending.hpp
 * \brief Headers of the CBSplineBlending class.
 *        Defines blending using uniform BSplines
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFreeFormBlending.hpp"
#include <vector>

using namespace std;

/*!
 * \class CBSplineBlending
 * \brief Class that defines the blending using uniform BSplines.
 * \author T. Albring
 */
class CBSplineBlending : public CFreeFormBlending {
 private:
  vector<su2double> U;          /*!< \brief The knot vector for uniform BSplines on the interval [0,1]. */
  vector<vector<su2double> > N; /*!< \brief The temporary matrix holding the j+p basis functions up to order p. */
  unsigned short KnotSize;      /*!< \brief The size of the knot vector. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CBSplineBlending(short val_order, short n_controlpoints);

  /*!
   * \brief Destructor of the class.
   */
  ~CBSplineBlending() override;

  /*!
   * \brief Returns the value of the i-th basis function and stores the values of the i+p basis functions in the matrix
   * N. \param[in] val_i - index of the basis function. \param[in] val_t - Point at which we want to evaluate the i-th
   * basis.
   */
  su2double GetBasis(short val_i, su2double val_t) override;

  /*!
   * \brief Returns the value of the derivative of the i-th basis function.
   * \param[in] val_i - index of the basis function.
   * \param[in] val_t - Point at which we want to evaluate the derivative of the i-th basis.
   * \param[in] val_order - Order of the derivative.
   */
  su2double GetDerivative(short val_i, su2double val_t, short val_order_der) override;

  /*!
   * \brief Set the order and number of control points.
   * \param[in] val_order - The new order of the function.
   * \param[in] n_controlpoints - the new number of control points.
   */
  void SetOrder(short val_order, short n_controlpoints) override;
};
