/*!
 * \file CBezierBlending.hpp
 * \brief Headers of the CBezierBlending class.
 *        Defines blending using Bernsteinpolynomials (Bezier Curves)
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
 * \class CBezierBlending
 * \brief Class that defines the blending using Bernsteinpolynomials (Bezier Curves).
 * \author F. Palacios, T. Albring
 */
class CBezierBlending : public CFreeFormBlending {
 private:
  vector<su2double> binomial; /*!< \brief Temporary vector for the Bernstein evaluation. */

  /*!
   * \brief Returns the value of the i-th Bernstein polynomial of order n.
   * \param[in] val_n - Order of the Bernstein polynomial.
   * \param[in] val_i - index of the basis function.
   * \param[in] val_t - Point at which we want to evaluate the i-th basis.
   */
  su2double GetBernstein(short val_n, short val_i, su2double val_t);

  /*!
   * \brief Returns the value of the derivative of the i-th Bernstein polynomial of order n.
   * \param[in] val_n - Order of the Bernstein polynomial.
   * \param[in] val_i - index of the basis function.
   * \param[in] val_t - Point at which we want to evaluate the i-th basis.
   * \param[in] val_order - Order of the derivative.
   */
  su2double GetBernsteinDerivative(short val_n, short val_i, su2double val_t, short val_order_der);

  /*!
   * \brief Get the binomial coefficient n over i, defined as n!/(m!(n-m)!)
   * \note If the denominator is 0, the value is 1.
   * \param[in] n - Upper coefficient.
   * \param[in] m - Lower coefficient.
   * \return Value of the binomial coefficient n over m.
   */
  su2double Binomial(unsigned short n, unsigned short m);

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_order - Max. order of the basis functions.
   * \param[in] n_controlpoints - Not used here.
   */
  CBezierBlending(short val_order, short n_controlpoints);

  /*!
   * \brief Destructor of the class.
   */
  ~CBezierBlending() override;

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
