/*!
 * \file CFreeFormBlending.hpp
 * \brief Headers of the CFreeFormBlending class.
 *        It is the parent class for the FFD blending function
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

/*!
 * \class CFreeFormBlending
 * \brief Class that defines the particular kind of blending function for the free form deformation.
 * \author T. Albring
 */

#pragma once

#include "../basic_types/datatype_structure.hpp"

class CFreeFormBlending {
 protected:
  unsigned short Order, /*!< \brief Order of the polynomial basis. */
      Degree,           /*!< \brief Degree (Order - 1) of the polynomial basis. */
      nControl;         /*!< \brief Number of control points. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CFreeFormBlending();

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFreeFormBlending();

  /*!
   * \brief A pure virtual member.
   * \param[in] val_i - index of the basis function.
   * \param[in] val_t - Point at which we want to evaluate the i-th basis.
   */
  inline virtual su2double GetBasis(short val_i, su2double val_t) { return 0.0; }

  /*!
   * \brief A pure virtual member.
   * \param[in] val_i - index of the basis function.
   * \param[in] val_t - Point at which we want to evaluate the derivative of the i-th basis.
   * \param[in] val_order - Order of the derivative.
   */
  inline virtual su2double GetDerivative(short val_i, su2double val_t, short val_order) { return 0.0; }

  /*!
   * \brief A pure virtual member.
   * \param[in] val_order - The new order of the function.
   * \param[in] n_controlpoints - the new number of control points.
   */
  inline virtual void SetOrder(short val_order, short n_controlpoints) {}

  /*!
   * \brief Returns the current order of the function.
   */
  inline su2double GetOrder() const { return Order; }

  /*!
   * \brief Returns the current degree of the function.
   */
  inline su2double GetDegree() const { return Degree; }
};
