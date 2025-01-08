/*!
 * \file CBezierBlending.cpp
 * \brief Subroutines for Bezier blending for FFDs
 * \author F. Palacios, T. Economon, S. Padron
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

#include "../../include/grid_movement/CBezierBlending.hpp"
#include "../../include/option_structure.hpp"

CBezierBlending::CBezierBlending(short val_order, short n_controlpoints) { SetOrder(val_order, n_controlpoints); }

CBezierBlending::~CBezierBlending() = default;

void CBezierBlending::SetOrder(short val_order, short n_controlpoints) {
  Order = val_order;
  Degree = Order - 1;
  binomial.resize(Order + 1, 0.0);
}

su2double CBezierBlending::GetBasis(short val_i, su2double val_t) { return GetBernstein(Degree, val_i, val_t); }

su2double CBezierBlending::GetBernstein(short val_n, short val_i, su2double val_t) {
  su2double value = 0.0;

  if (val_i > val_n) {
    value = 0.0;
    return value;
  }

  if (val_i == 0) {
    if (val_t == 0)
      value = 1.0;
    else if (val_t == 1)
      value = 0.0;
    else
      value = Binomial(val_n, val_i) * pow(val_t, val_i) * pow(1.0 - val_t, val_n - val_i);
  } else if (val_i == val_n) {
    if (val_t == 0)
      value = 0.0;
    else if (val_t == 1)
      value = 1.0;
    else
      value = pow(val_t, val_n);
  } else {
    if ((val_t == 0) || (val_t == 1)) value = 0.0;
    value = Binomial(val_n, val_i) * pow(val_t, val_i) * pow(1.0 - val_t, val_n - val_i);
  }

  return value;
}

su2double CBezierBlending::GetDerivative(short val_i, su2double val_t, short val_order_der) {
  return GetBernsteinDerivative(Degree, val_i, val_t, val_order_der);
}

su2double CBezierBlending::GetBernsteinDerivative(short val_n, short val_i, su2double val_t, short val_order_der) {
  su2double value = 0.0;

  /*--- Verify this subroutine, it provides negative val_n,
   which is a wrong value for GetBernstein ---*/

  if (val_order_der == 0) {
    value = GetBernstein(val_n, val_i, val_t);
    return value;
  }

  if (val_i == 0) {
    value = val_n * (-GetBernsteinDerivative(val_n - 1, val_i, val_t, val_order_der - 1));
    return value;
  }
  if (val_n == 0) {
    value = val_t;
    return value;
  } else {
    value = val_n * (GetBernsteinDerivative(val_n - 1, val_i - 1, val_t, val_order_der - 1) -
                     GetBernsteinDerivative(val_n - 1, val_i, val_t, val_order_der - 1));
    return value;
  }

  return value;
}

su2double CBezierBlending::Binomial(unsigned short n, unsigned short m) {
  unsigned short i, j;
  su2double result;

  binomial[0] = 1.0;
  for (i = 1; i <= n; ++i) {
    binomial[i] = 1.0;
    for (j = i - 1U; j > 0; --j) {
      binomial[j] += binomial[j - 1U];
    }
  }

  result = binomial[m];
  if (fabs(result) < EPS * EPS) {
    result = 0.0;
  }

  return result;
}
