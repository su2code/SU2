/*!
 * \file CBSplineBlending.cpp
 * \brief Subroutines for B-Spline blening for FFDs
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

#include "../../include/grid_movement/CBSplineBlending.hpp"
#include "../../include/parallelization/mpi_structure.hpp"

CBSplineBlending::CBSplineBlending(short val_order, short n_controlpoints) : CFreeFormBlending() {
  SetOrder(val_order, n_controlpoints);
}

CBSplineBlending::~CBSplineBlending() = default;

void CBSplineBlending::SetOrder(short val_order, short n_controlpoints) {
  unsigned short iKnot;

  Order = val_order;
  Degree = Order - 1;
  nControl = n_controlpoints;

  KnotSize = Order + nControl;

  U.resize(KnotSize, 0.0);

  /*--- Set up the knot vectors for open uniform B-Splines ---*/

  /*--- Note: the first knots are zero now.---*/

  /*--- The next knots are equidistantly distributed in [0,1] ---*/

  for (iKnot = 0; iKnot < nControl - Order; iKnot++) {
    U[Order + iKnot] = su2double(iKnot + 1) / su2double(nControl - Order + 1);
  }
  for (iKnot = nControl - Order; iKnot < nControl; iKnot++) {
    U[Order + iKnot] = 1.0;
  }

  /*--- Allocate the temporary vectors for the basis evaluation ---*/

  N.resize(Order, vector<su2double>(Order, 0.0));
}

su2double CBSplineBlending::GetBasis(short val_i, su2double val_t) {
  /*--- Evaluation is based on the algorithm from "The NURBS Book (Les Piegl and Wayne Tiller)" ---*/

  /*--- Special cases ---*/

  if ((val_i == 0 && val_t == U[0]) || (val_i == (short)U.size() - 1 && val_t == U.back())) {
    return 1.0;
  }

  /*--- Local property of BSplines ---*/

  if ((val_t < U[val_i]) || (val_t >= U[val_i + Order])) {
    return 0.0;
  }

  unsigned short j, k;
  su2double saved, temp;

  for (j = 0; j < Order; j++) {
    if ((val_t >= U[val_i + j]) && (val_t < U[val_i + j + 1]))
      N[j][0] = 1.0;
    else
      N[j][0] = 0;
  }

  for (k = 1; k < Order; k++) {
    if (N[0][k - 1] == 0.0)
      saved = 0.0;
    else
      saved = ((val_t - U[val_i]) * N[0][k - 1]) / (U[val_i + k] - U[val_i]);
    for (j = 0; j < Order - k; j++) {
      if (N[j + 1][k - 1] == 0.0) {
        N[j][k] = saved;
        saved = 0.0;
      } else {
        temp = N[j + 1][k - 1] / (U[val_i + j + k + 1] - U[val_i + j + 1]);
        N[j][k] = saved + (U[val_i + j + k + 1] - val_t) * temp;
        saved = (val_t - U[val_i + j + 1]) * temp;
      }
    }
  }
  return N[0][Order - 1];
}

su2double CBSplineBlending::GetDerivative(short val_i, su2double val_t, short val_order_der) {
  if ((val_t < U[val_i]) || (val_t >= U[val_i + Order])) {
    return 0.0;
  }

  /*--- Evaluate the i+p basis functions up to the order p (stored in the matrix N). ---*/

  GetBasis(val_i, val_t);

  /*--- Use the recursive definition for the derivative (hardcoded for 1st and 2nd derivative). ---*/

  if (val_order_der == 0) {
    return N[0][Order - 1];
  }

  if (val_order_der == 1) {
    return (Order - 1.0) / (1e-10 + U[val_i + Order - 1] - U[val_i]) * N[0][Order - 2] -
           (Order - 1.0) / (1e-10 + U[val_i + Order] - U[val_i + 1]) * N[1][Order - 2];
  }

  if (val_order_der == 2 && Order > 2) {
    const su2double left = (Order - 2.0) / (1e-10 + U[val_i + Order - 2] - U[val_i]) * N[0][Order - 3] -
                           (Order - 2.0) / (1e-10 + U[val_i + Order - 1] - U[val_i + 1]) * N[1][Order - 3];

    const su2double right = (Order - 2.0) / (1e-10 + U[val_i + Order - 1] - U[val_i + 1]) * N[1][Order - 3] -
                            (Order - 2.0) / (1e-10 + U[val_i + Order] - U[val_i + 2]) * N[2][Order - 3];

    return (Order - 1.0) / (1e-10 + U[val_i + Order - 1] - U[val_i]) * left -
           (Order - 1.0) / (1e-10 + U[val_i + Order] - U[val_i + 1]) * right;
  }

  /*--- Higher order derivatives are not implemented, so we exit if they are requested. ---*/

  if (val_order_der > 2) {
    SU2_MPI::Error("Higher order derivatives for BSplines are not implemented.", CURRENT_FUNCTION);
  }
  return 0.0;
}
