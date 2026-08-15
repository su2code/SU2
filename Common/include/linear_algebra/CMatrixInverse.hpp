/*!
 * \file CMatrixInverse.hpp
 * \brief Dense small-matrix inversion via Gauss-Jordan elimination, shared between the host
 *        (CSysMatrix::MatrixInverse) and device (CSysPreconditionerGPU.cu) implementations.
 * \author F. Palacios, A. Bueno, T. Economon, P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include <cmath>

#ifdef __CUDACC__
#define SU2_CUDA_HOST_DEVICE __host__ __device__
#else
#define SU2_CUDA_HOST_DEVICE
#endif

namespace SU2_LinAlg {

/*!
 * \brief Regularize a pivot that is too small to prevent divide-by-zero, on host and device
 *        this needs to clamp to the same value so that the two produce the same factors.
 */
template <class ScalarType>
SU2_CUDA_HOST_DEVICE inline void RegularizePivot(ScalarType& pivot) {
  const float eps = 1e-12;
#ifdef __CUDA_ARCH__
  if (fabs(pivot) < eps) pivot = copysign(ScalarType(eps), pivot);
#else
  if (std::abs(pivot) < eps) pivot = std::copysign(ScalarType(eps), pivot);
#endif
}

/*!
 * \brief Invert the \p nVar by \p nVar dense matrix \p matrix into \p inverse via Gauss-Jordan
 *        elimination with partial pivoting on the diagonal.
 * \note \p matrix is used as scratch space and destroyed, \p inverse must not alias it.
 */
template <class ScalarType>
SU2_CUDA_HOST_DEVICE inline void MatrixInverse(unsigned long nVar, ScalarType* matrix, ScalarType* inverse) {
#define A(I, J) matrix[(I)*nVar + (J)]
#define M(I, J) inverse[(I)*nVar + (J)]

  /*--- Initialize the inverse with the identity. ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++) M(iVar, jVar) = ScalarType(iVar == jVar);

  /*--- Transform system in Upper Matrix. ---*/
  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++) {
      RegularizePivot(A(jVar, jVar));

      const ScalarType weight = A(iVar, jVar) / A(jVar, jVar);
      for (auto kVar = jVar; kVar < nVar; kVar++) A(iVar, kVar) -= weight * A(jVar, kVar);

      /*--- At this stage M is lower triangular so not all cols need updating. ---*/
      for (auto kVar = 0ul; kVar <= jVar; kVar++) M(iVar, kVar) -= weight * M(jVar, kVar);
    }
  }

  /*--- Backwards substitution. ---*/
  for (auto iVar = nVar; iVar > 0ul;) {
    iVar--;  // unsigned type
    for (auto jVar = iVar + 1; jVar < nVar; jVar++)
      for (auto kVar = 0ul; kVar < nVar; kVar++) M(iVar, kVar) -= A(iVar, jVar) * M(jVar, kVar);

    RegularizePivot(A(iVar, iVar));

    for (auto kVar = 0ul; kVar < nVar; kVar++) M(iVar, kVar) /= A(iVar, iVar);
  }

#undef A
#undef M
}

}  // namespace SU2_LinAlg

#undef SU2_CUDA_HOST_DEVICE
