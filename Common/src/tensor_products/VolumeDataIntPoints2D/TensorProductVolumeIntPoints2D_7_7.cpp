/*!
 * \file TensorProductVolumeIntPoints2D_7_7.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (7,7)
 * \author Automatically generated file, do not change manually
 * \version 7.0.7 "Blackbird"
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

#include "../../../include/tensor_products/TensorProductVolumeIntPoints2D.hpp"

void TensorProductVolumeIntPoints2D_7_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C) {

  /*--- Cast the one dimensional input arrays for the A-tensor to 2D arrays.
        Note that C++ stores multi-dimensional arrays in row major order,
        hence the indices are reversed compared to the column major order
        storage of e.g. Fortran. ---*/
  const passivedouble (*ai)[8] = (const passivedouble (*)[8]) Ai;
  const passivedouble (*aj)[8] = (const passivedouble (*)[8]) Aj;

  /*--- Define the variables to store the intermediate results. ---*/
  su2double tmpJ[7][8];
  su2double tmpI[7][8];

  /*--- Outer loop over N. ---*/
  for(int l=0; l<N; ++l) {

    /*--- Cast the index l of B and C to multi-dimensional arrays. ---*/
    const su2double (*b)[7] = (const su2double (*)[7]) &B[l*ldb];
    su2double       (*c)[7] = (su2double (*)[7]) &C[l*ldc];

    /*--- Tensor product in j-direction to obtain the solution
          in the integration points in j-direction. ---*/
    for(int i=0; i<7; ++i) {
      SU2_OMP_SIMD
      for(int j=0; j<8; ++j) tmpJ[i][j] = 0.0;
      for(int jj=0; jj<7; ++jj) {
        SU2_OMP_SIMD_IF_NOT_AD
        for(int j=0; j<8; ++j)
          tmpJ[i][j] += aj[jj][j] * b[jj][i];
      }
    }

    /*--- Tensor product in i-direction to obtain the solution
          in the integration points in i-direction. This is
          the final result of the tensor product. ---*/
    for(int j=0; j<7; ++j) {
      SU2_OMP_SIMD
      for(int i=0; i<8; ++i) tmpI[j][i] = 0.0;
      for(int ii=0; ii<7; ++ii) {
        SU2_OMP_SIMD_IF_NOT_AD
        for(int i=0; i<8; ++i)
          tmpI[j][i] += ai[ii][i] * tmpJ[ii][j];
      }
    }

    /*--- Copy the values to the appropriate location in c. ---*/
    for(int j=0; j<7; ++j)
      for(int i=0; i<7; ++i)
        c[j][i] = tmpI[j][i];

  } /*--- End of the loop over N. ---*/
}
