/*!
 * \file TensorProductVolumeIntPoints3D_10_10.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (10,10)
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

#include "../../../include/tensor_products/TensorProductVolumeIntPoints3D.hpp"

void TensorProductVolumeIntPoints3D_10_10(const int           N,
                                          const int           ldb,
                                          const int           ldc,
                                          const passivedouble *Ai,
                                          const passivedouble *Aj,
                                          const passivedouble *Ak,
                                          const su2double     *B,
                                          su2double           *C) {

  /*--- Cast the one dimensional input arrays for the A-tensor to 2D arrays.
        Note that C++ stores multi-dimensional arrays in row major order,
        hence the indices are reversed compared to the column major order
        storage of e.g. Fortran. ---*/
  const passivedouble (*ai)[16] = (const passivedouble (*)[16]) Ai;
  const passivedouble (*aj)[16] = (const passivedouble (*)[16]) Aj;
  const passivedouble (*ak)[16] = (const passivedouble (*)[16]) Ak;

  /*--- Define the variables to store the intermediate results. ---*/
  su2double tmpK[10][10][16];
  su2double tmpJ[10][10][16];
  su2double tmpI[10][10][16];

  /*--- Outer loop over N. ---*/
  for(int l=0; l<N; ++l) {

    /*--- Cast the index l of B and C to multi-dimensional arrays. ---*/
    const su2double (*b)[10][10] = (const su2double (*)[10][10]) &B[l*ldb];
    su2double       (*c)[10][10] = (su2double (*)[10][10]) &C[l*ldc];

    /*--- Tensor product in k-direction to obtain the solution
          in the integration points in k-direction. ---*/
    for(int i=0; i<10; ++i) {
      for(int j=0; j<10; ++j) {
        SU2_OMP_SIMD
        for(int k=0; k<16; ++k) tmpK[i][j][k] = 0.0;
        for(int kk=0; kk<10; ++kk) {
          SU2_OMP_SIMD
          for(int k=0; k<16; ++k)
            tmpK[i][j][k] += ak[kk][k] * b[kk][j][i];
        }
      }
    }

    /*--- Tensor product in j-direction to obtain the solution
          in the integration points in j-direction. ---*/
    for(int k=0; k<10; ++k) {
      for(int i=0; i<10; ++i) {
        SU2_OMP_SIMD
        for(int j=0; j<16; ++j) tmpJ[k][i][j] = 0.0;
        for(int jj=0; jj<10; ++jj) {
          SU2_OMP_SIMD
          for(int j=0; j<16; ++j)
            tmpJ[k][i][j] += aj[jj][j] * tmpK[i][jj][k];
        }
      }
    }

    /*--- Tensor product in i-direction to obtain the solution
          in the integration points in i-direction. This is
          the final result of the tensor product. ---*/
    for(int k=0; k<10; ++k) {
      for(int j=0; j<10; ++j) {
        SU2_OMP_SIMD
        for(int i=0; i<16; ++i) tmpI[k][j][i] = 0.0;
        for(int ii=0; ii<10; ++ii) {
          SU2_OMP_SIMD
          for(int i=0; i<16; ++i)
            tmpI[k][j][i] += ai[ii][i] * tmpJ[k][ii][j];
        }
      }
    }

    /*--- Copy the values to the appropriate location in c. ---*/
    for(int k=0; k<10; ++k)
      for(int j=0; j<10; ++j)
        for(int i=0; i<10; ++i)
          c[k][j][i] = tmpI[k][j][i];

  } /*--- End of the loop over N. ---*/
}
