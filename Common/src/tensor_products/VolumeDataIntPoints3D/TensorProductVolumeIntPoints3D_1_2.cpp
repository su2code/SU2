/*!
 * \file TensorProductVolumeIntPoints3D_1_2.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (1,2)
 * \author Automatically generated file, do not change manually
 * \version 7.0.8 "Blackbird"
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
#include "../../../include/fem/CFEMStandardElementBase.hpp"

void TensorProductVolumeIntPoints3D_1_2(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const passivedouble *Ak,
                                        const su2double     *B,
                                        su2double           *C) {

  /*--- Compute the padded value of the number of integration points. ---*/
  const size_t baseVectorLen = CFEMStandardElementBase::baseVectorLen;
  const int MP = ((1+baseVectorLen)/baseVectorLen)*baseVectorLen;

  /*--- Cast the one dimensional input arrays for the A-tensor to 2D arrays.
        Note that C++ stores multi-dimensional arrays in row major order,
        hence the indices are reversed compared to the column major order
        storage of e.g. Fortran. ---*/
  const passivedouble (*ai)[MP] = (const passivedouble (*)[MP]) Ai;
  const passivedouble (*aj)[MP] = (const passivedouble (*)[MP]) Aj;
  const passivedouble (*ak)[MP] = (const passivedouble (*)[MP]) Ak;

  /*--- Define the variables to store the intermediate results. ---*/
  su2double tmpK[1][1][MP];
  su2double tmpJ[2][1][MP];
#if MP > 2
  su2double tmpI[2][2][MP];
#endif

  /*--- Outer loop over N. ---*/
  for(int l=0; l<N; ++l) {

    /*--- Cast the index l of B and C to multi-dimensional arrays. ---*/
    const su2double (*b)[1][1] = (const su2double (*)[1][1]) &B[l*ldb];
    su2double       (*c)[2][2] = (su2double (*)[2][2]) &C[l*ldc];

    /*--- Tensor product in k-direction to obtain the solution
          in the integration points in k-direction. ---*/
    for(int i=0; i<1; ++i) {
      for(int j=0; j<1; ++j) {
        SU2_OMP_SIMD
        for(int k=0; k<MP; ++k) tmpK[i][j][k] = 0.0;
        for(int kk=0; kk<1; ++kk) {
          SU2_OMP_SIMD_IF_NOT_AD
          for(int k=0; k<MP; ++k)
            tmpK[i][j][k] += ak[kk][k] * b[kk][j][i];
        }
      }
    }

    /*--- Tensor product in j-direction to obtain the solution
          in the integration points in j-direction. ---*/
    for(int k=0; k<2; ++k) {
      for(int i=0; i<1; ++i) {
        SU2_OMP_SIMD
        for(int j=0; j<MP; ++j) tmpJ[k][i][j] = 0.0;
        for(int jj=0; jj<1; ++jj) {
          SU2_OMP_SIMD_IF_NOT_AD
          for(int j=0; j<MP; ++j)
            tmpJ[k][i][j] += aj[jj][j] * tmpK[i][jj][k];
        }
      }
    }

    /*--- Tensor product in i-direction to obtain the solution
          in the integration points in i-direction. This is
          the final result of the tensor product. ---*/
    for(int k=0; k<2; ++k) {
      for(int j=0; j<2; ++j) {
#if MP > 2
        SU2_OMP_SIMD
        for(int i=0; i<MP; ++i) tmpI[k][j][i] = 0.0;
        for(int ii=0; ii<1; ++ii) {
          SU2_OMP_SIMD_IF_NOT_AD
          for(int i=0; i<MP; ++i)
            tmpI[k][j][i] += ai[ii][i] * tmpJ[k][ii][j];
#else
        SU2_OMP_SIMD
        for(int i=0; i<MP; ++i) c[k][j][i] = 0.0;
        for(int ii=0; ii<1; ++ii) {
          SU2_OMP_SIMD_IF_NOT_AD
          for(int i=0; i<MP; ++i)
            c[k][j][i] += ai[ii][i] * tmpJ[k][ii][j];
#endif
        }
      }
    }

#if MP > 2
    /*--- Copy the values to the appropriate location in c. ---*/
    for(int k=0; k<2; ++k)
      for(int j=0; j<2; ++j)
        for(int i=0; i<2; ++i)
          c[k][j][i] = tmpI[k][j][i];
#endif

  } /*--- End of the loop over N. ---*/
}
