/*!
 * \file TensorProductSurfaceIntPoints2D_1_3.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (1,3)
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

#include "../../../include/tensor_products/TensorProductSurfaceIntPoints2D.hpp"
#include "../../../include/fem/CFEMStandardElementBase.hpp"

void TensorProductSurfaceIntPoints2D_1_3(const int           N,
                                         const int           faceID,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *An,
                                         const passivedouble *At,
                                         const su2double     *B,
                                         su2double           *C) {

  /*--- Compute the padded value of the number of integration points. ---*/
  const int MP = CFEMStandardElementBase::PaddedValue(3);

  /*--- Cast the one dimensional input array for the tangential part of
        the A-tensor to a 2D array. The normal part is a 1D array.
        Note that C++ stores multi-dimensional arrays in row major order,
        hence the indices are reversed compared to the column major order
        storage of e.g. Fortran. ---*/
  const passivedouble *an       = An;
  const passivedouble (*at)[MP] = (const passivedouble (*)[MP]) At;

  /*--- Define the variable to store the intermediate results. ---*/
  su2double bFace[1];

  /*--- Outer loop over N. ---*/
  for(int l=0; l<N; ++l) {

    /*--- Cast the index l of B to a multi-dimensional array. ---*/
    const su2double (*b)[1] = (const su2double (*)[1]) &B[l*ldb];

    /*--- Cast the index l of C to a 1D array. ---*/
    su2double *c = &C[l*ldc];

    /*--- Tensor product in normal direction which generates the data on the face. ---*/
    bFace[0] = an[0]*b[0][0];

    /*--- Tensor product in tangential direction to obtain the data
          in the integration points of the line element. ---*/
    SU2_OMP_SIMD
    for(int i=0; i<MP; ++i) c[i] = 0.0;
    for(int ii=0; ii<1; ++ii) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(int i=0; i<MP; ++i)
        c[i] += at[ii][i]*bFace[ii];
    }

  } /*--- End of the loop over N. ---*/
}
