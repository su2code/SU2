/*!
 * \file TensorProductSurfaceResVolumeDOFs2D_1_3.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (1,3)
 * \author Automatically generated file, do not change manually
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/tensor_products/TensorProductSurfaceResVolumeDOFs2D.hpp"
#include "../../../include/fem/CFEMStandardElementBase.hpp"

void TensorProductSurfaceResVolumeDOFs2D_1_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C) {

  /*--- Compute the padded value of the number of 1D DOFs. ---*/
  const int KP = CFEMStandardElementBase::PaddedValue(1);

  /*--- Cast the one dimensional input array for the tangential part of
        the AT-tensor to a 2D array. The normal part is a 1D array.
        Note that C++ stores multi-dimensional arrays in row major order,
        hence the indices are reversed compared to the column major order
        storage of e.g. Fortran. ---*/
  const passivedouble *an        = An;
  const passivedouble (*aTt)[KP] = (const passivedouble (*)[KP]) ATt;

  /*--- Define the variable to store the intermediate results. ---*/
  su2double bFace[KP];

  /*--- Outer loop over N. ---*/
  for(int l=0; l<N; ++l) {

    /*--- Cast the index l of B to a 1D. ---*/
    const su2double *b = &B[l*ldb];

    /*--- Cast the index l of C to a multi-dimensional array. ---*/
    su2double (*c)[1] = (su2double (*)[1]) &C[l*ldc];

    /*--- Tensor product in tangential direction to obtain the data
          in the DOFs of the line element. ---*/
    SU2_OMP_SIMD
    for(int i=0; i<KP; ++i) bFace[i] = 0.0;
    for(int ii=0; ii<3; ++ii) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(int i=0; i<KP; ++i)
        bFace[i] += aTt[ii][i]*b[ii];
    }

    /*--- Tensor product in normal direction to update the residual. ---*/
    c[0][0] += an[0]*bFace[0];

  } /*--- End of the loop over N. ---*/
}
