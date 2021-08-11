/*!
 * \file TensorProductSurfaceResVolumeDOFs2D_3_4.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (3,4)
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

void TensorProductSurfaceResVolumeDOFs2D_3_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C) {

  /*--- Compute the padded value of the number of 1D DOFs. ---*/
  const int KP = CFEMStandardElementBase::PaddedValue(3);

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
    su2double (*c)[3] = (su2double (*)[3]) &C[l*ldc];

    /*--- Tensor product in tangential direction to obtain the data
          in the DOFs of the line element. ---*/
    SU2_OMP_SIMD
    for(int i=0; i<KP; ++i) bFace[i] = 0.0;
    for(int ii=0; ii<4; ++ii) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(int i=0; i<KP; ++i)
        bFace[i] += aTt[ii][i]*b[ii];
    }

    /*--- Tensor product in normal direction to update the residual.
          This depends on the face ID on which the face data resides. ---*/
    if((faceID == 0) || (faceID == 2)) {

      /*--- The normal direction corresponds to the j-direction.
            Carry out the multiplication accordingly. ----*/
      for(int j=0; j<3; ++j) {
        SU2_OMP_SIMD_IF_NOT_AD_(safelen(3))
        for(int i=0; i<3; ++i)
          c[j][i] += an[j]*bFace[i];
      }
    }
    else {

      /*--- The normal direction corresponds to the i-direction.
            Carry out the multiplication accordingly. ----*/
      for(int j=0; j<3; ++j) {
        SU2_OMP_SIMD_IF_NOT_AD_(safelen(3))
        for(int i=0; i<3; ++i)
          c[j][i] += an[i]*bFace[j];
      }
    }

  } /*--- End of the loop over N. ---*/
}
