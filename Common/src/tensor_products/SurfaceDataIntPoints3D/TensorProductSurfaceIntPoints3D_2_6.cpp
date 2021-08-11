/*!
 * \file TensorProductSurfaceIntPoints3D_2_6.cpp
 * \brief Function, which carries out the tensor product for (nDOFs1D,nInt1D) = (2,6)
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

#include "../../../include/tensor_products/TensorProductSurfaceIntPoints3D.hpp"
#include "../../../include/fem/CFEMStandardElementBase.hpp"

void TensorProductSurfaceIntPoints3D_2_6(const int           N,
                                         const int           faceID,
                                         const int           ldb,
                                         const int           ldc,
                                         const bool          swapTangDir,
                                         const passivedouble *An,
                                         const passivedouble *At0,
                                         const passivedouble *At1,
                                         const su2double     *B,
                                         su2double           *C) {

  /*--- Compute the padded value of the number of integration points. ---*/
  const int MP = CFEMStandardElementBase::PaddedValue(6);

  /*--- Cast the one dimensional input arrays for the tangential parts
        of the A-tensor to 2D arrays. The normal part is a 1D array.
        Note that C++ stores multi-dimensional arrays in row major order,
        hence the indices are reversed compared to the column major order
        storage of e.g. Fortran. ---*/
  const passivedouble *an        = An;
  const passivedouble (*at0)[MP] = (const passivedouble (*)[MP]) At0;
  const passivedouble (*at1)[MP] = (const passivedouble (*)[MP]) At1;

  /*--- Define the variables to store the intermediate results. ---*/
  su2double bFace[2][2];
  su2double tmpI[2][MP], tmpJ[6][MP];

  /*--- Outer loop over N. ---*/
  for(int l=0; l<N; ++l) {

    /*--- Cast the index l of B to a 3D array. ---*/
    const su2double (*b)[2][2] = (const su2double (*)[2][2]) &B[l*ldb];

    /*--- Cast the index l of C to a 2D array. ---*/
    su2double (*c)[6] = (su2double (*)[6]) &C[l*ldc];

    /*--- Tensor product in normal direction which generates the data on the face
          This depends on the face ID on which the data must be generated. ---*/
    for(int j=0; j<2; ++j) {
      SU2_OMP_SIMD_(safelen(2))
      for(int i=0; i<2; ++i) bFace[j][i] = 0.0;
    }

    if((faceID == 0) || (faceID == 1)) {

      /*--- The normal direction corresponds to the k-direction.
            Carry out the multiplication accordingly. ----*/
      for(int k=0; k<2; ++k) {
        for(int j=0; j<2; ++j) {
          SU2_OMP_SIMD_IF_NOT_AD_(safelen(2))
          for(int i=0; i<2; ++i)
            bFace[j][i] += an[k]*b[k][j][i];
        }
      }
    }
    else if((faceID == 2) || (faceID == 3)) {

      /*--- The normal direction corresponds to the j-direction.
            Carry out the multiplication accordingly. ----*/
      for(int k=0; k<2; ++k) {
        for(int j=0; j<2; ++j) {
          SU2_OMP_SIMD_IF_NOT_AD_(safelen(2))
          for(int i=0; i<2; ++i)
            bFace[k][i] += an[j]*b[k][j][i];
        }
      }
    }
    else {

      /*--- The normal direction corresponds to the i-direction.
            Carry out the multiplication accordingly. ----*/
      for(int k=0; k<2; ++k) {
        for(int i=0; i<2; ++i) {
          SU2_OMP_SIMD_IF_NOT_AD_(safelen(2))
          for(int j=0; j<2; ++j)
            bFace[k][j] += an[i]*b[k][j][i];
        }
      }
    }

    /*--- Tensor product in first tangential direction to obtain the data
          in the integration points in the i-direction of the face. ---*/
    for(int j=0; j<2; ++j) {
      SU2_OMP_SIMD
      for(int i=0; i<MP; ++i) tmpI[j][i] = 0.0;
      for(int ii=0; ii<2; ++ii) {
        SU2_OMP_SIMD_IF_NOT_AD
        for(int i=0; i<MP; ++i)
          tmpI[j][i] += at0[ii][i] * bFace[j][ii];
      }
    }

    /*--- Tensor product in second tangential direction to obtain the data
          in the integration points in the j-direction of the face. ---*/
    for(int i=0; i<6; ++i) {
      SU2_OMP_SIMD
      for(int j=0; j<MP; ++j) tmpJ[i][j] = 0.0;
      for(int jj=0; jj<2; ++jj) {
        SU2_OMP_SIMD_IF_NOT_AD
        for(int j=0; j<MP; ++j)
          tmpJ[i][j] += at1[jj][j] * tmpI[jj][i];
      }
    }

    /*--- Copy the value to the appropriate location in c.
          Take a possible swapping into account. ---*/
    if( swapTangDir ) {
      for(int j=0; j<6; ++j)
        for(int i=0; i<6; ++i)
          c[j][i] = tmpJ[j][i];
    }
    else {
      for(int j=0; j<6; ++j)
        for(int i=0; i<6; ++i)
          c[j][i] = tmpJ[i][j];
    }

  } /*--- End of the loop over N. ---*/
}
