/*!
 * \file TensorProductVolumeIntPoints2D.hpp
 * \brief Function prototypes for the tensor product to compute the data in the 2D integration points
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

#pragma once

#include <iostream>
#include <map>
#include "../basic_types/datatype_structure.hpp"
#include "../parallelization/omp_structure.hpp"
#include "../toolboxes/classes_multiple_integers.hpp"

using namespace std;

typedef void(*TPI2D)(const int           N,
                     const int           ldb,
                     const int           ldc,
                     const passivedouble *Ai,
                     const passivedouble *Aj,
                     const su2double     *B,
                     su2double           *C);
/*!
 * \brief Function, which stores the available function pointers for the tensor
 *        product for the 2D volume integration points in a map.
 * \param[out] mapFunctions - Map to store the function pointers to carry out the tensor product.
 */
void CreateMapTensorProductVolumeIntPoints2D(map<CUnsignedShort2T, TPI2D> &mapFunctions);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (1,2).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_1_2(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (1,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_1_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (1,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_1_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (1,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_1_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,1).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_1(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,2).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_2(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_9(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,10).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_10(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,11).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_11(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,12).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_12(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,13).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_13(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,14).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_14(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (2,15).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_2_15(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,1).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_1(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,2).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_2(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (3,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_3_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,1).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_1(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,2).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_2(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (4,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_4_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,1).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_1(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,2).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_2(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (5,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_5_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (6,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_6_9(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_9(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (7,3).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_7_3(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,4).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_4(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,5).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_5(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_8(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_9(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,10).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_10(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,11).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_11(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,12).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_12(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (8,13).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_8_13(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (9,6).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_9_6(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (9,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_9_7(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (9,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_9_9(const int           N,
                                        const int           ldb,
                                        const int           ldc,
                                        const passivedouble *Ai,
                                        const passivedouble *Aj,
                                        const su2double     *B,
                                        su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (9,13).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_9_13(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (9,14).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_9_14(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (10,10).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_10_10(const int           N,
                                          const int           ldb,
                                          const int           ldc,
                                          const passivedouble *Ai,
                                          const passivedouble *Aj,
                                          const su2double     *B,
                                          su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (10,14).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_10_14(const int           N,
                                          const int           ldb,
                                          const int           ldc,
                                          const passivedouble *Ai,
                                          const passivedouble *Aj,
                                          const su2double     *B,
                                          su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (11,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_11_8(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (12,7).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_12_7(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (12,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_12_8(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (13,8).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_13_8(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (13,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_13_9(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (14,9).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_14_9(const int           N,
                                         const int           ldb,
                                         const int           ldc,
                                         const passivedouble *Ai,
                                         const passivedouble *Aj,
                                         const su2double     *B,
                                         su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to obtain the data
 *        in the 2D integration points for (nDOFs1D,nInt1D) = (14,10).
 * \param[in]  N   - Number of variables to be determined in the integration points
 * \param[in]  ldb - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc - Leading dimension of C when stored as a matrix.
 * \param[in]  Ai  - I-componnent of the A tensor.
 * \param[in]  Aj  - J-componnent of the A tensor.
 * \param[in]  B   - Tensor, which contains the data to be interpolated.
 * \param[out] C   - Result of the tensor product C = A*B.
 */
void TensorProductVolumeIntPoints2D_14_10(const int           N,
                                          const int           ldb,
                                          const int           ldc,
                                          const passivedouble *Ai,
                                          const passivedouble *Aj,
                                          const su2double     *B,
                                          su2double           *C);
