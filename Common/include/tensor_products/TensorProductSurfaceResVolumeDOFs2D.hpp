/*!
 * \file TensorProductSurfaceResVolumeDOFs2D.hpp
 * \brief Function prototypes for the tensor product to compute the contribution to the residual of the surface data in the 1D integration points of a line
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

typedef void(*TPDR2D)(const int           N,
                      const int           faceID,
                      const int           ldb,
                      const int           ldc,
                      const passivedouble *An,
                      const passivedouble *ATt,
                      const su2double     *B,
                      su2double           *C);
/*!
 * \brief Function, which stores the available function pointers for the tensor
 *        product for the 1D line integration points adjacent to a quad in a map.
 * \param[out] mapFunctions - Map to store the function pointers to carry out the tensor product.
 */
void CreateMapTensorProductSurfaceResVolumeDOFs2D(map<CUnsignedShort2T, TPDR2D> &mapFunctions);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (2,1).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_2_1(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (3,1).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_3_1(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (4,1).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_4_1(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (5,1).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_5_1(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (2,2).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_2_2(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (3,2).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_3_2(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (4,2).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_4_2(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (5,2).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_5_2(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (3,3).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_3_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (4,3).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_4_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (5,3).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_5_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (6,3).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_6_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (7,3).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_7_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (8,3).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_8_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (4,4).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_4_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (5,4).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_5_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (6,4).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_6_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (7,4).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_7_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (8,4).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_8_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (5,5).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_5_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (6,5).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_6_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (7,5).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_7_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (8,5).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_8_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (6,6).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_6_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (7,6).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_7_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (8,6).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_8_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (9,6).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_9_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (7,7).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_7_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (8,7).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_8_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (9,7).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_9_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (12,7).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_12_7(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const passivedouble *An,
                                              const passivedouble *ATt,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (8,8).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_8_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (12,8).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_12_8(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const passivedouble *An,
                                              const passivedouble *ATt,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (13,8).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_13_8(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const passivedouble *An,
                                              const passivedouble *ATt,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (9,9).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_9_9(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const passivedouble *An,
                                             const passivedouble *ATt,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (13,9).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_13_9(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const passivedouble *An,
                                              const passivedouble *ATt,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (14,9).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_14_9(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const passivedouble *An,
                                              const passivedouble *ATt,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (10,10).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_10_10(const int           N,
                                               const int           faceID,
                                               const int           ldb,
                                               const int           ldc,
                                               const passivedouble *An,
                                               const passivedouble *ATt,
                                               const su2double     *B,
                                               su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 2D DOFs of a quad from the data in the 1D integration points
 *        of a line adjacent to the quad for (nInt1D,nDOFs1D) = (14,10).
 * \param[in]  N      - Number of variables to be determined in the integration points
 * \param[in]  faceID - Face ID of the quad on which the line data are defined
 * \param[in]  ldb    - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc    - Leading dimension of C when stored as a matrix.
 * \param[in]  An     - Component of the A tensor normal to the line.
 * \param[in]  ATt    - Component of the A transpose tensor tangential to the line.
 * \param[in]  B      - Tensor, which contains the residual in the integration points of the line.
 * \param[out] C      - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs2D_14_10(const int           N,
                                               const int           faceID,
                                               const int           ldb,
                                               const int           ldc,
                                               const passivedouble *An,
                                               const passivedouble *ATt,
                                               const su2double     *B,
                                               su2double           *C);
