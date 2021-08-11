/*!
 * \file TensorProductSurfaceResVolumeDOFs3D.hpp
 * \brief Function prototypes for the tensor product to compute the contribution to the residual of the surface data in the 2D integration points of a quad
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

typedef void(*TPDR3D)(const int           N,
                      const int           faceID,
                      const int           ldb,
                      const int           ldc,
                      const bool          swapTangDir,
                      const passivedouble *An,
                      const passivedouble *ATt0,
                      const passivedouble *ATt1,
                      const su2double     *B,
                      su2double           *C);
/*!
 * \brief Function, which stores the available function pointers for the tensor
 *        product for the 2D quad integration points adjacent to a hex in a map.
 * \param[out] mapFunctions - Map to store the function pointers to carry out the tensor product.
 */
void CreateMapTensorProductSurfaceResVolumeDOFs3D(map<CUnsignedShort2T, TPDR3D> &mapFunctions);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (1,2).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_1_2(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (1,3).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_1_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (1,4).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_1_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (1,5).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_1_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (2,2).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_2_2(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (2,3).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_2_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (2,4).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_2_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (2,5).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_2_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (3,3).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_3_3(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (3,4).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_3_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (3,5).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_3_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (3,6).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_3_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (3,7).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_3_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (3,8).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_3_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (4,4).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_4_4(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (4,5).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_4_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (4,6).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_4_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (4,7).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_4_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (4,8).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_4_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (5,5).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_5_5(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (5,6).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_5_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (5,7).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_5_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (5,8).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_5_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (6,6).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_6_6(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (6,7).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_6_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (6,8).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_6_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (6,9).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_6_9(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (7,7).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_7_7(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (7,8).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_7_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (7,9).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_7_9(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (7,12).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_7_12(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const bool          swapTanDir,
                                              const passivedouble *An,
                                              const passivedouble *ATt0,
                                              const passivedouble *ATt1,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (8,8).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_8_8(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (8,12).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_8_12(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const bool          swapTanDir,
                                              const passivedouble *An,
                                              const passivedouble *ATt0,
                                              const passivedouble *ATt1,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (8,13).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_8_13(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const bool          swapTanDir,
                                              const passivedouble *An,
                                              const passivedouble *ATt0,
                                              const passivedouble *ATt1,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (9,9).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_9_9(const int           N,
                                             const int           faceID,
                                             const int           ldb,
                                             const int           ldc,
                                             const bool          swapTanDir,
                                             const passivedouble *An,
                                             const passivedouble *ATt0,
                                             const passivedouble *ATt1,
                                             const su2double     *B,
                                             su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (9,13).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_9_13(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const bool          swapTanDir,
                                              const passivedouble *An,
                                              const passivedouble *ATt0,
                                              const passivedouble *ATt1,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (9,14).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_9_14(const int           N,
                                              const int           faceID,
                                              const int           ldb,
                                              const int           ldc,
                                              const bool          swapTanDir,
                                              const passivedouble *An,
                                              const passivedouble *ATt0,
                                              const passivedouble *ATt1,
                                              const su2double     *B,
                                              su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (10,10).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_10_10(const int           N,
                                               const int           faceID,
                                               const int           ldb,
                                               const int           ldc,
                                               const bool          swapTanDir,
                                               const passivedouble *An,
                                               const passivedouble *ATt0,
                                               const passivedouble *ATt1,
                                               const su2double     *B,
                                               su2double           *C);

/*!
 * \brief Function, which carries out the tensor product to update the residual
 *        in the 3D DOFs of a hex from the data in the 2D integration points
 *        of a quad adjacent to the hex for (nDOFs1D,nInt1D) = (10,14).
 * \param[in]  N          - Number of variables to be determined in the integration points
 * \param[in]  faceID     - Face ID of the hex on which the quad data are defined
 * \param[in]  ldb        - Leading dimension of B when stored as a matrix.
 * \param[in]  ldc        - Leading dimension of C when stored as a matrix.
 * \param[in]  swapTanDir - Whether or not to swap the tangential directions of the quad.
 * \param[in]  An         - Component of the A tensor normal to the face.
 * \param[in]  ATt0       - Component of the A transpose tensor in the first tangential direction.
 * \param[in]  ATt1       - Component of the A transpose tensor in the second tangential direction.
 * \param[in]  B          - Tensor, which contains the residual in the integration points of the quad.
 * \param[out] C          - Result of the tensor product C = A*B.
 */
void TensorProductSurfaceResVolumeDOFs3D_10_14(const int           N,
                                               const int           faceID,
                                               const int           ldb,
                                               const int           ldc,
                                               const bool          swapTanDir,
                                               const passivedouble *An,
                                               const passivedouble *ATt0,
                                               const passivedouble *ATt1,
                                               const su2double     *B,
                                               su2double           *C);
