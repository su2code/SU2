/*!
 * \file CreateMapTensorProductSurfaceResVolumeDOFs2D
 * \brief Function, which creates the map between the number of 1D DOFs and integration points and the function pointers.
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

void CreateMapTensorProductSurfaceResVolumeDOFs2D(map<CUnsignedShort2T, TPDR2D> &mapFunctions) {

  /*--- Make sure that the map is empty. ---*/
  mapFunctions.clear();

  /*--- Variable to store the number of DOFs and integration points as one entity. ---*/
  CUnsignedShort2T nDOFsAndInt;

  /*--- Insert the mappings from the CUnsignedShort2T to the function pointer. ---*/
  nDOFsAndInt.short0 = 2; nDOFsAndInt.short1 = 1;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_2_1);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 1;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_3_1);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 1;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_4_1);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 1;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_5_1);

  nDOFsAndInt.short0 = 2; nDOFsAndInt.short1 = 2;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_2_2);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 2;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_3_2);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 2;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_4_2);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 2;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_5_2);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_3_3);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_4_3);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_5_3);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_6_3);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_7_3);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_8_3);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_4_4);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_5_4);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_6_4);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_7_4);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_8_4);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_5_5);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_6_5);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_7_5);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_8_5);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_6_6);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_7_6);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_8_6);

  nDOFsAndInt.short0 = 9; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_9_6);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_7_7);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_8_7);

  nDOFsAndInt.short0 = 9; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_9_7);

  nDOFsAndInt.short0 = 12; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_12_7);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_8_8);

  nDOFsAndInt.short0 = 12; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_12_8);

  nDOFsAndInt.short0 = 13; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_13_8);

  nDOFsAndInt.short0 = 9; nDOFsAndInt.short1 = 9;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_9_9);

  nDOFsAndInt.short0 = 13; nDOFsAndInt.short1 = 9;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_13_9);

  nDOFsAndInt.short0 = 14; nDOFsAndInt.short1 = 9;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_14_9);

  nDOFsAndInt.short0 = 10; nDOFsAndInt.short1 = 10;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_10_10);

  nDOFsAndInt.short0 = 14; nDOFsAndInt.short1 = 10;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceResVolumeDOFs2D_14_10);
}
