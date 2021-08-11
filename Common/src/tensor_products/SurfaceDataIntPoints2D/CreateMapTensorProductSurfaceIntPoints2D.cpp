/*!
 * \file CreateMapTensorProductSurfaceIntPoints2D.cpp
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

#include "../../../include/tensor_products/TensorProductSurfaceIntPoints2D.hpp"

void CreateMapTensorProductSurfaceIntPoints2D(map<CUnsignedShort2T, TPIS2D> &mapFunctions) {

  /*--- Make sure that the map is empty. ---*/
  mapFunctions.clear();

  /*--- Variable to store the number of DOFs and integration points as one entity. ---*/
  CUnsignedShort2T nDOFsAndInt;

  /*--- Insert the mappings from the CUnsignedShort2T to the function pointer. ---*/
  nDOFsAndInt.short0 = 1; nDOFsAndInt.short1 = 2;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_1_2);

  nDOFsAndInt.short0 = 1; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_1_3);

  nDOFsAndInt.short0 = 1; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_1_4);

  nDOFsAndInt.short0 = 1; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_1_5);

  nDOFsAndInt.short0 = 2; nDOFsAndInt.short1 = 2;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_2_2);

  nDOFsAndInt.short0 = 2; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_2_3);

  nDOFsAndInt.short0 = 2; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_2_4);

  nDOFsAndInt.short0 = 2; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_2_5);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 3;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_3_3);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_3_4);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_3_5);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_3_6);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_3_7);

  nDOFsAndInt.short0 = 3; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_3_8);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 4;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_4_4);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_4_5);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_4_6);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_4_7);

  nDOFsAndInt.short0 = 4; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_4_8);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 5;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_5_5);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_5_6);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_5_7);

  nDOFsAndInt.short0 = 5; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_5_8);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 6;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_6_6);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_6_7);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_6_8);

  nDOFsAndInt.short0 = 6; nDOFsAndInt.short1 = 9;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_6_9);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 7;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_7_7);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_7_8);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 9;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_7_9);

  nDOFsAndInt.short0 = 7; nDOFsAndInt.short1 = 12;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_7_12);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 8;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_8_8);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 12;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_8_12);

  nDOFsAndInt.short0 = 8; nDOFsAndInt.short1 = 13;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_8_13);

  nDOFsAndInt.short0 = 9; nDOFsAndInt.short1 = 9;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_9_9);

  nDOFsAndInt.short0 = 9; nDOFsAndInt.short1 = 13;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_9_13);

  nDOFsAndInt.short0 = 9; nDOFsAndInt.short1 = 14;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_9_14);

  nDOFsAndInt.short0 = 10; nDOFsAndInt.short1 = 10;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_10_10);

  nDOFsAndInt.short0 = 10; nDOFsAndInt.short1 = 14;
  mapFunctions.emplace(nDOFsAndInt, &TensorProductSurfaceIntPoints2D_10_14);
}
