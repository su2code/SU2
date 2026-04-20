/*!
 * \file CLookUp_ANN_tests.cpp
 * \brief Unit tests for CLookUp_ANN and CIOMap classes.
 * \author E.C.Bunschoten
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "catch.hpp"
#include "../../../../Common/include/CConfig.hpp"
#if defined(HAVE_MLPCPP)
#define MLP_CUSTOM_TYPE su2double
#include "../../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif
#include <vector>

#ifdef USE_MLPCPP
TEST_CASE("LookUp ANN test", "[LookUpANN]") {
  MLPToolbox::CLookUp_ANN ANN;
  MLPToolbox::CNeuralNetwork mlp =
      MLPToolbox::CNeuralNetwork("src/SU2/UnitTests/Common/toolboxes/multilayer_perceptron/simple_mlp.mlp");

  ANN.AddNetwork(&mlp);
  su2double x, y, z, z_alt;

  /*--- Create a query where the value of z is calculated from x and y ---*/
  MLPToolbox::CIOMap iomap_ref;
  iomap_ref.AddQueryInput("x", &x);
  iomap_ref.AddQueryInput("y", &y);
  iomap_ref.AddQueryOutput("z", &z);

  ANN.PairVariableswithMLPs(iomap_ref);

  MLPToolbox::CIOMap iomap_vec;
  std::vector<std::string> input_names = {"x", "y"}, output_names = {"z"};
  std::vector<su2double*> output_refs = {&z_alt};
  std::vector<su2double> mlp_inputs;

  iomap_vec.SetQueryInput(input_names);
  iomap_vec.SetQueryOutput(output_names);
  ANN.PairVariableswithMLPs(iomap_vec);

  /*--- MLP evaluation on two points in the middle of the training data range ---*/
  x = 1.0;
  y = -0.5;

  bool inside = ANN.Predict(iomap_ref);
  CHECK(z == Approx(0.344829));
  CHECK(inside);

  mlp_inputs.resize(2);
  mlp_inputs[0] = x;
  mlp_inputs[1] = y;
  ANN.Predict(iomap_vec, mlp_inputs, output_refs);
  CHECK(z == z_alt);

  x = 0.5;
  y = -0.23;
  inside = ANN.Predict(iomap_ref);

  CHECK(z == Approx(0.224986));
  CHECK(inside);

  mlp_inputs[0] = x;
  mlp_inputs[1] = y;
  ANN.Predict(iomap_vec, mlp_inputs, output_refs);

  /*--- */
  x = 10.0;
  y = -20.0;
  inside = ANN.Predict(iomap_ref);
  CHECK(!inside);
}
#endif
