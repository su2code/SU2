/*!
 * \file CLookUp_ANN_tests.cpp
 * \brief Unit tests for CLookUp_ANN and CIOMap classes.
 * \author E.C.Bunschoten
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif
#include <vector>

#ifdef USE_MLPCPP
TEST_CASE("LookUp ANN test", "[LookUpANN]") {
  std::string MLP_input_files[] = {"src/SU2/UnitTests/Common/toolboxes/multilayer_perceptron/simple_mlp.mlp"};
  unsigned short n_MLPs = 1;
  MLPToolbox::CLookUp_ANN ANN(n_MLPs, MLP_input_files);
  std::vector<std::string> MLP_input_names, MLP_output_names;
  std::vector<double> MLP_inputs;
  std::vector<double*> MLP_outputs;
  su2double x, y, z;

  /*--- Define MLP inputs and outputs ---*/
  MLP_input_names.resize(2);
  MLP_input_names[0] = "x";
  MLP_input_names[1] = "y";
  MLP_inputs.resize(2);

  MLP_outputs.resize(1);
  MLP_output_names.resize(1);
  MLP_output_names[0] = "z";
  MLP_outputs[0] = &z;

  /*--- Generate input-output map ---*/
  MLPToolbox::CIOMap iomap(MLP_input_names, MLP_output_names);
  ANN.PairVariableswithMLPs(iomap);
  /*--- MLP evaluation on point in the middle of the training data range ---*/
  x = 1.0;
  y = -0.5;

  MLP_inputs[0] = x;
  MLP_inputs[1] = y;
  ANN.PredictANN(&iomap, MLP_inputs, MLP_outputs);
  CHECK(z == Approx(0.344829));

  /*--- MLP evaluation on point outside the training data range ---*/
  x = 3.0;
  y = -10;
  MLP_inputs[0] = x;
  MLP_inputs[1] = y;
  ANN.PredictANN(&iomap, MLP_inputs, MLP_outputs);
  CHECK(z == Approx(0.012737));
}
#endif
