/*!
 * \file CIOMap.cpp
 * \brief Implementation of the input-output mapping class for the
 *      use of multi-layer perceptrons in SU2
 * \author E. Bunschoten
 * \version 7.5.1 "Blackbird"
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
#include "../../../include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"
#include "../../../include/toolboxes/multilayer_perceptron/CIOMap.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

#include "../../../include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"

MLPToolbox::CIOMap::CIOMap(CLookUp_ANN* MLP_collection, su2vector<std::string>& inputs,
                           su2vector<std::string>& outputs) {

    /* Generate an input-output map given a set of call inputs and call outputs. These inputs 
  mapped to the inputs of the loaded MLPs in the MLP_collection object. The call outputs are
  then matched to the MLPs with matching inputs.
  */
  PairVariableswithMLPs(MLP_collection, inputs, outputs);

  // Perform checks on input-output validity
  if (outputs.size() > 0) {
    // Check wether all call inputs are used
    MLP_collection->CheckUseOfInputs(inputs, this);

    // Check wether all call outputs are present in the MLP collection
    MLP_collection->CheckUseOfOutputs(outputs, this);
  }
}
void MLPToolbox::CIOMap::PairVariableswithMLPs(CLookUp_ANN* MLP_collection, su2vector<std::string>& inputs,
                                               su2vector<std::string>& outputs) {
  /*
  In this function, the call inputs and outputs are matched to those within the MLP collection.
  */
  bool isInput, isOutput;

  // Looping over the loaded MLPs to check wether the MLP inputs match with the call inputs
  for (size_t iMLP = 0; iMLP < MLP_collection->GetNANNs(); iMLP++) {
    // Mapped call inputs to MLP inputs
    std::vector<pair<size_t, size_t>> Input_Indices = MLP_collection->FindVariableIndices(iMLP, inputs, true);
    isInput = Input_Indices.size() > 0;

    if (isInput) {
      // Only when the MLP inputs match with a portion of the call inputs are the output variable checks performed

      std::vector<pair<size_t, size_t>> Output_Indices = MLP_collection->FindVariableIndices(iMLP, outputs, false);
      isOutput = Output_Indices.size() > 0;

      if (isOutput) {
        // Update input and output mapping if both inputs and outputs match
        MLP_indices.push_back(iMLP);
        Input_Map.push_back(Input_Indices);
        Output_Map.push_back(Output_Indices);
      }
    }
  }
}