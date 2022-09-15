/*!
 * \file CIOMap.cpp
 * \brief Implementation of the input-output mapping class for the 
 *      use of multi-layer perceptrons in SU2
 * \author E. Bunschoten
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

MLPToolbox::CIOMap::CIOMap(CLookUp_ANN*ANN_collection, vector<string> &inputs, vector<string> &outputs){
    PairVariableswithANNs(ANN_collection, inputs, outputs);
    if(outputs.size() > 0){
        ANN_collection->Check_Use_of_Inputs(inputs, this);
        ANN_collection->Check_Use_of_Outputs(outputs, this);
        ANN_collection->Check_Duplicate_Outputs(outputs, this);
    }
}
void MLPToolbox::CIOMap::PairVariableswithANNs(CLookUp_ANN*ANN_collection, vector<string> &inputs, vector<string> &outputs){

    bool isInput, isOutput;
    for(size_t i_ANN=0; i_ANN<ANN_collection->GetNANNs(); i_ANN++){
        vector<pair<size_t, size_t>> Input_Indices = ANN_collection->FindVariable_Indices(i_ANN, inputs, true);
        isInput = Input_Indices.size() > 0;
        if(isInput){
            vector<pair<size_t, size_t>> Output_Indices = ANN_collection->FindVariable_Indices(i_ANN, outputs, false);
            isOutput = Output_Indices.size() > 0;
            if(isOutput){
                ANN_indices.push_back(i_ANN);
                Input_Map.push_back(Input_Indices);
                Output_Map.push_back(Output_Indices);
            }
        }
    }
}