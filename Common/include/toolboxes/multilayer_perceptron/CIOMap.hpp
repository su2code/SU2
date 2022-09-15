/*!
 * \file CIOMap.hpp
 * \brief Declaration of input-to-output mapping class for artifical neural networks
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
#pragma once

#include <cmath>
#include <string>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <vector>
#include "../../CConfig.hpp"
#include "../../linear_algebra/blas_structure.hpp"
#include "CLookUp_ANN.hpp"

using namespace std;
namespace MLPToolbox
{
        
    class CLookUp_ANN;

    class CIOMap
    {
        private:
        vector<size_t> ANN_indices;
        vector<vector<pair<size_t, size_t>>> Input_Map;
        vector<vector<pair<size_t, size_t>>> Output_Map;
        vector<vector<su2double>> input_indices;
        public:
        CIOMap(){};
        CIOMap(CLookUp_ANN*ANN_collection, vector<string> &inputs, vector<string> &outputs);
        ~CIOMap(){};
        void Push_ANN_index(size_t i_ANN){ANN_indices.push_back(i_ANN);};
        void PairVariableswithANNs(CLookUp_ANN * ANN_collection, vector<string> &inputs, vector<string> &outputs);
        size_t GetNANNs(){return ANN_indices.size();}
        size_t GetANNIndex(size_t i_Map){return ANN_indices[i_Map];}
        size_t GetInputIndex(size_t i_Map, size_t iInput){return Input_Map[i_Map][iInput].first;}
        size_t GetOutputIndex(size_t i_Map, size_t iOutput){return Output_Map[i_Map][iOutput].first;}
        size_t GetANNOutputIndex(size_t i_Map, size_t iOutput){return Output_Map[i_Map][iOutput].second;}
        size_t GetNMappedOutputs(size_t i_Map){return Output_Map[i_Map].size();}

        vector<pair<size_t, size_t>> GetOutputMapping(size_t i_map){return Output_Map[i_map];}
        vector<pair<size_t, size_t>> GetInputMapping(size_t i_map){return Input_Map[i_map];}
        
        vector<su2double> GetANN_Inputs(size_t i_Map, vector<su2double>&inputs){
            vector<su2double> ANN_input;
            ANN_input.resize(Input_Map[i_Map].size());

            for(size_t iInput=0; iInput<Input_Map[i_Map].size(); iInput++){
                ANN_input.at(iInput) = inputs[GetInputIndex(i_Map, iInput)];
            }
            return ANN_input;
        }
        vector<su2double*> GetANN_Outputs(size_t i_Map, vector<su2double*>&outputs){
            vector<su2double*> ANN_output;
            ANN_output.resize(Output_Map[i_Map].size());

            for(size_t iOutput=0; iOutput<Output_Map[i_Map].size(); iOutput++){
                ANN_output.at(iOutput) = outputs[GetOutputIndex(i_Map, iOutput)];
            }
            return ANN_output;
        }
    };
}
