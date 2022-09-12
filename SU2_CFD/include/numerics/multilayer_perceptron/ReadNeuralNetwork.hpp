/*!
 * \file ReadNeuralNetwork.hpp
 * \brief Declaration of MLP input file reader class
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
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/linear_algebra/blas_structure.hpp"

using namespace std;
class ReadNeuralNetwork {
    private:
    vector<string> input_names;
    vector<string> output_names;

    string filename;
    unsigned long n_layers;
    vector<unsigned long> n_neurons;
    vector<vector<vector<double long>>> weights;
    vector<vector<double long>> biases;
    vector<string> activation_functions;

    vector<pair<double long, double long>> input_norm;
    vector<pair<double long, double long>> output_norm;
    public:

    ReadNeuralNetwork(string filename_in);
    void ReadMLPFile();

    string SkipToFlag(ifstream *file_stream, string flag);
    
    unsigned long GetNInputs(){return n_neurons.at(0);}
    unsigned long GetNOutputs(){return n_neurons.at(n_layers - 1);}

    unsigned long GetNlayers(){return n_layers;}
    unsigned long GetNneurons(size_t iLayer){return n_neurons.at(iLayer);}
    double long GetWeight(size_t iLayer, size_t iNeuron, size_t jNeuron){return weights.at(iLayer).at(iNeuron).at(jNeuron);}
    double long GetBias(size_t iLayer, size_t iNeuron){return biases.at(iLayer).at(iNeuron);}
    pair<double long, double long> GetInputNorm(size_t iInput){return input_norm.at(iInput);}
    pair<double long, double long> GetOutputNorm(size_t iOutput){return output_norm.at(iOutput);}
    string GetActivationFunction(size_t iLayer){return activation_functions.at(iLayer);}

    string GetInputName(size_t iInput){return input_names.at(iInput);}
    string GetOutputName(size_t iOutput){return output_names.at(iOutput);}
    
    ~ReadNeuralNetwork(){};
};