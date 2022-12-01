/*!
 * \file CReadNeuralNetwork.hpp
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
#include <vector>
#include "../../CConfig.hpp"

namespace MLPToolbox{

class CReadNeuralNetwork {
    private:
    std::vector<std::string> input_names,   // Network input variable names.
                             output_names;  // Network output variable names.

    std::string filename;       // Network input filename.

    unsigned long n_layers;     // Number of layers in the network.
    
    std::vector<unsigned long> n_neurons;
    std::vector<std::vector<std::vector<su2double>>> weights;
    std::vector<std::vector<su2double>> biases;
    std::vector<std::string> activation_functions;

    std::vector<std::pair<su2double, su2double>> input_norm;
    std::vector<std::pair<su2double, su2double>> output_norm;
    public:

    CReadNeuralNetwork(std::string filename_in);
    void ReadMLPFile();

    std::string SkipToFlag(std::ifstream *file_stream, std::string flag);
    
    unsigned long GetNInputs(){return n_neurons.at(0);}
    unsigned long GetNOutputs(){return n_neurons.at(n_layers - 1);}

    unsigned long GetNlayers(){return n_layers;}
    unsigned long GetNneurons(std::size_t iLayer){return n_neurons.at(iLayer);}
    su2double GetWeight(std::size_t iLayer, std::size_t iNeuron, std::size_t jNeuron){return weights.at(iLayer).at(iNeuron).at(jNeuron);}
    su2double GetBias(std::size_t iLayer, std::size_t iNeuron){return biases.at(iLayer).at(iNeuron);}
    std::pair<su2double, su2double> GetInputNorm(std::size_t iInput){return input_norm.at(iInput);}
    std::pair<su2double, su2double> GetOutputNorm(std::size_t iOutput){return output_norm.at(iOutput);}
    std::string GetActivationFunction(std::size_t iLayer){return activation_functions.at(iLayer);}

    std::string GetInputName(std::size_t iInput){return input_names.at(iInput);}
    std::string GetOutputName(std::size_t iOutput){return output_names.at(iOutput);}
    
    ~CReadNeuralNetwork(){};
};
}
