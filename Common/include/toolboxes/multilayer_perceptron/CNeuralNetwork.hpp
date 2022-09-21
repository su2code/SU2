/*!
 * \file CNeuralNetwork.hpp
 * \brief Declaration of the neural network class
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
#include <cstring>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../CConfig.hpp"
#include "../../linear_algebra/blas_structure.hpp"
#include "CLayer.hpp"

using namespace std;
namespace MLPToolbox{
class CNeuralNetwork
{
private:
    vector<string> input_names,
                   output_names;

    unsigned long n_hidden_layers;

    CLayer *inputLayer,
           *outputLayer;

    vector<CLayer*> hiddenLayers,
                    total_layers;

    vector<vector<vector<su2double>>> weights;
    vector<su2activematrix> weights_mat;
    vector<vector<su2double>> values;

    vector<pair<su2double, su2double>> input_norm,
                                       output_norm;
    vector<su2double> last_inputs;

    su2double* ANN_outputs;

    enum ENUM_ACTIVATION_FUNCTION {
        NONE=0,
        LINEAR=1,
        RELU=2,
        SMOOTH_SLOPE=3,
        ELU = 4,
        SWISH = 5,
        SIGMOID = 6,
        TANH = 7,
        SELU = 8,
        GELU = 9
    };
    ENUM_ACTIVATION_FUNCTION * activation_function_types;

public:
    CNeuralNetwork();
    ~CNeuralNetwork(){
        for(size_t i=0; i<total_layers.size(); i++){
            delete total_layers.at(i);
        }
        delete [] ANN_outputs;
    };
    void defineInputLayer(unsigned long n_neurons);
    void defineOutputLayer(unsigned long n_neurons);
    void push_hidden_layer(unsigned long n_neurons);
    void setWeight(unsigned long i_layer, unsigned long i_neuron, unsigned long j_neuron, su2double value);
    void setBias(unsigned long i_layer, unsigned long i_neuron, su2double value){total_layers.at(i_layer)->setBias(i_neuron, value);}
    void setActivationFunction(unsigned long i_layer, string input);
    void displayNetwork();
    void sizeWeights();
    void sizeInputs(unsigned long n_inputs){last_inputs.resize(n_inputs); for(unsigned long iInput=0; iInput<n_inputs; iInput++) last_inputs.at(iInput) = 0.0;}
    unsigned long getNWeightLayers(){return weights.size();}
    unsigned long getNLayers(){return total_layers.size();}

    unsigned long getNNeurons(unsigned long iLayer){;
    return total_layers.at(iLayer)->getNNeurons();}
    unsigned long getNNeurons(unsigned long iLayer, unsigned long iNeuron){return weights.at(iLayer).at(iNeuron).size();}

    void predict(vector<su2double> &inputs);


    void SetInputNorm(unsigned long iInput, su2double input_min, su2double input_max){input_norm.at(iInput) = make_pair(input_min, input_max);}
    void SetOutputNorm(unsigned long iOutput, su2double output_min, su2double output_max){output_norm.at(iOutput) = make_pair(output_min, output_max);}
    
    void PushOutputName(string input){output_names.push_back(input);}
    void PushInputName(string input){input_names.push_back(input);}
   
    string GetInputName(size_t iInput){return input_names[iInput];}
    string GetOutputName(size_t iOutput){return output_names[iOutput];}

    size_t GetnInputs(){return input_names.size();}
    size_t GetnOutputs(){return output_names.size();}

    su2double GetANN_Output(size_t iOutput){return ANN_outputs[iOutput];}
    void SizeActivationFunctions(unsigned long n_layers){activation_function_types = new ENUM_ACTIVATION_FUNCTION[n_layers]; 
    }

    su2double ComputeX(size_t iLayer, size_t iNeuron){
        su2double x;
        x = total_layers[iLayer]->getBias(iNeuron);
        size_t nNeurons_previous = total_layers[iLayer - 1]->getNNeurons();
        for(size_t jNeuron=0; jNeuron<nNeurons_previous; jNeuron++){
            x += weights_mat[iLayer - 1][iNeuron][jNeuron] * total_layers[iLayer-1]->getOutput(jNeuron);
        }
        return x;
    }
};


}


