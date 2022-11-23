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

namespace MLPToolbox{
class CNeuralNetwork
{
private:
    std::vector<std::string> input_names,   // MLP input variable names.
                             output_names;  // MLP output variable names.

    unsigned long n_hidden_layers;  // Number of hidden layers (layers between input and output layer).

    CLayer *inputLayer,     // Pointer to network input layer.
           *outputLayer;    // Pointer to network output layer.

    std::vector<CLayer*> hiddenLayers,  // Hidden layer collection.
                         total_layers;  // Hidden layers plus in/output layers

    std::vector<su2activematrix> weights_mat;   // Weights of synapses connecting layers

    std::vector<std::pair<su2double, su2double>> input_norm,    // Normalization factors for network inputs
                                                 output_norm;   // Normalization factors for network outputs

    std::vector<su2double> last_inputs; // Inputs from previous lookup operation. Evaluation of the network 
                                        // is skipped if current inputs are the same as the last inputs.

    su2double* ANN_outputs; // Pointer to network outputs

    /*--- Activation function types that are currently supported ---*/
    enum ENUM_ACTIVATION_FUNCTION {
        NONE=0,
        LINEAR=1,
        RELU=2,
        SMOOTH_SLOPE=3,
        ELU = 4,
        GELU = 5,
        SELU = 6,
        SIGMOID = 7,
        SWISH = 8,
        TANH = 9
    };
    ENUM_ACTIVATION_FUNCTION * activation_function_types;

public:
    CNeuralNetwork();
    ~CNeuralNetwork(){
        delete inputLayer;
        delete outputLayer;
        for(std::size_t i=1; i<total_layers.size()-1; i++){
            delete total_layers[i];
        }
        delete [] ANN_outputs;
        delete [] activation_function_types;
    };
    /*!
    * \brief Set the input layer of the network.
    * \param[in] n_neurons - Number of inputs
    */
    void defineInputLayer(unsigned long n_neurons);

    /*!
    * \brief Set the output layer of the network.
    * \param[in] n_neurons - Number of outputs
    */
    void defineOutputLayer(unsigned long n_neurons);

    /*!
    * \brief Add a hidden layer to the network
    * \param[in] n_neurons - Hidden layer size.
    */
    void push_hidden_layer(unsigned long n_neurons);

    /*!
    * \brief Set the weight value of a specific synapse.
    * \param[in] i_layer - current layer.
    * \param[in] i_neuron - neuron index in current layer.
    * \param[in] j_neuron - neuron index of connecting neuron.
    * \param[in] value - weight value.
    */
    void setWeight(unsigned long i_layer, unsigned long i_neuron, unsigned long j_neuron, su2double value);

    /*!
    * \brief Set bias value at a specific neuron.
    * \param[in] i_layer - Layer index.
    * \param[in] i_neuron - Neuron index of current layer.
    * \param[in] value - Bias value.
    */
    void setBias(unsigned long i_layer, unsigned long i_neuron, su2double value){total_layers.at(i_layer)->setBias(i_neuron, value);}

    /*!
    * \brief Set layer activation function.
    * \param[in] i_layer - Layer index.
    * \param[in] input - Activation function name.
    */
    void setActivationFunction(unsigned long i_layer, std::string input);

    /*!
    * \brief Display the network architecture in the terminal.
    */
    void displayNetwork();

    /*!
    * \brief Size the weight layers in the network according to its architecture.
    */
    void sizeWeights();

    /*!
    * \brief Size the vector of previous inputs.
    * \param[in] n_inputs - Number of inputs.
    */
    void sizeInputs(unsigned long n_inputs) { last_inputs.resize(n_inputs); for(unsigned long iInput=0; iInput<n_inputs; iInput++) last_inputs.at(iInput) = 0.0; }
    
    /*!
    * \brief Get the number of connecting regions in the network.
    * \returns number of spaces in between layers.
    */
    unsigned long getNWeightLayers(){return total_layers.size()-1;}

    /*!
    * \brief Get the total number of layers in the network
    * \returns number of netowork layers.
    */
    unsigned long getNLayers(){return total_layers.size();}

    /*!
    * \brief Get neuron count in a layer.
    * \param[in] iLayer - Layer index.
    * \returns number of neurons in the layer.
    */
    unsigned long getNNeurons(unsigned long iLayer) { return total_layers.at(iLayer)->getNNeurons(); }

    /*!
    * \brief Evaluate the network.
    */
    void predict(std::vector<su2double> &inputs);

    /*!
    * \brief Set the normalization factors for the input layer
    * \param[in] iInput - Input index.
    * \param[in] input_min - Minimum input value.
    * \param[in] input_max - Maximum input value.
    */
    void SetInputNorm(unsigned long iInput, su2double input_min, su2double input_max) { input_norm.at(iInput) = make_pair(input_min, input_max); }

    /*!
    * \brief Set the normalization factors for the output layer
    * \param[in] iOutput - Input index.
    * \param[in] input_min - Minimum output value.
    * \param[in] input_max - Maximum output value.
    */
    void SetOutputNorm(unsigned long iOutput, su2double output_min, su2double output_max){ output_norm.at(iOutput) = make_pair(output_min, output_max); }
    
    /*!
    * \brief Add an output variable name to the network.
    * \param[in] input - Input variable name.
    */
    void PushOutputName(std::string input) { output_names.push_back(input); }

    /*!
    * \brief Add an input variable name to the network.
    * \param[in] output - Output variable name.
    */
    void PushInputName(std::string input) { input_names.push_back(input); }
   
    /*!
    * \brief Get network input variable name.
    * \param[in] iInput - Input variable index.
    * \returns input variable name.
    */
    std::string GetInputName(std::size_t iInput) { return input_names[iInput]; }

    /*!
    * \brief Get network output variable name.
    * \param[in] iOutput - Output variable index.
    * \returns output variable name.
    */
    std::string GetOutputName(std::size_t iOutput) { return output_names[iOutput]; }

    /*!
    * \brief Get network number of inputs.
    * \returns Number of network inputs
    */
    std::size_t GetnInputs() { return input_names.size(); }

    /*!
    * \brief Get network number of outputs.
    * \returns Number of network outputs
    */
    std::size_t GetnOutputs() { return output_names.size(); }

    /*!
    * \brief Get network evaluation output.
    * \param[in] iOutput - output index.
    * \returns Prediction value.
    */
    su2double GetANN_Output(std::size_t iOutput) { return ANN_outputs[iOutput]; }

    /*!
    * \brief Set the activation function array size.
    * \param[in] n_layers - network layer count.
    */
    void SizeActivationFunctions(unsigned long n_layers) { activation_function_types = new ENUM_ACTIVATION_FUNCTION[n_layers]; }

    /*!
    * \brief Compute neuron activation function input.
    * \param[in] iLayer - Network layer index.
    * \param[in] iNeuron - Layer neuron index.
    * \returns Neuron activation function input.
    */
    su2double ComputeX(std::size_t iLayer, std::size_t iNeuron){
        su2double x;
        x = total_layers[iLayer]->getBias(iNeuron);
        std::size_t nNeurons_previous = total_layers[iLayer - 1]->getNNeurons();
        for(std::size_t jNeuron=0; jNeuron<nNeurons_previous; jNeuron++){
            x += weights_mat[iLayer - 1][iNeuron][jNeuron] * total_layers[iLayer-1]->getOutput(jNeuron);
        }
        return x;
    }
};


}


