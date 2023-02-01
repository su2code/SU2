/*!
 * \file CNeuralNetwork.hpp
 * \brief Declaration of the neural network class
 * \author E.C.Bunschoten
 * \version 7.5.0 "Blackbird"
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
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

#include "../../CConfig.hpp"
#include "../../linear_algebra/blas_structure.hpp"
#include "CLayer.hpp"

namespace MLPToolbox {
class CNeuralNetwork {
 private:
  su2vector<std::string> input_names,  // MLP input variable names.
      output_names;                    // MLP output variable names.

  unsigned long n_hidden_layers = 0;  // Number of hidden layers (layers between input and output layer).

  CLayer *inputLayer = nullptr,  // Pointer to network input layer.
      *outputLayer = nullptr;    // Pointer to network output layer.

  std::vector<CLayer*> hiddenLayers;  // Hidden layer collection.
  su2vector<CLayer*> total_layers;    // Hidden layers plus in/output layers

  su2vector<su2activematrix> weights_mat;  // Weights of synapses connecting layers

  su2vector<std::pair<su2double, su2double>> input_norm,  // Normalization factors for network inputs
      output_norm;                                        // Normalization factors for network outputs

  su2vector<su2double> last_inputs;  // Inputs from previous lookup operation. Evaluation of the network
                                     // is skipped if current inputs are the same as the last inputs.

  su2double* ANN_outputs;                 // Pointer to network outputs
  su2matrix<su2double> dOutputs_dInputs;  // Network output derivatives w.r.t inputs

  /*--- Activation function types that are currently supported ---*/
  enum ENUM_ACTIVATION_FUNCTION {
    NONE = 0,
    LINEAR = 1,
    RELU = 2,
    ELU = 3,
    GELU = 4,
    SELU = 5,
    SIGMOID = 6,
    SWISH = 7,
    TANH = 8,
    EXPONENTIAL = 9
  };
  su2vector<ENUM_ACTIVATION_FUNCTION> activation_function_types;
  su2vector<string> activation_function_names;
 public:
  ~CNeuralNetwork() {
    delete inputLayer;
    delete outputLayer;
    for (std::size_t i = 1; i < total_layers.size() - 1; i++) {
      delete total_layers[i];
    }
    delete[] ANN_outputs;
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
  void setWeight(unsigned long i_layer, unsigned long i_neuron, unsigned long j_neuron, su2double value) {
    weights_mat[i_layer][j_neuron][i_neuron] = value;
  };

  /*!
   * \brief Set bias value at a specific neuron.
   * \param[in] i_layer - Layer index.
   * \param[in] i_neuron - Neuron index of current layer.
   * \param[in] value - Bias value.
   */
  void setBias(unsigned long i_layer, unsigned long i_neuron, su2double value) {
    total_layers[i_layer]->setBias(i_neuron, value);
  }

  /*!
   * \brief Set layer activation function.
   * \param[in] i_layer - Layer index.
   * \param[in] input - Activation function name.
   */
  void setActivationFunction(unsigned long i_layer, std::string input);

  /*!
   * \brief Display the network architecture in the terminal.
   */
  void displayNetwork() const;

  /*!
   * \brief Size the weight layers in the network according to its architecture.
   */
  void sizeWeights();

  /*!
   * \brief Size the vector of previous inputs.
   * \param[in] n_inputs - Number of inputs.
   */
  void sizeInputs(unsigned long n_inputs) {
    last_inputs.resize(n_inputs);
    for (unsigned long iInput = 0; iInput < n_inputs; iInput++) last_inputs[iInput] = 0.0;
  }

  /*!
   * \brief Get the number of connecting regions in the network.
   * \returns number of spaces in between layers.
   */
  unsigned long getNWeightLayers() const { return total_layers.size() - 1; }

  /*!
   * \brief Get the total number of layers in the network
   * \returns number of netowork layers.
   */
  unsigned long getNLayers() const { return total_layers.size(); }

  /*!
   * \brief Get neuron count in a layer.
   * \param[in] iLayer - Layer index.
   * \returns number of neurons in the layer.
   */
  unsigned long getNNeurons(unsigned long iLayer) const { return total_layers[iLayer]->getNNeurons(); }

  /*!
   * \brief Evaluate the network.
   * \param[in] inputs - Network input variable values.
   * \param[in] compute_gradient - Compute the derivatives of the outputs wrt inputs.
   */
  void predict(su2vector<su2double>& inputs, bool compute_gradient = false);

  /*!
   * \brief Set the normalization factors for the input layer
   * \param[in] iInput - Input index.
   * \param[in] input_min - Minimum input value.
   * \param[in] input_max - Maximum input value.
   */
  void SetInputNorm(unsigned long iInput, su2double input_min, su2double input_max) {
    input_norm[iInput] = make_pair(input_min, input_max);
  }

  /*!
   * \brief Set the normalization factors for the output layer
   * \param[in] iOutput - Input index.
   * \param[in] input_min - Minimum output value.
   * \param[in] input_max - Maximum output value.
   */
  void SetOutputNorm(unsigned long iOutput, su2double output_min, su2double output_max) {
    output_norm[iOutput] = make_pair(output_min, output_max);
  }

  std::pair<su2double, su2double> GetInputNorm(unsigned long iInput) const { return input_norm[iInput]; }

  std::pair<su2double, su2double> GetOutputNorm(unsigned long iOutput) const { return output_norm[iOutput]; }
  /*!
   * \brief Add an output variable name to the network.
   * \param[in] input - Input variable name.
   */
  void SetOutputName(size_t iOutput, std::string input) { output_names[iOutput] = input; }

  /*!
   * \brief Add an input variable name to the network.
   * \param[in] output - Output variable name.
   */
  void SetInputName(size_t iInput, std::string input) { input_names[iInput] = input; }

  /*!
   * \brief Get network input variable name.
   * \param[in] iInput - Input variable index.
   * \returns input variable name.
   */
  std::string GetInputName(std::size_t iInput) const { return input_names[iInput]; }

  /*!
   * \brief Get network output variable name.
   * \param[in] iOutput - Output variable index.
   * \returns output variable name.
   */
  std::string GetOutputName(std::size_t iOutput) const { return output_names[iOutput]; }

  /*!
   * \brief Get network number of inputs.
   * \returns Number of network inputs
   */
  std::size_t GetnInputs() const { return input_names.size(); }

  /*!
   * \brief Get network number of outputs.
   * \returns Number of network outputs
   */
  std::size_t GetnOutputs() const { return output_names.size(); }

  /*!
   * \brief Get network evaluation output.
   * \param[in] iOutput - output index.
   * \returns Prediction value.
   */
  su2double GetANN_Output(std::size_t iOutput) const { return ANN_outputs[iOutput]; }

  /*!
   * \brief Get network output derivative w.r.t specific input.
   * \param[in] iOutput - output variable index.
   * \param[in] iInput - input variable index.
   * \returns Output derivative w.r.t input.
   */
  su2double GetANN_Output_Input_Derivative(std::size_t iOutput, std::size_t iInput) const {
    return dOutputs_dInputs[iOutput][iInput];
  }

  /*!
   * \brief Set the activation function array size.
   * \param[in] n_layers - network layer count.
   */
  void SizeActivationFunctions(unsigned long n_layers) {
    activation_function_types.resize(n_layers);
    activation_function_names.resize(n_layers);
  }

  /*!
   * \brief Compute neuron activation function input.
   * \param[in] iLayer - Network layer index.
   * \param[in] iNeuron - Layer neuron index.
   * \returns Neuron activation function input.
   */
  su2double ComputeX(std::size_t iLayer, std::size_t iNeuron) const {
    su2double x;
    x = total_layers[iLayer]->getBias(iNeuron);
    std::size_t nNeurons_previous = total_layers[iLayer - 1]->getNNeurons();
    for (std::size_t jNeuron = 0; jNeuron < nNeurons_previous; jNeuron++) {
      x += weights_mat[iLayer - 1][iNeuron][jNeuron] * total_layers[iLayer - 1]->getOutput(jNeuron);
    }
    return x;
  }
  su2double ComputedOutputdInput(std::size_t iLayer, std::size_t iNeuron, std::size_t iInput) const {
    su2double doutput_dinput = 0;
    for (auto jNeuron = 0u; jNeuron < total_layers[iLayer - 1]->getNNeurons(); jNeuron++) {
      doutput_dinput += weights_mat[iLayer - 1][iNeuron][jNeuron] * total_layers[iLayer - 1]->getdYdX(jNeuron, iInput);
    }
    return doutput_dinput;
  }
};

}  // namespace MLPToolbox