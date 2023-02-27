/*!
 * \file CLayer.hpp
 * \brief Declaration of artificial neural network interpolation class
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
#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "../../CConfig.hpp"
#include "../../linear_algebra/blas_structure.hpp"
#include "CNeuron.hpp"

namespace MLPToolbox {
class CLayer {
  /*!
   *\class CLayer
   *\brief This class functions as one of the hidden, input, or output layers in the multi-layer perceptron
   * class. The CLayer class is used to communicate information (activation function inputs and outputs and gradients)
   * between the CNeuralNetwork class and the CNeuron class. Currently, only a single activation function
   * can be applied to the neuron inputs within the layer. 
   */
 private:
  unsigned long number_of_neurons; /*!< Neuron count in current layer */
  su2vector<CNeuron>neurons;                /*!< Array of neurons in current layer */
  bool is_input;                   /*!< Input layer identifyer */
  std::string activation_type;     /*!< Activation function type applied to the current layer*/
 public:
  CLayer();
  CLayer(unsigned long n_neurons);
  /*!
   * \brief Set current layer neuron count
   * \param[in] n_neurons - Number of neurons in this layer
   */
  void SetNNeurons(unsigned long n_neurons);

  /*!
   * \brief Get the current layer neuron count
   * \return Neuron count
   */
  unsigned long GetNNeurons() const { return number_of_neurons; }

  /*!
   * \brief Define current layer as input layer
   * \param[in] input - input layer identifyer
   */
  void SetInput(bool def) { is_input = def; }

  /*!
   * \brief Get input layer identifyer
   * \return input layer identifyer
   */
  bool IsInput() const { return is_input; }

  /*!
   * \brief Set the output value of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \param[in] output_value - Activation function output
   */
  void SetOutput(std::size_t i_neuron, su2double value) { neurons[i_neuron].SetOutput(value); }

  /*!
   * \brief Get the output value of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \return Neuron output value
   */
  su2double GetOutput(std::size_t i_neuron) const { return neurons[i_neuron].GetOutput(); }

  /*!
   * \brief Set the input value of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \param[in] input_value - Activation function input
   */
  void SetInput(std::size_t i_neuron, su2double value) { neurons[i_neuron].SetInput(value); }

  /*!
   * \brief Get the input value of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \return Neuron input value
   */
  su2double GetInput(std::size_t i_neuron) const { return neurons[i_neuron].GetInput(); }

  /*!
   * \brief Set the bias value of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \param[in] bias_value - Bias value
   */
  void SetBias(std::size_t i_neuron, su2double value) { neurons[i_neuron].SetBias(value); }

  /*!
   * \brief Get the bias value of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \return Neuron bias value
   */
  su2double GetBias(std::size_t i_neuron) const { return neurons[i_neuron].GetBias(); }

  /*!
   * \brief Get the output-input gradient of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \return Gradient of neuron output wrt input
   */
  su2double GetdYdX(std::size_t i_neuron, std::size_t iInput) const { return neurons[i_neuron].GetGradient(iInput); }

  /*!
   * \brief Get the output-input gradient of a neuron in the layer
   * \param[in] i_neuron - Neuron index
   * \return Gradient of neuron output wrt input
   */
  void SetdYdX(std::size_t i_neuron, std::size_t iInput, su2double dy_dx) {
    neurons[i_neuron].SetGradient(iInput, dy_dx);
  }

  /*!
   * \brief Size neuron output derivative wrt network inputs.
   * \param[in] nInputs - Number of network inputs.
   */
  void SizeGradients(std::size_t nInputs) {
    for (auto iNeuron = 0u; iNeuron < number_of_neurons; iNeuron++) neurons[iNeuron].SizeGradient(nInputs);
  }

  /*!
   * \brief Get the activation function name applied to this layer
   * \return name of the activation function
   */
  string GetActivationType() const { return activation_type; }
};

}  // namespace MLPToolbox
