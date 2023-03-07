/*!
 * \file CNeuron.hpp
 * \brief Declaration of artificial neural network perceptron class
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

namespace MLPToolbox {
class CNeuron {
  /*!
   *\class CNeuron
   *\brief This class functions as a neuron within the CLayer class, making up the
   * CNeuralNetwork class. The CNeuron class functions as a location to store activation
   * function inputs and outputs, as well as gradients and biases. These are accessed
   * through the CLayer class for network evalution operations.
   */
 private:
  unsigned long i_neuron; /*!< Neuron identification number */
  su2double output{0},    /*!< Output value of the current neuron */
      input{0},           /*!< Input value of the current neuron */
      doutput_dinput{0},  /*!< Gradient of output with respect to input */
      bias{0};            /*!< Bias value at current neuron */
  su2vector<su2double> doutput_dinputs;

 public:
  /*!
   * \brief Set neuron identification number
   * \param[in] input - Identification number
   */
  void SetNumber(unsigned long input) { i_neuron = input; }

  /*!
   * \brief Get neuron identification number
   * \return Identification number
   */
  unsigned long GetNumber() const { return i_neuron; }

  /*!
   * \brief Set neuron output value
   * \param[in] input - activation function output value
   */
  void SetOutput(su2double input) { output = input; }

  /*!
   * \brief Get neuron output value
   * \return Output value
   */
  su2double GetOutput() const { return output; }

  /*!
   * \brief Set neuron input value
   * \param[in] input - activation function input value
   */
  void SetInput(su2double x) { input = x; }

  /*!
   * \brief Get neuron input value
   * \return input value
   */
  su2double GetInput() const { return input; }

  /*!
   * \brief Set neuron bias
   * \param[in] input - bias value
   */
  void SetBias(su2double input) { bias = input; }

  /*!
   * \brief Get neuron bias value
   * \return bias value
   */
  su2double GetBias() const { return bias; }

  /*!
   * \brief Size the derivative of the neuron output wrt MLP inputs.
   * \param[in] nInputs - Number of MLP inputs.
   */
  void SizeGradient(std::size_t nInputs) { doutput_dinputs.resize(nInputs); }

  /*!
   * \brief Set neuron output gradient with respect to its input value
   * \param[in] input - Derivative of activation function with respect to input
   */
  void SetGradient(std::size_t iInput, su2double input) { doutput_dinputs[iInput] = input; }

  /*!
   * \brief Get neuron output gradient with respect to input value
   * \return output gradient wrt input value
   */
  su2double GetGradient(std::size_t iInput) const { return doutput_dinputs[iInput]; }
};

}  // namespace MLPToolbox