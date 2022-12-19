/*!
 * \file CNeuron.hpp
 * \brief Declaration of artificial neural network perceptron class
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

#include "../../CConfig.hpp"

namespace MLPToolbox{
class CNeuron
{
private:
    unsigned long i_neuron;         /*!< Neuron identification number */
    su2double output{0},            /*!< Output value of the current neuron */
              input{0},             /*!< Input value of the current neuron */
              doutput_dinput{0},    /*!< Gradient of output with respect to input */
              bias{0};              /*!< Bias value at current neuron */
    su2vector<su2double> doutput_dinputs;

public:

    /*!
    * \brief Set neuron identification number
    * \param[in] input - Identification number
    */
    void setNumber(unsigned long input){ i_neuron = input; }

    /*!
    * \brief Get neuron identification number
    * \return Identification number
    */
    unsigned long getNumber() const { return i_neuron; }

    /*!
    * \brief Set neuron output value
    * \param[in] input - activation function output value
    */
    void setOutput(su2double input){ output = input; }

    /*!
    * \brief Get neuron output value
    * \return Output value
    */
    su2double getOutput() const { return output; }

    /*!
    * \brief Set neuron input value
    * \param[in] input - activation function input value
    */
    void setInput(su2double x){ input = x; }

    /*!
    * \brief Get neuron input value
    * \return input value
    */
    su2double getInput() const { return input; }

    /*!
    * \brief Set neuron bias
    * \param[in] input - bias value
    */
    void setBias(su2double input){ bias = input; }

    /*!
    * \brief Get neuron bias value
    * \return bias value
    */
    su2double getBias() const { return bias; }

    void sizeGradient(std::size_t nInputs) { doutput_dinputs.resize(nInputs); }
    /*!
    * \brief Set neuron output gradient with respect to its input value
    * \param[in] input - Derivative of activation function with respect to input
    */
    void setGradient(std::size_t iInput, su2double input){ doutput_dinputs[iInput] = input; }

    /*!
    * \brief Get neuron output gradient with respect to input value
    * \return output gradient wrt input value
    */
    su2double getGradient(std::size_t iInput) const { return doutput_dinputs[iInput]; }
};

}

