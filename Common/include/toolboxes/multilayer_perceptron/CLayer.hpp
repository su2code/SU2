/*!
 * \file CLayer.hpp
 * \brief Declaration of artificial neural network interpolation class
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
#include "CNeuron.hpp"
#include "../../linear_algebra/blas_structure.hpp"

namespace MLPToolbox{
class CLayer
{
private:
    unsigned long number_of_neurons;    /*!< Neuron count in current layer */
    CNeuron * neurons;                   /*!< Array of neurons in current layer */
    bool is_input;                      /*!< Input layer identifyer */
    std::string activation_type;             /*!< Activation function type applied to the current layer*/
public:
    CLayer();
    CLayer(unsigned long n_neurons);
    ~CLayer(){delete neurons;}
    /*!
    * \brief Set current layer neuron count
    * \param[in] n_neurons - Number of neurons in this layer
    */
    void setNNeurons(unsigned long n_neurons);

    /*!
    * \brief Get the current layer neuron count
    * \return Neuron count
    */
    unsigned long getNNeurons() const {return number_of_neurons;};

    /*!
    * \brief Define current layer as input layer
    * \param[in] input - input layer identifyer
    */
    void setInput(bool def){is_input = def;};

    /*!
    * \brief Get input layer identifyer
    * \return input layer identifyer
    */
    bool isInput() const {return is_input;};

    /*!
    * \brief Set the output value of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \param[in] output_value - Activation function output
    */
    void setOutput(std::size_t i_neuron, su2double value){neurons[i_neuron].setOutput(value);}

    /*!
    * \brief Get the output value of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \return Neuron output value
    */
    su2double getOutput(std::size_t i_neuron) const {return neurons[i_neuron].getOutput();}

    /*!
    * \brief Set the input value of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \param[in] input_value - Activation function input
    */
    void setInput(std::size_t i_neuron, su2double value){neurons[i_neuron].setInput(value);}

    /*!
    * \brief Get the input value of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \return Neuron input value
    */
    su2double getInput(std::size_t i_neuron) const {return neurons[i_neuron].getInput();}

    /*!
    * \brief Set the bias value of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \param[in] bias_value - Bias value
    */
    void setBias(std::size_t i_neuron, su2double value){neurons[i_neuron].setBias(value);}

    /*!
    * \brief Get the bias value of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \return Neuron bias value
    */
    su2double getBias(std::size_t i_neuron){return neurons[i_neuron].getBias();}

    /*!
    * \brief Get the output-input gradient of a neuron in the layer
    * \param[in] i_neuron - Neuron index
    * \return Gradient of neuron output wrt input
    */
    su2double getdYdX(std::size_t i_neuron){return neurons[i_neuron].getGradient();}

    /*!
    * \brief Get the activation function name applied to this layer
    * \return name of the activation function
    */
    string getActivationType(){return activation_type;}
    
};

}
