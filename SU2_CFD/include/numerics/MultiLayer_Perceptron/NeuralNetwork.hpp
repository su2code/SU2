/*!
 * \file NeuralNetwork.hpp
 * \brief Delarations of numerics classes for definitition of a multilayer perceptron.
 * \version 7.3.1 "Blackbird"
 * \author E.C.Bunschoten
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

/*!
 * \class Neuron
 * \brief 
 * \details This class contains all the attributes for the Neuron class, to be used in a Layer of a NeuralNetwork class.
 * \author E.C.Bunschoten
 */

class Neuron
{
private:
    string activation_type{"linear"};   /*!< \brief Activation function type. */
    unsigned long i_neuron;             /*!< \brief Identification number within layer. */
    su2double value;                    /*!< \brief Output of the activation function. */
    su2double dy_dx;                    /*!< \brief Gradient of the activation function wrt the input. */
    su2double bias{0};                  /*!< \brief Bias of the neuron. */

public:
    /*!
     * \brief Constructor of the class.
    */
    Neuron();

    /*!
     * \brief Set the identification number of the neuron within the layer.
     * \param[in] input - Identification number.
    */
    void setNumber(unsigned long input){i_neuron = input;}

    /*!
     * \brief Get the identification number of the neuron within the layer.
     * \param[out] i_neuron - Identification number.
    */
    unsigned long getNumber(){return i_neuron;}

    /*!
     * \brief Set the activation function of the neuron.
     * \param[in] input - Activation function name.
    */
    void setFunctionType(string input){activation_type = input;}

    /*!
     * \brief Get the neuron activation function
     * \param[out] activation_type - Activation function name.
    */
    string getFunctionType(){return activation_type;}

    /*!
     * \brief Call the activation function.
     * \param[in] x - Neuron input.
    */
    su2double activation_function(su2double x);

    /*!
     * \brief Set the output value of the neuron.
     * \param[in] input - Neuron input.
    */
    void setValue(su2double input){value = input;}

    /*!
     * \brief Get the activation function output.
     * \param[out] value - Neuron output.
    */
    su2double getValue(){return value;}

    /*!
     * \brief Set the neuron bias.
     * \param[in] input - bias value.
    */
    void setBias(su2double input){bias = input;}

    /*!
     * \brief Get the neuron bias.
     * \param[out] bias - bias value.
    */
    su2double getBias(){return bias;}

    /*!
     * \brief Get the gradient of the activation function.
     * \param[out] dy_dx - derivative of the activation function.
    */
    su2double getGradient(){return dy_dx;}
    /*!
     * \brief Destructor of the class.
    */
    ~Neuron(){};
};

/*!
 * \class Layer
 * \brief .
 * \details This class contains all the attibutes regarding the Layer class, to be used within the NeuralNetwork class.
 * \author E.C.Bunschoten
 */

class Layer
{
private:
    unsigned long number_of_neurons;    /*!< \brief Number of neurons within the layer. */
    Neuron * neurons;                   /*!< \brief Array of neurons making up the layer.*/
    bool input;                         /*!< \brief Input layer flag.*/
    string activation_type;             /*!< \brief Activation function assigned to the layer. */

public:
    /*!
     * \brief Constructor of the class
    */
    Layer();

    /*!
     * \brief Constructor of the class.
     * \param[in] n_neurons - Number of neurons in this layer.
    */
    Layer(unsigned long n_neurons);

    /*!
     * \brief Set the number of neurons for this layer.
     * \param[in] n_neurons - Neuron count.
    */
    void setNNeurons(unsigned long n_neurons);

    /*!
     * \brief Get the neuron count for this layer.
     * \param[out] number_of_neurons - Neuron count.
    */
    unsigned long getNNeurons(){return number_of_neurons;};

    /*!
     * \brief Define the layer as an input layer.
     * \param[in] def - Input flag.
    */
    void setInput(bool def){input = def;};

    /*!
     * \brief Get input layer flag.
     * \param[out] flag - input layer flag.
    */
    bool isInput(){return input;};

    /*!
     * \brief Set the activation function of the layer.
     * \param[in] input - Name of the activation function. 
    */
    void setActivationFunction(string input);

    /*!
     * \brief Get the output from a neuron within the layer.
     * \param[in] i_neuron - Neuron index.
     * \param[in] x - Neuron input.
     * \param[out] y - Neuron output.
    */
    su2double activation_function(size_t i_neuron, su2double x){return neurons[i_neuron].activation_function(x);};
    
    /*!
     * \brief Set the neuron output value for a specific neuron.
     * \param[in] i_neuron - Neuron index.
     * \param[in] value - Input value.
    */
    void setValue(size_t i_neuron, su2double value){neurons[i_neuron].setValue(value);}
    
    /*!
     * \brief Get the output from a specific neuron.
     * \param[in] i_neuron - Neuron index.
     * \param[out] y - Activation function output.
    */
    su2double getValue(size_t i_neuron){return neurons[i_neuron].getValue();}
    
    /*!
     * \brief Set the bias for a neuron within the layer.
     * \param[in] i_neuron - Neuron index.
     * \param[in] bias - Bias value.
    */
    void setBias(size_t i_neuron, su2double value){neurons[i_neuron].setBias(value);}
    
    /*!
     * \brief Get the bias of a neuron within the layer.
     * \param[in] i_neuron - Neuron index.
     * \param[out] bias - Bias value.
    */
    su2double getBias(size_t i_neuron){return neurons[i_neuron].getBias();}
    
    /*!
     * \brief Get the activation function derivative of a neuron within the layer.
     * \param[in] i_neuron - Neuron index.
     * \param[out] dy_dx - Gradient of the neuron output wrt its input.
    */
    su2double getdYdX(size_t i_neuron){return neurons[i_neuron].getGradient();}
    
    /*!
     * \brief Get the activation function type assigned to the neurons in this layer.
     * \param[out] activation_type - Name of the activation function.
    */
    string getActivationType(){return activation_type;}
   
    /*!
     * \brief Destructor of the class
    */
    ~Layer(){
        delete [] neurons;
    };
};

/*!
 * \class NeuralNetwork
 * \brief Template class for a multilayer perceptron.
 * \details This class 
 * \author E.C.Bunschoten
 */

class NeuralNetwork
{
private:
    vector<string> input_names;     /*<! \brief Vector of input variable names. */
    vector<string> output_names;    /*<! \brief Vector of output variable names. */

    unsigned long n_hidden_layers;  /*<! \brief Number of hidden layers within the neural network. */
    
    Layer *inputLayer;              /*<! \brief Pointer to the input layer of the neural network.*/
    Layer *outputLayer;             /*<! \brief Poiner to the output layer of the neural network. */

    vector< Layer *> hiddenLayers; /*<! \brief Vector of hidden layer pointers. */
    vector<Layer *>total_layers;   /*<! \brief Vector containing pointers to all the sequential layers within the network.*/

    vector<su2activematrix> weights_mat;    /*<! \brief Vector of weight matrices defining the network. */
    vector<vector<su2double>> values;       /*<! \brief Vector containing the output values of the neurons within the network.*/
    
    vector<pair<su2double, su2double>> input_norm;  /*<! \brief Vector of input normalization values (minimum, maximum) */
    vector<pair<su2double, su2double>> output_norm; /*<! \brief Vector of output normalization values (minimum, maximum) */

public:

    /*!
     *\brief Constructor of the class.
    */
    NeuralNetwork();

    /*!
     *\brief Destructor of the class.
    */
    ~NeuralNetwork(){
    if(inputLayer != nullptr){
        delete inputLayer;
    }
    if(outputLayer != nullptr){
        delete outputLayer;
    }
    for(size_t i=0; i<hiddenLayers.size(); i++){
        delete hiddenLayers.at(i);
    }
    };

    /*!
     *\brief Generate the network input layer.
     *\param[in] n_neurons - Number of MLP inputs.
    */
    void defineInputLayer(unsigned long n_neurons);

    /*!
     *\brief Generate the MLP output layer.
     *\param[in] n_neurons - Number of MLP outputs.
    */
    void defineOutputLayer(unsigned long n_neurons);

    /*!
     *\brief Generate a new hidden layer in the MLP.
     *\param[in] n_neurons - Size of the hidden layer.
    */
    void push_hidden_layer(unsigned long n_neurons);

    /*!
     *\brief Set the weight value of a connection within the network.
     *\param[in] i_layer - Layer index.
     *\param[in] i_neuron - Neuron index of sending neuron.
     *\param[in] j_neuron - Neuron index of receiving neuron.
     *\param[in] value - Weight value.
    */
    void setWeight(unsigned long i_layer, unsigned long i_neuron, unsigned long j_neuron, su2double value);
    
    /*!
     *\brief Set the bias value of a neuron within a layer.
     *\param[in] i_layer - Layer index.
     *\param [in] i_neuron - Neuron index.
     *\param[in] bias - Bias value.
    */
    void setBias(unsigned long i_layer, unsigned long i_neuron, su2double value){total_layers.at(i_layer)->setBias(i_neuron, value);}
    
    /*!
     *\brief Set the activation function type of a neuron within the network.
     *\param[in] i_layer - Layer index.
     *\param[in] name - Activation function name.
    */
    void setActivationFunction(unsigned long i_layer, string input){total_layers.at(i_layer)->setActivationFunction(input);}
    
    /*!
     *\brief Display the layout and weights of the neural network on the command line.
    */
    void displayNetwork();
    
    /*!
     *\brief Size the connectivity of the neural network.
    */
    void sizeWeights();
    
    /*!
     *\brief Get the number of connectivity regions within the neural network.
     *\param[out] n_regions - Number of sections between layers.
    */
    unsigned long getNWeightLayers(){return weights_mat.size();}
    
    /*!
     *\brief Get the total layer count in the neural network.
     *\param[out] n_layers - Layer count.
    */
    unsigned long getNLayers(){return total_layers.size();}

    /*!
     *\brief Get the number of neurons of a specific layer.
     *\param[in] iLayer - Layer index.
     *\param[out] n_neurons - Number of neurons in the requested layer.
    */
    unsigned long getNNeurons(unsigned long iLayer){return total_layers.at(iLayer)->getNNeurons();}

    /*!
     *\brief Evaluate the output of the neural network.
     *\param[in] inputs - Vector of values fed to the input layer.
     *\param[in] outputs - Vector containing pointers to the output values.
     *\param[in] doutputs_dinputs - Optional output for the gradients of the outputs with respect to the inputs.
    */
    void predict(vector<su2double> &inputs, vector<su2double*> &outputs, su2double**doutputs_dinputs=nullptr);

    /*!
     *\brief Set the input normalization values
     *\param[in] iInput - Input index.
     *\param[in] input_min - Minimum input value.
     *\param[in] input_max - Maximum intput value.
    */
    void SetInputNorm(unsigned long iInput, su2double input_min, su2double input_max){input_norm.at(iInput) = make_pair(input_min, input_max);}
   
    /*!
     *\brief Get the input normalization values
     *\param[in] iInput - Input index.
     *\param[out] input_norm - Pair with input normalization values (minimum first, maximum second)
    */
    pair<su2double, su2double> GetInputNorm(unsigned long iInput){return input_norm.at(iInput);}
    
    /*!
     *\brief Set the output normalization values
     *\param[in] iOutput - Output index.
     *\param[in] output_min - Minimum output value.
     *\param[in] output_max - Maximum output value.
    */
    void SetOutputNorm(unsigned long iOutput, su2double output_min, su2double output_max){output_norm.at(iOutput) = make_pair(output_min, output_max);}
    
    /*!
     *\brief Get the output normalization values
     *\param[in] iOutput - Output index.
     *\param[out] output_norm - Pair with output normalization values (minimum first, maximum second)
    */
    pair<su2double, su2double> GetOutputNorm(unsigned long iOutput){return output_norm.at(iOutput);}

    /*!
     *\brief Append an output variable name.
     *\param[in] name - Output variable name
    */
    void PushOutputName(string name){output_names.push_back(name);}
    
    /*!
     *\brief Append an input variable name.
     *\param[in] name - Input variable name
    */
    void PushInputName(string name){input_names.push_back(name);}
   
    /*!
     *\brief Get MLP input name.
     *\param[in] iInput - Input variable index.
     *\param[out] name - Input variable name.
    */
    string GetInputName(size_t iInput){return input_names[iInput];}
    
    /*!
     *\brief Get MLP output name.
     *\param[in] iOutput - Output variable index.
     *\param[out] name - Output variable name.
    */
    string GetOutputName(size_t iOutput){return output_names[iOutput];}

    /*!
     *\brief Get the number of inputs of the MLP
     *\param[out] nInputs - number of inputs of the MLP
    */
    size_t GetnInputs(){return input_names.size();}
    
    /*!
     *\brief Get the number of outputs of the MLP
     *\param[out] nOutputs - number of outputs of the MLP
    */
    size_t GetnOutputs(){return output_names.size();}

};





