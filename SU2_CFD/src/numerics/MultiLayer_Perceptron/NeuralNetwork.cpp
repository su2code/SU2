/*!
 * \file NeuralNetwork.cpp
 * \brief Implementation of numerics classes for the evaluation of Multi-Layer Perceptrons.
 * \version 7.3.1 "Blackbird"
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
#include "../../../include/numerics/MultiLayer_Perceptron/NeuralNetwork.hpp"
#include <iostream>
#include "../../../include/numerics/MultiLayer_Perceptron/ReadNeuralNetwork.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


using namespace std;

NeuralNetwork::NeuralNetwork(){
    inputLayer = nullptr;
    outputLayer = nullptr;
    n_hidden_layers = 0;
}

void NeuralNetwork::predict(vector<su2double> &inputs, vector<su2double*>&outputs, su2double**doutputs_dinputs){
    // Evaluate the MLP based on the inputs. In case a pointer to the output gradient is provided, the partial
    // derivatives of the outputs with respect to the inputs are evaluated as well.

    su2double *y_values_old, // Perceptron values of the previous layer
              *y_values_new, // Perceptron values of the current layer
              **dy_dx_old,   // Derivatives of the previous layer
              **dy_dx_new;   // Derivatives of the current layer

    su2double x,    // Perceptron input
            x_norm, // Normalized MLP input
            y,      // Perceptron output
            y_norm; // Normalized MLP output

    size_t n_inputs{inputLayer->getNNeurons()}; // Number of inputs

    bool compute_gradient = (doutputs_dinputs != nullptr); // Gradient computation option

    // Stepping though the MLP layers, communicating the output from the previous layer to the next.
    // TODO: make this recursive to increase performance and code simplification
    for(size_t iLayer=0; iLayer < total_layers.size(); iLayer++){

        // In case the input layer is called, input normalization is applied.
        if(total_layers.at(iLayer)->isInput()){
            // Generating an array for the input layer preceptron values.
            y_values_old = new su2double[total_layers.at(iLayer)->getNNeurons()];

            if(compute_gradient)
                dy_dx_old = new su2double*[total_layers.at(iLayer)->getNNeurons()];

            for(size_t iNeuron=0; iNeuron<total_layers.at(iLayer)->getNNeurons(); iNeuron++){
                // Normalizing the inputs
                x_norm = (inputs[iNeuron] - input_norm.at(iNeuron).first)/(input_norm.at(iNeuron).second - input_norm.at(iNeuron).first);

                // Set perceptron values of input layer.
                total_layers.at(iLayer)->setValue(iNeuron, x_norm);
                y_values_old[iNeuron] = x_norm;

                if(compute_gradient){
                    dy_dx_old[iNeuron] = new su2double[n_inputs];
                    for(size_t jInput=0; jInput<n_inputs; jInput++){
                        if(iNeuron == jInput) dy_dx_old[iNeuron][jInput] = 1.0/(input_norm.at(iNeuron).second - input_norm.at(iNeuron).first);
                        else dy_dx_old[iNeuron][jInput] = 0.0;
                    }
                }
            }
        }else{
            // For the layers after the input layer, activation functions of each layer are called.
            y_values_new = new su2double[total_layers.at(iLayer)->getNNeurons()];

            if(compute_gradient)
                dy_dx_new = new su2double*[total_layers.at(iLayer)->getNNeurons()];

            for(size_t iNeuron=0; iNeuron<total_layers.at(iLayer)->getNNeurons(); iNeuron++){
                // Multiplying the weights and output values of the previous layer to generate the current perceptron input
                x = GeometryToolbox::DotProduct(total_layers.at(iLayer-1)->getNNeurons(), y_values_old, weights_mat[iLayer-1][iNeuron]);

                // Generating the perceptron output
                y = total_layers.at(iLayer)->activation_function(iNeuron, x);

                // Storing perceptron output
                total_layers.at(iLayer)->setValue(iNeuron, y);
                y_values_new[iNeuron] = y;

                if(compute_gradient){
                    dy_dx_new[iNeuron] = new su2double[n_inputs];
                    for(size_t jInput=0; jInput < n_inputs; jInput++) dy_dx_new[iNeuron][jInput] = 0.0;

                    for(size_t jNeuron=0; jNeuron<total_layers.at(iLayer-1)->getNNeurons(); jNeuron++){
                        for(size_t jInput=0; jInput<inputLayer->getNNeurons(); jInput++){
                            dy_dx_new[iNeuron][jInput] += total_layers.at(iLayer)->getdYdX(iNeuron) * weights_mat[iLayer-1][iNeuron][jNeuron] * dy_dx_old[jNeuron][jInput]; 
                        }
                    }
                }
                
            }
            // Resetting pointer to previous layer outputs
            delete y_values_old;
            y_values_old = y_values_new;

            if(compute_gradient){
                for(size_t jNeuron=0; jNeuron<total_layers.at(iLayer-1)->getNNeurons(); jNeuron++){
                    delete [] dy_dx_old[jNeuron];
                }
                delete [] dy_dx_old;
                dy_dx_old = dy_dx_new;
            }
        }
    }
    delete y_values_old;

    // Dimensionalize the MLP outputs and gradients
    for(size_t iOutput=0; iOutput<outputLayer->getNNeurons(); iOutput++){

        // Dimensionalize the outputs
        y_norm = outputLayer->getValue(iOutput);
        y = y_norm*(output_norm.at(iOutput).second - output_norm.at(iOutput).first) + output_norm.at(iOutput).first;
        
        // Setting output value
        *outputs.at(iOutput) = y;

        if(compute_gradient){
            for(size_t iInput=0; iInput<n_inputs; iInput++){
                doutputs_dinputs[iOutput][iInput] = (output_norm.at(iOutput).second - output_norm.at(iOutput).first) * dy_dx_old[iOutput][iInput];
            }
            delete [] dy_dx_old[iOutput];
        }
    }

    if(compute_gradient)
        delete [] dy_dx_old;
}



void NeuralNetwork::defineInputLayer(unsigned long n_neurons){
    inputLayer = new Layer(n_neurons);
    inputLayer->setInput(true);
    input_norm.resize(n_neurons);
}

void NeuralNetwork::defineOutputLayer(unsigned long n_neurons){
    outputLayer = new Layer(n_neurons);
    output_norm.resize(n_neurons);
}

void NeuralNetwork::push_hidden_layer(unsigned long n_neurons){
    Layer *newLayer = new Layer(n_neurons);
    hiddenLayers.push_back(newLayer);
    n_hidden_layers ++;
}

void NeuralNetwork::sizeWeights(){

    total_layers.resize(n_hidden_layers + 2);
    total_layers.at(0) = inputLayer;
    for(size_t iLayer=0; iLayer<n_hidden_layers; iLayer++){
        total_layers.at(iLayer + 1) = hiddenLayers.at(iLayer);
    }
    total_layers.at(total_layers.size()-1) = outputLayer;

    weights_mat.resize(n_hidden_layers+1);
    weights_mat[0].resize(hiddenLayers[0]->getNNeurons(), inputLayer->getNNeurons());


    for(size_t iLayer=1; iLayer<n_hidden_layers; iLayer++){
        weights_mat[iLayer].resize(hiddenLayers[iLayer]->getNNeurons(), hiddenLayers[iLayer-1]->getNNeurons());
    }
    weights_mat[n_hidden_layers].resize(outputLayer->getNNeurons(), hiddenLayers[n_hidden_layers-1]->getNNeurons());

}

void NeuralNetwork::setWeight(unsigned long i_layer, unsigned long i_neuron, unsigned long j_neuron, su2double value){
    weights_mat[i_layer][j_neuron][i_neuron] = value;
}
void NeuralNetwork::displayNetwork(){

    for(size_t iLayer=0; iLayer<total_layers.size()-1; iLayer++){
        if(total_layers.at(iLayer)->isInput()){
            cout << "Input Layer: ";
        }else{
            cout << "Hidden Layer: ";
        }
        cout << total_layers.at(iLayer)->getNNeurons() << " neurons, Activation Function: " << total_layers.at(iLayer)->getActivationType() << endl;
        for(size_t iNeuron=0; iNeuron<total_layers.at(iLayer)->getNNeurons(); iNeuron++){
            for(size_t jNeuron=0; jNeuron<total_layers.at(iLayer + 1)->getNNeurons(); jNeuron++){
                cout << weights_mat[iLayer][jNeuron][iNeuron] << " ";
            }
            cout << endl;
        }

    }
    cout << "Output Layer: " << outputLayer->getNNeurons() << " neurons, Activation Function: " << outputLayer->getActivationType() << endl;
}


Layer::Layer() : Layer(1) {};

Layer::Layer(unsigned long n_neurons) : number_of_neurons{n_neurons}, input{false}
{
    neurons = new Neuron[n_neurons];
    for(size_t i=0; i<number_of_neurons; i++){
        neurons[i].setNumber(i+1);
    }
}

void Layer::setNNeurons(unsigned long n_neurons){
    if(number_of_neurons != n_neurons){
        delete [] neurons;
        neurons = new Neuron[n_neurons];
        for(size_t i=0; i<number_of_neurons; i++){
            neurons[i].setNumber(i+1);
        }
    }
}

void Layer::setActivationFunction(string input){
    activation_type = input;
    for(size_t i=0; i<number_of_neurons; i++){
        neurons[i].setFunctionType(input);
    }
}

Neuron::Neuron()
{
}

su2double Neuron::activation_function(su2double x)
{
    x += bias;
    /* TODO: make this a switch-case */
    if(activation_type.compare("swish") == 0){
        su2double sigmoid = 1.0/(1 + exp(-x));
        su2double dsigmoid_dx = -(1/pow(1 + exp(-x), 2))*(-exp(-x));
        dy_dx = 1.0*sigmoid + x*dsigmoid_dx;
        return x * sigmoid;
    }
    if(activation_type.compare("sigmoid") == 0){
        dy_dx = -(1/pow(1 + exp(-x), 2))*(-exp(-x));
        return 1.0/(1 + exp(-x));
    }
    if(activation_type.compare("smoothSlope") == 0){
        su2double tanx = tanh(x);
        if(tanx > x){
            dy_dx = 4.0 / pow((exp(-x) + exp(x)), 2); 
            return tanx;
        }else{
            dy_dx = 1.0;
            return x;
        }
    }
    if(activation_type.compare("relu") == 0){
        su2double zero{0.0};
        if(x > zero){
            dy_dx = 1.0;
            return x;
        }else{
            dy_dx = 0.0;
            return zero;
        }
    }
    if(activation_type.compare("linear") == 0){
        dy_dx = 1.0;
        return x;
    }
    SU2_MPI::Error(string("Unknown activation function type: ")+activation_type,
                   CURRENT_FUNCTION);   
    
    return x;
}