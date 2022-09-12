/*!
 * \file NeuralNetwork.cpp
 * \brief Implementation of the NeuralNetwork class to be used
 *      for evaluation of multi-layer perceptrons.
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
#include "../../../include/numerics/multilayer_perceptron/NeuralNetwork.hpp"
#include "../../../include/numerics/multilayer_perceptron/Layer.hpp"
#include <iostream>
#include "../../../include/numerics/multilayer_perceptron/ReadNeuralNetwork.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
using namespace std;


void NeuralNetwork::predict(vector<su2double> &inputs, vector<su2double*>&outputs, su2double**doutputs_dinputs){
    su2double x, x_norm, y, y_norm, tanx;
    size_t iNeuron, jNeuron, iLayer, nNeurons_current, nNeurons_previous;
    bool same_point{true};
    for(iNeuron=0; iNeuron<inputLayer->getNNeurons(); iNeuron++){
        x_norm = (inputs[iNeuron] - input_norm[iNeuron].first)/(input_norm[iNeuron].second - input_norm[iNeuron].first);
        if(abs(x_norm - inputLayer->getValue(iNeuron)) > 0) same_point = false;
        inputLayer->setValue(iNeuron, x_norm);
    }
    if(!same_point){
        for(iLayer=1; iLayer<n_hidden_layers + 2; iLayer++){
            nNeurons_current = total_layers[iLayer]->getNNeurons();
            nNeurons_previous = total_layers[iLayer - 1]->getNNeurons();

            
            switch (activation_function_types[iLayer])
                {
                case ENUM_ACTIVATION_FUNCTION::SMOOTH_SLOPE:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = ComputeX(iLayer, iNeuron);
                        y = x > 0 ? x : tanh(x);
                        total_layers[iLayer]->setValue(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::ELU:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = ComputeX(iLayer, iNeuron);
                        y = x > 0 ? x : (exp(x) - 1);
                        total_layers[iLayer]->setValue(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::LINEAR:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = ComputeX(iLayer, iNeuron);
                        y = x;
                        total_layers[iLayer]->setValue(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::RELU:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = ComputeX(iLayer, iNeuron);
                        y = x > 0.0 ? x : 0.0;
                        total_layers[iLayer]->setValue(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::NONE:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        y = 0.0;
                        total_layers[iLayer]->setValue(iNeuron, y);
                    }
                    break;
                default:
                    break;
                }
                
        }
    }
    for(iNeuron=0; iNeuron<outputLayer->getNNeurons(); iNeuron++){
        y_norm = outputLayer->getValue(iNeuron);
        y = y_norm*(output_norm[iNeuron].second - output_norm[iNeuron].first) + output_norm[iNeuron].first;
        
        // Setting output value
        *outputs[iNeuron] = y;
    }
}
// void NeuralNetwork::predict(vector<su2double> &inputs, vector<su2double*>&outputs, su2double**doutputs_dinputs){
//     // Evaluate the MLP based on the inputs. In case a pointer to the output gradient is provided, the partial
//     // derivatives of the outputs with respect to the inputs are evaluated as well.

//     su2double *y_values_old, // Perceptron values of the previous layer
//               *y_values_new, // Perceptron values of the current layer
//               **dy_dx_old,   // Derivatives of the previous layer
//               **dy_dx_new;   // Derivatives of the current layer

//     su2double x,    // Perceptron input
//             x_norm, // Normalized MLP input
//             y,      // Perceptron output
//             y_norm; // Normalized MLP output

//     size_t n_inputs{inputLayer->getNNeurons()}; // Number of inputs

//     bool compute_gradient = (doutputs_dinputs != nullptr); // Gradient computation option

//     // Stepping though the MLP layers, communicating the output from the previous layer to the next.
//     // TODO: make this recursive to increase performance and code simplification
//     for(size_t iLayer=0; iLayer < total_layers.size(); iLayer++){

//         // In case the input layer is called, input normalization is applied.
//         if(total_layers.at(iLayer)->isInput()){
//             // Generating an array for the input layer preceptron values.
//             y_values_old = new su2double[total_layers.at(iLayer)->getNNeurons()];

//             if(compute_gradient)
//                 dy_dx_old = new su2double*[total_layers.at(iLayer)->getNNeurons()];

//             for(size_t iNeuron=0; iNeuron<total_layers.at(iLayer)->getNNeurons(); iNeuron++){
//                 // Normalizing the inputs
//                 x_norm = (inputs[iNeuron] - input_norm.at(iNeuron).first)/(input_norm.at(iNeuron).second - input_norm.at(iNeuron).first);

//                 // Set perceptron values of input layer.
//                 total_layers.at(iLayer)->setValue(iNeuron, x_norm);
//                 y_values_old[iNeuron] = x_norm;

//                 if(compute_gradient){
//                     dy_dx_old[iNeuron] = new su2double[n_inputs];
//                     for(size_t jInput=0; jInput<n_inputs; jInput++){
//                         if(iNeuron == jInput) dy_dx_old[iNeuron][jInput] = 1.0/(input_norm.at(iNeuron).second - input_norm.at(iNeuron).first);
//                         else dy_dx_old[iNeuron][jInput] = 0.0;
//                     }
//                 }
//             }
//         }else{
//             // For the layers after the input layer, activation functions of each layer are called.
//             y_values_new = new su2double[total_layers.at(iLayer)->getNNeurons()];

//             if(compute_gradient)
//                 dy_dx_new = new su2double*[total_layers.at(iLayer)->getNNeurons()];

//             for(size_t iNeuron=0; iNeuron<total_layers.at(iLayer)->getNNeurons(); iNeuron++){
//                 // Multiplying the weights and output values of the previous layer to generate the current perceptron input
//                 x = GeometryToolbox::DotProduct(total_layers.at(iLayer-1)->getNNeurons(), y_values_old, weights_mat[iLayer-1][iNeuron]);

//                 // Generating the perceptron output
//                 y = total_layers.at(iLayer)->activation_function(iNeuron, x);

//                 // Storing perceptron output
//                 total_layers.at(iLayer)->setValue(iNeuron, y);
//                 y_values_new[iNeuron] = y;

//                 if(compute_gradient){
//                     dy_dx_new[iNeuron] = new su2double[n_inputs];
//                     for(size_t jInput=0; jInput < n_inputs; jInput++) dy_dx_new[iNeuron][jInput] = 0.0;

//                     for(size_t jNeuron=0; jNeuron<total_layers.at(iLayer-1)->getNNeurons(); jNeuron++){
//                         for(size_t jInput=0; jInput<inputLayer->getNNeurons(); jInput++){
//                             dy_dx_new[iNeuron][jInput] += total_layers.at(iLayer)->getdYdX(iNeuron) * weights_mat[iLayer-1][iNeuron][jNeuron] * dy_dx_old[jNeuron][jInput]; 
//                         }
//                     }
//                 }
                
//             }
//             // Resetting pointer to previous layer outputs
//             delete y_values_old;
//             y_values_old = y_values_new;

//             if(compute_gradient){
//                 for(size_t jNeuron=0; jNeuron<total_layers.at(iLayer-1)->getNNeurons(); jNeuron++){
//                     delete [] dy_dx_old[jNeuron];
//                 }
//                 delete [] dy_dx_old;
//                 dy_dx_old = dy_dx_new;
//             }
//         }
//     }
//     delete y_values_old;

//     // Dimensionalize the MLP outputs and gradients
//     for(size_t iOutput=0; iOutput<outputLayer->getNNeurons(); iOutput++){

//         // Dimensionalize the outputs
//         y_norm = outputLayer->getValue(iOutput);
//         y = y_norm*(output_norm.at(iOutput).second - output_norm.at(iOutput).first) + output_norm.at(iOutput).first;
        
//         // Setting output value
//         *outputs.at(iOutput) = y;

//         if(compute_gradient){
//             for(size_t iInput=0; iInput<n_inputs; iInput++){
//                 doutputs_dinputs[iOutput][iInput] = (output_norm.at(iOutput).second - output_norm.at(iOutput).first) * dy_dx_old[iOutput][iInput];
//             }
//             delete [] dy_dx_old[iOutput];
//         }
//     }

//     if(compute_gradient)
//         delete [] dy_dx_old;
// }

NeuralNetwork::NeuralNetwork(){
    inputLayer = nullptr;
    outputLayer = nullptr;
    n_hidden_layers = 0;
}


void NeuralNetwork::defineInputLayer(unsigned long n_neurons){
    //cout << "Creating an input layer with " << n_neurons << " neurons. "<< endl;
    inputLayer = new Layer(n_neurons);
    inputLayer->setInput(true);
    input_norm.resize(n_neurons);
}

void NeuralNetwork::defineOutputLayer(unsigned long n_neurons){
    //cout << "Creating an output layer with " << n_neurons << " neurons. "<< endl;
    outputLayer = new Layer(n_neurons);
    output_norm.resize(n_neurons);
}

void NeuralNetwork::push_hidden_layer(unsigned long n_neurons){
    Layer *newLayer = new Layer(n_neurons);
    //cout << "Creating a hidden layer with " << n_neurons << " neurons. "<< endl;
    hiddenLayers.push_back(newLayer);
    n_hidden_layers ++;
}

void NeuralNetwork::sizeWeights(){
    unsigned long i_layer{0};
    weights.resize(n_hidden_layers + 1);
    weights.at(i_layer).resize(inputLayer->getNNeurons());
    Layer * previouslayer = inputLayer;

    if(n_hidden_layers != 0){
        for(size_t i_hidden_layer=0; i_hidden_layer < n_hidden_layers; i_hidden_layer++){
            for(size_t i_neuron=0; i_neuron < previouslayer->getNNeurons(); i_neuron++){
                weights.at(i_layer).at(i_neuron).resize(hiddenLayers.at(i_hidden_layer)->getNNeurons());
            }
            previouslayer = hiddenLayers.at(i_hidden_layer);
            i_layer ++;
            weights.at(i_layer).resize(previouslayer->getNNeurons());
        }
    }
    for(size_t i_neuron=0; i_neuron < previouslayer->getNNeurons(); i_neuron++){
        weights.at(i_layer).at(i_neuron).resize(outputLayer->getNNeurons());
    }

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
    //weights.at(i_layer).at(i_neuron).at(j_neuron) = value;
    weights_mat[i_layer][j_neuron][i_neuron] = value;
}
void NeuralNetwork::displayNetwork(){
    cout << "Input layer: " << inputLayer->getNNeurons() << " neurons, Activation Function: " <<  inputLayer->getActivationType() << endl;
    for(size_t i=0; i<total_layers[1]->getNNeurons(); i++){
        for(size_t j=0; j<inputLayer->getNNeurons(); j++){
            cout << weights_mat[0][i][j] << " ";
        }
        su2double thingy[] = {1, 1, 1};
        cout << GeometryToolbox::DotProduct(inputLayer->getNNeurons(), weights_mat[0][i], thingy) << endl;
        cout << endl;
    }
    for(size_t i_layer=0; i_layer < hiddenLayers.size(); i_layer++){
        cout << "Hidden layer " << i_layer + 1 << ": " << hiddenLayers.at(i_layer)->getNNeurons() << " neurons, Activation Function: " <<  hiddenLayers.at(i_layer) ->getActivationType() << endl;
        for(size_t i=0; i<weights.at(i_layer+1).size(); i++){
            for(size_t j=0; j<weights.at(i_layer+1).at(i).size(); j++){
                cout << weights_mat[i_layer][i][j] << " ";
            }
            cout << endl;
    }
    }
    cout << "Output layer: "<<outputLayer->getNNeurons() <<" neurons, Activation Function: " << outputLayer->getActivationType() << endl;
}

void NeuralNetwork::setActivationFunction(unsigned long i_layer, string input)
{
    if(input.compare("linear") == 0){
        //activation_functions[i_layer] = new Linear();
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::LINEAR;
        return;
    }
    if(input.compare("relu") == 0){
        //activation_functions[i_layer] = new Relu();
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::RELU;
        return;
    }
    if((input.compare("smooth_slope") == 0)){
        //activation_functions[i_layer] = new SmoothSlope();
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SMOOTH_SLOPE;
        return;
    }
    if((input.compare("elu") == 0)){
        //activation_functions[i_layer] = new SmoothSlope();
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::ELU;
        return;
    }
    if(input.compare("none") == 0){
        //activation_functions[i_layer] = new None();
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::NONE;
        return;
    }
    //total_layers.at(i_layer)->setActivationFunction(input);
    return;
}
