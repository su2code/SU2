/*!
 * \file Layer.cpp
 * \brief Implementation of the Layer class to be used in the NeuralNetwork
 *      class
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
#include "../../../include/numerics/multilayer_perceptron/Layer.hpp"
#include <cstring>
using namespace std;

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

// void Layer::setActivationFunction(string input){
//     activation_type = input;
//     for(size_t i=0; i<number_of_neurons; i++){
//         neurons[i].setFunctionType(input);
//     }
// }
// void Layer::sayhi(){
//     for(size_t i=0; i<number_of_neurons; i++){
//         neurons[i].sayhi();
//     }
// }
