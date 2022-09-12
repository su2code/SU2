/*!
 * \file Neuron.cpp
 * \brief Implementation of the neuron class to be used within the 
 *      Layer class as a part of the NeuralNetwork class.
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
#include "../../../include/numerics/multilayer_perceptron/Neuron.hpp"
#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;
// Neuron::Neuron()
// {
// }

// su2double Neuron::activation_function(su2double x)
// {
//     x += bias;
//     su2double tanx;
//     switch(ActivationFunction){
//         case ENUM_ACTIVATION_FUNCTION::LINEAR:
//             dy_dx = 1;
//             return x;
//             break;
//         case ENUM_ACTIVATION_FUNCTION::RELU:
//             if(x > 0.0){
//                 dy_dx = 1.0;
//                 return x;
//             }else{
//                 dy_dx = 0.0;
//                 return 0.0;
//             }
//             break;
//         case ENUM_ACTIVATION_FUNCTION::SMOOTH_SLOPE:
//             tanx = tanh(x);
//             if(tanx > x){
//                 dy_dx = 4.0 / pow((exp(-x) + exp(x)), 2); 
//                 return tanx;
//             }else{
//                 dy_dx = 1.0;
//                 return x;
//             }
//             break;
//         default:
//             SU2_MPI::Error(string("Unknown activation function type: ")+activation_type,
//                     CURRENT_FUNCTION);
//             return 0.0;
//         break;
//     }
//     // if(activation_type.compare("smooth_slope") == 0){
//     //     su2double tanx = tanh(x);
//     //     if(tanx > x){
//     //         dy_dx = 4.0 / pow((exp(-x) + exp(x)), 2); 
//     //         return tanx;
//     //     }else{
//     //         dy_dx = 1.0;
//     //         return x;
//     //     }

//     // }
//     // if(activation_type.compare("linear") == 0){
//     //     dy_dx = 1.0;
//     //     return x;
//     // }
//     // if(activation_type.compare("relu") == 0){
//     //     su2double zero{0.0};
//     //     if(x > zero){
//     //         dy_dx = 1.0;
//     //         return x;
//     //     }else{
//     //         dy_dx = 0.0;
//     //         return zero;
//     //     }
//     // }
    
//     // SU2_MPI::Error(string("Unknown activation function type: ")+activation_type,
//     //                CURRENT_FUNCTION);   
    
//     // return x;
// }

// void Neuron::setFunctionType(string input){
//     activation_type = input;

//     if(input.compare("smooth_slope") == 0) ActivationFunction = ENUM_ACTIVATION_FUNCTION::SMOOTH_SLOPE; return;
    
//     if(input.compare("linear") == 0) ActivationFunction = ENUM_ACTIVATION_FUNCTION::LINEAR; return;

//     if(input.compare("relu") == 0) ActivationFunction = ENUM_ACTIVATION_FUNCTION::RELU; return;

//     if(input.compare("none") == 0) ActivationFunction = ENUM_ACTIVATION_FUNCTION::NONE; return;
// }