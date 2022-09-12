/*!
 * \file Layer.hpp
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

#include "../../../../Common/include/CConfig.hpp"
#include "Neuron.hpp"
#include "../../../../Common/include/linear_algebra/blas_structure.hpp"
using namespace std;

class Layer
{
private:
    unsigned long number_of_neurons;
    Neuron * neurons;
    bool input;
    string activation_type;
public:
    Layer();
    Layer(unsigned long n_neurons);
    void setNNeurons(unsigned long n_neurons);
    unsigned long getNNeurons(){return number_of_neurons;};
    void sayhi();
    void setInput(bool def){input = def;};
    bool isInput(){return input;};
    void setValue(size_t i_neuron, su2double value){neurons[i_neuron].setValue(value);}
    su2double getValue(size_t i_neuron){return neurons[i_neuron].getValue();}
    void setBias(size_t i_neuron, su2double value){neurons[i_neuron].setBias(value);}
    su2double getBias(size_t i_neuron){return neurons[i_neuron].getBias();}
    su2double getdYdX(size_t i_neuron){return neurons[i_neuron].getGradient();}
    string getActivationType(){return activation_type;}
    ~Layer(){
        delete [] neurons;
    };
};
