/*!
 * \file LookUp_ANN.hpp
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

#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/linear_algebra/blas_structure.hpp"

using namespace std;
class Neuron
{
private:
    string activation_type="SmoothSlope";
    unsigned long ActivationFunction{ENUM_ACTIVATION_FUNCTION::LINEAR};
    unsigned long i_neuron;
    double value;
    su2double dy_dx;
    double bias{0};
    enum ENUM_ACTIVATION_FUNCTION {
        NONE=0,
        LINEAR=1,
        RELU=2,
        SMOOTH_SLOPE=3
    };

public:
    Neuron(){};
    void setNumber(unsigned long input){i_neuron = input;}
    unsigned long getNumber(){return i_neuron;}
    //void setFunctionType(string input);

    //string getFunctionType(){return activation_type;}
    //su2double activation_function(su2double x);
    void setValue(su2double input){value = input;}
    su2double getValue(){return value;}

    void setBias(su2double input){bias = input;}
    su2double getBias(){return bias;}

    su2double getGradient(){return dy_dx;}
    void setGradient(su2double input){dy_dx = input;}
    //void SetActivationFunction(unsigned long input){ActivationFunction = input;}

    ~Neuron(){};
};

