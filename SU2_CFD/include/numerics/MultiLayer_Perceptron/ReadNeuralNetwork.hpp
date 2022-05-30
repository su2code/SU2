#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/linear_algebra/blas_structure.hpp"

using namespace std;
class ReadNeuralNetwork {
    private:
    vector<string> input_names;
    vector<string> output_names;

    string filename;
    unsigned long n_layers;
    vector<unsigned long> n_neurons;
    vector<vector<vector<double long>>> weights;
    vector<vector<double long>> biases;
    vector<string> activation_functions;

    vector<pair<double long, double long>> input_norm;
    vector<pair<double long, double long>> output_norm;
    public:

    ReadNeuralNetwork(string filename_in);
    void ReadRawDRG();

    string SkipToFlag(ifstream *file_stream, string flag);
    
    unsigned long GetNInputs(){return n_neurons.at(0);}
    unsigned long GetNOutputs(){return n_neurons.at(n_layers - 1);}

    unsigned long GetNlayers(){return n_layers;}
    unsigned long GetNneurons(size_t iLayer){return n_neurons.at(iLayer);}
    double long GetWeight(size_t iLayer, size_t iNeuron, size_t jNeuron){return weights.at(iLayer).at(iNeuron).at(jNeuron);}
    double long GetBias(size_t iLayer, size_t iNeuron){return biases.at(iLayer).at(iNeuron);}
    pair<double long, double long> GetInputNorm(size_t iInput){return input_norm.at(iInput);}
    pair<double long, double long> GetOutputNorm(size_t iOutput){return output_norm.at(iOutput);}
    string GetActivationFunction(size_t iLayer){return activation_functions.at(iLayer);}

    string GetInputName(size_t iInput){return input_names.at(iInput);}
    string GetOutputName(size_t iOutput){return output_names.at(iOutput);}
    
    ~ReadNeuralNetwork(){};
};