/*!
 * \file LookUp_MLP.cpp
 * \brief Implementation of numerics classes for the evaluation of multiple Multi-Layer Perceptrons for regression purposes.
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
#include "../../../include/numerics/MultiLayer_Perceptron/LookUp_MLP.hpp"
#include "../../../include/numerics/MultiLayer_Perceptron/ReadNeuralNetwork.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
// #include "/home/evert/PhD/cppflow/include/cppflow/model.h"

using namespace std;

LookUp_MLP::LookUp_MLP(string inputFileName)
{
    
    // string model_dir;
    // unsigned short n_inputs, n_outputs;
    // cout << "Provide the directory to the model: ";
    // cin >> model_dir;

    //cppflow::model model("a");

    #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif
    if(rank == MASTER_NODE)
        cout << "Generating ANN collection" << endl;

    ReadANNInputFile(inputFileName);
    ANN_inputs.resize(number_of_variables);
    ANN_outputs.resize(number_of_variables);
    for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++){
        
        NeuralNetworks.at(i_ANN) = new NeuralNetwork;
        if(rank == MASTER_NODE)
            cout << "Generating neural network for " << ANN_filenames.at(i_ANN) << endl;
        GenerateANN(NeuralNetworks.at(i_ANN), ANN_filenames.at(i_ANN));
        ANN_inputs[i_ANN].resize(NeuralNetworks.at(i_ANN)->GetnInputs());
        ANN_outputs[i_ANN].resize(NeuralNetworks.at(i_ANN)->GetnOutputs());
    }
}

size_t LookUp_MLP::SelectANN(vector<string> input_names, vector<string> output_names){
    string total_inputs;
    string total_outputs;
    for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++){
        bool input_match{true}, output_match{true};
        if(input_names.size() != NeuralNetworks[i_ANN]->GetnInputs()) continue;
        if(output_names.size() != NeuralNetworks[i_ANN]->GetnOutputs()) continue;

        for(size_t iInput=0; iInput<NeuralNetworks[i_ANN]->GetnInputs(); iInput++){
            if(input_names.at(iInput).compare(NeuralNetworks[i_ANN]->GetInputName(iInput)) != 0){
                input_match = false;
                //cout << "Input error" << endl;
            }
        }
        for(size_t iOutput=0; iOutput<NeuralNetworks[i_ANN]->GetnOutputs(); iOutput++){
            if(output_names.at(iOutput).compare(NeuralNetworks[i_ANN]->GetOutputName(iOutput)) != 0){
                output_match = false;
                //cout << "Output error" << endl;
            }
        }
        
        if(input_match && output_match) return i_ANN;
        
    }
    
    for(size_t iInput=0; iInput<input_names.size(); iInput++){
        total_inputs += input_names.at(iInput) + " ";
    }
    for(size_t iOutput=0; iOutput<output_names.size(); iOutput++){
        total_outputs += output_names.at(iOutput) + " ";
    }
    SU2_MPI::Error(string("No ANN found with inputs ") + total_inputs + string("and outputs ") + total_outputs,
                   CURRENT_FUNCTION);
}

void LookUp_MLP::Predict_MLP(vector<string> input_names, vector<su2double> inputs, vector<string> output_names, vector<su2double*>outputs, su2double**doutputs_dinputs){
    size_t i_ANN = SelectANN(input_names, output_names);
    NeuralNetworks[i_ANN]->predict(inputs, outputs, doutputs_dinputs);
}
void LookUp_MLP::GenerateANN(NeuralNetwork * ANN, string fileName)
{
    ReadNeuralNetwork Reader = ReadNeuralNetwork(fileName);
    
    Reader.ReadRawDRG();
    ANN->defineInputLayer(Reader.GetNInputs());
    for(size_t iInput=0; iInput<Reader.GetNInputs(); iInput++){
        ANN->PushInputName(Reader.GetInputName(iInput));
    }
    for(size_t iLayer=1; iLayer<Reader.GetNlayers()-1; iLayer++){
        ANN->push_hidden_layer(Reader.GetNneurons(iLayer));
    }
    ANN->defineOutputLayer(Reader.GetNOutputs());
    for(size_t iOutput=0; iOutput<Reader.GetNOutputs(); iOutput++){
        ANN->PushOutputName(Reader.GetOutputName(iOutput));
    }
    ANN->sizeWeights();
    
    for(size_t i_layer = 0; i_layer < ANN->getNWeightLayers(); i_layer++){
        ANN->setActivationFunction(i_layer, Reader.GetActivationFunction(i_layer));
        for(size_t i_neuron=0; i_neuron < ANN->getNNeurons(i_layer); i_neuron++){
            for(size_t j_neuron=0; j_neuron<ANN->getNNeurons(i_layer+1); j_neuron++){
                ANN->setWeight(i_layer, i_neuron, j_neuron, Reader.GetWeight(i_layer, i_neuron, j_neuron));
            }
        }   
        
    }
    ANN->setActivationFunction(ANN->getNWeightLayers(), Reader.GetActivationFunction(ANN->getNWeightLayers()));
    
    for(size_t i_layer=0; i_layer<ANN->getNWeightLayers()+1; i_layer++){
        for(size_t i_neuron=0; i_neuron<ANN->getNNeurons(i_layer); i_neuron++){
            ANN->setBias(i_layer, i_neuron, Reader.GetBias(i_layer, i_neuron));
        }
    }
    
    for(unsigned long iInput=0; iInput<Reader.GetNInputs(); iInput++){
        ANN->SetInputNorm(iInput, Reader.GetInputNorm(iInput).first, Reader.GetInputNorm(iInput).second);
    }

    for(unsigned long iOutput=0; iOutput<Reader.GetNOutputs(); iOutput++){
        ANN->SetOutputNorm(iOutput, Reader.GetOutputNorm(iOutput).first, Reader.GetOutputNorm(iOutput).second);
    }
    
}
void LookUp_MLP::ReadANNInputFile(string inputFileName)
{
    ifstream file_stream;
    istringstream stream_names_var;
    istringstream stream_filenames;
    file_stream.open(inputFileName.c_str(), ifstream::in);
    string line;
    string word;

    line = SkipToFlag(&file_stream, "[number of variables]");
    getline(file_stream, line);
    number_of_variables = stoul(line);
    NeuralNetworks.resize(number_of_variables);

    line = SkipToFlag(&file_stream, "[variable names]");
    getline(file_stream, line);
    stream_names_var.str(line);
    while (stream_names_var) {
        stream_names_var >> word;
        int ixColon = (int)word.find(":");

        names_var.push_back(word.substr(ixColon + 1, word.size() - 1));
    }
    names_var.pop_back();  // removes last redundant element


    line = SkipToFlag(&file_stream, "[MLP file names]");
    getline(file_stream, line);
    stream_filenames.str(line);
    while(stream_filenames){
        stream_filenames >> word;
        ANN_filenames.push_back(word);
    }
    ANN_filenames.pop_back();
}
string LookUp_MLP::SkipToFlag(ifstream *file_stream, string flag) {
  string line;
  getline(*file_stream, line);

  while (line.compare(flag) != 0 && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof())
    SU2_MPI::Error("Flag not found in file", CURRENT_FUNCTION);

  return line;
}

pair<su2double, su2double> LookUp_MLP::GetOutputNorm(string outputName){
    for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++){
        for(size_t iOutput=0; iOutput<NeuralNetworks[i_ANN]->GetnOutputs(); iOutput++){
            if(outputName.compare(NeuralNetworks.at(i_ANN)->GetOutputName(iOutput)) == 0){
                return NeuralNetworks.at(i_ANN)->GetOutputNorm(iOutput);
            }
        }
    }
    SU2_MPI::Error("No MLP with output "+outputName, CURRENT_FUNCTION);
    return make_pair(0, 1);
}

pair<su2double, su2double> LookUp_MLP::GetInputNorm(string inputName){
    for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++){
        for(size_t iInput=0; iInput<NeuralNetworks[i_ANN]->GetnInputs(); iInput++){
            if(inputName.compare(NeuralNetworks.at(i_ANN)->GetInputName(iInput)) == 0){
                return NeuralNetworks.at(i_ANN)->GetInputNorm(iInput);
            }
        }
    }
    SU2_MPI::Error("No MLP with input "+inputName, CURRENT_FUNCTION);
    return make_pair(0, 1);
}