/*!
 * \file CLookUp_ANN.cpp
 * \brief Implementation of the multi-layer perceptron class to be 
 *      used for look-up operations.
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
#include "../../../include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"
#include "../../../include/toolboxes/multilayer_perceptron/CIOMap.hpp"
#include "../../../include/toolboxes/multilayer_perceptron/CReadNeuralNetwork.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

using namespace std;

MLPToolbox::CLookUp_ANN::CLookUp_ANN(string inputFileName)
{
    
    #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif
    if(rank == MASTER_NODE)
        cout << "Generating ANN collection" << endl;

    ReadANNInputFile(inputFileName);
    for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++){
        
        NeuralNetworks.at(i_ANN) = new CNeuralNetwork;
        if(rank == MASTER_NODE)
            cout << "Generating neural network for " << ANN_filenames.at(i_ANN) << endl;
        GenerateANN(NeuralNetworks.at(i_ANN), ANN_filenames.at(i_ANN));
    }
}


vector<pair<size_t, size_t>> MLPToolbox::CLookUp_ANN::FindVariable_Indices(size_t i_ANN, vector<string> variable_names, bool input) const {
    vector<pair<size_t, size_t>> variable_indices;
    size_t nVar = input ? NeuralNetworks[i_ANN]->GetnInputs() : NeuralNetworks[i_ANN]->GetnOutputs();

    for(size_t iVar=0; iVar<nVar; iVar++){
        for(size_t jVar=0; jVar<variable_names.size(); jVar++){
            string ANN_varname = input ? NeuralNetworks[i_ANN]->GetInputName(iVar) : NeuralNetworks[i_ANN]->GetOutputName(iVar);
            
            if(variable_names.at(jVar).compare(ANN_varname) == 0){
                variable_indices.push_back(make_pair(jVar, iVar));
                //cout << variable_names[jVar] << " : " << jVar << " | " << ANN_varname << " : " << iVar << endl;
            }
        }
    }
    return variable_indices;
}

void MLPToolbox::CLookUp_ANN::Predict_ANN(CIOMap *input_output_map, vector<su2double>& inputs, vector<su2double*>& outputs){
    for(size_t i_map=0; i_map<input_output_map->GetNANNs(); i_map++){
        size_t i_ANN = input_output_map->GetANNIndex(i_map);
        vector<su2double> ANN_inputs = input_output_map->GetANN_Inputs(i_map, inputs);
        NeuralNetworks[i_ANN]->predict(ANN_inputs);
        for(size_t i=0; i < input_output_map->GetNMappedOutputs(i_map); i++){
            *outputs[input_output_map->GetOutputIndex(i_map, i)] = NeuralNetworks[i_ANN]->GetANN_Output(input_output_map->GetANNOutputIndex(i_map, i));
        }
    }
}
void MLPToolbox::CLookUp_ANN::GenerateANN(CNeuralNetwork * ANN, string fileName)
{
    CReadNeuralNetwork Reader = CReadNeuralNetwork(fileName);
    
    // Read MLP input file
    Reader.ReadMLPFile();

    // Generate basic layer architectures
    ANN->defineInputLayer(Reader.GetNInputs());
    ANN->sizeInputs(Reader.GetNInputs());
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

    // Size weights of each layer
    ANN->sizeWeights();

    // Define weights and activation functions
    ANN->SizeActivationFunctions(ANN->getNWeightLayers()+1);
    for(size_t i_layer = 0; i_layer < ANN->getNWeightLayers(); i_layer++){
        ANN->setActivationFunction(i_layer, Reader.GetActivationFunction(i_layer));
        for(size_t i_neuron=0; i_neuron < ANN->getNNeurons(i_layer); i_neuron++){
            for(size_t j_neuron=0; j_neuron<ANN->getNNeurons(i_layer, i_neuron); j_neuron++){
                ANN->setWeight(i_layer, i_neuron, j_neuron, Reader.GetWeight(i_layer, i_neuron, j_neuron));
            }
        }

        
    }
    ANN->setActivationFunction(ANN->getNWeightLayers(), Reader.GetActivationFunction(ANN->getNWeightLayers()));
    
    // Set neuron biases
    for(size_t i_layer=0; i_layer<ANN->getNWeightLayers()+1; i_layer++){
        for(size_t i_neuron=0; i_neuron<ANN->getNNeurons(i_layer); i_neuron++){
            ANN->setBias(i_layer, i_neuron, Reader.GetBias(i_layer, i_neuron));
        }
    }
    
    // Define input and output layer normalization values
    for(unsigned long iInput=0; iInput<Reader.GetNInputs(); iInput++){
        ANN->SetInputNorm(iInput, Reader.GetInputNorm(iInput).first, Reader.GetInputNorm(iInput).second);
    }
    for(unsigned long iOutput=0; iOutput<Reader.GetNOutputs(); iOutput++){
        ANN->SetOutputNorm(iOutput, Reader.GetOutputNorm(iOutput).first, Reader.GetOutputNorm(iOutput).second);
    }
    
}

void MLPToolbox::CLookUp_ANN::ReadANNInputFile(string inputFileName)
{
    ifstream file_stream;
    istringstream stream_names_var;
    istringstream stream_filenames;
    file_stream.open(inputFileName.c_str(), ifstream::in);

    if (!file_stream.is_open()) {
      SU2_MPI::Error(string("There is no MLP collection file file called ") + inputFileName,
                    CURRENT_FUNCTION);
    }

    string line;
    string word;

    line = SkipToFlag(&file_stream, "[number of MLP files]");
    getline(file_stream, line);
    number_of_variables = stoul(line);
    NeuralNetworks.resize(number_of_variables);

    line = SkipToFlag(&file_stream, "[MLP file names]");
    getline(file_stream, line);
    stream_filenames.str(line);
    while(stream_filenames){
        stream_filenames >> word;
        ANN_filenames.push_back(word);
    }
    ANN_filenames.pop_back();
}
string MLPToolbox::CLookUp_ANN::SkipToFlag(ifstream *file_stream, string flag) {
  string line;
  getline(*file_stream, line);

  while (line.compare(flag) != 0 && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof())
    SU2_MPI::Error("Flag not found in file", CURRENT_FUNCTION);

  return line;
}


bool MLPToolbox::CLookUp_ANN::Check_Duplicate_Outputs(vector<string> &output_names, CIOMap *input_output_map) const {
    unsigned short n_occurances;
    bool duplicate{false};
    vector<string> duplicate_variables;
    for(size_t i_Output =0; i_Output < output_names.size(); i_Output++){
        n_occurances = 0;
        for(size_t i_map=0; i_map<input_output_map->GetNANNs(); i_map++){
            vector<pair<size_t, size_t>> output_map = input_output_map->GetOutputMapping(i_map);
            for(size_t j_Output=0; j_Output<output_map.size(); j_Output++){
                if(output_map[j_Output].first == i_Output) n_occurances++;
            }
        }
        if(n_occurances > 1){
            duplicate_variables.push_back(output_names[i_Output]);
            duplicate = true;
        }
    }
    if(duplicate){
        string message{"Variables "};
        for(size_t iVar=0; iVar<duplicate_variables.size(); iVar++) message += duplicate_variables[iVar] + " ";
        SU2_MPI::Error(message + "occur more than once in the loaded ANN outputs.", CURRENT_FUNCTION);
    }
    
    return duplicate;
}


bool MLPToolbox::CLookUp_ANN::Check_Use_of_Inputs(vector<string> &input_names, CIOMap *input_output_map) const {
vector<string> missing_inputs;
bool inputs_are_present{true};
for(size_t iInput=0; iInput<input_names.size(); iInput ++){
    bool found_input = false;
    for(size_t i_map=0; i_map < input_output_map->GetNANNs(); i_map++){
        vector<pair<size_t, size_t>> input_map = input_output_map->GetInputMapping(i_map);
        for(size_t jInput=0; jInput<input_map.size(); jInput++){
            if(input_map[jInput].first == iInput){
                found_input = true;
            }
        }
    }
    if(!found_input){
        missing_inputs.push_back(input_names[iInput]);
        inputs_are_present = false;
    };
}
if(missing_inputs.size() > 0){
    string message{"Inputs "};
    for(size_t iVar=0; iVar<missing_inputs.size(); iVar++) message += missing_inputs[iVar] + " ";
    SU2_MPI::Error(message + "are not present in any loaded ANN.", CURRENT_FUNCTION);
}
return inputs_are_present;
}

bool MLPToolbox::CLookUp_ANN::Check_Use_of_Outputs(vector<string> &output_names, CIOMap * input_output_map) const {
    /* Check wether all output variables are in the loaded MLPs */
    
    vector<string> missing_outputs;     
    bool outputs_are_present{true};
    /* Looping over the target outputs */
    for(size_t iOutput=0; iOutput<output_names.size(); iOutput ++){

        bool found_output{false};

        /* Looping over all the selected ANNs */
        for(size_t i_map=0; i_map < input_output_map->GetNANNs(); i_map++){
            vector<pair<size_t, size_t>> output_map = input_output_map->GetOutputMapping(i_map);

            /* Looping over the outputs of the output map of the current ANN */
            for(size_t jOutput=0; jOutput<output_map.size(); jOutput++){
                if(output_map[jOutput].first == iOutput) found_output = true;
            }
        }
        if(!found_output){
            missing_outputs.push_back(output_names[iOutput]);
            outputs_are_present = false;
        };
    }
    if(missing_outputs.size() > 0){
        string message{"Outputs "};
        for(size_t iVar=0; iVar<missing_outputs.size(); iVar++) message += missing_outputs[iVar] + " ";
        SU2_MPI::Error(message + "are not present in any loaded ANN.", CURRENT_FUNCTION);
    }
    return outputs_are_present;
}