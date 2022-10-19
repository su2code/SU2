/*!
 * \file ReadNeuralNetwork.cpp
 * \brief Implementation of the reader class to read .mlp input files
 *      used to set up multi-layer perceptrons.
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
#include "../../../include/toolboxes/multilayer_perceptron/CReadNeuralNetwork.hpp"
using namespace std;

MLPToolbox::CReadNeuralNetwork::CReadNeuralNetwork(string filename_in){
    filename = filename_in;
}

void MLPToolbox::CReadNeuralNetwork::ReadMLPFile(){
    ifstream file_stream;
    file_stream.open(filename.c_str(), ifstream::in);
    if (!file_stream.is_open()) {
      SU2_MPI::Error(string("There is no MLP file called ") + filename,
                    CURRENT_FUNCTION);
    }

    string line, word;
    double input_min, input_max, output_min, output_max;
    bool eoHeader{false}, found_layercount{false}, found_input_names{false}, found_output_names{false};

    /* Read general architecture information from file header */

    line = SkipToFlag(&file_stream, "<header>");

    while(getline(file_stream, line) && !eoHeader){

      /* Read layer count */
      if(line.compare("[number of layers]") == 0){
        getline(file_stream, line);
        n_layers = stoul(line);
        n_neurons.resize(n_layers);
        biases.resize(n_layers);
        weights.resize(n_layers-1);
        activation_functions.resize(n_layers);

        found_layercount = true;
      }

      /* Set number of neurons for each layer */
      if(line.compare("[neurons per layer]") == 0){

        /* In case layer count was not yet provided, return an error */
        if(!found_layercount){
          SU2_MPI::Error("No layer count provided before defining neuron count per layer", CURRENT_FUNCTION);
        }
        // Loop over layer count and size neuron count and bias count per layer accordingly
        for(size_t iLayer=0; iLayer<n_layers; iLayer++){
          getline(file_stream, line);
          n_neurons.at(iLayer) = stoul(line);
          biases.at(iLayer).resize(n_neurons.at(iLayer));
        }
        // Loop over spaces between layers and size the weight matrices accordingly
        for(size_t iLayer=0; iLayer<n_layers-1; iLayer++){
          weights.at(iLayer).resize(n_neurons.at(iLayer));
          for(size_t iNeuron=0; iNeuron<n_neurons.at(iLayer); iNeuron++){
            weights.at(iLayer).at(iNeuron).resize(n_neurons.at(iLayer+1));
          }
        }
        // Size input normalization according to number of inputs
        input_norm.resize(n_neurons.at(0));
        // Set default input normalization
        for(size_t iNeuron=0; iNeuron<n_neurons.at(0); iNeuron++) input_norm.at(iNeuron) = make_pair(0, 1);

        // Size output normalization according to number of outputs
        output_norm.resize(n_neurons.at(n_neurons.size()-1));
        // Set default output normalization
        for(size_t iNeuron=0; iNeuron<n_neurons.at(n_neurons.size()-1); iNeuron++) output_norm.at(iNeuron) = make_pair(0, 1);
      }

      /* Read layer activation function types */
      if(line.compare("[activation function]") == 0){
        if(!found_layercount){
          SU2_MPI::Error("No layer count provided before providing layer activation functions", CURRENT_FUNCTION);
        }
        for(size_t iLayer=0; iLayer<n_layers; iLayer++){
          getline(file_stream, line);
          istringstream activation_stream(line);
          activation_stream >> word;
          activation_functions.at(iLayer) = word;
        }    
      }

      /* Read MLP input variable names */
      if(line.compare("[input names]") == 0){
        found_input_names = true;
        getline(file_stream, line);
        while(line.compare("") != 0){
          input_names.push_back(line);
          getline(file_stream, line);
        }

        if(input_names.size() != n_neurons.at(0)){
          SU2_MPI::Error("Number of input variable names inconsistent with number of MLP inputs", CURRENT_FUNCTION);
        }
      }

      /* In case input normalization is applied, read upper and lower input bounds */
      if(line.compare("[input normalization]") == 0){
        for(size_t iInput=0; iInput<input_norm.size(); iInput++){
          getline(file_stream, line);
          if(line.compare("") != 0){
            istringstream input_norm_stream(line);
            input_norm_stream >> word;
            input_min = stold(word);
            input_norm_stream >> word;
            input_max = stold(word);
            input_norm.at(iInput) = make_pair(input_min, input_max);
          }
        }
      }

      /* Read MLP output variable names */
      if(line.compare("[output names]") == 0){
        found_output_names = true;
        getline(file_stream, line);
        while(line.compare("") != 0){
          output_names.push_back(line);
          getline(file_stream, line);
        }

        if(output_names.size() != (n_neurons.at(n_neurons.size()-1))){
          SU2_MPI::Error("Number of output variable names inconsistent with number of MLP outputs", CURRENT_FUNCTION);
        }
      }

      /* In case output normalization is applied, read upper and lower output bounds */
      if(line.compare("[output normalization]") == 0){
        for(size_t iOutput=0; iOutput<output_norm.size(); iOutput++){
          getline(file_stream, line);
          if(line.compare("") != 0){
            istringstream output_norm_stream(line);
            output_norm_stream >> word;
            output_min = stold(word);
            output_norm_stream >> word;
            output_max = stold(word);
            output_norm.at(iOutput) = make_pair(output_min, output_max);
          }
        }
      }

      if(line.compare("</header>") == 0){
        eoHeader = true;
      }
    }
    if(!found_input_names){
      SU2_MPI::Error("No MLP input variable names provided", CURRENT_FUNCTION);
    }
    if(!found_output_names){
      SU2_MPI::Error("No MLP output variable names provided", CURRENT_FUNCTION);
    }

    /* Read weights for each layer */
    line = SkipToFlag(&file_stream, "[weights per layer]");
    for(size_t iLayer=0; iLayer<n_layers-1; iLayer++){
      getline(file_stream, line);
      for(size_t iNeuron=0; iNeuron<n_neurons.at(iLayer); iNeuron++){
        getline(file_stream, line);
        istringstream weight_stream(line);
        for(size_t jNeuron=0; jNeuron<n_neurons.at(iLayer+1); jNeuron++){
          weight_stream >> word;
          weights.at(iLayer).at(iNeuron).at(jNeuron) = stold(word);
        }
      }
      getline(file_stream, line);
    }
    
    /* Read biases for each neuron */
    line = SkipToFlag(&file_stream, "[biases per layer]");
    for(size_t iLayer=0; iLayer<n_layers; iLayer++){
      getline(file_stream, line);
      istringstream bias_stream(line);
      for(size_t iNeuron=0; iNeuron<n_neurons.at(iLayer); iNeuron++){
        bias_stream >> word;
        biases.at(iLayer).at(iNeuron) = stold(word);
      }
    }



}

string MLPToolbox::CReadNeuralNetwork::SkipToFlag(ifstream *file_stream, string flag) {
  string line;
  getline(*file_stream, line);

  while (line.compare(flag) != 0 && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof())
    cout << "line not in file!" << endl;

  return line;
}