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
#include "../../../include/numerics/multilayer_perceptron/ReadNeuralNetwork.hpp"
using namespace std;

ReadNeuralNetwork::ReadNeuralNetwork(string filename_in){
    filename = filename_in;
};

void ReadNeuralNetwork::ReadMLPFile(){
    ifstream file_stream;
    file_stream.open(filename.c_str(), ifstream::in);
    if (!file_stream.is_open()) {
      SU2_MPI::Error(string("There is no MLP file file called ") + filename,
                    CURRENT_FUNCTION);
    }

    string line, word;
    double input_min, input_max, output_min, output_max;

    line = SkipToFlag(&file_stream, "[number of layers]");
    getline(file_stream, line);
    n_layers = stoul(line);
    n_neurons.resize(n_layers);
    biases.resize(n_layers);
    weights.resize(n_layers-1);
    activation_functions.resize(n_layers);
    line = SkipToFlag(&file_stream, "[neurons per layer]");
    for(size_t iLayer=0; iLayer<n_layers; iLayer++){
      getline(file_stream, line);
      n_neurons.at(iLayer) = stoul(line);
      biases.at(iLayer).resize(n_neurons.at(iLayer));
    }
    for(size_t iLayer=0; iLayer<n_layers-1; iLayer++){
      weights.at(iLayer).resize(n_neurons.at(iLayer));
      for(size_t iNeuron=0; iNeuron<n_neurons.at(iLayer); iNeuron++){
        weights.at(iLayer).at(iNeuron).resize(n_neurons.at(iLayer+1));
      }
    }

    input_norm.resize(n_neurons.at(0));
    output_norm.resize(n_neurons.at(n_neurons.size()-1));

    line = SkipToFlag(&file_stream, "[activation function]");
    for(size_t iLayer=0; iLayer<n_layers; iLayer++){
      getline(file_stream, line);
      istringstream activation_stream(line);
      activation_stream >> word;
      activation_functions.at(iLayer) = word;
    }    

    line = SkipToFlag(&file_stream, "[input names]");
    for(size_t iInput=0; iInput < n_neurons.at(0); iInput++){
      getline(file_stream, line);
      input_names.push_back(line);
    }

    line = SkipToFlag(&file_stream, "[input normalization]");
    for(size_t iInput=0; iInput<input_norm.size(); iInput++){
      getline(file_stream, line);
      istringstream input_norm_stream(line);
      input_norm_stream >> word;
      input_min = stold(word);
      input_norm_stream >> word;
      input_max = stold(word);
      input_norm.at(iInput) = make_pair(input_min, input_max);
    }

    line = SkipToFlag(&file_stream, "[output names]");
    for(size_t iOutput=0; iOutput < n_neurons.at(n_neurons.size()-1); iOutput++){
      getline(file_stream, line);
      output_names.push_back(line);
    }

    line = SkipToFlag(&file_stream, "[output normalization]");
    for(size_t iOutput=0; iOutput<output_norm.size(); iOutput++){
      getline(file_stream, line);
      istringstream output_norm_stream(line);
      output_norm_stream >> word;
      output_min = stold(word);
      output_norm_stream >> word;
      output_max = stold(word);
      output_norm.at(iOutput) = make_pair(output_min, output_max);
    }

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

string ReadNeuralNetwork::SkipToFlag(ifstream *file_stream, string flag) {
  string line;
  getline(*file_stream, line);

  while (line.compare(flag) != 0 && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof())
    cout << "line not in file!" << endl;

  return line;
}