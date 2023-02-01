/*!
 * \file CNeuralNetwork.cpp
 * \brief Implementation of the NeuralNetwork class to be used
 *      for evaluation of multi-layer perceptrons.
 * \author E.C.Bunschoten
 * \version 7.5.0 "Blackbird"
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
#include "../../../include/toolboxes/multilayer_perceptron/CNeuralNetwork.hpp"

#include <iomanip>
#include <iostream>

#include "../../../include/toolboxes/multilayer_perceptron/CLayer.hpp"
#include "../../../include/toolboxes/multilayer_perceptron/CReadNeuralNetwork.hpp"
using namespace std;

void MLPToolbox::CNeuralNetwork::predict(su2vector<su2double>& inputs, bool compute_gradient) {
  /*--- Evaluate MLP for given inputs ---*/

  su2double y = 0, dy_dx = 0;  // Activation function output.
  bool same_point = true;
  /* Normalize input and check if inputs are the same w.r.t last evaluation */
  for (auto iNeuron = 0u; iNeuron < inputLayer->getNNeurons(); iNeuron++) {
    su2double x_norm =
        (inputs[iNeuron] - input_norm[iNeuron].first) / (input_norm[iNeuron].second - input_norm[iNeuron].first);
    if (abs(x_norm - inputLayer->getOutput(iNeuron)) > 0) same_point = false;
    inputLayer->setOutput(iNeuron, x_norm);

    if (compute_gradient)
      inputLayer->setdYdX(iNeuron, iNeuron, 1 / (input_norm[iNeuron].second - input_norm[iNeuron].first));
  }
  /* Skip evaluation process if current point is the same as during the previous evaluation */
  if (!same_point) {
    su2double alpha = 1.67326324;
    su2double lambda = 1.05070098;
    /* Traverse MLP and compute inputs and outputs for the neurons in each layer */
    for (auto iLayer = 1u; iLayer < n_hidden_layers + 2; iLayer++) {
      auto nNeurons_current = total_layers[iLayer]->getNNeurons();  // Neuron count of current layer
      su2double x;                                                  // Neuron input value

      /* Compute and store input value for each neuron */
      for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
        x = ComputeX(iLayer, iNeuron);
        total_layers[iLayer]->setInput(iNeuron, x);
        if (compute_gradient) {
          for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
            dy_dx = ComputedOutputdInput(iLayer, iNeuron, iInput);
            total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx);
          }
        }
      }

      /* Compute and store neuron output based on activation function */
      switch (activation_function_types[iLayer]) {
        case ENUM_ACTIVATION_FUNCTION::ELU:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            if (x > 0) {
              y = x;
              if (compute_gradient) dy_dx = 1.0;
            } else {
              y = exp(x) - 1;
              if (compute_gradient) dy_dx = exp(x);
            }
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::LINEAR:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            y = x;
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = 1.0;
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::EXPONENTIAL:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            y = exp(x);
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = 1.0;
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::RELU:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            if (x > 0) {
              y = x;
              if (compute_gradient) dy_dx = 1.0;
            } else {
              y = 0.0;
              if (compute_gradient) dy_dx = 0.0;
            }
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::SWISH:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            y = x / (1 + exp(-x));
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = exp(x) * (x + exp(x) + 1) / pow(exp(x) + 1, 2);
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::TANH:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            y = tanh(x);
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = 1 / pow(cosh(x), 2);
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::SIGMOID:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            y = 1.0 / (1 + exp(-x));
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = exp(-x) / pow(exp(-x) + 1, 2);
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::SELU:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);
            if (x > 0.0) {
              y = lambda * x;
              if (compute_gradient) dy_dx = lambda;
            } else {
              y = lambda * alpha * (exp(x) - 1);
              if (compute_gradient) dy_dx = lambda * alpha * exp(x);
            }
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::GELU:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            x = total_layers[iLayer]->getInput(iNeuron);

            y = 0.5 * x * (1 + tanh(sqrt(2 / PI_NUMBER) * (x + 0.044715 * pow(x, 3))));
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = 0.5 * (tanh(0.0356774 * pow(x, 3) + 0.797885 * x) + (0.107032 * pow(x, 3) + 0.797885 * x) *
                                                                              pow(cosh(x), -2) *
                                                                              (0.0356774 * pow(x, 3) + 0.797885 * x));
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        case ENUM_ACTIVATION_FUNCTION::NONE:
          for (auto iNeuron = 0u; iNeuron < nNeurons_current; iNeuron++) {
            y = 0.0;
            total_layers[iLayer]->setOutput(iNeuron, y);
            if (compute_gradient) {
              dy_dx = 0.0;
              for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
                total_layers[iLayer]->setdYdX(iNeuron, iInput, dy_dx * total_layers[iLayer]->getdYdX(iNeuron, iInput));
              }
            }
          }
          break;
        default:
          break;
      }  // activation_function_types
    }
  }
  /* Compute and de-normalize MLP output */
  for (auto iNeuron = 0u; iNeuron < outputLayer->getNNeurons(); iNeuron++) {
    su2double y_norm = outputLayer->getOutput(iNeuron);
    y = y_norm * (output_norm[iNeuron].second - output_norm[iNeuron].first) + output_norm[iNeuron].first;
    if (compute_gradient) {
      dy_dx = (output_norm[iNeuron].second - output_norm[iNeuron].first);
      for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++) {
        outputLayer->setdYdX(iNeuron, iInput, dy_dx * outputLayer->getdYdX(iNeuron, iInput));
        dOutputs_dInputs[iNeuron][iInput] = outputLayer->getdYdX(iNeuron, iInput);
      }
    }
    /* Storing output value */
    ANN_outputs[iNeuron] = y;
  }
}

void MLPToolbox::CNeuralNetwork::defineInputLayer(unsigned long n_neurons) {
  /*--- Define the input layer of the network ---*/
  inputLayer = new CLayer(n_neurons);

  /* Mark layer as input layer */
  inputLayer->setInput(true);
  input_norm.resize(n_neurons);
  input_names.resize(n_neurons);
}

void MLPToolbox::CNeuralNetwork::defineOutputLayer(unsigned long n_neurons) {
  /*--- Define the output layer of the network ---*/
  outputLayer = new CLayer(n_neurons);
  output_norm.resize(n_neurons);
  output_names.resize(n_neurons);
}

void MLPToolbox::CNeuralNetwork::push_hidden_layer(unsigned long n_neurons) {
  /*--- Add a hidden layer to the network ---*/
  CLayer* newLayer = new CLayer(n_neurons);
  hiddenLayers.push_back(newLayer);
  n_hidden_layers++;
}

void MLPToolbox::CNeuralNetwork::sizeWeights() {
  /*--- Size weight matrices based on neuron counts in each layer ---*/

  /* Generate vector containing input, output, and hidden layer references */
  total_layers.resize(n_hidden_layers + 2);
  total_layers[0] = inputLayer;
  for (auto iLayer = 0u; iLayer < n_hidden_layers; iLayer++) {
    total_layers[iLayer + 1] = hiddenLayers[iLayer];
  }
  total_layers[total_layers.size() - 1] = outputLayer;

  weights_mat.resize(n_hidden_layers + 1);
  weights_mat[0].resize(hiddenLayers[0]->getNNeurons(), inputLayer->getNNeurons());
  for (auto iLayer = 1u; iLayer < n_hidden_layers; iLayer++) {
    weights_mat[iLayer].resize(hiddenLayers[iLayer]->getNNeurons(), hiddenLayers[iLayer - 1]->getNNeurons());
  }
  weights_mat[n_hidden_layers].resize(outputLayer->getNNeurons(), hiddenLayers[n_hidden_layers - 1]->getNNeurons());

  ANN_outputs = new su2double[outputLayer->getNNeurons()];
  dOutputs_dInputs.resize(outputLayer->getNNeurons(), inputLayer->getNNeurons());
  for (auto iLayer = 0u; iLayer < n_hidden_layers + 2; iLayer++) {
    total_layers[iLayer]->sizeGradients(inputLayer->getNNeurons());
  }
}

void MLPToolbox::CNeuralNetwork::displayNetwork() const {
  /*--- Display information on the MLP architecture ---*/
  int display_width = 54;
  int column_width = int(display_width / 3.0) - 1;

  /*--- Input layer information ---*/
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  cout << "|" << left << setw(display_width - 1) << "Input Layer Information:"
       << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  cout << "|" << left << setw(column_width) << "Input Variable:"
       << "|" << left << setw(column_width) << "Lower limit:"
       << "|" << left << setw(column_width) << "Upper limit:"
       << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');

  /*--- Hidden layer information ---*/
  for (auto iInput = 0u; iInput < inputLayer->getNNeurons(); iInput++)
    cout << "|" << left << setw(column_width) << to_string(iInput + 1) + ": " + input_names[iInput] << "|" << right
         << setw(column_width) << input_norm[iInput].first << "|" << right << setw(column_width)
         << input_norm[iInput].second << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  cout << "|" << left << setw(display_width - 1) << "Hidden Layers Information:"
       << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  cout << "|" << setw(column_width) << left << "Layer index"
       << "|" << setw(column_width) << left << "Neuron count"
       << "|" << setw(column_width) << left << "Function"
       << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  for (auto iLayer = 0u; iLayer < n_hidden_layers; iLayer++)
    cout << "|" << setw(column_width) << right << iLayer + 1 << "|" << setw(column_width) << right
         << hiddenLayers[iLayer]->getNNeurons() << "|" << setw(column_width) << right
         << activation_function_names[iLayer + 1] << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');

  /*--- Output layer information ---*/
  cout << "|" << left << setw(display_width - 1) << "Output Layer Information:"
       << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  cout << "|" << left << setw(column_width) << "Output Variable:"
       << "|" << left << setw(column_width) << "Lower limit:"
       << "|" << left << setw(column_width) << "Upper limit:"
       << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  for (auto iOutput = 0u; iOutput < outputLayer->getNNeurons(); iOutput++)
    cout << "|" << left << setw(column_width) << to_string(iOutput + 1) + ": " + output_names[iOutput] << "|" << right
         << setw(column_width) << output_norm[iOutput].first << "|" << right << setw(column_width)
         << output_norm[iOutput].second << "|" << endl;
  cout << "+" << setfill('-') << setw(display_width) << right << "+" << endl;
  cout << setfill(' ');
  cout << endl;
}

void MLPToolbox::CNeuralNetwork::setActivationFunction(unsigned long i_layer, string input) {
  /*--- Translate activation function name from input file to a number ---*/
  activation_function_names[i_layer] = input;
  if (input.compare("linear") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::LINEAR;
    return;
  }
  if (input.compare("relu") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::RELU;
    return;
  }
  if ((input.compare("elu") == 0)) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::ELU;
    return;
  }
  if (input.compare("swish") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SWISH;
    return;
  }
  if (input.compare("sigmoid") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SIGMOID;
    return;
  }
  if (input.compare("tanh") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::TANH;
    return;
  }
  if (input.compare("selu") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SELU;
    return;
  }
  if (input.compare("gelu") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::GELU;
    return;
  }
  if (input.compare("none") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::NONE;
    return;
  }
  if (input.compare("exponential") == 0) {
    activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::EXPONENTIAL;
    return;
  }
  SU2_MPI::Error("Unknown activation function type: " + input, CURRENT_FUNCTION);
  return;
}
