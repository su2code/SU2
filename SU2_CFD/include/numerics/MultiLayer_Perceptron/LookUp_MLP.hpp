/*!
 * \file LookUp_MLP.hpp
 * \brief Delarations of numerics class for a collection of MLP's used for look-up functionalities.
 * \version 7.3.1 "Blackbird"
 * \author E.C.Bunschoten
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
#include <string>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <vector>
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/linear_algebra/blas_structure.hpp"
#include "NeuralNetwork.hpp"

using namespace std;

class LookUp_MLP
/*!
 * \class LookUp_MLP
 * \brief 
 * \details This class contains all the attributes for the LookUp_MLP class, which allows for interpolation of parameters via MLP.
 * \author E.C.Bunschoten
 */
{
    private:

    int rank{0};                            /*!< \brief MPI rank */
    vector<NeuralNetwork*> NeuralNetworks;  /*!< \brief Vector containing MLP's for different variables. */
    vector<string> names_var;               /*!< \brief Vector containing the names for the MLP's */
    vector<string> ANN_filenames;           /*!< \brief Vector containing the .mlp filenames for the MLP's */
    vector<vector<string>> ANN_inputs;      /*!< \brief Vector containing the input names for each MLP */
    vector<vector<string>> ANN_outputs;     /*!< \brief Vector containing the output names for each MLP */
    unsigned long number_of_variables;      /*!< \brief Number of variables in the MLP collection */

    /*!
     *\brief Skip over to specific string in file.
     *\param[in] file_stream - file stream of the input file
     *\param[in] flag - String to skip to.
     *\param[out] line - Line at flag.
    */
    string SkipToFlag(ifstream *file_stream, string flag);

    /*!
    *\brief Get the index of a variable from the MLP list.
    *\param[in] nameVar - Variable name.
    *\param[out] index - Variable index.
    */
    inline int GetIndexOfVar(string nameVar) {
        int index;
        int endoflist;
        index =  (int)(find(names_var.begin(), names_var.end(), nameVar) - names_var.begin());
        endoflist = names_var.size();
        if (index == endoflist){
        index = -1;
        string error_msg = "Variable '";
        error_msg.append(nameVar);
        error_msg.append("' is not in the artificial neural network.");
        SU2_MPI::Error(error_msg, CURRENT_FUNCTION);
        }
        return index;
    }

    /*!
    *\brief Generate an MLP based on .mlp input file contents.
    *\param[in] ANN - Pointer to NeuralNetwork class to be modified.
    *\param[in] filename - .mlp file describing MLP architecture.
    */
    void GenerateANN(NeuralNetwork*, string);

    /*!
    *\brief Read .mlp input file.
    *\param[in] fileName - .mlp file.
    */
    void ReadANNInputFile(string fileName);

    /*!
    *\brief Get the index of the NeuralNetwork class with the corresponding inputs and outputs.
    *\param[in] input_names - Vector containing the desired input variable names.
    *\param[in] output_names - Vector containing the desired output names.
    *\param [out] idx - Index of the NeuralNetwork class with corresponding same input and output names.
    */
    size_t SelectANN(vector<string> input_names, vector<string>output_names);

    public: 
    /*!
    *\brief Class constructor
    */
    LookUp_MLP(string fileName);

    /*!
    *\brief Evaluate MLP within the MLP collection.
    *\param[in] input_names - Vector containing the desired input variable names.
    *\param[in] output_names - Vector containing the desired output variable names.
    *\param[in] output - Vector containing pointers to output variables.
    *\param[in] doutputs_dinputs - Pointer for the derivatives of the outputs with respect to the inputs.
    */
    void Predict_MLP(vector<string> input_names, vector<su2double> input, vector<string> output_names, vector<su2double*> output, su2double**doutputs_dinputs=nullptr);

    pair<su2double, su2double> GetInputNorm(string input_name);
    pair<su2double, su2double> GetOutputNorm(string output_name);
    
    /*!
    *\brief Class destructor
    */
    ~LookUp_MLP(){for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++)
        delete NeuralNetworks.at(i_ANN);
    };


};