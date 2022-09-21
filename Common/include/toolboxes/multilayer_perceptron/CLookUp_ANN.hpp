/*!
 * \file CLookUp_ANN.hpp
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
#include <string>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <vector>
#include "../../CConfig.hpp"
#include "../../linear_algebra/blas_structure.hpp"
#include "CNeuralNetwork.hpp"
#include "CIOMap.hpp"

using namespace std;

namespace MLPToolbox{
class CIOMap;

class CLookUp_ANN
{
    private:
    int rank{0};
    vector<CNeuralNetwork*> NeuralNetworks;  /*!< Vector containing all loaded neural networks. */

    vector<string> ANN_filenames;       /*!< Filenames containing ANN architecture information. */

    unsigned long number_of_variables;  /*!< Number of loaded ANNs. */

    /*!
   * \brief Skip to certain flag in file.
   * \param[in] file_stream - input stream of file
   * \param[in] flag - line to search for
   * \return line in file.
   */
    string SkipToFlag(ifstream *file_stream, string flag); 

    /*!
    * \brief Load ANN architecture
    * \param[in] ANN - pointer to target NeuralNetwork class
    * \param[in] filename - filename containing ANN architecture information
    */
    void GenerateANN(CNeuralNetwork*ANN, string filename);

    /*!
    * \brief Read ANN architecture input file
    * \param[in] filename - filename containing ANN architecture information
    */
    void ReadANNInputFile(string fileName); 
    
    public: 
    
    /*!
    * \brief ANN collection class constructor
    * \param[in] filename - filename containing list of ANN input files
    */
    CLookUp_ANN(string fileName);
    
    /*!
    * \brief Evaluate loaded ANNs for given inputs and outputs
    * \param[in] input_output_map - input-output map coupling desired inputs and outputs to loaded ANNs
    * \param[in] inputs - input values
    * \param[in] outputs - pointers to output variables
    */
    void Predict_ANN(CIOMap *input_output_map, vector<su2double> &inputs, vector<su2double*> &outputs);

    ~CLookUp_ANN(){for(size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++)
        delete NeuralNetworks.at(i_ANN);
    };

    /*!
    * \brief Get number of loaded ANNs
    * \return number of loaded ANNs
    */
    size_t GetNANNs() const {return NeuralNetworks.size();}

    /*!
    * \brief Check if variables occur more than once in ANN outputs
    * \param[in] output_names - output variable names to check
    * \param[in] input_output_map - pointer to input-output map to be checked
    */
    bool Check_Duplicate_Outputs(vector<string> &output_names, CIOMap *input_output_map) const;

    /*!
    * \brief Check if all output variables are present in the loaded ANNs
    * \param[in] output_names - output variable names to check
    * \param[in] input_output_map - pointer to input-output map to be checked
    */
    bool Check_Use_of_Outputs(vector<string> &output_names, CIOMap *input_output_map) const;

    /*!
    * \brief Check if all input variables are present in the loaded ANNs
    * \param[in] input_names - input variable names to check
    * \param[in] input_output_map - pointer to input-output map to be checked
    */
    bool Check_Use_of_Inputs(vector<string> &input_names, CIOMap *input_output_map) const;

    /*!
    * \brief Map variable names to ANN inputs or outputs
    * \param[in] i_ANN - loaded ANN index
    * \param[in] variable_names - variable names to map to ANN inputs or outputs
    * \param[in] input - map to inputs (true) or outputs (false)
    */
    vector<pair<size_t, size_t>> FindVariable_Indices(size_t i_ANN, vector<string> variable_names, bool input) const;
    
};  

}

