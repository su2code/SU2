/*!
 * \file CLookUp_ANN.hpp
 * \brief Declaration of artificial neural network interpolation class
 * \author E. Bunschoten
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


namespace MLPToolbox{

class CLookUp_ANN
{
    /*! 
    *\class CLookUp_ANN
    *\brief This class allows for the evaluation of one or more multi-layer perceptrons in for example
    * thermodynamic state look-up operations. The multi-layer perceptrons are loaded in the order listed
    * in the MLP collection file. Each multi-layer perceptron is generated based on the 
    * architecture described in its respective input file. When evaluating the MLP collection, an
    * input-output map is used to find the correct MLP corresponding to the call function inputs and 
    * outputs.  
    */

    private:
    int rank{0};
    su2vector<CNeuralNetwork*> NeuralNetworks;  /*!< std::vector containing all loaded neural networks. */

    unsigned short number_of_variables;  /*!< Number of loaded ANNs. */

    /*!
    * \brief Load ANN architecture
    * \param[in] ANN - pointer to target NeuralNetwork class
    * \param[in] filename - filename containing ANN architecture information
    */
    void GenerateANN(CNeuralNetwork*ANN, std::string filename);
    
    public: 
    
    /*!
    * \brief ANN collection class constructor
    * \param[in] n_inputs - Number of MLP files to be loaded.
    * \param[in] input_filenames - String array containing MLP input file names. 
    */
    CLookUp_ANN(const unsigned short n_inputs, const std::string*input_filenames);

    /*!
    * \brief Evaluate loaded ANNs for given inputs and outputs
    * \param[in] input_output_map - input-output map coupling desired inputs and outputs to loaded ANNs
    * \param[in] inputs - input values
    * \param[in] outputs - pointers to output variables
    * \returns Within output normalization range.
    */
    unsigned long Predict_ANN(CIOMap *input_output_map, su2vector<su2double> &inputs, su2vector<su2double*> &outputs);

    ~CLookUp_ANN() { 
        for(std::size_t i_ANN=0; i_ANN<number_of_variables; i_ANN++)
            delete NeuralNetworks[i_ANN];
    };

    /*!
    * \brief Get number of loaded ANNs
    * \return number of loaded ANNs
    */
    std::size_t GetNANNs() const { return NeuralNetworks.size(); }

    /*!
    * \brief Check if all output variables are present in the loaded ANNs
    * \param[in] output_names - output variable names to check
    * \param[in] input_output_map - pointer to input-output map to be checked
    */
    bool Check_Use_of_Outputs(su2vector<std::string> &output_names, CIOMap *input_output_map) const;

    /*!
    * \brief Check if all input variables are present in the loaded ANNs
    * \param[in] input_names - input variable names to check
    * \param[in] input_output_map - pointer to input-output map to be checked
    */
    bool Check_Use_of_Inputs(su2vector<std::string> &input_names, CIOMap *input_output_map) const;

    /*!
    * \brief Map variable names to ANN inputs or outputs
    * \param[in] i_ANN - loaded ANN index
    * \param[in] variable_names - variable names to map to ANN inputs or outputs
    * \param[in] input - map to inputs (true) or outputs (false)
    */
    std::vector<pair<std::size_t, std::size_t>> FindVariable_Indices(std::size_t i_ANN, su2vector<std::string> variable_names, bool input) const;
    
};  

}

