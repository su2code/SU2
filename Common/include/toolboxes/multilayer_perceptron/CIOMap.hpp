/*!
 * \file CIOMap.hpp
 * \brief Declaration of input-to-output mapping class for artifical neural networks
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
#include "../../linear_algebra/blas_structure.hpp"
#include "CLookUp_ANN.hpp"

namespace MLPToolbox
{
        
    class CLookUp_ANN;

    class CIOMap
    /*! 
    *\class CIOMap
    *\brief This class is used by the CLookUp_ANN class to assign user-defined inputs 
    * and outputs to loaded multi-layer perceptrons. When a look-up operation is called
    * with a specific CIOMap, the multi-layer perceptrons are evaluated with input and
    * output variables coinciding with the desired input and output variable names.
    * 
    *
    * For example, in a custom, data-driven fluid model, MLP's are used for thermodynamic state
    * definition. There are three MLP's loaded. MLP_1 predicts temperature and specific heat 
    * based on density and energy. MLP_2 predicts pressure and speed of sound based on density and
    * energy as well. MLP_3 predicts density and energy based on pressure and temperature.
    * During a certain look-up operation in the CFluidModel, temperature, speed of sound and pressure 
    * are needed for a given density and energy. What the CIOMap does is to point to MLP_1 for 
    * temperature evalutation, and to MLP_2 for pressure and speed of sound evaluation. MLP_3 is 
    * not considered, as the respective inputs and outputs don't match with the function call
    * inputs and outputs.
    * 
    *  call variables:      MLP inputs:                     MLP outputs:           call outputs:
    *   
    *                        2--> energy --|            |--> temperature --> 1       
    *                                      |--> MLP_1 --|       
    *  1:density            1--> density --|            |--> c_p                   1:temperature
    *  2:energy                                                                    2:speed of sound
    *                       1--> density --|            |--> pressure --> 3        3:pressure
    *                                      |--> MLP_2 --|
    *                        2--> energy --|            |--> speed of sound --> 2
    * 
    *                           pressure --|            |--> density
    *                                      |--> MLP_3 --|
    *                        temperature --|            |--> energy
    * 
    * 
    * \author E.Bunschoten
    */

    {
        private:
        std::vector<std::size_t> MLP_indices;       /*!< Loaded MLP index */

        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> 
            Input_Map,      /*!< Mapping of call variable inputs to matching MLP inputs */
            Output_Map;     /*!< Mapping of call variable outputs to matching MLP outputs */
        
        public:

        /*!
        * \brief CIOMap class constructor;
        * \param[in] MLP_collection - Pointer to CLookUp_ANN class for which the input-output amap is to be created
        * \param[in] inputs - Input names for the call function. These should match with at least one of the MLP inputs.
        * \param[in] outputs - Output names for the call function. These should match with at least one of the MLP outputs.
        */
        CIOMap(CLookUp_ANN*MLP_collection, su2vector<std::string> &inputs, su2vector<std::string> &outputs);

        /*!
        * \brief Set MLP index in IO map
        * \param[in] iMLP - loaded MLP index
        */
        void Push_MLP_index(std::size_t iMLP){ MLP_indices.push_back(iMLP); }

        /*!
        * \brief Pair call inputs and outputs with the inputs and outputs of 
                the loaded MLPs
        * \param[in] MLP_collection - pointer to MLP collection class
        * \param[in] inputs - vector with call input variable names
        * \param[in] outputs - vector with call output variable names
        */
        void PairVariableswithMLPs(CLookUp_ANN * MLP_collection, su2vector<std::string> &inputs, su2vector<std::string> &outputs);

        /*!
        * \brief Get the number of MLPs in the current IO map
        * \return number of MLPs with matching inputs and output(s)
        */
        std::size_t GetNMLPs() const { return MLP_indices.size(); }

        /*!
        * \brief Get the loaded MLP index
        * \return MLP index
        */
        std::size_t GetMLPIndex(std::size_t i_Map) const { return MLP_indices[i_Map]; }

        /*!
        * \brief Get the call input variable index 
        * \param[in] i_Map - input-output mapping index of the IO map
        * \param[in] iInput - input index of the call input variable
        * \return MLP input variable index
        */
        std::size_t GetInputIndex(std::size_t i_Map, std::size_t iInput) const { return Input_Map[i_Map][iInput].first; }

        /*!
        * \brief Get the call output variable index 
        * \param[in] i_Map - input-output mapping index of the IO map
        * \param[in] iOutput - output index of the call input variable
        * \return call variable output index
        */
        std::size_t GetOutputIndex(std::size_t i_Map, std::size_t iOutput) const { return Output_Map[i_Map][iOutput].first; }

        /*!
        * \brief Get the MLP output variable index 
        * \param[in] i_Map - input-output mapping index of the IO map
        * \param[in] iOutput - output index of the call input variable
        * \return MLP output variable index
        */
        std::size_t GetMLPOutputIndex(std::size_t i_Map, std::size_t iOutput) const { return Output_Map[i_Map][iOutput].second; }

        /*!
        * \brief Get the number of matching output variables between the call and MLP outputs
        * \param[in] i_Map - input-output mapping index of the IO map
        * \return Number of matching variables between the loaded MLP and call variables
        */
        std::size_t GetNMappedOutputs(std::size_t i_Map) const { return Output_Map[i_Map].size(); }

        /*!
        * \brief Get the mapping of MLP outputs matching to call outputs
        * \param[in] i_Map - input-output mapping index of the IO map
        * \return Mapping of MLP output variables to call variables
        */
        std::vector<std::pair<std::size_t, std::size_t>> GetOutputMapping(std::size_t i_map) const { return Output_Map[i_map]; }

        /*!
        * \brief Get the mapping of MLP inputs to call inputs
        * \param[in] i_Map - input-output mapping index of the IO map
        * \return Mapping of MLP input variables to call inputs
        */
        std::vector<std::pair<std::size_t, std::size_t>> GetInputMapping(std::size_t i_map) const { return Input_Map[i_map]; }
        
        /*!
        * \brief Get the mapped inputs for the MLP at i_Map
        * \param[in] i_Map - input-output mapping index of the IO map
        * \param[in] inputs - call inputs
        * \return Vector with call inputs in the correct order of the loaded MLP
        */
        su2vector<su2double> GetMLP_Inputs(std::size_t i_Map, su2vector<su2double>&inputs) const {
            su2vector<su2double> MLP_input;
            MLP_input.resize(Input_Map[i_Map].size());

            for(std::size_t iInput=0; iInput<Input_Map[i_Map].size(); iInput++){
                MLP_input[iInput] = inputs[GetInputIndex(i_Map, iInput)];
            }
            return MLP_input;
        }
    };
}
