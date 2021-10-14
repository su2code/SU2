/*!
 * \file ReadBFMInput.hpp
 * \brief Declaration and inlines of the classes used to read
 *        the blade geometry input file for body-force model
 *        simulations.
 * \author E.C.Bunschoten
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include <vector>
#include "../variables/CVariable.hpp"

class ReadBFMInput{
private:

  string version_input_file;                  // Blade geometric parameter input file version
  unsigned short n_blade_rows{1};             // Number of blade rows in input file
  vector<unsigned long> n_axial_points{};     // Number of axial data entries per blade row
  vector<unsigned long> n_radial_points{};    // Number of radial data entries per blade row
  vector<unsigned long> n_tangential_points{};// Number of tangential sections per blade row
  vector<unsigned short> n_blades{};          // Number of blades in each blade row
  vector<unsigned short> rotation_factor{};   // Rotation factor for each blade row
  vector<string> variable_names{};            // Blade geometric parameter names

  vector<vector<CVectorOfMatrix>> Geometric_Parameters;  // Geometric blade parameters read from input file
  vector<string> translated_names;                               // Variable name for each blade parameter
  vector<pair<unsigned short, unsigned short>> name_translation; // Pairing name with variable index
    
  /*!
   * \brief Size blade input data matrix according to header information
   */
  void AllocateMemory(){
    // Allocate memory from heap
    //Geometric_Parameters = new vector<vector<CVectorOfMatrix>>;
    
    // Size according to number of blade rows
    Geometric_Parameters.resize(n_blade_rows);
    
    // Looping through blade rows
    for(unsigned short i_row = 0; i_row<n_blade_rows; ++i_row){
      
      // Size according to number of geometric parameters
      Geometric_Parameters[i_row].resize(N_BFM_PARAMS);
      
      // Looping through parameters to size data structure according to number of spatial data entries
      for(unsigned short i_param=0; i_param<N_BFM_PARAMS; ++i_param){
        Geometric_Parameters[i_row][i_param].resize(n_tangential_points[i_row], n_radial_points[i_row], n_axial_points[i_row]);
      }
    }
  }

  /*!
   * \brief Read blade geometry parameter input file
   * \param [in] file_name - input file name
   */
  void ReadInputFile(string file_name);

  /*!
   * \brief Skip to flag in input file
   * \param [in] file_stream - pointer to file stream of input file
   * \param [in] flag - flag to which to skip
   * \returns line at flag
   */
  string SkipToFlag(ifstream &file_stream, string flag);

  /*!
   * \brief Set number of axial data entries at a blade row
   * \param [in] n_input - number of axial data entries in the blade row
   * \param [in] i_row - Blade row index
   */
  void SetNAxialPoints(unsigned long n_input, unsigned short i_row){n_axial_points[i_row] = n_input;}

  /*!
   * \brief Set number of radial data entries at a blade row
   * \param [in] n_input - number of radial data entries in the blade row
   * \param [in] i_row - Blade row index
   */
  void SetNRadialPoints(unsigned long n_input, unsigned short i_row){n_radial_points[i_row] = n_input;}

  /*!
   * \brief Set number of tangential data entries at a blade row
   * \param [in] n_input - number of tangential data entries in the blade row
   * \param [in] i_row - Blade row index
   */
  void SetNTangentialPoints(unsigned long n_input, unsigned short i_row){n_tangential_points[i_row] = n_input;}

  /*!
   * \brief Set number of blade rows
   * \param [in] n_input - number of blade rows
   */
  void SetNBladeRows(unsigned short n_input)
  {
    n_blade_rows = n_input;
  }
  
  void TranslateVariables();
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - pointer to config class
   * \param[in] file_input_name - Name of blade geometry input file.
   */
  ReadBFMInput(const CConfig *config, string file_input_name);

  /*!
   * \brief Destructor of the class.
   */
  ~ReadBFMInput() {};

  /*!
  * \brief Get number of blade rows provided in BFM input file
  * \returns Blade row count.
  */
  inline unsigned short GetNBladeRows() const {return n_blade_rows;}

  /*!
  * \brief Get number of data points in axial direction in BFM input file
  * \param[in] i_row - blade row index
  * \returns Number of data points in axial direction.
  */
  inline unsigned long GetNAxialPoints(unsigned short i_row) const {return n_axial_points[i_row];}

  /*!
  * \brief Get number of data points in radial direction in BFM input file
  * \param[in] i_row - blade row index
  * \returns Number of data points in radial direction.
  */
  inline unsigned long GetNRadialPoints(unsigned short i_row) const {return n_radial_points[i_row];}

  /*!
  * \brief Get number of data points in tangential direction in BFM input file
  * \param[in] i_row - blade row index
  * \returns Number of data points in tangential direction.
  */
  inline unsigned long GetNTangentialPoints(unsigned short i_row) const {return n_tangential_points[i_row];}

  /*!
  * \brief Get blade geometric parameter at a specific entry
  * \param[in] iRow - Blade row index
  * \param[in] iTang - Tangential section index
  * \param[in] iRad - Radial section index
  * \param[in] iAx - Axial point index
  * \param[in] iVar - Parameter index
  * \returns Blade geometric parameter value
  */
  inline su2double GetBFMParameter(unsigned short iRow, unsigned long iTang, unsigned long iRad, unsigned long iAx, unsigned short iVar) const {
    return Geometric_Parameters[iRow][iVar](iTang, iRad, iAx);
  }

  inline int GetIndexOfVar(string nameVar) {
    int index;
    int endoflist;
    index =  (int)(find(variable_names.begin(), variable_names.end(), nameVar) - variable_names.begin());
    endoflist = variable_names.size();
    if (index == endoflist){
      index = -1;
      string error_msg = "Variable '";
      error_msg.append(nameVar);
      error_msg.append("' is not in BFM input file");
      SU2_MPI::Error(error_msg, CURRENT_FUNCTION);
    }
    return index;
  }
  
  inline unsigned short GetBladeCount(unsigned short i_row) const {return n_blades[i_row];}
  inline unsigned short GetRotationFactor(unsigned short i_row) const {return rotation_factor[i_row];}

};
