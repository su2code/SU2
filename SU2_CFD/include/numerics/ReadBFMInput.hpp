/*!
 * \file ReadBFMInput.hpp
 * \brief Declarations of template (empty) numerics classes, these give
 *        an idea of the methods that need to be defined to implement
 *        new schemes in SU2, in practice you should look for a similar
 *        scheme and try to re-use functionality (not by copy-paste).
 * \author F. Palacios, T. Economon
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

/*!
 * \class CConvectiveTemplate
 * \brief Class for setting up new method for spatial discretization of convective terms in flow equations.
 * \ingroup ConvDiscr
 * \author A. Lonkar
 */
class ReadBFMInput{
private:
  string version_input_file;
  

  unsigned short n_blade_rows{1};
  vector<unsigned long> n_axial_points{};
  vector<unsigned long> n_radial_points{};
  vector<unsigned long> n_tangential_points{};
  vector<unsigned short> n_blades{};
  vector<unsigned short> rotation_factor{};
  vector<string> variable_names{};

  vector<CVectorOfMatrix> *axial_coordinate = NULL;
  vector<CVectorOfMatrix> *radial_coordinate = NULL;
  vector<CVectorOfMatrix> *tangential_angle = NULL;

  vector<vector<CVectorOfMatrix>> *Geometric_Parameters = NULL;
  vector<string> translated_names;
  vector<pair<unsigned short, unsigned short>> name_translation;
    
  void AllocateMemory(){
    Geometric_Parameters = new vector<vector<CVectorOfMatrix>>;
    
    Geometric_Parameters->resize(n_blade_rows);
    
    for(unsigned short i_row = 0; i_row<n_blade_rows; ++i_row){
      
      Geometric_Parameters->at(i_row).resize(N_BFM_PARAMS);
      
      for(unsigned short i_param=0; i_param<N_BFM_PARAMS; ++i_param){
        Geometric_Parameters->at(i_row).at(i_param).resize(n_tangential_points.at(i_row), n_radial_points.at(i_row), n_axial_points.at(i_row));
      }
    }
  }

  void ReadInputFile(string file_name);

  string SkipToFlag(ifstream *file_stream, string flag);
  void SetNAxialPoints(unsigned long n_input, unsigned short i_row){n_axial_points.at(i_row) = n_input;}

  void SetNRadialPoints(unsigned long n_input, unsigned short i_row){n_radial_points.at(i_row) = n_input;}

  void SetNTangentialPoints(unsigned long n_input, unsigned short i_row){n_tangential_points.at(i_row) = n_input;}

  void SetNBladeRows(unsigned short n_input)
  {
    n_blade_rows = n_input;
  }
  
  void TranslateVariables();
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  ReadBFMInput(CConfig *config, string file_input_name);

  /*!
   * \brief Destructor of the class.
   */
  ~ReadBFMInput();

  /*!
  * \brief Get number of blade rows provided in BFM input file
  * \returns Blade row count.
  */
  unsigned short GetNBladeRows(){return n_blade_rows;}

  /*!
  * \brief Get number of data points in axial direction in BFM input file
  * \returns Number of data points in axial direction.
  */
  unsigned long GetNAxialPoints(unsigned short i_row){return n_axial_points.at(i_row);}

  unsigned long GetNRadialPoints(unsigned short i_row){return n_radial_points.at(i_row);}

  unsigned long GetNTangentialPoints(unsigned short i_row){return n_tangential_points.at(i_row);}

  su2double GetBFMParameter(unsigned short iRow, unsigned long iTang, unsigned long iRad, unsigned long iAx, unsigned short iVar){
    return Geometric_Parameters->at(iRow).at(iVar)(iTang, iRad, iAx);
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
  
  unsigned short GetBladeCount(unsigned short i_row){return n_blades.at(i_row);}
  unsigned short GetRotationFactor(unsigned short i_row){return rotation_factor.at(i_row);}

};
