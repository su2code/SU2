/*!
 * \file ReadBFMInput.cpp
 * \brief 
 * \author E.C. Bunschoten
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

#include "../../include/numerics/ReadBFMInput.hpp"

ReadBFMInput::ReadBFMInput(CConfig *config, string file_inputname) 
{

}

void ReadBFMInput::ReadInputFile(string input_file)
{
    ifstream file_stream;
    string line, word;
    int ixColon;
    bool eoHeader{false};
    bool eoData{false};
    bool eoRow{false};
    bool eoTang{false};
    bool eoRad{false};

    file_stream.open(input_file.c_str(), ifstream::in);

    if (!file_stream.is_open()) {
        SU2_MPI::Error(string("There is no look-up-table file file called ") + input_file,
                CURRENT_FUNCTION);
    }

    /* Read header */
    line = SkipToFlag(&file_stream, "<header>");

    while (getline(file_stream, line) && !eoHeader) {
        if (line.compare("[version]") == 0) {
            getline(file_stream, line);
            version_input_file = line;
        }
        if (line.compare("[number of blade rows]") == 0){
            getline(file_stream, line);
            SetNBladeRows(stoi(line));
        }
        if (line.compare("[number of tangential locations]") == 0){
            getline(file_stream, line);
            istringstream stream_tang_locations(line);
            while(stream_tang_locations){
                stream_tang_locations >> word;
                unsigned long tmp = stoul(word);
                n_tangential_points.push_back(tmp);
            }
            n_tangential_points.pop_back();
        }

        if (line.compare("[number of data entries in chordwise direction]") == 0){
            getline(file_stream, line);
            istringstream stream_tang_locations(line);
            while(stream_tang_locations){
                stream_tang_locations >> word;
                unsigned long tmp = stoul(word);
                n_axial_points.push_back(tmp);
            }
            n_axial_points.pop_back();
        }

        if (line.compare("[number of data entries in spanwise direction]") == 0){
            getline(file_stream, line);
            istringstream stream_tang_locations(line);
            while(stream_tang_locations){
                stream_tang_locations >> word;
                unsigned long tmp = stoul(word);
                n_radial_points.push_back(tmp);
            }
            n_radial_points.pop_back();
        }

        if (line.compare("[variable names]") == 0) {
        getline(file_stream, line);
        istringstream stream_names_var(line);
        while (stream_names_var) {
            stream_names_var >> word;
            ixColon = (int)word.find(":");
            variable_names.push_back(word.substr(ixColon + 1, word.size() - 1));
        }
        variable_names.pop_back();  // removes last redundant element
        }

        if (line.compare("</header>") == 0) eoHeader = true;
    }

    if((n_axial_points.size() != n_blade_rows) && (n_radial_points.size() != n_blade_rows) && (n_tangential_points.size() != n_blade_rows)){
        SU2_MPI::Error("Number of chordwise, spanwise or tangential data entries not provided for all blade rows", CURRENT_FUNCTION);
    }

    if(variable_names.size() != N_BFM_PARAMS + 2){
        SU2_MPI::Error("Number of variables in input file are unequal to the number required for Blade geometry definition", CURRENT_FUNCTION);
    }

    AllocateMemory();

    line = SkipToFlag(&file_stream, "<data>");
    

}

string ReadBFMInput::SkipToFlag(ifstream *file_stream, string flag) {
  string line;
  getline(*file_stream, line);

  while (line.compare(flag) != 0 && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof())
    SU2_MPI::Error("Flag not found in file", CURRENT_FUNCTION);

  return line;
}

ReadBFMInput::~ReadBFMInput(void) 
{
    if(axial_coordinate != NULL){
        delete axial_coordinate;
        delete radial_coordinate;
        delete tangential_angle;
        delete Geometric_Parameters;
    }
}

