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
    //vector<string> translated_names;
    translated_names.resize(N_BFM_PARAMS);
    translated_names.at(I_AXIAL_COORDINATE) = "axial_coordinate";
    translated_names.at(I_RADIAL_COORDINATE) = "radial_coordinate";
    translated_names.at(I_BLOCKAGE_FACTOR) = "blockage_factor";
    translated_names.at(I_CAMBER_NORMAL_AXIAL) = "n_ax";
    translated_names.at(I_CAMBER_NORMAL_TANGENTIAL) = "n_tang";
    translated_names.at(I_CAMBER_NORMAL_RADIAL) = "n_rad";
    translated_names.at(I_LEADING_EDGE_AXIAL) = "x_LE";
    translated_names.at(I_ROTATION_FACTOR) = "rotation_factor";
    translated_names.at(I_BLADE_COUNT) = "blade_count";

    ReadInputFile(file_inputname);
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
        if (line.compare("[row blade count]") == 0){
            getline(file_stream, line);
            istringstream stream_tang_locations(line);
            while(stream_tang_locations){
                stream_tang_locations >> word;
                unsigned long tmp = stoul(word);
                n_blades.push_back(tmp);
            }
            n_blades.pop_back();
        }
        if (line.compare("[rotation factor]") == 0){
            getline(file_stream, line);
            istringstream stream_tang_locations(line);
            while(stream_tang_locations){
                stream_tang_locations >> word;
                unsigned long tmp = stoi(word);
                rotation_factor.push_back(tmp);
            }
            rotation_factor.pop_back();
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

    TranslateVariables();

    cout << "Number of blade rows: " << GetNBladeRows() << endl;
    cout << "Number of spanwise sections: "<< n_axial_points.at(0) << endl;
    cout << "Variables: ";
    for(size_t i=0; i<variable_names.size(); ++i){
        cout << variable_names.at(i) << " ";
    }
    cout << endl;

    //TranslateVariables();
    
    if((n_axial_points.size() != n_blade_rows) && (n_radial_points.size() != n_blade_rows) && (n_tangential_points.size() != n_blade_rows)){
        SU2_MPI::Error("Number of chordwise, spanwise or tangential data entries not provided for all blade rows", CURRENT_FUNCTION);
    }

    // if(variable_names.size() != N_BFM_PARAMS-2){
    //     SU2_MPI::Error("Number of variables in input file are unequal to the number required for Blade geometry definition", CURRENT_FUNCTION);
    // }

    AllocateMemory();
    
    line = SkipToFlag(&file_stream, "<data>");
    unsigned short rowCounter{0};
    unsigned long tangCounter{0};
    unsigned long radCounter{0};
    unsigned long pointCounter{0};
    su2double temp;
    while (getline(file_stream, line) && !eoData) {
        if (line.compare("</data>") == 0){
            eoData = true;
        }
        if(!eoData){
            if(line.compare("<blade row>") == 0){
                eoRow = false;
                tangCounter = 0;
                getline(file_stream, line);
            }
            if(line.compare("<tang section>") == 0){
                eoTang = false;
                radCounter = 0;
                getline(file_stream, line);
            }
            if(line.compare("<radial section>") == 0){
                eoRad = false;
                pointCounter = 0;
                getline(file_stream, line);
            }
            if(line.compare("</blade row>") == 0){
                eoRow = true;
                rowCounter++;
                //getline(file_stream, line);
            }
            if(line.compare("</tang section>") == 0){
                eoTang = true;
                tangCounter++;
                //getline(file_stream, line);
            }
            if(line.compare("</radial section>") == 0){
                eoRad = true;
                radCounter++;
                //getline(file_stream, line);
            }

            if(!eoRow && !eoTang && !eoRad){
                istringstream streamDataLine(line);
                for (unsigned long iVar = 0; iVar < variable_names.size(); iVar++) {
                    streamDataLine >> word;
                    temp = stod(word);
                    for(unsigned short iName=0; iName < translated_names.size(); ++iName){
                        if((translated_names.at(iName).compare(variable_names.at(iVar)) == 0)
                         && (translated_names.at(I_ROTATION_FACTOR).compare(variable_names.at(iVar)) != 0)
                         && (translated_names.at(I_BLADE_COUNT).compare(variable_names.at(iVar)) != 0)){
                            Geometric_Parameters->at(rowCounter).at(iName)(tangCounter, radCounter, pointCounter) = temp;
                        }
                    }
                }
                
                pointCounter ++;
            }
        }
    }

}

void ReadBFMInput::TranslateVariables(){
    name_translation.resize(variable_names.size());
    for(unsigned short iVar=0; iVar<variable_names.size(); ++iVar){
        if(variable_names.at(iVar) == "axial_coordinate"){
            name_translation.at(iVar) = make_pair(iVar, I_AXIAL_COORDINATE);
        }
        if(variable_names.at(iVar) == "radial coordinte"){
            name_translation.at(iVar) = make_pair(iVar, I_RADIAL_COORDINATE);
        }
        if(variable_names.at(iVar) == "blockage_factor"){
            name_translation.at(iVar) = make_pair(iVar, I_BLOCKAGE_FACTOR);
        }
        if(variable_names.at(iVar) == "n_ax"){
            name_translation.at(iVar) = make_pair(iVar, I_CAMBER_NORMAL_AXIAL);
        }
        if(variable_names.at(iVar) == "n_tang"){
            name_translation.at(iVar) = make_pair(iVar, I_CAMBER_NORMAL_TANGENTIAL);
        }
        if(variable_names.at(iVar) == "n_rad"){
            name_translation.at(iVar) = make_pair(iVar, I_CAMBER_NORMAL_RADIAL);
        }
        if(variable_names.at(iVar) == "x_LE"){
            name_translation.at(iVar) = make_pair(iVar, I_LEADING_EDGE_AXIAL);
        }
        if(variable_names.at(iVar) == "rotation_factor"){
            name_translation.at(iVar) = make_pair(iVar, I_ROTATION_FACTOR);
        }
        if(variable_names.at(iVar) == "blade_count"){
            name_translation.at(iVar) = make_pair(iVar, I_BLADE_COUNT);
        }
        
        
    }
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
    if(Geometric_Parameters != NULL){
        delete Geometric_Parameters;
    }
    if(axial_coordinate != NULL){
        delete axial_coordinate;
        delete radial_coordinate;
        delete tangential_angle;
        delete Geometric_Parameters;
    }
}

