/*!
 * \file output_driver.cpp
 * \brief Main subroutines for multizone output
 * \author R. Sanchez, T. Albring
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/output/CDriverOutput.hpp"

CDriverOutput::CDriverOutput(CConfig* driver_config, CConfig** config, unsigned short nDim) : COutput(driver_config, nDim) {

  unsigned short iZone = 0;
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  nZone = driver_config->GetnZone();

  field_width = 12;
  
  bgs_res_name = "BGS_RES";
  
  write_zone = false;

  /*--- If the zone output is disabled for every zone ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    write_zone = config[iZone]->GetWrt_ZoneConv();
  }

  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("ITER");
    for (iZone = 0; iZone < nZone; iZone++){
      RequestedHistoryFields.push_back(bgs_res_name + "[" + PrintingToolbox::to_string(iZone) + "]");      
      RequestedHistoryFields.push_back("AVG_RES[" + PrintingToolbox::to_string(iZone) + "]");
    }
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  
  if (nRequestedScreenFields == 0){
    if (config[ZONE_0]->GetTime_Domain()) RequestedScreenFields.push_back("TIME_ITER");    
    RequestedScreenFields.push_back("OUTER_ITER");
    for (iZone = 0; iZone < nZone; iZone++){
      RequestedScreenFields.push_back("AVG_" + bgs_res_name + "[" + PrintingToolbox::to_string(iZone) + "]"); 
    }
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  
  MultiZoneHeaderString = "Multizone Summary";
  
  HistoryFilename = "multizone_history";

  /*--- Set the default convergence field --- */

  if (Conv_Field.size() == 0 ) Conv_Field = "AVG_BGS_RES[0]";

}

CDriverOutput::~CDriverOutput() {


}


void CDriverOutput::LoadMultizoneHistoryData(COutput **output, CConfig **config) {

  unsigned short iZone, iField, nField;
  string name, header;

  if (config[ZONE_0]->GetTime_Domain()){
    SetHistoryOutputValue("TIME_ITER", curr_TimeIter);
  }
  SetHistoryOutputValue("OUTER_ITER", curr_OuterIter);
  
  
  for (iZone = 0; iZone < nZone; iZone++){
    
    map<string, HistoryOutputField> ZoneHistoryFields = output[iZone]->GetHistoryFields();
    vector<string>                  ZoneHistoryNames  = output[iZone]->GetHistoryOutput_List();
    
    nField = ZoneHistoryNames.size();
    
    
    /*-- For all the variables per solver --*/
    for (iField = 0; iField < nField; iField++){
      
      if (ZoneHistoryNames[iField] != "TIME_ITER" && ZoneHistoryNames[iField] != "OUTER_ITER"){
        
        name   = ZoneHistoryNames[iField]+ "[" + PrintingToolbox::to_string(iZone) + "]";
        
        SetHistoryOutputValue(name, ZoneHistoryFields[ZoneHistoryNames[iField]].Value);

      }
    }
  }
}

void CDriverOutput::SetMultizoneHistoryOutputFields(COutput **output, CConfig **config) {
  
  unsigned short iZone, iField, nField;
  string name, header, group;

  if (config[ZONE_0]->GetTime_Domain()){
    AddHistoryOutput("TIME_ITER", "Time_Iter", FORMAT_INTEGER,  "ITER", "Time iteration index");
  }
  AddHistoryOutput("OUTER_ITER", "Outer_Iter", FORMAT_INTEGER,  "ITER", "Outer iteration index");
  
  
  /*--- Set the fields ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    
    map<string, HistoryOutputField> ZoneHistoryFields = output[iZone]->GetHistoryFields();
    vector<string>                  ZoneHistoryNames  = output[iZone]->GetHistoryOutput_List();
    
    nField = ZoneHistoryNames.size();
    
    
    /*-- For all the variables per solver --*/
    for (iField = 0; iField < nField; iField++){
      
      if (ZoneHistoryNames[iField] != "TIME_ITER" && ZoneHistoryNames[iField] != "OUTER_ITER"){
        
        name   = ZoneHistoryNames[iField]+ "[" + PrintingToolbox::to_string(iZone) + "]";
        header = ZoneHistoryFields[ZoneHistoryNames[iField]].FieldName + "[" + PrintingToolbox::to_string(iZone) + "]";
        group  = ZoneHistoryFields[ZoneHistoryNames[iField]].OutputGroup + "[" + PrintingToolbox::to_string(iZone) + "]";
        
        AddHistoryOutput(name, header, ZoneHistoryFields[ZoneHistoryNames[iField]].ScreenFormat, group, "", ZoneHistoryFields[ZoneHistoryNames[iField]].FieldType );
      }
    }
  }
}

inline bool CDriverOutput::WriteScreen_Header(CConfig *config) {

  bool write_header = true;

  /*--- If the outer iteration is zero ---*/
  write_header = (write_header && (curr_OuterIter == 0)) || write_zone;

  return write_header;

}
inline bool CDriverOutput::WriteScreen_Output(CConfig *config) {
  
  su2double* ScreenWrt_Freq = config->GetScreen_Wrt_Freq();

  /*--- Check if screen output should be written --- */
  
  if (!PrintOutput(curr_TimeIter, SU2_TYPE::Int(ScreenWrt_Freq[0]))&& 
      !(curr_TimeIter == config->GetnTime_Iter() - 1)){
    
    return false;
    
  }
  
  if (Convergence) {return true;}
  
  if (!PrintOutput(curr_OuterIter, SU2_TYPE::Int(ScreenWrt_Freq[1])) && 
      !(curr_OuterIter == config->GetnOuter_Iter() - 1)){
    
    return false;
    
  }
  
 
  return true;
}

inline bool CDriverOutput::WriteHistoryFile_Output(CConfig *config){
  
  su2double* HistoryWrt_Freq = config->GetHistory_Wrt_Freq();
    
  /*--- Check if screen output should be written --- */
  
  if (!PrintOutput(curr_TimeIter, SU2_TYPE::Int(HistoryWrt_Freq[0]))&& 
      !(curr_TimeIter == config->GetnTime_Iter() - 1)){
    
    return false;
    
  }
  
  if (Convergence) {return true;}
  
  if (!PrintOutput(curr_OuterIter, SU2_TYPE::Int(HistoryWrt_Freq[1])) && 
      !(curr_OuterIter == config->GetnOuter_Iter() - 1)){
    
    return false;
    
  }
 
  return true;
}
