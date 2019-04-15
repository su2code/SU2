/*!
 * \file output_driver.cpp
 * \brief Main subroutines for multizone output
 * \author R. Sanchez, T. Albring
 * \version 6.1.0 "Falcon"
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

#include "../include/output/output_driver.hpp"

CDriverOutput::CDriverOutput(CConfig* driver_config, CConfig** config) : COutput(driver_config) {

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
    }
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  
  

 
  MultiZoneHeaderString = "Multizone Summary";
  
  HistoryFilename = "multizone_history";
      
}

CDriverOutput::~CDriverOutput() {


}


void CDriverOutput::LoadMultizoneHistoryData(COutput **output, CConfig **config) {

  unsigned short iZone, iField, nField;
  string name, header;
  su2double val, avgzone;

  if (config[ZONE_0]->GetTime_Domain()){
    SetHistoryOutputValue("TIME_ITER", curr_TimeIter);
  }
  SetHistoryOutputValue("OUTER_ITER", curr_OuterIter);
  
  
  for (iZone = 0; iZone < nZone; iZone++){
    
    map<string, HistoryOutputField> ZoneHistoryFields = output[iZone]->GetHistoryFields();
    vector<string>                  ZoneHistoryNames  = output[iZone]->GetHistoryOutput_List();
    
    nField = ZoneHistoryFields.size();
    
    
    /*-- For all the variables per solver --*/
    for (iField = 0; iField < nField; iField++){
      
      if (ZoneHistoryNames[iField] != "TIME_ITER" && ZoneHistoryNames[iField] != "OUTER_ITER"){
        
        name   = ZoneHistoryNames[iField]+ "[" + PrintingToolbox::to_string(iZone) + "]";
        
        SetHistoryOutputValue(name, ZoneHistoryFields[ZoneHistoryNames[iField]].Value);

      }
    }
  }

//  for (iZone = 0; iZone < nZone; iZone++){

//    /*--- Initialize values per zone ---*/
//    nVar_Zone = 0;
//    avgzone = 0.0;
    
//    InnerIter = output[iZone]->GetHistoryFieldValue("INNER_ITER");
    
//    name = string("INNER_ITER") + "_" + PrintingToolbox::to_string(iZone);
    
//    SetHistoryOutputValue(name, InnerIter);
    
//    vector<HistoryOutputField> BGSGroup = output[iZone]->GetHistoryGroup("BGS_RES");   
    
        
//    nField = BGSGroup.size();
    
//    for (iField = 0; iField < nField; iField++){
      
//      name = BGSGroup[iField].FieldName + "_" + PrintingToolbox::to_string(iZone);
      
//      val = BGSGroup[iField].Value;
      
//      SetHistoryOutputValue(name, val);
      
//      /*--- Add values to averages ---*/
//      avgzone += val;
//      nVar_Zone++;
      
//    }

//    /*--- Get an unique name for the averaged history data per zone ---*/
//    name = "ZONE" + PrintingToolbox::to_string(iZone);

//    /*--- Compute the average and set the value for the zone iZone---*/
//    avgzone = avgzone / nVar_Zone;
//    SetHistoryOutputValue(name, avgzone);
    
//    vector<HistoryOutputField> RMSGroup = output[iZone]->GetHistoryGroup("RMS_RES");        
        
//    nField = RMSGroup.size();
    
//    for (iField = 0; iField < nField; iField++){
      
//      name = RMSGroup[iField].FieldName + "_" + PrintingToolbox::to_string(iZone);
      
//      val = RMSGroup[iField].Value;
      
//      SetHistoryOutputValue(name, val);
      
//    }

//  }


}

void CDriverOutput::SetMultizoneHistoryOutputFields(COutput **output, CConfig **config) {
  
  unsigned short iZone, iField, nField;
  string name, header, group;
  
  if (nRequestedScreenFields == 0){
    if (config[ZONE_0]->GetTime_Domain()) RequestedScreenFields.push_back("TIME_ITER");
    RequestedScreenFields.push_back("OUTER_ITER");
    
    for (iZone = 0; iZone < nZone; iZone++){
      
      map<string, HistoryOutputField> ZoneHistoryFields = output[iZone]->GetHistoryFields();
      vector<string>                  ZoneHistoryNames  = output[iZone]->GetHistoryOutput_List();
      nField = ZoneHistoryFields.size();
      
      for (iField = 0; iField < nField; iField++){
        
        if (ZoneHistoryFields[ZoneHistoryNames[iField]].OutputGroup == bgs_res_name){
          
          name   = ZoneHistoryNames[iField]+ "[" + PrintingToolbox::to_string(iZone) + "]";
          
          RequestedScreenFields.push_back(name);
          
          break;
          
        }
      }
    }
    
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  
  if (config[ZONE_0]->GetTime_Domain()){
    AddHistoryOutput("TIME_ITER", "Time_Iter", FORMAT_INTEGER,  "ITER");
  }
  AddHistoryOutput("OUTER_ITER", "Outer_Iter", FORMAT_INTEGER,  "ITER");
  
  
  /*--- Set the fields ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    
    map<string, HistoryOutputField> ZoneHistoryFields = output[iZone]->GetHistoryFields();
    vector<string>                  ZoneHistoryNames  = output[iZone]->GetHistoryOutput_List();
    
    nField = ZoneHistoryFields.size();
    
    
    /*-- For all the variables per solver --*/
    for (iField = 0; iField < nField; iField++){
      
      if (ZoneHistoryNames[iField] != "TIME_ITER" && ZoneHistoryNames[iField] != "OUTER_ITER"){
        
        name   = ZoneHistoryNames[iField]+ "[" + PrintingToolbox::to_string(iZone) + "]";
        header = ZoneHistoryFields[ZoneHistoryNames[iField]].FieldName + "[" + PrintingToolbox::to_string(iZone) + "]";
        group  = ZoneHistoryFields[ZoneHistoryNames[iField]].OutputGroup + "[" + PrintingToolbox::to_string(iZone) + "]";
        
        AddHistoryOutput(name, header, ZoneHistoryFields[ZoneHistoryNames[iField]].ScreenFormat, group, ZoneHistoryFields[ZoneHistoryNames[iField]].FieldType );
      }
    }
  }

//  /*--- Set the fields ---*/
//  for (iZone = 0; iZone < nZone; iZone++){
    
//    InnerIter = output[iZone]->GetHistoryFieldValue("INNER_ITER");
    
//    name   = string("INNER_ITER") + "_" + PrintingToolbox::to_string(iZone);
//    header = string("Inner_Iter") + "[" + PrintingToolbox::to_string(iZone) + "]";
    
//    AddHistoryOutput(name, header, FORMAT_INTEGER, "ITER");
    
////    RequestedScreenFields.push_back(name);
    
//    vector<HistoryOutputField> BGSGroup = output[iZone]->GetHistoryGroup("BGS_RES");        
    
//    nField = BGSGroup.size();
 
//    /*-- For all the variables per solver --*/
//    for (iField = 0; iField < nField; iField++){
      
//      /*--- Set an unique name for the history data ---*/
//      name   = BGSGroup[iField].FieldName + "_" + PrintingToolbox::to_string(iZone);
      
//      header = BGSGroup[iField].FieldName + "[" + PrintingToolbox::to_string(iZone) + "]";

//      AddHistoryOutput(name, header, FORMAT_FIXED,  "BGS_RES", TYPE_RESIDUAL);
      
//      /*--- Request the variable residual for history output ---*/
//      RequestedScreenFields.push_back(name);
//    }

//    vector<HistoryOutputField> RMSGroup = output[iZone]->GetHistoryGroup("RMS_RES");
//    nField = RMSGroup.size();
    
//    /*-- For all the variables per solver --*/
//    for (iField = 0; iField < nField; iField++){
      
//      /*--- Set an unique name for the history data ---*/
//      name   = RMSGroup[iField].FieldName + "_" + PrintingToolbox::to_string(iZone);
      
//      header = RMSGroup[iField].FieldName + "[" + PrintingToolbox::to_string(iZone) + "]";

//      AddHistoryOutput(name, header, FORMAT_FIXED,  "RMS_RES", TYPE_RESIDUAL);
      
//      /*--- Request the variable residual for history output ---*/
////      RequestedScreenFields.push_back(name);
//    }
    

    
 
//    /*--- Set an unique name for the averaged history data ---*/
//    name = "ZONE" + PrintingToolbox::to_string(iZone);
//    /*--- Set an unique name for the history headers of the averaged data ---*/
//    header = "avgres[" + PrintingToolbox::to_string(iZone) + "]";

//    AddHistoryOutput(name, header, FORMAT_FIXED,  "ZONE_AVGRES", TYPE_RESIDUAL);
//  }
  
//  RequestedHistoryFields.push_back("BGS_RES");
//  RequestedHistoryFields.push_back("RMS_RES");
//  RequestedHistoryFields.push_back("ZONE_AVGRES");
  
//  /*--- Retrieve number of requested fields ---*/
//  nRequestedScreenFields = RequestedScreenFields.size();
//  nRequestedHistoryFields = RequestedHistoryFields.size();


}

inline bool CDriverOutput::WriteScreen_Header(CConfig *config) {

  bool write_header = true;

  /*--- If the outer iteration is zero ---*/
  write_header = (write_header && (curr_OuterIter == 0)) || write_zone;

  return write_header;

}
inline bool CDriverOutput::WriteScreen_Output(CConfig *config) {

  bool write_output = true;

  return write_output;
}

