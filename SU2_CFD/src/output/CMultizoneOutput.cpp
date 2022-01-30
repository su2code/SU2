/*!
 * \file CMultizoneOutput.cpp
 * \brief Main subroutines for multizone output
 * \author R. Sanchez, T. Albring
 * \version 7.3.0 "Blackbird"
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

#include "../../include/output/CMultizoneOutput.hpp"

CMultizoneOutput::CMultizoneOutput(const CConfig* driver_config, const CConfig* const* config, unsigned short nDim) :
  COutput(driver_config, nDim, false) {

  unsigned short iZone = 0;

  nZone = driver_config->GetnZone();

  bgs_res_name = "BGS_RES";

  write_zone = false;

  /*--- If the zone output is disabled for every zone ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    write_zone = config[iZone]->GetWrt_ZoneConv();
  }

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    for (iZone = 0; iZone < nZone; iZone++){
      requestedHistoryFields.emplace_back(bgs_res_name + "[" + PrintingToolbox::to_string(iZone) + "]");
      requestedHistoryFields.emplace_back("AVG_RES[" + PrintingToolbox::to_string(iZone) + "]");
    }
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0){
    if (config[ZONE_0]->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    requestedScreenFields.emplace_back("OUTER_ITER");
    for (iZone = 0; iZone < nZone; iZone++){
      requestedScreenFields.emplace_back("AVG_" + bgs_res_name + "[" + PrintingToolbox::to_string(iZone) + "]");
    }
    nRequestedScreenFields = requestedScreenFields.size();
  }

  multiZoneHeaderString = "Multizone Summary";

  historyFilename = driver_config->GetCaseName();

  /*--- Add the correct file extension depending on the file format ---*/

  string hist_ext = ".csv";
  if (driver_config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_TECPLOT) hist_ext = ".dat";

  historyFilename += hist_ext;

  /*--- Set the default convergence field --- */

  if (convFields.empty()) convFields.emplace_back("AVG_BGS_RES[0]");

}

void CMultizoneOutput::LoadMultizoneHistoryData(const COutput* const* output) {

  string nameMultizone, zoneIndex;
  su2double comboValue = 0;

  for (unsigned short iZone = 0; iZone < nZone; iZone++) {
    zoneIndex = "[" + PrintingToolbox::to_string(iZone) + "]";

    /*--- For all the variables per solver ---*/
    for (const auto& item : output[iZone]->GetHistoryFields()) {
      const auto& name = item.first;
      if (name != "TIME_ITER" && name != "OUTER_ITER") {
        nameMultizone = name + zoneIndex;
        SetHistoryOutputValue(nameMultizone, item.second.value);
      }
      if (name == "COMBO") {
        comboValue += item.second.value;
      }
    }
  }
  SetHistoryOutputValue("COMBO", comboValue);
}

void CMultizoneOutput::SetMultizoneHistoryOutputFields(const COutput* const* output, const CConfig* const* config) {

  string name, header, group, zoneIndex;

  /*--- Set the fields ---*/
  for (unsigned short iZone = 0; iZone < nZone; iZone++) {

    const auto& ZoneHistoryFields = output[iZone]->GetHistoryFields();
    zoneIndex = "[" + PrintingToolbox::to_string(iZone) + "]";

    /*--- For all the variables per solver ---*/
    for (const auto& nameSinglezone : output[iZone]->GetHistoryOutput_List()) {

      if (nameSinglezone != "TIME_ITER" && nameSinglezone != "OUTER_ITER") {

        const auto& field = ZoneHistoryFields.at(nameSinglezone);

        name   = nameSinglezone + zoneIndex;
        header = field.fieldName + zoneIndex;
        group  = field.outputGroup + zoneIndex;

        AddHistoryOutput(name, header, field.screenFormat, group, "", field.fieldType );
      }
    }
    /*--- Also add the PerSurface outputs ---*/
    const auto& ZoneHistoryPerSurfaceFields = output[iZone]->GetHistoryPerSurfaceFields();

    vector<string> Marker_Analyze;
    for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config[iZone]->GetnMarker_Analyze(); iMarker_Analyze++){
      Marker_Analyze.push_back(config[iZone]->GetMarker_Analyze_TagBound(iMarker_Analyze));
    }
    //Marker_Analyze.push_back("test");

    for (const auto& nameSinglezone : output[iZone]->GetHistoryOutputPerSurface_List()) {

      if (nameSinglezone != "TIME_ITER" && nameSinglezone != "OUTER_ITER") {

        const auto& field = ZoneHistoryPerSurfaceFields.at(nameSinglezone);

        name   = nameSinglezone + zoneIndex;

        // field[0].fieldName contains names like Avg_Massflow(fluid_outlet) but looking at CFlowOutput.cpp it is also:
        // AddHistoryOutputPerSurface("SURFACE_MASSFLOW",         "Avg_Massflow", ... without the markerTag already attached
        // So first the fieldName for one of the markers (here the 0th) is taken and the markerTag in brackets is removed
        string baseheader = field[0].fieldName;
        std::string::size_type pos = baseheader.find('(');
        if (pos != std::string::npos)
            baseheader = baseheader.substr(0, pos);
        else
            SU2_MPI::Error("Cannot proccess PerSurface *_SURF history output.", CURRENT_FUNCTION);

        header = baseheader + zoneIndex; // field[i] where i is the marker
        group  = field[0].outputGroup + zoneIndex;

        AddHistoryOutputPerSurface(name, header, field[0].screenFormat, group, Marker_Analyze, field[0].fieldType );
        //AddHistoryOutputPerSurface("SURFACE_MASSFLOW", "Avg_Massflow", ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
        //AddHistoryOutput           "SURFACE_MASSFLOW", "Avg_Massflow", ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF",      "desc",         HistoryFieldType::COEFFICIENT);
      }
    }
  }
  AddHistoryOutput("COMBO", "ComboObj", ScreenOutputFormat::SCIENTIFIC, "COMBO", "Combined obj. function value.", HistoryFieldType::COEFFICIENT);
}

bool CMultizoneOutput::WriteScreen_Header(const CConfig *config) {

  /*--- Print header if the outer iteration is zero or zonal convergence is printed ---*/

  return curOuterIter == 0 || write_zone;

  /*--- Always print header if it is forced ---*/

  if (headerNeeded){
    headerNeeded = false;
    return true;
  }

  return false;
}

bool CMultizoneOutput::WriteScreen_Output(const CConfig *config) {

  unsigned long ScreenWrt_Freq_Outer = config->GetScreen_Wrt_Freq(1);
  unsigned long ScreenWrt_Freq_Time  = config->GetScreen_Wrt_Freq(0);

  /*--- Check if screen output should be written --- */

  if (!PrintOutput(curTimeIter, ScreenWrt_Freq_Time)&&
      !(curTimeIter == config->GetnTime_Iter() - 1)){

    return false;

  }

  if (convergence) {return true;}

  if (!PrintOutput(curOuterIter, ScreenWrt_Freq_Outer) &&
      !(curOuterIter == config->GetnOuter_Iter() - 1)){

    return false;

  }

  return true;
}

bool CMultizoneOutput::WriteHistoryFile_Output(const CConfig *config){

  unsigned long HistoryWrt_Freq_Outer = config->GetHistory_Wrt_Freq(1);
  unsigned long HistoryWrt_Freq_Time  = config->GetHistory_Wrt_Freq(0);

  /*--- Check if screen output should be written --- */

  if (!PrintOutput(curTimeIter, HistoryWrt_Freq_Time)&&
      !(curTimeIter == config->GetnTime_Iter() - 1)){

    return false;

  }

  if (convergence) {return true;}

  if (!PrintOutput(curOuterIter, HistoryWrt_Freq_Outer) &&
      !(curOuterIter == config->GetnOuter_Iter() - 1)){

    return false;

  }

  return true;
}
