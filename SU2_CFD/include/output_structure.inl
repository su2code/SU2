/*!
 * \file output_structure.inl
 * \brief In-Line subroutines of the <i>output_structure.hpp</i> file.
 * \author J. Smith
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

#pragma once
#include "output_structure.hpp"

inline su2double COutput::GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan) { return EntropyGen[iMarkerTP][iSpan]; }

inline su2double COutput::GetFlowAngleOut(unsigned short iMarkerTP, unsigned short iSpan) { return FlowAngleOut[iMarkerTP][iSpan]*180.0/PI_NUMBER; }

inline su2double COutput::GetMassFlowIn(unsigned short iMarkerTP, unsigned short iSpan) { return MassFlowIn[iMarkerTP][iSpan]; }

inline bool COutput::PrintOutput(unsigned long iIter, unsigned long iFreq) { return (iIter % iFreq == 0); }

inline void COutput::SetHistoryOutputFields(CConfig *config){}

inline void COutput::SetConvHistory_Header(CConfig *config, unsigned short val_iZone, unsigned short val_iInst) { }
 
inline bool COutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { return true; }

inline bool COutput::WriteScreen_Header(CConfig *config) { return true; }

inline bool COutput::WriteScreen_Output(CConfig *config, bool write_dualtime) { return true; }

inline void COutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) { }

inline void COutput::PrintScreenFixed(stringstream& stream, su2double val) {
  stream.precision(6); stream.setf(ios::fixed, ios::floatfield); stream.width(field_width);
  stream << std::right << val;
  stream.unsetf(ios::fixed);
}

inline void COutput::PrintScreenScientific(stringstream& stream, su2double val) {
  stream.precision(4); stream.setf(ios::scientific, ios::floatfield); stream.width(field_width);
  stream << std::right << val;
  stream.unsetf(ios::scientific);  
}

inline void COutput::PrintScreenInteger(stringstream& stream, unsigned long val) {
  stream.width(field_width);
  stream << std::right << val;
}

inline void COutput::PrintScreenHeaderString(stringstream& stream, string header) {
  if (header.size() > field_width-1) header.resize(field_width-1);
  stream << std::right << std::setw(field_width) << header; 
}

inline void COutput::PrintHistorySep(stringstream& stream){
  stream << HistorySep;
}

inline void COutput::AddHistoryOutput(string name, string field_name, unsigned short format, string groupname , unsigned short field_type){
  HistoryOutput_Map[name] = HistoryOutputField(field_name, format, groupname, field_type);
  HistoryOutput_List.push_back(name);
}

inline void COutput::AddHistoryOutputPerSurface(string name, string field_name, unsigned short format, string groupname, vector<string> marker_names, unsigned short field_type){
  if (marker_names.size() != 0){
    HistoryOutputPerSurface_List.push_back(name);
    for (unsigned short i = 0; i < marker_names.size(); i++){
      HistoryOutputPerSurface_Map[name].push_back(HistoryOutputField(field_name+"("+marker_names[i]+")", format, groupname, field_type));
    }
  }
}

inline void COutput::SetHistoryOutputValue(string name, su2double value){
  if (HistoryOutput_Map.count(name) > 0){
    HistoryOutput_Map[name].Value = value;
  } else {
    SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
  }
}

inline void COutput::SetHistoryOutputPerSurfaceValue(string name, su2double value, unsigned short iMarker){
  if (HistoryOutputPerSurface_Map.count(name) > 0){
    HistoryOutputPerSurface_Map[name][iMarker].Value = value;
  } else {
    SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
  }
}

inline void COutput::AddVolumeOutput(string name, string field_name, string groupname){
  VolumeOutput_Map[name] = VolumeOutputField(field_name, -1, groupname);
  VolumeOutput_List.push_back(name);
}

inline void COutput::SetVolumeOutputValue(string name, unsigned long iPoint, su2double value){
  if (VolumeOutput_Map.count(name) > 0){
    if (VolumeOutput_Map[name].Offset != -1){
      Local_Data[iPoint][VolumeOutput_Map[name].Offset] = value;
    }
  } else {
    SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);    
  }
}
inline void COutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint) {}

inline void COutput::SetVolumeOutputFields(CConfig *config) {}

inline void COutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex) {}

inline bool COutput::SetInit_Residuals(CConfig *config) {return false;}

inline bool COutput::SetUpdate_Averages(CConfig *config, bool dualtime) {return false;}