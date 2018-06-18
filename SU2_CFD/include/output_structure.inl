/*!
 * \file output_structure.inl
 * \brief In-Line subroutines of the <i>output_structure.hpp</i> file.
 * \author J. Smith
 * \version 6.0.1 "Falcon"
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

inline su2double COutput::GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan) { return EntropyGen[iMarkerTP][iSpan]; }

inline su2double COutput::GetFlowAngleOut(unsigned short iMarkerTP, unsigned short iSpan) { return FlowAngleOut[iMarkerTP][iSpan]*180.0/PI_NUMBER; }

inline su2double COutput::GetMassFlowIn(unsigned short iMarkerTP, unsigned short iSpan) { return MassFlowIn[iMarkerTP][iSpan]; }

inline bool COutput::PrintOutput(unsigned long iIter, unsigned long iFreq) { return (iIter % iFreq == 0); }

inline void COutput::SetOutputFields(CConfig *config){}

inline void COutput::SetConvHistory_Header(CConfig *config, unsigned short val_iZone, unsigned short val_iInst) { }
 
inline bool COutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { }

inline bool COutput::WriteScreen_Header(CConfig *config) { }

inline bool COutput::WriteScreen_Output(CConfig *config, bool write_dualtime) { }

inline void COutput::LoadOutput_Data(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
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

inline void COutput::AddHistoryValue(su2double val) { 
  HistoryValues.push_back(val);
}

inline void COutput::AddHistoryHeaderString(string header) {
  HistoryHeader.push_back(header); 
}

inline void COutput::PrintHistorySep(stringstream& stream){
  stream << HistorySep;
}

inline void COutput::AddOutputField(string name, string field_name, unsigned short format, string groupname ){
  Output_Fields[name] = OutputField(field_name, format, groupname);
}

inline void COutput::SetOutputFieldValue(string name, su2double value){
  Output_Fields[name].Value = value;
}