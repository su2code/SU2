/*!
 * \file datatype_structure.inl
 * \brief In-Line subroutines of the <i>datatype_structure.hpp</i> file.
 * \author T. Albring
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

inline void CTypeWrapper::SetPrimary(su2double& data, const double &val){
  data = val;
}


inline double CTypeWrapper::GetPrimary(su2double& data){
  return data;
}

inline void CTypeWrapper::SetSecondary(su2double& data, const double &val){

}


inline double CTypeWrapper::GetSecondary(su2double& data){
  return 0.0;
}

#ifdef COMPLEX_TYPE
inline void CComplexTypeWrapper::SetPrimary(su2double& data, const double &val){
  data = su2double(val, data.imag());
}


inline double CComplexTypeWrapper::GetPrimary(su2double& data){
  return data.real();
}

inline void CComplexTypeWrapper::SetSecondary(su2double& data, const double &val){
  data = su2double(data.real(), val);
}

inline double CComplexTypeWrapper::GetSecondary(su2double& data){
  return data.imag();
}
#endif
