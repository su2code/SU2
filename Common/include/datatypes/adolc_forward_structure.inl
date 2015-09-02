/*!
 * \file adolc_forward_structure.inl
 * \brief Inline subroutines for <i>datatype_structure.hpp<i>.
 * \author T. Albring
 * \version 4.0.1 "Cardinal"
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

namespace SU2_TYPE{
  inline void SetPrimary(su2double& data, const double &val){data.setValue(val);}

  inline double GetPrimary(const su2double& data){return data.getValue();}

  inline void SetSecondary(su2double& data, const double &val){data.setADValue(&val);}

  inline double GetSecondary(const su2double& data){return *data.getADValue();}

  inline double GetDerivative(const su2double& data){return *data.getADValue();}

  inline void SetDerivative(su2double& data, const double &val){data.setADValue(&val);}
}

namespace adtl{
  inline std::ostream& operator<<(std::ostream& out, const su2double& data){
    out << data.getValue();
    return out;
  }
  inline std::istream& operator>>(std::istream& in, su2double& data){
    double val;
    in >> val;
    data.setValue(val);
    return in;
  }
}

/* --- We need additional functions that are not defined yet --- */

inline su2double min(const su2double& a, const su2double& b){ return fmin(a,b);}
inline su2double max(const su2double& a, const su2double& b){ return fmax(a,b);}
inline su2double abs(const su2double&a){ return fabs(a);}
inline su2double atanh(const su2double& a){return 0.5*log(1+a)/log(1-a);}
