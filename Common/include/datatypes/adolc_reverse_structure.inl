/*!
 * \file adolc_reverse_structure.inl
 * \brief Inline subroutines for <i>datatype_structure.hpp<i>.
 * \author T. Albring
 * \version 4.0.0 "Cardinal"
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

namespace AD{
  /* --- Stores a copy of the input variables (since they might be overwritten) --- */

  extern std::vector<double> inputVariables;

  /* --- Stores the seeding vector for the adjoint computation (adjoints of output variables) --- */

  extern std::vector<double> seedVector;

  /* --- The current position in the adjoint vector --- */

  extern int adjointVector_Position;

  /* --- Holds the adjoint values of the input variables after the adjoint computation --- */

  extern double* adjointVector;

}
namespace SU2_TYPE{
  inline void SetValue(su2double& data, const double &val){data = val;}

  inline double GetValue(const su2double& data){return data.value();}

  inline void SetSecondary(su2double& data, const double &val){AD::seedVector.push_back(val);}

  inline double GetSecondary(const su2double& data){return AD::adjointVector[AD::adjointVector_Position++];}

  inline double GetDerivative(const su2double& data){return AD::adjointVector[AD::adjointVector_Position++];}

  inline void SetDerivative(su2double& data, const double &val){AD::seedVector.push_back(val);}
}
namespace AD{

  inline void RegisterInput(su2double &data){inputVariables.push_back(data.value()); data <<= data.value();}

  inline void RegisterOutput(su2double& data){double temp; data >>= temp;}

  inline void StartRecording(){trace_on(1,1);}

  inline void StopRecording(){trace_off();}

}

/* --- We need additional functions that are not defined yet --- */

inline su2double min(const su2double& a, const su2double& b){ return fmin(a,b);}

inline su2double max(const su2double& a, const su2double& b){ return fmax(a,b);}

inline su2double abs(const su2double&a){ return fabs(a);}

inline su2double atanh(const su2double& a){return 0.5*log(1+a)/log(1-a);}

typedef decltype(std::cout) cout_type;

inline cout_type& operator << (cout_type& out, const adouble& data) {
  out << data.value();
  return out;
}

inline cout_type& operator << (cout_type& out, const adub& data) {
  out << data.value();
  return out;
}
template<> struct Impl_getValue<adub> {
  typedef double OUT; // Default implementation has the same output type as input type
  static inline OUT getValue(const adub &value) {
    return value.value();
  }
};
