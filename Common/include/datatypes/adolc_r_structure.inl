/*!
 * \file adolc_r_structure.inl
 * \brief Inline subroutines for <i>datatype_structure.hpp<i>.
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

namespace AD{

  extern std::vector<double> inputVariables;
  extern std::vector<double> seedVector;
  extern int adjointVectorPos;
  extern double* adjointVector;

}
namespace SU2_TYPE{
  inline void SetPrimary(su2double& data, const double &val){data = val;}

  inline double GetPrimary(const su2double& data){return data.value();}

  inline void SetSecondary(su2double& data, const double &val){AD::seedVector.push_back(val);}

  inline double GetSecondary(const su2double& data){return AD::adjointVector[AD::adjointVectorPos++];}

  inline double GetDerivative(const su2double& data){return AD::adjointVector[AD::adjointVectorPos++];}

  inline void SetDerivative(su2double& data, const double &val){AD::seedVector.push_back(val);}
}
namespace AD{

  inline void RegisterInputVariable(su2double &data){
    inputVariables.push_back(data.value());
    data <<= data.value();
  }

  inline void RegisterOutputVariable(su2double& data){
    double temp;
    data >>= temp;
  }

  inline void StartRecording(){trace_on(1,1);}

  inline void StopRecording(){trace_off();}

  inline void ClearAdjoints(){
    delete [] adjointVector;
    adjointVector = NULL;
    seedVector.clear();
  }
  inline void ComputeAdjoint(){
    size_t tape_stats[STAT_SIZE];
    tapestats(1, tape_stats);

    double* input = new double[tape_stats[0]];
    double* gradient = new double[tape_stats[0]];
    double* adjoint = new double[tape_stats[1]];

    for (int i = 0; i < tape_stats[0]; ++i) {
      input[i] = inputVariables[i];
      gradient[i] = 0.0;
    }

    for(int i = 0; i < tape_stats[1]; ++i) {
      adjoint[i] = seedVector[i];
    }

    fos_reverse(1,tape_stats[1], tape_stats[0],adjoint, gradient);

    adjointVector = gradient;
    adjointVectorPos = 0;

    inputVariables.clear();
    seedVector.clear();

    delete [] input;
    delete [] adjoint;
  }
}

/* --- We need additional functions that are not defined yet --- */

inline su2double min(const su2double& a, const su2double& b){ return fmin(a,b);}
inline su2double max(const su2double& a, const su2double& b){ return fmax(a,b);}
inline su2double abs(const su2double&a){ return fabs(a);}
inline su2double atanh(const su2double& a){return 0.5*log(1+a)/log(1-a);}


template<> struct Impl_getValue<adub> {
  typedef double OUT; // Default implementation has the same output type as input type
  static inline OUT getValue(const adub &value) {
    return value.value();
  }
};
