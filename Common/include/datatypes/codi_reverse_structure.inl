/*!
 * \file codi_reverse_structure.inl
 * \brief Inline subroutines for <i>datatype_structure.hpp<i>.
 * \author T. Albring
 * \version 7.0.4 "Blackbird"
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

#pragma once

namespace SU2_TYPE{
  inline void SetValue(su2double& data, const double &val) {data.setValue(val);}

  inline double GetValue(const su2double& data) { return data.getValue();}

  inline void SetSecondary(su2double& data, const double &val) {data.setGradient(val);}

  inline double GetSecondary(const su2double& data) { return AD::globalTape.getGradient(AD::inputValues[AD::adjointVectorPosition++]);}

  inline double GetDerivative(const su2double& data) { return AD::globalTape.getGradient(AD::inputValues[AD::adjointVectorPosition++]);}

  inline void SetDerivative(su2double& data, const double &val) {data.setGradient(val);}
}

/*--- Object for the definition of getValue used in the printfOver definition.
 * Necessary for cases where the argument of sprintfOver is an expression, e.g:
 * SPRINTF("Residual: %d", log10(Residual)) ---*/

template<class A> struct Impl_getValue<codi::Expression<double, A> > {
  typedef double OUT;
  static inline OUT getValue(const codi::Expression<double, A> &value) {
    return value.getValue();
  }
};

