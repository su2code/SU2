/*!
 * \file complex_structure.hpp
 * \brief Headers for complex datatype definition.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
#include <complex>

class CComplexType;
typedef CComplexType su2double;

/*--- So the real() statement can be used even with the double type.
   Eases the implementation of some operators. ---*/
inline double real(const double& r);

/*!
 * \class CComplexType
 * \brief Class for defining the complex datatype for complex step gradient computation.
 * Based on complexify.h by Peter Sturdza.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 */

class CComplexType : public std::complex<double> {
public:
  CComplexType() : std::complex<double>() {}

  CComplexType(const double& d) : std::complex<double>(d) {}

  CComplexType(const double& r, const double& i) : std::complex<double>(r,i) {}

  CComplexType(const std::complex<double>& z) : std::complex<double>(z) {}

  CComplexType(const std::complex<float>& z) : std::complex<double>(z) {}

  operator double() { return this->real();}

  operator int() { return int(this->real());}

  operator short() { return short(this->real());}

  /*--- To get rid of some ambiguities, we need to reimplement some operators
   * (although they are already implemented for std::complex). ---*/

  /*--- Comparison operators (they are defined by comparing only the real parts
   * as we assume a very small imag. step ( < 1e-50)---*/

  template<typename T, typename Z>
  friend bool operator==(const Z&, const T&);

  template<typename T, typename Z>
  friend bool operator!=(const Z&, const T&);

  template<typename T, typename Z>
  friend bool operator<(const Z&, const T&);

  template<typename T, typename Z>
  friend bool operator>(const Z&, const T&);

  template<typename T, typename Z>
  friend bool operator<=(const Z&, const T&);

  template<typename T, typename Z>
  friend bool operator>=(const Z&, const T&);

  /*--- Basic arithmetic (some are templated to work with double, int, long int etc.) ---*/

  CComplexType operator+() const;

  CComplexType operator+(const CComplexType&) const;

  template<typename T>
  friend CComplexType operator+(const CComplexType&, const T&);

  template<typename T>
  friend CComplexType operator+(const T&, const CComplexType&);

  CComplexType operator-() const;

  CComplexType operator-(const CComplexType&) const;

  template<typename T>
  friend CComplexType operator-(const CComplexType&, const T&);

  template<typename T>
  friend CComplexType operator-(const T&, const CComplexType&);

  CComplexType operator*(const CComplexType&) const;

  template<typename T>
  friend CComplexType operator*(const CComplexType&, const T&);

  template<typename T>
  friend CComplexType operator*(const T&, const CComplexType&);

  CComplexType operator/(const CComplexType&) const;

  template<typename T>
  friend CComplexType operator/(const CComplexType&, const T&);

  template<typename T>
  friend CComplexType operator/(const T&, const CComplexType&);


  /*--- From <math.h> ---*/

  friend CComplexType sin(const CComplexType&);
  friend CComplexType sinh(const CComplexType&);
  friend CComplexType cos(const CComplexType&);
  friend CComplexType cosh(const CComplexType&);
  friend CComplexType tan(const CComplexType&);
  friend CComplexType tanh(const CComplexType&);
  friend CComplexType log10(const CComplexType&);
  friend CComplexType log(const CComplexType&);
  friend CComplexType sqrt(const CComplexType&);
  friend CComplexType exp(const CComplexType&);
  friend CComplexType pow(const CComplexType&, const CComplexType&);
  friend CComplexType pow(const CComplexType&, const double&);
  friend CComplexType pow(const CComplexType&, const int&);
  friend CComplexType pow(const double&, const CComplexType&);
  friend CComplexType pow(const int&, const CComplexType&);
  friend CComplexType min(const CComplexType& t, const CComplexType& z);
  template<typename T>
  friend CComplexType min(const T& t, const CComplexType& z);
  template<typename T>
  friend CComplexType min(const CComplexType& z, const T& t);
  friend CComplexType max(const CComplexType& t, const CComplexType& z);
  template<typename T>
  friend CComplexType max(const T& t, const CComplexType& z);
  template<typename T>
  friend CComplexType max(const CComplexType& z, const T& t);

  /*--- std::complex versions of these are not in standard library
   or they need to be redefined: (frexp, modf, and fmod have not been dealt with) ---*/

  friend CComplexType fabs(const CComplexType&);
  friend CComplexType asin(const CComplexType&);
  friend CComplexType acos(const CComplexType&);
  friend CComplexType atan(const CComplexType&);
  friend CComplexType atan2(const CComplexType&, const CComplexType&);
  friend CComplexType ceil(const CComplexType&);
  friend CComplexType floor(const CComplexType&);
  friend CComplexType ldexp(const CComplexType&, const int&);
  friend std::ostream& operator<< (std::ostream &out, const CComplexType &);

};
