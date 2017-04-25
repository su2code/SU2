/*!
 * \file complex_structure.inl
 * \brief In-Line subroutines of the <i>datatype_structure.hpp</i> file.
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

namespace SU2_TYPE{
  inline void SetValue(su2double& data, const double &val) {
    data = su2double(val, data.imag());
  }

  inline double GetValue(const su2double& data) {
    return data.real();
  }

  inline void SetSecondary(su2double& data, const double &val) {
    data = su2double(data.real(), val);
  }

  inline double GetSecondary(const su2double& data) {
    return data.imag();
  }

  inline double GetDerivative(const su2double& data) {
    return data.imag()/1e-50;
  }

  inline void SetDerivative(su2double& data, const double &val) {
    data = su2double(data.real(), val*1e-50);
  }
}
inline double real(const double& r) {
  return r;
}

template<typename T, typename Z>
inline bool operator==(const Z& lhs, const T& rhs) {
  return real(lhs) == real(rhs);
}

template<typename T, typename Z>
inline bool operator!=(const Z& lhs, const T& rhs) {
  return real(lhs) != real(rhs);
}

template<typename T, typename Z>
inline bool operator<(const Z& lhs, const T& rhs) {
  return real(lhs) < real(rhs);
}

template<typename T, typename Z>
inline bool operator>(const Z& lhs, const T& rhs) {
  return real(lhs) > real(rhs);
}

template<typename T, typename Z>
inline bool operator<=(const Z& lhs, const T& rhs) {
  return real(lhs) <= real(rhs);
}

template<typename T, typename Z>
inline bool operator>=(const Z& lhs, const T& rhs) {
  return real(lhs) >= real(rhs);
}

inline CComplexType min(const CComplexType& t, const CComplexType& z) {
  if (real(t) < real(z)) return t;
  else return z;
}

template<typename T>
inline CComplexType min(const T& t, const CComplexType& z) {
  if (real(t) < real(z)) return CComplexType(t, 0.0);
  else return z;
}

template<typename T>
inline CComplexType min(const CComplexType& z, const T& t) {
  if (real(t) < real(z)) return CComplexType(t, 0.0);
  else return z;
}

inline CComplexType max(const CComplexType& t, const CComplexType& z) {
  if (real(t) > real(z)) return t;
  else return z;
}

template<typename T>
inline CComplexType max(const T& t, const CComplexType& z) {
  if (real(t) > real(z)) return CComplexType(t, 0.0);
  else return z;
}

template<typename T>
inline CComplexType max(const CComplexType& z, const T& t) {
  if (real(t) > real(z)) return CComplexType(t, 0.0);
  else return z;
}

inline CComplexType CComplexType::operator+() const{
  return +std::complex<double>(*this);
}

inline CComplexType CComplexType::operator+(const CComplexType& z) const{
  return std::complex<double>(*this)+std::complex<double>(z);
}

template<typename T>
inline CComplexType operator+(const CComplexType& z, const T& t) {
  return std::complex<double>(z)+std::complex<double>(t);
}

template<typename T>
inline CComplexType operator+(const T& t, const CComplexType& z) {
  return std::complex<double>(z)+std::complex<double>(t);
}

inline CComplexType CComplexType::operator-() const{
  return -std::complex<double>(*this);
}

inline CComplexType CComplexType::operator-(const CComplexType& z) const{
  return std::complex<double>(*this)-std::complex<double>(z);
}

template<typename T>
inline CComplexType operator-(const CComplexType& z, const T& t) {
  return std::complex<double>(z)-std::complex<double>(t);
}

template<typename T>
inline CComplexType operator-(const T& t, const CComplexType& z) {
  return std::complex<double>(t)-std::complex<double>(z);
}

inline CComplexType CComplexType::operator*(const CComplexType& z) const{
  return std::complex<double>(*this)*std::complex<double>(z);
}

template<typename T>
inline CComplexType operator*(const CComplexType& z, const T& t) {
  return std::complex<double>(z)*std::complex<double>(t);
}

template<typename T>
inline CComplexType operator*(const T& t, const CComplexType& z) {
  return std::complex<double>(z)*std::complex<double>(t);
}

inline CComplexType CComplexType::operator/(const CComplexType& z) const{
  return std::complex<double>(*this)/std::complex<double>(z);
}

template<typename T>
inline CComplexType operator/(const CComplexType& z, const T& t) {
  return std::complex<double>(z)/std::complex<double>(t);
}

template<typename T>
inline CComplexType operator/(const T& t, const CComplexType& z) {
  return std::complex<double>(t)/std::complex<double>(z);
}

inline CComplexType sin(const CComplexType& z) {
  return sin(std::complex<double>(z));
}

inline CComplexType sinh(const CComplexType& z) {
  return sinh(std::complex<double>(z));
}

inline CComplexType cos(const CComplexType& z) {
  return cos(std::complex<double>(z));
}

inline CComplexType cosh(const CComplexType& z) {
  return cosh(std::complex<double>(z));
}

#ifdef __GNUC__ // bug in gcc ?? get segv w/egcs-2.91.66 and 2.95.2
inline CComplexType tan(const CComplexType& z) {
  return sin(std::complex<double>(z))/cos(std::complex<double>(z));
}

inline CComplexType tanh(const CComplexType& z) {
  return sinh(std::complex<double>(z))/cosh(std::complex<double>(z));
}

inline CComplexType log10(const CComplexType& z) {
  return log(std::complex<double>(z))/log(10.);
}
#else
inline CComplexType tan(const CComplexType& z) {
  return tan(std::complex<double>(z));
}

inline CComplexType tanh(const CComplexType& z) {
  return tanh(std::complex<double>(z));
}

inline CComplexType log10(const CComplexType& z) {
  return log10(std::complex<double>(z));
}
#endif

inline CComplexType log(const CComplexType& z) {
  return log(std::complex<double>(z));
}

inline CComplexType sqrt(const CComplexType& z) {
  return sqrt(std::complex<double>(z));
}

inline CComplexType exp(const CComplexType& z) {
  return exp(std::complex<double>(z));
}

inline CComplexType pow(const CComplexType& a, const CComplexType& b) {
  return pow(std::complex<double>(a),std::complex<double>(b));
}

inline CComplexType pow(const CComplexType& a, const double& b) {
  return pow(std::complex<double>(a),b);
}

inline CComplexType pow(const CComplexType& a, const int& b) {
  return pow(std::complex<double>(a),double(b));
}

inline CComplexType pow(const double& a, const CComplexType& b) {
  return pow(a,std::complex<double>(b));
}

inline CComplexType pow(const int& a, const CComplexType& b) {
  return pow(double(a),std::complex<double>(b));
}

inline CComplexType fabs(const CComplexType& z) {
  return (real(z)<0.0) ? -z:z;
}

#define surr_TEENY (1.e-24) /* machine zero compared to nominal magnitude of
                               the real part */

inline CComplexType asin(const CComplexType& z) {
  // derivative trouble if imag(z) = +/- 1.0
  return CComplexType(asin(real(z)),imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

inline CComplexType acos(const CComplexType& z) {
  // derivative trouble if imag(z) = +/- 1.0
  return CComplexType(acos(real(z)),-imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

#undef surr_TEENY

inline CComplexType atan(const CComplexType& z) {
  return CComplexType(atan(real(z)),imag(z)/(1.0+real(z)*real(z)));
}

inline CComplexType atan2(const CComplexType& z1, const CComplexType& z2) {
  return CComplexType(atan2(real(z1),real(z2)),
              (real(z2)*imag(z1)-real(z1)*imag(z2))
              /(real(z1)*real(z1)+real(z2)*real(z2)));
}

inline CComplexType ceil(const CComplexType& z) {
  return CComplexType(ceil(real(z)),0.);
}

inline CComplexType floor(const CComplexType& z) {
  return CComplexType(floor(real(z)),0.);
}

inline CComplexType ldexp(const CComplexType& z, const int& i) {
  return CComplexType(ldexp(real(z),i),ldexp(imag(z),i));
}

inline std::ostream& operator<< (std::ostream &out, const CComplexType & a) {
  out << real(a);
  return out;
}

