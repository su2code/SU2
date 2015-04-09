/*!
 * \file complex_structure.inl
 * \brief In-Line subroutines of the <i>complex_structure.hpp</i> file.
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

inline double real(const double& r){
  return r;
}

template<typename T, typename Z>
inline bool operator==(const Z& lhs, const T& rhs){
  return real(lhs) == real(rhs);
}

template<typename T, typename Z>
inline bool operator!=(const Z& lhs, const T& rhs){
  return real(lhs) != real(rhs);
}

template<typename T, typename Z>
inline bool operator<(const Z& lhs, const T& rhs){
  return real(lhs) < real(rhs);
}

template<typename T, typename Z>
inline bool operator>(const Z& lhs, const T& rhs){
  return real(lhs) > real(rhs);
}

template<typename T, typename Z>
inline bool operator<=(const Z& lhs, const T& rhs){
  return real(lhs) <= real(rhs);
}

template<typename T, typename Z>
inline bool operator>=(const Z& lhs, const T& rhs){
  return real(lhs) >= real(rhs);
}

inline CComplexType min(const CComplexType& t, const CComplexType& z){
  if (real(t) < real(z)) return t;
  else return z;
}

template<typename T>
inline CComplexType min(const T& t, const CComplexType& z){
  if (real(t) < real(z)) return CComplexType(t, 0.0);
  else return z;
}

template<typename T>
inline CComplexType min(const CComplexType& z, const T& t){
  if (real(t) < real(z)) return CComplexType(t, 0.0);
  else return z;
}

inline CComplexType max(const CComplexType& t, const CComplexType& z){
  if (real(t) > real(z)) return t;
  else return z;
}

template<typename T>
inline CComplexType max(const T& t, const CComplexType& z){
  if (real(t) > real(z)) return CComplexType(t, 0.0);
  else return z;
}

template<typename T>
inline CComplexType max(const CComplexType& z, const T& t){
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
inline CComplexType operator+(const CComplexType& z, const T& t){
  return std::complex<double>(z)+std::complex<double>(t);
}

template<typename T>
inline CComplexType operator+(const T& t, const CComplexType& z){
  return std::complex<double>(z)+std::complex<double>(t);
}

inline CComplexType CComplexType::operator-() const{
  return -std::complex<double>(*this);
}

inline CComplexType CComplexType::operator-(const CComplexType& z) const{
  return std::complex<double>(*this)-std::complex<double>(z);
}

template<typename T>
inline CComplexType operator-(const CComplexType& z, const T& t){
  return std::complex<double>(z)-std::complex<double>(t);
}

template<typename T>
inline CComplexType operator-(const T& t, const CComplexType& z){
  return std::complex<double>(t)-std::complex<double>(z);
}

inline CComplexType CComplexType::operator*(const CComplexType& z) const{
  return std::complex<double>(*this)*std::complex<double>(z);
}

template<typename T>
inline CComplexType operator*(const CComplexType& z, const T& t){
  return std::complex<double>(z)*std::complex<double>(t);
}

template<typename T>
inline CComplexType operator*(const T& t, const CComplexType& z){
  return std::complex<double>(z)*std::complex<double>(t);
}

inline CComplexType CComplexType::operator/(const CComplexType& z) const{
  return std::complex<double>(*this)/std::complex<double>(z);
}

template<typename T>
inline CComplexType operator/(const CComplexType& z, const T& t){
  return std::complex<double>(z)/std::complex<double>(t);
}

template<typename T>
inline CComplexType operator/(const T& t, const CComplexType& z){
  return std::complex<double>(t)/std::complex<double>(z);
}

inline CComplexType sin(const CComplexType& z){
  return sin(std::complex<double>(z));
}

inline CComplexType sinh(const CComplexType& z){
  return sinh(std::complex<double>(z));
}

inline CComplexType cos(const CComplexType& z){
  return cos(std::complex<double>(z));
}

inline CComplexType cosh(const CComplexType& z){
  return cosh(std::complex<double>(z));
}

#ifdef __GNUC__ // bug in gcc ?? get segv w/egcs-2.91.66 and 2.95.2
inline CComplexType tan(const CComplexType& z){
  return sin(std::complex<double>(z))/cos(std::complex<double>(z));
}

inline CComplexType tanh(const CComplexType& z){
  return sinh(std::complex<double>(z))/cosh(std::complex<double>(z));
}

inline CComplexType log10(const CComplexType& z){
  return log(std::complex<double>(z))/log(10.);
}
#else
inline CComplexType tan(const CComplexType& z){
  return tan(std::complex<double>(z));
}

inline CComplexType tanh(const CComplexType& z){
  return tanh(std::complex<double>(z));
}

inline CComplexType log10(const CComplexType& z){
  return log10(std::complex<double>(z));
}
#endif

inline CComplexType log(const CComplexType& z){
  return log(std::complex<double>(z));
}

inline CComplexType sqrt(const CComplexType& z){
  return sqrt(std::complex<double>(z));
}

inline CComplexType exp(const CComplexType& z){
  return exp(std::complex<double>(z));
}

inline CComplexType pow(const CComplexType& a, const CComplexType& b){
  return pow(std::complex<double>(a),std::complex<double>(b));
}

inline CComplexType pow(const CComplexType& a, const double& b){
  return pow(std::complex<double>(a),b);
}

inline CComplexType pow(const CComplexType& a, const int& b){
  return pow(std::complex<double>(a),double(b));
}

inline CComplexType pow(const double& a, const CComplexType& b){
  return pow(a,std::complex<double>(b));
}

inline CComplexType pow(const int& a, const CComplexType& b){
  return pow(double(a),std::complex<double>(b));
}

inline CComplexType fabs(const CComplexType& z){
  return (real(z)<0.0) ? -z:z;
}

#define surr_TEENY (1.e-24) /* machine zero compared to nominal magnitude of
                               the real part */

inline CComplexType asin(const CComplexType& z){
  // derivative trouble if imag(z) = +/- 1.0
  return CComplexType(asin(real(z)),imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

inline CComplexType acos(const CComplexType& z){
  // derivative trouble if imag(z) = +/- 1.0
  return CComplexType(acos(real(z)),-imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

#undef surr_TEENY

inline CComplexType atan(const CComplexType& z){
  return CComplexType(atan(real(z)),imag(z)/(1.0+real(z)*real(z)));
}

inline CComplexType atan2(const CComplexType& z1, const CComplexType& z2){
  return CComplexType(atan2(real(z1),real(z2)),
              (real(z2)*imag(z1)-real(z1)*imag(z2))
              /(real(z1)*real(z1)+real(z2)*real(z2)));
}

inline CComplexType ceil(const CComplexType& z){
  return CComplexType(ceil(real(z)),0.);
}

inline CComplexType floor(const CComplexType& z){
  return CComplexType(floor(real(z)),0.);
}

inline CComplexType ldexp(const CComplexType& z, const int& i){
  return CComplexType(ldexp(real(z),i),ldexp(imag(z),i));
}

inline std::ostream& operator<< (std::ostream &out, const CComplexType & a){
  out << real(a);
  return out;
}

/* --- The following functions are necessary for the handling of the sprintf function --- */

template< typename IN > struct Impl_getValue {
  typedef IN OUT; // Default implementation has the same output type as input type
  static inline const OUT& getValue(const IN &value) {
    return value;
  }
};

template< typename IN, int n> struct Impl_getValue<IN[n]> {
  typedef IN OUT; // Default implementation has the same output type as input type
  static inline OUT* getValue(IN* value) {
    return value;
  }

  static inline const OUT* getValue(const IN* value) {
      return value;
    }
};

static inline void setValue(double &x, const double& value) {
  x = value;
}

template <typename IN>
static inline typename Impl_getValue<IN>::OUT getValue(const IN& x) {
  return Impl_getValue<IN>::getValue(x);
}

template<> struct Impl_getValue<CComplexType> {
  typedef double OUT;
  static inline OUT getValue(const CComplexType& value) {
    return real(value);
  }
};
inline void sprintfOver(char * str, const char * format) {
  sprintf(str, format);
}
template <class A001>
inline void sprintfOver(char * str, const char * format, const A001 &a001) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001));
}
template <class A001, class A002>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002));
}
template <class A001, class A002, class A003>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003));
}
template <class A001, class A002, class A003, class A004>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004));
}
template <class A001, class A002, class A003, class A004, class A005>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005));
}
template <class A001, class A002, class A003, class A004, class A005, class A006>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019, class A020>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019, const A020 &a020) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019), Impl_getValue<A020 >::getValue(a020));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019, class A020, class A021>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019, const A020 &a020, const A021 &a021) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019), Impl_getValue<A020 >::getValue(a020), Impl_getValue<A021 >::getValue(a021));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019, class A020, class A021, class A022>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019, const A020 &a020, const A021 &a021, const A022 &a022) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019), Impl_getValue<A020 >::getValue(a020), Impl_getValue<A021 >::getValue(a021), Impl_getValue<A022 >::getValue(a022));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019, class A020, class A021, class A022, class A023>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019, const A020 &a020, const A021 &a021, const A022 &a022, const A023 &a023) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019), Impl_getValue<A020 >::getValue(a020), Impl_getValue<A021 >::getValue(a021), Impl_getValue<A022 >::getValue(a022), Impl_getValue<A023 >::getValue(a023));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019, class A020, class A021, class A022, class A023, class A024>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019, const A020 &a020, const A021 &a021, const A022 &a022, const A023 &a023, const A024 &a024) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019), Impl_getValue<A020 >::getValue(a020), Impl_getValue<A021 >::getValue(a021), Impl_getValue<A022 >::getValue(a022), Impl_getValue<A023 >::getValue(a023), Impl_getValue<A024 >::getValue(a024));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018, class A019, class A020, class A021, class A022, class A023, class A024, class A025>
inline void sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018, const A019 &a019, const A020 &a020, const A021 &a021, const A022 &a022, const A023 &a023, const A024 &a024, const A025 &a025) {
  sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018), Impl_getValue<A019 >::getValue(a019), Impl_getValue<A020 >::getValue(a020), Impl_getValue<A021 >::getValue(a021), Impl_getValue<A022 >::getValue(a022), Impl_getValue<A023 >::getValue(a023), Impl_getValue<A024 >::getValue(a024), Impl_getValue<A025 >::getValue(a025));
}
