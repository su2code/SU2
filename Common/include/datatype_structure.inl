/*!
 * \file datatype_structure.inl
 * \brief In-Line subroutines of the <i>datatype_structure.hpp</i> file.
 * \author T. Albring
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

/*--- Explicit cast functions ---*/

namespace SU2_TYPE{
  inline int Int(const su2double& data) {
    return int(SU2_TYPE::GetValue(data));
  }

  inline short Short(const su2double& data) {
    return short(SU2_TYPE::GetValue(data));
  }
}

/*--- Special handling of the sprint routine for non-primitive types. ---*/

#if  defined CODI_REVERSE_TYPE  || \
     defined CODI_FORWARD_TYPE

/*--- This objective is used for primitive types,
 where the output type of the getValue coincides with the input type. ---*/

template< typename IN > struct Impl_getValue {
  typedef IN OUT;
  static inline const OUT& getValue(const IN &value) {
    return value;
  }
};

/*--- This objective is used for non-primitive types,
 where the output type is double and the input type is su2double. ---*/

template<> struct Impl_getValue<su2double> {
  typedef double OUT;
  static inline OUT getValue(const su2double& value) {
    return SU2_TYPE::GetValue(value);
  }
};

/*--- Other objects are implemented in the corresponding header files of the datatypes.
 For example there may be an expression in the argument. ---*/


/*--- Terminating definition of sprintfOver ---*/

inline void sprintfOver(char * str, const char * format) {
  sprintf(str, "%s", format);
}


template <class A001>
inline int sprintfOver(char * str, const char * format, const A001 &a001) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001));
}

template <class A001, class A002>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002));
}

template <class A001, class A002, class A003>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003));
}

template <class A001, class A002, class A003, class A004>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004));
}

template <class A001, class A002, class A003, class A004, class A005>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005));
}

template <class A001, class A002, class A003, class A004, class A005, class A006>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015));
}

template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017));
}
template <class A001, class A002, class A003, class A004, class A005, class A006, class A007, class A008, class A009, class A010, class A011, class A012, class A013, class A014, class A015, class A016, class A017, class A018>
inline int sprintfOver(char * str, const char * format, const A001 &a001, const A002 &a002, const A003 &a003, const A004 &a004, const A005 &a005, const A006 &a006, const A007 &a007, const A008 &a008, const A009 &a009, const A010 &a010, const A011 &a011, const A012 &a012, const A013 &a013, const A014 &a014, const A015 &a015, const A016 &a016, const A017 &a017, const A018 &a018) {
  return sprintf(str, format, Impl_getValue<A001 >::getValue(a001), Impl_getValue<A002 >::getValue(a002), Impl_getValue<A003 >::getValue(a003), Impl_getValue<A004 >::getValue(a004), Impl_getValue<A005 >::getValue(a005), Impl_getValue<A006 >::getValue(a006), Impl_getValue<A007 >::getValue(a007), Impl_getValue<A008 >::getValue(a008), Impl_getValue<A009 >::getValue(a009), Impl_getValue<A010 >::getValue(a010), Impl_getValue<A011 >::getValue(a011), Impl_getValue<A012 >::getValue(a012), Impl_getValue<A013 >::getValue(a013), Impl_getValue<A014 >::getValue(a014), Impl_getValue<A015 >::getValue(a015), Impl_getValue<A016 >::getValue(a016), Impl_getValue<A017 >::getValue(a017), Impl_getValue<A018 >::getValue(a018));
}
#endif

#if defined CODI_REVERSE_TYPE
#include "datatypes/codi_reverse_structure.inl"
#elif defined CODI_FORWARD_TYPE
#include "datatypes/codi_forward_structure.inl"
#else
#include "datatypes/primitive_structure.inl"
#endif
