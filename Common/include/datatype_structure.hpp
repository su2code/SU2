/*!
 * \file datatype_structure.hpp
 * \brief Headers for generalized datatypes.
 *        The subroutines and functions are in the <i>datatype_structure.cpp</i> file.
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

#include <ostream>
#include <cstdio>

#ifdef COMPLEX_TYPE
#include "complex_structure.hpp"

/* --- Define the complex datatype to be the base type used in SU2 --- */

typedef CComplexType su2double;

class CComplexTypeWrapper;

typedef CComplexTypeWrapper SU2_TYPE;

#else

#define SPRINTF sprintf

/* --- Define the double datatype to be the base type used in SU2 --- */

typedef double su2double;

class CTypeWrapper;

typedef CTypeWrapper SU2_TYPE;

#endif

/*!
 * \class CTypeWrapper
 * \brief Class for defining the datatype wrapper routines; this class features as a base class for
 * type interfaces for non-primitive dataypes e.g. used by AD, complex etc.
 * \author T. Albring
 * \version 3.2.9 "eagle"
 */
class CTypeWrapper{
public:
  static void SetPrimary(su2double& data, const double &val);

  static void SetSecondary(su2double& data, const double &val);

  static double GetPrimary(su2double &data);

  static double GetSecondary(su2double &data);
};

#ifdef COMPLEX_TYPE
class CComplexTypeWrapper : CTypeWrapper{
public:
  static void SetPrimary(su2double& data, const double &val);

  static void SetSecondary(su2double& data, const double &val);

  static double GetPrimary(su2double &data);

  static double GetSecondary(su2double &data);
};
#endif

#include "datatype_structure.inl"
