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

/* --- Define the complex datatype interface --- */

class CComplexTypeWrapper;

typedef CComplexTypeWrapper SU2_TYPE;

#else

#define SPRINTF sprintf

/* --- Define the double datatype to be the base type used in SU2 --- */

typedef double su2double;

/* --- Define the default datatype interface --- */

class CTypeWrapper;

typedef CTypeWrapper SU2_TYPE;

#endif


/*!
 * \class CTypeWrapper
 * \brief Class for defining the datatype wrapper routines; this class features as a base class for
 * type interfaces for non-primitive dataypes e.g. used by AD, complex etc. However, if we use the normal double datatype
 * then this base class is used which does essentially nothing.
 * \author T. Albring
 * \version 3.2.9 "eagle"
 */
class CTypeWrapper{
public:

  /*!
   * \brief Set the primary (primitive) value of the datatype.
   * \param[in] data - The non-primitive datatype.
   * \param[in] val - The primitive value.
   */
  static void SetPrimary(su2double& data, const double &val);

  /*!
   * \brief Set the secondary (primitive) value of the datatype (in this particular case it does nothing).
   * \param[in] data - The non-primitive datatype.
   * \param[in] val - The primitive value.
   */
  static void SetSecondary(su2double& data, const double &val);

  /*!
   * \brief Get the primary (primitive) value of the datatype.
   * \param[in] data - The non-primitive datatype.
   * \return The primitive value.
   */
  static double GetPrimary(su2double &data);

  /*!
   * \brief Get the secondary (primitive) value of the datatype (in this case it just return zero).
   * \param[in] data - The non-primitive datatype.
   * \return The primitive value.
   */
  static double GetSecondary(su2double &data);

  /*!
   * \brief Get the derivative (primitive) value of the datatype (in this case it just return zero).
   * \param[in] data - The non-primitive datatype.
   * \return The derivative value.
   */
  static double GetDerivative(su2double &data);

  /*!
   * \brief Set the derivative (primitive) value of the datatype (in this case it does nothing).
   * \param[in] data - The non-primitive datatype.
   * \param[in] val - The value of the derivative.
   */
  static void SetDerivative(su2double &data, const double &val);

};
/*!
 * \class CComplexTypeWrapper
 * \brief Class for defining the datatype wrapper routines for the complex datatype. It enables the extraction and setting of real and imag
 * values as well as the derivative values.
 * \author T. Albring
 * \version 3.2.9 "eagle"
 */
#ifdef COMPLEX_TYPE
class CComplexTypeWrapper : CTypeWrapper{
public:
  /*!
   * \brief Set the primary (primitive) value of the datatype (aka the real value).
   * \param[in] data - The non-primitive datatype.
   * \param[in] val - The primitive value.
   */
  static void SetPrimary(su2double& data, const double &val);

  /*!
   * \brief Set the secondary (primitive) value of the datatype (aka the imag value).
   * \param[in] data - The non-primitive datatype.
   * \param[in] val - The imag value.
   */
  static void SetSecondary(su2double& data, const double &val);

  /*!
   * \brief Get the primary (primitive) value of the datatype (aka the real value).
   * \param[in] data - The non-primitive datatype.
   * \return The primitive value.
   */
  static double GetPrimary(su2double &data);

  /*!
   * \brief Get the secondary (primitive) value of the datatype (aka the imag value).
   * \param[in] data - The non-primitive datatype.
   * \return The imag value.
   */
  static double GetSecondary(su2double &data);

  /*!
   * \brief Get the derivative (primitive) value of the datatype (imag divided by the stepsize).
   * \param[in] data - The non-primitive datatype.
   * \return The derivative value.
   */
  static double GetDerivative(su2double &data);

  /*!
   * \brief Set the derivative (primitive) value of the datatype (sets imag value to val times the stepsize).
   * \param[in] data - The non-primitive datatype.
   * \param[in] val - The value of the derivative.
   */
  static void SetDerivative(su2double &data, const double &val);

};
#endif

#include "datatype_structure.inl"
