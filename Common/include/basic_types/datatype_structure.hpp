/*!
 * \file datatype_structure.hpp
 * \brief Headers for generalized datatypes, defines an interface for AD types.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include <iostream>
#include <complex>
#include <cstdio>

#include "../code_config.hpp"
#include "ad_structure.hpp"

/*!
 * \namespace SU2_TYPE
 * \brief Namespace for defining the datatype wrapper routines, this acts as a base
 * class for type interfaces for non-primitive dataypes e.g. used by AD, complex etc.
 * \author T. Albring
 */
namespace SU2_TYPE {
/*!
 * \brief Set the (primitive) value of the datatype (needs to be implemented for each new type).
 * \param[in] data - The non-primitive datatype.
 * \param[in] val - The primitive value.
 */
void SetValue(su2double& data, const passivedouble& val);

/*!
 * \brief Set the secondary value of the datatype (needs to be implemented for each new type).
 * \param[in] data - The non-primitive datatype.
 * \param[in] val - The primitive value.
 */
void SetSecondary(su2double& data, const passivedouble& val);

/*!
 * \brief Get the (primitive) value of the datatype (needs to be specialized for active types).
 * \param[in] data - The non-primitive datatype.
 * \return The primitive value.
 */
passivedouble GetValue(const su2double& data);

/*!
 * \brief Get the secondary value of the datatype (needs to be implemented for each new type).
 * \param[in] data - The non-primitive datatype.
 * \return The primitive value.
 */
passivedouble GetSecondary(const su2double& data);

/*!
 * \brief Get the derivative value of the datatype (needs to be implemented for each new type).
 * \param[in] data - The non-primitive datatype.
 * \return The derivative value.
 */
passivedouble GetDerivative(const su2double& data);

/*!
 * \brief Set the derivative value of the datatype (needs to be implemented for each new type).
 * \param[in] data - The non-primitive datatype.
 * \param[in] val - The value of the derivative.
 */
void SetDerivative(su2double& data, const passivedouble& val);

/*--- Implementation of the above for the different types. ---*/

#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)

FORCEINLINE void SetValue(su2double& data, const passivedouble& val) { data.setValue(val); }

FORCEINLINE passivedouble GetValue(const su2double& data) { return data.getValue(); }

FORCEINLINE void SetSecondary(su2double& data, const passivedouble& val) { data.setGradient(val); }

FORCEINLINE void SetDerivative(su2double& data, const passivedouble& val) { data.setGradient(val); }

FORCEINLINE passivedouble GetSecondary(const su2double& data) { return data.getGradient(); }

FORCEINLINE passivedouble GetDerivative(const su2double& data) { return data.getGradient(); }

#else  // passive type, no AD

FORCEINLINE void SetValue(su2double& data, const passivedouble& val) { data = val; }

FORCEINLINE passivedouble GetValue(const su2double& data) { return data; }

FORCEINLINE void SetSecondary(su2double&, const passivedouble&) {}

FORCEINLINE passivedouble GetDerivative(const su2double&) { return 0.0; }

FORCEINLINE passivedouble GetSecondary(const su2double&) { return 0.0; }

FORCEINLINE void SetDerivative(su2double&, const passivedouble&) {}
#endif

/*!
 * \brief Get the passive value of any variable. For most types return directly,
 * specialize for su2double to call GetValue.
 * \note This is a struct instead of a function because the return type of the
 * su2double specialization changes.
 */
template <class T>
struct Passive {
  FORCEINLINE static T Value(const T& val) { return val; }
};
template <>
struct Passive<su2double> {
  FORCEINLINE static passivedouble Value(const su2double& val) { return GetValue(val); }
};

/*!
 * \brief Casts the primitive value to int (uses GetValue, already implemented for each type).
 * \param[in] data - The non-primitive datatype.
 * \return - The primary value casted to int.
 */
FORCEINLINE int Int(const su2double& data) { return static_cast<int>(SU2_TYPE::GetValue(data)); }

/*!
 * \brief Casts the primitive value to short (uses GetValue, already implemented for each type).
 * \param[in] data - The non-primitive datatype.
 * \return - The primary value casted to short.
 */
FORCEINLINE short Short(const su2double& data) { return static_cast<short>(SU2_TYPE::GetValue(data)); }

/*--- Special handling of the sprintf routine for non-primitive types. ---*/
/*--- Pass-through for built-in types. ---*/
template <class T, su2enable_if<std::is_trivial<T>::value> = 0>
FORCEINLINE const T& _printGetValue(const T& val) {
  return val;
}
/*--- Overload for expressions of active types. ---*/
template <class T, su2enable_if<!std::is_trivial<T>::value> = 0>
FORCEINLINE passivedouble _printGetValue(const T& val) {
  return val.getValue();
}

/*!
 * \brief Wrapper to sprintf to be able to print active types and AD expressions.
 * \note This is for compatibility with old code, stringstreams should be the preferred way to build strings.
 * \param[in] str - Target char buffer.
 * \param[in] format - Format string.
 * \param[in] args - Values to be printed to the string.
 */
template <class... Ts>
FORCEINLINE void sprintf(char* str, const char* format, Ts&&... args) {
  ::sprintf(str, format, SU2_TYPE::_printGetValue(args)...);
}
FORCEINLINE void sprintf(char* str, const char* literal) { ::sprintf(str, "%s", literal); }
#define SPRINTF SU2_TYPE::sprintf

}  // namespace SU2_TYPE
