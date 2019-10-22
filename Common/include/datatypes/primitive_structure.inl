/*!
 * \file primitive_structure.inl
 * \brief Inline subroutines for <i>datatype_structure.hpp<i>.
 * \author T. Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
  inline void SetValue(double& data, const double &val) {data = val;}

  inline double GetValue(const double& data) { return data;}

  inline void SetSecondary(double& data, const double &val) {}

  inline double GetDerivative(const double& data) { return 0.0;}

  inline double GetSecondary(const double& data) { return 0.0;}

  inline void SetDerivative(double &data, const double &val) {}
}
