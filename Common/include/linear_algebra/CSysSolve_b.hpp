/*!
 * \file CSysSolve_b.hpp
 * \brief Routines for the linear solver used in the reverse sweep of AD.
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

#include "../basic_types/datatype_structure.hpp"

#ifdef CODI_REVERSE_TYPE
template <class ScalarType>
struct CSysSolve_b {
  static void Solve_b(const su2double::Real* x, su2double::Real* x_b, size_t m, const su2double::Real* y,
                      const su2double::Real* y_b, size_t n, codi::ExternalFunctionUserData* d);
};
#endif
