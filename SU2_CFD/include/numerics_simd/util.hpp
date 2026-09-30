/*!
 * \file util.hpp
 * \brief Vector, matrix and index types bound to the SIMD Double/Int of CNumericsSIMD.hpp.
 * \author P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "CNumericsSIMD.hpp"
#include "../numerics/util.hpp"

template <size_t Size> using VectorInt = Vector<Int, Size>;
template <size_t Size> using VectorDbl = Vector<Double, Size>;

template <size_t Rows, size_t Cols = Rows> using MatrixInt = Matrix<Int, Rows, Cols>;
template <size_t Rows, size_t Cols = Rows> using MatrixDbl = Matrix<Double, Rows, Cols>;
