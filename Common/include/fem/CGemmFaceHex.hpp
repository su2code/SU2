/*!
 * \file CGemmFaceHex.hpp
 * \brief Class for carrying out a GEMM multiplication for a face
 *        when the adjacent element is a hexahedron. In this
 *        case a tensor product is used.
 *        The functions are in the <i>CGemmFaceHex.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CGemmBase.hpp"
#include "../CConfig.hpp"

using namespace std;

/*!
 * \class CGemmFaceHex
 * \brief Class to carry out a GEMM multiplication for a face when the adjacent
 *        element is a hexahedron. In this case a tensor product is used.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CGemmFaceHex final : public CGemmBase {

private:

  int M;   /*!< \brief First tensor dimensions of A and C in the gemm call. */
  int N;   /*!< \brief Last tensor dimension of B and C in the gemm call. */
  int K;   /*!< \brief First tensor dimensions of B and last tensor dimensions
                       of A in the gemm call. */

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CGemmFaceHex() = delete;

  /*!
   * \overload
   * \param[in] val_M - First tensor dimensions of A and C in the gemm call.
   * \param[in] val_N - Last tensor dimension of B and C in the gemm call.
   * \param[in] val_K - First tensor dimensions of B and last tensor dimensions
   *                    of A in the gemm call.
   */
  CGemmFaceHex(const int val_M, const int val_N, const int val_K);

  /*!
   * \brief Destructor, nothing to be done.
   */
  ~CGemmFaceHex() = default;
};
