/*!
 * \file CGemmStandard.hpp
 * \brief Standard class for carrying out a GEMM multiplication.
 *        The functions are in the <i>CGemmStandard.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
#include "mkl.h"
#else
#include "../blas_structure.hpp"
#endif

using namespace std;

/*!
 * \class CGemmStandard
 * \brief Standard class to carry out a GEMM multiplication.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CGemmStandard final : public CGemmBase {

private:

  int M;   /*!< \brief First matrix dimension of A and C in the gemm call. */
  int N;   /*!< \brief Second matrix dimension of B and C in the gemm call. */
  int K;   /*!< \brief First matrix dimension of B and second matrix dimension
                       of A in the gemm call. */

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  void *jitter;                 /*!< \brief Pointer to the data for the jitted gemm function. */
  dgemm_jit_kernel_t my_dgemm;  /*!< \brief Pointer to the function to carry out the jitted gemm call. */
#else
  CBlasStructure blasFunctions; /*!< \brief  The object to carry out the BLAS functionalities. */
#endif

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CGemmStandard() = delete;

  /*!
   * \overload
   * \param[in] val_M - First matrix dimension of A and C in the gemm call.
   * \param[in] val_N - Second matrix dimension of B and C in the gemm call.
   * \param[in] val_K - First matrix dimension of B and second matrix dimension
   *                    of A in the gemm call.
   */
  CGemmStandard(const int val_M, const int val_N, const int val_K);

  /*!
   * \brief Destructor.
   */
  ~CGemmStandard();
};
