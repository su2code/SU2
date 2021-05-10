/*!
 * \file CGemmStandard.cpp
 * \brief Functions for the class CGemmStandard.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#include "../../include/fem/CGemmStandard.hpp"
#include "../../include/fem/CFEMStandardElementBase.hpp"

/*----------------------------------------------------------------------------------*/
/*                  Public member functions of CGemmStandard.                       */
/*----------------------------------------------------------------------------------*/

CGemmStandard::CGemmStandard(const int val_M,
                             const int val_N,
                             const int val_K,
                             const int val_gemm_type)
  : CGemmBase() {

  /*--- Copy the arguments into the member variables. ---*/
  M = val_M;
  N = val_N;
  K = val_K;

  /*--- Create the padded value of M. ---*/
  MP = CFEMStandardElementBase::PaddedValue(M);

  /*--- Determine the gemm type and set some values accordingly. ---*/
  switch( val_gemm_type ) {
    case DOFS_TO_INT: {
      initZero = true;
      LDB      = K;
      break;
    }

    case INT_TO_DOFS: {
      initZero = false;
      LDB      = CFEMStandardElementBase::PaddedValue(K);
      break;
    }
  }

  /*--- Check if the jitted GEMM call must be created. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- Create the GEMM Kernel and check if it went okay. ---*/
  passivedouble alpha = 1.0;
  passivedouble beta  = initZero ? 0.0 : 1.0;
  mkl_jit_status_t status = mkl_jit_create_dgemm(&jitterFace, MKL_COL_MAJOR, MKL_NOTRANS,
                                                 MKL_NOTRANS, MP, N, K, alpha, MP, LDB, beta, MP);

  if(status == MKL_JIT_ERROR)
    SU2_MPI::Error(string("Jitted gemm kernel could not be created"), CURRENT_FUNCTION);

  /*--- Retrieve the function pointer to the DGEMM kernel
        void gemmFace(void*, double*, double*, double*). ---*/
  gemmFace = mkl_jit_get_dgemm_ptr(jitterFace);

#endif
}

CGemmStandard::~CGemmStandard() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterFace ) {
    mkl_jit_destroy(jitterFace);
    jitterFace = nullptr;
  }
#endif
}

void CGemmStandard::gemm(ColMajorMatrix<passivedouble> &A,
                         const int                     nVar,
                         ColMajorMatrix<su2double>     &B,
                         ColMajorMatrix<su2double>     &C,
                         const CConfig                 *config) {

  /*--- A couple of checks to see if this function is used correctly. ---*/
  assert(MP   == A.rows());
  assert(MP   == C.rows());
  assert(LDB  == B.rows());
  assert(K    == A.cols());
  assert(nVar == C.cols());
  assert(nVar == B.cols());

  /*--- Check if the jitted gemm call can be used. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- An additional check when the jitted gemm is used. ---*/
  assert(nVar == N);

  /*--- Carry out the timing, if desired and call the jitted gemm function. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

  gemmFace(jitterFace, A.data(), B.data(), C.data());

#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif

#else

  /*--- Use the interface to the more standard BLAS functionality. ---*/
  blasFunctions.gemm(MP, nVar, K, A.rows(), B.rows(), C.rows(), initZero,
                     A.data(), B.data(), C.data(), config);
#endif
}
