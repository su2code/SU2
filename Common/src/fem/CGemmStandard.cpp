/*!
 * \file CGemmStandard.cpp
 * \brief Functions for the class CGemmStandard.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

CGemmStandard::CGemmStandard(const int val_M, const int val_N, const int val_K)
  : CGemmBase() {

  /*--- Copy the arguments into the member variables. ---*/
  M = val_M;
  N = val_N;
  K = val_K;

  /*--- Create the padded value of M. ---*/
  MP = CFEMStandardElementBase::PaddedValue(M);

  /*--- Check if the jitted GEMM call must be created. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- Create the GEMM Kernel and check if it went okay. ---*/
  passivedouble alpha = 1.0;
  passivedouble beta  = 0.0;
  mkl_jit_status_t status = mkl_jit_create_dgemm(&jitterFace, MKL_COL_MAJOR, MKL_NOTRANS,
                                                 MKL_NOTRANS, MP, N, K, alpha, MP, K, beta, MP);

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

void CGemmStandard::DOFs2Int(ColMajorMatrix<passivedouble> &basis,
                             const int                     nVar,
                             ColMajorMatrix<su2double>     &dataDOFs,
                             ColMajorMatrix<su2double>     &dataInt,
                             const CConfig                 *config) {

  /*--- Check if the jitted gemm call can be used. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- A couple of checks to see if this function is used correctly. ---*/
  assert(nVar == N);
  assert(MP   == basis.rows());
  assert(MP   == dataInt.rows());
  assert(K    == dataDOFs.rows());
  assert(K    == basis.cols());
  assert(nVar == dataInt.cols());
  assert(nVar == dataDOFs.cols());

/*--- Carry out the timing, if desired and call the jitted gemm function. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

  gemmFace(jitterFace, basis.data(), dataDOFs.data(), dataInt.data());

#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif

#else

  /*--- Use the interface to the more standard BLAS functionality. ---*/
  blasFunctions.gemm(MP, nVar, K, basis.rows(), dataDOFs.rows(), dataInt.rows(),
                     basis.data(), dataDOFs.data(), dataInt.data(), config);
#endif
}
