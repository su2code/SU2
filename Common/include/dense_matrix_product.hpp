/*!
 * \file dense_matrix_product.hpp
 * \brief Include files and headers of the functions to carry out a dense
 *        matrix product. The functions are in the
 *        <i>dense_matrix_product.cpp</i> file.
 * \author E. van der Weide
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "datatype_structure.hpp"

/* LIBXSMM include files, if supported. */
#ifdef HAVE_LIBXSMM
#include "libxsmm.h"
#endif

/* MKL or BLAS include files, if supported. */
#ifdef HAVE_MKL
#include "mkl.h"
#elif HAVE_CBLAS
#include "cblas.h"
#endif

/*-----------------------------------------------------------------------------*/
/*---                          Function prototype.                          ---*/
/*-----------------------------------------------------------------------------*/

/*!
 * \brief Function, which carries out a dense matrix product.
 * \param[in]  M  - Number of rows of A and C.
 * \param[in]  N  - Number of columns of B and C.
 * \param[in]  K  - Number of columns of A and number of rows of B.
 * \param[in]  A  - Input matrix in the multiplication.
 * \param[in]  B  - Input matrix in the multiplication.
 * \param[out] C  - Result of the matrix product A*B.
 */
void DenseMatrixProduct(const int M,        const int N,        const int K,
                        const su2double *A, const su2double *B, su2double *C);
