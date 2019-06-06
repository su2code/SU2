/*!
 * \file matrix_structure.inl
 * \brief In-Line subroutines of the <i>matrix_structure.hpp</i> file.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

template<class ScalarType>
inline void CSysMatrix<ScalarType>::SetValZero(void) { 
  if(NULL != matrix) {
	  for (unsigned long index = 0; index < nnz*nVar*nEqn; index++)
		matrix[index] = 0.0;
  }
}

template<class ScalarType>
template<class DstType, class SrcType>
inline DstType CSysMatrix<ScalarType>::ActiveAssign(const SrcType & val) const { return val; }

#ifdef CODI_REVERSE_TYPE
template<> template<>
inline passivedouble CSysMatrix<passivedouble>::ActiveAssign(const su2double & val) const { return SU2_TYPE::GetValue(val); }

template<> template<>
inline passivedouble CSysMatrix<su2double>::ActiveAssign(const su2double & val) const { return SU2_TYPE::GetValue(val); }
#endif

template<class ScalarType>
template<class DstType, class SrcType>
inline DstType CSysMatrix<ScalarType>::PassiveAssign(const SrcType & val) const {
#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
  return SU2_TYPE::GetValue(val);
#else
  return val;
#endif
}

template<class ScalarType>
inline ScalarType *CSysMatrix<ScalarType>::GetBlock(unsigned long block_i, unsigned long block_j) {

  for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
    if (col_ind[index] == block_j)
      return &(matrix[index*nVar*nEqn]);

  return NULL;
}

template<class ScalarType>
inline ScalarType CSysMatrix<ScalarType>::GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar, unsigned short jVar) {

  for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
    if (col_ind[index] == block_j)
      return matrix[index*nVar*nEqn+iVar*nEqn+jVar];

  return 0;
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SetBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[index*nVar*nEqn+iVar*nEqn+jVar] = PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SetBlock(unsigned long block_i, unsigned long block_j, OtherType *val_block) {
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[index*nVar*nEqn+iVar*nEqn+jVar] = PassiveAssign<ScalarType,OtherType>(val_block[iVar*nVar+jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::AddBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[index*nVar*nEqn+iVar*nEqn+jVar] += PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SubtractBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[index*nVar*nEqn+iVar*nEqn+jVar] -= PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::AddVal2Diag(unsigned long block_i, OtherType val_matrix) {
  
  unsigned long iVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    if (col_ind[index] == block_i) {	// Only elements on the diagonal
      for (iVar = 0; iVar < nVar; iVar++)
        matrix[index*nVar*nVar+iVar*nVar+iVar] += PassiveAssign<ScalarType,OtherType>(val_matrix);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SetVal2Diag(unsigned long block_i, OtherType val_matrix) {
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    if (col_ind[index] == block_i) {	// Only elements on the diagonal
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          matrix[index*nVar*nVar+iVar*nVar+jVar] = 0.0;
      
      for (iVar = 0; iVar < nVar; iVar++)
        matrix[index*nVar*nVar+iVar*nVar+iVar] = PassiveAssign<ScalarType,OtherType>(val_matrix);
      
      break;
    }
  }
  
}

template<class ScalarType>
inline ScalarType *CSysMatrix<ScalarType>::GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j) {

  for (unsigned long index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++)
    if (col_ind_ilu[index] == block_j)
      return &(ILU_matrix[index*nVar*nEqn]);

  return NULL;
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block) {
  
  unsigned long iVar, index;
  
  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar*nEqn; iVar++)
        ILU_matrix[index*nVar*nEqn+iVar] = val_block[iVar];
      break;
    }
  }
  
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::SetBlockTransposed_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block) {

  unsigned long iVar, jVar, index;

  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[index*nVar*nEqn+iVar*nEqn+jVar] = val_block[jVar*nVar+iVar];
      break;
    }
  }

}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block) {
  
  for (unsigned long index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    if (col_ind_ilu[index] == block_j) {
      MatrixSubtraction(&ILU_matrix[index*nVar*nEqn], val_block, &ILU_matrix[index*nVar*nEqn]);
      break;
    }
  }
  
}

template<class T, bool alpha, bool beta, bool transp>
inline void gemv_impl(const unsigned long n, const T *a, const T *b, T *c) {
  /*---
   This is a templated version of GEMV with the constants as boolean
   template parameters so that they can be optimized away at compilation.
   This is still the traditional "row dot vector" method.
  ---*/
  unsigned long i, j;
  for (i = 0; i < n; i++) {
    if (!beta) c[i] = 0.0;
    for (j = 0; j < n; j++) {
      if (alpha) c[transp? j:i] += a[i*n+j] * b[transp? i:j];
      else       c[transp? j:i] -= a[i*n+j] * b[transp? i:j];
    }
  }
}

template<class T>
inline void gemm_impl(const unsigned long n, const T *a, const T *b, T *c) {
  /*--- Same deal as for GEMV but here only the type is templated. ---*/
  unsigned long i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      c[i*n+j] = 0.0;
      for (k = 0; k < n; k++)
        c[i*n+j] += a[i*n+k] * b[k*n+j];
    }
  }
}

#define __MATVECPROD_SIGNATURE__(TYPE,NAME) \
inline void CSysMatrix<TYPE>::NAME(const TYPE *matrix, const TYPE *vector, TYPE *product)

#define MATVECPROD_SIGNATURE(NAME) template<class ScalarType> __MATVECPROD_SIGNATURE__(ScalarType,NAME)

#if !defined(USE_MKL)
MATVECPROD_SIGNATURE( MatrixVectorProduct ) {
  /*---
   Without MKL (default) picture copying the body of gemv_impl
   here and resolving the conditionals at compilation.
  ---*/
  gemv_impl<ScalarType,true,false,false>(nVar, matrix, vector, product);
}

MATVECPROD_SIGNATURE( MatrixVectorProductAdd ) {
  gemv_impl<ScalarType,true,true,false>(nVar, matrix, vector, product);
}

MATVECPROD_SIGNATURE( MatrixVectorProductSub ) {
  gemv_impl<ScalarType,false,true,false>(nVar, matrix, vector, product);
}

MATVECPROD_SIGNATURE( MatrixVectorProductTransp ) {
  gemv_impl<ScalarType,true,true,true>(nVar, matrix, vector, product);
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::MatrixMatrixProduct(const ScalarType *matrix_a, const ScalarType *matrix_b, ScalarType *product) {
  gemm_impl<ScalarType>(nVar, matrix_a, matrix_b, product);
}
#else
MATVECPROD_SIGNATURE( MatrixVectorProduct ) {
  /*--- With MKL we use the just-in-time kernels instead of the naive implementation. ---*/
  MatrixVectorProductKernelBetaZero(MatrixVectorProductJitterBetaZero, const_cast<ScalarType*>(vector),
                                    const_cast<ScalarType*>(matrix), product );
}

MATVECPROD_SIGNATURE( MatrixVectorProductAdd ) {
  MatrixVectorProductKernelBetaOne(MatrixVectorProductJitterBetaOne, const_cast<ScalarType*>(vector),
                                   const_cast<ScalarType*>(matrix), product );
}

MATVECPROD_SIGNATURE( MatrixVectorProductSub ) {
  MatrixVectorProductKernelAlphaMinusOne(MatrixVectorProductJitterAlphaMinusOne, const_cast<ScalarType*>(vector),
                                         const_cast<ScalarType*>(matrix), product );
}

MATVECPROD_SIGNATURE( MatrixVectorProductTransp ) {
  MatrixVectorProductTranspKernelBetaOne(MatrixVectorProductTranspJitterBetaOne, const_cast<ScalarType*>(matrix),
                                         const_cast<ScalarType*>(vector), product );
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::MatrixMatrixProduct(const ScalarType *matrix_a, const ScalarType *matrix_b, ScalarType *product) {
  MatrixMatrixProductKernel(MatrixMatrixProductJitter, const_cast<ScalarType*>(matrix_a),
                            const_cast<ScalarType*>(matrix_b), product );
}
#ifdef CODI_REVERSE_TYPE
/*--- WHEN using MKL, AND compiling for AD, we need to specialize for su2double to avoid mixing incompatible types. ---*/
#define MATVECPROD_SPECIALIZATION(NAME) template<> __MATVECPROD_SIGNATURE__(su2double,NAME)
MATVECPROD_SPECIALIZATION( MatrixVectorProduct ) {
  gemv_impl<su2double,true,false,false>(nVar, matrix, vector, product);
}

MATVECPROD_SPECIALIZATION( MatrixVectorProductAdd ) {
  gemv_impl<su2double,true,true,false>(nVar, matrix, vector, product);
}

MATVECPROD_SPECIALIZATION( MatrixVectorProductSub ) {
  gemv_impl<su2double,false,true,false>(nVar, matrix, vector, product);
}

MATVECPROD_SPECIALIZATION( MatrixVectorProductTransp ) {
  gemv_impl<su2double,true,true,true>(nVar, matrix, vector, product);
}

template<>
inline void CSysMatrix<su2double>::MatrixMatrixProduct(const su2double *matrix_a, const su2double *matrix_b, su2double *product) {
  gemm_impl<su2double>(nVar, matrix_a, matrix_b, product);
}
#undef MATVECPROD_SPECIALIZATION
#endif // CODI_REVERSE_TYPE
#endif // USE_MKL

#undef MATVECPROD_SIGNATURE
#undef __MATVECPROD_SIGNATURE__

template<class ScalarType>
inline void CSysMatrix<ScalarType>::VectorSubtraction(const ScalarType *a, const ScalarType *b, ScalarType *c) {

  for(unsigned long iVar = 0; iVar < nVar; iVar++)
    c[iVar] = a[iVar] - b[iVar];
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::MatrixSubtraction(const ScalarType *a, const ScalarType *b, ScalarType *c) {

  for(unsigned long iVar = 0; iVar < nVar*nEqn; iVar++)
    c[iVar] = a[iVar] - b[iVar];
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::Gauss_Elimination(ScalarType* matrix, ScalarType* vec) {

  /*---
   This is a relatively large method to inline but maybe better
   code will be generated for the special case nVar=1 this way.
  ---*/

  if (nVar==1) {vec[0] /= matrix[0]; return;}

#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  LAPACKE_dgetrf( LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, mkl_ipiv );
  LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', nVar, 1, matrix, nVar, mkl_ipiv, vec, 1 );
#else
  int iVar, jVar, kVar, nvar = int(nVar);
  ScalarType weight;

  /*--- Transform system in Upper Matrix ---*/
  for (iVar = 1; iVar < nvar; iVar++) {
    for (jVar = 0; jVar < iVar; jVar++) {
      weight = matrix[iVar*nvar+jVar] / matrix[jVar*nvar+jVar];
      for (kVar = jVar; kVar < nvar; kVar++)
        matrix[iVar*nvar+kVar] -= weight*matrix[jVar*nvar+kVar];
      vec[iVar] -= weight*vec[jVar];
    }
  }

  /*--- Backwards substitution ---*/
  for (iVar = nvar-1; iVar >= 0; iVar--) {
    for (jVar = iVar+1; jVar < nvar; jVar++)
      vec[iVar] -= matrix[iVar*nvar+jVar]*vec[jVar];
    vec[iVar] /= matrix[iVar*nvar+iVar];
  }
#endif
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::MatrixInverse(const ScalarType *matrix, ScalarType *inverse) {

  /*---
   This is a generalization of Gaussian elimination for multiple rhs' (the basis vectors).
   We could call "Gauss_Elimination" multiple times or fully generalize it for multiple rhs,
   the performance of both routines would suffer in both cases without the use of exotic templating.
   And so it feels reasonable to have some duplication here.
  ---*/

  if (nVar==1) {inverse[0] = 1.0/matrix[0]; return;}

  int iVar, jVar, nvar = int(nVar);

  /*--- Initialize the inverse and make a copy of the matrix ---*/
  for (iVar = 0; iVar < nvar; iVar++) {
    for (jVar = 0; jVar < nvar; jVar++) {
      block[iVar*nvar+jVar] = matrix[iVar*nvar+jVar];
      inverse[iVar*nvar+jVar] = ScalarType(iVar==jVar); // identity
    }
  }

  /*--- Inversion ---*/
#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  LAPACKE_dgetrf( LAPACK_ROW_MAJOR, nVar, nVar, block, nVar, mkl_ipiv );
  LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', nVar, nVar, block, nVar, mkl_ipiv, inverse, nVar );
#else
  int kVar;
  ScalarType weight;

  /*--- Transform system in Upper Matrix ---*/
  for (iVar = 1; iVar < nvar; iVar++) {
    for (jVar = 0; jVar < iVar; jVar++)
    {
      weight = block[iVar*nvar+jVar] / block[jVar*nvar+jVar];

      for (kVar = jVar; kVar < nvar; kVar++)
        block[iVar*nvar+kVar] -= weight*block[jVar*nvar+kVar];

      /*--- at this stage "inverse" is lower triangular so not all cols need updating ---*/
      for (kVar = 0; kVar <= jVar; kVar++)
        inverse[iVar*nvar+kVar] -= weight*inverse[jVar*nvar+kVar];
    }
  }

  /*--- Backwards substitution ---*/
  for (iVar = nvar-1; iVar >= 0; iVar--)
  {
    for (jVar = iVar+1; jVar < nvar; jVar++)
      for (kVar = 0; kVar < nvar; kVar++)
        inverse[iVar*nvar+kVar] -= block[iVar*nvar+jVar] * inverse[jVar*nvar+kVar];

    for (kVar = 0; kVar < nvar; kVar++)
      inverse[iVar*nvar+kVar] /= block[iVar*nvar+iVar];
  }
#endif
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::InverseDiagonalBlock(unsigned long block_i, ScalarType *invBlock, bool transpose) {
  
  const ScalarType* mat = GetBlock(block_i, block_i);
  MatrixInverse(mat, invBlock);
  
  if (transpose) // swap off-diag
    for (unsigned long iVar = 0; iVar < nVar-1; ++iVar)
      for (unsigned long jVar = iVar+1; jVar < nVar; ++jVar) {
        ScalarType tmp = invBlock[iVar*nVar+jVar];
        invBlock[iVar*nVar+jVar] = invBlock[jVar*nVar+iVar];
        invBlock[jVar*nVar+iVar] = tmp;
      }
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::InverseDiagonalBlock_ILUMatrix(unsigned long block_i, ScalarType *invBlock) {

  const ScalarType* mat = GetBlock_ILUMatrix(block_i, block_i);
  MatrixInverse(mat, invBlock);
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::EnforceSolutionAtNode(const unsigned long node_i, const OtherType *x_i, CSysVector<OtherType> & b) {

  /*--- Both row and column associated with node i are eliminated (Block_ii = I and all else 0) to preserve eventual symmetry. ---*/
  /*--- The vector is updated with the product of column i by the known (enforced) solution at node i. ---*/

  unsigned long iPoint, iVar, jVar, index, mat_begin;

  /*--- Delete whole row first. ---*/
  for (index = row_ptr[node_i]*nVar*nVar; index < row_ptr[node_i+1]*nVar*nVar; ++index)
    matrix[index] = 0.0;

  /*--- Update b with the column product and delete column. ---*/
  for (iPoint = 0; iPoint < nPoint; ++iPoint) {
    for (index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; ++index) {
      if (col_ind[index] == node_i)
      {
        mat_begin = index*nVar*nVar;

        for(iVar = 0; iVar < nVar; ++iVar)
          for(jVar = 0; jVar < nVar; ++jVar)
            b[iPoint*nVar+iVar] -= matrix[mat_begin+iVar*nVar+jVar] * x_i[jVar];

        /*--- If on diagonal, set diagonal of block to 1, else delete block. ---*/
        if (iPoint == node_i)
          for (iVar = 0; iVar < nVar; ++iVar) matrix[mat_begin+iVar*(nVar+1)] = 1.0;
        else
          for (iVar = 0; iVar < nVar*nVar; iVar++) matrix[mat_begin+iVar] = 0.0;
      }
    }
  }

  /*--- Set know solution in rhs vector. ---*/
  for (iVar = 0; iVar < nVar; iVar++) b[node_i*nVar+iVar] = x_i[iVar];

}

template<class ScalarType>
inline CSysMatrixVectorProduct<ScalarType>::CSysMatrixVectorProduct(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;  
}

template<class ScalarType>
inline void CSysMatrixVectorProduct<ScalarType>::operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CSysMatrixVectorProduct::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->MatrixVectorProduct(u, v, geometry, config);
}

template<class ScalarType>
inline CSysMatrixVectorProductTransposed<ScalarType>::CSysMatrixVectorProductTransposed(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;
}

template<class ScalarType>
inline void CSysMatrixVectorProductTransposed<ScalarType>::operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CSysMatrixVectorProduct::operator()(const CSysVector &, CSysVector &): " << endl;
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->MatrixVectorProductTransposed(u, v, geometry, config);
}

template<class ScalarType>
inline CJacobiPreconditioner<ScalarType>::CJacobiPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;  
}

template<class ScalarType>
inline void CJacobiPreconditioner<ScalarType>::operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CJacobiPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeJacobiPreconditioner(u, v, geometry, config);
}

template<class ScalarType>
inline CILUPreconditioner<ScalarType>::CILUPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;
}

template<class ScalarType>
inline void CILUPreconditioner<ScalarType>::operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CILUPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl;
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeILUPreconditioner(u, v, geometry, config);
}

template<class ScalarType>
inline CLU_SGSPreconditioner<ScalarType>::CLU_SGSPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
      geometry = geometry_ref;
  config = config_ref;
}

template<class ScalarType>
inline void CLU_SGSPreconditioner<ScalarType>::operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CLU_SGSPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeLU_SGSPreconditioner(u, v, geometry, config);
}

template<class ScalarType>
inline CLineletPreconditioner<ScalarType>::CLineletPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;  
}

template<class ScalarType>
inline void CLineletPreconditioner<ScalarType>::operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CLineletPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeLineletPreconditioner(u, v, geometry, config);
}
