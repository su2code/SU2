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
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[index*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar*nVar+jVar];
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
  
  unsigned long iVar, jVar, index;
  
  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[index*nVar*nEqn+iVar*nEqn+jVar] -= val_block[iVar*nVar+jVar];
      break;
    }
  }
  
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::MatrixVectorProduct(ScalarType *matrix, ScalarType *vector, ScalarType *product) {

#if defined(HAVE_MKL) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
  // NOTE: matrix/vector swapped due to column major kernel -- manual "CBLAS" setup.
  if (useMKL) 
  {
    MatrixVectorProductKernelBetaZero( MatrixVectorProductJitterBetaZero, vector, matrix, product );
    return;
  }
#endif
  
  unsigned short iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    product[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      product[iVar] += matrix[iVar*nVar+jVar] * vector[jVar];
    }
  }
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::MatrixMatrixProduct(ScalarType *matrix_a, ScalarType *matrix_b, ScalarType *product) {

#if defined(HAVE_MKL) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
  if (useMKL)
  {
    MatrixMatrixProductKernel( MatrixMatrixProductJitter, matrix_a, matrix_b, product );
    return;
  }
#endif
  
  unsigned short iVar, jVar, kVar;

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      product[iVar*nVar+jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++) {
        product[iVar*nVar+jVar] += matrix_a[iVar*nVar+kVar]*matrix_b[kVar*nVar+jVar];
      }
    }
  }
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::GetSubsBlock(ScalarType *c, ScalarType *a, ScalarType *b) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      c[iVar*nVar+jVar] = a[iVar*nVar+jVar] - b[iVar*nVar+jVar];
  
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::GetSubsVector(ScalarType *c, ScalarType *a, ScalarType *b) {
  
  unsigned long iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    c[iVar] = a[iVar] - b[iVar];
  
}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::Gauss_Elimination(ScalarType* matrix, ScalarType* vec) {

  short iVar, jVar, kVar;
  ScalarType weight, sum;

  if (nVar == 1) {
    vec[0] /= matrix[0];
  }
  else {

#if defined(HAVE_MKL) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
    if (useMKL) {
      // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
      lapack_int * ipiv = new lapack_int [ nVar ];
      LAPACKE_dgetrf( LAPACK_ROW_MAJOR, nVar, nVar, (double *)&block[0], nVar, ipiv );
      LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', nVar, 1, (double *)&block[0], nVar, ipiv, rhs, 1 );

      delete [] ipiv;
      return;
    }
#endif

    /*--- Transform system in Upper Matrix ---*/
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = matrix[iVar*nVar+jVar] / matrix[jVar*nVar+jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          matrix[iVar*nVar+kVar] -= weight*matrix[jVar*nVar+kVar];
        vec[iVar] -= weight*vec[jVar];
      }
    }

    /*--- Backwards substitution ---*/
    vec[nVar-1] = vec[nVar-1] / matrix[nVar*nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      sum = 0.0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        sum += matrix[iVar*nVar+jVar]*vec[jVar];
      vec[iVar] = (vec[iVar]-sum) / matrix[iVar*nVar+iVar];
    }
  }

}

template<class ScalarType>
inline void CSysMatrix<ScalarType>::InverseBlock(ScalarType *Block, ScalarType *invBlock) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++)
      aux_vector[jVar] = 0.0;
    aux_vector[iVar] = 1.0;
    
    /*--- Compute the i-th column of the inverse matrix ---*/
    Gauss_Elimination(Block, aux_vector);
    
    for (jVar = 0; jVar < nVar; jVar++)
      invBlock[jVar*nVar+iVar] = aux_vector[jVar];
  }
  
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
