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
template<class OtherType>
inline void CSysMatrix<ScalarType>::SetBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SetBlock(unsigned long block_i, unsigned long block_j, OtherType *val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = PassiveAssign<ScalarType,OtherType>(val_block[iVar*nVar+jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::AddBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] += PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SubtractBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] -= PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::AddVal2Diag(unsigned long block_i, OtherType val_matrix) {
  
  unsigned long step = 0, iVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_i) {	// Only elements on the diagonal
      for (iVar = 0; iVar < nVar; iVar++)
        matrix[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+iVar] += PassiveAssign<ScalarType,OtherType>(val_matrix);
      break;
    }
  }
  
}

template<class ScalarType>
template<class OtherType>
inline void CSysMatrix<ScalarType>::SetVal2Diag(unsigned long block_i, OtherType val_matrix) {
  
  unsigned long step = 0, iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_i) {	// Only elements on the diagonal
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+jVar] = 0.0;
      
      for (iVar = 0; iVar < nVar; iVar++)
        matrix[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+iVar] = PassiveAssign<ScalarType,OtherType>(val_matrix);
      
      break;
    }
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
