/*!
 * \file sparse_structure.inl
 * \brief In-Line subroutines of the <i>sparse_structure.hpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.1
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#pragma once

inline void CSparseMatrix::SetValZero(void) { 
	for (unsigned long index = 0; index < nnz*nVar*nVar; index++) 
		val[index] = 0.0; 
}

inline void CSparseMatrix::ScaleVals(double val_scale) { 
	for (unsigned long index = 0; index < nnz*nVar*nVar; index++) 
		val[index] *= val_scale; 
}


inline CSparseMatrixVectorProduct::CSparseMatrixVectorProduct(CSparseMatrix & matrix_ref) {
  sparse_matrix = &matrix_ref;
}

inline void CSparseMatrixVectorProduct::operator()(const CSysVector & u, CSysVector & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CSparseMatrixVectorProduct::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->MatrixVectorProduct(u, v);
}

inline CJacobiPreconditioner::CJacobiPreconditioner(CSparseMatrix & matrix_ref) {
  sparse_matrix = &matrix_ref;
}

inline void CJacobiPreconditioner::operator()(const CSysVector & u, CSysVector & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CJacobiPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeJacobiPreconditioner(u, v);
}

inline CLineletPreconditioner::CLineletPreconditioner(CSparseMatrix & matrix_ref) {
  sparse_matrix = &matrix_ref;
}

inline void CLineletPreconditioner::operator()(const CSysVector & u, CSysVector & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CJacobiPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeLineletPreconditioner(u, v);
}

inline void CIdentityPreconditioner::operator()(const CSysVector & u, CSysVector & v) const {
  v = u;
}

inline CSparseMatrixSolMPI::CSparseMatrixSolMPI(CSparseMatrix & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;
}

inline void CSparseMatrixSolMPI::operator()(CSysVector & u) const {
  if (sparse_matrix == NULL) {
    cerr << "CSparseMatrixSolMPI::operator()(CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->SendReceive_Solution(u, geometry, config);
}
