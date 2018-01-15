/*!
 * \file matrix_structure.inl
 * \brief In-Line subroutines of the <i>matrix_structure.hpp</i> file.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

template<class CalcType>
inline void TCSysMatrix<CalcType>::SetValZero(void) { 
  if(NULL != matrix) {
	  for (unsigned long index = 0; index < nnz*nVar*nEqn; index++)
		matrix[index] = 0.0;
  }
}
template<class CalcType>
inline TCSysMatrixVectorProduct<CalcType>::TCSysMatrixVectorProduct(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;  
}

template<class CalcType>
inline void TCSysMatrixVectorProduct<CalcType>::operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CSysMatrixVectorProduct::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->MatrixVectorProduct(u, v, geometry, config);
}

template<class CalcType>
inline TCSysMatrixVectorProductTransposed<CalcType>::TCSysMatrixVectorProductTransposed(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;
}

template<class CalcType>
inline void TCSysMatrixVectorProductTransposed<CalcType>::operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CSysMatrixVectorProduct::operator()(const CSysVector &, CSysVector &): " << endl;
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->MatrixVectorProductTransposed(u, v, geometry, config);
}


template<class CalcType>
inline TCJacobiPreconditioner<CalcType>::TCJacobiPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;  
}

template<class CalcType>
inline void TCJacobiPreconditioner<CalcType>::operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CJacobiPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeJacobiPreconditioner(u, v, geometry, config);
}

template<class CalcType>
inline TCILUPreconditioner<CalcType>::TCILUPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;
}

template<class CalcType>
inline void TCILUPreconditioner<CalcType>::operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CILUPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl;
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeILUPreconditioner(u, v, geometry, config);
}

template<class CalcType>
inline TCLU_SGSPreconditioner<CalcType>::TCLU_SGSPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
      geometry = geometry_ref;
  config = config_ref;
}

template<class CalcType>
inline void TCLU_SGSPreconditioner<CalcType>::operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CLU_SGSPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeLU_SGSPreconditioner(u, v, geometry, config);
}

template<class CalcType>
inline TCLineletPreconditioner<CalcType>::TCLineletPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref) {
  sparse_matrix = &matrix_ref;
  geometry = geometry_ref;
  config = config_ref;  
}

template<class CalcType>
inline void TCLineletPreconditioner<CalcType>::operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const {
  if (sparse_matrix == NULL) {
    cerr << "CLineletPreconditioner::operator()(const CSysVector &, CSysVector &): " << endl; 
    cerr << "pointer to sparse matrix is NULL." << endl;
    throw(-1);
  }
  sparse_matrix->ComputeLineletPreconditioner(u, v, geometry, config);
}
