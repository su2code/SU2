/*!
 * \file matrix_structure.cpp
 * \brief Main subroutines for doing the sparse structures
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 6.0.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/matrix_structure.hpp"

template <class CalcType>
TCSysMatrix<CalcType>::TCSysMatrix(void) {
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();
  
  /*--- Array initialization ---*/

  matrix            = NULL;
  ILU_matrix        = NULL;
  row_ptr           = NULL;
  col_ind           = NULL;
  row_ptr_ilu       = NULL;
  col_ind_ilu       = NULL;
  block             = NULL;
  prod_block_vector = NULL;
  prod_row_vector   = NULL;
  aux_vector        = NULL;
  sum_vector        = NULL;
  invM              = NULL;
  block_weight      = NULL;
  block_inverse     = NULL;

  /*--- Linelet preconditioner ---*/
  
  LineletBool     = NULL;
  LineletPoint    = NULL;
  UBlock          = NULL;
  invUBlock       = NULL;
  LBlock          = NULL;
  yVector         = NULL;
  zVector         = NULL;
  rVector         = NULL;
  LFBlock         = NULL;
  LyVector        = NULL;
  FzVector        = NULL;
  max_nElem       = 0;
  ilu_fill_in = 0;
  
}

template <class CalcType>
TCSysMatrix<CalcType>::~TCSysMatrix(void) {
  
  unsigned long iElem;

  /*--- Memory deallocation ---*/
  
  if (matrix != NULL)             delete [] matrix;
  if (ILU_matrix != NULL)         delete [] ILU_matrix;
  if (row_ptr != NULL)            delete [] row_ptr;
  if (col_ind != NULL)            delete [] col_ind;

  if (ilu_fill_in != 0) {
    if (row_ptr_ilu != NULL) delete [] row_ptr_ilu;
    if (col_ind_ilu != NULL) delete [] col_ind_ilu;
  }
  
  if (block != NULL)              delete [] block;
  if (block_weight != NULL)       delete [] block_weight;
  if (block_inverse != NULL)      delete [] block_inverse;
  
  if (prod_block_vector != NULL)  delete [] prod_block_vector;
  if (prod_row_vector != NULL)    delete [] prod_row_vector;
  if (aux_vector != NULL)         delete [] aux_vector;
  if (sum_vector != NULL)         delete [] sum_vector;
  if (invM != NULL)               delete [] invM;
  if (LineletBool != NULL)        delete [] LineletBool;
  if (LineletPoint != NULL)       delete [] LineletPoint;
  
  for (iElem = 0; iElem < max_nElem; iElem++) {
    if (UBlock[iElem] != NULL)      delete [] UBlock[iElem];
    if (invUBlock[iElem] != NULL)   delete [] invUBlock[iElem];
    if (LBlock[iElem] != NULL)      delete [] LBlock[iElem];
    if (yVector[iElem] != NULL)     delete [] yVector[iElem];
    if (zVector[iElem] != NULL)     delete [] zVector[iElem];
    if (rVector[iElem] != NULL)     delete [] rVector[iElem];
  }
  if (UBlock != NULL)     delete [] UBlock;
  if (invUBlock != NULL)  delete [] invUBlock;
  if (LBlock != NULL)     delete [] LBlock;
  if (yVector != NULL)    delete [] yVector;
  if (zVector != NULL)    delete [] zVector;
  if (rVector != NULL)    delete [] rVector;

  if (LFBlock != NULL)    delete [] LFBlock;
  if (LyVector != NULL)   delete [] LyVector;
  if (FzVector != NULL)   delete [] FzVector;
  
}

template <class CalcType>
void TCSysMatrix<CalcType>::Initialize(unsigned long nPoint, unsigned long nPointDomain,
                            unsigned short nVar, unsigned short nEqn,
                            CSparsityPattern *pattern) {

  unsigned long iVar;
  
  sparsity_pattern = pattern;
  
  nnz = sparsity_pattern->GetnNonZero();
  
  this->nPoint       = nPoint;        // Assign number of points in the mesh
  this->nPointDomain = nPointDomain;  // Assign number of points in the mesh
  this->nVar         = nVar;          // Assign number of vars in each block system
  this->nEqn         = nEqn;           // Assign number of eqns in each block system
  
  matrix            = new CalcType [nnz*nVar*nEqn];  // Reserve memory for the values of the matrix
  block             = new CalcType [nVar*nEqn];
  block_weight      = new CalcType [nVar*nEqn];
  block_inverse     = new CalcType [nVar*nEqn];

  prod_block_vector = new CalcType [nEqn];
  prod_row_vector   = new CalcType [nVar];
  aux_vector        = new CalcType [nVar];
  sum_vector        = new CalcType [nVar];
  
  /*--- Memory initialization ---*/
  
  for (iVar = 0; iVar < nnz*nVar*nEqn; iVar++) matrix[iVar]        = 0.0;
  for (iVar = 0; iVar < nVar*nEqn; iVar++)     block[iVar]         = 0.0;
  for (iVar = 0; iVar < nVar*nEqn; iVar++)     block_weight[iVar]  = 0.0;
  for (iVar = 0; iVar < nVar*nEqn; iVar++)     block_inverse[iVar] = 0.0;

  for (iVar = 0; iVar < nEqn; iVar++)          prod_block_vector[iVar] = 0.0;
  for (iVar = 0; iVar < nVar; iVar++)          prod_row_vector[iVar]   = 0.0;
  for (iVar = 0; iVar < nVar; iVar++)          aux_vector[iVar]        = 0.0;
  for (iVar = 0; iVar < nVar; iVar++)          sum_vector[iVar]        = 0.0;

}


template<class CalcType>
CalcType *TCSysMatrix<CalcType>::GetBlock(unsigned long block_i, unsigned long block_j) {
  return &matrix[sparsity_pattern->GetIndex(block_i, block_j)*nVar*nEqn];
}

template<class CalcType>
CalcType TCSysMatrix<CalcType>::GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar, unsigned short jVar) {
  return matrix[sparsity_pattern->GetIndex(block_i, block_j)*nVar*nEqn+iVar*nEqn+jVar];
}

template<class CalcType>
void TCSysMatrix<CalcType>::SetBlock(unsigned long block_i, unsigned long block_j, CalcType **val_block) {
  
  unsigned short iVar, jVar;
  
  unsigned long index = sparsity_pattern->GetIndex(block_i, block_j);
  
  for (iVar = 0; iVar < nVar; iVar++){
    for (jVar = 0; jVar < nEqn; jVar++){
      matrix[index*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar][jVar];
    }
  }
}

template<class CalcType>  
void TCSysMatrix<CalcType>::SetBlock(unsigned long block_i, unsigned long block_j, CalcType *val_block) {
  
  unsigned long iVar, jVar;
  
  unsigned long index = sparsity_pattern->GetIndex(block_i, block_j);
  
  for (iVar = 0; iVar < nVar; iVar++){
    for (jVar = 0; jVar < nEqn; jVar++){
      matrix[index*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar*nVar+jVar];
    }
  } 
}

template<class CalcType>
void TCSysMatrix<CalcType>::AddBlock(unsigned long block_i, unsigned long block_j, CalcType **val_block) {
  
  unsigned long iVar, jVar;
  
  unsigned long index = sparsity_pattern->GetIndex(block_i, block_j);
  
  for (iVar = 0; iVar < nVar; iVar++){
    for (jVar = 0; jVar < nEqn; jVar++){
      matrix[index*nVar*nEqn+iVar*nEqn+jVar] += val_block[iVar][jVar];
    }
  }
}

template<class CalcType>
void TCSysMatrix<CalcType>::SubtractBlock(unsigned long block_i, unsigned long block_j, CalcType **val_block) {
  
  unsigned long iVar, jVar;
  
  unsigned long index = sparsity_pattern->GetIndex(block_i, block_j);

  for (iVar = 0; iVar < nVar; iVar++){
    for (jVar = 0; jVar < nEqn; jVar++){
      matrix[index*nVar*nEqn+iVar*nEqn+jVar] -= val_block[iVar][jVar];
    }
  }
}

template<class CalcType>
CalcType *TCSysMatrix<CalcType>::GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j) {
  
  unsigned long step = 0, index;
  
  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    step++;
    if (col_ind_ilu[index] == block_j) { return &(ILU_matrix[(row_ptr_ilu[block_i]+step-1)*nVar*nEqn]); }
  }
  return NULL;
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, CalcType *val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    step++;
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[(row_ptr_ilu[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar*nVar+jVar];
      break;
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::SetBlockTransposed_ILUMatrix(unsigned long block_i, unsigned long block_j, CalcType *val_block) {

  unsigned long iVar, jVar, index, step = 0;

  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    step++;
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[(row_ptr_ilu[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = val_block[jVar*nVar+iVar];
      break;
    }
  }

}

template<class CalcType>
void TCSysMatrix<CalcType>::SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, CalcType *val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr_ilu[block_i]; index < row_ptr_ilu[block_i+1]; index++) {
    step++;
    if (col_ind_ilu[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[(row_ptr_ilu[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] -= val_block[iVar*nVar+jVar];
      break;
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixVectorProduct(CalcType *matrix, CalcType *vector, CalcType *product) {
  
  unsigned short iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    product[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      product[iVar] += matrix[iVar*nVar+jVar] * vector[jVar];
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixMatrixProduct(CalcType *matrix_a, CalcType *matrix_b, CalcType *product) {
  
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

template<class CalcType>
void TCSysMatrix<CalcType>::AddVal2Diag(unsigned long block_i, CalcType val_matrix) {
  
  unsigned long iVar;
  
  unsigned long index = sparsity_pattern->GetIndex(block_i, block_i);
  
  for (iVar = 0; iVar < nVar; iVar++){
    matrix[index*nVar*nVar+iVar*nVar+iVar] += val_matrix;
  }
}

template<class CalcType>
void TCSysMatrix<CalcType>::SetVal2Diag(unsigned long block_i, CalcType val_matrix) {
  
  unsigned long iVar, jVar;
  
  unsigned long index = sparsity_pattern->GetIndex(block_i, block_i);  
  
  for (iVar = 0; iVar < nVar; iVar++){
    for (jVar = 0; jVar < nVar; jVar++){
      matrix[index*nVar*nVar+iVar*nVar+jVar] = 0.0;
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar++){
    matrix[index*nVar*nVar+iVar*nVar+iVar] = val_matrix;
  }
}

template<class CalcType>
void TCSysMatrix<CalcType>::DeleteValsRowi(unsigned long i) {
  
  unsigned long block_i = i/nVar;
  unsigned long row = i - block_i*nVar;
  unsigned long index, iVar;
  
  for (index = sparsity_pattern->GetRowPointer(block_i); index < sparsity_pattern->GetRowPointer(block_i+1); index++) {
    for (iVar = 0; iVar < nVar; iVar++)
      matrix[index*nVar*nVar+row*nVar+iVar] = 0.0; // Delete row values in the block
    if (sparsity_pattern->GetColumnIndex(index) == block_i)
      matrix[index*nVar*nVar+row*nVar+row] = 1.0; // Set 1 to the diagonal element
  }
  
}


template<class CalcType>
CalcType TCSysMatrix<CalcType>::MatrixDeterminant(CalcType **a, unsigned long n) {
  
  unsigned long i, j, j1, j2;
  CalcType det = 0;
  CalcType **m = NULL;
  
  if (n < 1) { }
  else if (n == 1) { det = a[0][0]; }
  else if (n == 2) { det = a[0][0] * a[1][1] - a[1][0] * a[0][1]; }
  else {
    det = 0.0;

    for (j1=0;j1<n;j1++) {
      m = new CalcType*[n-1];
      for (i=0;i<n-1;i++)
        m[i] = new CalcType[n-1];
      
      for (i=1;i<n;i++) {
        j2 = 0;
        for (j=0;j<n;j++) {
          if (j == j1)
          continue;
          m[i-1][j2] = a[i][j];
          j2++;
        }
      }
      
      det += pow(-1.0,j1+2.0) * a[0][j1] * MatrixDeterminant(m,n-1);
      for (i=0;i<n-1;i++)
      delete [] m[i];
      delete [] m;
    }
    
  }
  
  return(det);
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixCoFactor(CalcType **a, unsigned long n, CalcType **b) {
  
  unsigned long i,j,ii,jj,i1,j1;
  CalcType det;
  CalcType **c;
  
  c = new CalcType*[n-1];
  for (i=0;i<n-1;i++)
    c[i] = new CalcType[n-1];
  
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      
      /*--- Form the adjoint a_ij ---*/
      i1 = 0;
      for (ii=0;ii<n;ii++) {
        if (ii == i)
        continue;
        j1 = 0;
        for (jj=0;jj<n;jj++) {
          if (jj == j)
          continue;
          c[i1][j1] = a[ii][jj];
          j1++;
        }
        i1++;
      }
      
      /*--- Calculate the determinate ---*/
      det = MatrixDeterminant(c,n-1);
      
      /*--- Fill in the elements of the cofactor ---*/
      b[i][j] = pow(-1.0,i+j+2.0) * det;
    }
  }
  for (i=0;i<n-1;i++)
    delete [] c[i];
  delete [] c;
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixTranspose(CalcType **a, unsigned long n) {
  
  unsigned long i, j;
  CalcType tmp;
  
  for (i=1;i<n;i++) {
    for (j=0;j<i;j++) {
      tmp = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = tmp;
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::Gauss_Elimination(unsigned long block_i, CalcType* rhs, bool transposed) {
  
  short iVar, jVar, kVar; // This is important, otherwise some compilers optimizations will fail
  CalcType weight, aux;
  
  CalcType *Block = GetBlock(block_i, block_i);
  
  /*--- Copy block matrix, note that the original matrix
   is modified by the algorithm---*/
  
  if (!transposed) {
    for (iVar = 0; iVar < (short)nVar; iVar++)
      for (jVar = 0; jVar < (short)nVar; jVar++)
        block[iVar*nVar+jVar] = Block[iVar*nVar+jVar];
  } else {
    for (iVar = 0; iVar < (short)nVar; iVar++)
      for (jVar = 0; jVar < (short)nVar; jVar++)
        block[iVar*nVar+jVar] = Block[jVar*nVar+iVar];
  }
  /*--- Gauss elimination ---*/
  
  if (nVar == 1) {
    rhs[0] /= block[0];
  }
  else {
    
    /*--- Transform system in Upper Matrix ---*/
    
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = block[iVar*nVar+jVar] / block[jVar*nVar+jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          block[iVar*nVar+kVar] -= weight*block[jVar*nVar+kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
    }
    
    /*--- Backwards substitution ---*/
    
    rhs[nVar-1] = rhs[nVar-1] / block[nVar*nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0.0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += block[iVar*nVar+jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux) / block[iVar*nVar+iVar];
      if (iVar == 0) break;
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::Gauss_Elimination_ILUMatrix(unsigned long block_i, CalcType* rhs) {
  
  short iVar, jVar, kVar; // This is important, otherwise some compilers optimizations will fail
  CalcType weight, aux;
  
  CalcType *Block = GetBlock_ILUMatrix(block_i, block_i);
  
  /*--- Copy block matrix, note that the original matrix
   is modified by the algorithm---*/
  
  for (iVar = 0; iVar < (short)nVar; iVar++)
    for (jVar = 0; jVar < (short)nVar; jVar++)
      block[iVar*nVar+jVar] = Block[iVar*nVar+jVar];
  
  /*--- Gauss elimination ---*/

  if (nVar == 1) {
    rhs[0] /= block[0];
  }
  else {
    
    /*--- Transform system in Upper Matrix ---*/
    
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = block[iVar*nVar+jVar] / block[jVar*nVar+jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          block[iVar*nVar+kVar] -= weight*block[jVar*nVar+kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
    }
    
    /*--- Backwards substitution ---*/
    
    rhs[nVar-1] = rhs[nVar-1] / block[nVar*nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0.0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += block[iVar*nVar+jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux) / block[iVar*nVar+iVar];
      if (iVar == 0) break;
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::Gauss_Elimination(CalcType* Block, CalcType* rhs) {
  
  short iVar, jVar, kVar; // This is important, otherwise some compilers optimizations will fail
  CalcType weight, aux;
  
  /*--- Copy block matrix, note that the original matrix
   is modified by the algorithm---*/
  
  for (iVar = 0; iVar < (short)nVar; iVar++)
    for (jVar = 0; jVar < (short)nVar; jVar++)
      block[iVar*nVar+jVar] = Block[iVar*nVar+jVar];
  
  
  if (nVar == 1) {
    rhs[0] /= block[0];
  }
  else {
    /*--- Transform system in Upper Matrix ---*/
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = block[iVar*nVar+jVar] / block[jVar*nVar+jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          block[iVar*nVar+kVar] -= weight*block[jVar*nVar+kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
    }
    
    /*--- Backwards substitution ---*/
    rhs[nVar-1] = rhs[nVar-1] / block[nVar*nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0.0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += block[iVar*nVar+jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux) / block[iVar*nVar+iVar];
      if (iVar == 0) break;
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::ProdBlockVector(unsigned long block_i, unsigned long block_j, const TCSysVector<CalcType> & vec) {
  
  unsigned long j = block_j*nVar;
  unsigned short iVar, jVar;
  
  CalcType *block = GetBlock(block_i, block_j);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    prod_block_vector[iVar] = 0;
    for (jVar = 0; jVar < nVar; jVar++)
      prod_block_vector[iVar] += block[iVar*nVar+jVar]*vec[j+jVar];
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::UpperProduct(TCSysVector<CalcType> & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = sparsity_pattern->GetRowPointer(row_i); index < sparsity_pattern->GetRowPointer(row_i+1); index++) {
    if (sparsity_pattern->GetColumnIndex(index) > row_i) {
      ProdBlockVector(row_i, sparsity_pattern->GetColumnIndex(index), vec);
      for (iVar = 0; iVar < nVar; iVar++)
        prod_row_vector[iVar] += prod_block_vector[iVar];
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::LowerProduct(TCSysVector<CalcType> & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = sparsity_pattern->GetRowPointer(row_i); index < sparsity_pattern->GetRowPointer(row_i+1); index++) {
    if (sparsity_pattern->GetColumnIndex(index) < row_i) {
      ProdBlockVector(row_i, sparsity_pattern->GetColumnIndex(index), vec);
      for (iVar = 0; iVar < nVar; iVar++)
        prod_row_vector[iVar] += prod_block_vector[iVar];
    }
  }

}

template<class CalcType>
void TCSysMatrix<CalcType>::DiagonalProduct(TCSysVector<CalcType> & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = sparsity_pattern->GetRowPointer(row_i); index < sparsity_pattern->GetRowPointer(row_i+1); index++) {
    if (sparsity_pattern->GetColumnIndex(index) == row_i) {
      ProdBlockVector(row_i, sparsity_pattern->GetColumnIndex(index), vec);
      for (iVar = 0; iVar < nVar; iVar++)
        prod_row_vector[iVar] += prod_block_vector[iVar];
    }
  }
  
}


template<class CalcType>
void TCSysMatrix<CalcType>::RowProduct(const TCSysVector<CalcType> & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = sparsity_pattern->GetRowPointer(row_i); index < sparsity_pattern->GetRowPointer(row_i+1); index++) {
    ProdBlockVector(row_i, sparsity_pattern->GetColumnIndex(index), vec);
    for (iVar = 0; iVar < nVar; iVar++)
      prod_row_vector[iVar] += prod_block_vector[iVar];
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixVectorProduct(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod) {
  
  unsigned long iPoint, iVar;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    RowProduct(vec, iPoint);
    for (iVar = 0; iVar < nVar; iVar++)
      prod[iPoint*nVar+iVar] = prod_row_vector[iVar];
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixVectorProduct(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config) {
  
  unsigned long prod_begin, vec_begin, mat_begin, index, iVar, jVar, row_i;
  
  assert (nVar   == vec.GetNVar());
  assert (nVar   == prod.GetNVar());
  assert (nPoint == vec.GetNBlk());
  assert (nPoint == prod.GetNBlk());
  
  prod = CalcType(0.0); // set all entries of prod to zero
  for (row_i = 0; row_i < nPointDomain; row_i++) {
    prod_begin = row_i*nVar; // offset to beginning of block row_i
    for (index = sparsity_pattern->GetRowPointer(row_i); index < sparsity_pattern->GetRowPointer(row_i+1); index++) {
      vec_begin = sparsity_pattern->GetColumnIndex(index)*nVar; // offset to beginning of block col_ind[index]
      mat_begin = (index*nVar*nVar); // offset to beginning of matrix block[row_i][col_ind[indx]]
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          prod[(unsigned long)(prod_begin+iVar)] += matrix[(unsigned long)(mat_begin+iVar*nVar+jVar)]*vec[(unsigned long)(vec_begin+jVar)];
        }
      }
    }
  }
  
  /*--- MPI Parallelization ---*/
  prod.SendReceive(geometry, config);
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::MatrixVectorProductTransposed(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config) {

  unsigned long prod_begin, vec_begin, mat_begin, index, iVar, jVar , row_i;

  assert (nVar   == vec.GetNVar());
  assert (nVar   == prod.GetNVar());
  assert (nPoint == vec.GetNBlk());
  assert (nPoint == prod.GetNBlk());

  prod = CalcType(0.0); // set all entries of prod to zero
  for (row_i = 0; row_i < nPointDomain; row_i++) {
    vec_begin = row_i*nVar; // offset to beginning of block col_ind[index]
    for (index = sparsity_pattern->GetRowPointer(row_i); index < sparsity_pattern->GetRowPointer(row_i+1); index++) {
      prod_begin = sparsity_pattern->GetColumnIndex(index)*nVar; // offset to beginning of block row_i
      mat_begin = (index*nVar*nVar); // offset to beginning of matrix block[row_i][col_ind[indx]]
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
            prod[(unsigned long)(prod_begin+jVar)] += matrix[(unsigned long)(mat_begin+iVar*nVar+jVar)]*vec[(unsigned long)(vec_begin+iVar)];
        }
      }
    }
  }

  /*--- MPI Parallelization ---*/
  
  prod.SendReceive_Reverse(geometry, config);

}

template<class CalcType>
void TCSysMatrix<CalcType>::GetMultBlockBlock(CalcType *c, CalcType *a, CalcType *b) {
  
  unsigned long iVar, jVar, kVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++) {
      c[iVar*nVar+jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        c[iVar*nVar+jVar] += a[iVar*nVar+kVar] * b[kVar*nVar+jVar];
    }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::GetMultBlockVector(CalcType *c, CalcType *a, CalcType *b) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    c[iVar] =  0.0;
    for (jVar = 0; jVar < nVar; jVar++)
      c[iVar] += a[iVar*nVar+jVar] * b[jVar];
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::GetSubsBlock(CalcType *c, CalcType *a, CalcType *b) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      c[iVar*nVar+jVar] = a[iVar*nVar+jVar] - b[iVar*nVar+jVar];
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::GetSubsVector(CalcType *c, CalcType *a, CalcType *b) {
  
  unsigned long iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    c[iVar] = a[iVar] - b[iVar];
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::InverseBlock(CalcType *Block, CalcType *invBlock) {
  
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

template<class CalcType>
void TCSysMatrix<CalcType>::InverseDiagonalBlock(unsigned long block_i, CalcType *invBlock, bool transpose) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++)
      aux_vector[jVar] = 0.0;
    aux_vector[iVar] = 1.0;
    
    /*--- Compute the i-th column of the inverse matrix ---*/
    
    Gauss_Elimination(block_i, aux_vector, transpose);
    for (jVar = 0; jVar < nVar; jVar++)
      invBlock[jVar*nVar+iVar] = aux_vector[jVar];
  }
  
  //  CalcType Det, **Matrix, **CoFactor;
  //  CalcType *Block = GetBlock(block_i, block_i);
  //
  //  Matrix = new CalcType*[nVar];
  //  CoFactor = new CalcType*[nVar];
  //  for (iVar=0;iVar<nVar;iVar++) {
  //    Matrix[iVar] = new CalcType[nVar];
  //    CoFactor[iVar] = new CalcType[nVar];
  //  }
  //
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    for (jVar = 0; jVar < nVar; jVar++)
  //    Matrix[iVar][jVar] = Block[jVar*nVar+iVar];
  //  }
  //
  //  Det =  MatrixDeterminant(Matrix, nVar);
  //  MatrixCoFactor(Matrix, nVar, CoFactor);
  //  MatrixTranspose(CoFactor, nVar);
  //
  //
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    for (jVar = 0; jVar < nVar; jVar++)
  //    invBlock[jVar*nVar+iVar] = CoFactor[iVar][jVar]/Det;
  //  }
  //
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    delete [] Matrix[iVar];
  //    delete [] CoFactor[iVar];
  //  }
  //  delete [] Matrix;
  //  delete [] CoFactor;
  
}


template<class CalcType>
void TCSysMatrix<CalcType>::InverseDiagonalBlock_ILUMatrix(unsigned long block_i, CalcType *invBlock) {
  
  unsigned long iVar, jVar;

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++)
      aux_vector[jVar] = 0.0;
    aux_vector[iVar] = 1.0;
    
    /*--- Compute the i-th column of the inverse matrix ---*/
    
    Gauss_Elimination_ILUMatrix(block_i, aux_vector);
    for (jVar = 0; jVar < nVar; jVar++)
      invBlock[jVar*nVar+iVar] = aux_vector[jVar];
  }
  
  //  CalcType Det, **Matrix, **CoFactor;
  //  CalcType *Block = GetBlock_ILUMatrix(block_i, block_i);
  //
  //  Matrix = new CalcType*[nVar];
  //  CoFactor = new CalcType*[nVar];
  //  for (iVar=0;iVar<nVar;iVar++) {
  //    Matrix[iVar] = new CalcType[nVar];
  //    CoFactor[iVar] = new CalcType[nVar];
  //  }
  //
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    for (jVar = 0; jVar < nVar; jVar++)
  //    Matrix[iVar][jVar] = Block[jVar*nVar+iVar];
  //  }
  //
  //  Det =  MatrixDeterminant(Matrix, nVar);
  //  MatrixCoFactor(Matrix, nVar, CoFactor);
  //  MatrixTranspose(CoFactor, nVar);
  //
  //
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    for (jVar = 0; jVar < nVar; jVar++)
  //    invBlock[jVar*nVar+iVar] = CoFactor[iVar][jVar]/Det;
  //  }
  //
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    delete [] Matrix[iVar];
  //    delete [] CoFactor[iVar];
  //  }
  //  delete [] Matrix;
  //  delete [] CoFactor;
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::BuildJacobiPreconditioner(bool transpose) {

  unsigned long iPoint, iVar, jVar;

  /*--- Compute Jacobi Preconditioner ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Compute the inverse of the diagonal block ---*/
    InverseDiagonalBlock(iPoint, block_inverse, transpose);

    /*--- Set the inverse of the matrix to the invM structure (which is a vector) ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        invM[iPoint*nVar*nVar+iVar*nVar+jVar] = block_inverse[iVar*nVar+jVar];
  }

}


template<class CalcType>
void TCSysMatrix<CalcType>::ComputeJacobiPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, iVar, jVar;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      prod[(unsigned long)(iPoint*nVar+iVar)] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++)
        prod[(unsigned long)(iPoint*nVar+iVar)] +=
        invM[(unsigned long)(iPoint*nVar*nVar+iVar*nVar+jVar)]*vec[(unsigned long)(iPoint*nVar+jVar)];
    }
  }
  
  /*--- MPI Parallelization ---*/
  
  prod.SendReceive(geometry, config);
  
}

template<class CalcType>
unsigned long TCSysMatrix<CalcType>::Jacobi_Smoother(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec, CalcType tol, unsigned long m, CalcType *residual, bool monitoring, CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, iVar, jVar;
  
  /*---  Check the number of iterations requested ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for smoothing iterations, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  /*--- Create vectors to hold the residual and the Matrix-Vector product
   of the Jacobian matrix with the current solution (x^k). These must be
   stored in order to perform multiple iterations of the smoother. ---*/
  
  TCSysVector<CalcType> r(b);
  TCSysVector<CalcType> A_x(b);
  
  /*--- Calculate the initial residual, compute norm, and check
   if system is already solved. Recall, r holds b initially. ---*/
  
  mat_vec(x, A_x);
  r -= A_x;
  CalcType norm_r = r.norm();
  CalcType norm0  = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysMatrix<CalcType>::Jacobi_Smoother(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "\n# " << "Jacobi Smoother" << " residual history" << endl;
    cout << "# Residual tolerance target = " << tol << endl;
    cout << "# Initial residual norm     = " << norm_r << endl;
    cout << "     " << i << "     " << norm_r/norm0 << endl;
  }
  
  /*---  Loop over all smoothing iterations ---*/
  
  for (i = 0; i < (int)m; i++) {
    
    /*--- Apply the Jacobi smoother, i.e., multiply by the inverse of the
     diagonal matrix of A, which was built in the preprocessing phase. Note
     that we are directly updating the solution (x^k+1) during the loop. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++)
          x[(unsigned long)(iPoint*nVar+iVar)] +=
          invM[(unsigned long)(iPoint*nVar*nVar+iVar*nVar+jVar)]*r[(unsigned long)(iPoint*nVar+jVar)];
      }
    }
    
    /*--- MPI Parallelization ---*/
    
    x.SendReceive(geometry, config);
    
    /*--- Update the residual (r^k+1 = b - A*x^k+1) with the new solution ---*/
    
    r = b;
    mat_vec(x, A_x);
    r -= A_x;
    
    /*--- Check if solution has converged, else output the relative
     residual if necessary. ---*/
    
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 5 == 0))
      cout << "     " << i << "     " << norm_r/norm0 << endl;
    
  }
  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# Jacobi smoother final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
  return (unsigned long) i;
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::BuildILUPreconditioner(bool transposed) {
  
  unsigned long index, index_, iVar;
  CalcType *Block_ij, *Block_jk;
  long iPoint, jPoint, kPoint;
  

  /*--- Copy block matrix, note that the original matrix
   is modified by the algorithm, so that we have the factorization stored
   in the ILUMatrix at the end of this preprocessing. ---*/

  for (iVar = 0; iVar < nnz_ilu*nVar*nEqn; iVar++) ILU_matrix[iVar] = 0.0;

  for (iPoint = 0; iPoint < (long)nPointDomain; iPoint++) {
    for (index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; index++) {
      jPoint = col_ind[index];
      if (transposed) {
        Block_ij = GetBlock(jPoint, iPoint);
        SetBlockTransposed_ILUMatrix(iPoint, jPoint, Block_ij);
      } else {
        Block_ij = GetBlock(iPoint, jPoint);
        SetBlock_ILUMatrix(iPoint, jPoint, Block_ij);
      }
    }
  }
  
  /*--- Transform system in Upper Matrix ---*/
  
  for (iPoint = 1; iPoint < (long)nPointDomain; iPoint++) {
    
    /*--- For each row (unknown), loop over all entries in A on this row
     row_ptr_ilu[iPoint+1] will have the index for the first entry on the next
     row. ---*/
    
    for (index = row_ptr_ilu[iPoint]; index < row_ptr_ilu[iPoint+1]; index++) {
      
      /*--- jPoint here is the column for each entry on this row ---*/
      
      jPoint = col_ind_ilu[index];
      
      /*--- Check that this column is in the lower triangular portion ---*/
      
      if ((jPoint < iPoint) && (jPoint < (long)nPointDomain)) {
        
        /*--- If we're in the lower triangle, get the pointer to this block,
         invert it, and then right multiply against the original block ---*/
        
        Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
        InverseDiagonalBlock_ILUMatrix(jPoint, block_inverse);
        MatrixMatrixProduct(Block_ij, block_inverse, block_weight);
        
        /*--- block_weight holds Aij*inv(Ajj). Jump to the row for jPoint ---*/
        
        for (index_ = row_ptr_ilu[jPoint]; index_ < row_ptr_ilu[jPoint+1]; index_++) {
          
          /*--- Get the column of the entry ---*/
          
          kPoint = col_ind_ilu[index_];
          
          /*--- If the column is greater than or equal to jPoint, i.e., the
           upper triangular part, then multiply and modify the matrix.
           Here, Aik' = Aik - Aij*inv(Ajj)*Ajk. ---*/
          
          if ((kPoint >= jPoint) && (jPoint < (long)nPointDomain)) {
            
            Block_jk = GetBlock_ILUMatrix(jPoint, kPoint);
            MatrixMatrixProduct(block_weight, Block_jk, block);
            SubtractBlock_ILUMatrix(iPoint, kPoint, block);
            
          }
        }
        
        /*--- Lastly, store block_weight in the lower triangular part, which
         will be reused during the forward solve in the precon/smoother. ---*/
        
        SetBlock_ILUMatrix(iPoint, jPoint, block_weight);
        
      }
    }
  }
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::ComputeILUPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config) {
  
  unsigned long index;
  CalcType *Block_ij;
  long iPoint, jPoint;
  unsigned short iVar;
  
  /*--- Copy block matrix, note that the original matrix
   is modified by the algorithm---*/
  
  for (iPoint = 0; iPoint < (long)nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      prod[iPoint*nVar+iVar] = vec[iPoint*nVar+iVar];
    }
  }
  
  /*--- Forward solve the system using the lower matrix entries that
   were computed and stored during the ILU preprocessing. Note
   that we are overwriting the residual vector as we go. ---*/
  
  for (iPoint = 1; iPoint < (long)nPointDomain; iPoint++) {
    for (index = row_ptr_ilu[iPoint]; index < row_ptr_ilu[iPoint+1]; index++) {
      jPoint = col_ind_ilu[index];
      if ((jPoint < iPoint) && (jPoint < (long)nPointDomain)) {
        Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
        MatrixVectorProduct(Block_ij, &prod[jPoint*nVar], aux_vector);
        for (iVar = 0; iVar < nVar; iVar++)
          prod[iPoint*nVar+iVar] -= aux_vector[iVar];
        
      }
    }
  }
  
  /*--- Backwards substitution (starts at the last row) ---*/
  
  InverseDiagonalBlock_ILUMatrix((nPointDomain-1), block_inverse);
  MatrixVectorProduct(block_inverse, &prod[(nPointDomain-1)*nVar], aux_vector);
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod[ (nPointDomain-1)*nVar + iVar] = aux_vector[iVar];
  
  for (iPoint = nPointDomain-2; iPoint >= 0; iPoint--) {
    for (iVar = 0; iVar < nVar; iVar++) sum_vector[iVar] = 0.0;
    for (index = row_ptr_ilu[iPoint]; index < row_ptr_ilu[iPoint+1]; index++) {
      jPoint = col_ind_ilu[index];
      if ((jPoint >= iPoint+1) && (jPoint < (long)nPointDomain)) {
        Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
        MatrixVectorProduct(Block_ij, &prod[jPoint*nVar], aux_vector);
        for (iVar = 0; iVar < nVar; iVar++) sum_vector[iVar] += aux_vector[iVar];
      }
    }
    for (iVar = 0; iVar < nVar; iVar++) prod[iPoint*nVar+iVar] = (prod[iPoint*nVar+iVar]-sum_vector[iVar]);
    InverseDiagonalBlock_ILUMatrix(iPoint, block_inverse);
    MatrixVectorProduct(block_inverse, &prod[iPoint*nVar], aux_vector);
    for (iVar = 0; iVar < nVar; iVar++) prod[iPoint*nVar+iVar] = aux_vector[iVar];
    if (iPoint == 0) break;
  }
  
  /*--- MPI Parallelization ---*/
  
  prod.SendReceive(geometry, config);
  
}

template<class CalcType>
unsigned long TCSysMatrix<CalcType>::ILU_Smoother(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec, CalcType tol, unsigned long m, CalcType *residual, bool monitoring, CGeometry *geometry, CConfig *config) {
  
  unsigned long index;
  CalcType *Block_ij, omega = 1.0;
  long iPoint, jPoint;
  unsigned short iVar;
  
  /*---  Check the number of iterations requested ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for smoothing iterations, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  /*--- Create vectors to hold the residual and the Matrix-Vector product
   of the Jacobian matrix with the current solution (x^k). These must be
   stored in order to perform multiple iterations of the smoother. ---*/
  
  TCSysVector<CalcType> r(b);
  TCSysVector<CalcType> A_x(b);
  
  /*--- Calculate the initial residual, compute norm, and check
   if system is already solved. Recall, r holds b initially. ---*/
  
  mat_vec(x, A_x);
  r -= A_x;
  CalcType norm_r = r.norm();
  CalcType norm0  = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysMatrix<CalcType>::ILU_Smoother(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "\n# " << "ILU Smoother" << " residual history" << endl;
    cout << "# Residual tolerance target = " << tol << endl;
    cout << "# Initial residual norm     = " << norm_r << endl;
    cout << "     " << i << "     " << norm_r/norm0 << endl;
  }
  
  /*---  Loop over all smoothing iterations ---*/
  
  for (i = 0; i < (int)m; i++) {
    
    /*--- Forward solve the system using the lower matrix entries that
     were computed and stored during the ILU preprocessing. Note
     that we are overwriting the residual vector as we go. ---*/
    
    for (iPoint = 1; iPoint < (long)nPointDomain; iPoint++) {
      
      /*--- For each row (unknown), loop over all entries in A on this row
       row_ptr_ilu[iPoint+1] will have the index for the first entry on the next
       row. ---*/
      
      for (index = row_ptr_ilu[iPoint]; index < row_ptr_ilu[iPoint+1]; index++) {
        
        /*--- jPoint here is the column for each entry on this row ---*/
        
        jPoint = col_ind_ilu[index];
        
        /*--- Check that this column is in the lower triangular portion ---*/
        
        if ((jPoint < iPoint) && (jPoint < (long)nPointDomain)) {
          
          /*--- Lastly, get Aij*inv(Ajj) from the lower triangular part, which
           was calculated in the preprocessing, and apply to r. ---*/
          
          Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
          MatrixVectorProduct(Block_ij, &r[jPoint*nVar], aux_vector);
          for (iVar = 0; iVar < nVar; iVar++)
            r[iPoint*nVar+iVar] -= aux_vector[iVar];
          
        }
      }
    }
    
    /*--- Backwards substitution (starts at the last row) ---*/
    
    InverseDiagonalBlock_ILUMatrix((nPointDomain-1), block_inverse);
    MatrixVectorProduct(block_inverse, &r[(nPointDomain-1)*nVar], aux_vector);
    
    for (iVar = 0; iVar < nVar; iVar++)
      r[(nPointDomain-1)*nVar + iVar] = aux_vector[iVar];
    
    for (iPoint = nPointDomain-2; iPoint >= 0; iPoint--) {
      for (iVar = 0; iVar < nVar; iVar++) sum_vector[iVar] = 0.0;
      for (index = row_ptr_ilu[iPoint]; index < row_ptr_ilu[iPoint+1]; index++) {
        jPoint = col_ind_ilu[index];
        if ((jPoint >= iPoint+1) && (jPoint < (long)nPointDomain)) {
          Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
          MatrixVectorProduct(Block_ij, &r[jPoint*nVar], aux_vector);
          for (iVar = 0; iVar < nVar; iVar++) sum_vector[iVar] += aux_vector[iVar];
        }
      }
      for (iVar = 0; iVar < nVar; iVar++) r[iPoint*nVar+iVar] = (r[iPoint*nVar+iVar]-sum_vector[iVar]);
      InverseDiagonalBlock_ILUMatrix(iPoint, block_inverse);
      MatrixVectorProduct(block_inverse, &r[iPoint*nVar], aux_vector);
      for (iVar = 0; iVar < nVar; iVar++) r[iPoint*nVar+iVar] = aux_vector[iVar];
      if (iPoint == 0) break;
    }
    
    /*--- Update solution (x^k+1 = x^k + w*M^-1*r^k) using the residual vector,
     which holds the update after applying the ILU smoother, i.e., M^-1*r^k.
     Omega is a relaxation factor that we have currently set to 1.0. ---*/
    
    x.Plus_AX(omega, r);
    
    /*--- MPI Parallelization ---*/
    
    x.SendReceive(geometry, config);
    
    /*--- Update the residual (r^k+1 = b - A*x^k+1) with the new solution ---*/
    
    r = b;
    mat_vec(x, A_x);
    r -= A_x;
    
    /*--- Check if solution has converged, else output the relative 
     residual if necessary. ---*/
    
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 5 == 0))
      cout << "     " << i << "     " << norm_r/norm0 << endl;
    
  }
  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# ILU smoother final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
  return (unsigned int) i;
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::ComputeLU_SGSPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, iVar;
  
  /*--- First part of the symmetric iteration: (D+L).x* = b ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    LowerProduct(prod, iPoint);                                        // Compute L.x*
    for (iVar = 0; iVar < nVar; iVar++)
      aux_vector[iVar] = vec[iPoint*nVar+iVar] - prod_row_vector[iVar]; // Compute aux_vector = b - L.x*
    Gauss_Elimination(iPoint, aux_vector);                            // Solve D.x* = aux_vector
    for (iVar = 0; iVar < nVar; iVar++)
      prod[iPoint*nVar+iVar] = aux_vector[iVar];                       // Assesing x* = solution
  }
  
  /*--- MPI Parallelization ---*/
  
  prod.SendReceive(geometry, config);
  
  /*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/
  
  for (iPoint = nPointDomain-1; (int)iPoint >= 0; iPoint--) {
    DiagonalProduct(prod, iPoint);                 // Compute D.x*
    for (iVar = 0; iVar < nVar; iVar++)
      aux_vector[iVar] = prod_row_vector[iVar];   // Compute aux_vector = D.x*
    UpperProduct(prod, iPoint);                    // Compute U.x_(n+1)
    for (iVar = 0; iVar < nVar; iVar++)
      aux_vector[iVar] -= prod_row_vector[iVar];  // Compute aux_vector = D.x*-U.x_(n+1)
    Gauss_Elimination(iPoint, aux_vector);        // Solve D.x* = aux_vector
    for (iVar = 0; iVar < nVar; iVar++)
      prod[iPoint*nVar + iVar] = aux_vector[iVar]; // Assesing x_(1) = solution
  }
  
  /*--- MPI Parallelization ---*/
  
  prod.SendReceive(geometry, config);
  
}

template<class CalcType>
unsigned long TCSysMatrix<CalcType>::LU_SGS_Smoother(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec, CalcType tol, unsigned long m, CalcType *residual, bool monitoring, CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, iVar;
  CalcType omega = 1.0;
  
  /*---  Check the number of iterations requested ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for smoothing iterations, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  /*--- Create vectors to hold the residual and the Matrix-Vector product
   of the Jacobian matrix with the current solution (x^k). These must be
   stored in order to perform multiple iterations of the smoother. ---*/
  
  TCSysVector<CalcType> r(b);
  TCSysVector<CalcType> A_x(b);
  TCSysVector<CalcType> xStar(x);
  
  /*--- Calculate the initial residual, compute norm, and check
   if system is already solved. Recall, r holds b initially. ---*/
  
  mat_vec(x, A_x);
  r -= A_x;
  CalcType norm_r = r.norm();
  CalcType norm0  = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysMatrix<CalcType>::LU_SGS_Smoother(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "\n# " << "LU_SGS Smoother" << " residual history" << endl;
    cout << "# Residual tolerance target = " << tol << endl;
    cout << "# Initial residual norm     = " << norm_r << endl;
    cout << "     " << i << "     " << norm_r/norm0 << endl;
  }
  
  /*---  Loop over all smoothing iterations ---*/
  
  for (i = 0; i < (int)m; i++) {

    /*--- First part of the symmetric iteration: (D+L).x* = b ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      LowerProduct(xStar, iPoint);                                      // Compute L.x*
      for (iVar = 0; iVar < nVar; iVar++)
        aux_vector[iVar] = r[iPoint*nVar+iVar] - prod_row_vector[iVar]; // Compute aux_vector = b - L.x*
      Gauss_Elimination(iPoint, aux_vector);                            // Solve D.x* = aux_vector
      for (iVar = 0; iVar < nVar; iVar++)
        xStar[iPoint*nVar+iVar] = aux_vector[iVar];                     // Assesing x* = solution, stored in r
    }
    
    /*--- MPI Parallelization ---*/
    
    xStar.SendReceive(geometry, config);
    
    /*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/
    
    for (iPoint = nPointDomain-1; (int)iPoint >= 0; iPoint--) {
      DiagonalProduct(xStar, iPoint);               // Compute D.x*
      for (iVar = 0; iVar < nVar; iVar++)
        aux_vector[iVar] = prod_row_vector[iVar];   // Compute aux_vector = D.x*
      UpperProduct(xStar, iPoint);                  // Compute U.x_(n+1)
      for (iVar = 0; iVar < nVar; iVar++)
        aux_vector[iVar] -= prod_row_vector[iVar];  // Compute aux_vector = D.x*-U.x_(n+1)
      Gauss_Elimination(iPoint, aux_vector);        // Solve D.x* = aux_vector
      for (iVar = 0; iVar < nVar; iVar++)
        xStar[iPoint*nVar+iVar] = aux_vector[iVar]; // Assesing x_(1) = solution
    }
    
    /*--- Update solution (x^k+1 = x^k + w*M^-1*r^k) using the xStar vector,
     which holds the update after applying the LU_SGS smoother, i.e., M^-1*r^k.
     Omega is a relaxation factor that we have currently set to 1.0. ---*/
    
    x.Plus_AX(omega, xStar);
    
    /*--- MPI Parallelization ---*/
    
    x.SendReceive(geometry, config);
    
    /*--- Update the residual (r^k+1 = b - A*x^k+1) with the new solution ---*/
    
    r = b;
    mat_vec(x, A_x);
    r -= A_x;
    xStar = x;
    
    /*--- Check if solution has converged, else output the relative
     residual if necessary. ---*/
    
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 5 == 0))
      cout << "     " << i << "     " << norm_r/norm0 << endl;
    
  }
  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# LU_SGS smoother final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
  return (unsigned int) i;
  
}

template<class CalcType>
void TCSysMatrix<CalcType>::ComputeResidual(const TCSysVector<CalcType> & sol, const TCSysVector<CalcType> & f, TCSysVector<CalcType> & res) {
  
  unsigned long iPoint, iVar;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    RowProduct(sol, iPoint);
    for (iVar = 0; iVar < nVar; iVar++) {
      res[iPoint*nVar+iVar] = prod_row_vector[iVar] - f[iPoint*nVar+iVar];
    }
  }
  
}


template class TCSysMatrix<su2double>;
#ifdef CODI_REVERSE_TYPE
template class TCSysMatrix<passivedouble>;
#endif