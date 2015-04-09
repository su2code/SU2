/*!
 * \file matrix_structure.cpp
 * \brief Main subroutines for doing the sparse structures
 * \author F. Palacios, A. Bueno
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

CSysMatrix::CSysMatrix(void) {
  
  /*--- Array initialization ---*/

  matrix            = NULL;
  row_ptr           = NULL;
  col_ind           = NULL;
  block             = NULL;
  prod_block_vector = NULL;
  prod_row_vector   = NULL;
  aux_vector        = NULL;
  sum_vector        = NULL;
  invM              = NULL;
  
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
  
}

CSysMatrix::~CSysMatrix(void) {
  
  unsigned long iElem;
  
  /*--- Memory deallocation ---*/
  
  if (matrix != NULL)             delete [] matrix;
  if (row_ptr != NULL)            delete [] row_ptr;
  if (col_ind != NULL)            delete [] col_ind;
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

void CSysMatrix::Initialize(unsigned long nPoint, unsigned long nPointDomain,
                            unsigned short nVar, unsigned short nEqn,
                            bool EdgeConnect, CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, *row_ptr, *col_ind, index, nnz, Elem;
  unsigned short iNeigh, iElem, iNode, *nNeigh;
  vector<unsigned long>::iterator it;
  vector<unsigned long> vneighs;
  
  /*--- Don't delete *row_ptr, *col_ind because they are
   asigned to the Jacobian structure. ---*/
  
  /*--- Compute the number of neighbors ---*/
  
  nNeigh = new unsigned short [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    if (EdgeConnect) {
      nNeigh[iPoint] = (geometry->node[iPoint]->GetnPoint()+1);  // +1 -> to include diagonal element
    }
    else {
      vneighs.clear();
      for (iElem = 0; iElem < geometry->node[iPoint]->GetnElem(); iElem++) {
        Elem =  geometry->node[iPoint]->GetElem(iElem);
        for (iNode = 0; iNode < geometry->elem[Elem]->GetnNodes(); iNode++)
          vneighs.push_back(geometry->elem[Elem]->GetNode(iNode));
      }
      vneighs.push_back(iPoint);
      
      sort(vneighs.begin(), vneighs.end());
      it = unique(vneighs.begin(), vneighs.end());
      vneighs.resize(it - vneighs.begin());
      nNeigh[iPoint] = vneighs.size();
    }
    
  }
  
  /*--- Create row_ptr structure, using the number of neighbors ---*/
  
  row_ptr = new unsigned long [nPoint+1];
  row_ptr[0] = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    row_ptr[iPoint+1] = row_ptr[iPoint] + nNeigh[iPoint];
  nnz = row_ptr[nPoint];
  
  /*--- Create col_ind structure ---*/
  
  col_ind = new unsigned long [nnz];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    vneighs.clear();
    
    if (EdgeConnect) {
      for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++)
        vneighs.push_back(geometry->node[iPoint]->GetPoint(iNeigh));
      vneighs.push_back(iPoint);
    }
    else {
      for (iElem = 0; iElem < geometry->node[iPoint]->GetnElem(); iElem++) {
        Elem =  geometry->node[iPoint]->GetElem(iElem);
        for (iNode = 0; iNode < geometry->elem[Elem]->GetnNodes(); iNode++)
          vneighs.push_back(geometry->elem[Elem]->GetNode(iNode));
      }
      vneighs.push_back(iPoint);
    }
    
    sort(vneighs.begin(), vneighs.end());
    it = unique(vneighs.begin(), vneighs.end());
    vneighs.resize( it - vneighs.begin() );
    
    index = row_ptr[iPoint];
    for (iNeigh = 0; iNeigh < vneighs.size(); iNeigh++) {
      col_ind[index] = vneighs[iNeigh];
      index++;
    }
    
  }
  
  /*--- Set the indices in the in the sparce matrix structure, and memory allocation ---*/
  
  SetIndexes(nPoint, nPointDomain, nVar, nEqn, row_ptr, col_ind, nnz, config);
  
  /*--- Initialization matrix to zero ---*/
  
  SetValZero();
  
  delete [] nNeigh;
  
}

void CSysMatrix::SetIndexes(unsigned long val_nPoint, unsigned long val_nPointDomain, unsigned short val_nVar, unsigned short val_nEq, unsigned long* val_row_ptr, unsigned long* val_col_ind, unsigned long val_nnz, CConfig *config) {
  
  unsigned long iVar;
  
  nPoint = val_nPoint;              // Assign number of points in the mesh
  nPointDomain = val_nPointDomain;  // Assign number of points in the mesh
  nVar = val_nVar;                  // Assign number of vars in each block system
  nEqn = val_nEq;                   // Assign number of eqns in each block system
  nnz = val_nnz;                    // Assign number of possible non zero blocks
  row_ptr = val_row_ptr;
  col_ind = val_col_ind;
  
  matrix            = new su2double [nnz*nVar*nEqn];	// Reserve memory for the values of the matrix
  block             = new su2double [nVar*nEqn];
  block_weight      = new su2double [nVar*nEqn];
  block_inverse     = new su2double [nVar*nEqn];

  prod_block_vector = new su2double [nEqn];
  prod_row_vector   = new su2double [nVar];
  aux_vector        = new su2double [nVar];
  sum_vector        = new su2double [nVar];
  
  /*--- Memory initialization ---*/
  
  for (iVar = 0; iVar < nnz*nVar*nEqn; iVar++)    matrix[iVar] = 0.0;
  for (iVar = 0; iVar < nVar*nEqn; iVar++)        block[iVar] = 0.0;
  for (iVar = 0; iVar < nVar*nEqn; iVar++)        block_weight[iVar] = 0.0;
  for (iVar = 0; iVar < nVar*nEqn; iVar++)        block_inverse[iVar] = 0.0;

  for (iVar = 0; iVar < nEqn; iVar++)             prod_block_vector[iVar] = 0.0;
  for (iVar = 0; iVar < nVar; iVar++)             prod_row_vector[iVar] = 0.0;
  for (iVar = 0; iVar < nVar; iVar++)             aux_vector[iVar] = 0.0;
  for (iVar = 0; iVar < nVar; iVar++)             sum_vector[iVar] = 0.0;
  
  
  /*--- Set specific preconditioner matrices (ILU) ---*/
  
  if ((config->GetKind_Linear_Solver_Prec() == ILU) ||
    (config->GetKind_Linear_Solver() == SMOOTHER_ILU)) {
    ILU_matrix = new su2double [nnz*nVar*nEqn];	// Reserve memory for the ILU matrix
    for (iVar = 0; iVar < nnz*nVar*nEqn; iVar++)    ILU_matrix[iVar] = 0.0;
  }
  
  /*--- Set specific preconditioner matrices (Jacobi and Linelet) ---*/
  
  if ((config->GetKind_Linear_Solver_Prec() == JACOBI) ||
      (config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_JACOBI) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET))   {
    invM = new su2double [nPoint*nVar*nEqn];	// Reserve memory for the values of the inverse of the preconditioner
    for (iVar = 0; iVar < nPoint*nVar*nEqn; iVar++) invM[iVar] = 0.0;
  }

}

su2double *CSysMatrix::GetBlock(unsigned long block_i, unsigned long block_j) {
  
  unsigned long step = 0, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) { return &(matrix[(row_ptr[block_i]+step-1)*nVar*nEqn]); }
  }
  return NULL;
  
}

su2double CSysMatrix::GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar, unsigned short jVar) {
  
  unsigned long step = 0, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) { return matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar]; }
  }
  return 0;
  
}

void CSysMatrix::SetBlock(unsigned long block_i, unsigned long block_j, su2double **val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar][jVar];
      break;
    }
  }
  
}
  
void CSysMatrix::SetBlock(unsigned long block_i, unsigned long block_j, su2double *val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar*nVar+jVar];
      break;
    }
  }
  
}

void CSysMatrix::AddBlock(unsigned long block_i, unsigned long block_j, su2double **val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] += val_block[iVar][jVar];
      break;
    }
  }
  
}

void CSysMatrix::SubtractBlock(unsigned long block_i, unsigned long block_j, su2double **val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] -= val_block[iVar][jVar];
      break;
    }
  }
  
}

su2double *CSysMatrix::GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j) {
  
  unsigned long step = 0, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) { return &(ILU_matrix[(row_ptr[block_i]+step-1)*nVar*nEqn]); }
  }
  return NULL;
  
}

void CSysMatrix::SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, su2double *val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] = val_block[iVar*nVar+jVar];
      break;
    }
  }
  
}

void CSysMatrix::SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, su2double *val_block) {
  
  unsigned long iVar, jVar, index, step = 0;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_j) {
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++)
          ILU_matrix[(row_ptr[block_i]+step-1)*nVar*nEqn+iVar*nEqn+jVar] -= val_block[iVar*nVar+jVar];
      break;
    }
  }
  
}

void CSysMatrix::MatrixVectorProduct(su2double *matrix, su2double *vector, su2double *product) {
  
  unsigned short iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    product[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      product[iVar] += matrix[iVar*nVar+jVar] * vector[jVar];
    }
  }
  
}

void CSysMatrix::MatrixMatrixProduct(su2double *matrix_a, su2double *matrix_b, su2double *product) {
  
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

void CSysMatrix::AddVal2Diag(unsigned long block_i, su2double val_matrix) {
  
  unsigned long step = 0, iVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_i) {	// Only elements on the diagonal
      for (iVar = 0; iVar < nVar; iVar++)
        matrix[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+iVar] += val_matrix;
      break;
    }
  }
  
}

void CSysMatrix::SetVal2Diag(unsigned long block_i, su2double val_matrix) {
  
  unsigned long step = 0, iVar, jVar, index;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    step++;
    if (col_ind[index] == block_i) {	// Only elements on the diagonal
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          matrix[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+jVar] = 0.0;
      
      for (iVar = 0; iVar < nVar; iVar++)
        matrix[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+iVar] = val_matrix;
      
      break;
    }
  }
  
}

void CSysMatrix::DeleteValsRowi(unsigned long i) {
  
  unsigned long block_i = i/nVar;
  unsigned long row = i - block_i*nVar;
  unsigned long index, iVar;
  
  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    for (iVar = 0; iVar < nVar; iVar++)
      matrix[index*nVar*nVar+row*nVar+iVar] = 0.0; // Delete row values in the block
    if (col_ind[index] == block_i)
      matrix[index*nVar*nVar+row*nVar+row] = 1.0; // Set 1 to the diagonal element
  }
  
}

void CSysMatrix::Gauss_Elimination(unsigned long block_i, su2double* rhs) {
  
  short iVar, jVar, kVar; // This is important, otherwise some compilers optimizations will fail
  su2double weight, aux;
  
  su2double *Block = GetBlock(block_i, block_i);
  
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
        for (kVar = jVar; kVar < nVar; kVar++)
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

void CSysMatrix::Gauss_Elimination_ILUMatrix(unsigned long block_i, su2double* rhs) {
  
  short iVar, jVar, kVar; // This is important, otherwise some compilers optimizations will fail
  su2double weight, aux;
  
  su2double *Block = GetBlock_ILUMatrix(block_i, block_i);
  
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
        for (kVar = jVar; kVar < nVar; kVar++)
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

void CSysMatrix::Gauss_Elimination(su2double* Block, su2double* rhs) {
  
  short iVar, jVar, kVar; // This is important, otherwise some compilers optimizations will fail
  su2double weight, aux;
  
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
        for (kVar = jVar; kVar < nVar; kVar++)
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

void CSysMatrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, const CSysVector & vec) {
  
  unsigned long j = block_j*nVar;
  unsigned short iVar, jVar;
  
  su2double *block = GetBlock(block_i, block_j);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    prod_block_vector[iVar] = 0;
    for (jVar = 0; jVar < nVar; jVar++)
      prod_block_vector[iVar] += block[iVar*nVar+jVar]*vec[j+jVar];
  }
  
}

void CSysMatrix::UpperProduct(CSysVector & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] > row_i) {
      ProdBlockVector(row_i, col_ind[index], vec);
      for (iVar = 0; iVar < nVar; iVar++)
        prod_row_vector[iVar] += prod_block_vector[iVar];
    }
  }
  
}

void CSysMatrix::LowerProduct(CSysVector & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] < row_i) {
      ProdBlockVector(row_i, col_ind[index], vec);
      for (iVar = 0; iVar < nVar; iVar++)
        prod_row_vector[iVar] += prod_block_vector[iVar];
    }
  }

}

void CSysMatrix::DiagonalProduct(CSysVector & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] == row_i) {
      ProdBlockVector(row_i, col_ind[index], vec);
      for (iVar = 0; iVar < nVar; iVar++)
        prod_row_vector[iVar] += prod_block_vector[iVar];
    }
  }
  
}

void CSysMatrix::SendReceive_Solution(CSysVector & x, CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive = NULL, *Buffer_Send = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive = new su2double [nBufferR_Vector];
      Buffer_Send = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send[iVertex*nVar+iVar] = x[iPoint*nVar+iVar];
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      
      SU2_MPI::Sendrecv(Buffer_Send, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive[iVar*nVertexR+iVertex] = Buffer_Send[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
          x[iPoint*nVar+iVar] = Buffer_Receive[iVertex*nVar+iVar];
        
      }
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive;
      
    }
    
  }
  
}

void CSysMatrix::RowProduct(const CSysVector & vec, unsigned long row_i) {
  
  unsigned long iVar, index;
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod_row_vector[iVar] = 0;
  
  for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    ProdBlockVector(row_i, col_ind[index], vec);
    for (iVar = 0; iVar < nVar; iVar++)
      prod_row_vector[iVar] += prod_block_vector[iVar];
  }
  
}

void CSysMatrix::MatrixVectorProduct(const CSysVector & vec, CSysVector & prod) {
  
  unsigned long iPoint, iVar;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    RowProduct(vec, iPoint);
    for (iVar = 0; iVar < nVar; iVar++)
      prod[iPoint*nVar+iVar] = prod_row_vector[iVar];
  }
  
}

void CSysMatrix::MatrixVectorProduct(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config) {
  
  unsigned long prod_begin, vec_begin, mat_begin, index, iVar, jVar, row_i;
  
  /*--- Some checks for consistency between CSysMatrix and the CSysVectors ---*/
  if ( (nVar != vec.GetNVar()) || (nVar != prod.GetNVar()) ) {
    cerr << "CSysMatrix::MatrixVectorProduct(const CSysVector&, CSysVector): "
    << "nVar values incompatible." << endl;
    throw(-1);
  }
  if ( (nPoint != vec.GetNBlk()) || (nPoint != prod.GetNBlk()) ) {
    cerr << "CSysMatrix::MatrixVectorProduct(const CSysVector&, CSysVector): "
    << "nPoint and nBlk values incompatible." << endl;
    throw(-1);
  }
  
  prod = su2double(0.0); // set all entries of prod to zero
  for (row_i = 0; row_i < nPointDomain; row_i++) {
    prod_begin = row_i*nVar; // offset to beginning of block row_i
    for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
      vec_begin = col_ind[index]*nVar; // offset to beginning of block col_ind[index]
      mat_begin = (index*nVar*nVar); // offset to beginning of matrix block[row_i][col_ind[indx]]
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          prod[(unsigned long)(prod_begin+iVar)] += matrix[(unsigned long)(mat_begin+iVar*nVar+jVar)]*vec[(unsigned long)(vec_begin+jVar)];
        }
      }
    }
  }
  
  /*--- MPI Parallelization ---*/
  SendReceive_Solution(prod, geometry, config);
  
}

void CSysMatrix::GetMultBlockBlock(su2double *c, su2double *a, su2double *b) {
  
  unsigned long iVar, jVar, kVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++) {
      c[iVar*nVar+jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        c[iVar*nVar+jVar] += a[iVar*nVar+kVar] * b[kVar*nVar+jVar];
    }
  
}

void CSysMatrix::GetMultBlockVector(su2double *c, su2double *a, su2double *b) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    c[iVar] =  0.0;
    for (jVar = 0; jVar < nVar; jVar++)
      c[iVar] += a[iVar*nVar+jVar] * b[jVar];
  }
  
}

void CSysMatrix::GetSubsBlock(su2double *c, su2double *a, su2double *b) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      c[iVar*nVar+jVar] = a[iVar*nVar+jVar] - b[iVar*nVar+jVar];
  
}

void CSysMatrix::GetSubsVector(su2double *c, su2double *a, su2double *b) {
  
  unsigned long iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    c[iVar] = a[iVar] - b[iVar];
  
}

void CSysMatrix::InverseBlock(su2double *Block, su2double *invBlock) {
  
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

void CSysMatrix::InverseDiagonalBlock(unsigned long block_i, su2double *invBlock) {
  
  unsigned long iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++)
      aux_vector[jVar] = 0.0;
    aux_vector[iVar] = 1.0;
    
    /*--- Compute the i-th column of the inverse matrix ---*/
    
    Gauss_Elimination(block_i, aux_vector);
    for (jVar = 0; jVar < nVar; jVar++)
      invBlock[jVar*nVar+iVar] = aux_vector[jVar];
  }
  
}

void CSysMatrix::InverseDiagonalBlock_ILUMatrix(unsigned long block_i, su2double *invBlock) {
  
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
  
}

void CSysMatrix::BuildJacobiPreconditioner(void) {
  
  unsigned long iPoint, iVar, jVar;
  
  /*--- Compute Jacobi Preconditioner ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Compute the inverse of the diagonal block ---*/
    InverseDiagonalBlock(iPoint, block_inverse);
    
    /*--- Set the inverse of the matrix to the invM structure (which is a vector) ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        invM[iPoint*nVar*nVar+iVar*nVar+jVar] = block_inverse[iVar*nVar+jVar];
  }
  
}

void CSysMatrix::BuildILUPreconditioner(void) {

/*--- Reimplement is such a way the LU is a preprocessing ---*/

}

unsigned short CSysMatrix::BuildLineletPreconditioner(CGeometry *geometry, CConfig *config) {
  
  bool *check_Point, add_point;
  unsigned long iEdge, iPoint, jPoint, index_Point, iLinelet, iVertex, next_Point, counter, iElem;
  unsigned short iMarker, iNode, ExtraLines = 100, MeanPoints;
  su2double alpha = 0.9, weight, max_weight, *normal, area, volume_iPoint, volume_jPoint;
  unsigned long Local_nPoints, Local_nLineLets, Global_nPoints, Global_nLineLets;

  /*--- Memory allocation --*/

  check_Point = new bool [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    check_Point[iPoint] = true;
  
  LineletBool = new bool[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
    LineletBool[iPoint] = false;
  
  nLinelet = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC    ) ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC   ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC) ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == DISPLACEMENT_BOUNDARY)) {
      nLinelet += geometry->nVertex[iMarker];
    }
  }
  
  /*--- If the domain contains well defined Linelets ---*/
  
  if (nLinelet != 0) {
    
    /*--- Basic initial allocation ---*/
    
    LineletPoint = new vector<unsigned long>[nLinelet + ExtraLines];
    
    /*--- Define the basic linelets, starting from each vertex ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
          (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC    ) ||
          (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC ) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC   ) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC) ||
          (config->GetMarker_All_KindBC(iMarker) == EULER_WALL             ) ||
          (config->GetMarker_All_KindBC(iMarker) == DISPLACEMENT_BOUNDARY)) {
        iLinelet = 0;
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          LineletPoint[iLinelet].push_back(iPoint);
          check_Point[iPoint] = false;
          iLinelet++;
        }
      }
    }
    
    /*--- Create the linelet structure ---*/
    
    iLinelet = 0;
    
    do {
      
      index_Point = 0;
      
      do {
        
        /*--- Compute the value of the max weight ---*/
        
        iPoint = LineletPoint[iLinelet][index_Point];
        max_weight = 0.0;
        for (iNode = 0; iNode < geometry->node[iPoint]->GetnPoint(); iNode++) {
          jPoint = geometry->node[iPoint]->GetPoint(iNode);
          if ((check_Point[jPoint]) && geometry->node[jPoint]->GetDomain()) {
            iEdge = geometry->FindEdge(iPoint, jPoint);
            normal = geometry->edge[iEdge]->GetNormal();
            if (geometry->GetnDim() == 3) area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
            else area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
            volume_iPoint = geometry->node[iPoint]->GetVolume();
            volume_jPoint = geometry->node[jPoint]->GetVolume();
            weight = 0.5*area*((1.0/volume_iPoint)+(1.0/volume_jPoint));
            max_weight = max(max_weight, weight);
          }
        }
        
        /*--- Verify if any face of the control volume must be added ---*/
        
        add_point = false;
        counter = 0;
        next_Point = geometry->node[iPoint]->GetPoint(0);
        for (iNode = 0; iNode < geometry->node[iPoint]->GetnPoint(); iNode++) {
          jPoint = geometry->node[iPoint]->GetPoint(iNode);
          iEdge = geometry->FindEdge(iPoint, jPoint);
          normal = geometry->edge[iEdge]->GetNormal();
          if (geometry->GetnDim() == 3) area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
          else area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
          volume_iPoint = geometry->node[iPoint]->GetVolume();
          volume_jPoint = geometry->node[jPoint]->GetVolume();
          weight = 0.5*area*((1.0/volume_iPoint)+(1.0/volume_jPoint));
          if (((check_Point[jPoint]) && (weight/max_weight > alpha) && (geometry->node[jPoint]->GetDomain())) &&
              ((index_Point == 0) || ((index_Point > 0) && (jPoint != LineletPoint[iLinelet][index_Point-1])))) {
            add_point = true;
            next_Point = jPoint;
            counter++;
          }
        }
        
        /*--- We have arrived to an isotropic zone ---*/
        
        if (counter > 1) add_point = false;
        
        /*--- Add a typical point to the linelet, no leading edge ---*/
        
        if (add_point) {
          LineletPoint[iLinelet].push_back(next_Point);
          check_Point[next_Point] = false;
          index_Point++;
        }
        
      } while (add_point);
      iLinelet++;
    } while (iLinelet < nLinelet);
    
    /*--- Identify the points that belong to a Linelet ---*/
    
    for (iLinelet = 0; iLinelet < nLinelet; iLinelet++) {
      for (iElem = 0; iElem < LineletPoint[iLinelet].size(); iElem++) {
        iPoint = LineletPoint[iLinelet][iElem];
        LineletBool[iPoint] = true;
      }
    }
    
    /*--- Identify the maximum number of elements in a Linelet ---*/
    
    max_nElem = LineletPoint[0].size();
    for (iLinelet = 1; iLinelet < nLinelet; iLinelet++)
      if (LineletPoint[iLinelet].size() > max_nElem)
        max_nElem = LineletPoint[iLinelet].size();
    
  }
  
  /*--- The domain doesn't have well defined linelets ---*/
  
  else {
    
    max_nElem = 0;
  
  }
  
  /*--- Screen output ---*/
  
  Local_nPoints = 0;
  for (iLinelet = 0; iLinelet < nLinelet; iLinelet++) {
    Local_nPoints += LineletPoint[iLinelet].size();
  }
  Local_nLineLets = nLinelet;

#ifndef HAVE_MPI
  Global_nPoints = Local_nPoints;
  Global_nLineLets = Local_nLineLets;
#else
  SU2_MPI::Allreduce(&Local_nPoints, &Global_nPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nLineLets, &Global_nLineLets, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  MeanPoints = int(su2double(Global_nPoints)/su2double(Global_nLineLets));
  
  /*--- Memory allocation --*/

  UBlock = new su2double* [max_nElem];
  invUBlock = new su2double* [max_nElem];
  LBlock = new su2double* [max_nElem];
  yVector = new su2double* [max_nElem];
  zVector = new su2double* [max_nElem];
  rVector = new su2double* [max_nElem];
  for (iElem = 0; iElem < max_nElem; iElem++) {
    UBlock[iElem] = new su2double [nVar*nVar];
    invUBlock[iElem] = new su2double [nVar*nVar];
    LBlock[iElem] = new su2double [nVar*nVar];
    yVector[iElem] = new su2double [nVar];
    zVector[iElem] = new su2double [nVar];
    rVector[iElem] = new su2double [nVar];
  }
  
  LFBlock = new su2double [nVar*nVar];
  LyVector = new su2double [nVar];
  FzVector = new su2double [nVar];
  
  /*--- Memory deallocation --*/

  delete [] check_Point;

  return MeanPoints;
  
}

void CSysMatrix::ComputeJacobiPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config) {
  
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
  
  SendReceive_Solution(prod, geometry, config);
  
}

void CSysMatrix::ComputeILUPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config) {
  
  unsigned long index, index_;
  su2double *Block_ij, *Block_jk;
  long iPoint, jPoint, kPoint;
  unsigned short iVar;
  
  /*--- Copy block matrix, note that the original matrix
   is modified by the algorithm---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; index++) {
      jPoint = col_ind[index];
      Block_ij = GetBlock(iPoint, jPoint);
      SetBlock_ILUMatrix(iPoint, jPoint, Block_ij);
    }
    for (iVar = 0; iVar < nVar; iVar++) {
      prod[iPoint*nVar+iVar] = vec[iPoint*nVar+iVar];
    }
  }
  
  /*--- Transform system in Upper Matrix ---*/
  
  for (iPoint = 1; iPoint < nPointDomain; iPoint++) {
    for (index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; index++) {
      jPoint = col_ind[index];
      if ((jPoint < iPoint) && (jPoint < nPointDomain)) {
        Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
        InverseDiagonalBlock_ILUMatrix(jPoint, block_inverse);
        MatrixMatrixProduct(Block_ij, block_inverse, block_weight);
        for (index_ = row_ptr[jPoint]; index_ < row_ptr[jPoint+1]; index_++) {
          kPoint = col_ind[index_];
          if (kPoint < nPointDomain) {
            Block_jk = GetBlock_ILUMatrix(jPoint, kPoint);
            if (kPoint >= jPoint) {
              MatrixMatrixProduct(Block_jk, block_weight, block);
              SubtractBlock_ILUMatrix(iPoint, kPoint, block);
            }
          }
        }
        MatrixVectorProduct(block_weight, &prod[jPoint*nVar], aux_vector);
        for (iVar = 0; iVar < nVar; iVar++)
          prod[iPoint*nVar+iVar] -= aux_vector[iVar];
        
      }
    }
  }
  
  /*--- Backwards substitution ---*/
  
  InverseDiagonalBlock_ILUMatrix((nPointDomain-1), block_inverse);
  MatrixVectorProduct(block_inverse, &prod[(nPointDomain-1)*nVar], aux_vector);
  
  for (iVar = 0; iVar < nVar; iVar++)
    prod[ (nPointDomain-1)*nVar + iVar] = aux_vector[iVar];
  
  for (iPoint = nPointDomain-2; iPoint >= 0; iPoint--) {
    for (iVar = 0; iVar < nVar; iVar++) sum_vector[iVar] = 0.0;
    for (index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; index++) {
      jPoint = col_ind[index];
      if (jPoint < nPointDomain) {
        Block_ij = GetBlock_ILUMatrix(iPoint, jPoint);
        if ((jPoint >= iPoint+1) && (jPoint < nPointDomain)) {
          MatrixVectorProduct(Block_ij, &prod[jPoint*nVar], aux_vector);
          for (iVar = 0; iVar < nVar; iVar++) sum_vector[iVar] += aux_vector[iVar];
        }
      }
    }
    for (iVar = 0; iVar < nVar; iVar++) prod[iPoint*nVar+iVar] = (prod[iPoint*nVar+iVar]-sum_vector[iVar]);
    InverseDiagonalBlock_ILUMatrix(iPoint, block_inverse);
    MatrixVectorProduct(block_inverse, &prod[iPoint*nVar], aux_vector);
    for (iVar = 0; iVar < nVar; iVar++) prod[iPoint*nVar+iVar] = aux_vector[iVar];
    if (iPoint == 0) break;
  }
  
  /*--- MPI Parallelization ---*/
  
  SendReceive_Solution(prod, geometry, config);
  
}

void CSysMatrix::ComputeLU_SGSPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config) {
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
  
  SendReceive_Solution(prod, geometry, config);
  
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
  
  SendReceive_Solution(prod, geometry, config);
  
}

void CSysMatrix::ComputeLineletPreconditioner(const CSysVector & vec, CSysVector & prod,
                                              CGeometry *geometry, CConfig *config) {
  
  unsigned long iVar, jVar, nElem = 0, iLinelet, im1Point, iPoint, ip1Point, iElem;
  long iElemLoop;
  su2double *block;
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  if (size == SINGLE_NODE) {
    
    /*--- Jacobi preconditioning if there is no linelet ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!LineletBool[iPoint]) {
        for (iVar = 0; iVar < nVar; iVar++) {
          prod[(unsigned long)(iPoint*nVar+iVar)] = 0.0;
          for (jVar = 0; jVar < nVar; jVar++)
            prod[(unsigned long)(iPoint*nVar+iVar)] +=
            invM[(unsigned long)(iPoint*nVar*nVar+iVar*nVar+jVar)]*vec[(unsigned long)(iPoint*nVar+jVar)];
        }
      }
    }
    
    /*--- MPI Parallelization ---*/
    
    SendReceive_Solution(prod, geometry, config);
    
    /*--- Solve linelet using a Thomas' algorithm ---*/
    
    for (iLinelet = 0; iLinelet < nLinelet; iLinelet++) {
      
      nElem = LineletPoint[iLinelet].size();
      
      /*--- Copy vec vector to the new structure ---*/
      
      for (iElem = 0; iElem < nElem; iElem++) {
        iPoint = LineletPoint[iLinelet][iElem];
        for (iVar = 0; iVar < nVar; iVar++)
          rVector[iElem][iVar] = vec[(unsigned long)(iPoint*nVar+iVar)];
      }
      
      /*--- Initialization (iElem = 0) ---*/
      
      iPoint = LineletPoint[iLinelet][0];
      block = GetBlock(iPoint, iPoint);
      for (iVar = 0; iVar < nVar; iVar++) {
        yVector[0][iVar] = rVector[0][iVar];
        for (jVar = 0; jVar < nVar; jVar++)
          UBlock[0][iVar*nVar+jVar] = block[iVar*nVar+jVar];
      }
      
      /*--- Main loop (without iElem = 0) ---*/
      
      for (iElem = 1; iElem < nElem; iElem++) {
        
        im1Point = LineletPoint[iLinelet][iElem-1];
        iPoint = LineletPoint[iLinelet][iElem];
        
        InverseBlock(UBlock[iElem-1], invUBlock[iElem-1]);
        block = GetBlock(iPoint, im1Point); GetMultBlockBlock(LBlock[iElem], block, invUBlock[iElem-1]);
        block = GetBlock(im1Point, iPoint); GetMultBlockBlock(LFBlock, LBlock[iElem], block);
        block = GetBlock(iPoint, iPoint); GetSubsBlock(UBlock[iElem], block, LFBlock);
        
        /*--- Forward substituton ---*/
        
        GetMultBlockVector(LyVector, LBlock[iElem], yVector[iElem-1]);
        GetSubsVector(yVector[iElem], rVector[iElem], LyVector);
        
      }
      
      /*--- Backward substituton ---*/
      
      InverseBlock(UBlock[nElem-1], invUBlock[nElem-1]);
      GetMultBlockVector(zVector[nElem-1], invUBlock[nElem-1], yVector[nElem-1]);
      
      for (iElemLoop = nElem-2; iElemLoop >= 0; iElemLoop--) {
        iPoint = LineletPoint[iLinelet][iElemLoop];
        ip1Point = LineletPoint[iLinelet][iElemLoop+1];
        block = GetBlock(iPoint, ip1Point); GetMultBlockVector(FzVector, block, zVector[iElemLoop+1]);
        GetSubsVector(aux_vector, yVector[iElemLoop], FzVector);
        GetMultBlockVector(zVector[iElemLoop], invUBlock[iElemLoop], aux_vector);
      }
      
      /*--- Copy zVector to the prod vector ---*/
      
      for (iElem = 0; iElem < nElem; iElem++) {
        iPoint = LineletPoint[iLinelet][iElem];
        for (iVar = 0; iVar < nVar; iVar++)
          prod[(unsigned long)(iPoint*nVar+iVar)] = zVector[iElem][iVar];
      }
      
    }
    
    /*--- MPI Parallelization ---*/
    
    SendReceive_Solution(prod, geometry, config);
    
  }
  else {
    
    if (rank == MASTER_NODE) cout << "ERROR: Linelet not implemented in parallel." << endl;
    
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
    
  }
  
}

void CSysMatrix::ComputeResidual(const CSysVector & sol, const CSysVector & f, CSysVector & res) {
  
  unsigned long iPoint, iVar;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    RowProduct(sol, iPoint);
    for (iVar = 0; iVar < nVar; iVar++) {
      res[iPoint*nVar+iVar] = prod_row_vector[iVar] - f[iPoint*nVar+iVar];
    }
  }
  
}
