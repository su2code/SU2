/*!
 * \file sparse_structure.cpp
 * \brief Main subroutines for doing the sparse structures.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

#include "../include/sparse_structure.hpp"

CSparseMatrix::CSparseMatrix(void) { }

CSparseMatrix::~CSparseMatrix(void) {
	delete [] val;
	delete [] row_ptr;
	delete [] col_ind;
	delete [] block;
	delete [] prod_block_vector;
	delete [] prod_row_vector;
	delete [] aux_vector;

	if (invM != NULL)
		delete [] invM;
}

void CSparseMatrix::SetIndexes(unsigned long val_nPoint, unsigned short val_nVar, unsigned long* val_row_ptr, 
		unsigned long* val_col_ind, unsigned long val_nnz, bool preconditioner) {

	nPoint = val_nPoint;	// Assign number of points in the mesh
	nVar = val_nVar;		// Assign number of vars in each block system
	nnz = val_nnz;			// Assign number of possible non zero blocks
	row_ptr = val_row_ptr;
	col_ind = val_col_ind;
	val = new double [nnz*nVar*nVar];	// Reserve memory for the values of the matrix
	block = new double [nVar*nVar]; 
	prod_block_vector = new double [nVar];
	prod_row_vector = new double [nVar];
	aux_vector = new double [nVar];

	if (preconditioner)
		invM = new double [nPoint*nVar*nVar];	// Reserve memory for the values of the inverse of the preconditioner
}

void CSparseMatrix::GetBlock(unsigned long block_i, unsigned long block_j) {
	unsigned long step = 0;
	for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
		step++;
		if (col_ind[index] == block_j) {
			for (unsigned long iVar = 0; iVar < nVar*nVar; iVar++)
				block[iVar] = val[(row_ptr[block_i]+step-1)*nVar*nVar+iVar];
			break;
		}
	}
}

void CSparseMatrix::DisplayBlock(void) {
	unsigned short iVar, jVar;
	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++)
			cout << block[iVar*nVar+jVar] << "  ";
		cout << endl;
	}
}

void CSparseMatrix::AddBlock(unsigned long block_i, unsigned long block_j, double **val_block) {
	unsigned long iVar, jVar, index, step = 0;
	for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
		step++;
		if (col_ind[index] == block_j) {
			for (iVar = 0; iVar < nVar; iVar++)
				for (jVar = 0; jVar < nVar; jVar++)
					val[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+jVar] += val_block[iVar][jVar];
			break;
		}
	}
}

void CSparseMatrix::SubtractBlock(unsigned long block_i, unsigned long block_j, double **val_block) {
	unsigned long iVar, jVar, index, step = 0;
	for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
		step++;
		if (col_ind[index] == block_j) {
			for (iVar = 0; iVar < nVar; iVar++)
				for (jVar = 0; jVar < nVar; jVar++)
					val[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+jVar] -= val_block[iVar][jVar];
			break;
		}
	}
}

void CSparseMatrix::AddVal2Diag(unsigned long block_i, double val_val) {
	unsigned long step = 0;
	unsigned long iVar;
	for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
		step++;
		if (col_ind[index] == block_i) {	// Only elements on the diagonal
			for (iVar = 0; iVar < nVar; iVar++)
				val[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+iVar] += val_val;
			break;
		}
	}
}

void CSparseMatrix::DeleteValsRowi(unsigned long i) {
	unsigned long block_i = i/nVar;
	unsigned long row = i - block_i*nVar;
	unsigned long index, iVar;

	for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
		for (iVar = 0; iVar < nVar; iVar++)
			val[index*nVar*nVar+row*nVar+iVar] = 0.0; // Delete row values in the block
		if (col_ind[index] == block_i)
			val[index*nVar*nVar+row*nVar+row] = 1.0; // Set 1 to the diagonal element
	}
}

double CSparseMatrix::SumAbsRowi(unsigned long i) {
	unsigned long block_i = i/nVar;
	unsigned long row = i - block_i*nVar;

	double sum = 0;
	for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
		for (unsigned long iVar = 0; iVar < nVar; iVar ++)
			sum += fabs(val[index*nVar*nVar+row*nVar+iVar]);

	return sum;
}

void CSparseMatrix::Gauss_Elimination(unsigned long block_i, double* rhs) {
	unsigned long iVar, jVar, kVar;
	double weight;

	GetBlock(block_i,block_i);

	if (nVar == 1)
		rhs[0] /= (block[0]+EPS*EPS);
	else {
		/*--- Transform system in Upper Matrix ---*/
		for (iVar = 1; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < iVar; jVar++) {
				weight = block[iVar*nVar+jVar]/(block[jVar*nVar+jVar]+EPS*EPS);
				for (kVar = jVar; kVar < nVar; kVar++)
					block[iVar*nVar+kVar] -= weight*block[jVar*nVar+kVar];
				rhs[iVar] -= weight*rhs[jVar];
			}
		}
		/*--- Backwards substitution ---*/
		double aux;
		rhs[nVar-1] = rhs[nVar-1]/(block[nVar*nVar-1]+EPS*EPS);
		for (iVar = nVar-2; iVar >= 0; iVar--) {
			aux = 0;
			for (jVar = iVar+1; jVar < nVar; jVar++)
				aux += block[iVar*nVar+jVar]*rhs[jVar];
			rhs[iVar] = (rhs[iVar]-aux)/(block[iVar*nVar+iVar]+EPS*EPS);
			if (iVar == 0) break;
		}
	}
}

void CSparseMatrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, double* vec) {
	unsigned long j = block_j*nVar;
	unsigned short iVar, jVar;

	GetBlock(block_i,block_j);

	for (iVar = 0; iVar < nVar; iVar++) {
		prod_block_vector[iVar] = 0;
		for (jVar = 0; jVar < nVar; jVar++)
			prod_block_vector[iVar] += block[iVar*nVar+jVar]*vec[j+jVar];
	}
}

void CSparseMatrix::UpperProduct(double* vec, unsigned long row_i) {
	unsigned long iVar, index;
	
	for (iVar = 0; iVar < nVar; iVar++)
		prod_row_vector[iVar] = 0;

	for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
		if (col_ind[index] > row_i) {
			ProdBlockVector(row_i,col_ind[index],vec);
			for (iVar = 0; iVar < nVar; iVar++)
				prod_row_vector[iVar] += prod_block_vector[iVar];
		}
	}
}

void CSparseMatrix::LowerProduct(double* vec, unsigned long row_i) {
	unsigned long iVar, index;
	
	for (iVar = 0; iVar < nVar; iVar++)
		prod_row_vector[iVar] = 0;

	for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
		if (col_ind[index] < row_i) {
			ProdBlockVector(row_i,col_ind[index],vec);
			for (iVar = 0; iVar < nVar; iVar++)
				prod_row_vector[iVar] += prod_block_vector[iVar];
		}
	}
}

void CSparseMatrix::DiagonalProduct(double* vec, unsigned long row_i) {
	unsigned long iVar, index;
	
	for (iVar = 0; iVar < nVar; iVar++)
		prod_row_vector[iVar] = 0;

	for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
		if (col_ind[index] == row_i) {
			ProdBlockVector(row_i,col_ind[index],vec);
			for (iVar = 0; iVar < nVar; iVar++)
				prod_row_vector[iVar] += prod_block_vector[iVar];
		}
	}
}

void CSparseMatrix::LU_SGSIteration(double* b, double* x_n) {
	unsigned long iPoint, iVar;

	/*--- First part of the symmetric iteration: (D+L).x* = b ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		LowerProduct(x_n,iPoint);  // compute L.x*
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] = b[iPoint*nVar+iVar] - prod_row_vector[iVar]; // compute aux_vector = b-L.x*
		Gauss_Elimination(iPoint, aux_vector); // solve D.x* = aux_vector
		for (iVar = 0; iVar < nVar; iVar++)
			x_n[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x* = solution
	}

	/*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/
	for (iPoint = nPoint-1; iPoint >= 0; iPoint--) {
		DiagonalProduct(x_n,iPoint);  // compute D.x*
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] = prod_row_vector[iVar]; // compute aux_vector = D.x*
		UpperProduct(x_n,iPoint);  // compute U.x_(n+1)
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] -= prod_row_vector[iVar]; // compute aux_vector = D.x*-U·x_(n+1)
		Gauss_Elimination(iPoint, aux_vector); // solve D.x* = aux_vector
		for (iVar = 0; iVar < nVar; iVar++)
			x_n[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x_(1) = solution
		if (iPoint == 0) break;
	}
}

double CSparseMatrix::SGSIteration(double* b, double* x_n) {
	unsigned long iPoint, iVar;
	double norm_j = 0.0, norm_2 = 0.0;

	/*--- First part of the symmetric iteration: (D+L).x* = b-U.x_n ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		UpperProduct(x_n,iPoint);  // compute U.x_n (result written in prod_row_vector)
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] = b[iPoint*nVar+iVar] - prod_row_vector[iVar]; // compute aux_vector = b-U.x_n
		LowerProduct(x_n,iPoint);  // compute L.x* (idem)
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] -= prod_row_vector[iVar]; // compute aux_vector = b-U.x_n-L.x*
		Gauss_Elimination(iPoint, aux_vector); // solve D.x* = aux_vector
		for (iVar = 0; iVar < nVar; iVar++) {
			norm_j += (x_n[iPoint*nVar+iVar]-aux_vector[iVar])*(x_n[iPoint*nVar+iVar]-aux_vector[iVar]);
			x_n[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x* = solution
		}
	}
	norm_j = sqrt(norm_j); 

	/*--- Second part of the symmetric iteration: (D+U)·x_(n+1) = b-L·x* ---*/
	for (iPoint = nPoint-1; iPoint >= 0; iPoint--) {
		LowerProduct(x_n,iPoint);  // compute L.x* (result written in prod_row_vector)
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] = b[iPoint*nVar+iVar] - prod_row_vector[iVar]; // compute aux_vector = b-L.x*
		UpperProduct(x_n,iPoint);  // compute U.x_(n+1) (idem)
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] -= prod_row_vector[iVar]; // compute aux_vector = b-L.x*-U.x_(n+1)
		Gauss_Elimination(iPoint, aux_vector); // solve D.x_(n+1) = aux_vector
		for (iVar = 0; iVar < nVar; iVar++) {
			norm_2 += (x_n[iPoint*nVar+iVar]-aux_vector[iVar])*(x_n[iPoint*nVar+iVar]-aux_vector[iVar]);
			x_n[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x* = solution
		}
		if (iPoint == 0) break;
	}
	norm_2 = sqrt(norm_2);

	return (norm_j+norm_2);
}

void CSparseMatrix::SGSSolution(double* b, double* x_i, double tol, int max_it, bool monitoring) {
	int it = 0;
	double norm = tol + 1;
	
	while (norm > tol) {
		it++;
		norm = SGSIteration(b,x_i);
		
		if ((monitoring == true) && (it%100 == 0))
			cout << "SGS-Solution:: Iteration = " << it << " ; Norm of the residual: " << norm << endl;
		if (it == max_it)
			return;
	}
}

void CSparseMatrix::RowProduct(double* vec, unsigned long row_i) {
	unsigned long iVar, index;
	
	for (iVar = 0; iVar < nVar; iVar++)
		prod_row_vector[iVar] = 0;

	for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
		ProdBlockVector(row_i,col_ind[index],vec);
		for (iVar = 0; iVar < nVar; iVar++)
			prod_row_vector[iVar] += prod_block_vector[iVar];
	}
}

void CSparseMatrix::MatrixVectorProduct(double* vec, double* prod) {
	unsigned long iPoint, iVar;

	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		RowProduct(vec,iPoint);
		for (iVar = 0; iVar < nVar; iVar++)
			prod[iPoint*nVar+iVar] = prod_row_vector[iVar];
	}
}

void CSparseMatrix::MatrixVectorProduct(const CSysVector & vec, CSysVector & prod) {

  unsigned long prod_begin, vec_begin, mat_begin, index, iVar, jVar, row_i;
  // some checks for consistency between CSparseMatrix and the CSysVectors
  if ( (nVar != vec.GetNVar()) || (nVar != prod.GetNVar()) ) {
    cerr << "CSparseMatrix::MatrixVectorProduct(const CSysVector&, CSysVector): " 
	 << "nVar values incompatible." << endl;
    throw(-1);
  }
  if ( (nPoint != vec.GetNBlk()) || (nPoint != prod.GetNBlk()) ) {
    cerr << "CSparseMatrix::MatrixVectorProduct(const CSysVector&, CSysVector): " 
	 << "nPoint and nBlk values incompatible." << endl;
    throw(-1);
  }

  prod = 0.0; // set all entries of prod to zero
  for (row_i = 0; row_i < nPoint; row_i++) {
    prod_begin = row_i*nVar; // offset to beginning of block row_i
    for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
      vec_begin = col_ind[index]*nVar; // offset to beginning of block col_ind[index]
      mat_begin = (index*nVar*nVar); // offset to beginning of matrix block[row_i][col_ind[indx]]
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          prod[prod_begin+iVar] += val[mat_begin+iVar*nVar+jVar]*vec[vec_begin+jVar];
        }
      }
    }
  }
}

void CSparseMatrix::CGSolution(double* b, double* x_i, double tol, int max_it, bool monitoring) {
	unsigned long iPoint, iVar, total_index;
	double *d, *r, *Atimesd;
	d = new double [nPoint * nVar];
	r = new double [nPoint * nVar];
	Atimesd = new double [nPoint * nVar];

	/*--- CG-algorithm ---*/
	unsigned short it = 0;
	double norm; // = tol + 1;
	double alpha, alpha_den, beta, norm_r, norm_r_new;

	/*--- First iteration: determine residual (r) and search direction (d) 
          d(0) = r(0) = b-A.x(0) ---*/
	MatrixVectorProduct(x_i, r);
        norm = 0.0;
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			r[total_index] = b[total_index] - r[total_index];
			d[total_index] = r[total_index];
                        norm += r[total_index]*r[total_index];                        
		}
	}

	while (norm > tol) {
		it++;

		/*--- Compute descend step (alpha) ---*/
		MatrixVectorProduct(d, Atimesd);
		alpha_den = 0.0; norm_r = 0;
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				norm_r += r[total_index]*r[total_index];
				alpha_den += d[total_index]*Atimesd[total_index];
			}
		}

		alpha = norm_r / alpha_den;

		/*--- Update solution and residual ---*/
		norm_r_new = 0.0; // we need the norm of the updated residual
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				x_i[total_index] += alpha*d[total_index];
				r[total_index] -= alpha*Atimesd[total_index];
				norm_r_new += r[total_index]*r[total_index];
			}
		}

		/*--- Compute Gram-Schmidt constant (beta) ---*/
		beta = norm_r_new / norm_r;

		/*--- Update search direction (d) ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				d[total_index] = r[total_index] + beta*d[total_index];
			}
		}

		/*--- Monitoring of the result ---*/
		norm = sqrt(norm_r_new);
		if ((monitoring == true) && (it%500 == 0))
			cout << "CG-Solution:: Iteration = " << it << " ; Norm of the residual: " << norm << endl;
		if (it == max_it)
			return;
	}

	delete [] Atimesd;
	delete [] r;
	delete [] d;
}

void CSparseMatrix::InverseDiagonalBlock(unsigned long block_i, double **invBlock) {

	unsigned long iVar, jVar;

	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++)
			aux_vector[jVar] = 0.0;
		aux_vector[iVar] = 1.0;
		
		/*--- Compute the i-th column of the inverse matrix ---*/
		Gauss_Elimination(block_i, aux_vector);
		for (jVar = 0; jVar < nVar; jVar++)
			invBlock[jVar][iVar] = aux_vector[jVar];
	}

}

void CSparseMatrix::BuildJacobiPreconditioner(void) {
	unsigned long iPoint, iVar, jVar;
	double **invBlock;

	/*--- Small nVar x nVar matrix for intermediate computations ---*/
	invBlock = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		invBlock[iVar] = new double [nVar];

	/*--- Compute Jacobi Preconditioner ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		
		/*--- Compute the inverse of the diagonal block ---*/
		InverseDiagonalBlock(iPoint, invBlock);

		/*--- Set the inverse of the matrix to the invM structure (which is a vector) ---*/
		for (iVar = 0; iVar < nVar; iVar++) 
			for (jVar = 0; jVar < nVar; jVar++) 
				invM[iPoint*nVar*nVar+iVar*nVar+jVar] = invBlock[iVar][jVar];
	}

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] invBlock[iVar];
	delete [] invBlock;
}

void CSparseMatrix::PrecondVectorProduct(double* vec, double* prod) {
	unsigned long iPoint, iVar, jVar;

	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			prod[iPoint*nVar+iVar] = 0;
			for (jVar = 0; jVar < nVar; jVar++)
				prod[iPoint*nVar+iVar] += invM[iPoint*nVar*nVar+iVar*nVar+jVar]*vec[iPoint*nVar+jVar];
		}
	}
}

void CSparseMatrix::PrecondVectorProduct(const CSysVector & vec, CSysVector & prod) {
	unsigned long iPoint, iVar, jVar;

	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			prod[iPoint*nVar+iVar] = 0;
			for (jVar = 0; jVar < nVar; jVar++)
				prod[iPoint*nVar+iVar] += invM[iPoint*nVar*nVar+iVar*nVar+jVar]*vec[iPoint*nVar+jVar];
		}
	}
}

void CSparseMatrix::PreconditionedCGSolution(double* b, double* x_i, double tol, int max_it, bool monitoring) {
	unsigned long iPoint, iVar, total_index;
	double *d, *r, *Atimesd, *invMtimesr;
	d = new double [nPoint * nVar];
	r = new double [nPoint * nVar];
	Atimesd = new double [nPoint * nVar];
	invMtimesr = new double [nPoint * nVar];
	double rinvMr;

	BuildJacobiPreconditioner();

	/*--- First iteration: determine residual (r) and search direction (d) ---*/
	/*--- r(0) = b-Ax(0) ---*/
	MatrixVectorProduct(x_i, r); // r = Ax(0)
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			r[total_index] = b[total_index] - r[total_index]; // r(0) = b-A·x(0);
		}
	}
	
	/*--- d(0) = invM*r(0) ---*/
	PrecondVectorProduct(r, invMtimesr);
	rinvMr = 0;
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			d[total_index] = invMtimesr[total_index];
			rinvMr += r[total_index]*invMtimesr[total_index];
		}
	}

	/*--- CG-algorithm ---*/
	unsigned short it = 0;
	double norm = tol + 1;
	double alpha, alpha_den, beta, norm_r_new, beta_den;
	
	while (norm > tol) {

		it++;

		/*--- Compute descend step (alpha) ---*/
		MatrixVectorProduct(d, Atimesd);
		alpha_den = 0.0; 
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				alpha_den += d[total_index]*Atimesd[total_index];
			}
		}
		alpha = rinvMr / (alpha_den + EPS*EPS);
		
		/*--- Update solution and residual ---*/
		norm_r_new = 0.0; // we need the norm of the updated residual
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				x_i[total_index] += alpha*d[total_index];
				r[total_index] -= alpha*Atimesd[total_index];
				norm_r_new += r[total_index]*r[total_index];
			}
		}

		/*--- Compute the preconditioner multiplication (just one per iteration) ---*/
		PrecondVectorProduct(r, invMtimesr);
		
		/*--- Update r^T x invM x r ---*/
		beta_den = rinvMr; // old value necessary to compute beta
		rinvMr = 0;
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				rinvMr += r[total_index]*invMtimesr[total_index];
			}
		}

		/*--- Compute Gram-Schmidt constant (beta) ---*/
		beta = rinvMr / beta_den;

		/*--- Update search direction (d) ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				d[total_index] = invMtimesr[total_index] + beta*d[total_index];
			}
		}		

		/*--- Monitoring of the result ---*/
		norm = sqrt(norm_r_new);
		if ((monitoring == true) && (it%500 == 0))
			cout << "Prec. CG-Solution:: Iteration = " << it << " ; Norm of the residual: " << norm << endl;
		if (it == max_it) break;

	}

	cout << "Prec. CG-Solution:: Iteration = " << it << " ; Norm of the residual: " << norm << endl;

	delete [] invMtimesr;
	delete [] Atimesd;
	delete [] r;
	delete [] d;
}
