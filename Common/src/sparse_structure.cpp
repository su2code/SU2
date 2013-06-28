/*!
 * \file sparse_structure.cpp
 * \brief Main subroutines for doing the sparse structures.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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
		unsigned long* val_col_ind, unsigned long val_nnz, bool preconditioner, bool blockDiagonal, unsigned short val_nSub_blocks,  unsigned short *val_Sub_block_sizes) {

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
	nSub_blocks = val_nSub_blocks;
	Sub_block_sizes = new  unsigned short[nSub_blocks];
	for (unsigned short isize = 0; isize < nSub_blocks; isize++)
		Sub_block_sizes[isize] = val_Sub_block_sizes[isize];
	if (preconditioner)
		invM = new double [nPoint*nVar*nVar];	// Reserve memory for the values of the inverse of the preconditioner

	blockDiagonalJacobian = blockDiagonal;

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

void CSparseMatrix::AddVal2Diag(unsigned long block_i,  double* val_val, unsigned short num_dim) {
	unsigned long step = 0;
	unsigned long iVar, iSpecies;
	for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
		step++;
		if (col_ind[index] == block_i) {	// Only elements on the diagonal
			for (iVar = 0; iVar < nVar; iVar++) {
				iSpecies = iVar/(num_dim + 2);
				val[(row_ptr[block_i]+step-1)*nVar*nVar+iVar*nVar+iVar] += val_val[iSpecies];
			}
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
	unsigned short jVar, kVar, iBlock_start_pos = 0;
	short iVar;
	double weight;

	GetBlock(block_i,block_i);

	if (nVar == 1)
		rhs[0] /= (block[0]+EPS*EPS);
	else {
		if (blockDiagonalJacobian) {
			for (short iBlock = 0; iBlock < nSub_blocks; iBlock ++ ) {
				iBlock_start_pos = iBlock_start_pos + iBlock*Sub_block_sizes[iBlock];
				for (iVar = iBlock_start_pos +1; iVar < iBlock_start_pos + Sub_block_sizes[iBlock]; iVar++) {
					for (jVar = iBlock_start_pos; jVar < iVar; jVar++) {
						weight = block[iVar*nVar+jVar]/(block[jVar*nVar+jVar]+EPS*EPS);
						for (kVar = jVar; kVar < iBlock_start_pos + Sub_block_sizes[iBlock]; kVar++)
							block[iVar*nVar+kVar] -= weight*block[jVar*nVar+kVar];
						rhs[iVar] -= weight*rhs[jVar];
					}
				}
			}

			double aux;
			unsigned short iBlock_end_pos;
			for (unsigned short iBlock = 0; iBlock < nSub_blocks; iBlock ++ ) {
				iBlock_end_pos = (iBlock+1)*Sub_block_sizes[iBlock]-1;
				rhs[iBlock_end_pos] = rhs[iBlock_end_pos]/(block[iBlock_end_pos*nVar+iBlock_end_pos]+EPS*EPS);
				for (iVar = iBlock_end_pos-1; iVar >= iBlock*Sub_block_sizes[iBlock]; iVar--) {
					aux = 0;
					for (jVar = iVar+1; jVar <= iBlock_end_pos; jVar++)
						aux += block[iVar*nVar+jVar]*rhs[jVar];

					rhs[iVar] = (rhs[iVar]-aux)/(block[iVar*nVar+iVar]+EPS*EPS);
					if (iVar == 0) break;
				}
			}
		}
		else {
			/*--- Transform system in Upper Matrix ---*/
			for (iVar = 1; iVar < (short)nVar; iVar++) {
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

void CSparseMatrix::SendReceive_Solution(double* x, CGeometry *geometry, CConfig *config) {

#ifndef NO_MPI
	unsigned long iPoint, iVertex;
	unsigned short iVar, iMarker;
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				
				/*--- Allocate buffer vector ---*/
				unsigned long nBuffer_Vector = nVertex*nVar;
				int send_to = SendRecv-1;
				double *Buffer_Send_X = new double[nBuffer_Vector];
				
				/*--- Copy the solution to the buffer vector ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Send_X[iVar*nVertex+iVertex] = x[iPoint*nVar+iVar];
				}
				
				/*--- Send the buffer information ---*/
				MPI::COMM_WORLD.Bsend(Buffer_Send_X, nBuffer_Vector, MPI::DOUBLE, send_to, 0);

				/*--- Deallocate buffer vector ---*/
				delete [] Buffer_Send_X;
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				
				/*--- Allocate buffer vector ---*/
				int receive_from = abs(SendRecv)-1;
				unsigned long nBuffer_Vector = nVertex*nVar;
				double *Buffer_Receive_X = new double [nBuffer_Vector];
				
				/*--- Receive the information ---*/
				MPI::COMM_WORLD.Recv(Buffer_Receive_X, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);
				
				/*--- Store the received information ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					for (iVar = 0; iVar < nVar; iVar++)
						x[iPoint*nVar+iVar] = Buffer_Receive_X[iVar*nVertex+iVertex];
					
				}
				
				/*--- Deallocate buffer vector ---*/
				delete [] Buffer_Receive_X;
			}		
		}
	
	MPI::COMM_WORLD.Barrier();
	
#endif
}

void CSparseMatrix::LU_SGSIteration(double* b, double* x_n, CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, iVar;
	
	/*--- First part of the symmetric iteration: (D+L).x* = b ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		LowerProduct(x_n,iPoint);  // compute L.x*
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] = b[iPoint*nVar+iVar] - prod_row_vector[iVar]; // compute aux_vector = b-L.x*
		Gauss_Elimination(iPoint, aux_vector); // solve D.x* = aux_vector
		for (iVar = 0; iVar < nVar; iVar++)
			x_n[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x* = solution
	}	
	
	/*--- Send-Receive the solution vector (MPI) ---*/
	SendReceive_Solution(x_n, geometry, config);

	/*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/
	for (iPoint = geometry->GetnPointDomain()-1; iPoint >= 0; iPoint--) {
		DiagonalProduct(x_n,iPoint);  // compute D.x*
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] = prod_row_vector[iVar]; // compute aux_vector = D.x*
		UpperProduct(x_n,iPoint);  // compute U.x_(n+1)
		for (iVar = 0; iVar < nVar; iVar++)
			aux_vector[iVar] -= prod_row_vector[iVar]; // compute aux_vector = D.x*-U.x_(n+1)
		Gauss_Elimination(iPoint, aux_vector); // solve D.x* = aux_vector
		for (iVar = 0; iVar < nVar; iVar++)
			x_n[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x_(1) = solution
		if (iPoint == 0) break;
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

void CSparseMatrix::MatrixVectorProduct(double* vec, double* prod, double nPoint) {
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

void CSparseMatrix::PrecondVectorProduct(double* vec, double* prod, double nPoint) {
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

void CSparseMatrix::SGSSolution(double* b, double* x_i, double tol, int max_it, 
																bool monitoring, CGeometry *geometry, CConfig *config) {
	int iter;
	unsigned long iPoint, iVar;
	double norm_1 = 0.0, norm_2 = 0.0, norm = 0.0, my_norm;
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
	for (iter = 0; iter < max_it; iter++) {
		
		norm_1 = 0.0; norm_2 = 0.0;
		
		/*--- First part of the symmetric iteration: (D+L).x* = b-U.x_n ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			UpperProduct(x_i,iPoint);  // compute U.x_n (result written in prod_row_vector)
			for (iVar = 0; iVar < nVar; iVar++)
				aux_vector[iVar] = b[iPoint*nVar+iVar] - prod_row_vector[iVar]; // compute aux_vector = b-U.x_n
			LowerProduct(x_i,iPoint);  // compute L.x* (idem)
			for (iVar = 0; iVar < nVar; iVar++)
				aux_vector[iVar] -= prod_row_vector[iVar]; // compute aux_vector = b-U.x_n-L.x*
			Gauss_Elimination(iPoint, aux_vector); // solve D.x* = aux_vector
			for (iVar = 0; iVar < nVar; iVar++) {
				norm_1 += (x_i[iPoint*nVar+iVar]-aux_vector[iVar])*(x_i[iPoint*nVar+iVar]-aux_vector[iVar]);
				x_i[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x* = solution
			}
		}
		norm_1 = sqrt(norm_1);
		
		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(x_i, geometry, config);
		
		/*--- Second part of the symmetric iteration: (D+U).x_(n+1) = b-L·x* ---*/
		for (iPoint = geometry->GetnPointDomain()-1; iPoint >= 0; iPoint--) {
			LowerProduct(x_i,iPoint);  // compute L.x* (result written in prod_row_vector)
			for (iVar = 0; iVar < nVar; iVar++)
				aux_vector[iVar] = b[iPoint*nVar+iVar] - prod_row_vector[iVar]; // compute aux_vector = b-L.x*
			UpperProduct(x_i,iPoint);  // compute U.x_(n+1) (idem)
			for (iVar = 0; iVar < nVar; iVar++)
				aux_vector[iVar] -= prod_row_vector[iVar]; // compute aux_vector = b-L.x*-U.x_(n+1)
			Gauss_Elimination(iPoint, aux_vector); // solve D.x_(n+1) = aux_vector
			for (iVar = 0; iVar < nVar; iVar++) {
				norm_2 += (x_i[iPoint*nVar+iVar]-aux_vector[iVar])*(x_i[iPoint*nVar+iVar]-aux_vector[iVar]);
				x_i[iPoint*nVar+iVar] = aux_vector[iVar]; // assesing x* = solution
			}
			if (iPoint == 0) break;
		}
		norm_2 = sqrt(norm_2);
		
		my_norm = (norm_1+norm_2);
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_norm, &norm, 1, MPI::DOUBLE, MPI::MIN); 
#else
		norm = my_norm;
#endif
		
//		if ((monitoring == true) && (iter%10 == 0) && (rank == MASTER_NODE))
//			cout << "SGS-Solution:: Iteration = " << iter << " ; Norm of the residual: " << norm << endl;
		
		if (norm < tol) break;
		
	}
	
	if ((monitoring == true) && (rank == MASTER_NODE))
		cout << "SGS-Solution => Iter: " << iter << " | Norm Res: " << norm << endl;
	return;
	
}

void CSparseMatrix::CGSolution(double* b, double* x_i, double tol, int max_it, bool monitoring, 
															 CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, iVar, total_index;
	unsigned short iter = 0;
	double *d, *r, *Atimesd, *invMtimesr, rinvMr, my_rinvMr, norm = 0.0, alpha, alpha_den, beta, norm_r_new, beta_den, 
	my_alpha_den, my_norm_r_new;
	
	d = new double [nPoint * nVar];
	r = new double [nPoint * nVar];
	Atimesd = new double [nPoint * nVar];
	invMtimesr = new double [nPoint * nVar];
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
	BuildJacobiPreconditioner();
	
	/*--- First iteration: determine residual (r) and search direction (d) ---*/
	/*--- r(0) = b-Ax(0) ---*/
	MatrixVectorProduct(x_i, r, geometry->GetnPointDomain()); // r = Ax(0)
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			r[total_index] = b[total_index] - r[total_index]; // r(0) = b-A·x(0);
		}
	}
	
	/*--- Send-Receive the solution vector (MPI) ---*/
	SendReceive_Solution(r, geometry, config);

	/*--- d(0) = invM*r(0) ---*/
	PrecondVectorProduct(r, invMtimesr, geometry->GetnPointDomain());
	
	/*--- Send-Receive the solution vector (MPI) ---*/
	SendReceive_Solution(invMtimesr, geometry, config);

	my_rinvMr = 0;
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			d[total_index] = invMtimesr[total_index];
			my_rinvMr += r[total_index]*invMtimesr[total_index];
		}
	}
	
#ifndef NO_MPI
	MPI::COMM_WORLD.Allreduce(&my_rinvMr, &rinvMr, 1, MPI::DOUBLE, MPI::SUM); 	
#else
	rinvMr = my_rinvMr;
#endif
	
	/*--- Send-Receive the solution vector (MPI) ---*/
	SendReceive_Solution(d, geometry, config);

	/*--- CG-algorithm ---*/
	for (iter = 0; iter < max_it; iter++) {
		
		/*--- Compute descend step (alpha) ---*/
		MatrixVectorProduct(d, Atimesd, geometry->GetnPointDomain());
		my_alpha_den = 0.0; 
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				my_alpha_den += d[total_index]*Atimesd[total_index];
			}
		}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_alpha_den, &alpha_den, 1, MPI::DOUBLE, MPI::SUM); 	
#else
		alpha_den = my_alpha_den;
#endif
		
		if (fabs(alpha_den) < 1E-30) alpha_den = 1E-30;

		alpha = rinvMr / alpha_den;
		
		/*--- Update solution and residual ---*/
		my_norm_r_new = 0.0; // we need the norm of the updated residual
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				x_i[total_index] += alpha*d[total_index];
				r[total_index] -= alpha*Atimesd[total_index];
				my_norm_r_new += r[total_index]*r[total_index];
			}
		}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_norm_r_new, &norm_r_new, 1, MPI::DOUBLE, MPI::SUM); 	
#else
		norm_r_new = my_norm_r_new;
#endif
		
		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(x_i, geometry, config);
		SendReceive_Solution(r, geometry, config);

		/*--- Compute the preconditioner multiplication (just one per iteration) ---*/
		PrecondVectorProduct(r, invMtimesr, geometry->GetnPointDomain());
		
		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(invMtimesr, geometry, config);

		/*--- Update r^T x invM x r ---*/
		beta_den = rinvMr; // old value necessary to compute beta
		my_rinvMr = 0;
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				my_rinvMr += r[total_index]*invMtimesr[total_index];
			}
		}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_rinvMr, &rinvMr, 1, MPI::DOUBLE, MPI::SUM); 	
#else
		rinvMr = my_rinvMr;
#endif
		
		if (fabs(beta_den)< 1E-30) beta_den = 1E-30;

		/*--- Compute Gram-Schmidt constant (beta) ---*/
		beta = rinvMr / beta_den;
		
		/*--- Update search direction (d) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				d[total_index] = invMtimesr[total_index] + beta*d[total_index];
			}
		}		
		
		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(d, geometry, config);

		/*--- Monitoring of the result ---*/
		norm = sqrt(norm_r_new);
		
//		if ((monitoring == true) && (iter%10 == 0) && (rank == MASTER_NODE))
//			cout << "Prec. CG-Solution:: Iteration = " << iter << " ; Norm of the residual: " << norm << endl;
		
		if (norm < tol) break;
		
	}
	
	delete [] invMtimesr; delete [] Atimesd;
	delete [] r; delete [] d;
	
	if ((monitoring == true) && (rank == MASTER_NODE))
		cout << "Prec. CG-Solution => Iter: " << iter << " | Norm Res: " << norm << endl;
	return;	
	
}



void CSparseMatrix::BCGSTABSolution(double* b, double* x_i, double tol, int max_it, 
																	 bool monitoring, CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, iVar, total_index;
	double *p, *v, *phat, *shat, *res, *res0, *t, *s, *tmp, rho = 0.0, alpha, omega, norm = 0.0, my_norm, beta,
	my_buf[2], buf[2] = {0.0, 0.0}, my_gamma, gamma = 0.0, rho_prime, my_rho;
  unsigned short iter = 0;
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif

	p = new double [geometry->GetnPointDomain() * nVar];
	v = new double [geometry->GetnPointDomain() * nVar];
	phat = new double [geometry->GetnPoint() * nVar];
	shat = new double [geometry->GetnPoint() * nVar];
	res = new double [geometry->GetnPointDomain() * nVar];
	res0 = new double [geometry->GetnPointDomain() * nVar];
	t = new double [geometry->GetnPointDomain() * nVar];
	s = new double [geometry->GetnPointDomain() * nVar];
	tmp = new double [nVar];
	
	BuildJacobiPreconditioner();

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
      p[total_index] = 0.0;	phat[total_index] = 0.0;
      v[total_index] = 0.0; shat[total_index] = 0.0;
    }
		
  rho = 1.0; alpha = 1.0; omega = 1.0;
	
	MatrixVectorProduct(x_i, res, geometry->GetnPointDomain());
	
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
      res[total_index] = b[total_index]-res[total_index];
      res0[total_index] = res[total_index];
    }
	
  my_norm= 0.0;
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
      tmp[iVar] = res[total_index];
		}
		Gauss_Elimination(iPoint, tmp);
    for (iVar = 0; iVar < nVar; iVar++)
      my_norm += tmp[iVar]*tmp[iVar];
  }
		
	for (iter = 0; iter < max_it; iter++) {
		
    rho_prime = -omega*rho;
		
    my_rho = 0.0;
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
        my_rho += res[total_index]*res0[total_index];
			}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_rho, &rho, 1, MPI::DOUBLE, MPI::SUM); 	
#else
		rho = my_rho;
#endif
		
    if (fabs(rho_prime) < 1E-30) rho_prime = 1E-30;
		
    beta = alpha*rho/rho_prime;
		
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
        p[total_index] = res[total_index]-beta*(p[total_index]-omega*v[total_index]);
			}
		
    /*--- Block-diagonal preconditioning (solve the system) ---*/
		PrecondVectorProduct(p, phat, geometry->GetnPointDomain());

		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(phat, geometry, config);
		
		MatrixVectorProduct(phat, v, geometry->GetnPointDomain());
		
    my_gamma = 0.0;
		
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
        my_gamma += v[total_index]*res0[total_index];
			}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_gamma, &gamma, 1, MPI::DOUBLE, MPI::SUM); 	
#else
		gamma = my_gamma;
#endif

    if (fabs(gamma)< 1E-30) gamma = 1E-30;
		
    alpha = rho/gamma;
		
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;				
        s[total_index] = res[total_index]-alpha*v[total_index];
			}
		
    /*--- Block-diagonal preconditioning (solve the system) ---*/
		PrecondVectorProduct(s, shat, geometry->GetnPointDomain());

		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(shat, geometry, config);
		
		MatrixVectorProduct(shat, t, geometry->GetnPointDomain());
		
		my_buf[0] = 0.0;
		my_buf[1] = 0.0;
		
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;		
        my_buf[0] += s[total_index]*t[total_index];
        my_buf[1] += t[total_index]*t[total_index];
      }
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(my_buf, buf, 2, MPI::DOUBLE, MPI::SUM); 	
#else
		buf[0] = my_buf[0]; buf[1] = my_buf[1];
#endif
				
		if (fabs(buf[1]) < 1E-30) buf[1] = 1E-30;

    omega = buf[0]/buf[1];
		
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
        x_i[total_index] += alpha*phat[total_index]+omega*shat[total_index];
			}
		
		/*--- Send-Receive the solution vector (MPI) ---*/
		SendReceive_Solution(x_i, geometry, config);
				
		MatrixVectorProduct(x_i, res, geometry->GetnPointDomain());
		
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;				
        res[total_index] = b[total_index]-res[total_index];
			}
		
		my_norm = 0.0;
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {			
      for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
        tmp[iVar] = res[total_index];
			}
			Gauss_Elimination(iPoint, tmp);
      for (iVar = 0; iVar < nVar; iVar++)
        my_norm += tmp[iVar]*tmp[iVar];
		}
			
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_norm, &norm, 1, MPI::DOUBLE, MPI::SUM); 
#else
		norm = my_norm;
#endif
		
		norm = sqrt(norm);
		
//		if ((monitoring == true) && (iter%10 == 0) && (rank == MASTER_NODE))
//			cout << "Prec. BCGSTAB-Solution:: Iteration = " << iter << " ; Norm of the residual: " << norm << endl;

		if (norm < tol) break;
		
  }
	
	delete[] res; delete[] res0;
	delete[] p; delete[] v; delete[] t; delete[] s;
	delete[] phat; delete[] shat; delete[] tmp;

	if ((monitoring == true) && (rank == MASTER_NODE))
		cout << "Prec. BCGSTAB-Solution => Iter: " << iter << " | Norm Res: " << norm << endl;
	return;
}
