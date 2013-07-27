/*!
 * \file vector_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author Current Development: Stanford University.
 * \version 2.0.6
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

#include "../include/vector_structure.hpp"

CSysVector::CSysVector(void) {
  
  vec_val = NULL;
  
}

CSysVector::CSysVector(const unsigned int & size, const double & val) {
  
  nElm = size; nElmDomain = size;
  nBlk = nElm; nBlkDomain = nElmDomain;
  nVar = 1;
  
  /*--- Check for invalid size, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,double): "
    << "invalid input: size = " << size << endl;
    throw(-1);
  }

  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  
#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  unsigned long nElmLocal = (unsigned long)nElm;
  MPI::COMM_WORLD.Allreduce(&nElmLocal, &nElmGlobal, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
  
}

CSysVector::CSysVector(const unsigned int & numBlk, const unsigned int & numBlkDomain, const unsigned short & numVar,
                       const double & val) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;
  
  /*--- Check for invalid input, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,unsigned int,double): "
    << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << endl;
    throw(-1);
  }
	
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  
#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  unsigned long nElmLocal = (unsigned long)nElm;
  MPI::COMM_WORLD.Allreduce(&nElmLocal, &nElmGlobal, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
  
}

CSysVector::CSysVector(const CSysVector & u) {
  
  /*--- Copy size information, allocate memory, and initialize values ---*/
  nElm = u.nElm; nElmDomain = u.nElmDomain;
  nBlk = u.nBlk; nBlkDomain = u.nBlkDomain;
  nVar = u.nVar;
  
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];
  
#ifndef NO_MPI
  myrank = u.myrank;
  nElmGlobal = u.nElmGlobal;
#endif
  
}

CSysVector::CSysVector(const unsigned int & size, const double* u_array) {
  
  nElm = size; nElmDomain = size;
  nBlk = nElm; nBlkDomain = nElmDomain;
  nVar = 1;
  
  /*--- Check for invalid size, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int, double*): "
    << "invalid input: size = " << size << endl;
    throw(-1);
  }

  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];
  
#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  unsigned long nElmLocal = (unsigned long)nElm;
  MPI::COMM_WORLD.Allreduce(&nElmLocal, &nElmGlobal, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
  
}

CSysVector::CSysVector(const unsigned int & numBlk, const unsigned int & numBlkDomain, const unsigned short & numVar,
                       const double* u_array) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;
  
  /*--- check for invalid input, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,unsigned int,double*): "
    << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << endl;
    throw(-1);
  }

  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];
  
#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  unsigned long nElmLocal = (unsigned long)nElm;
  MPI::COMM_WORLD.Allreduce(&nElmLocal, &nElmGlobal, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
  
}

CSysVector::~CSysVector() {
  delete [] vec_val;
  nElm = -1;
	nElmDomain = -1;
  nBlk = -1;
  nBlkDomain = -1;
  nVar = -1;
#ifndef NO_MPI
  myrank = -1;
#endif
}

void CSysVector::Initialize(const unsigned int & numBlk, const unsigned int & numBlkDomain, const unsigned short & numVar, const double & val) {
  
  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;
  
  /*--- Check for invalid input, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,unsigned int,double): "
    << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << endl;
    throw(-1);
  }
	
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  
#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  unsigned long nElmLocal = (unsigned long)nElm;
  MPI::COMM_WORLD.Allreduce(&nElmLocal, &nElmGlobal, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
  
}

void CSysVector::Equals_AX(const double & a, CSysVector & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    cerr << "CSysVector::Equals_AX(): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i];
}

void CSysVector::Plus_AX(const double & a, CSysVector & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    cerr << "CSysVector::Plus_AX(): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] += a * x.vec_val[i];
}

void CSysVector::Equals_AX_Plus_BY(const double & a, CSysVector & x, const double & b, CSysVector & y) {
  /*--- check that *this, x and y are compatible ---*/
  if ((nElm != x.nElm) || (nElm != y.nElm)) {
    cerr << "CSysVector::Equals_AX_Plus_BY(): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i] + b * y.vec_val[i];
}

CSysVector & CSysVector::operator=(const CSysVector & u) {
  
  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if (this == &u) return *this;
  
  delete [] vec_val; // in case the size is different
  nElm = u.nElm;
  nElmDomain = u.nElmDomain;
  
  nBlk = u.nBlk;
	nBlkDomain = u.nBlkDomain;
  
  nVar = u.nVar;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];
  
#ifndef NO_MPI
  myrank = u.myrank;
  nElmGlobal = u.nElmGlobal;
#endif
  
  return *this;
}

CSysVector & CSysVector::operator=(const double & val) {
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  return *this;
}

CSysVector CSysVector::operator+(const CSysVector & u) const {
  
  /*--- Use copy constructor and compound addition-assignment ---*/
  CSysVector sum(*this);
  sum += u;
  return sum;
}

CSysVector & CSysVector::operator+=(const CSysVector & u) {
  
  /*--- Check for consistent sizes, then add elements ---*/
  if (nElm != u.nElm) {
    cerr << "CSysVector::operator+=(CSysVector): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] += u.vec_val[i];
  return *this;
}

CSysVector CSysVector::operator-(const CSysVector & u) const {
  
  /*--- Use copy constructor and compound subtraction-assignment ---*/
  CSysVector diff(*this);
  diff -= u;
  return diff;
}

CSysVector & CSysVector::operator-=(const CSysVector & u) {
  
  /*--- Check for consistent sizes, then subtract elements ---*/
  if (nElm != u.nElm) {
    cerr << "CSysVector::operator-=(CSysVector): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] -= u.vec_val[i];
  return *this;
}

CSysVector CSysVector::operator*(const double & val) const {
  
  /*--- use copy constructor and compound scalar
   multiplication-assignment ---*/
  CSysVector prod(*this);
  prod *= val;
  return prod;
}

CSysVector operator*(const double & val, const CSysVector & u) {
  
  /*--- use copy constructor and compound scalar
   multiplication-assignment ---*/
  CSysVector prod(u);
  prod *= val;
  return prod;
}

CSysVector & CSysVector::operator*=(const double & val) {
  
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] *= val;
  return *this;
}

CSysVector CSysVector::operator/(const double & val) const {
  
  /*--- use copy constructor and compound scalar
   division-assignment ---*/
  CSysVector quotient(*this);
  quotient /= val;
  return quotient;
}

CSysVector & CSysVector::operator/=(const double & val) {
  
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] /= val;
  return *this;
}

double CSysVector::norm() const {
  
  /*--- just call dotProd on this*, then sqrt ---*/
  double val = dotProd(*this,*this);
  if (val < 0.0) {
    cerr << "CSysVector::norm(): " << "inner product of CSysVector is negative";
    throw(-1);
  }
  return sqrt(val);
}

void CSysVector::CopyToArray(double* u_array) {
  
  for (unsigned int i = 0; i < nElm; i++)
    u_array[i] = vec_val[i];
}

void CSysVector::AddBlock(unsigned long val_ipoint, double *val_residual) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] += val_residual[iVar];
}

void CSysVector::SubtractBlock(unsigned long val_ipoint, double *val_residual) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] -= val_residual[iVar];
}

void CSysVector::SetBlock(unsigned long val_ipoint, double *val_residual) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] = val_residual[iVar];
}

void CSysVector::SetBlock(unsigned long val_ipoint, unsigned short val_var, double val_residual) {

  vec_val[val_ipoint*nVar+val_var] = val_residual;
}

void CSysVector::SetBlock_Zero(unsigned long val_ipoint) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] = 0.0;
}

void CSysVector::SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var) {
    vec_val[val_ipoint*nVar+val_var] = 0.0;
}

double CSysVector::GetBlock(unsigned long val_ipoint, unsigned short val_var) {
  return vec_val[val_ipoint*nVar + val_var];
}

double *CSysVector::GetBlock(unsigned long val_ipoint) {
  return &vec_val[val_ipoint*nVar];
}

double dotProd(const CSysVector & u, const CSysVector & v) {
  
  /*--- check for consistent sizes ---*/
  if (u.nElm != v.nElm) {
    cerr << "CSysVector friend dotProd(CSysVector,CSysVector): "
    << "CSysVector sizes do not match";
    throw(-1);
  }
  
  /*--- find local inner product and, if a parallel run, sum over all
   processors (we use nElemDomain instead of nElem) ---*/
  double loc_prod = 0.0;
  for (unsigned int i = 0; i < u.nElmDomain; i++)
    loc_prod += u.vec_val[i]*v.vec_val[i];
  double prod = 0.0;
  
#ifndef NO_MPI
  MPI::COMM_WORLD.Allreduce(&loc_prod, &prod, 1, MPI::DOUBLE, MPI::SUM);
#else
  prod = loc_prod;
#endif
  
  return prod;
}