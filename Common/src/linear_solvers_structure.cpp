/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author Current Development: Stanford University.
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

#include "../include/linear_solvers_structure.hpp"

CSysVector::CSysVector(const unsigned int & size, const double & val) {

  /*--- check for invalid size, then allocate memory and initialize
        values ---*/
  if ( (size <= 0) || (size >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,double): " 
	 << "invalid input: size = " << size << endl;
    throw(-1);
  }
  nElm = size;
  nBlk = nElm;
  nVar = 1;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;

#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  // TODO(J. Hicken): the following should probably be a static
  // element of the class, otherwise it would be too expensive to call
  // each time a new CSysVector is created
  //unsigned long nElmLocal = (unsigned long)nElm;
  //MPI::COMM_WORLD.Allreduce(nElmLocal,nElmGlobal,1,MPI_UNSIGNED_LONG,MPI_SUM);
#endif

}

CSysVector::CSysVector(const unsigned int & numBlk, const unsigned short & numVar, 
		 const double & val) {

  /*--- check for invalid input, then allocate memory and initialize
        values ---*/
  nElm = numBlk*numVar;
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,unsigned int,double): " 
	 << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << endl;
    throw(-1);
  }
  nBlk = numBlk;
  nVar = numVar;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;

#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  //unsigned long nElmLocal = (unsigned long)nElm;
  //MPI::COMM_WORLD.Allreduce(nElmLocal,nElmGlobal,1,MPI_UNSIGNED_LONG,MPI_SUM);
#endif
}

CSysVector::~CSysVector() {
  delete [] vec_val;
  nElm = -1;
  nBlk = -1;
  nVar = -1;
#ifndef NO_MPI
  myrank = -1;
#endif
}

CSysVector::CSysVector(const CSysVector & u) {
  
  /*--- copy size information, allocate memory, and initialize
        values ---*/
  nElm = u.nElm;
  nBlk = u.nBlk;
  nVar = u.nVar;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];

#ifndef NO_MPI
  myrank = u.myrank;
  //nElmGlobal = u.nElmGlobal;
#endif
}

CSysVector::CSysVector(const unsigned int & size, const double* u_array) {

  /*--- check for invalid size, then allocate memory and initialize
        values ---*/
  if ( (size <= 0) || (size >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int, double*): " 
	 << "invalid input: size = " << size << endl;
    throw(-1);
  }
  nElm = size;
  nBlk = nElm;
  nVar = 1;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];

#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  //unsigned long nElmLocal = (unsigned long)nElm;
  //MPI::COMM_WORLD.Allreduce(nElmLocal,nElmGlobal,1,MPI_UNSIGNED_LONG,MPI_SUM);
#endif
  
}

CSysVector::CSysVector(const unsigned int & numBlk, const unsigned short & numVar, 
		 const double* u_array) {

  /*--- check for invalid input, then allocate memory and initialize
        values ---*/
  nElm = numBlk*numVar;
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    cerr << "CSysVector::CSysVector(unsigned int,unsigned int,double*): " 
	 << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << endl;
    throw(-1);
  }
  nBlk = numBlk;
  nVar = numVar;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];

#ifndef NO_MPI
  myrank = MPI::COMM_WORLD.Get_rank();
  //unsigned long nElmLocal = (unsigned long)nElm;
  //MPI::COMM_WORLD.Allreduce(nElmLocal,nElmGlobal,1,MPI_UNSIGNED_LONG,MPI_SUM);
#endif
}

void CSysVector::Equals_AX(const double & a, CSysVector & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    cerr << "CSysVector::Equals_AX(): " 
	 << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i];
}

void CSysVector::Plus_AX(const double & a, CSysVector & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    cerr << "CSysVector::Plus_AX(): " 
	 << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] += a * x.vec_val[i];
}

void CSysVector::Equals_AX_Plus_BY(const double & a, CSysVector & x, const double & b, CSysVector & y) {
  /*--- check that *this, x and y are compatible ---*/
  if ((nElm != x.nElm) || (nElm != y.nElm)) {
    cerr << "CSysVector::Equals_AX_Plus_BY(): " 
	 << "sizes do not match";
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
  nBlk = u.nBlk;
  nVar = u.nVar;
  vec_val = new double[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];

#ifndef NO_MPI
  myrank = u.myrank;
  //nElmGlobal = u.nElmGlobal;
#endif

  return *this;
}

CSysVector & CSysVector::operator=(const double & val) {
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  return *this;
}

CSysVector CSysVector::operator+(const CSysVector & u) const {

  /*--- use copy constructor and compound addition-assignment ---*/
  CSysVector sum(*this);
  sum += u;
  return sum;
}

CSysVector & CSysVector::operator+=(const CSysVector & u) {
  
  /*--- check for consistent sizes, then add elements ---*/
  if (nElm != u.nElm) {
    cerr << "CSysVector::operator+=(CSysVector): " 
	 << "sizes do not match";
    throw(-1);
  }
  for (unsigned int i = 0; i < nElm; i++) 
    vec_val[i] += u.vec_val[i];
  return *this;
}

CSysVector CSysVector::operator-(const CSysVector & u) const {

  /*--- use copy constructor and compound subtraction-assignment ---*/
  CSysVector diff(*this);
  diff -= u;
  return diff;
}

CSysVector & CSysVector::operator-=(const CSysVector & u) {
  
  /*--- check for consistent sizes, then subtract elements ---*/
  if (nElm != u.nElm) {
    cerr << "CSysVector::operator-=(CSysVector): " 
	 << "sizes do not match";
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
    cerr << "CSysVector::norm(): " 
	 << "inner product of CSysVector is negative";
    throw(-1);
  }
  return sqrt(val);
}

void CSysVector::CopyToArray(double* u_array) {

  for (unsigned int i = 0; i < nElm; i++)
    u_array[i] = vec_val[i];
}

double dotProd(const CSysVector & u, const CSysVector & v) {

  /*--- check for consistent sizes ---*/
  if (u.nElm != v.nElm) {
    cerr << "CSysVector friend dotProd(CSysVector,CSysVector): " 
	 << "CSysVector sizes do not match";
    throw(-1);
  }

  /*--- find local inner product and, if a parallel run, sum over all
        processors ---*/
  double loc_prod = 0.0;
  for (unsigned int i = 0; i < u.nElm; i++) 
    loc_prod += u.vec_val[i]*v.vec_val[i];
  double prod = 0.0;
#ifndef NO_MPI
//  MPI::COMM_WORLD.Allreduce(loc_prod,prod,1,MPI_DOUBLE,MPI_SUM);
#else
  prod = loc_prod;
#endif
  return prod;
}

// ======================================================================
// start of routines for CSysSolve class

void CSysSolve::applyGivens(const double & s, const double & c, double & h1, double & h2) {

  double temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

void CSysSolve::generateGivens(double & dx, double & dy, double & s, double & c) {

  if ( (dx == 0.0) && (dy == 0.0) ) {
    c = 1.0;
    s = 0.0;
  }
  else if ( fabs(dy) > fabs(dx) ) {
    double tmp = dx/dy;
    dx = sqrt(1.0 + tmp*tmp);
    s = sign(1.0/dx,dy);
    c = tmp*s;
  }
  else if ( fabs(dy) <= fabs(dx) ) {
    double tmp = dy/dx;
    dy = sqrt(1.0 + tmp*tmp);
    c = sign(1.0/dy,dx);
    s = tmp*c;
  }
  else {
    // dx and/or dy must be invalid
    dx = 0.0;
    dy = 0.0;
    c = 1.0;
    s = 0.0;
  }
  dx = fabs(dx*dy);
  dy = 0.0;
}

void CSysSolve::solveReduced(const int & n, const vector<vector<double> > & Hsbg, 
			     const vector<double> & rhs, vector<double> & x) {
  // initialize... 
  for (int i = 0; i < n; i++)
    x[i] = rhs[i];
  // ... and backsolve
  for (int i = n-1; i >= 0; i--) {
    x[i] /= Hsbg[i][i];
    for (int j = i-1; j >= 0; j--) {
      x[j] -= Hsbg[j][i]*x[i];
    }
  }
}

void CSysSolve::modGramSchmidt(int i, vector<vector<double> > & Hsbg, vector<CSysVector> & w) {

  // parameter for reorthonormalization
  static const double reorth = 0.98; 

  // get the norm of the vector being orthogonalized, and find the 
  // threshold for re-orthogonalization
  double nrm = dotProd(w[i+1],w[i+1]);
  double thr = nrm*reorth;
  if (nrm <= 0.0) {
    // the norm of w[i+1] < 0.0
    cerr << "CSysSolve::modGramSchmidt: dotProd(w[i+1],w[i+1]) < 0.0" << endl;
    throw(-1);
  }
  else if (nrm != nrm) { 
    // this is intended to catch if nrm = NaN, but some optimizations
    // may mess it up (according to posts on stackoverflow.com)
    cerr << "CSysSolve::modGramSchmidt: w[i+1] = NaN" << endl;
    throw(-1);
  }

  // begin main Gram-Schmidt loop
  for (int k = 0; k < i+1; k++) {
    double prod = dotProd(w[i+1],w[k]);
    Hsbg[k][i] = prod;
    w[i+1].Plus_AX(-prod, w[k]);
    
    // check if reorthogonalization is necessary
    if (prod*prod > thr) {
      prod = dotProd(w[i+1],w[k]);
      Hsbg[k][i] += prod;
      w[i+1].Plus_AX(-prod, w[k]);
    }

    // update the norm and check its size
    nrm -= Hsbg[k][i]*Hsbg[k][i];
    if (nrm < 0.0) nrm = 0.0;
    thr = nrm*reorth;
  }

  // test the resulting vector
  nrm = w[i+1].norm();
  Hsbg[i+1][i] = nrm;
  if (nrm <= 0.0) {
    // w[i+1] is a linear combination of the w[0:i]
    cerr << "CSysSolve::modGramSchmidt: w[i+1] linearly dependent on w[0:i]" << endl;
    throw(-1);
  }

  // scale the resulting vector
  w[i+1] /= nrm;
}

void CSysSolve::writeHeader(ostream & os, const std::string & solver, 
			    const double & restol, const double & resinit) {
  if (!(os.good())) {
    cerr << "CSysSolve::writeHeader(): "
	 << "ostream is not good for i/o operations." << endl;
    throw(-1);
  }

  std::ios_base::fmtflags old_fmt;
  old_fmt = os.setf(ios::scientific, ios::floatfield);

  /* os << "# " << solver << " residual history" << endl;
  os << "# " << "residual tolerance target = " 
     << std::setw(10) << std::setprecision(4) << restol << endl;
  os << "# " << "initial residual norm     = " 
     << std::setw(10) << std::setprecision(4) << resinit << endl;
  os << std::setw(6) << "# iter" << std::setw(12) << "rel. res." << endl; */
  
  os.setf(old_fmt, ios::floatfield);
}

void CSysSolve::writeHistory(ostream & os, const int & iter, 
			     const double & res, const double & resinit) {
  if (!(os.good())) {
    cerr << "CSysSolve::writeHistory: "
	 << "ostream is not good for i/o operations." << endl;
    throw(-1);
  }
 
  std::ios_base::fmtflags old_fmt;
  old_fmt = os.setf(ios::scientific, ios::floatfield);

  os << std::setw(5) << iter
     << std::setw(12) << std::setprecision(4) << res/resinit 
     << endl;

  os.setf(old_fmt, ios::floatfield);
}

void CSysSolve::ConjugateGradient(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                  CPreconditioner & precond, double tol, int m, bool monitoring) {
  // check the subspace size
  if (m < 1) {
    cerr << "CSysSolve::ConjugateGradient: illegal value for subspace size, m = " << m << endl;
    throw(-1);
  }

  CSysVector r(b);
  CSysVector A_p(b);

  // calculate the initial residual, compute norm, and check if system is already solved
  mat_vec(x,A_p);
  r -= A_p; // recall, r holds b initially
  double norm_r = r.norm();  
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
    return;
  }  

  double alpha, beta, r_dot_z;
  CSysVector z(r);
  precond(r, z);
  CSysVector p(z);

  // output header information including initial residual
  int i = 0;
  if (monitoring) {
    writeHeader(cout, "CG", tol, norm_r);
    writeHistory(cout, i, norm_r, norm0);
  }

  // loop over all serach directions
  for (i = 0; i < m; i++) {

    // apply matrix to p to build Krylov subspace
    mat_vec(p, A_p);

    // calculate step-length alpha
    r_dot_z = dotProd(r, z);
    alpha = dotProd(A_p, p);
    alpha = r_dot_z / alpha;

    // update solution and residual: 
    x.Plus_AX(alpha, p);
    r.Plus_AX(-alpha, A_p);

    // check if solution has converged, else output the relative residual if necessary
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if ((monitoring) && ((i+1) % 100 == 0)) writeHistory(cout, i+1, norm_r, norm0);    
 
    precond(r, z);

    // calculate Gram-Schmidt coefficient beta
    // beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i})
    beta = 1.0 / r_dot_z;
    r_dot_z = dotProd(r, z);
    beta *= r_dot_z;

    // Gram-Schmidt orthogonalization; p = beta *p + z
    p.Equals_AX_Plus_BY(beta, p, 1.0, z);
  }

  // recalculate final residual (this should be optional)
  mat_vec(x, A_p);
  r = b;
  r -= A_p;

  double true_res = r.norm();
  
  if (monitoring) {
    cout << "# Conjugate Gradient final (true) residual:" << endl;
    cout << "# iteration = " << i << ": |res|/|res0| = "  << true_res/norm0 << endl;
  }

  if (fabs(true_res - norm_r) > tol*10.0) {
    cout << "# WARNING in CSysSolve::ConjugateGradient(): " << endl;
    cout << "# true residual norm and calculated residual norm do not agree." << endl;
    cout << "# true_res - calc_res = " << true_res - norm_r << endl;
  }
}

void CSysSolve::FlexibleGMRES(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                              CPreconditioner & precond, double tol, int m, bool monitoring) {

  // check the subspace size
  if (m < 1) {
    cerr << "CSysSolve::FlexibleGMRES: illegal value for subspace size, m = " << m << endl;
    throw(-1);
  }
  
  // define various arrays
  // Note: elements in w and z are initialized to x to avoid creating
  // a temporary CSysVector object for the copy constructor
  vector<CSysVector> w(m+1,x);
  vector<CSysVector> z(m+1,x);
  vector<double> g(m+1,0.0);
  vector<double> sn(m+1,0.0);
  vector<double> cs(m+1,0.0);
  vector<double> y(m,0.0);
  vector<vector<double> > H(m+1,vector<double>(m,0.0));
  
  // calculate the norm of the rhs vector
  double norm0 = b.norm();

  // calculate the initial residual (actually the negative residual)
  // and compute its norm
  mat_vec(x,w[0]);
  w[0] -= b;

  double beta = w[0].norm();  
  if ( (beta < tol*norm0) || (beta < eps) ) {
    // system is already solved
    cout << "CSysSolve::FlexibleGMRES(): system solved by initial guess." << endl;
    return;
  }

  cout << "beta " << beta << endl;
  cout << "eps*norm0 " << eps*norm0 << endl;

  // normalize residual to get w_{0} (the negative sign is because w[0]
  // holds the negative residual, as mentioned above)
  w[0] /= -beta;
  
  // initialize the RHS of the reduced system
  g[0] = beta;

  // output header information including initial residual
  int i = 0;
  if (monitoring) {
    writeHeader(cout, "FGMRES", tol, beta);
    writeHistory(cout, i, beta, norm0);
  }

  // loop over all serach directions
  for (i = 0; i < m; i++) {
    
    // check if solution has converged
    if (beta < tol*norm0) break;

    // precondition the CSysVector w[i] and store result in z[i]
    precond(w[i],z[i]);

    // add to Krylov subspace
    mat_vec(z[i],w[i+1]);

    // modified Gram-Schmidt orthogonalization
    modGramSchmidt(i,H,w);

    // apply old Givens rotations to new column of the Hessenberg matrix
    // then generate the new Givens rotation matrix and apply it to
    // the last two elements of H[:][i] and g
    for (int k = 0; k < i; k++)
      applyGivens(sn[k],cs[k],H[k][i],H[k+1][i]);
    generateGivens(H[i][i],H[i+1][i],sn[i],cs[i]);
    applyGivens(sn[i],cs[i],g[i],g[i+1]);

    // set L2 norm of residual and check if solution has converged
    beta = fabs(g[i+1]);

    // output the relative residual if necessary
    if ((monitoring) && ((i+1) % 100 == 0)) writeHistory(cout,i+1,beta,norm0);
  }

  // solve the least-squares system and update solution
  solveReduced(i,H,g,y);
  for (int k = 0; k < i; k++) {
    x.Plus_AX(y[k], z[k]);
  }
  
  // recalculate final (neg.) residual (this should be optional)
  mat_vec(x,w[0]);
  // w[0] = w[0] - b = A*x -b = -r
  w[0] -= b;

  double res = w[0].norm();
  
  if (monitoring) {
    cout << "# FGMRES final (true) residual:" << endl;
    cout << "# iteration = " << i << ": |res|/|res0| = " << res/norm0 << endl;
  }

  if (fabs(res - beta) > tol*10.0) {
    cout << "# WARNING in CSysSolve::FlexibleGMRES(): " << endl;
    cout << "# true residual norm and calculated residual norm do not agree." << endl;
    cout << "# res - beta = " << res - beta << endl;
  }
}
