/*!
 * \file linear_solvers_structure.cpp
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

#include "../include/linear_solvers_structure.hpp"

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

void CSysSolve::writeHeader(const string & solver, const double & restol, const double & resinit) {
  
  cout << "# " << solver << " residual history" << endl;
  cout << "# Residual tolerance target = " << restol << endl;
  cout << "# Initial residual norm     = " << resinit << endl;
  
}

void CSysSolve::writeHistory(const int & iter, const double & res, const double & resinit) {
  
  cout << "     " << iter << "     " << res/resinit << endl;
  
}

unsigned long CSysSolve::ConjugateGradient(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                           CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
  int rank = 0;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Check the subspace size ---*/
  if (m < 1) {
    if (rank == 0) cerr << "CSysSolve::ConjugateGradient: illegal value for subspace size, m = " << m << endl;
#ifdef NO_MPI
    exit(1);
#else
    MPI::COMM_WORLD.Abort(1);
    MPI::Finalize();
#endif
  }
  
  CSysVector r(b);
  CSysVector A_p(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  mat_vec(x,A_p);
  
  r -= A_p; // recall, r holds b initially
  double norm_r = r.norm();
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == 0) cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
    return 0;
  }
  
  double alpha, beta, r_dot_z;
  CSysVector z(r);
  precond(r, z);
  CSysVector p(z);
  
  /*--- Set the norm to the initial initial residual value ---*/
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  int i = 0;
  if ((monitoring) && (rank == 0))  {
    writeHeader("CG", tol, norm_r);
    writeHistory(i, norm_r, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  for (i = 0; i < m; i++) {
    
    /*--- Apply matrix to p to build Krylov subspace ---*/
    mat_vec(p, A_p);
    
    /*--- Calculate step-length alpha ---*/
    r_dot_z = dotProd(r, z);
    alpha = dotProd(A_p, p);
    alpha = r_dot_z / alpha;
    
    /*--- Update solution and residual: ---*/
    x.Plus_AX(alpha, p);
    r.Plus_AX(-alpha, A_p);
    
    /*--- Check if solution has converged, else output the relative residual if necessary ---*/
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == 0)) && ((i+1) % 5 == 0)) writeHistory(i+1, norm_r, norm0);
    
    precond(r, z);
    
    /*--- Calculate Gram-Schmidt coefficient beta,
		 beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/
    beta = 1.0 / r_dot_z;
    r_dot_z = dotProd(r, z);
    beta *= r_dot_z;
    
    /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/
    p.Equals_AX_Plus_BY(beta, p, 1.0, z);
  }
  

  
  if ((monitoring) && (rank == 0))  {
    cout << "# Conjugate Gradient final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << endl;
  }
  
//  /*--- Recalculate final residual (this should be optional) ---*/
//  mat_vec(x, A_p);
//  r = b;
//  r -= A_p;
//  double true_res = r.norm();
//  
//  if (fabs(true_res - norm_r) > tol*10.0) {
//    if (rank == 0) {
//      cout << "# WARNING in CSysSolve::ConjugateGradient(): " << endl;
//      cout << "# true residual norm and calculated residual norm do not agree." << endl;
//      cout << "# true_res - calc_res = " << true_res - norm_r << endl;
//    }
//  }
	
	return i;
  
}

unsigned long CSysSolve::FGMRES(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                               CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
  int rank = 0;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*---  Check the subspace size ---*/
  if (m < 1) {
    if (rank == 0) cerr << "CSysSolve::FGMRES: illegal value for subspace size, m = " << m << endl;
#ifdef NO_MPI
    exit(1);
#else
    MPI::COMM_WORLD.Abort(1);
    MPI::Finalize();
#endif
  }

  /*---  Check the subspace size ---*/
  if (m > 1000) {
    if (rank == 0) cerr << "CSysSolve::FGMRES: illegal value for subspace size (too high), m = " << m << endl;
#ifdef NO_MPI
    exit(1);
#else
    MPI::COMM_WORLD.Abort(1);
    MPI::Finalize();
#endif
  }
  
  /*---  Define various arrays
	 Note: elements in w and z are initialized to x to avoid creating
	 a temporary CSysVector object for the copy constructor ---*/
  vector<CSysVector> w(m+1, x);
  vector<CSysVector> z(m+1, x);
  vector<double> g(m+1, 0.0);
  vector<double> sn(m+1, 0.0);
  vector<double> cs(m+1, 0.0);
  vector<double> y(m, 0.0);
  vector<vector<double> > H(m+1, vector<double>(m, 0.0));
  
  /*---  Calculate the norm of the rhs vector ---*/
  double norm0 = b.norm();
  
  /*---  Calculate the initial residual (actually the negative residual)
	 and compute its norm ---*/
  mat_vec(x,w[0]);
  w[0] -= b;
  
  double beta = w[0].norm();
  
  if ( (beta < tol*norm0) || (beta < eps) ) {
    /*---  System is already solved ---*/
    if (rank == 0) cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*---  Normalize residual to get w_{0} (the negative sign is because w[0]
	 holds the negative residual, as mentioned above) ---*/
  w[0] /= -beta;
  
  /*---  Initialize the RHS of the reduced system ---*/
  g[0] = beta;
  
  /*--- Set the norm to the initial initial residual value ---*/
  norm0 = beta;

  /*---  Output header information including initial residual ---*/
  int i = 0;
  if ((monitoring) && (rank == 0)) {
    writeHeader("FGMRES", tol, beta);
    writeHistory(i, beta, norm0);
  }
  
  /*---  Loop over all serach directions ---*/
  for (i = 0; i < m; i++) {
    
    /*---  Check if solution has converged ---*/
    if (beta < tol*norm0) break;
    
    /*---  Precondition the CSysVector w[i] and store result in z[i] ---*/
    precond(w[i], z[i]);
    
    /*---  Add to Krylov subspace ---*/
    mat_vec(z[i], w[i+1]);
    
    /*---  Modified Gram-Schmidt orthogonalization ---*/
    modGramSchmidt(i, H, w);
    
    /*---  Apply old Givens rotations to new column of the Hessenberg matrix
		 then generate the new Givens rotation matrix and apply it to
		 the last two elements of H[:][i] and g ---*/
    for (int k = 0; k < i; k++)
      applyGivens(sn[k], cs[k], H[k][i], H[k+1][i]);
    generateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
    applyGivens(sn[i], cs[i], g[i], g[i+1]);
    
    /*---  Set L2 norm of residual and check if solution has converged ---*/
    beta = fabs(g[i+1]);
    
    /*---  Output the relative residual if necessary ---*/
    if ((((monitoring) && (rank == 0)) && ((i+1) % 5 == 0)) && (rank == 0)) writeHistory(i+1, beta, norm0);
  }

  /*---  Solve the least-squares system and update solution ---*/
  solveReduced(i, H, g, y);
  for (int k = 0; k < i; k++) {
    x.Plus_AX(y[k], z[k]);
  }
  
  if ((monitoring) && (rank == 0)) {
    cout << "# FGMRES final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = " << beta/norm0 << endl;
  }
  
//  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
//  mat_vec(x, w[0]);
//  w[0] -= b;
//  double res = w[0].norm();
//  
//  if (fabs(res - beta) > tol*10) {
//    if (rank == 0) {
//      cout << "# WARNING in CSysSolve::FGMRES(): " << endl;
//      cout << "# true residual norm and calculated residual norm do not agree." << endl;
//      cout << "# res - beta = " << res - beta << endl;
//    }
//  }
	
	return i;
  
}

unsigned long CSysSolve::BCGSTAB(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                 CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
  int rank = 0;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Check the subspace size ---*/
  if (m < 1) {
    if (rank == 0) cerr << "CSysSolve::BCGSTAB: illegal value for subspace size, m = " << m << endl;
#ifdef NO_MPI
    exit(1);
#else
    MPI::COMM_WORLD.Abort(1);
    MPI::Finalize();
#endif
  }
	
  CSysVector r(b);
  CSysVector r_0(b);
  CSysVector p(b);
	CSysVector v(b);
  CSysVector s(b);
	CSysVector t(b);
	CSysVector phat(b);
	CSysVector shat(b);
  CSysVector A_x(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
	mat_vec(x,A_x);
  r -= A_x; r_0 = r; // recall, r holds b initially
  double norm_r = r.norm();
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == 0) cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
    return 0;
  }
	
	/*--- Initialization ---*/
  double alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
	
  /*--- Set the norm to the initial initial residual value ---*/
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  int i = 0;
  if ((monitoring) && (rank == 0)) {
    writeHeader("BCGSTAB", tol, norm_r);
    writeHistory(i, norm_r, norm0);
  }
	
  /*---  Loop over all search directions ---*/
  for (i = 0; i < m; i++) {
		
		/*--- Compute rho_prime ---*/
		rho_prime = rho;
		
		/*--- Compute rho_i ---*/
		rho = dotProd(r, r_0);
		
		/*--- Compute beta ---*/
		beta = (rho / rho_prime) * (alpha /omega);
		
		/*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/
		double beta_omega = -beta*omega;
		p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
		p.Plus_AX(1.0, r);
		
		/*--- Preconditioning step ---*/
		precond(p, phat);
		mat_vec(phat, v);
    
		/*--- Calculate step-length alpha ---*/
    double r_0_v = dotProd(r_0, v);
    alpha = rho / r_0_v;
    
		/*--- s_{i} = r_{i-1} - alpha * v_{i} ---*/
		s.Equals_AX_Plus_BY(1.0, r, -alpha, v);
		
		/*--- Preconditioning step ---*/
		precond(s, shat);
		mat_vec(shat, t);
    
		/*--- Calculate step-length omega ---*/
    omega = dotProd(t, s) / dotProd(t, t);
    
		/*--- Update solution and residual: ---*/
    x.Plus_AX(alpha, phat); x.Plus_AX(omega, shat);
		r.Equals_AX_Plus_BY(1.0, s, -omega, t);
    
    /*--- Check if solution has converged, else output the relative residual if necessary ---*/
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == 0)) && ((i+1) % 5 == 0) && (rank == 0)) writeHistory(i+1, norm_r, norm0);
    
  }
	  
  if ((monitoring) && (rank == 0)) {
    cout << "# BCGSTAB final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << endl;
  }
	
//  /*--- Recalculate final residual (this should be optional) ---*/
//	mat_vec(x, A_x);
//  r = b; r -= A_x;
//  double true_res = r.norm();
//  
//  if ((fabs(true_res - norm_r) > tol*10.0) && (rank == 0)) {
//    cout << "# WARNING in CSysSolve::BCGSTAB(): " << endl;
//    cout << "# true residual norm and calculated residual norm do not agree." << endl;
//    cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << endl;
//  }
	
	return i;
}
