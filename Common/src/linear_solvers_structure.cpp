/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios
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

#include "../include/linear_solvers_structure.hpp"

void CSysSolve::ApplyGivens(const double & s, const double & c, double & h1, double & h2) {
  
  double temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

void CSysSolve::GenerateGivens(double & dx, double & dy, double & s, double & c) {
  
  if ( (dx == 0.0) && (dy == 0.0) ) {
    c = 1.0;
    s = 0.0;
  }
  else if ( fabs(dy) > fabs(dx) ) {
    double tmp = dx/dy;
    dx = sqrt(1.0 + tmp*tmp);
    s = Sign(1.0/dx, dy);
    c = tmp*s;
  }
  else if ( fabs(dy) <= fabs(dx) ) {
    double tmp = dy/dx;
    dy = sqrt(1.0 + tmp*tmp);
    c = Sign(1.0/dy, dx);
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

void CSysSolve::SolveReduced(const int & n, const vector<vector<double> > & Hsbg,
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

void CSysSolve::ModGramSchmidt(int i, vector<vector<double> > & Hsbg, vector<CSysVector> & w) {
  
  bool Convergence = true;
  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  /*--- Parameter for reorthonormalization ---*/
  
  static const double reorth = 0.98;
  
  /*--- Get the norm of the vector being orthogonalized, and find the
  threshold for re-orthogonalization ---*/
  
  double nrm = dotProd(w[i+1], w[i+1]);
  double thr = nrm*reorth;
  
  /*--- The norm of w[i+1] < 0.0 or w[i+1] = NaN ---*/

  if ((nrm <= 0.0) || (nrm != nrm)) Convergence = false;
  
  /*--- Synchronization point to check the convergence of the solver ---*/

#ifdef HAVE_MPI
  
  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
  
  /*--- Convergence criteria ---*/
  
  sbuf_conv[0] = Convergence;
  MPI_Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  
  /*-- Compute global convergence criteria in the master node --*/
  
  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }
  
  MPI_Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
  
  if (sbuf_conv[0] == 1) Convergence = true;
  else Convergence = false;
  
  delete [] sbuf_conv;
  delete [] rbuf_conv;
  
#endif
  
  if (!Convergence) {
    if (rank == MASTER_NODE)
      cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
#ifndef HAVE_MPI
		exit(EXIT_DIVERGENCE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
  }
  
  /*--- Begin main Gram-Schmidt loop ---*/
  
  for (int k = 0; k < i+1; k++) {
    double prod = dotProd(w[i+1], w[k]);
    Hsbg[k][i] = prod;
    w[i+1].Plus_AX(-prod, w[k]);
    
    /*--- Check if reorthogonalization is necessary ---*/
    
    if (prod*prod > thr) {
      prod = dotProd(w[i+1], w[k]);
      Hsbg[k][i] += prod;
      w[i+1].Plus_AX(-prod, w[k]);
    }
    
    /*--- Update the norm and check its size ---*/
    
    nrm -= Hsbg[k][i]*Hsbg[k][i];
    if (nrm < 0.0) nrm = 0.0;
    thr = nrm*reorth;
  }
  
  /*--- Test the resulting vector ---*/
  
  nrm = w[i+1].norm();
  Hsbg[i+1][i] = nrm;
  
//  if (nrm <= 0.0) {
//    
//    /*--- w[i+1] is a linear combination of the w[0:i] ---*/
//    
//    cerr << "The FGMRES linear solver has diverged" << endl;
//#ifndef HAVE_MPI
//    exit(EXIT_DIVERGENCE);
//#else
//    MPI_Abort(MPI_COMM_WORLD,1);
//    MPI_Finalize();
//#endif
//    
//  }
  
  /*--- Scale the resulting vector ---*/
  
  w[i+1] /= nrm;
}

void CSysSolve::WriteHeader(const string & solver, const double & restol, const double & resinit) {
  
  cout << "\n# " << solver << " residual history" << endl;
  cout << "# Residual tolerance target = " << restol << endl;
  cout << "# Initial residual norm     = " << resinit << endl;
  
}

void CSysSolve::WriteHistory(const int & iter, const double & res, const double & resinit) {
  
  cout << "     " << iter << "     " << res/resinit << endl;
  
}

unsigned long CSysSolve::CG_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                           CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
int rank = 0;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Check the subspace size ---*/
  if (m < 1) {
    if (rank == MASTER_NODE) cerr << "CSysSolve::ConjugateGradient: illegal value for subspace size, m = " << m << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
	MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  CSysVector r(b);
  CSysVector A_p(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  mat_vec(x, A_p);
  
  r -= A_p; // recall, r holds b initially
  double norm_r = r.norm();
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
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
  if ((monitoring) && (rank == MASTER_NODE)) {
    WriteHeader("CG", tol, norm_r);
    WriteHistory(i, norm_r, norm0);
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
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 5 == 0)) WriteHistory(i+1, norm_r, norm0);
    
    precond(r, z);
    
    /*--- Calculate Gram-Schmidt coefficient beta,
		 beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/
    beta = 1.0 / r_dot_z;
    r_dot_z = dotProd(r, z);
    beta *= r_dot_z;
    
    /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/
    p.Equals_AX_Plus_BY(beta, p, 1.0, z);
  }
  

  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# Conjugate Gradient final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
//  /*--- Recalculate final residual (this should be optional) ---*/
//  mat_vec(x, A_p);
//  r = b;
//  r -= A_p;
//  double true_res = r.norm();
//  
//  if (fabs(true_res - norm_r) > tol*10.0) {
//    if (rank == MASTER_NODE) {
//      cout << "# WARNING in CSysSolve::ConjugateGradient(): " << endl;
//      cout << "# true residual norm and calculated residual norm do not agree." << endl;
//      cout << "# true_res - calc_res = " << true_res - norm_r << endl;
//    }
//  }
	
	return i;
  
}

unsigned long CSysSolve::FGMRES_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                               CPreconditioner & precond, double tol, unsigned long m, double *residual, bool monitoring) {
	
int rank = 0;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*---  Check the subspace size ---*/
  
  if (m < 1) {
    if (rank == MASTER_NODE) cerr << "CSysSolve::FGMRES: illegal value for subspace size, m = " << m << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
	MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*---  Check the subspace size ---*/
  
  if (m > 1000) {
    if (rank == MASTER_NODE) cerr << "CSysSolve::FGMRES: illegal value for subspace size (too high), m = " << m << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
	MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
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
  
  mat_vec(x, w[0]);
  w[0] -= b;
  
  double beta = w[0].norm();
  
  if ( (beta < tol*norm0) || (beta < eps) ) {
    
    /*---  System is already solved ---*/
    
    if (rank == MASTER_NODE) cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*---  Normalize residual to get w_{0} (the negative sign is because w[0]
	 holds the negative residual, as mentioned above) ---*/
  
  w[0] /= -beta;
  
  /*---  Initialize the RHS of the reduced system ---*/
  
  g[0] = beta;
  
  /*--- Set the norm to the initial residual value ---*/
  
  norm0 = beta;

  /*---  Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    WriteHeader("FGMRES", tol, beta);
    WriteHistory(i, beta, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < m; i++) {
    
    /*---  Check if solution has converged ---*/
    
    if (beta < tol*norm0) break;
    
    /*---  Precondition the CSysVector w[i] and store result in z[i] ---*/
    
    precond(w[i], z[i]);
    
    /*---  Add to Krylov subspace ---*/
    
    mat_vec(z[i], w[i+1]);
    
    /*---  Modified Gram-Schmidt orthogonalization ---*/
    
    ModGramSchmidt(i, H, w);
    
    /*---  Apply old Givens rotations to new column of the Hessenberg matrix
		 then generate the new Givens rotation matrix and apply it to
		 the last two elements of H[:][i] and g ---*/
    
    for (int k = 0; k < i; k++)
      ApplyGivens(sn[k], cs[k], H[k][i], H[k+1][i]);
    GenerateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
    ApplyGivens(sn[i], cs[i], g[i], g[i+1]);
    
    /*---  Set L2 norm of residual and check if solution has converged ---*/
    
    beta = fabs(g[i+1]);
    
    /*---  Output the relative residual if necessary ---*/
    
    if ((((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 50 == 0)) && (rank == MASTER_NODE)) WriteHistory(i+1, beta, norm0);
    
  }

  /*---  Solve the least-squares system and update solution ---*/
  
  SolveReduced(i, H, g, y);
  for (int k = 0; k < i; k++) {
    x.Plus_AX(y[k], z[k]);
  }
  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# FGMRES final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = " << beta/norm0 << ".\n" << endl;
  }
  
//  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
//  mat_vec(x, w[0]);
//  w[0] -= b;
//  double res = w[0].norm();
//
//  if (fabs(res - beta) > tol*10) {
//    if (rank == MASTER_NODE) {
//      cout << "# WARNING in CSysSolve::FGMRES(): " << endl;
//      cout << "# true residual norm and calculated residual norm do not agree." << endl;
//      cout << "# res - beta = " << res - beta << endl;
//    }
//  }
	
  (*residual) = beta;
	return i;
  
}

unsigned long CSysSolve::BCGSTAB_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                 CPreconditioner & precond, double tol, unsigned long m, double *residual, bool monitoring) {
	
  int rank = 0;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Check the subspace size ---*/
  
  if (m < 1) {
    if (rank == MASTER_NODE) cerr << "CSysSolve::BCGSTAB: illegal value for subspace size, m = " << m << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
	MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
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
  
	mat_vec(x, A_x);
  r -= A_x; r_0 = r; // recall, r holds b initially
  double norm_r = r.norm();
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
    return 0;
  }
	
	/*--- Initialization ---*/
  
  double alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
	
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    WriteHeader("BCGSTAB", tol, norm_r);
    WriteHistory(i, norm_r, norm0);
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
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 50 == 0) && (rank == MASTER_NODE)) WriteHistory(i+1, norm_r, norm0);
    
  }
	  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# BCGSTAB final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
	
//  /*--- Recalculate final residual (this should be optional) ---*/
//	mat_vec(x, A_x);
//  r = b; r -= A_x;
//  double true_res = r.norm();
//  
//  if ((fabs(true_res - norm_r) > tol*10.0) && (rank == MASTER_NODE)) {
//    cout << "# WARNING in CSysSolve::BCGSTAB(): " << endl;
//    cout << "# true residual norm and calculated residual norm do not agree." << endl;
//    cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << endl;
//  }
	
  (*residual) = norm_r;
	return i;
}

unsigned long CSysSolve::Solve(CSysMatrix & Jacobian, CSysVector & LinSysRes, CSysVector & LinSysSol, CGeometry *geometry, CConfig *config) {
  
  double SolverTol = config->GetLinear_Solver_Error(), Residual;
  unsigned long MaxIter = config->GetLinear_Solver_Iter();
  unsigned long IterLinSol = 0;
  
  /*--- Solve the linear system using a Krylov subspace method ---*/
  
  if (config->GetKind_Linear_Solver() == BCGSTAB || config->GetKind_Linear_Solver() == FGMRES
      || config->GetKind_Linear_Solver() == RESTARTED_FGMRES) {
    
    CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
    
    CPreconditioner* precond = NULL;
    
    switch (config->GetKind_Linear_Solver_Prec()) {
      case JACOBI:
        Jacobian.BuildJacobiPreconditioner();
        precond = new CJacobiPreconditioner(Jacobian, geometry, config);
        break;
      case ILU:
        Jacobian.BuildILUPreconditioner();
        precond = new CILUPreconditioner(Jacobian, geometry, config);
        break;
      case LU_SGS:
        precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
        break;
      case LINELET:
        Jacobian.BuildJacobiPreconditioner();
        precond = new CLineletPreconditioner(Jacobian, geometry, config);
        break;
      default:
        Jacobian.BuildJacobiPreconditioner();
        precond = new CJacobiPreconditioner(Jacobian, geometry, config);
        break;
    }
    
    switch (config->GetKind_Linear_Solver()) {
      case BCGSTAB:
        IterLinSol = BCGSTAB_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
        break;
      case FGMRES:
        IterLinSol = FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
        break;
      case RESTARTED_FGMRES:
        IterLinSol = 0;
        while (IterLinSol < config->GetLinear_Solver_Iter()) {
          if (IterLinSol + config->GetLinear_Solver_Restart_Frequency() > config->GetLinear_Solver_Iter())
            MaxIter = config->GetLinear_Solver_Iter() - IterLinSol;
          IterLinSol += FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
          if (LinSysRes.norm() < SolverTol) break;
          SolverTol = SolverTol*(1.0/LinSysRes.norm());
        }
        break;
    }
    
    /*--- Dealocate memory of the Krylov subspace method ---*/
    
    delete mat_vec;
    delete precond;
    
  }
  
  /*--- Smooth the linear system. ---*/
  
  else {
    switch (config->GetKind_Linear_Solver()) {
      case SMOOTHER_LUSGS:
        Jacobian.ComputeLU_SGSPreconditioner(LinSysRes, LinSysSol, geometry, config);
        break;
      case SMOOTHER_JACOBI:
        Jacobian.BuildJacobiPreconditioner();
        Jacobian.ComputeJacobiPreconditioner(LinSysRes, LinSysSol, geometry, config);
        break;
      case SMOOTHER_ILU:
        Jacobian.BuildILUPreconditioner();
        Jacobian.ComputeILUPreconditioner(LinSysRes, LinSysSol, geometry, config);
        break;
      case SMOOTHER_LINELET:
        Jacobian.BuildJacobiPreconditioner();
        Jacobian.ComputeLineletPreconditioner(LinSysRes, LinSysSol, geometry, config);
        break;
        IterLinSol = 1;
    }
  }
  
  return IterLinSol;
  
}
