/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios, T. Economon
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

#include "../include/linear_solvers_structure.hpp"
#include "../include/linear_solvers_structure_b.hpp"

template<class CalcType>
void TCSysSolve<CalcType>::ApplyGivens(const CalcType & s, const CalcType & c, CalcType & h1, CalcType & h2) {
  
  CalcType temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

template<class CalcType>
void TCSysSolve<CalcType>::GenerateGivens(CalcType & dx, CalcType & dy, CalcType & s, CalcType & c) {
  
  if ( (dx == 0.0) && (dy == 0.0) ) {
    c = 1.0;
    s = 0.0;
  }
  else if ( fabs(dy) > fabs(dx) ) {
    CalcType tmp = dx/dy;
    dx = sqrt(1.0 + tmp*tmp);
    s = Sign(1.0/dx, dy);
    c = tmp*s;
  }
  else if ( fabs(dy) <= fabs(dx) ) {
    CalcType tmp = dy/dx;
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

template<class CalcType>
void TCSysSolve<CalcType>::SolveReduced(const int & n, const vector<vector<CalcType> > & Hsbg,
                             const vector<CalcType> & rhs, vector<CalcType> & x) {
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

template<class CalcType>
void TCSysSolve<CalcType>::ModGramSchmidt(int i, vector<vector<CalcType> > & Hsbg, vector<TCSysVector<CalcType> > & w) {
  
  bool Convergence = true;
  int rank = SU2_MPI::GetRank();
  int size = SU2_MPI::GetSize();
  
  /*--- Parameter for reorthonormalization ---*/
  
  static const CalcType reorth = 0.98;
  
  /*--- Get the norm of the vector being orthogonalized, and find the
  threshold for re-orthogonalization ---*/
  
  CalcType nrm = dotProd(w[i+1], w[i+1]);
  CalcType thr = nrm*reorth;
  
  /*--- The norm of w[i+1] < 0.0 or w[i+1] = NaN ---*/

  if ((nrm <= 0.0) || (nrm != nrm)) Convergence = false;
  
  /*--- Synchronization point to check the convergence of the solver ---*/

#ifdef HAVE_MPI
  
  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
  
  /*--- Convergence criteria ---*/
  
  sbuf_conv[0] = Convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  
  /*-- Compute global convergence criteria in the master node --*/
  
  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }
  
  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
  
  if (sbuf_conv[0] == 1) Convergence = true;
  else Convergence = false;
  
  delete [] sbuf_conv;
  delete [] rbuf_conv;
  
#endif
  
  if (!Convergence) {
    SU2_MPI::Error("SU2 has diverged.", CURRENT_FUNCTION);
  }
  
  /*--- Begin main Gram-Schmidt loop ---*/
  
  for (int k = 0; k < i+1; k++) {
    CalcType prod = dotProd(w[i+1], w[k]);
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

  /*--- Scale the resulting vector ---*/
  
  w[i+1] /= nrm;

}

template<class CalcType>
void TCSysSolve<CalcType>::WriteHeader(const string & solver, const CalcType & restol, const CalcType & resinit) {
  
  cout << "\n# " << solver << " residual history" << endl;
  cout << "# Residual tolerance target = " << restol << endl;
  cout << "# Initial residual norm     = " << resinit << endl;
  
}

template<class CalcType>
void TCSysSolve<CalcType>::WriteHistory(const int & iter, const CalcType & res, const CalcType & resinit) {
  
  cout << "     " << iter << "     " << res/resinit << endl;
  
}

template<class CalcType>
unsigned long TCSysSolve<CalcType>::CG_LinSolver(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec,
                                           TCPreconditioner<CalcType> & precond, CalcType tol, unsigned long m, CalcType *residual, bool monitoring) {

  int rank = SU2_MPI::GetRank();

  /*--- Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  TCSysVector<CalcType> r(b);
  TCSysVector<CalcType> A_p(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  
  mat_vec(x, A_p);
  
  r -= A_p; // recall, r holds b initially
  CalcType norm_r = r.norm();
  CalcType norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysSolve<CalcType>::ConjugateGradient(): system solved by initial guess." << endl;
    return 0;
  }
  
  CalcType alpha, beta, r_dot_z;
  TCSysVector<CalcType> z(r);
  precond(r, z);
  TCSysVector<CalcType> p(z);
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    WriteHeader("CG", tol, norm_r);
    WriteHistory(i, norm_r, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)m; i++) {
    
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
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0)) WriteHistory(i+1, norm_r, norm0);
    
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
  
  /*--- Recalculate final residual (this should be optional) ---*/
  
  if (monitoring) {
    
    mat_vec(x, A_p);
    r = b;
    r -= A_p;
    CalcType true_res = r.norm();
    
    if (fabs(true_res - norm_r) > tol*10.0) {
      if (rank == MASTER_NODE) {
        cout << "# WARNING in TCSysSolve<CalcType>::CG_LinSolver(): " << endl;
        cout << "# true residual norm and calculated residual norm do not agree." << endl;
        cout << "# true_res = " << true_res <<", calc_res = " << norm_r <<", tol = " << tol*10 <<"."<< endl;
        cout << "# true_res - calc_res = " << true_res - norm_r << endl;
      }
    }
    
  }

  
  (*residual) = norm_r;
	return (unsigned long) i;
  
}

template<class CalcType>
unsigned long TCSysSolve<CalcType>::FGMRES_LinSolver(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec,
                               TCPreconditioner<CalcType> & precond, CalcType tol, unsigned long m, CalcType *residual, bool monitoring) {
	
  int rank = SU2_MPI::GetRank();
  
  /*---  Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  /*---  Check the subspace size ---*/
  
  if (m > 1000) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size (too high), m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  /*---  Define various arrays
	 Note: elements in w and z are initialized to x to avoid creating
	 a temporary TCSysVector object for the copy constructor ---*/
  
  vector<TCSysVector<CalcType> > w(m+1, x);
  vector<TCSysVector<CalcType >> z(m+1, x);
  vector<CalcType> g(m+1, 0.0);
  vector<CalcType> sn(m+1, 0.0);
  vector<CalcType> cs(m+1, 0.0);
  vector<CalcType> y(m, 0.0);
  vector<vector<CalcType> > H(m+1, vector<CalcType>(m, 0.0));
  
  /*---  Calculate the norm of the rhs vector ---*/
  
  CalcType norm0 = b.norm();
  
  /*---  Calculate the initial residual (actually the negative residual)
	 and compute its norm ---*/
  
  mat_vec(x, w[0]);
  w[0] -= b;
  
  CalcType beta = w[0].norm();
  
  if ( (beta < tol*norm0) || (beta < eps) ) {
    
    /*---  System is already solved ---*/
    
    if (rank == MASTER_NODE) cout << "TCSysSolve<CalcType>::FGMRES(): system solved by initial guess." << endl;
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
  
  for (i = 0; i < (int)m; i++) {
    
    /*---  Check if solution has converged ---*/
    
    if (beta < tol*norm0) break;
    
    /*---  Precondition the TCSysVector w[i] and store result in z[i] ---*/
    
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
    
    if ((((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0)) && (rank == MASTER_NODE)) WriteHistory(i+1, beta, norm0);
    
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
  
  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
  
  if (monitoring) {
    mat_vec(x, w[0]);
    w[0] -= b;
    CalcType res = w[0].norm();
    
    if (fabs(res - beta) > tol*10) {
      if (rank == MASTER_NODE) {
        cout << "# WARNING in TCSysSolve<CalcType>::FGMRES_LinSolver(): " << endl;
        cout << "# true residual norm and calculated residual norm do not agree." << endl;
        cout << "# res = " << res <<", beta = " << beta <<", tol = " << tol*10 <<"."<< endl;
        cout << "# res - beta = " << res - beta << endl << endl;
      }
    }
  }
  
  (*residual) = beta;
  return (unsigned long) i;
  
}

template<class CalcType>
unsigned long TCSysSolve<CalcType>::BCGSTAB_LinSolver(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec,
                                           TCPreconditioner<CalcType> & precond, CalcType tol, unsigned long m, CalcType *residual, bool monitoring) {
  
  int rank = SU2_MPI::GetRank();
  
  /*--- Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  TCSysVector<CalcType> r(b);
  TCSysVector<CalcType> r_0(b);
  TCSysVector<CalcType> p(b);
  TCSysVector<CalcType> v(b);
  TCSysVector<CalcType> s(b);
  TCSysVector<CalcType> t(b);
  TCSysVector<CalcType> phat(b);
  TCSysVector<CalcType> shat(b);
  TCSysVector<CalcType> A_x(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  
  mat_vec(x, A_x);
  r -= A_x; r_0 = r; // recall, r holds b initially
  CalcType norm_r = r.norm();
  CalcType norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysSolve<CalcType>::BCGSTAB(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*--- Initialization ---*/
  
  CalcType alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == MASTER_NODE)) {
    WriteHeader("BCGSTAB", tol, norm_r);
    WriteHistory(i, norm_r, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)m; i++) {
    
    /*--- Compute rho_prime ---*/
    
    rho_prime = rho;
    
    /*--- Compute rho_i ---*/
    
    rho = dotProd(r, r_0);
    
    /*--- Compute beta ---*/
    
    beta = (rho / rho_prime) * (alpha /omega);
    
    /*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/
    
    CalcType beta_omega = -beta*omega;
    p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
    p.Plus_AX(1.0, r);
    
    /*--- Preconditioning step ---*/
    
    precond(p, phat);
    mat_vec(phat, v);
    
    /*--- Calculate step-length alpha ---*/
    
    CalcType r_0_v = dotProd(r_0, v);
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
    if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0) && (rank == MASTER_NODE)) WriteHistory(i+1, norm_r, norm0);
    
  }
  
  if ((monitoring) && (rank == MASTER_NODE)) {
    cout << "# BCGSTAB final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
    /*--- Recalculate final residual (this should be optional) ---*/
  if (monitoring) {
    mat_vec(x, A_x);
    r = b; r -= A_x;
    CalcType true_res = r.norm();
    
    if ((fabs(true_res - norm_r) > tol*10.0) && (rank == MASTER_NODE)) {
      cout << "# WARNING in TCSysSolve<CalcType>::BCGSTAB_LinSolver(): " << endl;
      cout << "# true residual norm and calculated residual norm do not agree." << endl;
      cout << "# true_res = " << true_res <<", calc_res = " << norm_r <<", tol = " << tol*10 <<"."<< endl;
      cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << endl;
    }
  }
  
  (*residual) = norm_r;
  return (unsigned long) i;
}

template<>
unsigned long TCSysSolve<su2double>::Solve(TCSysMatrix<su2double> & Jacobian, TCSysVector<su2double> & LinSysRes, TCSysVector<su2double> & LinSysSol, CGeometry *geometry, CConfig *config) {
  
  su2double SolverTol = config->GetLinear_Solver_Error(), Residual;
  unsigned long MaxIter = config->GetLinear_Solver_Iter();
  unsigned long IterLinSol = 0;
  TCMatrixVectorProduct<su2double> *mat_vec;

  bool TapeActive = NO;
  
  if (config->GetDiscrete_Adjoint()) {
    
    
    TapeActive = AD::globalTape.isActive();

    /*--- Initialize the external function helper with storage of primal input and output disabled  ---*/
    
    AD::InitExtFunc(false, false);
    
    /*--- Set the right-hand side of the linear system as input of the external function ---*/
    
    AD::SetExtFuncIn(LinSysRes.vec_val, LinSysRes.GetLocSize());


    AD::StopRecording();
  }

  /*--- Solve the linear system using a Krylov subspace method ---*/
  
  if (config->GetKind_Linear_Solver() == BCGSTAB ||
      config->GetKind_Linear_Solver() == FGMRES ||
      config->GetKind_Linear_Solver() == RESTARTED_FGMRES ||
      config->GetKind_Linear_Solver() == CONJUGATE_GRADIENT) {
    
    mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
    TCPreconditioner<su2double>* precond = NULL;
    
    switch (config->GetKind_Linear_Solver_Prec()) {
      case JACOBI:
        Jacobian.BuildJacobiPreconditioner();
        precond = new TCJacobiPreconditioner<su2double>(Jacobian, geometry, config);
        break;
      case ILU:
        Jacobian.BuildILUPreconditioner();
        precond = new TCILUPreconditioner<su2double>(Jacobian, geometry, config);
        break;
      case LU_SGS:
        precond = new TCLU_SGSPreconditioner<su2double>(Jacobian, geometry, config);
        break;
      case LINELET:
        Jacobian.BuildJacobiPreconditioner();
        precond = new TCLineletPreconditioner<su2double>(Jacobian, geometry, config);
        break;
      default:
        Jacobian.BuildJacobiPreconditioner();
        precond = new TCJacobiPreconditioner<su2double>(Jacobian, geometry, config);
        break;
    }
    
    switch (config->GetKind_Linear_Solver()) {
      case BCGSTAB:
        IterLinSol = BCGSTAB_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
        break;
      case FGMRES:
        IterLinSol = FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
        break;
      case CONJUGATE_GRADIENT:
        IterLinSol = CG_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
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
        mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
        IterLinSol = Jacobian.LU_SGS_Smoother(LinSysRes, LinSysSol, *mat_vec, SolverTol, MaxIter, &Residual, false, geometry, config);
        delete mat_vec;
        break;
      case SMOOTHER_JACOBI:
        mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
        Jacobian.BuildJacobiPreconditioner();
        IterLinSol = Jacobian.Jacobi_Smoother(LinSysRes, LinSysSol, *mat_vec, SolverTol, MaxIter, &Residual, false, geometry, config);
        delete mat_vec;
        break;
      case SMOOTHER_ILU:
        mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
        Jacobian.BuildILUPreconditioner();
        IterLinSol = Jacobian.ILU_Smoother(LinSysRes, LinSysSol, *mat_vec, SolverTol, MaxIter, &Residual, false, geometry, config);
        delete mat_vec;
        break;
      case SMOOTHER_LINELET:
        Jacobian.BuildJacobiPreconditioner();
        Jacobian.ComputeLineletPreconditioner(LinSysRes, LinSysSol, geometry, config);
        IterLinSol = 1;
        break;
    }
  }


  if(TapeActive) {
    
    AD::StartRecording();

    AD::SetExtFuncOut(LinSysSol.vec_val, LinSysSol.GetLocSize());
    
    TCSysVector<passivedouble>* LinSysRes_b = new TCSysVector<passivedouble>(LinSysRes.GetNBlk(), LinSysRes.GetNBlkDomain(), LinSysRes.GetNVar(), 0.0);
    TCSysVector<passivedouble>* LinSysSol_b = new TCSysVector<passivedouble>(LinSysSol.GetNBlk(), LinSysSol.GetNBlkDomain(), LinSysSol.GetNVar(), 0.0);
    
    TCSysMatrix<passivedouble>* Jacobian_b = new TCSysMatrix<passivedouble>();
    Jacobian_b->Initialize(geometry->GetnPoint(), geometry->GetnPointDomain(), LinSysRes.GetNVar(), LinSysRes.GetNVar(), true, geometry, config );
    unsigned short nVar = LinSysRes.GetNVar();
    
    unsigned long iPoint, jPoint;
    su2double* Block_i, *Block_j, *Block_ii, *Block_jj;
    passivedouble* Block_i_p, *Block_j_p, *Block_ii_p, *Block_jj_p;
    Block_i_p = new passivedouble[nVar*nVar];
    Block_j_p = new passivedouble[nVar*nVar];
    Block_ii_p = new passivedouble[nVar*nVar];
    Block_jj_p = new passivedouble[nVar*nVar];
    
    for (unsigned long iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
      
      Block_i = Jacobian.GetBlock(iPoint,jPoint);
      Block_j = Jacobian.GetBlock(jPoint, iPoint);
      Block_ii = Jacobian.GetBlock(iPoint, iPoint);
      Block_jj = Jacobian.GetBlock(jPoint, jPoint);
      
      for (unsigned short iVar = 0; iVar < nVar; iVar++){
        for (unsigned short jVar = 0; jVar < nVar; jVar++){
          
          Block_i_p[iVar*nVar+jVar] = SU2_TYPE::GetValue(Block_i[iVar*nVar+jVar]);
          Block_j_p[iVar*nVar+jVar] = SU2_TYPE::GetValue(Block_j[iVar*nVar+jVar]);
          Block_ii_p[iVar*nVar+jVar] = SU2_TYPE::GetValue(Block_ii[iVar*nVar+jVar]);
          Block_jj_p[iVar*nVar+jVar] = SU2_TYPE::GetValue(Block_jj[iVar*nVar+jVar]);
          
        }
      }      
      Jacobian_b->SetBlock(iPoint, jPoint, Block_i_p);
      Jacobian_b->SetBlock(jPoint, iPoint, Block_j_p);
      Jacobian_b->SetBlock(iPoint, iPoint, Block_ii_p);
      Jacobian_b->SetBlock(iPoint, iPoint, Block_jj_p);
    }
    
#ifdef CODI_REVERSE_TYPE
    AD::FuncHelper->addUserData(LinSysRes_b);
    AD::FuncHelper->addUserData(LinSysSol_b);
    AD::FuncHelper->addUserData(Jacobian_b);
    AD::FuncHelper->addUserData(geometry);
    AD::FuncHelper->addUserData(config);
    
    AD::FuncHelper->addToTape(CSysSolve_b::Solve_b);
    
    delete AD::FuncHelper;
#endif

  }

  return IterLinSol;
  
}

template<>
unsigned long TCSysSolve<passivedouble>::Solve(TCSysMatrix<passivedouble> & Jacobian, TCSysVector<passivedouble> & LinSysRes, TCSysVector<passivedouble> & LinSysSol, CGeometry *geometry, CConfig *config) {
    
  }

template<class CalcType>
void TCSysSolve<CalcType>::SetExternalSolve(TCSysMatrix<CalcType> & Jacobian, TCSysVector<CalcType> & LinSysRes, TCSysVector<CalcType> & LinSysSol, CGeometry *geometry, CConfig *config) {

//#ifdef CODI_REVERSE_TYPE
  
//  unsigned long size = LinSysRes.GetLocSize();
//  unsigned long i, nBlk = LinSysRes.GetNBlk(),
//                nVar = LinSysRes.GetNVar(),
//                nBlkDomain = LinSysRes.GetNBlkDomain();

//  /*--- Arrays to store the indices of the input/output of the linear solver.
//     * Note: They will be deleted in the TCSysSolve<CalcType>_b::Delete_b routine. ---*/

//  CalcType::GradientData *LinSysRes_Indices = new CalcType::GradientData[size];
//  CalcType::GradientData *LinSysSol_Indices = new CalcType::GradientData[size];
//#if CODI_PRIMAL_INDEX_TAPE
//  CalcType::Real *oldValues = new CalcType::Real[size];
//#endif

//  for (i = 0; i < size; i++) {

//    /*--- Register the solution of the linear system (could already be registered when using multigrid) ---*/

//    if (!LinSysSol[i].isActive()) {
//#if CODI_PRIMAL_INDEX_TAPE
//      oldValues[i] = AD::globalTape.registerExtFunctionOutput(LinSysSol[i]);
//#else
//      AD::globalTape.registerInput(LinSysSol[i]);
//#endif
//    }

//    /*--- Store the indices ---*/

//    LinSysRes_Indices[i] = LinSysRes[i].getGradientData();
//    LinSysSol_Indices[i] = LinSysSol[i].getGradientData();
//  }

//  /*--- Push the data to the checkpoint handler for access in the reverse sweep ---*/

//  AD::CheckpointHandler* dataHandler = new AD::CheckpointHandler;

//  dataHandler->addData(LinSysRes_Indices);
//  dataHandler->addData(LinSysSol_Indices);
//#if CODI_PRIMAL_INDEX_TAPE
//  dataHandler->addData(oldValues);
//#endif
//  dataHandler->addData(size);
//  dataHandler->addData(nBlk);
//  dataHandler->addData(nVar);
//  dataHandler->addData(nBlkDomain);
//  dataHandler->addData(&Jacobian);
//  dataHandler->addData(geometry);
//  dataHandler->addData(config);

//  /*--- Build preconditioner for the transposed Jacobian ---*/

//  switch(config->GetKind_DiscAdj_Linear_Prec()) {
//    case ILU:
//      Jacobian.BuildILUPreconditioner(true);
//      break;
//    case JACOBI:
//      Jacobian.BuildJacobiPreconditioner(true);
//      break;
//    default:
//      SU2_MPI::Error("The specified preconditioner is not yet implemented for the discrete adjoint method.", CURRENT_FUNCTION);
//      break;
//  }

//  /*--- Push the external function to the AD tape ---*/

//  AD::globalTape.pushExternalFunctionHandle(&TCSysSolve<CalcType>_b::Solve_b, dataHandler, &TCSysSolve<CalcType>_b::Delete_b);

//#endif
}

template<class CalcType>
void TCSysSolve<CalcType>::SetExternalSolve_Mesh(TCSysMatrix<CalcType> & Jacobian, TCSysVector<CalcType> & LinSysRes, TCSysVector<CalcType> & LinSysSol, CGeometry *geometry, CConfig *config){

//#ifdef CODI_REVERSE_TYPE

//  unsigned long size = LinSysRes.GetLocSize();
//  unsigned long i, nBlk = LinSysRes.GetNBlk(),
//                nVar = LinSysRes.GetNVar(),
//                nBlkDomain = LinSysRes.GetNBlkDomain();

//  /*--- Arrays to store the indices of the input/output of the linear solver.
//     * Note: They will be deleted in the TCSysSolve<CalcType>_b::Delete_b routine. ---*/

//  CalcType::GradientData *LinSysRes_Indices = new CalcType::GradientData[size];
//  CalcType::GradientData *LinSysSol_Indices = new CalcType::GradientData[size];

//  for (i = 0; i < size; i++){

//    /*--- Register the solution of the linear system (could already be registered when using multigrid) ---*/

//    if (!LinSysSol[i].isActive()){
//      AD::globalTape.registerInput(LinSysSol[i]);
//    }

//    /*--- Store the indices ---*/

//    LinSysRes_Indices[i] = LinSysRes[i].getGradientData();
//    LinSysSol_Indices[i] = LinSysSol[i].getGradientData();
//  }

//  /*--- Push the data to the checkpoint handler for access in the reverse sweep ---*/

//  AD::CheckpointHandler* dataHandler = new AD::CheckpointHandler;

//  dataHandler->addData(LinSysRes_Indices);
//  dataHandler->addData(LinSysSol_Indices);
//  dataHandler->addData(size);
//  dataHandler->addData(nBlk);
//  dataHandler->addData(nVar);
//  dataHandler->addData(nBlkDomain);
//  dataHandler->addData(&Jacobian);
//  dataHandler->addData(geometry);
//  dataHandler->addData(config);

//  /*--- Build preconditioner for the transposed Jacobian ---*/

//  switch(config->GetKind_DiscAdj_Linear_Prec()){
//    case ILU:
//      Jacobian.BuildILUPreconditioner(false);
//      break;
//    case JACOBI:
//      Jacobian.BuildJacobiPreconditioner(false);
//      break;
//    default:
//      SU2_MPI::Error("The specified preconditioner is not yet implemented for the discrete adjoint method.", CURRENT_FUNCTION);
//      break;
//  }

//  /*--- Push the external function to the AD tape ---*/

//  AD::globalTape.pushExternalFunctionHandle(&TCSysSolve<CalcType>_b::Solve_g, dataHandler, &TCSysSolve<CalcType>_b::Delete_b);

//#endif
}
