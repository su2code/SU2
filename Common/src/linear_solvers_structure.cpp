/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios, T. Economon
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

#include "../include/linear_solvers_structure.hpp"
#include "../include/linear_solvers_structure_b.hpp"

template<class CalcType, class BaseType>
Convert<CalcType, BaseType>::Convert(){
  

}

template<class CalcType, class BaseType>
Convert<CalcType, BaseType>::~Convert(){
  
}

template<class CalcType, class BaseType>
void Convert<CalcType, BaseType>::Initialize(unsigned short MaxSize){

  Block_CalcType = new CalcType*[MaxSize];
  Block_BaseType = new BaseType*[MaxSize];

  for (unsigned short i=0; i < MaxSize; i++){
    Block_CalcType[i] = new CalcType[MaxSize];
    Block_BaseType[i] = new BaseType[MaxSize];
  } 
  BlockLin_CalcType = new CalcType[MaxSize];
  BlockLin_BaseType = new BaseType[MaxSize];
  
}


template<class CalcType>
TCLinSolver<CalcType>::TCLinSolver(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring){
  
  SubSpaceSize = m;
  Output       = monitoring;
  Tolerance    = tol;
  nElem        = elem;
  nElemDomain  = elemdomain;
  Residual     = 0.0;
  BlockSize = blocksize;
  
}

template<class CalcType>
TCLinSolver<CalcType>::~TCLinSolver(){}

template<class CalcType>
void TCLinSolver<CalcType>::WriteHeader(const string & solver, const CalcType & restol, const CalcType & resinit) {
  
  cout << "\n# " << solver << " residual history" << endl;
  cout << "# Residual TCLinSolverBase::Tolerance target = " << restol << endl;
  cout << "# Initial residual norm     = " << resinit << endl;
  
}

template<class CalcType>
void TCLinSolver<CalcType>::WriteHistory(const int & iter, const CalcType & res, const CalcType & resinit) {
  
  cout << "     " << iter << "     " << res/resinit << endl;
  
}

template<class CalcType>
unsigned long TCLinSolver<CalcType>::Solve(const TCSysVector<CalcType> &b, TCSysVector<CalcType> &x, TCMatrixVectorProduct<CalcType> &mat_vec, TCPreconditioner<CalcType> &precond){
  
  return 0;
}

template<class CalcType>
TCLinSolver_CG<CalcType>::TCLinSolver_CG(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring) : TCLinSolver<CalcType>(blocksize, elem, elemdomain, tol, m, monitoring){
    
  /*--- Check the subspace size ---*/
  
  if (TCLinSolverBase::SubSpaceSize < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", TCLinSolverBase::SubSpaceSize );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  /*--- Initialize vector structures ---*/
  
  r.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  A_p.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  z.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  p.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  
  /*--- Initialize variables ---*/
  
  norm_r = 0.0;
  norm0  = 0.0;
  alpha  = 0.0;
  beta   = 0.0;
  r_dot_z = 0.0;
  
}

template<class CalcType>
TCLinSolver_CG<CalcType>::~TCLinSolver_CG(){}

template<class CalcType>
unsigned long TCLinSolver_CG<CalcType>::Solve(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec,
                                              TCPreconditioner<CalcType> & precond) {

  int rank = SU2_MPI::GetRank();

  r   = b;
  A_p = b;
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  
  mat_vec(x, A_p);
  
  r -= A_p; // recall, r holds b initially
  norm_r = r.norm();
  norm0 = b.norm();
  if ( (norm_r < TCLinSolverBase::Tolerance*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysSolve<CalcType>::ConjugateGradient(): system solved by initial guess." << endl;
    return 0;
  }
  
  z = r;
  precond(r, z);
  p = z;
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- TCLinSolverBase::Output header information including initial residual ---*/
  
  int i = 0;
  if ((TCLinSolverBase::Output) && (rank == MASTER_NODE)) {
    TCLinSolverBase::WriteHeader("CG", TCLinSolverBase::Tolerance, norm_r);
    TCLinSolverBase::WriteHistory(i, norm_r, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)TCLinSolverBase::SubSpaceSize; i++) {
    
    /*--- Apply matrix to p to build Krylov subspace ---*/
    
    mat_vec(p, A_p);
    
    /*--- Calculate step-length alpha ---*/
    
    r_dot_z = dotProd(r, z);
    alpha = dotProd(A_p, p);
    alpha = r_dot_z / alpha;
    
    /*--- Update solution and residual: ---*/
    
    x.Plus_AX(alpha, p);
    r.Plus_AX(-alpha, A_p);
    
    /*--- Check if solution has converged, else TCLinSolverBase::Output the relative residual if necessary ---*/
    
    norm_r = r.norm();
    if (norm_r < TCLinSolverBase::Tolerance*norm0) break;
    if (((TCLinSolverBase::Output) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0)) TCLinSolverBase::WriteHistory(i+1, norm_r, norm0);
    
    precond(r, z);
    
    /*--- Calculate Gram-Schmidt coefficient beta,
		 beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/
    
    beta = 1.0 / r_dot_z;
    r_dot_z = dotProd(r, z);
    beta *= r_dot_z;
    
    /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/
    
    p.Equals_AX_Plus_BY(beta, p, 1.0, z);
    
  }
  
  if ((TCLinSolverBase::Output) && (rank == MASTER_NODE)) {
    cout << "# Conjugate Gradient final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
  /*--- Recalculate final residual (this should be optional) ---*/
  
  if (TCLinSolverBase::Output) {
    
    mat_vec(x, A_p);
    r = b;
    r -= A_p;
    CalcType true_res = r.norm();
    
    if (fabs(true_res - norm_r) > TCLinSolverBase::Tolerance*10.0) {
      if (rank == MASTER_NODE) {
        cout << "# WARNING in TCSysSolve<CalcType>::CG_LinSolver(): " << endl;
        cout << "# true residual norm and calculated residual norm do not agree." << endl;
        cout << "# true_res = " << true_res <<", calc_res = " << norm_r <<", tol = " << TCLinSolverBase::Tolerance*10 <<"."<< endl;
        cout << "# true_res - calc_res = " << true_res - norm_r << endl;
      }
    }
    
  }

  
  TCLinSolverBase::Residual = norm_r;
	return (unsigned long) i;
  
}

template<class CalcType>
TCLinSolver_FGMRES<CalcType>::TCLinSolver_FGMRES(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring) : TCLinSolver<CalcType>(blocksize, elem, elemdomain, tol, m, monitoring){
    
  /*---  Check the subspace size ---*/
  
  if (TCLinSolverBase::SubSpaceSize < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", TCLinSolverBase::SubSpaceSize );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  /*---  Check the subspace size ---*/
  
  if (TCLinSolverBase::SubSpaceSize > 5000) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size (too high), m = %lu", TCLinSolverBase::SubSpaceSize );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  /*--- Initialize arrays ---*/
  
  w.resize(TCLinSolverBase::SubSpaceSize+1, TCSysVector<CalcType>(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain,TCLinSolverBase::BlockSize, 0.0));
  z.resize(TCLinSolverBase::SubSpaceSize,   TCSysVector<CalcType>(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain,TCLinSolverBase::BlockSize, 0.0));
  g.resize(TCLinSolverBase::SubSpaceSize+1, 0.0);
  sn.resize(TCLinSolverBase::SubSpaceSize+1, 0.0);
  cs.resize(TCLinSolverBase::SubSpaceSize+1, 0.0);
  y.resize(TCLinSolverBase::SubSpaceSize, 0.0);
  H.resize(TCLinSolverBase::SubSpaceSize+1, vector<CalcType>(TCLinSolverBase::SubSpaceSize, 0.0));
  
  /*--- Initialize variables ---*/
  
  norm0 = 0.0;
  beta  = 0.0;
}

template<class CalcType>
TCLinSolver_FGMRES<CalcType>::~TCLinSolver_FGMRES(){}


template<class CalcType>
unsigned long TCLinSolver_FGMRES<CalcType>::Solve(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec,
                                                  TCPreconditioner<CalcType> & precond) {
	
  int rank = SU2_MPI::GetRank();
  int i = 0;
 
  /*--- Initialize all vectors with 0 ---*/
  
  std::fill(g.begin(), g.end(), 0.0);
  std::fill(sn.begin(), sn.end(), 0.0);
  std::fill(cs.begin(), cs.end(), 0.0);
  std::fill(y.begin(), y.end(), 0.0);

  /*---  Calculate the norm of the rhs vector ---*/
  
  norm0 = b.norm();
  
  /*---  Calculate the initial residual (actually the negative residual)
	 and compute its norm ---*/
  
  mat_vec(x, w[0]);
  w[0] -= b;
  
  beta = w[0].norm();
  
  if ( (beta < TCLinSolverBase::Tolerance*norm0) || (beta < eps) ) {
    
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

  /*---  TCLinSolverBase::Output header information including initial residual ---*/
  
  if ((TCLinSolverBase::Output) && (rank == MASTER_NODE)) {
    TCLinSolverBase::WriteHeader("FGMRES", TCLinSolverBase::Tolerance, beta);
    TCLinSolverBase::WriteHistory(i, beta, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)TCLinSolverBase::SubSpaceSize; i++) {
    
    /*---  Check if solution has converged ---*/
    
    if (beta < TCLinSolverBase::Tolerance*norm0) break;
    
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
    
    /*---  TCLinSolverBase::Output the relative residual if necessary ---*/
    
    if ((((TCLinSolverBase::Output) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0)) && (rank == MASTER_NODE)) TCLinSolverBase::WriteHistory(i+1, beta, norm0);
    
  }

  /*---  Solve the least-squares system and update solution ---*/
  
  SolveReduced(i, H, g, y);
  for (int k = 0; k < i; k++) {
    x.Plus_AX(y[k], z[k]);
  }
  
  if ((TCLinSolverBase::Output) && (rank == MASTER_NODE)) {
    cout << "# FGMRES final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = " << beta/norm0 << ".\n" << endl;
  }
  
  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
  
  if (TCLinSolverBase::Output) {
    mat_vec(x, w[0]);
    w[0] -= b;
    CalcType res = w[0].norm();
    
    if (fabs(res - beta) > TCLinSolverBase::Tolerance*10) {
      if (rank == MASTER_NODE) {
        cout << "# WARNING in TCSysSolve<CalcType>::FGMRES_LinSolver(): " << endl;
        cout << "# true residual norm and calculated residual norm do not agree." << endl;
        cout << "# res = " << res <<", beta = " << beta <<", tol = " << TCLinSolverBase::Tolerance*10 <<"."<< endl;
        cout << "# res - beta = " << res - beta << endl << endl;
      }
    }
  }
  
 TCLinSolverBase::Residual = beta;
  return (unsigned long) i;
  
}


template<class CalcType>
void TCLinSolver_FGMRES<CalcType>::ApplyGivens(const CalcType & s, const CalcType & c, CalcType & h1, CalcType & h2) {
  
  CalcType temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

template<class CalcType>
void TCLinSolver_FGMRES<CalcType>::GenerateGivens(CalcType & dx, CalcType & dy, CalcType & s, CalcType & c) {
  
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
void TCLinSolver_FGMRES<CalcType>::SolveReduced(const int & n, const vector<vector<CalcType> > & Hsbg,
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
void TCLinSolver_FGMRES<CalcType>::ModGramSchmidt(int i, vector<vector<CalcType> > & Hsbg, vector<TCSysVector<CalcType> > & w) {
  
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
    SU2_MPI::Error("Divergence of linear solver.", CURRENT_FUNCTION);
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
TCLinSolver_BCGSTAB<CalcType>::TCLinSolver_BCGSTAB(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring): TCLinSolver<CalcType>(blocksize, elem, elemdomain, tol, m, monitoring){
  
  /*--- Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
  
  r.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  r_0.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  p.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  v.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  s.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  t.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  phat.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  shat.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  A_x.Initialize(TCLinSolverBase::nElem, TCLinSolverBase::nElemDomain, TCLinSolverBase::BlockSize, 0.0);
  
  norm_r = 0.0;
  norm0 = 0.0;
  alpha = 0.0;
  beta = 0.0;
  omega = 0.0;
  rho = 0.0;
  rho_prime = 0.0;
  beta_omega = 0.0;
  r_0_v = 0.0;

}

template<class CalcType>
TCLinSolver_BCGSTAB<CalcType>::~TCLinSolver_BCGSTAB(){}

template<class CalcType>
unsigned long TCLinSolver_BCGSTAB<CalcType>::Solve(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec,
                                                   TCPreconditioner<CalcType> & precond) {
  
  int rank = SU2_MPI::GetRank();
  
  r = b;
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  
  mat_vec(x, A_x);
  r -= A_x; r_0 = r; // recall, r holds b initially
  norm_r = r.norm();
  norm0 = b.norm();
  if ( (norm_r < TCLinSolverBase::Tolerance*norm0) || (norm_r < eps) ) {
    if (rank == MASTER_NODE) cout << "TCSysSolve<CalcType>::BCGSTAB(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*--- Initialization ---*/
  
  alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
  
  /*--- Set the norm to the initial initial residual value ---*/
  
  norm0 = norm_r;
  
  /*--- TCLinSolverBase::Output header information including initial residual ---*/
  
  int i = 0;
  if ((TCLinSolverBase::Output) && (rank == MASTER_NODE)) {
    TCLinSolverBase::WriteHeader("BCGSTAB", TCLinSolverBase::Tolerance, norm_r);
    TCLinSolverBase::WriteHistory(i, norm_r, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)TCLinSolverBase::SubSpaceSize; i++) {
    
    /*--- Compute rho_prime ---*/
    
    rho_prime = rho;
    
    /*--- Compute rho_i ---*/
    
    rho = dotProd(r, r_0);
    
    /*--- Compute beta ---*/
    
    beta = (rho / rho_prime) * (alpha /omega);
    
    /*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/
    
    beta_omega = -beta*omega;
    p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
    p.Plus_AX(1.0, r);
    
    /*--- Preconditioning step ---*/
    
    precond(p, phat);
    mat_vec(phat, v);
    
    /*--- Calculate step-length alpha ---*/
    
    r_0_v = dotProd(r_0, v);
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
    
    /*--- Check if solution has converged, else TCLinSolverBase::Output the relative residual if necessary ---*/
    
    norm_r = r.norm();
    if (norm_r < TCLinSolverBase::Tolerance*norm0) break;
    if (((TCLinSolverBase::Output) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0) && (rank == MASTER_NODE)) TCLinSolverBase::WriteHistory(i+1, norm_r, norm0);
    
  }
  
  if ((TCLinSolverBase::Output) && (rank == MASTER_NODE)) {
    cout << "# BCGSTAB final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
  }
  
    /*--- Recalculate final residual (this should be optional) ---*/
  if (TCLinSolverBase::Output) {
    mat_vec(x, A_x);
    r = b; r -= A_x;
    CalcType true_res = r.norm();
    
    if ((fabs(true_res - norm_r) > TCLinSolverBase::Tolerance*10.0) && (rank == MASTER_NODE)) {
      cout << "# WARNING in TCSysSolve<CalcType>::BCGSTAB_LinSolver(): " << endl;
      cout << "# true residual norm and calculated residual norm do not agree." << endl;
      cout << "# true_res = " << true_res <<", calc_res = " << norm_r <<", tol = " << TCLinSolverBase::Tolerance*10 <<"."<< endl;
      cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << endl;
    }
  }
  
  TCLinSolverBase::Residual = norm_r;
  return (unsigned long) i;
}
template<class CalcType, class BaseType>
TCSysSolve<CalcType, BaseType>::TCSysSolve(){
  
}

template<class CalcType, class BaseType>
TCSysSolve<CalcType, BaseType>::~TCSysSolve(){
  
//  delete LinSolver;
//  delete MatVec;
//  delete MatVec_b;
//  delete Precond;
//  delete Precond_b;
  
}
template<>
void TCSysSolve<su2double, su2double>::Initialize_System_Adjoint(unsigned short BlockSize, CGeometry *geometry, CConfig *config){
  
}
template<class CalcType, class BaseType>
void TCSysSolve<CalcType, BaseType>::Initialize_System(unsigned short blocksize, bool edgeconnect, CGeometry *geometry, CConfig *config){
  
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  BlockSize    = blocksize;
  
  /*--- Create Matrix structure ---*/
  
  Matrix.Initialize(nPoint, nPointDomain, BlockSize, BlockSize, edgeconnect, geometry, config);
  
  convert.Initialize(BlockSize);
  
  Initialize_System_Adjoint(blocksize, geometry, config);
  
}


template<>
void TCSysSolve<su2double, su2double>::Initialize_Linear_Solver_Adjoint(unsigned short blocksize,
                                                                        unsigned short kind_solver, 
                                                                        unsigned short kind_preconditioner, 
                                                                        unsigned long max_iter,
                                                                        su2double solver_error, 
                                                                        CGeometry *geometry, 
                                                                        CConfig *config){

}

template<class CalcType, class BaseType>
void TCSysSolve<CalcType, BaseType>::Initialize_Linear_Solver(unsigned short blocksize,
                                                                unsigned short kind_solver, 
                                                                unsigned short kind_preconditioner, 
                                                                unsigned long max_iter,
                                                                BaseType solver_error, 
                                                                CGeometry *geometry, 
                                                                CConfig *config){
  
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  BlockSize    = blocksize;
  
  CalcType SolverTol    = convert.ToCalcType(solver_error);
  
  kind_prec = kind_preconditioner;
  
  /*--- Solve the linear system using a Krylov subspace method ---*/
  
  MatVec = new TCSysMatrixVectorProduct<CalcType>(Matrix, geometry, config);
  
  switch (kind_preconditioner) {
  case JACOBI:
    Precond = new TCJacobiPreconditioner<CalcType>(Matrix, geometry, config);
    break;
  case ILU:
    Precond = new TCILUPreconditioner<CalcType>(Matrix, geometry, config);
    break;
  case LU_SGS:
    Precond = new TCLU_SGSPreconditioner<CalcType>(Matrix, geometry, config);
    break;
  case LINELET:
    Precond = new TCLineletPreconditioner<CalcType>(Matrix, geometry, config);
    break;
  default:
    Precond = new TCJacobiPreconditioner<CalcType>(Matrix, geometry, config);
    break;
  }
  switch (kind_solver) {
  case BCGSTAB:
    LinSolver = new TCLinSolver_BCGSTAB<CalcType>(blocksize, nPoint, nPointDomain, SolverTol, max_iter, false);
    break;
  case FGMRES: case RESTARTED_FGMRES:
    LinSolver = new TCLinSolver_FGMRES<CalcType>(blocksize, nPoint, nPointDomain, SolverTol, max_iter, true);
    break;
  case CONJUGATE_GRADIENT:
    LinSolver = new TCLinSolver_CG<CalcType>(blocksize, nPoint, nPointDomain, SolverTol, max_iter, false);
    break;
  default:
    LinSolver = new TCLinSolver_FGMRES<CalcType>(blocksize, nPoint, nPointDomain, SolverTol, max_iter, false);
    break;
  }
  

  
  Initialize_Linear_Solver_Adjoint(blocksize, kind_solver, kind_preconditioner, max_iter, solver_error, geometry, config);
  
}


template<class CalcType, class BaseType>
void TCSysSolve<CalcType, BaseType>::Build_Preconditioner(unsigned short kind_prec, bool transpose){
  switch (kind_prec) {
    case JACOBI:
      Matrix.BuildJacobiPreconditioner(transpose);
      break;
    case ILU:
      Matrix.BuildILUPreconditioner(transpose);
      break;
    case LU_SGS:
      break;
    default:
      Matrix.BuildJacobiPreconditioner(transpose);
      break;
  }
}

template<>
void TCSysSolve<su2double, su2double>::Solve_System(TCSysVector<su2double>& Rhs, TCSysVector<su2double>& Sol){
  Build_Preconditioner(kind_prec, false);  
  LinSolver->Solve(Rhs, Sol, *MatVec, *Precond);
}
#ifdef CODI_REVERSE_TYPE
template<>
void TCSysSolve<passivedouble, su2double>::Initialize_System_Adjoint(unsigned short blocksize, CGeometry *geometry, CConfig *config){
  
  LinSysRes_calc.Initialize(nPoint, nPointDomain, BlockSize, 0.0);
  LinSysSol_calc.Initialize(nPoint, nPointDomain, BlockSize, 0.0);
  
}
template<>
void TCSysSolve<passivedouble, su2double>::Initialize_Linear_Solver_Adjoint(unsigned short blocksize,
                                                                        unsigned short kind_solver, 
                                                                        unsigned short kind_preconditioner, 
                                                                        unsigned long max_iter,
                                                                        su2double solver_error, 
                                                                        CGeometry *geometry, 
                                                                        CConfig *config){

  kind_prec_b = ILU;
  
  MatVec_b = new TCSysMatrixVectorProductTransposed<passivedouble>(Matrix, geometry, config);
  Precond_b = new TCILUPreconditioner<passivedouble>(Matrix, geometry, config);
  
}

template<>
void TCSysSolve<passivedouble, su2double>::Solve_System(TCSysVector<su2double> &Rhs, TCSysVector<su2double> &Sol){
  
  bool TapeActive = NO;
  unsigned long iPoint;
  
  TapeActive = AD::globalTape.isActive();
  
  if (TapeActive){
    /*--- Initialize the external function helper with storage of primal input and output disabled  ---*/
    
    AD::InitExtFunc(false, false);
    
    /*--- Set the right-hand side of the linear system as input of the external function ---*/
    
    AD::SetExtFuncIn(Rhs.vec_val, Rhs.GetLocSize());
    
    /*--- Stop the recording ---*/
    
    AD::StopRecording();
    
  }
  
  /*--- Convert data from the basetype to the calc type---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    LinSysSol_calc.SetBlock(iPoint, convert.ToCalcType(Rhs.GetBlock(iPoint), BlockSize));
    LinSysRes_calc.SetBlock(iPoint, convert.ToCalcType(Sol.GetBlock(iPoint), BlockSize));
  }
  
  Build_Preconditioner(kind_prec, false);
  LinSolver->Solve(LinSysRes_calc, LinSysSol_calc, *MatVec, *Precond);
  
  /*--- Convert data back from the calc type to the base type---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    Rhs.SetBlock(iPoint, convert.ToBaseType(LinSysSol_calc.GetBlock(iPoint), BlockSize));
    Sol.SetBlock(iPoint, convert.ToBaseType(LinSysRes_calc.GetBlock(iPoint), BlockSize));
  }
  
  if(TapeActive) {
    
    AD::StartRecording();
    
    /*--- Set the solution of the linear system as input of the external function ---*/
    
    AD::SetExtFuncOut(Sol.vec_val, Sol.GetLocSize());

    /*--- Add some pointers to additional data that is required in the reverse sweep ---*/
    
    AD::FuncHelper->addUserData(&LinSysRes_calc);
    AD::FuncHelper->addUserData(&LinSysSol_calc);
    AD::FuncHelper->addUserData(MatVec_b);
    AD::FuncHelper->addUserData(Precond_b);
    AD::FuncHelper->addUserData(LinSolver);
    AD::FuncHelper->addToTape(CSysSolve_b::Solve_b<passivedouble>);    
    
    /*--- Build the preconditioner for the transposed system ---*/
    
    Build_Preconditioner(kind_prec_b, true);
    
    delete AD::FuncHelper;
  }
  
}

#endif
