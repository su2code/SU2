/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

template<class ScalarType>
CSysSolve<ScalarType>::CSysSolve(const bool mesh_deform_mode) : cg_ready(false), bcg_ready(false), gmres_ready(false) {

  mesh_deform = mesh_deform_mode;
  LinSysRes_ptr = NULL;
  LinSysSol_ptr = NULL;
}

template<class ScalarType>
void CSysSolve<ScalarType>::ApplyGivens(const ScalarType & s, const ScalarType & c, ScalarType & h1, ScalarType & h2) {
  
  ScalarType temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

template<class ScalarType>
void CSysSolve<ScalarType>::GenerateGivens(ScalarType & dx, ScalarType & dy, ScalarType & s, ScalarType & c) {
  
  if ( (dx == 0.0) && (dy == 0.0) ) {
    c = 1.0;
    s = 0.0;
  }
  else if ( fabs(dy) > fabs(dx) ) {
    ScalarType tmp = dx/dy;
    dx = sqrt(1.0 + tmp*tmp);
    s = Sign(1.0/dx, dy);
    c = tmp*s;
  }
  else if ( fabs(dy) <= fabs(dx) ) {
    ScalarType tmp = dy/dx;
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

template<class ScalarType>
void CSysSolve<ScalarType>::SolveReduced(const int & n, const vector<vector<ScalarType> > & Hsbg,
                             const vector<ScalarType> & rhs, vector<ScalarType> & x) {
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

template<class ScalarType>
void CSysSolve<ScalarType>::ModGramSchmidt(int i, vector<vector<ScalarType> > & Hsbg, vector<CSysVector<ScalarType> > & w) {

  bool Convergence = true;

  /*--- Parameter for reorthonormalization ---*/
  
  static const ScalarType reorth = 0.98;
  
  /*--- Get the norm of the vector being orthogonalized, and find the
  threshold for re-orthogonalization ---*/
  
  ScalarType nrm = dotProd(w[i+1], w[i+1]);
  ScalarType thr = nrm*reorth;
  
  /*--- The norm of w[i+1] < 0.0 or w[i+1] = NaN ---*/

  if ((nrm <= 0.0) || (nrm != nrm)) Convergence = false;
  
  /*--- Synchronization point to check the convergence of the solver ---*/

#ifdef HAVE_MPI

  int rank = SU2_MPI::GetRank();
  int size = SU2_MPI::GetSize();
  
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
    ScalarType prod = dotProd(w[i+1], w[k]);
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

template<class ScalarType>
void CSysSolve<ScalarType>::WriteHeader(const string & solver, const ScalarType & restol, const ScalarType & resinit) {
  
  cout << "\n# " << solver << " residual history" << endl;
  cout << "# Residual tolerance target = " << restol << endl;
  cout << "# Initial residual norm     = " << resinit << endl;
  
}

template<class ScalarType>
void CSysSolve<ScalarType>::WriteHistory(const int & iter, const ScalarType & res, const ScalarType & resinit) {
  
  cout << "     " << iter << "     " << res/resinit << endl;
  
}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::CG_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                  CMatrixVectorProduct<ScalarType> & mat_vec, CPreconditioner<ScalarType> & precond,
                                                  ScalarType tol, unsigned long m, ScalarType *residual, bool monitoring, CConfig *config) {

  int rank = SU2_MPI::GetRank();
  ScalarType norm_r = 0.0, norm0 = 0.0;
  int i = 0;
  
  /*--- Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet ---*/

  if (!cg_ready) {
    A_x = b;
    z = b;
    cg_ready = true;
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  
  mat_vec(x, A_x);
  r = b; r -= A_x;

  /*--- Only compute the residuals in full communication mode. ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
    
    norm_r = r.norm();
    norm0  = b.norm();
    if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
      if (rank == MASTER_NODE) cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
      return 0;
    }
    
    /*--- Set the norm to the initial initial residual value ---*/
    
    norm0 = norm_r;
    
    /*--- Output header information including initial residual ---*/
    
    if ((monitoring) && (rank == MASTER_NODE)) {
      WriteHeader("CG", tol, norm_r);
      WriteHistory(i, norm_r, norm0);
    }
    
  }
  
  ScalarType alpha, beta, r_dot_z;
  precond(r, z);
  p = z;

  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)m; i++) {
    
    /*--- Apply matrix to p to build Krylov subspace ---*/
    
    mat_vec(p, A_x);
    
    /*--- Calculate step-length alpha ---*/
    
    r_dot_z = dotProd(r, z);
    alpha = dotProd(A_x, p);
    alpha = r_dot_z / alpha;
    
    /*--- Update solution and residual: ---*/
    
    x.Plus_AX(alpha, p);
    r.Plus_AX(-alpha, A_x);
    
    /*--- Only compute the residuals in full communication mode. ---*/
    
    if (config->GetComm_Level() == COMM_FULL) {
      
      /*--- Check if solution has converged, else output the relative residual if necessary ---*/
      
      norm_r = r.norm();
      if (norm_r < tol*norm0) break;
      if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0)) WriteHistory(i+1, norm_r, norm0);
      
    }
    
    precond(r, z);
    
    /*--- Calculate Gram-Schmidt coefficient beta,
		 beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/
    
    beta = 1.0 / r_dot_z;
    r_dot_z = dotProd(r, z);
    beta *= r_dot_z;
    
    /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/
    
    p.Equals_AX_Plus_BY(beta, p, 1.0, z);
    
  }
  
  /*--- Recalculate final residual (this should be optional) ---*/
  
  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {
    
    if (rank == MASTER_NODE) {
      cout << "# Conjugate Gradient final (true) residual:" << endl;
      cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
    }
    
    mat_vec(x, A_x);
    r = b; r -= A_x;
    ScalarType true_res = r.norm();
    
    if (fabs(true_res - norm_r) > tol*10.0) {
      if (rank == MASTER_NODE) {
        cout << "# WARNING in CSysSolve::CG_LinSolver(): " << endl;
        cout << "# true residual norm and calculated residual norm do not agree." << endl;
        cout << "# true_res = " << true_res <<", calc_res = " << norm_r <<", tol = " << tol*10 <<"."<< endl;
        cout << "# true_res - calc_res = " << true_res - norm_r << endl;
      }
    }
    
  }
  
  (*residual) = norm_r;
	return (unsigned long) i;
  
}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                      CMatrixVectorProduct<ScalarType> & mat_vec, CPreconditioner<ScalarType> & precond,
                                                      ScalarType tol, unsigned long m, ScalarType *residual, bool monitoring, CConfig *config) {
	
  int rank = SU2_MPI::GetRank();
  
  /*---  Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  /*---  Check the subspace size ---*/
  
  if (m > 5000) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size (too high), m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet
   Note: elements in w and z are initialized to x to avoid creating
	 a temporary CSysVector object for the copy constructor ---*/

  if (!gmres_ready) {
    W.resize(m+1, x);
    Z.resize(m+1, x);
    gmres_ready = true;
  }

  /*---  Define various arrays ---*/

  vector<ScalarType> g(m+1, 0.0);
  vector<ScalarType> sn(m+1, 0.0);
  vector<ScalarType> cs(m+1, 0.0);
  vector<ScalarType> y(m, 0.0);
  vector<vector<ScalarType> > H(m+1, vector<ScalarType>(m, 0.0));
  
  /*---  Calculate the norm of the rhs vector ---*/
  
  ScalarType norm0 = b.norm();
  
  /*---  Calculate the initial residual (actually the negative residual)
	 and compute its norm ---*/
  
  mat_vec(x, W[0]);
  W[0] -= b;
  
  ScalarType beta = W[0].norm();
  
  if ( (beta < tol*norm0) || (beta < eps) ) {
    
    /*---  System is already solved ---*/
    
    if (rank == MASTER_NODE) cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*---  Normalize residual to get w_{0} (the negative sign is because w[0]
	 holds the negative residual, as mentioned above) ---*/
  
  W[0] /= -beta;
  
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
    
    /*---  Precondition the CSysVector w[i] and store result in z[i] ---*/
    
    precond(W[i], Z[i]);
    
    /*---  Add to Krylov subspace ---*/
    
    mat_vec(Z[i], W[i+1]);
    
    /*---  Modified Gram-Schmidt orthogonalization ---*/
    
    ModGramSchmidt(i, H, W);
    
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
    x.Plus_AX(y[k], Z[k]);
  }
  
  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
  
  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {
    
    if (rank == MASTER_NODE) {
      cout << "# FGMRES final (true) residual:" << endl;
      cout << "# Iteration = " << i << ": |res|/|res0| = " << beta/norm0 << ".\n" << endl;
    }
    
    mat_vec(x, W[0]);
    W[0] -= b;
    ScalarType res = W[0].norm();
    
    if (fabs(res - beta) > tol*10) {
      if (rank == MASTER_NODE) {
        cout << "# WARNING in CSysSolve::FGMRES_LinSolver(): " << endl;
        cout << "# true residual norm and calculated residual norm do not agree." << endl;
        cout << "# res = " << res <<", beta = " << beta <<", tol = " << tol*10 <<"."<< endl;
        cout << "# res - beta = " << res - beta << endl << endl;
      }
    }
    
  }
  
  (*residual) = beta;
  return (unsigned long) i;
  
}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::BCGSTAB_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                       CMatrixVectorProduct<ScalarType> & mat_vec, CPreconditioner<ScalarType> & precond,
                                                       ScalarType tol, unsigned long m, ScalarType *residual, bool monitoring, CConfig *config) {
  
  int rank = SU2_MPI::GetRank();
  ScalarType norm_r = 0.0, norm0 = 0.0;
  int i = 0;
  
  /*--- Check the subspace size ---*/
  
  if (m < 1) {
    char buf[100];
    SPRINTF(buf, "Illegal value for subspace size, m = %lu", m );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet ---*/

  if (!bcg_ready) {
    A_x = b;
    p = b;
    z = b;
    v = b;
    bcg_ready = true;
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

  mat_vec(x, A_x);
  r = b; r -= A_x;

  /*--- Only compute the residuals in full communication mode. ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
    
    norm_r = r.norm();
    norm0  = b.norm();
    if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
      if (rank == MASTER_NODE) cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
      return 0;
    }
    
    /*--- Set the norm to the initial initial residual value ---*/
    
    norm0 = norm_r;
    
    /*--- Output header information including initial residual ---*/
    
    if ((monitoring) && (rank == MASTER_NODE)) {
      WriteHeader("BCGSTAB", tol, norm_r);
      WriteHistory(i, norm_r, norm0);
    }
    
  }
  
  /*--- Initialization ---*/
  
  ScalarType alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
  p = ScalarType(0.0); v = ScalarType(0.0); r_0 = r;
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < (int)m; i++) {
    
    /*--- Compute rho_prime ---*/
    
    rho_prime = rho;
    
    /*--- Compute rho_i ---*/
    
    rho = dotProd(r, r_0);
    
    /*--- Compute beta ---*/
    
    beta = (rho / rho_prime) * (alpha /omega);
    
    /*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/
    
    ScalarType beta_omega = -beta*omega;
    p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
    p.Plus_AX(1.0, r);
    
    /*--- Preconditioning step ---*/
    
    precond(p, z);
    mat_vec(z, v);
    
    /*--- Calculate step-length alpha ---*/
    
    ScalarType r_0_v = dotProd(r_0, v);
    alpha = rho / r_0_v;

    /*--- Update solution and residual: ---*/

    /*--- x_{i-1/2} = x_{i-1} + alpha * z ---*/
    x.Plus_AX(alpha, z);
    /*--- r_{i-1/2} = r_{i-1} - alpha * v_{i} ---*/
    r.Plus_AX(-alpha, v);

    /*--- Preconditioning step ---*/
    
    precond(r, z);
    mat_vec(z, A_x);
    
    /*--- Calculate step-length omega ---*/
    
    omega = dotProd(A_x, r) / dotProd(A_x, A_x);
    
    /*--- Update solution and residual: ---*/
    
    /*--- x_{i} = x_{i-1/2} + omega * z ---*/
    x.Plus_AX(omega, z);
    /*--- r_{i} = r_{i-1/2} - omega * A * z ---*/
    r.Plus_AX(-omega, A_x);
    
    /*--- Only compute the residuals in full communication mode. ---*/
    
    if (config->GetComm_Level() == COMM_FULL) {
      
      /*--- Check if solution has converged, else output the relative residual if necessary ---*/
      
      norm_r = r.norm();
      if (norm_r < tol*norm0) break;
      if (((monitoring) && (rank == MASTER_NODE)) && ((i+1) % 10 == 0) && (rank == MASTER_NODE)) WriteHistory(i+1, norm_r, norm0);
      
    }
    
  }
  
  /*--- Recalculate final residual (this should be optional) ---*/
  
  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {
    
    if (rank == MASTER_NODE) {
      cout << "# BCGSTAB final (true) residual:" << endl;
      cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << ".\n" << endl;
    }
    
    mat_vec(x, A_x);
    r = b; r -= A_x;
    ScalarType true_res = r.norm();
    
    if ((fabs(true_res - norm_r) > tol*10.0) && (rank == MASTER_NODE)) {
      cout << "# WARNING in CSysSolve::BCGSTAB_LinSolver(): " << endl;
      cout << "# true residual norm and calculated residual norm do not agree." << endl;
      cout << "# true_res = " << true_res <<", calc_res = " << norm_r <<", tol = " << tol*10 <<"."<< endl;
      cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << endl;
    }
    
  }
  
  (*residual) = norm_r;
  return (unsigned long) i;
}

template<>
void CSysSolve<su2double>::HandleTemporariesIn(CSysVector<su2double> & LinSysRes, CSysVector<su2double> & LinSysSol) {

  /*--- When the type is the same the temporaties are not required ---*/
  /*--- Set the pointers ---*/
  LinSysRes_ptr = &LinSysRes;
  LinSysSol_ptr = &LinSysSol;
}

template<>
void CSysSolve<su2double>::HandleTemporariesOut(CSysVector<su2double> & LinSysSol) {

  /*--- When the type is the same the temporaties are not required ---*/
  /*--- Reset the pointers ---*/
  LinSysRes_ptr = NULL;
  LinSysSol_ptr = NULL;
}

#ifdef CODI_REVERSE_TYPE
template<>
void CSysSolve<passivedouble>::HandleTemporariesIn(CSysVector<su2double> & LinSysRes, CSysVector<su2double> & LinSysSol) {

  /*--- When the type is different we need to copy data to the temporaries ---*/
  /*--- Copy data, the solution is also copied because it serves as initial conditions ---*/
  LinSysRes_tmp.PassiveCopy(LinSysRes);
  LinSysSol_tmp.PassiveCopy(LinSysSol);

  /*--- Set the pointers ---*/
  LinSysRes_ptr = &LinSysRes_tmp;
  LinSysSol_ptr = &LinSysSol_tmp;
}

template<>
void CSysSolve<passivedouble>::HandleTemporariesOut(CSysVector<su2double> & LinSysSol) {

  /*--- When the type is different we need to copy data from the temporaries ---*/
  /*--- Copy data, only the solution needs to be copied ---*/
  LinSysSol.PassiveCopy(LinSysSol_tmp);

  /*--- Reset the pointers ---*/
  LinSysRes_ptr = NULL;
  LinSysSol_ptr = NULL;
}
#endif

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::Solve(CSysMatrix<ScalarType> & Jacobian, CSysVector<su2double> & LinSysRes,
                                           CSysVector<su2double> & LinSysSol, CGeometry *geometry, CConfig *config) {
  /*---
   A word about the templated types. It is assumed that the residual and solution vectors are always of su2doubles,
   meaning that they are active in the discrete adjoint. The same assumption is made in SetExternalSolve.
   When the Jacobian is passive (and therefore not compatible with the vectors) we go through the "HandleTemporaries"
   mechanisms. Note that CG, BCGSTAB, and FGMRES, all expect the vector to be compatible with the Product and
   Preconditioner (and therefore with the Matrix). Likewise for Solve_b (which is used by CSysSolve_b).
   There are no provisions here for active Matrix and passive Vectors as that makes no sense since we only handle the
   derivatives of the residual in CSysSolve_b.
  ---*/

  unsigned short KindSolver, KindPrecond;
  unsigned long MaxIter, RestartIter, IterLinSol = 0;
  ScalarType SolverTol, Norm0 = 0.0;
  bool ScreenOutput;

  /*--- Normal mode ---*/

  if(!mesh_deform) {

    KindSolver   = config->GetKind_Linear_Solver();
    KindPrecond  = config->GetKind_Linear_Solver_Prec();
    MaxIter      = config->GetLinear_Solver_Iter();
    RestartIter  = config->GetLinear_Solver_Restart_Frequency();
    SolverTol    = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
    ScreenOutput = false;
  }

  /*--- Mesh Deformation mode ---*/

  else {

    KindSolver   = config->GetKind_Deform_Linear_Solver();
    KindPrecond  = config->GetKind_Deform_Linear_Solver_Prec();
    MaxIter      = config->GetDeform_Linear_Solver_Iter();
    RestartIter  = config->GetLinear_Solver_Restart_Frequency();
    SolverTol    = SU2_TYPE::GetValue(config->GetDeform_Linear_Solver_Error());
    ScreenOutput = config->GetDeform_Output();
  }

  CMatrixVectorProduct<ScalarType> *mat_vec = NULL;

  bool TapeActive = NO;

  if (config->GetDiscrete_Adjoint()) {
#ifdef CODI_REVERSE_TYPE
    
    TapeActive = AD::globalTape.isActive();

    AD::StartExtFunc(false, false);

    AD::SetExtFuncIn(&LinSysRes[0], LinSysRes.GetLocSize());

    /*--- Stop the recording for the linear solver ---*/

    AD::StopRecording();
#endif
  }

  /*--- Solve the linear system using a Krylov subspace method ---*/

  HandleTemporariesIn(LinSysRes, LinSysSol);

  if (KindSolver == BCGSTAB || KindSolver == CONJUGATE_GRADIENT ||
      KindSolver == FGMRES  || KindSolver == RESTARTED_FGMRES ) {
    
    mat_vec = new CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);
    CPreconditioner<ScalarType>* precond = NULL;
    
    switch (KindPrecond) {
      case JACOBI:
        Jacobian.BuildJacobiPreconditioner();
        precond = new CJacobiPreconditioner<ScalarType>(Jacobian, geometry, config);
        break;
      case ILU:
        Jacobian.BuildILUPreconditioner();
        precond = new CILUPreconditioner<ScalarType>(Jacobian, geometry, config);
        break;
      case LU_SGS:
        precond = new CLU_SGSPreconditioner<ScalarType>(Jacobian, geometry, config);
        break;
      case LINELET:
        Jacobian.BuildJacobiPreconditioner();
        precond = new CLineletPreconditioner<ScalarType>(Jacobian, geometry, config);
        break;
      default:
        Jacobian.BuildJacobiPreconditioner();
        precond = new CJacobiPreconditioner<ScalarType>(Jacobian, geometry, config);
        break;
    }
    
    switch (KindSolver) {
      case BCGSTAB:
        IterLinSol = BCGSTAB_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol, MaxIter, &Residual, ScreenOutput, config);
        break;
      case FGMRES:
        IterLinSol = FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol, MaxIter, &Residual, ScreenOutput, config);
        break;
      case CONJUGATE_GRADIENT:
        IterLinSol = CG_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol, MaxIter, &Residual, ScreenOutput, config);
        break;
      case RESTARTED_FGMRES:
        IterLinSol = 0;
        Norm0 = LinSysRes_ptr->norm();
        while (IterLinSol < MaxIter) {
          /*--- Enforce a hard limit on total number of iterations ---*/
          unsigned long IterLimit = min(RestartIter, MaxIter-IterLinSol);
          IterLinSol += FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol, IterLimit, &Residual, ScreenOutput, config);
          if ( Residual < SolverTol*Norm0 ) break;
        }
        break;
    }
    
    /*--- Dealocate memory of the Krylov subspace method ---*/
    
    delete mat_vec;
    delete precond;
    
  }
  
  /*--- Smooth the linear system. ---*/
  
  else {
    switch (KindSolver) {
      case SMOOTHER_LUSGS:
        mat_vec = new CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);
        IterLinSol = Jacobian.LU_SGS_Smoother(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, SolverTol, MaxIter, &Residual, ScreenOutput, geometry, config);
        delete mat_vec;
        break;
      case SMOOTHER_JACOBI:
        mat_vec = new CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);
        Jacobian.BuildJacobiPreconditioner();
        IterLinSol = Jacobian.Jacobi_Smoother(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, SolverTol, MaxIter, &Residual, ScreenOutput, geometry, config);
        delete mat_vec;
        break;
      case SMOOTHER_ILU:
        mat_vec = new CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);
        Jacobian.BuildILUPreconditioner();
        IterLinSol = Jacobian.ILU_Smoother(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, SolverTol, MaxIter, &Residual, ScreenOutput, geometry, config);
        delete mat_vec;
        break;
      case SMOOTHER_LINELET:
        Jacobian.BuildJacobiPreconditioner();
        Jacobian.ComputeLineletPreconditioner(*LinSysRes_ptr, *LinSysSol_ptr, geometry, config);
        IterLinSol = 1;
        break;
    }
  }

  HandleTemporariesOut(LinSysSol);

  if(TapeActive) {

    bool RequiresTranspose = !mesh_deform; // jacobian is symmetric
    if (!mesh_deform) KindPrecond = config->GetKind_DiscAdj_Linear_Prec();
    else              KindPrecond = config->GetKind_Deform_Linear_Solver_Prec();
	
    /*--- Start recording if it was stopped for the linear solver ---*/

    AD::StartRecording();
    
    AD::SetExtFuncOut(&LinSysSol[0], (int)LinSysSol.GetLocSize());

#ifdef CODI_REVERSE_TYPE
    AD::FuncHelper->addUserData(&LinSysRes);
    AD::FuncHelper->addUserData(&LinSysSol);
    AD::FuncHelper->addUserData(&Jacobian);
    AD::FuncHelper->addUserData(geometry);
    AD::FuncHelper->addUserData(config);
    AD::FuncHelper->addUserData(this);
    AD::FuncHelper->addToTape(CSysSolve_b<ScalarType>::Solve_b);
#endif   

    /*--- Build preconditioner for the transposed Jacobian ---*/

    switch(KindPrecond) {
      case ILU:
        Jacobian.BuildILUPreconditioner(RequiresTranspose);
        break;
      case JACOBI:
        Jacobian.BuildJacobiPreconditioner(RequiresTranspose);
        break;
      default:
        SU2_MPI::Error("The specified preconditioner is not yet implemented for the discrete adjoint method.", CURRENT_FUNCTION);
        break;
    }

    AD::EndExtFunc();
  }

  return IterLinSol;
}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::Solve_b(CSysMatrix<ScalarType> & Jacobian, CSysVector<su2double> & LinSysRes,
                                             CSysVector<su2double> & LinSysSol, CGeometry *geometry, CConfig *config) {
#ifdef CODI_REVERSE_TYPE

  unsigned short KindSolver, KindPrecond;
  unsigned long MaxIter, RestartIter, IterLinSol = 0;
  ScalarType SolverTol, Norm0 = 0.0;
  bool ScreenOutput;

  /*--- Normal mode ---*/

  if(!mesh_deform) {

    KindSolver   = config->GetKind_DiscAdj_Linear_Solver();
    KindPrecond  = config->GetKind_DiscAdj_Linear_Prec();
    MaxIter      = config->GetLinear_Solver_Iter();
    RestartIter  = config->GetLinear_Solver_Restart_Frequency();
    SolverTol    = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
    ScreenOutput = false;
  }

  /*--- Mesh Deformation mode ---*/

  else {

    KindSolver   = config->GetKind_Deform_Linear_Solver();
    KindPrecond  = config->GetKind_Deform_Linear_Solver_Prec();
    MaxIter      = config->GetDeform_Linear_Solver_Iter();
    RestartIter  = config->GetLinear_Solver_Restart_Frequency();
    SolverTol    = SU2_TYPE::GetValue(config->GetDeform_Linear_Solver_Error());
    ScreenOutput = config->GetDeform_Output();
  }

  /*--- Set up preconditioner and matrix-vector product ---*/

  CPreconditioner<ScalarType>* precond  = NULL;

  switch(KindPrecond) {
    case ILU:
      precond = new CILUPreconditioner<ScalarType>(Jacobian, geometry, config);
      break;
    case JACOBI:
      precond = new CJacobiPreconditioner<ScalarType>(Jacobian, geometry, config);
      break;
  }

  CMatrixVectorProduct<ScalarType>* mat_vec = new CSysMatrixVectorProductTransposed<ScalarType>(Jacobian, geometry, config);

  /*--- Solve the system ---*/

  HandleTemporariesIn(LinSysRes, LinSysSol);

  switch(KindSolver) {
    case FGMRES:
      IterLinSol = FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol , MaxIter, &Residual, ScreenOutput, config);
      break;
    case BCGSTAB:
      IterLinSol = BCGSTAB_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol , MaxIter, &Residual, ScreenOutput, config);
      break;
    case CONJUGATE_GRADIENT:
      IterLinSol = CG_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol, MaxIter, &Residual, ScreenOutput, config);
      break;
    case RESTARTED_FGMRES:
      IterLinSol = 0;
      Norm0 = LinSysRes_ptr->norm();
      while (IterLinSol < MaxIter) {
        /*--- Enforce a hard limit on total number of iterations ---*/
        unsigned long IterLimit = min(RestartIter, MaxIter-IterLinSol);
        IterLinSol += FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, *mat_vec, *precond, SolverTol , IterLimit, &Residual, ScreenOutput, config);
        if ( Residual < SolverTol*Norm0 ) break;
      }
      break;
    default:
      SU2_MPI::Error("The specified linear solver is not yet implemented for the discrete adjoint method.", CURRENT_FUNCTION);
      break;
  }

  HandleTemporariesOut(LinSysSol);

  delete mat_vec;
  delete precond;

  return IterLinSol;
#else
  return 0;
#endif
}

/*--- Explicit instantiations ---*/
template class CSysSolve<su2double>;

#ifdef CODI_REVERSE_TYPE
template class CSysSolve<passivedouble>;
#endif
