/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/linear_algebra/CSysSolve.hpp"
#include "../../include/linear_algebra/CSysSolve_b.hpp"
#include "../../include/omp_structure.hpp"
#include "../../include/option_structure.hpp"
#include "../../include/config_structure.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../include/linear_algebra/CPreconditioner.hpp"

#include <limits>

const su2double eps = numeric_limits<passivedouble>::epsilon(); /*!< \brief machine epsilon */

template<class ScalarType>
CSysSolve<ScalarType>::CSysSolve(const bool mesh_deform_mode) : cg_ready(false), bcg_ready(false),
                                                                gmres_ready(false), smooth_ready(false) {
  mesh_deform = mesh_deform_mode;
  LinSysRes_ptr = nullptr;
  LinSysSol_ptr = nullptr;
  Residual = 0.0;
}

template<class ScalarType>
void CSysSolve<ScalarType>::ApplyGivens(ScalarType s, ScalarType c, ScalarType & h1, ScalarType & h2) const {

  ScalarType temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

template<class ScalarType>
void CSysSolve<ScalarType>::GenerateGivens(ScalarType & dx, ScalarType & dy, ScalarType & s, ScalarType & c) const {

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
void CSysSolve<ScalarType>::SolveReduced(int n, const vector<vector<ScalarType> > & Hsbg,
                                         const vector<ScalarType> & rhs, vector<ScalarType> & x) const {
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
void CSysSolve<ScalarType>::ModGramSchmidt(int i, vector<vector<ScalarType> > & Hsbg,
                                           vector<CSysVector<ScalarType> > & w) const {

  /*--- Parameter for reorthonormalization ---*/

  const ScalarType reorth = 0.98;

  /*--- Get the norm of the vector being orthogonalized, and find the
  threshold for re-orthogonalization ---*/

  ScalarType nrm = w[i+1].squaredNorm();
  ScalarType thr = nrm*reorth;

  /*--- The norm of w[i+1] < 0.0 or w[i+1] = NaN ---*/

  if ((nrm <= 0.0) || (nrm != nrm)) {
    /*--- nrm is the result of a dot product, communications are implicitly handled. ---*/
    SU2_OMP_MASTER
    SU2_MPI::Error("FGMRES orthogonalization failed, linear solver diverged.", CURRENT_FUNCTION);
  }

  /*--- Begin main Gram-Schmidt loop ---*/

  for (int k = 0; k < i+1; k++) {
    ScalarType prod = w[i+1].dot(w[k]);
    Hsbg[k][i] = prod;
    w[i+1].Plus_AX(-prod, w[k]);

    /*--- Check if reorthogonalization is necessary ---*/

    if (prod*prod > thr) {
      prod = w[i+1].dot(w[k]);
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
void CSysSolve<ScalarType>::WriteHeader(string solver, ScalarType restol, ScalarType resinit) const {

  cout << "\n# " << solver << " residual history\n";
  cout << "# Residual tolerance target = " << restol << "\n";
  cout << "# Initial residual norm     = " << resinit << endl;
}

template<class ScalarType>
void CSysSolve<ScalarType>::WriteHistory(unsigned long iter, ScalarType res) const {

  cout << "     " << iter << "     " << res << endl;
}

template<class ScalarType>
void CSysSolve<ScalarType>::WriteFinalResidual(string solver, unsigned long iter, ScalarType res) const {

  cout << "# " << solver << " final (true) residual:\n";
  cout << "# Iteration = " << iter << ": |res|/|res0| = " << res << ".\n" << endl;
}

template<class ScalarType>
void CSysSolve<ScalarType>::WriteWarning(ScalarType res_calc, ScalarType res_true, ScalarType tol) const {

  cout << "# WARNING:\n";
  cout << "# true residual norm and calculated residual norm do not agree.\n";
  cout << "# true_res = " << res_true << ", calc_res = " << res_calc << ", tol = " << tol*10 << ".\n";
  cout << "# true_res - calc_res = " << res_true - res_calc << endl;
}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::CG_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                  const CMatrixVectorProduct<ScalarType> & mat_vec, const CPreconditioner<ScalarType> & precond,
                                                  ScalarType tol, unsigned long m, ScalarType & residual, bool monitoring, CConfig *config) const {

  const bool master = (SU2_MPI::GetRank() == MASTER_NODE) && (omp_get_thread_num() == 0);
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Check the subspace size ---*/

  if (m < 1) {
    SU2_OMP_MASTER
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet, only one thread can
   *    do this since the working vectors are shared. ---*/

  if (!cg_ready) {
    SU2_OMP_MASTER
    {
      auto nVar = b.GetNVar();
      auto nBlk = b.GetNBlk();
      auto nBlkDomain = b.GetNBlkDomain();

      A_x.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      z.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      p.Initialize(nBlk, nBlkDomain, nVar, nullptr);

      cg_ready = true;
    }
    SU2_OMP_BARRIER
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

  mat_vec(x, A_x);
  r = b; r -= A_x;

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    norm_r = r.norm();
    norm0  = b.norm();
    if ((norm_r < tol*norm0) || (norm_r < eps)) {
      if (master) cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
      return 0;
    }

    /*--- Set the norm to the initial initial residual value ---*/

    norm0 = norm_r;

    /*--- Output header information including initial residual ---*/

    if ((monitoring) && (master)) {
      WriteHeader("CG", tol, norm_r);
      WriteHistory(i, norm_r/norm0);
    }

  }

  ScalarType alpha, beta, r_dot_z, r_dot_z_old;
  precond(r, z);
  p = z;
  r_dot_z = r.dot(z);

  /*---  Loop over all search directions ---*/

  for (i = 0; i < m; i++) {

    /*--- Apply matrix to p to build Krylov subspace ---*/

    mat_vec(p, A_x);

    /*--- Calculate step-length alpha ---*/

    alpha = r_dot_z / A_x.dot(p);

    /*--- Update solution and residual: ---*/

    x.Plus_AX(alpha, p);
    r.Plus_AX(-alpha, A_x);

    /*--- Only compute the residuals in full communication mode. ---*/

    if (config->GetComm_Level() == COMM_FULL) {

      /*--- Check if solution has converged, else output the relative residual if necessary ---*/

      norm_r = r.norm();
      if (norm_r < tol*norm0) break;
      if (((monitoring) && (master)) && ((i+1) % 10 == 0))
        WriteHistory(i+1, norm_r/norm0);

    }

    precond(r, z);

    /*--- Calculate Gram-Schmidt coefficient beta,
     beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/

    r_dot_z_old = r_dot_z;
    r_dot_z = r.dot(z);
    beta = r_dot_z / r_dot_z_old;

    /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/

    p.Equals_AX_Plus_BY(beta, p, 1.0, z);

  }

  /*--- Recalculate final residual (this should be optional) ---*/

  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {

    if (master) WriteFinalResidual("CG", i, norm_r/norm0);

    mat_vec(x, A_x);
    r = b; r -= A_x;
    ScalarType true_res = r.norm();

    if (fabs(true_res - norm_r) > tol*10.0) {
      if (master) {
        WriteWarning(norm_r, true_res, tol);
      }
    }

  }

  residual = norm_r/norm0;
  return i;

}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                      const CMatrixVectorProduct<ScalarType> & mat_vec, const CPreconditioner<ScalarType> & precond,
                                                      ScalarType tol, unsigned long m, ScalarType & residual, bool monitoring, CConfig *config) const {

  const bool master = (SU2_MPI::GetRank() == MASTER_NODE) && (omp_get_thread_num() == 0);

  /*---  Check the subspace size ---*/

  if (m < 1) {
    SU2_OMP_MASTER
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  if (m > 5000) {
    SU2_OMP_MASTER
    SU2_MPI::Error("FGMRES subspace is too large.", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet
   Note: elements in w and z are initialized to x to avoid creating
   a temporary CSysVector object for the copy constructor ---*/

  if (!gmres_ready) {
    SU2_OMP_MASTER
    {
      W.resize(m+1, x);
      Z.resize(m+1, x);
      gmres_ready = true;
    }
    SU2_OMP_BARRIER
  }

  /*--- Define various arrays. In parallel, each thread of each rank has and works
   on its own thread, since calculations on these arrays are based on dot products
   (reduced across all threads and ranks) all threads do the same computations. ---*/

  vector<ScalarType> g(m+1, 0.0);
  vector<ScalarType> sn(m+1, 0.0);
  vector<ScalarType> cs(m+1, 0.0);
  vector<ScalarType> y(m, 0.0);
  vector<vector<ScalarType> > H(m+1, vector<ScalarType>(m, 0.0));

  /*--- Calculate the norm of the rhs vector. ---*/

  ScalarType norm0 = b.norm();

  /*--- Calculate the initial residual (actually the negative residual) and compute its norm. ---*/

  mat_vec(x, W[0]);
  W[0] -= b;

  ScalarType beta = W[0].norm();

  if ((beta < tol*norm0) || (beta < eps)) {

    /*--- System is already solved ---*/

    if (master) cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
    residual = beta;
    return 0;
  }

  /*--- Normalize residual to get w_{0} (the negative sign is because w[0]
        holds the negative residual, as mentioned above). ---*/

  W[0] /= -beta;

  /*--- Initialize the RHS of the reduced system ---*/

  g[0] = beta;

  /*--- Set the norm to the initial residual value ---*/

  norm0 = beta;

  /*--- Output header information including initial residual ---*/

  unsigned long i = 0;
  if ((monitoring) && (master)) {
    WriteHeader("FGMRES", tol, beta);
    WriteHistory(i, beta/norm0);
  }

  /*---  Loop over all search directions ---*/

  for (i = 0; i < m; i++) {

    /*---  Check if solution has converged ---*/

    if (beta < tol*norm0) break;

    /*---  Precondition the CSysVector w[i] and store result in z[i] ---*/

    precond(W[i], Z[i]);

    /*---  Add to Krylov subspace ---*/

    mat_vec(Z[i], W[i+1]);

    /*---  Modified Gram-Schmidt orthogonalization ---*/

    ModGramSchmidt(i, H, W);

    /*---  Apply old Givens rotations to new column of the Hessenberg matrix then generate the
     new Givens rotation matrix and apply it to the last two elements of H[:][i] and g ---*/

    for (unsigned long k = 0; k < i; k++)
      ApplyGivens(sn[k], cs[k], H[k][i], H[k+1][i]);
    GenerateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
    ApplyGivens(sn[i], cs[i], g[i], g[i+1]);

    /*---  Set L2 norm of residual and check if solution has converged ---*/

    beta = fabs(g[i+1]);

    /*---  Output the relative residual if necessary ---*/

    if ((((monitoring) && (master)) && ((i+1) % 10 == 0)) && (master))
      WriteHistory(i+1, beta/norm0);
  }

  /*---  Solve the least-squares system and update solution ---*/

  SolveReduced(i, H, g, y);
  for (unsigned long k = 0; k < i; k++) {
    x.Plus_AX(y[k], Z[k]);
  }

  /*---  Recalculate final (neg.) residual (this should be optional) ---*/

  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {

    if (master) WriteFinalResidual("FGMRES", i, beta/norm0);

    mat_vec(x, W[0]);
    W[0] -= b;
    ScalarType res = W[0].norm();

    if (fabs(res - beta) > tol*10) {
      if (master) {
        WriteWarning(beta, res, tol);
      }
    }

  }

  residual = beta/norm0;
  return i;

}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::BCGSTAB_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                       const CMatrixVectorProduct<ScalarType> & mat_vec, const CPreconditioner<ScalarType> & precond,
                                                       ScalarType tol, unsigned long m, ScalarType & residual, bool monitoring, CConfig *config) const {

  const bool master = (SU2_MPI::GetRank() == MASTER_NODE) && (omp_get_thread_num() == 0);
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Check the subspace size ---*/

  if (m < 1) {
    SU2_OMP_MASTER
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet ---*/

  if (!bcg_ready) {
    SU2_OMP_MASTER
    {
      auto nVar = b.GetNVar();
      auto nBlk = b.GetNBlk();
      auto nBlkDomain = b.GetNBlkDomain();

      A_x.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r_0.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      p.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      v.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      z.Initialize(nBlk, nBlkDomain, nVar, nullptr);

      bcg_ready = true;
    }
    SU2_OMP_BARRIER
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

  mat_vec(x, A_x);
  r = b; r -= A_x;

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    norm_r = r.norm();
    norm0  = b.norm();
    if ((norm_r < tol*norm0) || (norm_r < eps)) {
      if (master) cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
      return 0;
    }

    /*--- Set the norm to the initial initial residual value ---*/

    norm0 = norm_r;

    /*--- Output header information including initial residual ---*/

    if ((monitoring) && (master)) {
      WriteHeader("BCGSTAB", tol, norm_r);
      WriteHistory(i, norm_r/norm0);
    }

  }

  /*--- Initialization ---*/

  ScalarType alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
  p = ScalarType(0.0); v = ScalarType(0.0); r_0 = r;

  /*--- Loop over all search directions ---*/

  for (i = 0; i < m; i++) {

    /*--- Compute rho_prime ---*/

    rho_prime = rho;

    /*--- Compute rho_i ---*/

    rho = r.dot(r_0);

    /*--- Compute beta ---*/

    beta = (rho / rho_prime) * (alpha /omega);

    /*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/

    ScalarType beta_omega = -beta*omega;
    p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
    p += r;

    /*--- Preconditioning step ---*/

    precond(p, z);
    mat_vec(z, v);

    /*--- Calculate step-length alpha ---*/

    ScalarType r_0_v = r_0.dot(v);
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

    omega = A_x.dot(r) / A_x.squaredNorm();

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
      if (((monitoring) && (master)) && ((i+1) % 10 == 0) && (master))
        WriteHistory(i+1, norm_r/norm0);

    }

  }

  /*--- Recalculate final residual (this should be optional) ---*/

  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {

    if (master) WriteFinalResidual("BCGSTAB", i, norm_r/norm0);

    mat_vec(x, A_x);
    r = b; r -= A_x;
    ScalarType true_res = r.norm();

    if ((fabs(true_res - norm_r) > tol*10.0) && (master)) {
      WriteWarning(norm_r, true_res, tol);
    }

  }

  residual = norm_r/norm0;
  return i;
}

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::Smoother_LinSolver(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x,
                                                        const CMatrixVectorProduct<ScalarType> & mat_vec, const CPreconditioner<ScalarType> & precond,
                                                        ScalarType tol, unsigned long m, ScalarType & residual, bool monitoring, CConfig *config) const {

  const bool master = (SU2_MPI::GetRank() == MASTER_NODE) && (omp_get_thread_num() == 0);
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Relaxation factor, see comments inside the loop over the smoothing iterations. ---*/
  ScalarType omega = SU2_TYPE::GetValue(config->GetLinear_Solver_Smoother_Relaxation());

  if (m < 1) {
    SU2_OMP_MASTER
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate vectors for residual (r), solution increment (z), and matrix-vector
   product (A_x), for the latter two this is done only on the first call to the method. ---*/

  if (!smooth_ready) {
    SU2_OMP_MASTER
    {
      auto nVar = b.GetNVar();
      auto nBlk = b.GetNBlk();
      auto nBlkDomain = b.GetNBlkDomain();

      A_x.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      z.Initialize(nBlk, nBlkDomain, nVar, nullptr);

      smooth_ready = true;
    }
    SU2_OMP_BARRIER
  }

  /*--- Compute the initial residual and check if the system is already solved (if in COMM_FULL mode). ---*/

  mat_vec(x, A_x);
  r = b; r -= A_x;

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    norm_r = r.norm();
    norm0  = b.norm();
    if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
      if (master) cout << "CSysSolve::Smoother_LinSolver(): system solved by initial guess." << endl;
      return 0;
    }

    /*--- Set the norm to the initial initial residual value. ---*/

    norm0 = norm_r;

    /*--- Output header information including initial residual. ---*/

    if ((monitoring) && (master)) {
      WriteHeader("Smoother", tol, norm_r);
      WriteHistory(i, norm_r/norm0);
    }

  }

  /*--- Smoothing Iterations ---*/

  for (i=0; i<m; i++) {

    /*--- Compute the solution increment (z), or "search" direction, by applying
     the preconditioner to the residual, i.e. z = M^{-1} * r. ---*/

    precond(r, z);

    /*--- The increment will be added to the solution (with relaxation omega)
     to update the residual, needed to compute the next increment, we get the
     product of the matrix with the direction (into A_x) and subtract from the
     current residual, the system is linear so this saves some computation
     compared to re-evaluating r = b-A*x. ---*/

    mat_vec(z, A_x);

    /*--- Update solution and residual with relaxation omega. Mathematically this
     is a modified Richardson iteration for the left-preconditioned system
     M^{-1}(b-A*x) which converges if ||I-w*M^{-1}*A|| < 1. Combining this method
     with a Gauss-Seidel preconditioner and w>1 is NOT equivalent to SOR. ---*/

    x.Plus_AX(omega, z);
    r.Plus_AX(-omega, A_x);

    /*--- Only compute the residuals in full communication mode. ---*/
    /*--- Check if solution has converged, else output the relative residual if necessary. ---*/

    if (config->GetComm_Level() == COMM_FULL) {
      norm_r = r.norm();
      if (norm_r < tol*norm0) break;
      if (((monitoring) && (master)) && ((i+1) % 5 == 0))
        WriteHistory(i+1, norm_r/norm0);
    }
  }

  if ((monitoring) && (master) && (config->GetComm_Level() == COMM_FULL)) {
    WriteFinalResidual("Smoother", i, norm_r/norm0);
  }

  residual = norm_r/norm0;
  return i;
}

template<>
void CSysSolve<su2double>::HandleTemporariesIn(const CSysVector<su2double> & LinSysRes, CSysVector<su2double> & LinSysSol) {

  /*--- When the type is the same the temporaties are not required ---*/
  /*--- Set the pointers ---*/
  LinSysRes_ptr = &LinSysRes;
  LinSysSol_ptr = &LinSysSol;
}

template<>
void CSysSolve<su2double>::HandleTemporariesOut(CSysVector<su2double> & LinSysSol) {

  /*--- When the type is the same the temporaties are not required ---*/
  /*--- Reset the pointers ---*/
  LinSysRes_ptr = nullptr;
  LinSysSol_ptr = nullptr;
}

#ifdef CODI_REVERSE_TYPE
template<>
void CSysSolve<passivedouble>::HandleTemporariesIn(const CSysVector<su2double> & LinSysRes, CSysVector<su2double> & LinSysSol) {

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
  LinSysRes_ptr = nullptr;
  LinSysSol_ptr = nullptr;
}
#endif

template<class ScalarType>
unsigned long CSysSolve<ScalarType>::Solve(CSysMatrix<ScalarType> & Jacobian, const CSysVector<su2double> & LinSysRes,
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
  unsigned long MaxIter, RestartIter;
  ScalarType SolverTol;
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

  /*--- Stop the recording for the linear solver ---*/

  bool TapeActive = NO;

  if (config->GetDiscrete_Adjoint()) {
#ifdef CODI_REVERSE_TYPE

    TapeActive = AD::globalTape.isActive();

    AD::StartExtFunc(false, false);

    AD::SetExtFuncIn(&LinSysRes[0], LinSysRes.GetLocSize());

    AD::StopRecording();
#endif
  }

  /*--- Create matrix-vector product, preconditioner, and solve the linear system ---*/

  HandleTemporariesIn(LinSysRes, LinSysSol);

  auto mat_vec = CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);
  CPreconditioner<ScalarType>* precond = nullptr;

  switch (KindPrecond) {
    case JACOBI:
      precond = new CJacobiPreconditioner<ScalarType>(Jacobian, geometry, config, false);
      break;
    case ILU:
      precond = new CILUPreconditioner<ScalarType>(Jacobian, geometry, config, false);
      break;
    case LU_SGS:
      precond = new CLU_SGSPreconditioner<ScalarType>(Jacobian, geometry, config);
      break;
    case LINELET:
      precond = new CLineletPreconditioner<ScalarType>(Jacobian, geometry, config);
      break;
    case PASTIX_ILU: case PASTIX_LU_P: case PASTIX_LDLT_P:
      precond = new CPastixPreconditioner<ScalarType>(Jacobian, geometry, config, KindPrecond, false);
      break;
    default:
      precond = new CJacobiPreconditioner<ScalarType>(Jacobian, geometry, config, false);
      break;
  }

  /*--- Start a thread-parallel section covering the preparation of the
   *    preconditioner and the solution of the linear solver.
   *    Beware of shared variables, i.e. defined outside the section or
   *    members of ANY class used therein, they should be treated as
   *    read-only or explicitly synchronized if written to. ---*/

  unsigned long IterLinSol = 0;

  SU2_OMP_PARALLEL
  {
    /*--- Build preconditioner in parallel. ---*/
    precond->Build();

    /*--- Thread-local variables. ---*/
    unsigned long iter = 0;
    ScalarType residual = 0.0, norm0 = 0.0;

    switch (KindSolver) {
      case BCGSTAB:
        iter = BCGSTAB_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual, ScreenOutput, config);
        break;
      case FGMRES:
        iter = FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual, ScreenOutput, config);
        break;
      case CONJUGATE_GRADIENT:
        iter = CG_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual, ScreenOutput, config);
        break;
      case RESTARTED_FGMRES:
        norm0 = LinSysRes_ptr->norm();
        while (iter < MaxIter) {
          /*--- Enforce a hard limit on total number of iterations ---*/
          unsigned long IterLimit = min(RestartIter, MaxIter-iter);
          iter += FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, IterLimit, residual, ScreenOutput, config);
          if ( residual < SolverTol*norm0 ) break;
        }
        break;
      case SMOOTHER:
        iter = Smoother_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual, ScreenOutput, config);
        break;
      case PASTIX_LDLT : case PASTIX_LU:
        Jacobian.BuildPastixPreconditioner(geometry, config, KindSolver);
        Jacobian.ComputePastixPreconditioner(*LinSysRes_ptr, *LinSysSol_ptr, geometry, config);
        iter = 1;
        break;
      default:
        SU2_MPI::Error("Unknown type of linear solver.",CURRENT_FUNCTION);
    }

    /*--- Only one thread modifies shared variables, synchronization
     *    is not required as we are exiting the parallel section. ---*/
    SU2_OMP_MASTER
    {
      IterLinSol = iter;
      Residual = residual;
    }

  } // end SU2_OMP_PARALLEL

  delete precond;

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
      case PASTIX_ILU: case PASTIX_LU_P: case PASTIX_LDLT_P:
        Jacobian.BuildPastixPreconditioner(geometry, config, KindPrecond, RequiresTranspose);
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
unsigned long CSysSolve<ScalarType>::Solve_b(CSysMatrix<ScalarType> & Jacobian, const CSysVector<su2double> & LinSysRes,
                                             CSysVector<su2double> & LinSysSol, CGeometry *geometry, CConfig *config) {
#ifdef CODI_REVERSE_TYPE

  unsigned short KindSolver, KindPrecond;
  unsigned long MaxIter, RestartIter, IterLinSol = 0;
  ScalarType SolverTol, Norm0 = 0.0;
  bool ScreenOutput, RequiresTranspose = !mesh_deform; // jacobian is symmetric

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

  CPreconditioner<ScalarType>* precond  = nullptr;

  switch(KindPrecond) {
    case ILU:
      precond = new CILUPreconditioner<ScalarType>(Jacobian, geometry, config, RequiresTranspose);
      break;
    case JACOBI:
      precond = new CJacobiPreconditioner<ScalarType>(Jacobian, geometry, config, RequiresTranspose);
      break;
    case PASTIX_ILU: case PASTIX_LU_P: case PASTIX_LDLT_P:
      precond = new CPastixPreconditioner<ScalarType>(Jacobian, geometry, config, KindPrecond, RequiresTranspose);
      break;
  }

  auto mat_vec = CSysMatrixVectorProductTransposed<ScalarType>(Jacobian, geometry, config);

  /*--- Solve the system ---*/

  HandleTemporariesIn(LinSysRes, LinSysSol);

  switch(KindSolver) {
    case FGMRES:
      IterLinSol = FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol , MaxIter, Residual, ScreenOutput, config);
      break;
    case BCGSTAB:
      IterLinSol = BCGSTAB_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol , MaxIter, Residual, ScreenOutput, config);
      break;
    case CONJUGATE_GRADIENT:
      IterLinSol = CG_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, Residual, ScreenOutput, config);
      break;
    case RESTARTED_FGMRES:
      IterLinSol = 0;
      Norm0 = LinSysRes_ptr->norm();
      while (IterLinSol < MaxIter) {
        /*--- Enforce a hard limit on total number of iterations ---*/
        unsigned long IterLimit = min(RestartIter, MaxIter-IterLinSol);
        IterLinSol += FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol , IterLimit, Residual, ScreenOutput, config);
        if ( Residual < SolverTol*Norm0 ) break;
      }
      break;
    case PASTIX_LDLT : case PASTIX_LU:
      Jacobian.BuildPastixPreconditioner(geometry, config, KindSolver, RequiresTranspose);
      Jacobian.ComputePastixPreconditioner(*LinSysRes_ptr, *LinSysSol_ptr, geometry, config);
      IterLinSol = 1;
      break;
    default:
      SU2_MPI::Error("The specified linear solver is not yet implemented for the discrete adjoint method.", CURRENT_FUNCTION);
      break;
  }

  HandleTemporariesOut(LinSysSol);

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
