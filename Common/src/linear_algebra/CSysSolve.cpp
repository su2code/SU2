/*!
 * \file CSysSolve.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios, T. Economon, P. Gomes
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

#include "../../include/linear_algebra/CSysSolve.hpp"
#include "../../include/linear_algebra/CSysSolve_b.hpp"
#include "../../include/parallelization/omp_structure.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../include/linear_algebra/CPreconditioner.hpp"

#include <complex>
#include <iostream>
#include <limits>
#include <type_traits>
#include <memory>

namespace {
/*!
 * \brief Epsilon used in CSysSolve depending on datatype to
 * decide if the linear system is already solved.
 */
template <class T>
constexpr T linSolEpsilon() {
  return numeric_limits<passivedouble>::epsilon();
}
template <>
constexpr float linSolEpsilon<float>() {
  return 1e-14;
}

/*--- Computes v = vs * ws or v += vs * ws with unrolling of up to 4 iterations. ---*/
template <class ScalarType, class Weights, class Vectors>
void LinearCombination(const unsigned long n, const Vectors& vs, const Weights& ws, CSysVector<ScalarType>& v,
                       bool inc = false) {
  if (n == 0) {
    if (!inc) v = ScalarType{};
    return;
  }
  for (unsigned long i = 0; i < n; i += 4) {
    switch (n - i) {
      case 1: {
        if (inc) {
          v += vs(i) * ws(i);
        } else {
          v = vs(i) * ws(i);
        }
      } break;
      case 2: {
        if (inc) {
          v += vs(i) * ws(i) + vs(i + 1) * ws(i + 1);
        } else {
          v = vs(i) * ws(i) + vs(i + 1) * ws(i + 1);
        }
      } break;
      case 3: {
        if (inc) {
          v += vs(i) * ws(i) + vs(i + 1) * ws(i + 1) + vs(i + 2) * ws(i + 2);
        } else {
          v = vs(i) * ws(i) + vs(i + 1) * ws(i + 1) + vs(i + 2) * ws(i + 2);
        }
      } break;
      default: {
        if (inc) {
          v += vs(i) * ws(i) + vs(i + 1) * ws(i + 1) + vs(i + 2) * ws(i + 2) + vs(i + 3) * ws(i + 3);
        } else {
          v = vs(i) * ws(i) + vs(i + 1) * ws(i + 1) + vs(i + 2) * ws(i + 2) + vs(i + 3) * ws(i + 3);
        }
      }
    }
    inc = true;
  }
}

/*--- Overload to handle a vector of CSysVector directly. ---*/
template <class ScalarType, class Weights>
void LinearCombination(const unsigned long n, const std::vector<CSysVector<ScalarType>>& vs, const Weights& ws,
                       CSysVector<ScalarType>& v, bool inc = false) {
  LinearCombination(
      n, [&vs](auto i) -> auto& { return vs[i]; }, ws, v, inc);
}

/*--- Overload to handle a std::vector<T> of weights directly. ---*/
template <class ScalarType, class Vectors>
void LinearCombination(const unsigned long n, const Vectors& vs, const std::vector<ScalarType>& ws,
                       CSysVector<ScalarType>& v, bool inc = false) {
  LinearCombination(
      n, vs, [&ws](auto i) { return ws[i]; }, v, inc);
}
}  // namespace

template <class ScalarType>
CSysSolve<ScalarType>::CSysSolve(LINEAR_SOLVER_MODE linear_solver_mode)
    : eps(linSolEpsilon<ScalarType>()),
      lin_sol_mode(linear_solver_mode),
      cg_ready(false),
      bcg_ready(false),
      smooth_ready(false),
      LinSysSol_ptr(nullptr),
      LinSysRes_ptr(nullptr) {}

template <class ScalarType>
void CSysSolve<ScalarType>::ApplyGivens(ScalarType s, ScalarType c, ScalarType& h1, ScalarType& h2) const {
  ScalarType temp = c * h1 + s * h2;
  h2 = c * h2 - s * h1;
  h1 = temp;
}

template <class ScalarType>
void CSysSolve<ScalarType>::GenerateGivens(ScalarType& dx, ScalarType& dy, ScalarType& s, ScalarType& c) const {
  if ((dx == 0.0) && (dy == 0.0)) {
    c = 1.0;
    s = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    ScalarType tmp = dx / dy;
    dx = sqrt(1.0 + tmp * tmp);
    s = Sign(1.0 / dx, dy);
    c = tmp * s;
  } else if (fabs(dy) <= fabs(dx)) {
    ScalarType tmp = dy / dx;
    dy = sqrt(1.0 + tmp * tmp);
    c = Sign(1.0 / dy, dx);
    s = tmp * c;
  } else {
    // dx and/or dy must be invalid
    dx = 0.0;
    dy = 0.0;
    c = 1.0;
    s = 0.0;
  }
  dx = fabs(dx * dy);
  dy = 0.0;
}

template <class ScalarType>
void CSysSolve<ScalarType>::SolveReduced(int n, const su2matrix<ScalarType>& Hsbg, const su2vector<ScalarType>& rhs,
                                         su2vector<ScalarType>& x) const {
  // initialize...
  for (int i = 0; i < n; i++) x[i] = rhs[i];
  // ... and backsolve
  for (int i = n - 1; i >= 0; i--) {
    x[i] /= Hsbg(i, i);
    for (int j = i - 1; j >= 0; j--) {
      x[j] -= Hsbg(j, i) * x[i];
    }
  }
}

template <class ScalarType>
void CSysSolve<ScalarType>::ModGramSchmidt(bool shared_hsbg, int i, su2matrix<ScalarType>& Hsbg,
                                           vector<CSysVector<ScalarType>>& w) const {
  const auto thread = omp_get_thread_num();

  /*--- If Hsbg is shared by multiple threads calling this function, only one
   * thread can write into it. If Hsbg is private, all threads need to write. ---*/

  auto SetHsbg = [&](int row, int col, const ScalarType& value) {
    if (!shared_hsbg || thread == 0) {
      Hsbg(row, col) = value;
    }
  };

  /*--- Parameter for reorthonormalization ---*/

  const ScalarType reorth = 0.98;

  /*--- Get the norm of the vector being orthogonalized, and find the
  threshold for re-orthogonalization ---*/

  ScalarType nrm = w[i + 1].squaredNorm();
  ScalarType thr = nrm * reorth;

  /*--- The norm of w[i+1] < 0.0 or w[i+1] = NaN ---*/

  if ((nrm <= 0.0) || (nrm != nrm)) {
    /*--- nrm is the result of a dot product, communications are implicitly handled. ---*/
    SU2_MPI::Error("FGMRES orthogonalization failed, linear solver diverged.", CURRENT_FUNCTION);
  }

  /*--- Begin main Gram-Schmidt loop ---*/

  for (int k = 0; k < i + 1; k++) {
    ScalarType prod = w[i + 1].dot(w[k]);
    ScalarType h_ki = prod;
    w[i + 1] -= prod * w[k];

    /*--- Check if reorthogonalization is necessary ---*/

    if (prod * prod > thr) {
      prod = w[i + 1].dot(w[k]);
      h_ki += prod;
      w[i + 1] -= prod * w[k];
    }
    SetHsbg(k, i, h_ki);

    /*--- Update the norm and check its size ---*/

    nrm -= pow(h_ki, 2);
    nrm = max<ScalarType>(nrm, 0.0);
    thr = nrm * reorth;
  }

  /*--- Test the resulting vector ---*/

  nrm = w[i + 1].norm();
  SetHsbg(i + 1, i, nrm);

  /*--- Scale the resulting vector ---*/

  w[i + 1] /= nrm;
}

template <class ScalarType>
void CSysSolve<ScalarType>::WriteHeader(const string& solver, ScalarType restol, ScalarType resinit) const {
  cout << "\n# " << solver << " residual history\n";
  cout << "# Residual tolerance target = " << restol << "\n";
  cout << "# Initial residual norm     = " << resinit << endl;
}

template <class ScalarType>
void CSysSolve<ScalarType>::WriteHistory(unsigned long iter, ScalarType res) const {
  cout << "     " << iter << "     " << res << endl;
}

template <class ScalarType>
void CSysSolve<ScalarType>::WriteFinalResidual(const string& solver, unsigned long iter, ScalarType res) const {
  cout << "# " << solver << " final (true) residual:\n";
  cout << "# Iteration = " << iter << ": |res|/|res0| = " << res << ".\n" << endl;
}

template <class ScalarType>
void CSysSolve<ScalarType>::WriteWarning(ScalarType res_calc, ScalarType res_true, ScalarType tol) const {
  cout << "# WARNING:\n";
  cout << "# true residual norm and calculated residual norm do not agree.\n";
  cout << "# true_res = " << res_true << ", calc_res = " << res_calc << ", tol = " << tol * 10 << ".\n";
  cout << "# true_res - calc_res = " << res_true - res_calc << endl;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::CG_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                  const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                  const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                  unsigned long m, ScalarType& residual, bool monitoring,
                                                  const CConfig* config) const {
  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Check the subspace size ---*/

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet, only one thread can
   *    do this since the working vectors are shared. ---*/

  if (!cg_ready) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      auto nVar = b.GetNVar();
      auto nBlk = b.GetNBlk();
      auto nBlkDomain = b.GetNBlkDomain();

      A_x.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      z.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      p.Initialize(nBlk, nBlkDomain, nVar, nullptr);

      cg_ready = true;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

  if (!xIsZero) {
    mat_vec(x, A_x);
    r = b - A_x;
  } else {
    r = b;
  }

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
    norm_r = r.norm();
    norm0 = b.norm();

    /*--- Set the norm to the initial initial residual value ---*/

    if (tol_type == LinearToleranceType::RELATIVE) norm0 = norm_r;

    if ((norm_r < tol * norm0) || (norm_r < eps)) {
      if (masterRank && (lin_sol_mode != LINEAR_SOLVER_MODE::MESH_DEFORM)) {
        SU2_OMP_MASTER
        cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
        END_SU2_OMP_MASTER
      }
      return 0;
    }

    /*--- Output header information including initial residual ---*/

    if (monitoring && masterRank) {
      SU2_OMP_MASTER {
        WriteHeader("CG", tol, norm_r);
        WriteHistory(i, norm_r / norm0);
      }
      END_SU2_OMP_MASTER
    }
  }

  precond(r, z);
  p = z;
  ScalarType r_dot_z = r.dot(z);

  /*---  Loop over all search directions ---*/

  for (i = 0; i < m; i++) {
    /*--- Apply matrix to p to build Krylov subspace ---*/

    mat_vec(p, A_x);

    /*--- Calculate step-length alpha ---*/

    ScalarType alpha = r_dot_z / A_x.dot(p);

    /*--- Update solution and residual: ---*/

    x += alpha * p;
    r -= alpha * A_x;

    /*--- Only compute the residuals in full communication mode. ---*/

    if (config->GetComm_Level() == COMM_FULL) {
      /*--- Check if solution has converged, else output the relative residual if necessary ---*/

      norm_r = r.norm();
      if (norm_r < tol * norm0) break;
      if (((monitoring) && (masterRank)) && ((i + 1) % monitorFreq == 0)) {
        SU2_OMP_MASTER
        WriteHistory(i + 1, norm_r / norm0);
        END_SU2_OMP_MASTER
      }
    }

    precond(r, z);

    /*--- Calculate Gram-Schmidt coefficient, beta = (r_{i+1}, z_{i+1}) / (r_{i}, z_{i}) ---*/

    ScalarType beta = r_dot_z;
    r_dot_z = r.dot(z);
    beta = r_dot_z / beta;

    /*--- Gram-Schmidt orthogonalization. ---*/

    p = beta * p + z;
  }

  /*--- Recalculate final residual (this should be optional) ---*/

  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("CG", i, norm_r / norm0);
      END_SU2_OMP_MASTER
    }

    if (recomputeRes) {
      mat_vec(x, A_x);
      r = b - A_x;
      ScalarType true_res = r.norm();

      if (fabs(true_res - norm_r) > tol * 10.0) {
        if (masterRank) {
          SU2_OMP_MASTER
          WriteWarning(norm_r, true_res, tol);
          END_SU2_OMP_MASTER
        }
      }
    }
  }

  residual = norm_r / norm0;
  return i;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                      const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                      const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                      unsigned long m, ScalarType& residual, bool monitoring,
                                                      const CConfig* config) const {
  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  const bool flexible = !precond.IsIdentity();
  /*--- If we call the solver outside of a parallel region, but the number of threads allows,
   * we still want to parallelize some of the expensive operations. ---*/
  const bool nestedParallel = !omp_in_parallel() && omp_get_max_threads() > 1;

  /*---  Check the subspace size ---*/

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  if (m > 1000) {
    SU2_MPI::Error("FGMRES subspace is too large (>1000).", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet ---*/

  if (V.size() <= m || (flexible && Z.size() <= m)) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      V.resize(m + 1);
      for (auto& w : V) w.Initialize(x.GetNBlk(), x.GetNBlkDomain(), x.GetNVar(), nullptr);
      if (flexible) {
        Z.resize(m + 1);
        for (auto& z : Z) z.Initialize(x.GetNBlk(), x.GetNBlkDomain(), x.GetNVar(), nullptr);
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Define various arrays. In parallel, each thread of each rank has and works
   on its own thread, since calculations on these arrays are based on dot products
   (reduced across all threads and ranks) all threads do the same computations. ---*/

  su2vector<ScalarType> g(m + 1), sn(m + 1), cs(m + 1), y(m);
  g = ScalarType(0);
  sn = ScalarType(0);
  cs = ScalarType(0);
  y = ScalarType(0);
  su2matrix<ScalarType> H(m + 1, m);
  H = ScalarType(0);

  /*--- Calculate the norm of the rhs vector. ---*/

  ScalarType norm0 = b.norm();

  /*--- Calculate the initial residual (actually the negative residual) and compute its norm. ---*/

  if (!xIsZero) {
    mat_vec(x, V[0]);
    V[0] -= b;
  } else {
    V[0] = -b;
  }

  ScalarType beta = xIsZero ? norm0 : V[0].norm();

  /*--- FGMRES shares V and Z with FGCRODR, reset deflation because we're breaking
   * relations between subspaces that FGCRODR relies on for recycling. ---*/

  SU2_OMP_MASTER
  ResetDeflation();
  END_SU2_OMP_MASTER

  /*--- Set the norm to the initial initial residual value ---*/

  if (tol_type == LinearToleranceType::RELATIVE) norm0 = beta;

  if (beta < tol * norm0 || beta < eps) {
    /*--- System is already solved ---*/

    if (masterRank) {
      SU2_OMP_MASTER
      cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
      END_SU2_OMP_MASTER
    }
    residual = beta;
    return 0;
  }

  /*--- Normalize residual to get w_{0} (the negative sign is because w[0]
        holds the negative residual, as mentioned above). ---*/

  V[0] /= -beta;

  /*--- Initialize the RHS of the reduced system ---*/

  g[0] = beta;

  /*--- Output header information including initial residual ---*/

  unsigned long i = 0;
  if (monitoring && masterRank) {
    SU2_OMP_MASTER {
      WriteHeader("FGMRES", tol, beta);
      WriteHistory(i, beta / norm0);
    }
    END_SU2_OMP_MASTER
  }

  /*---  Loop over all search directions ---*/

  for (i = 0; i < m; i++) {
    /*---  Check if solution has converged ---*/

    if (beta < tol * norm0) break;

    if (flexible) {
      /*---  Precondition the CSysVector v[i] and store result in z[i] ---*/

      precond(V[i], Z[i]);

      /*---  Add to Krylov subspace ---*/

      mat_vec(Z[i], V[i + 1]);
    } else {
      mat_vec(V[i], V[i + 1]);
    }

    /*---  Modified Gram-Schmidt orthogonalization ---*/

    if (nestedParallel) {
      /*--- "omp parallel if" does not work well here ---*/
      SU2_OMP_PARALLEL
      ModGramSchmidt(true, i, H, V);
      END_SU2_OMP_PARALLEL
    } else {
      ModGramSchmidt(false, i, H, V);
    }

    /*---  Apply old Givens rotations to new column of the Hessenberg matrix then generate the
     new Givens rotation matrix and apply it to the last two elements of H[:][i] and g ---*/

    for (unsigned long k = 0; k < i; k++) ApplyGivens(sn[k], cs[k], H(k, i), H(k + 1, i));
    GenerateGivens(H(i, i), H(i + 1, i), sn[i], cs[i]);
    ApplyGivens(sn[i], cs[i], g[i], g[i + 1]);

    /*---  Set L2 norm of residual and check if solution has converged ---*/

    beta = fabs(g[i + 1]);

    /*---  Output the relative residual if necessary ---*/

    if (monitoring && masterRank && ((i + 1) % monitorFreq == 0)) {
      SU2_OMP_MASTER
      WriteHistory(i + 1, beta / norm0);
      END_SU2_OMP_MASTER
    }
  }

  /*---  Solve the least-squares system and update solution ---*/

  SolveReduced(i, H, g, y);

  const auto& basis = flexible ? Z : V;

  if (nestedParallel) {
    SU2_OMP_PARALLEL
    LinearCombination(i, basis, y, x, true);
    END_SU2_OMP_PARALLEL
  } else {
    LinearCombination(i, basis, y, x, true);
  }

  /*---  Recalculate final (neg.) residual (this should be optional) ---*/

  if (monitoring && config->GetComm_Level() == COMM_FULL) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("FGMRES", i, beta / norm0);
      END_SU2_OMP_MASTER
    }

    if (recomputeRes) {
      mat_vec(x, V[0]);
      V[0] -= b;
      ScalarType res = V[0].norm();

      if (fabs(res - beta) > tol * 10) {
        if (masterRank) {
          SU2_OMP_MASTER
          WriteWarning(beta, res, tol);
          END_SU2_OMP_MASTER
        }
      }
    }
  }

  residual = beta / norm0;
  return i;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::RFGMRES_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                       const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                       const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                       unsigned long MaxIter, ScalarType& residual, bool monitoring,
                                                       const CConfig* config) const {
  const auto restartIter = config->GetLinear_Solver_Restart_Frequency();

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    xIsZero = false;
    tol_type = LinearToleranceType::ABSOLUTE;
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  for (auto totalIter = 0ul; totalIter < MaxIter;) {
    /*--- Enforce a hard limit on total number of iterations ---*/
    auto iterLimit = min(restartIter, MaxIter - totalIter);
    auto iter = FGMRES_LinSolver(b, x, mat_vec, precond, tol, iterLimit, residual, monitoring, config);
    totalIter += iter;
    if (residual <= tol || iter < iterLimit) return totalIter;
  }
  return 0;
}

template <class ScalarType>
template <class Dummy>
unsigned long CSysSolve<ScalarType>::FGCRODR_LinSolverImpl(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                           const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                           const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                           unsigned long max_iter, ScalarType& residual,
                                                           bool monitoring, const CConfig* config,
                                                           FgcrodrMode mode) const {
  using EigenMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
  using EigenVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

  auto m = min(config->GetLinear_Solver_Restart_Frequency(), max_iter);
  const auto deflation = min(config->GetLinear_Solver_Restart_Deflation(), m - 1);

  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  /*--- If we call the solver outside of a parallel region, but the number of threads allows,
   * we still want to parallelize some of the expensive operations. ---*/
  const bool nestedParallel = !omp_in_parallel() && omp_get_max_threads() > 1;

  /*--- Check the subspace size. ---*/

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver restart iterations must be greater than 0.", CURRENT_FUNCTION);
  }
  if (m > 200) {
    SU2_MPI::Error("FGCRODR subspace is too large (> 200).", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet. ---*/

  if (V.size() < m + 1 || W.size() < deflation + 1) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      const auto nVar = b.GetNVar();
      const auto nBlk = b.GetNBlk();
      const auto nBlkDomain = b.GetNBlkDomain();

      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      V.resize(m + 1);
      for (auto& v : V) v.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      Z.resize(m);
      for (auto& z : Z) z.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      W.resize(deflation + 1);
      for (auto& w : W) w.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      T.resize(deflation + 1);
      for (auto& t : T) t.Initialize(nBlk, nBlkDomain, nVar, nullptr);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Hessenberg matrix. See FGMRES for parallelization notes. ---*/

  su2matrix<ScalarType> H(m + 1, m);
  EigenMatrix Heigen(m + 1, m);

  /*--- Calculate the norm of the rhs vector. ---*/

  ScalarType norm0 = b.norm();

  /*--- Calculate the initial residual and compute its norm. ---*/

  if (!xIsZero) {
    mat_vec(x, Z.back());
    r = b - Z.back();
  } else {
    r = b;
  }

  /*--- Rebuild Z, V, and W for the new matrix if necessary.
   * Q * R = A * Z
   * V = Q = A * (Z * R^-1)
   * Such that the first part of the Arnoldi remains:
   * A * (Zk * R^-1) = Vk * I ---*/
  if (k > 0) {
    EigenMatrix R = EigenMatrix::Zero(k, k);
    EigenVector vr = EigenVector::Zero(k);
    for (auto j = 0ul; j < k; ++j) {
      if (mode != FgcrodrMode::SAME_MAT) {
        /*--- When k = 0, Z = M(W), we could keep that property but it is not
         * critical and so we choose to save the cost of precond(W[j], Z[j]); ---*/
        mat_vec(Z[j], V[j]);

        for (auto i = 0ul; i < j; ++i) {
          R(i, j) = V[i].dot(V[j]);
          V[j] -= R(i, j) * V[i];
        }
        R(j, j) = V[j].norm();
        V[j] /= R(j, j);
      }
      /*--- Make r orthogonal to the rebuilt V so we can proceed with the usual Arnoldi process. ---*/
      vr(j) = r.dot(V[j]);
    }
    LinearCombination(k, V, -vr, r, true);

    /*--- Apply R^-1 to Z and W and update x accordingly. R is uppper triangular,
     * so we loop backwards to compute the products in-place. ---*/
    if (mode != FgcrodrMode::SAME_MAT) {
      EigenMatrix invR = R.template triangularView<Eigen::Upper>().solve(EigenMatrix::Identity(k, k));
      for (auto j = k - 1;; --j) {
        for (auto* basis : {&W, &Z}) {
          auto reversed = [&](auto i) -> const auto& { return (*basis)[j - i]; };
          LinearCombination(
              j + 1, reversed, [&](auto i) { return invR(j - i, j); }, (*basis)[j]);
        }
        if (j == 0) break;  // j is unsigned, avoid underflow.
      }
    }
    LinearCombination(k, Z, vr, x, true);
  }
  ScalarType rNorm = r.norm();
  auto iter = k;

  /*--- Set the norm to the initial initial residual value for relative tolerance. ---*/

  if (tol_type == LinearToleranceType::RELATIVE) norm0 = rNorm;
  residual = rNorm / norm0;

  if (residual < tol || rNorm < eps) {
    /*--- The system is already solved. ---*/
    if (masterRank) {
      SU2_OMP_MASTER
      cout << "CSysSolve::FGCRODR(): system solved by initial guess." << endl;
      END_SU2_OMP_MASTER
    }
    return iter;
  }

  /*--- Output header information including initial residual. ---*/

  if (monitoring && masterRank) {
    SU2_OMP_MASTER {
      WriteHeader("FGCRODR", tol, rNorm);
      WriteHistory(iter, residual);
    }
    END_SU2_OMP_MASTER
  }

  while (iter < max_iter) {
    /*--- Initial direction. ---*/
    V[k] = r / rNorm;

    /*--- Initialize the H matrix and RHS of the reduced system. ---*/
    H = ScalarType(0);
    Heigen = EigenMatrix::Zero(m + 1, m);
    for (auto j = 0ul; j < k; ++j) {
      Heigen(j, j) = H(j, j) = 1;
    }
    EigenVector c = EigenVector::Zero(m + 1), rls = EigenVector::Zero(m + 1), y;
    /*--- c = r_norm * e1 when k=0 is a special case of c = V' * r. ---*/
    c(k) = rNorm;

    /*--- Inner loop to build the generalized Arnoldi relation, AZ = VH. ---*/
    bool converged = false;
    for (auto j = k; j < m; ++j) {
      ++iter;
      precond(V[j], Z[j]);
      mat_vec(Z[j], V[j + 1]);

      if (nestedParallel) {
        /*--- "omp parallel if" does not work well here ---*/
        SU2_OMP_PARALLEL
        ModGramSchmidt(true, j, H, V);
        END_SU2_OMP_PARALLEL
      } else {
        ModGramSchmidt(false, j, H, V);
      }

      /*--- Solve the reduced system. We do not use Given's rotations to factor H in place
       * because we need the Arnoldi relation to hold to compute the Ritz values. ---*/
      const auto mj = j + 1;
      auto Hj = Heigen.leftCols(mj);

      for (auto i = 0ul; i <= mj; ++i) {
        Hj(i, j) = H(i, j);
      }
      y = Hj.colPivHouseholderQr().solve(c);

      /*--- Residual of the least squares problem, its norm is equivalent to ||b - Ax||. ---*/
      rls = c - Hj * y;
      rNorm = rls.norm();

      residual = rNorm / norm0;

      if (monitoring && masterRank) {
        SU2_OMP_MASTER
        WriteHistory(iter, residual);
        END_SU2_OMP_MASTER
      }
      if ((iter >= max_iter || residual < tol) && j >= deflation) {
        /*--- Adjust m if we stop early. ---*/
        converged = true;
        m = mj;
        break;
      }
    }

    /*--- We don't store the part of W that is equal to V explicitly,
     * W(:, k:m) = V(:, k:m). ---*/

    /*--- Update the solution and residual. The latter is only required if we restart. ---*/

    LinearCombination(m, Z, y, x, true);
    if (!converged) LinearCombination(m + 1, V, rls, r);

    /*--- Update deflation vectors. ---*/

    if (deflation == 0) {
      if (converged) break;
      continue;
    }

    /*--- Compute Ritz values and keep the ones with the smallest real part. ---*/

    EigenMatrix VW = EigenMatrix::Identity(m + 1, m);
    for (auto i = 0ul; i <= m; ++i) {
      for (auto j = 0ul; j < k; ++j) {
        // TODO(pedro): There are clever ways to avoid this multiplication, or at least use BLAS.
        VW(i, j) = V[i].dot(W[j]);
      }
    }
    const auto Hm = Heigen.topLeftCorner(m + 1, m);
    EigenMatrix HTVW = Hm.transpose() * VW;

    /*--- If the "B" matrix in the generalized eigenvalue problem is not invertible we reset. ---*/
    if (Eigen::ColPivHouseholderQR<EigenMatrix> qr(HTVW); !qr.isInvertible()) {
      if (masterRank) {
        SU2_OMP_MASTER
        cout << "WARNING: (VH)^T W in FGCRODR is not invertible.\n";
        END_SU2_OMP_MASTER
      }
      SU2_OMP_SAFE_GLOBAL_ACCESS(ResetDeflation();)
      if (converged) break;
      continue;
    }
    EigenMatrix HTH = Hm.transpose() * Hm;
    Eigen::GeneralizedEigenSolver<EigenMatrix> ges(HTH, HTVW);
    const auto& lambda = ges.eigenvalues();

    std::vector<int> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&lambda](int i, int j) { return fabs(std::real(lambda(i))) < fabs(std::real(lambda(j))); });

    EigenMatrix P(m, deflation + 1);
    auto k_new = 0ul;
    for (auto i = 0ul; i < m; ++i) {
      const auto j = order[i];

      /*--- Skip conjugate pairs because we split complex vectors into real and imag. ---*/
      if (i > 0 && abs(lambda(j) - std::conj(lambda(order[i - 1]))) / abs(lambda(j)) < 1e-3) {
        continue;
      }
      if (monitoring && masterRank && config->GetComm_Level() == COMM_FULL) {
        SU2_OMP_MASTER
        cout << "     FGCRODR Ritz value #" << i << ": " << lambda(j) << "\n";
        END_SU2_OMP_MASTER
      }
      P.col(k_new++) = ges.eigenvectors().col(j).real();

      if (fabs(std::imag(lambda(j))) > 1e-2 * fabs(std::real(lambda(j)))) {
        P.col(k_new++) = ges.eigenvectors().col(j).imag();
      }
      if (k_new >= deflation) break;
    }
    P.conservativeResize(m, k_new);

    /*--- Modify the Krylov basis vectors using P.
     * A Z P = V H P <=> A Z P = V Q R <=> A Z P R^-1 = V Q <=> A Zk = Vk.
     * W is updated the same way as Z since Z = M(W). ---*/

    EigenMatrix HP = Hm * P;
    Eigen::HouseholderQR<EigenMatrix> qr(HP);
    EigenMatrix Q = qr.householderQ() * EigenMatrix::Identity(m + 1, k_new);
    auto R = qr.matrixQR().topRows(k_new).template triangularView<Eigen::Upper>();

    EigenMatrix PinvR = P * (R.solve(EigenMatrix::Identity(k_new, k_new)));

    auto modify = [&](const EigenMatrix& mod, const auto& basis) {
      for (auto j = 0ul; j < k_new; ++j) {
        LinearCombination(mod.rows(), basis, mod.col(j), T[j]);
      }
    };
    modify(
        PinvR, [&](auto i) -> auto& { return i < k ? W[i] : V[i]; });

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      /*--- T and W are the same size, so we can swap them. ---*/
      std::swap(T, W);
      k = k_new;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- For other basis we need to swap vectors one by one. ---*/
    auto update = [&](auto& basis) {
      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
        for (auto j = 0ul; j < k_new; ++j) T[j].swap(basis[j]);
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS
    };

    modify(PinvR, Z);
    update(Z);

    /*--- Update V only if necessary. ---*/
    if (!converged || mode == FgcrodrMode::SAME_MAT) {
      modify(Q, V);
      update(V);
    }
    if (converged) break;
  }

  /*--- Recalculate final (neg.) residual (this should be optional). ---*/

  if (monitoring && config->GetComm_Level() == COMM_FULL) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("FGCRODR", iter, residual);
      END_SU2_OMP_MASTER
    }
    if (recomputeRes) {
      mat_vec(x, Z.back());
      Z.back() -= b;
      ScalarType res = Z.back().norm();

      if (fabs(res - rNorm) > tol * 10) {
        if (masterRank) {
          SU2_OMP_MASTER
          WriteWarning(rNorm, res, tol);
          END_SU2_OMP_MASTER
        }
      }
    }
  }
  return iter;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::FGCRODR_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                       const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                       const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                       unsigned long max_iter, ScalarType& residual, bool monitoring,
                                                       const CConfig* config, FgcrodrMode mode) const {
  if constexpr (std::is_same_v<ScalarType, float> || std::is_same_v<ScalarType, double>) {
    return FGCRODR_LinSolverImpl<>(b, x, mat_vec, precond, tol, max_iter, residual, monitoring, config, mode);
  } else {
    return RFGMRES_LinSolver(b, x, mat_vec, precond, tol, max_iter, residual, monitoring, config);
  }
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::BCGSTAB_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                       const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                       const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                       unsigned long m, ScalarType& residual, bool monitoring,
                                                       const CConfig* config) const {
  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Check the subspace size ---*/

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet ---*/

  if (!bcg_ready) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
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
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

  if (!xIsZero) {
    mat_vec(x, A_x);
    r = b - A_x;
  } else {
    r = b;
  }

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
    norm_r = r.norm();
    norm0 = b.norm();

    /*--- Set the norm to the initial initial residual value ---*/

    if (tol_type == LinearToleranceType::RELATIVE) norm0 = norm_r;

    if ((norm_r < tol * norm0) || (norm_r < eps)) {
      if (masterRank) {
        SU2_OMP_MASTER
        cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
        END_SU2_OMP_MASTER
      }
      return 0;
    }

    /*--- Output header information including initial residual ---*/

    if ((monitoring) && (masterRank)) {
      SU2_OMP_MASTER {
        WriteHeader("BCGSTAB", tol, norm_r);
        WriteHistory(i, norm_r / norm0);
      }
      END_SU2_OMP_MASTER
    }
  }

  /*--- Initialization ---*/

  ScalarType alpha = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
  p = ScalarType(0.0);
  v = ScalarType(0.0);
  r_0 = r;

  /*--- Loop over all search directions ---*/

  for (i = 0; i < m; i++) {
    /*--- Compute rho_prime ---*/

    rho_prime = rho;

    /*--- Compute rho_i ---*/

    rho = r.dot(r_0);

    /*--- Compute beta ---*/

    ScalarType beta = (rho / rho_prime) * (alpha / omega);

    /*--- Update p ---*/

    p = beta * (p - omega * v) + r;

    /*--- Preconditioning step ---*/

    precond(p, z);
    mat_vec(z, v);

    /*--- Calculate step-length alpha ---*/

    ScalarType r_0_v = r_0.dot(v);
    alpha = rho / r_0_v;

    /*--- Update solution and residual ---*/

    x += alpha * z;
    r -= alpha * v;

    /*--- Preconditioning step ---*/

    precond(r, z);
    mat_vec(z, A_x);

    /*--- Calculate step-length omega, avoid division by 0. ---*/

    omega = A_x.squaredNorm();
    if (omega == ScalarType(0)) break;
    omega = A_x.dot(r) / omega;

    /*--- Update solution and residual ---*/

    x += omega * z;
    r -= omega * A_x;

    /*--- Only compute the residuals in full communication mode. ---*/

    if (config->GetComm_Level() == COMM_FULL) {
      /*--- Check if solution has converged, else output the relative residual if necessary ---*/

      norm_r = r.norm();
      if (norm_r < tol * norm0) break;
      if (((monitoring) && (masterRank)) && ((i + 1) % monitorFreq == 0)) {
        SU2_OMP_MASTER
        WriteHistory(i + 1, norm_r / norm0);
        END_SU2_OMP_MASTER
      }
    }
  }

  /*--- Recalculate final residual (this should be optional) ---*/

  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("BCGSTAB", i, norm_r / norm0);
      END_SU2_OMP_MASTER
    }

    if (recomputeRes) {
      mat_vec(x, A_x);
      r = b - A_x;
      ScalarType true_res = r.norm();

      if ((fabs(true_res - norm_r) > tol * 10.0) && (masterRank)) {
        SU2_OMP_MASTER
        WriteWarning(norm_r, true_res, tol);
        END_SU2_OMP_MASTER
      }
    }
  }

  residual = norm_r / norm0;
  return i;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::Smoother_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                        const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                        const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                        unsigned long m, ScalarType& residual, bool monitoring,
                                                        const CConfig* config) const {
  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  const bool fix_iter_mode = tol < eps;
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Relaxation factor, see comments inside the loop over the smoothing iterations. ---*/
  const ScalarType omega = SU2_TYPE::GetValue(config->GetLinear_Solver_Smoother_Relaxation());

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate vectors for residual (r), solution increment (z), and matrix-vector
   product (A_x), for the latter two this is done only on the first call to the method. ---*/

  if (!smooth_ready) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      auto nVar = b.GetNVar();
      auto nBlk = b.GetNBlk();
      auto nBlkDomain = b.GetNBlkDomain();

      A_x.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      z.Initialize(nBlk, nBlkDomain, nVar, nullptr);

      smooth_ready = true;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Compute the initial residual and check if the system is already solved (if in COMM_FULL mode). ---*/

  if (!xIsZero) {
    mat_vec(x, A_x);
    r = b - A_x;
  } else {
    r = b;
  }

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
    norm_r = r.norm();
    norm0 = b.norm();

    /*--- Set the norm to the initial initial residual value ---*/

    if (tol_type == LinearToleranceType::RELATIVE) norm0 = norm_r;

    if ((norm_r < tol * norm0) || (norm_r < eps)) {
      if (masterRank) {
        SU2_OMP_MASTER
        cout << "CSysSolve::Smoother_LinSolver(): system solved by initial guess." << endl;
        END_SU2_OMP_MASTER
      }
      return 0;
    }

    /*--- Output header information including initial residual. ---*/

    if ((monitoring) && (masterRank)) {
      SU2_OMP_MASTER {
        WriteHeader("Smoother", tol, norm_r);
        WriteHistory(i, norm_r / norm0);
      }
      END_SU2_OMP_MASTER
    }
  }

  /*--- Smoothing Iterations ---*/

  for (i = 0; i < m; i++) {
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

    x += omega * z;
    r -= omega * A_x;

    /*--- Only compute the residuals in full communication mode. ---*/
    /*--- Check if solution has converged, else output the relative residual if necessary. ---*/

    if (!fix_iter_mode && config->GetComm_Level() == COMM_FULL) {
      norm_r = r.norm();
      if (norm_r < tol * norm0) break;
      if (((monitoring) && (masterRank)) && ((i + 1) % monitorFreq == 0)) {
        SU2_OMP_MASTER
        WriteHistory(i + 1, norm_r / norm0);
        END_SU2_OMP_MASTER
      }
    }
  }

  if (fix_iter_mode) norm_r = r.norm();

  if ((monitoring) && (masterRank) && (config->GetComm_Level() == COMM_FULL)) {
    SU2_OMP_MASTER
    WriteFinalResidual("Smoother", i, norm_r / norm0);
    END_SU2_OMP_MASTER
  }

  residual = norm_r / norm0;
  return i;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::Solve(CSysMatrix<ScalarType>& Jacobian, const CSysVector<su2double>& LinSysRes,
                                           CSysVector<su2double>& LinSysSol, CGeometry* geometry,
                                           const CConfig* config) {
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
  unsigned long MaxIter;
  ScalarType SolverTol;
  bool ScreenOutput;

  switch (lin_sol_mode) {
    /*--- Mesh Deformation mode ---*/
    case LINEAR_SOLVER_MODE::MESH_DEFORM: {
      KindSolver = config->GetKind_Deform_Linear_Solver();
      KindPrecond = config->GetKind_Deform_Linear_Solver_Prec();
      MaxIter = config->GetDeform_Linear_Solver_Iter();
      SolverTol = SU2_TYPE::GetValue(config->GetDeform_Linear_Solver_Error());
      ScreenOutput = config->GetDeform_Output();
      break;
    }

    /*--- Gradient Smoothing mode ---*/
    case LINEAR_SOLVER_MODE::GRADIENT_MODE: {
      KindSolver = config->GetKind_Grad_Linear_Solver();
      KindPrecond = config->GetKind_Grad_Linear_Solver_Prec();
      MaxIter = config->GetGrad_Linear_Solver_Iter();
      SolverTol = SU2_TYPE::GetValue(config->GetGrad_Linear_Solver_Error());
      ScreenOutput = true;
      break;
    }

    /*--- Normal mode assumes that 'lin_sol_mode==LINEAR_SOLVER_MODE::STANDARD',
     * but does not enforce it to avoid compiler warning.  ---*/
    default: {
      KindSolver = config->GetKind_Linear_Solver();
      KindPrecond = config->GetKind_Linear_Solver_Prec();
      MaxIter = config->GetLinear_Solver_Iter();
      SolverTol = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
      ScreenOutput = false;
      break;
    }
  }

  const bool nested = (KindSolver == FGMRES || KindSolver == RESTARTED_FGMRES || KindSolver == SMOOTHER) &&
                      config->GetKind_Linear_Solver_Inner() != LINEAR_SOLVER_INNER::NONE;

  if (nested && !inner_solver) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      inner_solver = std::make_unique<CSysSolve<ScalarType>>(LINEAR_SOLVER_MODE::STANDARD);
      inner_solver->SetxIsZero(true);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Stop the recording for the linear solver ---*/
  bool TapeActive = NO;

  if (config->GetDiscrete_Adjoint()) {
#ifdef CODI_REVERSE_TYPE

    TapeActive = AD::TapeActive();

    /*--- Declare external function inputs, outputs, and data ---*/
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      AD::SetExtFuncIn(&LinSysRes[0], LinSysRes.GetLocSize());
      AD::SetExtFuncOut(&LinSysSol[0], LinSysSol.GetLocSize());
      AD::FuncHelper.addUserData(&LinSysRes);
      AD::FuncHelper.addUserData(&LinSysSol);
      AD::FuncHelper.addUserData(&Jacobian);
      AD::FuncHelper.addUserData(geometry);
      AD::FuncHelper.addUserData(config);
      AD::FuncHelper.addUserData(this);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
#endif
  }

  unsigned long IterLinSol = 0;

  /*--- Declaration of the external function ---*/
  auto externalFunction = [&]() {
    /*--- Create matrix-vector product, preconditioner, and solve the linear system ---*/

    HandleTemporariesIn(LinSysRes, LinSysSol);

    auto mat_vec = CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);

    /*--- Build preconditioner. ---*/

    const auto kindPrec = static_cast<ENUM_LINEAR_SOLVER_PREC>(KindPrecond);
    auto* normal_prec = CPreconditioner<ScalarType>::Create(kindPrec, Jacobian, geometry, config);
    normal_prec->Build();

    CPreconditioner<ScalarType>* nested_prec = nullptr;
    if (nested) {
      auto f = [&](const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) {
        /*--- Initialize to 0 to be safe. ---*/
        v = ScalarType{};
        ScalarType unused{};
        /*--- Handle other types here if desired but do not call Solve because
         * that will create issues with the AD external function. ---*/
        (void)inner_solver->BCGSTAB_LinSolver(u, v, mat_vec, *normal_prec, sqrt(SolverTol), MaxIter, unused, false,
                                              config);
      };
      nested_prec = new CAbstractPreconditioner<ScalarType>(f);
    }
    const auto* precond = nested ? nested_prec : normal_prec;

    /*--- Solve system. ---*/

    ScalarType residual = 0.0;

    switch (KindSolver) {
      case BCGSTAB:
        IterLinSol = BCGSTAB_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                       ScreenOutput, config);
        break;
      case FGMRES:
        IterLinSol = FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                      ScreenOutput, config);
        break;
      case FGCRODR:
        IterLinSol = FGCRODR_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                       ScreenOutput, config);
        break;
      case RESTARTED_FGMRES:
        IterLinSol = RFGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                       ScreenOutput, config);
        break;
      case CONJUGATE_GRADIENT:
        IterLinSol = CG_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                  ScreenOutput, config);
        break;
      case SMOOTHER:
        IterLinSol = Smoother_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                        ScreenOutput, config);
        break;
      case PASTIX_LDLT:
      case PASTIX_LU:
        Jacobian.BuildPastixPreconditioner(geometry, config, KindSolver);
        Jacobian.ComputePastixPreconditioner(*LinSysRes_ptr, *LinSysSol_ptr, geometry, config);
        IterLinSol = 1;
        residual = 1e-20;
        break;
      default:
        SU2_MPI::Error("Unknown type of linear solver.", CURRENT_FUNCTION);
    }

    SU2_OMP_MASTER {
      Residual = residual;
      Iterations = IterLinSol;
    }
    END_SU2_OMP_MASTER

    HandleTemporariesOut(LinSysSol);

    delete normal_prec;
    delete nested_prec;

    if (TapeActive) {
      /*--- To keep the behavior of SU2_DOT, but not strictly required since jacobian is symmetric(?). ---*/
      const bool RequiresTranspose =
          ((lin_sol_mode != LINEAR_SOLVER_MODE::MESH_DEFORM) || (config->GetKind_SU2() == SU2_COMPONENT::SU2_DOT));

      if (lin_sol_mode == LINEAR_SOLVER_MODE::MESH_DEFORM)
        KindPrecond = config->GetKind_Deform_Linear_Solver_Prec();
      else if (lin_sol_mode == LINEAR_SOLVER_MODE::GRADIENT_MODE)
        KindPrecond = config->GetKind_Grad_Linear_Solver_Prec();
      else
        KindPrecond = config->GetKind_DiscAdj_Linear_Prec();

      /*--- Build preconditioner for the transposed Jacobian ---*/

      if (RequiresTranspose) Jacobian.TransposeInPlace();

      switch (KindPrecond) {
        case ILU:
          if (RequiresTranspose) Jacobian.BuildILUPreconditioner();
          break;
        case JACOBI:
        case LINELET:
          if (RequiresTranspose) Jacobian.BuildJacobiPreconditioner();
          break;
        case LU_SGS:
          /*--- Nothing to build. ---*/
          break;
        case PASTIX_ILU:
        case PASTIX_LU_P:
        case PASTIX_LDLT_P:
          /*--- It was already built. ---*/
          break;
        default:
          SU2_MPI::Error("The specified preconditioner is not yet implemented for the discrete adjoint method.",
                         CURRENT_FUNCTION);
          break;
      }
    }
  }; /*--- Finish declaration of the external function ---*/

#ifdef CODI_REVERSE_TYPE
  /*--- Call the external function with appropriate AD handling ---*/
  AD::FuncHelper.callPrimalFuncWithADType(externalFunction);

  AD::FuncHelper.addToTape(CSysSolve_b<ScalarType>::Solve_b);
#else
  /*--- Without reverse AD, call the external function directly ---*/
  externalFunction();
#endif

  return IterLinSol;
}

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::Solve_b(CSysMatrix<ScalarType>& Jacobian, const CSysVector<su2double>& LinSysRes,
                                             CSysVector<su2double>& LinSysSol, CGeometry* geometry,
                                             const CConfig* config, const bool directCall) {
  unsigned short KindSolver, KindPrecond;
  unsigned long MaxIter, IterLinSol = 0;
  ScalarType SolverTol;
  bool ScreenOutput;

  switch (lin_sol_mode) {
    /*--- Mesh Deformation mode ---*/
    case LINEAR_SOLVER_MODE::MESH_DEFORM: {
      KindSolver = config->GetKind_Deform_Linear_Solver();
      KindPrecond = config->GetKind_Deform_Linear_Solver_Prec();
      MaxIter = config->GetDeform_Linear_Solver_Iter();
      SolverTol = SU2_TYPE::GetValue(config->GetDeform_Linear_Solver_Error());
      ScreenOutput = config->GetDeform_Output();
      break;
    }

    /*--- Gradient Smoothing mode ---*/
    case LINEAR_SOLVER_MODE::GRADIENT_MODE: {
      KindSolver = config->GetKind_Grad_Linear_Solver();
      KindPrecond = config->GetKind_Grad_Linear_Solver_Prec();
      MaxIter = config->GetGrad_Linear_Solver_Iter();
      SolverTol = SU2_TYPE::GetValue(config->GetGrad_Linear_Solver_Error());
      ScreenOutput = true;
      break;
    }

    /*--- Normal mode assumes that 'lin_sol_mode==LINEAR_SOLVER_MODE::STANDARD',
     * but does not enforce it to avoid compiler warning. ---*/
    default: {
      KindSolver = config->GetKind_Linear_Solver();
      KindPrecond = config->GetKind_Linear_Solver_Prec();
      MaxIter = config->GetLinear_Solver_Iter();
      SolverTol = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
      ScreenOutput = false;
      break;
    }
  }

  const bool nested = (KindSolver == FGMRES || KindSolver == RESTARTED_FGMRES || KindSolver == SMOOTHER) &&
                      config->GetKind_Linear_Solver_Inner() != LINEAR_SOLVER_INNER::NONE;

  if (nested && !inner_solver) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      inner_solver = std::make_unique<CSysSolve<ScalarType>>(LINEAR_SOLVER_MODE::STANDARD);
      inner_solver->SetxIsZero(true);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Set up preconditioner and matrix-vector product ---*/

  auto mat_vec = CSysMatrixVectorProduct<ScalarType>(Jacobian, geometry, config);

  const auto kindPrec = static_cast<ENUM_LINEAR_SOLVER_PREC>(KindPrecond);
  auto* normal_prec = CPreconditioner<ScalarType>::Create(kindPrec, Jacobian, geometry, config);

  /*--- If there was no call to solve first the preconditioner needs to be built here. ---*/
  if (directCall) {
    Jacobian.TransposeInPlace();
    normal_prec->Build();
  }

  CPreconditioner<ScalarType>* nested_prec = nullptr;
  if (nested) {
    auto f = [&](const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) {
      /*--- Initialize to 0 to be safe. ---*/
      v = ScalarType{};
      ScalarType unused{};
      /*--- Handle other types here if desired but do not call Solve because
       * that will create issues with the AD external function. ---*/
      (void)inner_solver->BCGSTAB_LinSolver(u, v, mat_vec, *normal_prec, sqrt(SolverTol), MaxIter, unused, false,
                                            config);
    };
    nested_prec = new CAbstractPreconditioner<ScalarType>(f);
  }
  const auto* precond = nested ? nested_prec : normal_prec;

  /*--- Solve the system ---*/

  /*--- Local variable to prevent all threads from writing to a shared location (this->Residual). ---*/
  ScalarType residual = 0.0;

  HandleTemporariesIn(LinSysRes, LinSysSol);

  switch (KindSolver) {
    case FGMRES:
      IterLinSol = FGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                    ScreenOutput, config);
      break;
    case FGCRODR:
      IterLinSol = FGCRODR_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                     ScreenOutput, config);
      break;
    case RESTARTED_FGMRES:
      IterLinSol = RFGMRES_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                     ScreenOutput, config);
      break;
    case BCGSTAB:
      IterLinSol = BCGSTAB_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                     ScreenOutput, config);
      break;
    case CONJUGATE_GRADIENT:
      IterLinSol = CG_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                ScreenOutput, config);
      break;
    case SMOOTHER:
      IterLinSol = Smoother_LinSolver(*LinSysRes_ptr, *LinSysSol_ptr, mat_vec, *precond, SolverTol, MaxIter, residual,
                                      ScreenOutput, config);
      break;
    case PASTIX_LDLT:
    case PASTIX_LU:
      if (directCall) Jacobian.BuildPastixPreconditioner(geometry, config, KindSolver);
      Jacobian.ComputePastixPreconditioner(*LinSysRes_ptr, *LinSysSol_ptr, geometry, config);
      IterLinSol = 1;
      residual = 1e-20;
      break;
    default:
      SU2_MPI::Error("Unknown type of linear solver.", CURRENT_FUNCTION);
      break;
  }

  HandleTemporariesOut(LinSysSol);

  delete normal_prec;
  delete nested_prec;

  SU2_OMP_MASTER {
    Residual = residual;
    Iterations = IterLinSol;
  }
  END_SU2_OMP_MASTER

  return IterLinSol;
}

/*--- Explicit instantiations ---*/

#ifdef CODI_FORWARD_TYPE
template class CSysSolve<su2double>;
#else
template class CSysSolve<su2mixedfloat>;
#ifdef USE_MIXED_PRECISION
template class CSysSolve<passivedouble>;
#endif
#endif
