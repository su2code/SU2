/*!
 * \file CSysSolveGPU.cu
 * \brief CUDA/GPU skeleton implementations for Krylov solvers.
 * \author Jesse Li
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../include/linear_algebra/CPreconditioner.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"

#include <cmath>
#include <iostream>

void SU2_GPU_BeginSolverBLASContext();
void SU2_GPU_EndSolverBLASContext();

namespace {

class CGPUSolverBLASContextGuard {
 public:
  CGPUSolverBLASContextGuard() { SU2_GPU_BeginSolverBLASContext(); }

  ~CGPUSolverBLASContextGuard() { SU2_GPU_EndSolverBLASContext(); }

  CGPUSolverBLASContextGuard(const CGPUSolverBLASContextGuard&) = delete;
  CGPUSolverBLASContextGuard& operator=(const CGPUSolverBLASContextGuard&) = delete;
};

}  // namespace

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolver_GPU(const VectorType& b, VectorType& x,
                                                          const ProductType& mat_vec, const PrecondType& precond,
                                                          ScalarType tol, unsigned long m, ScalarType& residual,
                                                          bool monitoring, const CConfig* config) const {
  SU2_ZONE_SCOPED

  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  const bool flexible = !precond.IsIdentity();

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  if (m > 1000) {
    SU2_MPI::Error("FGMRES subspace is too large (>1000).", CURRENT_FUNCTION);
  }

  CGPUSolverBLASContextGuard blas_context;
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

  su2vector<ScalarType> g(m + 1), sn(m + 1), cs(m + 1), y(m);
  g = ScalarType(0);
  sn = ScalarType(0);
  cs = ScalarType(0);
  y = ScalarType(0);

  su2matrix<ScalarType> H(m + 1, m);
  H = ScalarType(0);

  /*--- b and x enter from the host side. The solver owns all subsequent
   * synchronization for vectors it updates with GPU primitives. ---*/
  b.HtDTransfer();
  if (xIsZero) {
    x.GPUSetVal(ScalarType(0));
  } else {
    x.HtDTransfer();
  }

  ScalarType norm0 = b.GPUNorm();

  if (!xIsZero) {
    mat_vec(x, V[0]);
    V[0].GPUAxpy(ScalarType(-1), b);
  } else {
    V[0].GPUCopy(b);
    V[0].GPUScale(ScalarType(-1));
  }

  ScalarType beta = xIsZero ? norm0 : V[0].GPUNorm();

  SU2_OMP_MASTER
  ResetDeflation();
  END_SU2_OMP_MASTER

  if (tol_type == LinearToleranceType::RELATIVE) norm0 = beta;

  if (beta < tol * norm0 || beta < eps) {
    if (masterRank) {
      SU2_OMP_MASTER
      cout << "CSysSolve::FGMRES_GPU(): system solved by initial guess." << endl;
      END_SU2_OMP_MASTER
    }
    residual = beta;
    x.DtHTransfer();
    return 0;
  }

  V[0].GPUScale(ScalarType(-1) / beta);
  g[0] = beta;

  unsigned long i = 0;
  if (monitoring && masterRank) {
    SU2_OMP_MASTER {
      WriteHeader("FGMRES_GPU", tol, beta);
      WriteHistory(i, beta / norm0);
    }
    END_SU2_OMP_MASTER
  }

  for (i = 0; i < m; i++) {
    if (beta < tol * norm0) break;

    /*--- mat_vec consumes and produces device-resident vectors on the CUDA path. ---*/
    if (flexible) {
      precond(V[i], Z[i]);
      mat_vec(Z[i], V[i + 1]);
    } else {
      mat_vec(V[i], V[i + 1]);
    }

    /*--- Classical Gram-Schmidt twice, matching the CPU implementation's
     * numerical structure while using GPU dot/axpy/norm primitives. ---*/
    for (unsigned long k = 0; k <= i; k++) {
      H(k, i) = V[k].GPUDot(V[i + 1]);
      V[i + 1].GPUAxpy(-H(k, i), V[k]);
    }

    for (unsigned long k = 0; k <= i; k++) {
      const ScalarType dh = V[k].GPUDot(V[i + 1]);
      H(k, i) += dh;
      V[i + 1].GPUAxpy(-dh, V[k]);
    }

    const ScalarType nrm = V[i + 1].GPUNorm();
    if (nrm <= ScalarType(0) || nrm != nrm) {
      H(i + 1, i) = ScalarType(0);
      if (masterRank) {
        SU2_OMP_MASTER
        cout << "WARNING: FGMRES GPU orthogonalization failed, linear solver diverged." << endl;
        END_SU2_OMP_MASTER
      }
      break;
    }

    H(i + 1, i) = nrm;
    V[i + 1].GPUScale(ScalarType(1) / nrm);

    for (unsigned long k = 0; k < i; k++) ApplyGivens(sn[k], cs[k], H(k, i), H(k + 1, i));
    GenerateGivens(H(i, i), H(i + 1, i), sn[i], cs[i]);
    ApplyGivens(sn[i], cs[i], g[i], g[i + 1]);

    beta = fabs(g[i + 1]);

    if (monitoring && masterRank && ((i + 1) % monitorFreq == 0)) {
      SU2_OMP_MASTER
      WriteHistory(i + 1, beta / norm0);
      END_SU2_OMP_MASTER
    }
  }

  SolveReduced(i, H, g, y);

  const auto& basis = flexible ? Z : V;

  for (unsigned long k = 0; k < i; k++) {
    x.GPUAxpy(y[k], basis[k]);
  }

  x.DtHTransfer();

  if (monitoring && config->GetComm_Level() == COMM_FULL) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("FGMRES_GPU", i, beta / norm0);
      END_SU2_OMP_MASTER
    }

    if (recomputeRes) {
      mat_vec(x, V[0]);
      V[0].GPUAxpy(ScalarType(-1), b);
      const ScalarType res = V[0].GPUNorm();

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

/*--- Explicit instantiations for the GPU solver helper.
 * Keep these aligned with the scalar types instantiated for CSysSolve. ---*/
template unsigned long CSysSolve<su2mixedfloat>::FGMRES_LinSolver_GPU(
    const CSysVector<su2mixedfloat>& b, CSysVector<su2mixedfloat>& x,
    const CMatrixVectorProduct<su2mixedfloat>& mat_vec, const CPreconditioner<su2mixedfloat>& precond,
    su2mixedfloat tol, unsigned long m, su2mixedfloat& residual, bool monitoring, const CConfig* config) const;

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template unsigned long CSysSolve<passivedouble>::FGMRES_LinSolver_GPU(
    const CSysVector<passivedouble>& b, CSysVector<passivedouble>& x,
    const CMatrixVectorProduct<passivedouble>& mat_vec, const CPreconditioner<passivedouble>& precond,
    passivedouble tol, unsigned long m, passivedouble& residual, bool monitoring, const CConfig* config) const;
#endif
