/*!
 * \file CSysSolveFGMRES.inl
 * \brief Shared implementation of the Flexible GMRES linear solver.
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

namespace fgmres_impl {

template <class ScalarType, class FGMRESOps>
bool Orthogonalize(unsigned long i, su2matrix<ScalarType>& H, std::vector<CSysVector<ScalarType>>& V,
                   const FGMRESOps& ops) {
  for (unsigned long k = 0; k <= i; k++) {
    H(k, i) = ops.Dot(V[k], V[i + 1]);
    ops.AddScaled(V[i + 1], -H(k, i), V[k]);
  }

  for (unsigned long k = 0; k <= i; k++) {
    const ScalarType dh = ops.Dot(V[k], V[i + 1]);
    H(k, i) += dh;
    ops.AddScaled(V[i + 1], -dh, V[k]);
  }

  const ScalarType nrm = ops.Norm(V[i + 1]);
  if (nrm <= ScalarType(0) || nrm != nrm) {
    H(i + 1, i) = ScalarType(0);
    return false;
  }

  H(i + 1, i) = nrm;
  ops.Divide(V[i + 1], nrm);
  return true;
}

}  // namespace fgmres_impl

template <class ScalarType>
template <class FGMRESOps>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolverImpl(const VectorType& b, VectorType& x,
                                                          const ProductType& mat_vec, const PrecondType& precond,
                                                          ScalarType tol, unsigned long m, ScalarType& residual,
                                                          bool monitoring, const CConfig* config,
                                                          const FGMRESOps& ops) const {
  SU2_ZONE_SCOPED

  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  const bool flexible = !precond.IsIdentity();
  const bool nestedParallel = ops.NestedParallel();

  /*--- Check the subspace size. ---*/

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  if (m > 1000) {
    SU2_MPI::Error("FGMRES subspace is too large (>1000).", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet. ---*/

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

  ops.PrepareInput(b, x, xIsZero);

  /*--- Calculate the norm of the rhs vector. ---*/

  ScalarType norm0 = ops.Norm(b);

  /*--- Calculate the initial residual (actually the negative residual) and compute its norm. ---*/

  if (!xIsZero) {
    mat_vec(x, V[0]);
    ops.Subtract(V[0], b);
  } else {
    ops.AssignNegative(V[0], b);
  }

  ScalarType beta = xIsZero ? norm0 : ops.Norm(V[0]);

  /*--- FGMRES shares V and Z with FGCRODR, reset deflation because we're breaking
   * relations between subspaces that FGCRODR relies on for recycling. ---*/

  SU2_OMP_MASTER
  ResetDeflation();
  END_SU2_OMP_MASTER

  /*--- Set the norm to the initial initial residual value. ---*/

  if (tol_type == LinearToleranceType::RELATIVE) norm0 = beta;

  if (beta < tol * norm0 || beta < eps) {
    /*--- System is already solved. ---*/

    if (masterRank) {
      SU2_OMP_MASTER
      cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
      END_SU2_OMP_MASTER
    }
    residual = beta;
    ops.FinalizeOutput(x);
    return 0;
  }

  /*--- Normalize residual to get w_{0} (the negative sign is because w[0]
        holds the negative residual, as mentioned above). ---*/

  ops.Divide(V[0], -beta);

  /*--- Initialize the RHS of the reduced system. ---*/

  g[0] = beta;

  /*--- Output header information including initial residual. ---*/

  unsigned long i = 0;
  if (monitoring && masterRank) {
    SU2_OMP_MASTER {
      WriteHeader("FGMRES", tol, beta);
      WriteHistory(i, beta / norm0);
    }
    END_SU2_OMP_MASTER
  }

  /*--- Loop over all search directions. ---*/

  for (i = 0; i < m; i++) {
    /*--- Check if solution has converged. ---*/

    if (beta < tol * norm0) break;

    if (flexible) {
      /*--- Precondition the CSysVector v[i] and store result in z[i]. ---*/

      precond(V[i], Z[i]);

      /*--- Add to Krylov subspace. ---*/

      mat_vec(Z[i], V[i + 1]);
    } else {
      mat_vec(V[i], V[i + 1]);
    }

    /*--- Modified Gram-Schmidt orthogonalization. ---*/

    bool orthog_ok;
    if constexpr (FGMRESOps::UseSolverModGramSchmidt) {
      if (nestedParallel) {
        /*--- "omp parallel if" does not work well here. ---*/
        SU2_OMP_PARALLEL
        orthog_ok = ModGramSchmidt(true, i, H, V);
        END_SU2_OMP_PARALLEL
      } else {
        orthog_ok = ModGramSchmidt(false, i, H, V);
      }
    } else {
      orthog_ok = fgmres_impl::Orthogonalize(i, H, V, ops);
    }

    if (!orthog_ok) {
      if (masterRank) {
        SU2_OMP_MASTER
        cout << "WARNING: FGMRES orthogonalization failed, linear solver diverged." << endl;
        END_SU2_OMP_MASTER
      }
      break;
    }

    /*--- Apply old Givens rotations to new column of the Hessenberg matrix then generate the
     new Givens rotation matrix and apply it to the last two elements of H[:][i] and g. ---*/

    for (unsigned long k = 0; k < i; k++) ApplyGivens(sn[k], cs[k], H(k, i), H(k + 1, i));
    GenerateGivens(H(i, i), H(i + 1, i), sn[i], cs[i]);
    ApplyGivens(sn[i], cs[i], g[i], g[i + 1]);

    /*--- Set L2 norm of residual and check if solution has converged. ---*/

    beta = fabs(g[i + 1]);

    /*--- Output the relative residual if necessary. ---*/

    if (monitoring && masterRank && ((i + 1) % monitorFreq == 0)) {
      SU2_OMP_MASTER
      WriteHistory(i + 1, beta / norm0);
      END_SU2_OMP_MASTER
    }
  }

  /*--- Solve the least-squares system and update solution. ---*/

  SolveReduced(i, H, g, y);

  const auto& basis = flexible ? Z : V;

  ops.UpdateSolution(i, basis, y, x, nestedParallel);
  ops.FinalizeOutput(x);

  /*--- Recalculate final (neg.) residual (this should be optional). ---*/

  if (monitoring && config->GetComm_Level() == COMM_FULL) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("FGMRES", i, beta / norm0);
      END_SU2_OMP_MASTER
    }

    if (recomputeRes) {
      mat_vec(x, V[0]);
      ops.Subtract(V[0], b);
      ScalarType res = ops.Norm(V[0]);

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
