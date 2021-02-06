/*!
 * \file CNewtonIntegration.cpp
 * \brief Newton-Krylov integration.
 * \author P. Gomes
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/integration/CNewtonIntegration.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"

using Scalar = CNewtonIntegration::Scalar;

namespace {

class CMatrixFreeProductWrapper final : public CMatrixVectorProduct<Scalar> {
  CNewtonIntegration* integration;
public:
  CMatrixFreeProductWrapper(CNewtonIntegration* i) : integration(i) {}

  /*!
   * \brief Operator for the product operation.
   */
  inline void operator()(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const override {
    integration->MatrixFreeProduct(u, v);
  }
};

class CPreconditionerWrapper final : public CPreconditioner<Scalar> {
  const CNewtonIntegration* integration;
public:
  CPreconditionerWrapper(const CNewtonIntegration* i) : integration(i) {}

  /*!
   * \brief Operator for the preconditioning operation.
   */
  inline void operator()(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const override {
    integration->Preconditioner(u, v);
  }
};
}

CNewtonIntegration::~CNewtonIntegration() { delete preconditioner; }

void CNewtonIntegration::Setup() {

  const auto nVar = solvers[FLOW_SOL]->GetnVar();
  const auto nPoint = geometry->GetnPoint();
  const auto nPointDomain = geometry->GetnPointDomain();

  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), 1024);

  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  if (!std::is_same<Scalar,su2double>::value) {
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, nullptr);
  }

  /*--- Check if the solver is able to provide a linear preconditioner. ---*/
  if (config->GetKind_TimeIntScheme() != EULER_IMPLICIT) return;

  switch (config->GetKind_Linear_Solver_Prec()) {
    case JACOBI:
      preconditioner = new CJacobiPreconditioner<MixedScalar>(solvers[FLOW_SOL]->Jacobian, geometry, config, false);
      break;
    case LINELET:
      preconditioner = new CLineletPreconditioner<MixedScalar>(solvers[FLOW_SOL]->Jacobian, geometry, config);
      break;
    case LU_SGS:
      preconditioner = new CLU_SGSPreconditioner<MixedScalar>(solvers[FLOW_SOL]->Jacobian, geometry, config);
      break;
    case ILU:
      preconditioner = new CILUPreconditioner<MixedScalar>(solvers[FLOW_SOL]->Jacobian, geometry, config, false);
      break;
    case PASTIX_ILU: case PASTIX_LU_P: case PASTIX_LDLT_P:
      preconditioner = new CPastixPreconditioner<MixedScalar>(solvers[FLOW_SOL]->Jacobian, geometry, config,
                                                              config->GetKind_Linear_Solver_Prec(), false);
      break;
  }

  if (!std::is_same<Scalar,MixedScalar>::value) {
    precondIn.Initialize(nPoint, nPointDomain, nVar, nullptr);
    precondOut.Initialize(nPoint, nPointDomain, nVar, nullptr);
  }
}

void CNewtonIntegration::PerturbSolution(const CSysVector<Scalar>& dir, Scalar mag) {

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < solvers[FLOW_SOL]->GetnVar(); ++iVar)
      solvers[FLOW_SOL]->GetNodes()->AddSolution(iPoint,iVar, mag*dir(iPoint,iVar));
  }
}

void CNewtonIntegration::ComputeResiduals(ResEvalType type) {

  /*--- Save the default integration scheme, and force to explicit if required. ---*/
  auto TimeIntScheme = config->GetKind_TimeIntScheme();
  if (type == EXPLICIT) {
    SU2_OMP_MASTER
    config->SetKind_TimeIntScheme(EULER_EXPLICIT);
    SU2_OMP_BARRIER
  }

  solvers[FLOW_SOL]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  Space_Integration(geometry, solvers, numerics[FLOW_SOL], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS);

  /*--- Restore default. ---*/
  if (type == EXPLICIT) {
    SU2_OMP_MASTER
    config->SetKind_TimeIntScheme(TimeIntScheme);
    SU2_OMP_BARRIER
  }

}

void CNewtonIntegration::MultiGrid_Iteration(CGeometry ****geometry_, CSolver *****solvers_, CNumerics ******numerics_,
                                             CConfig **config_, unsigned short EqSystem, unsigned short iZone,
                                             unsigned short iInst) {
  config = config_[iZone];
  solvers = solvers_[iZone][iInst][MESH_0];
  geometry = geometry_[iZone][iInst][MESH_0];
  numerics = numerics_[iZone][iInst][MESH_0];

  if (!setup) { Setup(); setup = true; }

  /*--- The step for finite-difference-based matrix-free product depends on the RMS of the solution. ---*/
  su2double rmsSol = 0.0;

  SU2_OMP_PARALLEL_(if(solvers[FLOW_SOL]->GetHasHybridParallel())) {

  /*--- Compute the current residual and the approximate Jacobian for preconditioning. ---*/

  ComputeResiduals(DEFAULT);

  solvers[FLOW_SOL]->SetTime_Step(geometry, solvers, config, MESH_0, config->GetTimeIter());

  solvers[FLOW_SOL]->PrepareImplicitIteration(geometry, solvers, config);

  if (preconditioner) preconditioner->Build();

  /*--- Save current residuals and the solution to be able to perturb it. ---*/

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto i = 0ul; i < LinSysRes.GetNElmDomain(); ++i)
    LinSysRes[i] = SU2_TYPE::GetValue(solvers[FLOW_SOL]->LinSysRes[i]);

  solvers[FLOW_SOL]->Set_OldSolution();

  /*--- Compute RMS(solution). ---*/

  su2double rmsSol_loc = 0.0;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint)
    for (auto iVar = 0ul; iVar < solvers[FLOW_SOL]->GetnVar(); ++iVar)
      rmsSol_loc += pow(solvers[FLOW_SOL]->GetNodes()->GetSolution(iPoint,iVar), 2);

  atomicAdd(rmsSol_loc, rmsSol);

  SU2_OMP_BARRIER
  SU2_OMP_MASTER {
    su2double t = rmsSol;
    SU2_MPI::Allreduce(&t, &rmsSol, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    /// TODO: Customize the step size (1e-4).
    finDiffStep = 1e-4 * max(1.0, sqrt(SU2_TYPE::GetValue(rmsSol) / geometry->GetGlobal_nPointDomain()));
  }
  SU2_OMP_BARRIER

  /*--- Solve for the solution update. ---*/

  CMatrixFreeProductWrapper product(this);
  CPreconditionerWrapper precond(this);
  auto& linSysSol = GetSolutionVec(solvers[FLOW_SOL]->LinSysSol);

  Scalar eps = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
  auto iter = config->GetLinear_Solver_Iter();

  iter = LinSolver.FGMRES_LinSolver(LinSysRes, linSysSol, product, precond, eps, iter, eps, false, config, true);

  SetSolutionResult(solvers[FLOW_SOL]->LinSysSol);

  SU2_OMP_MASTER {
    solvers[FLOW_SOL]->SetIterLinSolver(iter);
    solvers[FLOW_SOL]->SetResLinSolver(eps);
  }
  SU2_OMP_BARRIER

  /// TODO: Clever back-tracking and CFL adaptation based on residual reduction.

  /*--- Update solution. ---*/

  solvers[FLOW_SOL]->CompleteImplicitIteration(geometry, solvers, config);

  /*--- Call the various post processings. ---*/

  solvers[FLOW_SOL]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  solvers[FLOW_SOL]->Postprocessing(geometry, solvers, config, MESH_0);

  SU2_OMP_MASTER {
    solvers[FLOW_SOL]->Pressure_Forces(geometry, config);
    solvers[FLOW_SOL]->Momentum_Forces(geometry, config);
    solvers[FLOW_SOL]->Friction_Forces(geometry, config);
  }

  } // end SU2_OMP_PARALLEL
}

void CNewtonIntegration::MatrixFreeProduct(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) {

  Scalar factor = finDiffStep / u.norm();

  PerturbSolution(u, factor);

  ComputeResiduals(EXPLICIT);

  /*--- Finalize product. ---*/
  factor = 1.0 / factor;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {
    su2double delta = (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint)) /
                      max(EPS, solvers[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < LinSysRes.GetNVar(); ++iVar) {
      Scalar perturbRes = SU2_TYPE::GetValue(solvers[FLOW_SOL]->LinSysRes(iPoint,iVar));

      /*--- The global residual had its sign flipped, so we add to get the difference. ---*/
      v(iPoint,iVar) = (perturbRes + LinSysRes(iPoint,iVar)) * factor;

      /*--- Pseudotime term of the true Jacobian. ---*/
      v(iPoint,iVar) += SU2_TYPE::GetValue(delta) * u(iPoint,iVar);
    }
  }

  solvers[FLOW_SOL]->Jacobian.InitiateComms(v, geometry, config, SOLUTION_MATRIX);
  solvers[FLOW_SOL]->Jacobian.CompleteComms(v, geometry, config, SOLUTION_MATRIX);
}

void CNewtonIntegration::Preconditioner(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const {

  if (preconditioner) {
    Preconditioner_impl(u, v);
  }
  else {
    /*--- Approximate diagonal preconditioner. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {
      su2double delta = solvers[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) /
                        (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint));
      SU2_OMP_SIMD
      for (auto iVar = 0ul; iVar < u.GetNVar(); ++iVar)
        v(iPoint,iVar) = SU2_TYPE::GetValue(delta) * u(iPoint,iVar);
    }

    solvers[FLOW_SOL]->Jacobian.InitiateComms(v, geometry, config, SOLUTION_MATRIX);
    solvers[FLOW_SOL]->Jacobian.CompleteComms(v, geometry, config, SOLUTION_MATRIX);
  }
}
