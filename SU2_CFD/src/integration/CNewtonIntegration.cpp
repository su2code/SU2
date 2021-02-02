/*!
 * \file CNewtonIntegration.cpp
 * \brief Coupled Newton integration.
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
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/linear_algebra/CPreconditioner.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../../Common/include/linear_algebra/CSysSolve.hpp"

#define PARALLEL_FOR SU2_OMP(for schedule(static,omp_chunk_size) nowait)

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

class CBlockJacobiPrecondWrapper final : public CPreconditioner<Scalar> {
  const CNewtonIntegration* integration;
public:
  CBlockJacobiPrecondWrapper(const CNewtonIntegration* i) : integration(i) {}

  /*!
   * \brief Operator for the preconditioning operation.
   */
  inline void operator()(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const override {
    integration->BlockJacobiPrecond(u, v);
  }
};
}

CNewtonIntegration::CNewtonIntegration() {
  KindSol2EqSys[FLOW_SOL] = RUNTIME_FLOW_SYS;
  KindSol2EqSys[TURB_SOL] = RUNTIME_TURB_SYS;
/// TODO: These solvers need Prepare/CompleteImplicitIteration methods.
//  KindSol2EqSys[HEAT_SOL] = RUNTIME_HEAT_SYS;
//  KindSol2EqSys[RAD_SOL] = RUNTIME_RADIATION_SYS;
}

CNewtonIntegration::~CNewtonIntegration() {
  for (auto p : preconditioners) delete p;
}

void CNewtonIntegration::Setup() {

  if (!kindSol.empty()) return;

  KindFlowSol = EULER;
  if (config->GetViscous()) KindFlowSol = NAVIER_STOKES;
  if (solvers[TURB_SOL]) KindFlowSol = RANS;

  const auto nPoint = geometry->GetnPoint();
  const auto nPointDomain = geometry->GetnPointDomain();
  unsigned long nVarTot = 0;

  thread_safe = true;
  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), 512);

  for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
    if (solvers[iSol] && iSol != MESH_SOL && iSol != ADJMESH_SOL) {

      if (KindSol2EqSys[iSol] == 0)
        SU2_MPI::Error("Some of the solvers do not support the coupled Newton method.", CURRENT_FUNCTION);

      auto nVar = solvers[iSol]->GetnVar();
      kindSol.push_back(iSol);
      nVars.push_back(nVar);
      nVarTot += nVar;

      thread_safe &= solvers[iSol]->GetHasHybridParallel();

      preconditioners.push_back(nullptr);
      precondIn.push_back(CSysVector<MixedScalar>());
      precondOut.push_back(CSysVector<MixedScalar>());

      /*--- Check if the solver is able to provide a linear preconditioner. ---*/
      if (config->GetKind_TimeIntScheme() != EULER_IMPLICIT) continue;

      auto& p = preconditioners.back();

      switch (config->GetKind_Linear_Solver_Prec()) {
        case JACOBI:
          p = new CJacobiPreconditioner<MixedScalar>(solvers[iSol]->Jacobian, geometry, config, false);
          break;
        case LINELET:
          p = new CLineletPreconditioner<MixedScalar>(solvers[iSol]->Jacobian, geometry, config);
          break;
        case LU_SGS:
          p = new CLU_SGSPreconditioner<MixedScalar>(solvers[iSol]->Jacobian, geometry, config);
          break;
        case ILU:
          p = new CILUPreconditioner<MixedScalar>(solvers[iSol]->Jacobian, geometry, config, false);
          break;
        case PASTIX_ILU: case PASTIX_LU_P: case PASTIX_LDLT_P:
          p = new CPastixPreconditioner<MixedScalar>(solvers[iSol]->Jacobian, geometry, config,
                                                     config->GetKind_Linear_Solver_Prec(), false);
          break;
      }

      if (!std::is_same<Scalar,MixedScalar>::value) {
        precondIn.back().Initialize(nPoint, nPointDomain, nVar, nullptr);
        precondOut.back().Initialize(nPoint, nPointDomain, nVar, nullptr);
      }
    }
  }

  LinSysRes.Initialize(nPoint, nPointDomain, nVarTot, nullptr);
  LinSysSol.Initialize(nPoint, nPointDomain, nVarTot, nullptr);

}

void CNewtonIntegration::ComputeResiduals(ResEvalType type) {

  /*--- Running all the pre and post processings first captures e.g. the
   * dependency of the flow residuals on the eddy viscosity. Sounds good
   * but the linear solver does not like it. If not done that way we get
   * a more triangular Jacobian, i.e. better conditioned. ---*/

  constexpr bool reallyCoupled = false;

  for (int step=0; step<1+reallyCoupled; ++step) {
    for (auto pos : kindSol) {
      /*--- Set global config parameters for the solver. ---*/
      const auto eqSys = KindSol2EqSys[pos];
      SU2_OMP_MASTER
      config->SetGlobalParam(KindFlowSol, eqSys);
      SU2_OMP_BARRIER

      /*--- Save the default integration scheme, and force to explicit if required. ---*/
      auto TimeIntScheme = config->GetKind_TimeIntScheme();
      if (type == EXPLICIT) {
        SU2_OMP_MASTER
        config->SetKind_TimeIntScheme(EULER_EXPLICIT);
        SU2_OMP_BARRIER
      }

      if (step==0) {
        solvers[pos]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, eqSys, false);
        if (reallyCoupled) solvers[pos]->Postprocessing(geometry, solvers, config, MESH_0);
      }
      if (step>0 || !reallyCoupled) {
        Space_Integration(geometry, solvers, numerics[pos], config, MESH_0, NO_RK_ITER, eqSys);
      }

      /*--- Restore default. ---*/
      if (type == EXPLICIT) {
        SU2_OMP_MASTER
        config->SetKind_TimeIntScheme(TimeIntScheme);
        SU2_OMP_BARRIER
      }
    }
  }
}

void CNewtonIntegration::MultiGrid_Iteration(CGeometry ****geometry_, CSolver *****solvers_, CNumerics ******numerics_,
                                             CConfig **config_, unsigned short EqSystem, unsigned short iZone,
                                             unsigned short iInst) {
  config = config_[iZone];
  solvers = solvers_[iZone][iInst][MESH_0];
  geometry = geometry_[iZone][iInst][MESH_0];
  numerics = numerics_[iZone][iInst][MESH_0];

  Setup();

  /*--- The step for finite-difference-based matrix-free product depends on the RMS of the solution. ---*/
  su2double rmsSol = 0.0;

  SU2_OMP_PARALLEL_(if(thread_safe)) {

  /*--- Compute the current residual and the approximate Jacobians for preconditioning. ---*/

  ComputeResiduals(DEFAULT);

  for (unsigned long iSol = 0, offset = 0; iSol < kindSol.size(); ++iSol) {
    const auto pos = kindSol[iSol];
    const auto nVar = nVars[iSol];

    const auto eqSys = KindSol2EqSys[pos];
    SU2_OMP_MASTER
    config->SetGlobalParam(KindFlowSol, eqSys);
    SU2_OMP_BARRIER

    solvers[pos]->SetTime_Step(geometry, solvers, config, MESH_0, config->GetTimeIter());

    solvers[pos]->PrepareImplicitIteration(geometry, solvers, config);

    /*--- Save current solution to be able to perturb it. ---*/
    solvers[pos]->GetNodes()->Set_OldSolution();

    /*--- Aggregate residuals. ---*/
    su2double rmsSol_loc = 0.0;

    PARALLEL_FOR
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {
      for (auto iVar = 0ul; iVar < nVar; ++iVar) {
        rmsSol_loc += pow(solvers[pos]->GetNodes()->GetSolution(iPoint,iVar), 2);
        LinSysRes(iPoint, offset+iVar) = SU2_TYPE::GetValue(solvers[pos]->LinSysRes(iPoint, iVar));
      }
    }
    atomicAdd(rmsSol_loc, rmsSol);

    offset += nVar;

    /*--- Build preconditioner. ---*/
    if (preconditioners[iSol]) preconditioners[iSol]->Build();

  }
  SU2_OMP_BARRIER

  SU2_OMP_MASTER {
    SU2_MPI::Allreduce(&rmsSol, &finDiffStep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /// TODO: Customize the step size (1e-4).
    /// TODO: Can we have one step size per variable? Probably not...
    finDiffStep = 1e-4 * max(1.0, sqrt(finDiffStep / geometry->GetGlobal_nPointDomain()));
  }
  SU2_OMP_BARRIER

  /*--- Solve for the search direction. ---*/

  CMatrixFreeProductWrapper product(this);
  CBlockJacobiPrecondWrapper precond(this);

  LinSysSol = Scalar(0.0);

  auto iter = config->GetLinear_Solver_Iter();
  Scalar tol = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
  Scalar eps = 0.0;
  auto nIter = LinSolver.FGMRES_LinSolver(LinSysRes, LinSysSol, product, precond,
                                          tol, iter, eps, false, config, true);
  SU2_OMP_MASTER
  for (auto pos : kindSol) {
    solvers[pos]->SetIterLinSolver(nIter);
    solvers[pos]->SetResLinSolver(eps);
  }
  SU2_OMP_BARRIER

  /*--- Update solution. ---*/
  /* For now we let each solver handle the search direction, using its own form of under-relaxation
   * to then adapt the CFL. However we should also check global descent, and use a back-tracking
   * strategy that informs the CFL adaptation to increase/decrease/freeze. */

  for (unsigned long iSol = 0, offset = 0; iSol < kindSol.size(); ++iSol) {
    const auto pos = kindSol[iSol];
    const auto nVar = nVars[iSol];

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint)
      for (auto iVar = 0ul; iVar < nVar; ++iVar)
        solvers[pos]->LinSysSol(iPoint,iVar) = LinSysSol(iPoint,offset+iVar);

    solvers[pos]->CompleteImplicitIteration(geometry, solvers, config);

    /*--- Call the various post processings to replicate what happens in Single and MG iterations,
     * it does not seem be important to run the flow solver preprocessing in output-mode. ---*/

    solvers[pos]->Postprocessing(geometry, solvers, config, MESH_0);

    SU2_OMP_MASTER
    switch (KindSol2EqSys[pos]) {
      case RUNTIME_FLOW_SYS:
        solvers[pos]->Pressure_Forces(geometry, config);
        solvers[pos]->Momentum_Forces(geometry, config);
        solvers[pos]->Friction_Forces(geometry, config);
        break;
      case RUNTIME_HEAT_SYS:
        solvers[pos]->Heat_Fluxes(geometry, solvers, config);
        break;
    }
    SU2_OMP_BARRIER

    offset += nVar;
  }

  } // end SU2_OMP_PARALLEL
}

void CNewtonIntegration::MatrixFreeProduct(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) {

  /*--- Perturb the solution. ---*/
  Scalar factor = finDiffStep / u.norm();

  for (unsigned long iSol = 0, offset = 0; iSol < kindSol.size(); ++iSol) {
    const auto pos = kindSol[iSol];
    const auto nVar = nVars[iSol];

    PARALLEL_FOR
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint)
      for (auto iVar = 0ul; iVar < nVar; ++iVar)
        solvers[pos]->GetNodes()->Add_DeltaSolution(iPoint,iVar, u(iPoint,offset+iVar)*factor);

    offset += nVar;
  }
  SU2_OMP_BARRIER


  /*--- Compute residuals after perturbation. ---*/

  ComputeResiduals(EXPLICIT);


  /*--- Finalize product and restore the old solution. ---*/
  factor = 1.0 / factor;

  for (unsigned long iSol = 0, offset = 0; iSol < kindSol.size(); ++iSol) {
    const auto pos = kindSol[iSol];
    const auto nVar = nVars[iSol];

    solvers[pos]->GetNodes()->Set_Solution();

    PARALLEL_FOR
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {
      su2double delta = (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint)) /
                        max(EPS, solvers[pos]->GetNodes()->GetDelta_Time(iPoint));

      for (auto iVar = 0ul; iVar < nVar; ++iVar) {
        Scalar perturbRes = SU2_TYPE::GetValue(solvers[pos]->LinSysRes(iPoint,iVar));

        /*--- The global residual had its sign flipped, so we add to get the difference. ---*/
        v(iPoint,offset+iVar) = (perturbRes + LinSysRes(iPoint,offset+iVar)) * factor;

        /*--- Pseudotime term of the true Jacobian. ---*/
        v(iPoint,offset+iVar) += SU2_TYPE::GetValue(delta) * u(iPoint,offset+iVar);
      }
    }
    offset += nVar;
  }
  SU2_OMP_BARRIER
}

void CNewtonIntegration::BlockJacobiPrecond(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const {

  for (unsigned long iSol = 0, offset = 0; iSol < kindSol.size(); ++iSol) {

    const auto nVar = nVars[iSol];
    const auto nPoint = geometry->GetnPoint();

    if (preconditioners[iSol]) {
      /*--- Get work vectors with nVar compatible with the preconditioner. ---*/
      auto& uLoc = GetPrecVecIn<su2double>(iSol);
      auto& vLoc = GetPrecVecOut<su2double>(iSol);

      PARALLEL_FOR
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < nVar; ++iVar)
          uLoc(iPoint, iVar) = u(iPoint, offset+iVar);

      /*--- Apply the preconditioner of this solver. ---*/
      (*preconditioners[iSol])(uLoc, vLoc);

      PARALLEL_FOR
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < nVar; ++iVar)
          v(iPoint, offset+iVar) = vLoc(iPoint, iVar);
    }
    else {
      /*--- Identity. ---*/
      /// TODO: Probably even some type of volume/timestep scaling would be better?
      PARALLEL_FOR
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < nVar; ++iVar)
          v(iPoint, offset+iVar) = u(iPoint, offset+iVar);
    }

    offset += nVar;
  }
  SU2_OMP_BARRIER
}
