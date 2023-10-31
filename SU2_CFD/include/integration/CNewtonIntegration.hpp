/*!
 * \file CNewtonIntegration.hpp
 * \brief Newton-Krylov integration.
 * \author P. Gomes
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CIntegration.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/linear_algebra/CPreconditioner.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../../Common/include/linear_algebra/CSysSolve.hpp"

#ifdef HAVE_OMP
#ifdef HAVE_OMP_SIMD
#define CNEWTON_PARFOR SU2_OMP_FOR_(simd schedule(static,omp_chunk_size) SU2_NOWAIT)
#else
#define CNEWTON_PARFOR SU2_OMP_FOR_(schedule(static,omp_chunk_size) SU2_NOWAIT)
#endif
#define END_CNEWTON_PARFOR END_SU2_OMP_FOR
#else
#define CNEWTON_PARFOR SU2_OMP_SIMD
#define END_CNEWTON_PARFOR
#endif

/*!
 * \class CNewtonIntegration
 * \ingroup Drivers
 * \brief Class for time integration using a Newton-Krylov method, based
 * on matrix-free products with the true Jacobian via finite differences.
 * \author P. Gomes
 */
class CNewtonIntegration final : public CIntegration {
public:
#ifdef CODI_FORWARD_TYPE
  using Scalar = su2double;
  using MixedScalar = su2double;
#else
  /*--- No point having single precision matrix-free products. ---*/
  using Scalar = passivedouble;
  /*--- The block preconditioners may still use single precision. ---*/
  using MixedScalar = su2mixedfloat;
#endif

private:
  /*--- Residual evaluation modes, explicit for products, default to allow preconditioners to be built. ---*/
  enum class ResEvalType {EXPLICIT, DEFAULT};

  bool setup = false;
  Scalar finDiffStepND = 0.0;
  Scalar finDiffStep = 0.0; /*!< \brief Based on RMS(solution), used in matrix-free products. */
  unsigned long omp_chunk_size; /*!< \brief Chunk size used in light point loops. */

  /*--- Number of iterations and tolerance for the linear preconditioner,
   * 0 iterations forces "weak" preconditioning, i.e. not iterative. ---*/
  unsigned short precondIters = 0;
  su2double precondTol = 0.0;

  /*--- For a number of iterations, or before a certain residual drop,
   * use the quasi-Newton approach instead of Newton-Krylov. If both
   * criteria are zero, or the solver does not provide a linear
   * preconditioner, there is no startup phase. ---*/
  bool startupPeriod = false;
  unsigned short startupIters = 0;
  su2double startupResidual = 0.0;
  su2double firstResidual = -20.0;

  /*--- Relax (increase) the tolerance for NK solves by a factor, until a
   * certain drop in residuals, to reduce the cost of early iterations. ---*/
  unsigned short tolRelaxFactor = 0;
  su2double fullTolResidual = 0.0;

  CConfig* config = nullptr;
  CSolver** solvers = nullptr;
  CGeometry* geometry = nullptr;
  CNumerics*** numerics = nullptr;

  /*--- Residual and linear solver. ---*/
  CSysVector<Scalar> LinSysRes;
  CSysSolve<Scalar> LinSolver;

  /*--- If possible the solution vector of the solver is re-used, otherwise this temporary is used. ---*/
  CSysVector<Scalar> LinSysSol;

  template<class T, su2enable_if<std::is_same<T,Scalar>::value> = 0>
  inline CSysVector<Scalar>& GetSolutionVec(CSysVector<T>& x) { return x; }

  template<class T, su2enable_if<std::is_same<T,Scalar>::value> = 0>
  inline void SetSolutionResult(CSysVector<T>&) const { }

  template<class T, su2enable_if<!std::is_same<T,Scalar>::value> = 0>
  inline CSysVector<Scalar>& GetSolutionVec(CSysVector<T>&) {
    LinSysSol = Scalar(0.0);
    return LinSysSol;
  }

  template<class T, su2enable_if<!std::is_same<T,Scalar>::value> = 0>
  inline void SetSolutionResult(CSysVector<T>& x) const {
    CNEWTON_PARFOR
    for (auto i = 0ul; i < x.GetLocSize(); ++i) x[i] = LinSysSol[i];
    END_CNEWTON_PARFOR
  }

  /*--- Preconditioner objects for each active solver. ---*/
  CPreconditioner<MixedScalar>* preconditioner = nullptr;

  /*--- If mixed precision is used, these temporaries are used to interface with the preconditioner. ---*/
  mutable CSysVector<MixedScalar> precondIn, precondOut;

  template<class T, su2enable_if<!std::is_same<T,MixedScalar>::value> = 0>
  inline unsigned long Preconditioner_impl(const CSysVector<T>& u, CSysVector<T>& v,
                                           unsigned long iters, Scalar& eps) const {
    CNEWTON_PARFOR
    for (auto i = 0ul; i < u.GetLocSize(); ++i) precondIn[i] = u[i];
    END_CNEWTON_PARFOR

    iters = Preconditioner_impl(precondIn, precondOut, iters, eps);

    CNEWTON_PARFOR
    for (auto i = 0ul; i < u.GetLocSize(); ++i) v[i] = precondOut[i];
    END_CNEWTON_PARFOR
    SU2_OMP_BARRIER

    return iters;
  }

  /*--- Otherwise they are not needed. ---*/
  template<class T, su2enable_if<std::is_same<T,MixedScalar>::value> = 0>
  inline unsigned long Preconditioner_impl(const CSysVector<T>& u, CSysVector<T>& v,
                                           unsigned long iters, Scalar& eps) const {
    if (iters == 0) {
      (*preconditioner)(u, v);
      return 0;
    }
    auto product = CSysMatrixVectorProduct<MixedScalar>(solvers[FLOW_SOL]->Jacobian, geometry, config);
    v = MixedScalar(0.0);
    MixedScalar eps_t = eps;
    iters = solvers[FLOW_SOL]->System.FGMRES_LinSolver(u, v, product, *preconditioner, eps, iters, eps_t, false, config);
    eps = eps_t;
    return iters;
  }

  /*!
   * \brief Gather solver info, etc..
   */
  void Setup();

  /*!
   * \brief Increment the solution, x := x+mag*dir.
   */
  void PerturbSolution(const CSysVector<Scalar>& direction, Scalar magnitude);

  /*!
   * \brief Evaluate the nonlinear residual of the solver, which should be capable of alternating
   * between implicit and explicit iterations to save time during matrix-free products.
   */
  void ComputeResiduals(ResEvalType type);

  /*!
   * \brief Compute the step size for finite differences.
   */
  void ComputeFinDiffStep();

public:
  /*!
   * \brief Constructor.
   */
  CNewtonIntegration() = default;

  /*!
   * \brief Destructor.
   */
  ~CNewtonIntegration();

  /*!
   * \brief This class overrides this method to make it a drop-in replacement for CMultigridIntegration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] iZone - Current zone.
   * \param[in] iInst - Current instance.
   */
  void MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                           CNumerics ******numerics_container, CConfig **config,
                           unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) override;

  /*!
   * \brief Implementation of matrix-vector product with the real Jacobian of the nonlinear residuals.
   */
  void MatrixFreeProduct(const CSysVector<Scalar>& u, CSysVector<Scalar>& v);

  /*!
   * \brief Wrapper for the preconditioner.
   */
  void Preconditioner(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const;

};

#undef CNEWTON_PARFOR
#undef END_CNEWTON_PARFOR
