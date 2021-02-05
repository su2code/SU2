/*!
 * \file CNewtonIntegration.hpp
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

#include "CIntegration.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/linear_algebra/CPreconditioner.hpp"
#include "../../../Common/include/linear_algebra/CSysSolve.hpp"

#ifdef HAVE_OMP
#ifdef HAVE_OMP_SIMD
#define CNEWTON_PARFOR SU2_OMP(for simd schedule(static,omp_chunk_size) nowait)
#else
#define CNEWTON_PARFOR SU2_OMP(for schedule(static,omp_chunk_size) nowait)
#endif
#else
#define CNEWTON_PARFOR SU2_OMP_SIMD
#endif

/*!
 * \class CNewtonIntegration
 * \brief Class for time integration using a Newton-Krylov method, based
 * on matrix-free products with the true Jacobian via finite differences.
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
  enum ResEvalType {EXPLICIT, DEFAULT};

  bool setup = false;
  Scalar finDiffStep = 0.0; /*!< \brief Based on RMS(solution), used in matrix-free products. */
  unsigned long omp_chunk_size; /*!< \brief Chunk size used in light point loops. */

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
  }

  /*--- Preconditioner objects for each active solver. ---*/
  CPreconditioner<MixedScalar>* preconditioner = nullptr;

  /*--- If mixed precision is used, these temporaries are used to interface with the preconditioner. ---*/
  mutable CSysVector<MixedScalar> precondIn, precondOut;

  template<class T, su2enable_if<!std::is_same<T,MixedScalar>::value> = 0>
  inline void Preconditioner_impl(const CSysVector<T>& u, CSysVector<T>& v) const {
    CNEWTON_PARFOR
    for (auto i = 0ul; i < u.GetLocSize(); ++i) precondIn[i] = u[i];

    (*preconditioner)(precondIn, precondOut);

    CNEWTON_PARFOR
    for (auto i = 0ul; i < u.GetLocSize(); ++i) v[i] = precondOut[i];
    SU2_OMP_BARRIER
  }

  /*--- Otherwise they are not needed. ---*/
  template<class T, su2enable_if<std::is_same<T,MixedScalar>::value> = 0>
  inline void Preconditioner_impl(const CSysVector<T>& u, CSysVector<T>& v) const { (*preconditioner)(u, v); }

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
