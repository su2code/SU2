/*!
 * \file CNewtonIntegration.hpp
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

#include "CIntegration.hpp"
#include "../../../Common/include/linear_algebra/CSysSolve.hpp"

/*!
 * \class CNewtonIntegration
 * \brief Class for time integration using a coupled Newton method, based
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

  bool thread_safe;             /*!< \brief If all target solvers support OpenMP. */
  unsigned long omp_chunk_size; /*!< \brief Chunk size used in light point loops. */

  unsigned short KindFlowSol = 0;
  unsigned short KindSol2EqSys[MAX_SOLS] = {0}; /*!< \brief Deduce runtime equations from solver position. */

  Scalar finDiffStep = 0.0;     /*!< \brief Based on RMS(solution), used in matrix-free products. */

  CConfig* config = nullptr;
  CSolver** solvers = nullptr;
  CGeometry* geometry = nullptr;
  CNumerics*** numerics = nullptr;

  std::vector<unsigned> kindSol; /*!< \brief Positions of the target solvers. */
  std::vector<unsigned> nVars;   /*!< \brief Number of variables for each target solver. */

  /*--- Residual, solution, and linear solver for the coupled problem. ---*/
  CSysVector<Scalar> LinSysRes;
  CSysVector<Scalar> LinSysSol;
  CSysSolve<Scalar> LinSolver;

  /*--- Preconditioner objects for each active solver. ---*/
  std::vector<CPreconditioner<MixedScalar>*> preconditioners;

  /*--- If mixed precision is used these temporaries are
   *    used to interface with the preconditioners. ---*/
  mutable std::vector<CSysVector<MixedScalar> > precondIn, precondOut;

  template<class T, su2enable_if<!std::is_same<T,MixedScalar>::value> = 0>
  inline CSysVector<MixedScalar>& GetPrecVecIn(size_t i) const { return precondIn[i]; }

  template<class T, su2enable_if<!std::is_same<T,MixedScalar>::value> = 0>
  inline CSysVector<MixedScalar>& GetPrecVecOut(size_t i) const { return precondOut[i]; }

  /*--- Otherwise we borrow the memory of the solvers. ---*/
  template<class T, su2enable_if<std::is_same<T,MixedScalar>::value> = 0>
  inline CSysVector<MixedScalar>& GetPrecVecIn(size_t i) const { return solvers[kindSol[i]]->LinSysRes; }

  template<class T, su2enable_if<std::is_same<T,MixedScalar>::value> = 0>
  inline CSysVector<MixedScalar>& GetPrecVecOut(size_t i) const { return solvers[kindSol[i]]->LinSysSol; }

  /*!
   * \brief List target solvers, gather their info, etc..
   */
  void Setup();

  /*!
   * \brief Evaluate the nonlinear residual of target solvers, which should be capable of alternating
   * between implicit and explicit iterations to save time during matrix-free products.
   */
  void ComputeResiduals(ResEvalType type);

public:
  /*!
   * \brief Constructor.
   */
  CNewtonIntegration();

  /*!
   * \brief Destructor.
   */
  ~CNewtonIntegration();

  /*!
   * \brief Return true if the integration already considers all solvers.
   */
  inline bool IsFullyCoupled(void) const override { return true; }

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
   * \brief Implementation of the block-Jacobi preconditioner.
   */
  void BlockJacobiPrecond(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const;

};
