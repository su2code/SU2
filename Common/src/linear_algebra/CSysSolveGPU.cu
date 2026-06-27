/*!
 * \file CSysSolveGPU.cu
 * \brief CUDA/GPU backend dispatch for Krylov solvers.
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

class CDeviceExpressionContextGuard {
 private:
  bool previous = false;

 public:
  CDeviceExpressionContextGuard() : previous(VecExpr::DeviceExpressionsEnabled()) {
    VecExpr::SetDeviceExpressionsEnabled(true);
  }

  ~CDeviceExpressionContextGuard() { VecExpr::SetDeviceExpressionsEnabled(previous); }

  CDeviceExpressionContextGuard(const CDeviceExpressionContextGuard&) = delete;
  CDeviceExpressionContextGuard& operator=(const CDeviceExpressionContextGuard&) = delete;
};

template <class ScalarType>
class CDeviceFGMRESOps {
 public:
  static constexpr bool UseSolverModGramSchmidt = false;

  bool NestedParallel() const { return false; }

  void PrepareInput(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x, bool x_is_zero) const {
    b.HtDTransfer();
    if (x_is_zero) {
      x = ScalarType(0);
    } else {
      x.HtDTransfer();
    }
  }

  void FinalizeOutput(CSysVector<ScalarType>& x) const { x.DtHTransfer(); }

  ScalarType Norm(const CSysVector<ScalarType>& v) const { return v.GPUNorm(); }

  void AssignNegative(CSysVector<ScalarType>& dst, const CSysVector<ScalarType>& src) const { dst = -src; }

  void Subtract(CSysVector<ScalarType>& dst, const CSysVector<ScalarType>& src) const { dst -= src; }

  void Divide(CSysVector<ScalarType>& dst, ScalarType val) const { dst /= val; }

  ScalarType Dot(const CSysVector<ScalarType>& lhs, const CSysVector<ScalarType>& rhs) const {
    return lhs.GPUDot(rhs);
  }

  void AddScaled(CSysVector<ScalarType>& dst, ScalarType alpha, const CSysVector<ScalarType>& src) const {
    dst += alpha * src;
  }

  void UpdateSolution(unsigned long n, const std::vector<CSysVector<ScalarType>>& basis,
                      const su2vector<ScalarType>& y, CSysVector<ScalarType>& x, bool) const {
    for (unsigned long k = 0; k < n; k++) {
      AddScaled(x, y[k], basis[k]);
    }
  }
};

}  // namespace

#include "../../include/linear_algebra/CSysSolveFGMRES.inl"

template <class ScalarType>
unsigned long CSysSolve<ScalarType>::FGMRES_LinSolver_GPU(const VectorType& b, VectorType& x,
                                                          const ProductType& mat_vec, const PrecondType& precond,
                                                          ScalarType tol, unsigned long m, ScalarType& residual,
                                                          bool monitoring, const CConfig* config) const {
  CGPUSolverBLASContextGuard blas_context;
  CDeviceExpressionContextGuard device_expression_context;
  return FGMRES_LinSolverImpl(b, x, mat_vec, precond, tol, m, residual, monitoring, config,
                              CDeviceFGMRESOps<ScalarType>());
}

/*--- Explicit instantiations for the GPU backend wrapper.
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
