/*!
 * \file CQuasiNewtonInvLeastSquares.hpp
 * \brief Implements a method to accelerate and stabilize the convergence
 * of fixed point iterations, the history of past iterates is used to compute
 * a least squares approximation of the inverse of the Jacobian, which is then
 * used to correct the natural solution of the fixed point iteration.
 * \note Based on the IQN-ILS method, see DOI 10.1007/s11831-013-9085-5 and
 * references therein.
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

#pragma once

#include "../parallelization/omp_structure.hpp"
#include "../parallelization/mpi_structure.hpp"
#include "CSymmetricMatrix.hpp"

/*!
 * \brief A quasi-Newton fixed-point (FP) accelerator based on IQN-ILS.
 * \note The implementation prioritizes storage, the LS problem is solved
 * via the normal equations as that is easy to make parallel over MPI, it
 * may however be unstable with large sample sizes (compared to QR decomp).
 * Usage: Allocate, store the initial solution (operator (i,j), default is 0),
 * run the FP, store its result ("FPresult"), compute new solution, use it
 * as the new input of the FP, run the FP, etc.
 * \ingroup BLAS
 */
template <class Scalar_t, bool WithMPI = true>
class CQuasiNewtonInvLeastSquares {
 public:
  using Scalar = Scalar_t;
  using Index = typename su2matrix<Scalar>::Index;
  static_assert(std::is_floating_point<Scalar>::value, "");

 private:
  using MPI_Wrapper = typename SelectMPIWrapper<Scalar>::W;

  enum : size_t { BLOCK_SIZE = 1024 };  /*!< \brief Loop tiling parameter. */
  std::vector<su2matrix<Scalar> > X, R; /*!< \brief Input and residual history of the FP. */
  su2matrix<Scalar> work;               /*!< \brief Work matrix (FP result, correction, and approx solution). */
  su2vector<Scalar> mat, rhs, sol;      /*!< \brief Matrix, rhs, and solution of the normal equations. */
  Index iSample = 0;                    /*!< \brief Current sample index. */
  Index nPtDomain = 0;                  /*!< \brief Local size of the history, considered in dot products. */

  void shiftHistoryLeft() {
    for (Index i = 1; i < X.size(); ++i) {
      /*--- Swap instead of moving to re-use the memory of the first sample.
       * This is why X and R are not stored as contiguous blocks of mem. ---*/
      std::swap(X[i - 1], X[i]);
      std::swap(R[i - 1], R[i]);
    }
  }

  void computeNormalEquations() {
    /*--- Size for the dot products. ---*/
    const auto end = std::min<Index>(nPtDomain, work.rows()) * work.cols();

    mat = Scalar(0);
    rhs = Scalar(0);

    /*--- Tiled part of the loop. ---*/
    Index begin = 0;
    while (end - begin >= BLOCK_SIZE) {
      computeNormalEquations<BLOCK_SIZE>(mat, rhs, begin);
      begin += BLOCK_SIZE;
    }

    /*--- Remainder of the loop. ---*/
    if (begin != end) {
      computeNormalEquations<0>(mat, rhs, begin, end - begin);
    }

    /*--- MPI reduction of the dot products. ---*/
    if (WithMPI) {
      const auto type = (sizeof(Scalar) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;

      su2vector<Scalar> tmp(mat.size());
      MPI_Wrapper::Allreduce(mat.data(), tmp.data(), iSample * (iSample + 1) / 2, type, MPI_SUM, SU2_MPI::GetComm());
      mat = std::move(tmp);

      MPI_Wrapper::Allreduce(rhs.data(), sol.data(), iSample, type, MPI_SUM, SU2_MPI::GetComm());
      std::swap(rhs, sol);
    }
  }

  template <Index StaticSize>
  void computeNormalEquations(su2vector<Scalar>& mat, su2vector<Scalar>& vec, Index start, Index dynSize = 0) const {
    /*--- Select either the static or dynamic size, optimizes inner loop. ---*/
    const auto blkSize = StaticSize ? StaticSize : dynSize;

    for (Index i = 0; i < iSample; ++i) {
      const auto ri1 = R[i + 1].data() + start;
      const auto ri0 = R[i].data() + start;
      /*--- Off-diagonal coefficients. ---*/
      for (Index j = 0; j < i; ++j) {
        const auto rj1 = R[j + 1].data() + start;
        const auto rj0 = R[j].data() + start;
        /*--- Sum of partial sums to reduce trunc. error. ---*/
        Scalar sum = 0;
        SU2_OMP_SIMD
        for (Index k = 0; k < blkSize; ++k) {
          sum += (ri1[k] - ri0[k]) * (rj1[k] - rj0[k]);
        }
        /*--- 1D index of (i,j) in lower triangular storage. ---*/
        const auto iCoeff = i * (i + 1) / 2 + j;
        mat(iCoeff) += sum;
      }
      /*--- Diagonal coeff and residual fused. ---*/
      Scalar diag = 0, res = 0;
      const auto r = R[iSample].data() + start;
      SU2_OMP_SIMD
      for (Index k = 0; k < blkSize; ++k) {
        diag += pow(ri1[k] - ri0[k], 2);
        res += (ri1[k] - ri0[k]) * r[k];
      }
      mat(i * (i + 3) / 2) += diag;
      vec(i) -= res;
    }
  }

 public:
  /*! \brief Default construction without allocation. */
  CQuasiNewtonInvLeastSquares() = default;

  /*! \brief Construction with allocation, see "resize". */
  CQuasiNewtonInvLeastSquares(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {
    resize(nsample, npt, nvar, nptdomain);
  }

  /*!
   * \brief Resize the object.
   * \param[in] nsample - Number of samples used to build the FP history.
   * \param[in] npt - Size of the solution including any halos.
   * \param[in] nvar - Number of solution variables.
   * \param[in] nptdomain - Local size (<= npt), if 0 it defaults to npt.
   */
  void resize(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {
    if (nptdomain > npt || nsample < 2) SU2_MPI::Error("Invalid quasi-Newton parameters", CURRENT_FUNCTION);
    iSample = 0;
    nPtDomain = nptdomain ? nptdomain : npt;
    work.resize(npt, nvar);
    X.clear();
    R.clear();
    for (Index i = 0; i < nsample; ++i) {
      X.emplace_back(npt, nvar);
      R.emplace_back(npt, nvar);
    }
    X[0] = Scalar(0);
    /*--- Lower triangular packed storage. ---*/
    mat.resize(nsample * (nsample - 1) / 2);
    rhs.resize(nsample - 1);
    sol.resize(nsample - 1);
  }

  /*! \brief Size of the object, the number of samples. */
  Index size() const { return X.size(); }

  /*! \brief Discard all history, keeping the current solution. */
  void reset() {
    std::swap(X[0], X[iSample]);
    iSample = 0;
  }

  /*!
   * \brief Access the current fixed-point result.
   * \note Use these to STORE the result of running the FP.
   */
  su2matrix<Scalar>& FPresult() { return work; }
  const su2matrix<Scalar>& FPresult() const { return work; }
  Scalar& FPresult(Index iPt, Index iVar) { return work(iPt, iVar); }
  const Scalar& FPresult(Index iPt, Index iVar) const { return work(iPt, iVar); }

  /*!
   * \brief Access the current solution approximation.
   * \note Use these, after calling compute, to GET the new FP solution estimate.
   */
  su2matrix<Scalar>& solution() { return X[iSample]; }
  const su2matrix<Scalar>& solution() const { return X[iSample]; }
  Scalar& operator()(Index iPt, Index iVar) { return X[iSample](iPt, iVar); }
  const Scalar& operator()(Index iPt, Index iVar) const { return X[iSample](iPt, iVar); }

  /*!
   * \brief Compute and return a new approximation.
   * \note To be used after storing the FP result.
   */
  const su2matrix<Scalar>& compute() {
    /*--- Compute FP residual, clear correction. ---*/
    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) {
      R[iSample].data()[i] = work.data()[i] - X[iSample].data()[i];
      work.data()[i] = Scalar(0);
    }

    if (iSample > 0) {
      /*--- Solve the normal equations. ---*/
      computeNormalEquations();
      CSymmetricMatrix pseudoInv(iSample);
      for (Index i = 0, k = 0; i < iSample; ++i)
        for (Index j = 0; j <= i; ++j) pseudoInv(i, j) = mat(k++);
      pseudoInv.Invert(true);
      pseudoInv.MatVecMult(rhs.data(), sol.data());

      /*--- Compute correction, cleared before for less trunc. error. ---*/
      for (Index k = 0; k < iSample; ++k) {
        const auto x1 = X[k + 1].data();
        const auto r1 = R[k + 1].data();
        const auto x0 = X[k].data();
        const auto r0 = R[k].data();
        SU2_OMP_SIMD
        for (Index i = 0; i < work.size(); ++i) {
          Scalar dy = r1[i] - r0[i] + x1[i] - x0[i];
          work.data()[i] += sol(k) * dy;
        }
      }
    }

    /*--- Check for need to shift left. ---*/
    if (iSample + 1 == X.size()) {
      shiftHistoryLeft();
      iSample--;
    }

    /*--- Set new solution. ---*/
    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) work.data()[i] += R[iSample].data()[i] + X[iSample].data()[i];
    std::swap(X[++iSample], work);

    return solution();
  }
};
