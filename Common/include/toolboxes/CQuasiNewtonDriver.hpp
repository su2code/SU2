/*!
 * \file CQuasiNewtonDriver.hpp
 * \brief Implements a method to accelerate and stabilize the convergence
 * of fixed point iterations, the history of past iterates is used to compute
 * a least squares approximation of the inverse of the Jacobian, which is then
 * used to correct the natural solution of the fixed point iteration.
 * \note Based on the IQN-ILS method, see DOI 10.1007/s11831-013-9085-5 and
 * references therein.
 * \author P. Gomes
 * \version 7.0.5 "Blackbird"
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

#pragma once

#include "../omp_structure.hpp"
#include "../mpi_structure.hpp"
#include "CSymmetricMatrix.hpp"

/*!
 * \class A quasi-Newton fixed-point (FP) accelerator based on IQN-ILS.
 * \note The implementation prioritizes storage, the LS problem is solved
 * via the normal equations as that is easy to make parallel over MPI, it
 * may however be unstable with large sample sizes (compared to QR decomp).
 * Usage: Allocate, set the initial solution (operator (i,j), default is 0),
 * run the FP, log its result ("logFPresult"), compute new solution, use it
 * as the input of the FP, run the FP, etc.
 */
template<class Scalar_t>
class CQuasiNewtonDriver {
public:
  using Scalar = Scalar_t;
  using Index = typename su2matrix<Scalar>::Index;
  static_assert(std::is_floating_point<Scalar>::value,"");

private:
  enum: size_t {BLOCK_SIZE = 1024};     /*!< \brief Loop tiling parameter. */
  std::vector<su2matrix<Scalar> > X, R; /*!< \brief Input and residual history of the FP. */
  su2matrix<Scalar> corr;               /*!< \brief Correction to the natural FP result. */
  su2vector<Scalar> mat, rhs, sol;      /*!< \brief Matrix, rhs, and solution of the normal equations. */
  Index nSample = 0, iSample = 0;       /*!< \brief Max num, and current sample. */
  Index nPtDomain = 0;                  /*!< \brief Local size of the history, considered in dot products. */

  void shiftHistoryLeft() {
    for (Index i = 1; i < nSample; ++i) {
      /*--- Swap instead of moving to re-use the memory of the first sample.
       * This is why X and R are not stored as contiguous blocks of mem. ---*/
      std::swap(X[i-1], X[i]);
      std::swap(R[i-1], R[i]);
    }
  }

  void computeNormalEquations() {
    /*--- Lower triangular packed storage. ---*/
    mat.resize(iSample*(iSample+1)/2) = Scalar(0);
    rhs.resize(iSample) = Scalar(0);

    /*--- Size for the dot products. ---*/
    const auto end = std::min<Index>(nPtDomain,corr.rows())*corr.cols();

    /*--- Tiled part of the loop. ---*/
    Index begin = 0;
    while (end-begin >= BLOCK_SIZE) {
      computeNormalEquations<BLOCK_SIZE>(mat, rhs, begin);
      begin += BLOCK_SIZE;
    }

    /*--- Remainder of the loop. ---*/
    if (begin != end) {
      computeNormalEquations<0>(mat, rhs, begin, end-begin);
    }

    /*--- MPI reduction of the dot products. ---*/
    if (nPtDomain < corr.rows()) {
      const auto type = (sizeof(Scalar) < sizeof(double))? MPI_FLOAT : MPI_DOUBLE;

      su2vector<Scalar> tmp(mat.size());
      SelectMPIWrapper<Scalar>::W::Allreduce(mat.data(), tmp.data(), mat.size(),
                                             type, MPI_SUM, MPI_COMM_WORLD);
      mat = std::move(tmp);

      sol = rhs;
      SelectMPIWrapper<Scalar>::W::Allreduce(sol.data(), rhs.data(), rhs.size(),
                                             type, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  template<Index StaticSize>
  void computeNormalEquations(su2vector<Scalar>& mat,
                              su2vector<Scalar>& vec,
                              Index start,
                              Index dynSize = 0) const {
    /*--- Select either the static or dynamic size, optimizes inner loop. ---*/
    const auto blkSize = StaticSize? StaticSize : dynSize;

    for (Index i = 0; i < iSample; ++i) {
      const auto ri1 = R[i+1].data() + start;
      const auto ri0 = R[i].data() + start;
      for (Index j = 0; j <= i; ++j) {
        const auto rj1 = R[j+1].data() + start;
        const auto rj0 = R[j].data() + start;
        Scalar sum = 0;
        SU2_OMP_SIMD
        for (Index k = 0; k < blkSize; ++k) {
          sum += (ri1[k]-ri0[k]) * (rj1[k]-rj0[k]);
        }
        const auto iCoeff = i*(i+1)/2 + j;
        mat(iCoeff) += sum;
      }
      const auto r = R[iSample].data() + start;
      Scalar sum = 0;
      SU2_OMP_SIMD
      for (Index k = 0; k < blkSize; ++k) {
        sum += (ri1[k]-ri0[k]) * r[k];
      }
      vec(i) -= sum;
    }
  }

public:
  /*!
   * \brief Default construction without allocation.
   */
  CQuasiNewtonDriver() = default;

  /*!
   * \brief Construction with allocation, see "resize".
   */
  CQuasiNewtonDriver(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {
    resize(nsample, npt, nvar, nptdomain);
  }

  /*!
   * \brief Resize the object.
   * \param[in] nsample - Number of samples used to build the FP history.
   * \param[in] npt - Size of the solution including any halos.
   * \param[in] nvar - Number of solution variables.
   * \param[in] nptdomain - Local size (< npt), if 0 (default), MPI parallelization is skipped.
   */
  void resize(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {
    assert(nptdomain <= npt);
    assert(nsample > 1);
    iSample = 0;
    nSample = nsample;
    nPtDomain = nptdomain? nptdomain : npt;
    corr.resize(npt,nvar);
    X.clear();
    R.clear();
    for (Index i = 0; i < nSample; ++i) {
      X.emplace_back(npt,nvar);
      R.emplace_back(npt,nvar);
    }
    X[0] = Scalar(0);
  }

  /*!
   * \brief Discard all history, keeping the current solution.
   */
  void reset() { std::swap(X[0], X[iSample]); iSample = 0; }

  /*!
   * \brief Store the result of the fixed point iteration for given point and variable index.
   */
  void logFPresult(Index iPt, Index iVar, Scalar val) { corr(iPt,iVar) = val; }

  /*!
   * \brief Access entire structure that stores the current FP result.
   */
  su2matrix<Scalar>& getFPresult() { return corr; }
  const su2matrix<Scalar>& getFPresult() const { return corr; }

  /*!
   * \brief Access the current solution approximation.
   */
  Scalar& operator() (Index iPt, Index iVar) { return X[iSample](iPt,iVar); }
  const Scalar& operator() (Index iPt, Index iVar) const { return X[iSample](iPt,iVar); }

  /*!
   * \brief Compute a new approximation.
   * \note To be used after logging the FP result for all points and variables.
   */
  void compute() {
    /*--- Compute FP residual, clear correction. ---*/
    SU2_OMP_SIMD
    for (Index i = 0; i < corr.size(); ++i) {
      R[iSample].data()[i] = corr.data()[i] - X[iSample].data()[i];
      corr.data()[i] = Scalar(0);
    }

    if (iSample > 0) {
      /*--- Solve the normal equations. ---*/
      computeNormalEquations();
      CSymmetricMatrix pseudoInv(iSample);
      for (Index i = 0, k = 0; i < iSample; ++i)
        for (Index j = 0; j <= i; ++j)
          pseudoInv(i,j) = mat(k++);
      pseudoInv.Invert(true);
      sol.resize(iSample);
      pseudoInv.MatVecMult(rhs.data(), sol.data());

      /*--- Compute correction, cleared before for less trunc. error. ---*/
      for (Index k = 0; k < sol.size(); ++k) {
        const auto x1 = X[k+1].data();
        const auto r1 = R[k+1].data();
        const auto x0 = X[k].data();
        const auto r0 = R[k].data();
        SU2_OMP_SIMD
        for (Index i = 0; i < corr.size(); ++i) {
          Scalar dy = r1[i]-r0[i] + x1[i]-x0[i];
          corr.data()[i] += sol(k) * dy;
        }
      }
    }

    /*--- Check for need to shift left. ---*/
    if (iSample == nSample-1) {
      shiftHistoryLeft();
      iSample--;
    }

    /*--- Set new solution. ---*/
    SU2_OMP_SIMD
    for (Index i = 0; i < corr.size(); ++i)
      corr.data()[i] += R[iSample].data()[i] + X[iSample].data()[i];
    std::swap(X[++iSample], corr);
  }
};
