/*!
 * \file CQuasiNewtonInvLeastSquares.hpp
 * \brief Implements a method to accelerate and stabilize the convergence
 * of fixed point iterations, the history of past iterates is used to compute
 * a least squares approximation of the inverse of the Jacobian, which is then
 * used to correct the natural solution of the fixed point iteration.
 * \note Based on the IQN-ILS method, see DOI 10.1007/s11831-013-9085-5 and
 * references therein.
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

#pragma once

#include "../parallelization/omp_structure.hpp"
#include "../parallelization/mpi_structure.hpp"
#include "CSymmetricMatrix.hpp"

//#ifdef HAVE_EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
typedef Eigen::Matrix<su2double, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
typedef Eigen::Matrix<su2double, Eigen::Dynamic, 1> EigenVector;
//#endif
//#include <Eigen/LU>
//#include <Eigen/SVD>
//typedef Eigen::JacobiSVD<EigenMatrix> EigenSVD;

/*!
 * \class A quasi-Newton fixed-point (FP) accelerator based on IQN-ILS.
 * \note The implementation prioritizes storage, the LS problem is solved
 * via the normal equations as that is easy to make parallel over MPI, it
 * may however be unstable with large sample sizes (compared to QR decomp).
 * Usage: Allocate, store the initial solution (operator (i,j), default is 0),
 * run the FP, store its result ("FPresult"), compute new solution, use it
 * as the new input of the FP, run the FP, etc.
 */
template<class Scalar_t>
class CQuasiNewtonInvLeastSquares {
public:
  using Scalar = Scalar_t;
  using Index = typename su2matrix<Scalar>::Index;
  static_assert(std::is_floating_point<Scalar>::value,"");

protected:                                                        // TODO: revert this back to private
  using MPI_Wrapper = typename SelectMPIWrapper<Scalar>::W;

  enum: size_t {BLOCK_SIZE = 1024};     /*!< \brief Loop tiling parameter. */
  std::vector<su2matrix<Scalar> > X, R; /*!< \brief Input and residual history of the FP. */
  su2matrix<Scalar> work;               /*!< \brief Work matrix (FP result, correction, and approx solution). */
  su2vector<Scalar> mat, rhs, sol;      /*!< \brief Matrix, rhs, and solution of the normal equations. */
  Index iSample = 0;                    /*!< \brief Current sample index. */
  Index nPtDomain = 0;                  /*!< \brief Local size of the history, considered in dot products. */

  void shiftHistoryLeft() {
    for (Index i = 1; i < X.size(); ++i) {
      /*--- Swap instead of moving to re-use the memory of the first sample.
       *    This is why X and R are not stored as contiguous blocks of mem. ---*/
      std::swap(X[i-1], X[i]);
      std::swap(R[i-1], R[i]);
    }
  }

  void computeNormalEquations() {
    /*--- Size for the dot products. ---*/
    const auto end = std::min<Index>(nPtDomain,work.rows())*work.cols();

    mat = Scalar(0); rhs = Scalar(0);

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
    if (nPtDomain < work.rows()) {
      const auto type = (sizeof(Scalar) < sizeof(double))? MPI_FLOAT : MPI_DOUBLE;

      su2vector<Scalar> tmp(mat.size());
      MPI_Wrapper::Allreduce(mat.data(), tmp.data(), iSample*(iSample+1)/2,
                             type, MPI_SUM, MPI_COMM_WORLD);
      mat = std::move(tmp);

      MPI_Wrapper::Allreduce(rhs.data(), sol.data(), iSample,
                             type, MPI_SUM, MPI_COMM_WORLD);
      std::swap(rhs, sol);
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
      /*--- Off-diagonal coefficients. ---*/
      for (Index j = 0; j < i; ++j) {
        const auto rj1 = R[j+1].data() + start;
        const auto rj0 = R[j].data() + start;
        /*--- Sum of partial sums to reduce trunc. error. ---*/
        Scalar sum = 0;
        SU2_OMP_SIMD
        for (Index k = 0; k < blkSize; ++k) {
          sum += (ri1[k]-ri0[k]) * (rj1[k]-rj0[k]);
        }
        /*--- 1D index of (i,j) in lower triangular storage. ---*/
        const auto iCoeff = i*(i+1)/2 + j;
        mat(iCoeff) += sum;
      }
      /*--- Diagonal coeff and residual fused. ---*/
      Scalar diag = 0, res = 0;
      const auto r = R[iSample].data() + start;
      SU2_OMP_SIMD
      for (Index k = 0; k < blkSize; ++k) {
        diag += pow(ri1[k]-ri0[k], 2);
        res += (ri1[k]-ri0[k]) * r[k];
      }
      mat(i*(i+3)/2) += diag;
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
   * \param[in] nptdomain - Local size (< npt), if 0 (default), MPI parallelization is skipped.
   */
  void resize(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {
    if (nptdomain > npt || nsample < 2)
      SU2_MPI::Error("Invalid quasi-Newton parameters", CURRENT_FUNCTION);
    iSample = 0;
    nPtDomain = nptdomain? nptdomain : npt;
    work.resize(npt,nvar);
    X.clear();
    R.clear();
    for (Index i = 0; i < nsample; ++i) {
      X.emplace_back(npt,nvar);
      R.emplace_back(npt,nvar);
    }
    X[0] = Scalar(0);
    /*--- Lower triangular packed storage. ---*/
    mat.resize(nsample*(nsample-1)/2);
    rhs.resize(nsample-1);
    sol.resize(nsample-1);
  }

  /*! \brief Size of the object, the number of samples. */
  Index size() const { return X.size(); }

  /*! \brief Discard all history, keeping the current solution. */
  void reset() { std::swap(X[0], X[iSample]); iSample = 0; }

  /*!
   * \brief Access the current fixed-point result.
   * \note Use these to STORE the result of running the FP.
   */
  su2matrix<Scalar>& FPresult() { return work; }
  const su2matrix<Scalar>& FPresult() const { return work; }
  Scalar& FPresult(Index iPt, Index iVar) { return work(iPt,iVar); }
  const Scalar& FPresult(Index iPt, Index iVar) const { return work(iPt,iVar); }

  /*!
   * \brief Access the current solution approximation.
   * \note Use these, after calling compute, to GET the new FP solution estimate.
   */
  su2matrix<Scalar>& solution() { return X[iSample]; }
  const su2matrix<Scalar>& solution() const { return X[iSample]; }
  Scalar& operator() (Index iPt, Index iVar) { return X[iSample](iPt,iVar); }
  const Scalar& operator() (Index iPt, Index iVar) const { return X[iSample](iPt,iVar); }

  /*!
   * \brief Compute and return a new approximation.
   * \note To be used after storing the FP result.
   */
  const su2matrix<Scalar>& compute() {
    /*--- Compute FP residual, clear correction. ---*/
    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) {
      /*--- Compute difference between the uncorrected and the last corrected solution. ---*/
      R[iSample].data()[i] = work.data()[i] - X[iSample].data()[i];
      work.data()[i] = Scalar(0);
    }

    if (iSample > 0) {
      /*--- Solve the normal equations. ---*/
      computeNormalEquations();
      CSymmetricMatrix pseudoInv(iSample);
      for (Index i = 0, k = 0; i < iSample; ++i)
        for (Index j = 0; j <= i; ++j)
          pseudoInv(i,j) = mat(k++);
      pseudoInv.Invert(true);
      pseudoInv.MatVecMult(rhs.data(), sol.data());

      /*--- Compute correction, cleared before for less trunc. error. ---*/
      for (Index k = 0; k < iSample; ++k) {
        const auto x1 = X[k+1].data();
        const auto r1 = R[k+1].data();
        const auto x0 = X[k].data();
        const auto r0 = R[k].data();
        SU2_OMP_SIMD
        for (Index i = 0; i < work.size(); ++i) {
          Scalar dy = r1[i]-r0[i] + x1[i]-x0[i];
          work.data()[i] += sol(k) * dy;
        }
      }
    }

    /*--- Check for need to shift left. ---*/
    if (iSample+1 == X.size()) {
      shiftHistoryLeft();
      iSample--;
    }

    /*--- Set new solution. ---*/
    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i)
      work.data()[i] += R[iSample].data()[i] + X[iSample].data()[i];
    std::swap(X[++iSample], work);

    return solution();
  }
};


/*--- Both classes operate on a window of past corrected solutions (X) and an
 *    data structure (R) of similar size to construct a (quasi-) Newton scheme.
 *    Have to identify a proper common base class. ---*/
template<class Scalar_t>
class CNewtonUpdateOnSubspace : public CQuasiNewtonInvLeastSquares<Scalar_t> {

public:

  using Scalar = Scalar_t;
  using Index = typename su2matrix<Scalar>::Index;

protected:

  using MPI_Wrapper = typename SelectMPIWrapper<Scalar>::W;

  using CQuasiNewtonInvLeastSquares<Scalar_t>::iSample;
  using CQuasiNewtonInvLeastSquares<Scalar_t>::nPtDomain;
  using CQuasiNewtonInvLeastSquares<Scalar_t>::work;        // the uncorrected (input) solutions and intermediate data
  using CQuasiNewtonInvLeastSquares<Scalar_t>::X;           // deltas of corrected solutions, nsample many
  using CQuasiNewtonInvLeastSquares<Scalar_t>::R;           // basis of unstable space, nbasis many

  /*--- Some more storage is needed because of the handling of two separated solution update strategies ---*/
  su2matrix<Scalar> work2;
  su2matrix<Scalar> p;                  // projected solution
  su2matrix<Scalar> p_n;                // old projected solution
  Eigen::VectorXd p_R;                  // projected solution in terms of basis R
  Eigen::VectorXd pn_R;                 // old projected solution in terms of basis R

  Eigen::MatrixXd EigenX;               // X in contiguous memory (Eigen object)
  Eigen::MatrixXd EigenR;               // R in contiguous memory (Eigen object)
  Eigen::MatrixXd DR;                   // Derivatives of basis vectors (needed for projected Jacobian)
  Eigen::MatrixXd ProjectedJacobian;    // Projected Jacobian to construct Newton step matrix
  Eigen::MatrixXd NewtonInverseMatrix;  // p = p_n + NewtonInverseMatrix*(p-p_n)

  Index iBasis = 0;

  void shiftHistoryLeft(std::vector<su2matrix<Scalar>> &history) {
    for (Index i = 1; i < history.size(); ++i) {
      /*--- Swap instead of moving to re-use the memory of the first sample.
       *    This is why X and R are not stored as contiguous blocks of mem. ---*/
      std::swap(history[i-1], history[i]);
    }
  }

  void projectOntoSubspace() {

    /*--- save p to p_n as needed for the Newton step ---*/
    std::swap(p,p_n);                                   // p_n addresses: old projected solution

    /*--- Get references ---*/
    Eigen::Map<Eigen::VectorXd> Eigen_work(work.data(),work.size());
    Eigen::Map<Eigen::VectorXd> Eigen_p(p.data(),p.size());

    /* --- Compute projection onto subspace of unstable/slow modes ---*/
    p_R = EigenR.transpose()*Eigen_work;                // \xi in original paper
    Eigen_p = EigenR*p_R;                               // p addresses: projected solution
  }

  void updateProjectedSolution() {

    /*--- Get references ---*/
    Eigen::Map<Eigen::VectorXd> Eigen_p(p.data(),p.size());
    Eigen::Map<Eigen::VectorXd> Eigen_pn(p_n.data(),p_n.size());

    pn_R = EigenR.transpose()*Eigen_pn;                 // z in original paper

    /*--- Compute update w.r.t. subspace basis (Eq. (5.6) in original paper of Shroff & Keller). ---*/
    p_R = pn_R + NewtonInverseMatrix * (p_R - pn_R);
//    p_R = pn_R + (p_R - pn_R);                        // linear algebra sanity check
    Eigen_p = EigenR*p_R;                               // updated projected solution
  }

public:

  /*! \brief Default construction without allocation. */
  CNewtonUpdateOnSubspace() = default;

  /*! \brief Construction with allocation, see "resize". */
  CNewtonUpdateOnSubspace(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {
    resize(nsample, npt, nvar, nptdomain);
  }

  /*!
   * \brief Resize the object.
   * \param[in] nsample - Number of samples used to build the FP history.
   * \param[in] nbasis - Dimension of basis the unstable space on which we apply the Newton update scheme.
   * \param[in] npt - Size of the solution including any halos.
   * \param[in] nvar - Number of solution variables.
   * \param[in] nptdomain - Local size (< npt), if 0 (default), MPI parallelization is skipped.
   */
  void resize(Index nsample, Index nbasis, Index npt, Index nvar, Index nptdomain = 0) {

    if (nptdomain > npt || nsample < 2)
      SU2_MPI::Error("Invalid Newton update parameters", CURRENT_FUNCTION);

    iSample = 0; iBasis = 0;
    nPtDomain = nptdomain? nptdomain : npt;
    work.resize(npt,nvar);
    work2.resize(npt,nvar);
    p.resize(npt,nvar);
    p_n.resize(npt,nvar);
    X.clear();                              // role here: store history of delta solutions in stable space
    R.clear();                              // role here: store basis of unstable subspace
    for (Index i = 0; i < nsample; ++i) {
      X.emplace_back(npt,nvar);
    }
    for (Index i = 0; i < nbasis; ++i) {
      R.emplace_back(npt,nvar);
    }
    X[0] = Scalar(0);
    R[0] = Scalar(0);

    EigenX.resize(npt*nvar,nsample);
    EigenR.resize(npt*nvar,1);
    DR.resize(npt*nvar,1);
  }

  /*! \brief Size of the object, the size of the subspace basis. */
  Index size() const { return R.size(); }

  /*! \brief Discard all history, keeping the current sample. */
  void reset() { std::swap(X[0], X[iSample]); iSample = 0; iBasis = 0; }

  /*!
   * \brief Check for new basis vector and eventually append to basis.
   */
  bool checkBasis(su2double KrylovCriterionValue) {

    /*--- Check whether we have collected enough samples, if not, return directly. ---*/
    if (iSample+1 < X.size())
      return false;

    if (iBasis < R.size()) {

      /*--- Create Eigen data structures for QR decomposition via Eigen ---*/
      for (Index i = 0; i < X.size(); ++i) {

        /*--- X is not stored in contiguous memory, copying it to an Eigen object is one
         *    alternative... (If it was, something like Map<MatrixXd> X(data, rows, cols)
         *    should be possible.) ---*/
        EigenX.col(i) = Eigen::VectorXd::Map(X[i].data(),X[0].size());
      }

      /* --- Compute QR decomposition and QR criterion ---*/
      Eigen::HouseholderQR<Eigen::MatrixXd> QR(EigenX.rows(),EigenX.cols());
      QR.compute(EigenX);
      auto Rdiag = QR.matrixQR().diagonal();

      std::vector<su2double> Krylov_Criterion_Quotients;
      Krylov_Criterion_Quotients.resize(X.size()-1);

      for (Index i = 0; i < X.size()-1; i++)
        if (Rdiag(i+1) != 0.0)
          Krylov_Criterion_Quotients[i] = Rdiag(i)/Rdiag(i+1);

      if ((abs(Krylov_Criterion_Quotients[0]) > KrylovCriterionValue) &&
          !(abs(Krylov_Criterion_Quotients[0])!=abs(Krylov_Criterion_Quotients[0]))
          ) {

        cout << "Krylov criterion fulfilled (" << Krylov_Criterion_Quotients[0] << "), appending new basis vector ... ";
        iBasis++;

        /*--- Get reference to new basis vector, extract first column from Q ---*/
        Eigen::Map<Eigen::VectorXd> Eigen_NewR(R[iBasis-1].data(),R[0].size());
        Eigen_NewR = QR.householderQ()*Eigen_NewR.setIdentity();

        for (Index i = 0; i < iBasis-1; ++i) {
          Eigen::Map<Eigen::VectorXd> PrecedingR(R[i].data(),R[0].size());
          Eigen_NewR = Eigen_NewR - Eigen_NewR.dot(PrecedingR)*PrecedingR;      // CHECK: this might be obsolete
        }
        Eigen_NewR.normalize();

        /*--- Update Eigen basis object ---*/
        EigenR.resize(EigenR.rows(),iBasis);
        for (Index i = 0; i < iBasis; ++i) {
          EigenR.col(i) = Eigen::VectorXd::Map(R[i].data(),R[0].size());
        }

        cout << "done." << endl;
        return true;
      }
    }
    else {
      // TODO: Maintain the basis
    }
    return false;
  }

  /*!
   * \brief Compute new projected subspace Jacobian and the inverse matrix for Newton steps.
   * \note To be used directly after basis dimension has been increased.
   */
  void computeProjectedJacobian(unsigned short iZone, su2matrix<int>& InputIndices, su2matrix<int>& OutputIndices) {

    ProjectedJacobian.resize(iBasis,iBasis);
    NewtonInverseMatrix.resize(iBasis,iBasis);
    DR.resize(DR.rows(),iBasis);

    cout << "Evaluate R^T (dG/du)^T R[i] for i = ";
    for (Index j = 0; j < iBasis; ++j) {

      AD::ClearAdjoints();
      for (Index it = 0; it < R[j].size(); ++it) {
        AD::SetDerivative(OutputIndices.data()[it], SU2_TYPE::GetValue(R[j].data()[it]));
      }
      AD::ComputeAdjoint();                                             // TODO: make this more efficient

      for (Index it = 0; it < R[j].size(); ++it)
        DR.col(j)(it) = AD::GetDerivative(InputIndices.data()[it]);     // extract DR = (dG/du)^T*R[j]
      for (Index i = 0; i < iBasis; ++i)
        ProjectedJacobian(i,j) = EigenR.col(i).transpose()*DR.col(j);   // R^T*DR

      cout << j+1 << ", ";
    }
    cout << "...";
    ProjectedJacobian = NewtonInverseMatrix.setIdentity() - ProjectedJacobian;
    NewtonInverseMatrix = ProjectedJacobian.inverse();
    cout << " done." << endl;
  }

  /*!
   * \brief Compute and return a new approximation.
   * \note To be used after storing the FP result.
   */
  const su2matrix<Scalar>& compute() {

    if (iBasis > 0) {

    /*--- Project solution update (loaded at work), store at p. ---*/
      projectOntoSubspace();
//      SU2_OMP_SIMD
      for (Index i = 0; i < work.size(); ++i) {
        work.data()[i] = work.data()[i] - p.data()[i];              // work addresses: q
      }
    } else {
      for (Index i = 0; i < p.size(); ++i)
        p.data()[i] = 0.0;
    }

    /*--- Check for need to shift left. ---*/
    if (iSample+1 == X.size()) {
      shiftHistoryLeft(X);
      iSample--;                                                    // X[0].data not needed anymore
    }
    /*--- Keep X updated to check for new basis elements. ---*/
//    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) {
      work2.data()[i] = work.data()[i] - work2.data()[i];           // work2 addresses: delta q
    }
    std::swap(X[++iSample], work2);                                 // X[iSample] addresses: delta q, address under work2 is free
    std::swap(work2, work);                                         // work2 addresses: q

    /*--- Newton correction for the slow/unstable part of the solution update. ---*/
    if (iBasis > 0)
      updateProjectedSolution();

    /*--- Set the corrected new solution in work2---*/
//    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) {
      work.data()[i] = work2.data()[i] + p.data()[i];               // work addresses: corrected solution
    }

    return CQuasiNewtonInvLeastSquares<Scalar_t>::FPresult();
  }
};
