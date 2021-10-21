/*!
 * \file CNewtonUpdateOnSubspace.hpp
 * \brief
 * \note Adopts some general layout from CQuasiNewtonInvLeastSquares.
 * \author O. Burghardt
 * \version 7.1.1 "Blackbird"
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

#include "CQuasiNewtonInvLeastSquares.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
//typedef Eigen::Matrix<su2double, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
//typedef Eigen::Matrix<su2double, Eigen::Dynamic, 1> EigenVector;
//#include <Eigen/LU>
//#include <Eigen/SVD>
//typedef Eigen::JacobiSVD<EigenMatrix> EigenSVD;

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
  su2matrix<Scalar> q,q_outer;
  su2matrix<Scalar> p;                  // projected solution
  Eigen::VectorXd p_R;                  // projected solution in terms of basis R
  Eigen::VectorXd pn_R;                 // old projected solution in terms of basis R

  Eigen::MatrixXd EigenX;               // X in contiguous memory (Eigen object)
  Eigen::MatrixXd EigenR;               // R in contiguous memory (Eigen object)
  Eigen::MatrixXd DR;                   // Derivatives of basis vectors (needed for projected Jacobian)
  Eigen::MatrixXd ProjectedJacobian;    // Projected Jacobian to construct Newton step matrix
  Eigen::MatrixXd NewtonInverseMatrix;  // p = p_n + NewtonInverseMatrix*(p-p_n)

  Index iBasis = 0;
  Index BasisSize_n = 0;

  void GlobalReduce(Scalar* valuesIn, Scalar* valuesOut, const int size) {

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(valuesIn, valuesOut, size, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
#else
    for (int i = 0; i < size; ++i)
      valuesOut[i] = valuesIn[i];
#endif
  }

  Scalar DotEigen (const Eigen::VectorXd& left, const Eigen::VectorXd& right) {

    Scalar sum = left.dot(right);
    Scalar sumTotal = 0.0;

    GlobalReduce(&sum, &sumTotal, 1);

    return sumTotal;
  }

  void NormalizeEigen(Eigen::MatrixXd& matrix, int i) {

    double norm = matrix.col(i).dot(matrix.col(i));
    double normTotal = 0.0;

    GlobalReduce(&norm, &normTotal, 1);
    norm = sqrt(normTotal);

    double oneOverNorm = 1.0/norm;
    matrix.col(i) = oneOverNorm*matrix.col(i);
  }

  void shiftHistoryLeft(std::vector<su2matrix<Scalar>> &history) {
    for (Index i = 1; i < history.size(); ++i) {
      /*--- Swap instead of moving to re-use the memory of the first sample.
       *    This is why X and R are not stored as contiguous blocks of mem. ---*/
      std::swap(history[i-1], history[i]);
    }
  }

  bool GramSchmidtQR(const Eigen::MatrixXd& Samples, std::vector<su2matrix<Scalar>>& Basis, su2double KrylovCriterionValue) {

    Eigen::MatrixXd Q;
    Eigen::VectorXd R;
    Q.resize(Samples.rows(),Samples.cols());
    R.resize(Samples.cols());

    /*--- QR decomposition of Samples ---*/

    for (Index i = 0; i < Samples.cols(); ++i) {

      Q.col(i) = Samples.col(i);

      for (Index k = 0; k < i; ++k) {
          double proj = - DotEigen(Q.col(i),Q.col(k));
          Q.col(i) = Q.col(i) + proj*Q.col(k);
      }

      NormalizeEigen(Q,i);
      R(i) = DotEigen(Q.col(i), Samples.col(i));
    }

    /* --- Compute quotients in diagonal of R. ---*/

    std::vector<su2double> Krylov_Criterion_Quotients;
    Krylov_Criterion_Quotients.resize(Samples.cols()-1);

    for (Index i = 0; i < Samples.cols()-1; i++)
      if (R(i+1) != 0.0)
        Krylov_Criterion_Quotients[i] = R(i)/R(i+1);

    /*--- Eventually append first column of Q to the basis. ---*/

    su2double FirstQuotient = abs(Krylov_Criterion_Quotients[0]);

    cout << " current criterion value: " << FirstQuotient << endl;
    if (FirstQuotient > KrylovCriterionValue) {

      cout << "Krylov criterion fulfilled (" << Krylov_Criterion_Quotients[0] << "), appending new basis vector ... ";
      iBasis++;

      /*--- Get reference to new basis vector, extract first column from Q ---*/
      Eigen::Map<Eigen::VectorXd> Eigen_NewR(Basis[iBasis-1].data(),Basis[0].size());
      Eigen_NewR = Q.col(0);

      /*--- Maintain an orthogonal Basis ---*/
//      for (Index i = 0; i < iBasis-1; ++i) {
//        Eigen::Map<Eigen::VectorXd> PrecedingR(R[i].data(),R[0].size());
//        Eigen_NewR = Eigen_NewR - Eigen_NewR.dot(PrecedingR)*PrecedingR;      // CHECK: this might be obsolete
//      }
//      Eigen_NewR.normalize();

      cout << "done." << endl;

      return true;
    }
    else return false;
  }

  void projectWork() {

    /*--- Get source reference ---*/
    Eigen::Map<Eigen::VectorXd> Eigen_work(work.data(),work.size());
    /*--- Get target reference ---*/
    Eigen::Map<Eigen::VectorXd> Eigen_p(p.data(),p.size());

    /* --- Compute projection onto subspace of unstable/slow modes ---*/

    for (Index i = 0; i < EigenR.cols(); ++i)
      p_R(i) = DotEigen(EigenR.col(i), Eigen_work);     // p_R: \xi in original paper

    Eigen_p = EigenR*p_R;                               // p addresses: (uncorrected) projected solution in standard basis (needed for computation of q)
  }

  // Also updates current q (so don't touch outside)
  void computeProjections(bool outer) {

    if (iBasis > 0) {

    /*--- Project solution-to-be-corrected update (loaded at address in work), store at p, coefficients at p_R. ---*/
      projectWork();
//      SU2_OMP_SIMD
      for (Index i = 0; i < work.size(); ++i) {
        work.data()[i] = work.data()[i] - p.data()[i];      // work addresses: q (subtract uncorrected projected solution)
      }
    } else {
      for (Index i = 0; i < p.size(); ++i)
        p.data()[i] = 0.0;
    }

    /*--- Update X to check for new basis elements. ---*/

    /*--- Check for need to shift left. ---*/
    if (iSample+1 == X.size()) {
      shiftHistoryLeft(X);
      iSample--;                                      // X[0].data not needed anymore
    }

    if(outer) std::swap(q,q_outer);                   // q addresses: outer q
//    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) {
      q.data()[i] = work.data()[i] - q.data()[i];     // q addresses: delta q
    }
    std::swap(X[++iSample], q);                       // X[iSample] addresses: delta q, address under q is unused
    std::swap(q, work);                               // q addresses: q ;-)
    if(outer) std::swap(q,q_outer);                   // q_outer addresses: outer q
  }

  void correctUnstableProjection() {

    if (EigenR.cols() > BasisSize_n) {
      pn_R = p_R;
      BasisSize_n = EigenR.cols();
    }

    /*--- Compute update w.r.t. subspace basis (Eq. (5.6) in original paper of Shroff & Keller). ---*/
    p_R = pn_R + NewtonInverseMatrix * (p_R - pn_R);
//    p_R = pn_R + (p_R - pn_R);                        // linear algebra sanity check

    /*--- Compute unstable part w.r.t. standard basis ---*/
    Eigen::Map<Eigen::VectorXd> Eigen_p(p.data(),p.size());
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
    q.resize(npt,nvar);
    q_outer.resize(npt,nvar);
    p.resize(npt,nvar);
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

      if(GramSchmidtQR(EigenX, R, KrylovCriterionValue)) {

        /*--- Update Eigen basis object as iBasis was increased ---*/
        EigenR.resize(EigenR.rows(),iBasis);
        for (Index i = 0; i < iBasis; ++i) {
          EigenR.col(i) = Eigen::VectorXd::Map(R[i].data(),R[0].size());
        }
        p_R.resize(iBasis);
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
        ProjectedJacobian(i,j) = DotEigen(EigenR.col(i),DR.col(j));     // R^T*DR
      cout << j+1 << ", ";
    }
    cout << "...";
    ProjectedJacobian = NewtonInverseMatrix.setIdentity() - ProjectedJacobian;
    NewtonInverseMatrix = ProjectedJacobian.inverse();
    cout << " done." << endl;
  }

  /*!
   * \brief Compute a new approximation.
   * \note To be used after storing the FP result.
   */
  void compute() {

    /*--- save p_R to pn_R to prepare next Newton step ---*/
    pn_R = p_R;                                           // pn_R: z in original paper

    /*--- Compute new projection solution p, its coefficients p_R, and new delta solution q ---*/
    computeProjections(false);

    /*--- Newton correction for the slow/unstable part of the solution update. ---*/
    if (iBasis > 0)
      correctUnstableProjection();

    /*--- Set the corrected new solution at address work, use work2 ---*/
//    SU2_OMP_SIMD
    for (Index i = 0; i < work.size(); ++i) {
      work.data()[i] = q.data()[i] + p.data()[i];         // work addresses: corrected solution (add corrected projected solution)
    }
  }
};
