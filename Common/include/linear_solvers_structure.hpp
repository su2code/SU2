/*!
 * \file linear_solvers_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "vector_structure.hpp"

#ifndef NO_MPI
#include <mpi.h>
#endif
#include <climits>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

/*!
 * \class CSysSolve
 * \brief Class for solving linear systems using classical and Krylov-subspace iterative methods
 * \author J. Hicken.
 * \version 2.0.6
 *
 * The individual solvers could be stand-alone subroutines, but by
 * creating CSysSolve objects we can more easily assign different
 * matrix-vector products and preconditioners to different problems
 * that may arise in a hierarchical solver (i.e. multigrid).
 */
class CSysSolve {
  
private:
  
  /*!
   * \brief sign transfer function
   * \param[in] x - value having sign prescribed
   * \param[in] y - value that defined the sign
   *
   * this may already be defined as a global function somewhere, if
   * so, feel free to delete this and replace it as needed with the
   * appropriate global function
   */
  double sign(const double & x, const double & y) const;
  
  /*!
   * \brief applys a Givens rotation to a 2-vector
   * \param[in] s - sine of the Givens rotation angle
   * \param[in] c - cosine of the Givens rotation angle
   * \param[in,out] h1 - first element of 2x1 vector being transformed
   * \param[in,out] h2 - second element of 2x1 vector being transformed
   */
  void applyGivens(const double & s, const double & c, double & h1, double & h2);
  
  /*!
   * \brief generates the Givens rotation matrix for a given 2-vector
   * \param[in,out] dx - element of 2x1 vector being transformed
   * \param[in,out] dy - element of 2x1 vector being set to zero
   * \param[in,out] s - sine of the Givens rotation angle
   * \param[in,out] c - cosine of the Givens rotation angle
   *
   * Based on givens() of SPARSKIT, which is based on p.202 of
   * "Matrix Computations" by Golub and van Loan.
   */
  void generateGivens(double & dx, double & dy, double & s, double & c);
  
  /*!
   * \brief finds the solution of the upper triangular system Hsbg*x = rhs
   *
   * \param[in] n - size of the reduced system
   * \param[in] Hsbg - upper triangular matrix
   * \param[in] rhs - right-hand side of the reduced system
   * \param[out] x - solution of the reduced system
   *
   * \pre the upper Hessenberg matrix has been transformed into a
   * triangular matrix.
   */
  void solveReduced(const int & n, const vector<vector<double> > & Hsbg,
                    const vector<double> & rhs, vector<double> & x);
  
  /*!
   * \brief Modified Gram-Schmidt orthogonalization
   * \author Based on Kesheng John Wu's mgsro subroutine in Saad's SPARSKIT
   *
   * \tparam Vec - a generic vector class
   * \param[in] i - index indicating which vector in w is being orthogonalized
   * \param[in,out] Hsbg - the upper Hessenberg begin updated
   * \param[in,out] w - the (i+1)th vector of w is orthogonalized against the
   *                    previous vectors in w
   *
   * \pre the vectors w[0:i] are orthonormal
   * \post the vectors w[0:i+1] are orthonormal
   *
   * Reothogonalization is performed if the cosine of the angle between
   * w[i+1] and w[k], k < i+1, is greater than 0.98.  The norm of the "new"
   * vector is kept in nrm0 and updated after operating with each vector
   *
   */
  void modGramSchmidt(int i, vector<vector<double> > & Hsbg, vector<CSysVector> & w);
  
  /*!
   * \brief writes header information for a CSysSolve residual history
   * \param[in,out] os - ostream class object for output
   * \param[in] solver - string describing the solver
   * \param[in] restol - the target tolerance to solve to
   * \param[in] resinit - the initial residual norm (absolute)
   *
   * \pre the ostream object os should be open
   */
  void writeHeader(const string & solver, const double & restol, const double & resinit);
  
  /*!
   * \brief writes residual convergence data for one iteration to a stream
   * \param[in,out] os - ostream class object for output
   * \param[in] iter - current iteration
   * \param[in] res - the (absolute) residual norm value
   * \param[in] resinit - the initial residual norm
   *
   * \pre the ostream object os should be open
   */
  void writeHistory(const int & iter, const double & res, const double & resinit);
  
public:
  
  /*! \brief Conjugate Gradient method
   * \param[in] b - the right hand size vector
   * \param[in,out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  unsigned long ConjugateGradient(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                  CPreconditioner & precond, double tol,
                                  unsigned long m, bool monitoring);
	
  /*!
   * \brief Flexible Generalized Minimal Residual method
   * \param[in] b - the right hand size vector
   * \param[in,out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  unsigned long FGMRES(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                      CPreconditioner & precond, double tol,
                      unsigned long m, bool monitoring);
	
	/*!
   * \brief Biconjugate Gradient Stabilized Method (BCGSTAB)
   * \param[in] b - the right hand size vector
   * \param[in,out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  unsigned long BCGSTAB(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                        CPreconditioner & precond, double tol,
                        unsigned long m, bool monitoring);
  
};

#include "linear_solvers_structure.inl"
