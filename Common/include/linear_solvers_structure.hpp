/*!
 * \file linear_solvers_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 * \author J. Hicken, F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "./mpi_structure.hpp"

#include <climits>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <string>

#include "option_structure.hpp"
#include "vector_structure.hpp"
#include "matrix_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"

using namespace std;

/*!
 * \class CSysSolve
 * \brief Class for solving linear systems using classical and Krylov-subspace iterative methods
 * \author J. Hicken.
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
  su2double Sign(const su2double & x, const su2double & y) const;
  
  /*!
   * \brief applys a Givens rotation to a 2-vector
   * \param[in] s - sine of the Givens rotation angle
   * \param[in] c - cosine of the Givens rotation angle
   * \param[in, out] h1 - first element of 2x1 vector being transformed
   * \param[in, out] h2 - second element of 2x1 vector being transformed
   */
  void ApplyGivens(const su2double & s, const su2double & c, su2double & h1, su2double & h2);
  
  /*!
   * \brief generates the Givens rotation matrix for a given 2-vector
   * \param[in, out] dx - element of 2x1 vector being transformed
   * \param[in, out] dy - element of 2x1 vector being set to zero
   * \param[in, out] s - sine of the Givens rotation angle
   * \param[in, out] c - cosine of the Givens rotation angle
   *
   * Based on givens() of SPARSKIT, which is based on p.202 of
   * "Matrix Computations" by Golub and van Loan.
   */
  void GenerateGivens(su2double & dx, su2double & dy, su2double & s, su2double & c);
  
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
  void SolveReduced(const int & n, const vector<vector<su2double> > & Hsbg,
                    const vector<su2double> & rhs, vector<su2double> & x);
  
  /*!
   * \brief Modified Gram-Schmidt orthogonalization
   * \author Based on Kesheng John Wu's mgsro subroutine in Saad's SPARSKIT
   *
   * \tparam Vec - a generic vector class
   * \param[in] i - index indicating which vector in w is being orthogonalized
   * \param[in, out] Hsbg - the upper Hessenberg begin updated
   * \param[in, out] w - the (i+1)th vector of w is orthogonalized against the
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
  void ModGramSchmidt(int i, vector<vector<su2double> > & Hsbg, vector<CSysVector> & w);
  
  /*!
   * \brief writes header information for a CSysSolve residual history
   * \param[in] solver - string describing the solver
   * \param[in] restol - the target tolerance to solve to
   * \param[in] resinit - the initial residual norm (absolute)
   *
   * \pre the ostream object os should be open
   */
  void WriteHeader(const string & solver, const su2double & restol, const su2double & resinit);
  
  /*!
   * \brief writes residual convergence data for one iteration to a stream
   * \param[in] iter - current iteration
   * \param[in] res - the (absolute) residual norm value
   * \param[in] resinit - the initial residual norm
   *
   * \pre the ostream object os should be open
   */
  void WriteHistory(const int & iter, const su2double & res, const su2double & resinit);
  
public:
  
  /*! \brief Conjugate Gradient method
   * \param[in] b - the right hand size vector
   * \param[in, out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  unsigned long CG_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                  CPreconditioner & precond, su2double tol,
                                  unsigned long m, su2double *residual, bool monitoring);
	
  /*!
   * \brief Flexible Generalized Minimal Residual method
   * \param[in] b - the right hand size vector
   * \param[in, out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  unsigned long FGMRES_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                      CPreconditioner & precond, su2double tol,
                      unsigned long m, su2double *residual, bool monitoring);
	
	/*!
   * \brief Biconjugate Gradient Stabilized Method (BCGSTAB)
   * \param[in] b - the right hand size vector
   * \param[in, out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  unsigned long BCGSTAB_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                        CPreconditioner & precond, su2double tol,
                        unsigned long m, su2double *residual, bool monitoring);
  
  /*!
   * \brief Solve the linear system using a Krylov subspace method
   * \param[in] Jacobian - Jacobian Matrix for the linear system
   * \param[in] LinSysRes - Linear system residual
   * \param[in] LinSysSol - Linear system solution
   * \param[in] geometry -  Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  unsigned long Solve(CSysMatrix & Jacobian, CSysVector & LinSysRes, CSysVector & LinSysSol, CGeometry *geometry, CConfig *config);
  

  /*!
   * \brief Prepare the linear solve during the reverse interpretation of the AD tape.
   * \param[in] Jacobian - Jacobian Matrix for the linear system
   * \param[in] LinSysRes - Linear system residual
   * \param[in] LinSysSol - Linear system solution
   * \param[in] geometry -  Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetExternalSolve(CSysMatrix & Jacobian, CSysVector & LinSysRes, CSysVector & LinSysSol, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Prepare the linear solve during the reverse interpretation of the AD tape.
   * \param[in] Jacobian - Jacobian Matrix for the linear system
   * \param[in] LinSysRes - Linear system residual
   * \param[in] LinSysSol - Linear system solution
   * \param[in] geometry -  Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetExternalSolve_Mesh(CSysMatrix & Jacobian, CSysVector & LinSysRes, CSysVector & LinSysSol, CGeometry *geometry, CConfig *config);


};

#include "linear_solvers_structure.inl"
