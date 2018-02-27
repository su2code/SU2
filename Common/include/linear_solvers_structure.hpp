/*!
 * \file linear_solvers_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 * \author J. Hicken, F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
 * \class TCLinSolver
 * \brief Class for solving linear systems using classical and Krylov-subspace iterative methods
 * \author J. Hicken, T. Albring
 * \version 5.0.0 "Raven"
 *
 * The individual solvers could be stand-alone subroutines, but by
 * creating CSysSolve objects we can more easily assign different
 * matrix-vector products and preconditioners to different problems
 * that may arise in a hierarchical solver (i.e. multigrid).
 */
template<class CalcType>
class TCLinSolver {
  
protected:  
  CalcType Residual, Tolerance;
  unsigned long SubSpaceSize, nElem, nElemDomain, BlockSize;
  bool Output;
  
public:
  
  TCLinSolver(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring);
   
  virtual ~TCLinSolver();
  
  virtual unsigned long Solve(const TCSysVector<CalcType>  & b, TCSysVector<CalcType>  & x,
                              TCMatrixVectorProduct<CalcType> & mat_vec, TCPreconditioner<CalcType> & precond);
  
 
protected:
  /*!
   * \brief writes header information for a CSysSolve residual history
   * \param[in] solver - string describing the solver
   * \param[in] restol - the target tolerance to solve to
   * \param[in] resinit - the initial residual norm (absolute)
   *
   */
  void WriteHeader(const string & solver, const CalcType & restol, const CalcType & resinit);
  
  /*!
   * \brief writes residual convergence data for one iteration to a stream
   * \param[in] iter - current iteration
   * \param[in] res - the (absolute) residual norm value
   * \param[in] resinit - the initial residual norm
   *
   */  
  void WriteHistory(const int & iter, const CalcType & res, const CalcType & resinit);
};

template<class CalcType>
class TCLinSolver_CG : public TCLinSolver<CalcType>{

  typedef TCLinSolver<CalcType> TCLinSolverBase;
  
private:

  TCSysVector<CalcType> r, A_p, z, p;
  
  CalcType norm_r, norm0, alpha, beta, r_dot_z;
  
public:
  
   TCLinSolver_CG(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring);
   
   ~TCLinSolver_CG();
  
   unsigned long Solve(const TCSysVector<CalcType>  & b, TCSysVector<CalcType>  & x,
                       TCMatrixVectorProduct<CalcType> & mat_vec, TCPreconditioner<CalcType> & precond);
};

template<class CalcType>
class TCLinSolver_FGMRES : public TCLinSolver<CalcType>{
  
  typedef TCLinSolver<CalcType> TCLinSolverBase;
  
private:
  
  vector<TCSysVector<CalcType> > w, z;
  vector<CalcType> g, sn, cs, y;
  vector<vector<CalcType> > H;
  
  CalcType norm0, beta;
  
public:
  
   TCLinSolver_FGMRES(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring);
  
   ~TCLinSolver_FGMRES();   
   
   unsigned long Solve(const TCSysVector<CalcType>  & b, TCSysVector<CalcType>  & x,
                       TCMatrixVectorProduct<CalcType> & mat_vec, TCPreconditioner<CalcType> & precond);
   
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
   CalcType Sign(const CalcType & x, const CalcType & y) const;   
   
  /*!
  * \brief applys a Givens rotation to a 2-vector
  * \param[in] s - sine of the Givens rotation angle
  * \param[in] c - cosine of the Givens rotation angle
  * \param[in, out] h1 - first element of 2x1 vector being transformed
  * \param[in, out] h2 - second element of 2x1 vector being transformed
  */
  void ApplyGivens(const CalcType & s, const CalcType & c, CalcType & h1, CalcType & h2);
  
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
  void GenerateGivens(CalcType & dx, CalcType & dy, CalcType & s, CalcType & c);
  
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
  void SolveReduced(const int & n, const vector<vector<CalcType> > & Hsbg,
                    const vector<CalcType> & rhs, vector<CalcType> & x);
  
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
  void ModGramSchmidt(int i, vector<vector<CalcType> > & Hsbg, vector<TCSysVector<CalcType> > & w);
};

template<class CalcType>
class TCLinSolver_BCGSTAB : public TCLinSolver<CalcType>{
  
  typedef TCLinSolver<CalcType> TCLinSolverBase;
  
private:

  TCSysVector<CalcType> r, r_0, p, v, s, t, phat, shat, A_x;
  
  CalcType norm_r, norm0, alpha, beta, omega, rho, rho_prime, beta_omega, r_0_v;

public:
  
  TCLinSolver_BCGSTAB(unsigned short blocksize, unsigned long elem, unsigned long elemdomain, CalcType tol, unsigned long m, bool monitoring);
  
  ~TCLinSolver_BCGSTAB();
  
   unsigned long Solve(const TCSysVector<CalcType>  & b, TCSysVector<CalcType>  & x,
                       TCMatrixVectorProduct<CalcType> & mat_vec, TCPreconditioner<CalcType> & precond);
   
   
};

template <class CalcType, class BaseType>
class Convert{
  
  unsigned short BlockSize;
  
  CalcType* BlockLin_CalcType;
  BaseType* BlockLin_BaseType;
  CalcType** Block_CalcType;
  BaseType** Block_BaseType;


public:
  Convert();
  
  ~Convert();
  
  void Initialize(unsigned short MaxSize);
  
  BaseType** ToBaseType(CalcType** block, unsigned short BlockSizeI, unsigned short BlockSizeJ);
  CalcType** ToCalcType(BaseType** block, unsigned short BlockSizeI, unsigned short BlockSizeJ);

  BaseType* ToBaseType(CalcType* block, unsigned short BlockSize);
  CalcType* ToCalcType(BaseType* block, unsigned short BlockSize);
  
  BaseType ToBaseType(CalcType val);
  CalcType ToCalcType(BaseType val);
};

/*!
 * \class TCSysSolve
 * \brief Class for solving linear systems using classical and Krylov-subspace iterative methods
 * \author J. Hicken.
 * \version 5.0.0 "Raven"
 *
 * The individual solvers could be stand-alone subroutines, but by
 * creating CSysSolve objects we can more easily assign different
 * matrix-vector products and preconditioners to different problems
 * that may arise in a hierarchical solver (i.e. multigrid).
 */
template<class CalcType, class BaseType>
class TCSysSolve {

private:
  
  unsigned long nPoint, nPointDomain;
  unsigned short BlockSize;
  
  TCLinSolver<CalcType>* LinSolver;
  TCLinSolver<CalcType>* LinSolver_b;
  
  TCPreconditioner<CalcType>*      Precond;
  TCMatrixVectorProduct<CalcType>* MatVec;
  
  TCPreconditioner<CalcType>*      Precond_b;
  TCMatrixVectorProduct<CalcType>* MatVec_b;
  
  TCSysMatrix<CalcType> Matrix;
  TCSysVector<BaseType> LinSysSol,
                        LinSysRes;
  
  TCSysVector<CalcType> LinSysRes_b,
                        LinSysSol_b;

  TCSysVector<CalcType> LinSysRes_calc,
                        LinSysSol_calc;
  
  Convert<CalcType, BaseType> convert;
  
  unsigned short kind_prec, kind_prec_b;
public:

  TCSysSolve();
  
  ~TCSysSolve();
  
  void Initialize_System(unsigned short BlockSize,  bool edgeconnect, CGeometry *geometry, CConfig *config);
  
  void Initialize_System_Adjoint(unsigned short BlockSize, CGeometry *geometry, CConfig* config);
  
  void Initialize_Linear_Solver(unsigned short blocksize, unsigned short kind_solver, unsigned short kind_preconditioner, unsigned long max_iter, BaseType error, CGeometry *geometry, CConfig *config);
  
  void Initialize_Linear_Solver_Adjoint(unsigned short blocksize, unsigned short kind_solver, unsigned short kind_preconditioner, unsigned long max_iter, BaseType error, CGeometry *geometry, CConfig *config);
  
  void SetValZero_Matrix();
  void SetValZero_Rhs();
  void SetValZero_Sol();
  void DeleteValsRowi(unsigned long i);
  
  void SetVal2Diag_Matrix(unsigned long iPoint, BaseType val);
  void AddVal2Diag_Matrix(unsigned long iPoint, BaseType val);
  
  void AddBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType** block);
  void SubtractBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType** block);
  
  BaseType* GetBlock_Matrix(unsigned long iPoint, unsigned long jPoint);
  BaseType GetBlock_Matrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar);
  void SetBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType* block);
  void SetBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType** block);
  
  void MatrixVectorProduct();    
  
  void Build_Preconditioner(unsigned short kind_prec, bool transpose);
  
  void Solve_System(TCSysVector<BaseType> &Rhs, TCSysVector<BaseType> &Sol);
    
};



//typedef TCSysSolve<su2double> CSysSolve;

template class TCSysSolve<passivedouble, su2double>;
template class TCLinSolver<passivedouble>;
template class TCLinSolver_BCGSTAB<passivedouble>;
template class TCLinSolver_CG<passivedouble>;
template class TCLinSolver_FGMRES<passivedouble>;
#ifdef CODI_REVERSE_TYPE
template class TCSysSolve<su2double, su2double>;
#endif
#include "linear_solvers_structure.inl"
