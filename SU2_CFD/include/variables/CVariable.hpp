/*!
 * \file CVariable.hpp
 * \brief Declaration and inlines of the parent class for defining problem
          variables, function definitions in file <i>CVariable.cpp</i>.
          All variables are children of at least this class.
 * \author F. Palacios, T. Economon
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

#include "../../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../../Common/include/config_structure.hpp"
#include "../fluid_model.hpp"
#include "../../../Common/include/toolboxes/C2DContainer.hpp"


using namespace std;

/*!
 * \class CVariable
 * \brief Main class for defining the variables.
 * \author F. Palacios
 */
class CVariable {
protected:
  using Idx_t = unsigned long; // ToDo: replace all before merging with develop
  using VectorType = C2DContainer<Idx_t, su2double, StorageType::ColumnMajor, 64, DynamicSize, 1>;
  using MatrixType = C2DContainer<Idx_t, su2double, StorageType::RowMajor,    64, DynamicSize, DynamicSize>;

  /*--- This contrived container is used to store matrices in a contiguous manner but still present the
   "su2double**" interface to the outside world, it will be replaced by something more efficient. ---*/
  struct VectorOfMatrix {
    su2activevector storage;
    su2matrix<su2double*> interface;
    Idx_t M, N;

    void resize(Idx_t length, Idx_t rows, Idx_t cols, su2double value) {
      M = rows;
      N = cols;
      storage.resize(length*rows*cols) = value;
      interface.resize(length,rows);

      for(Idx_t i=0; i<length; ++i)
        for(Idx_t j=0; j<rows; ++j)
          interface(i,j) = &(*this)(i,j,0);
    }

    su2double& operator() (Idx_t i, Idx_t j, Idx_t k) { return storage(i*M*N + j*N + k); }
    const su2double& operator() (Idx_t i, Idx_t j, Idx_t k) const { return storage(i*M*N + j*N + k); }

    su2double** operator[] (Idx_t i) { return interface[i]; }
  };

  MatrixType Solution;       /*!< \brief Solution of the problem. */
  MatrixType Solution_Old;   /*!< \brief Old solution of the problem R-K. */

  MatrixType External;       /*!< \brief External (outer) contribution in discrete adjoint multizone problems. */
  MatrixType External_Old;   /*!< \brief Old external (outer) contribution in discrete adjoint multizone problems. */

  su2vector<bool> Non_Physical;  /*!< \brief Non-physical points in the solution (force first order). */

  MatrixType Solution_time_n;    /*!< \brief Solution of the problem at time n for dual-time stepping technique. */
  MatrixType Solution_time_n1;   /*!< \brief Solution of the problem at time n-1 for dual-time stepping technique. */
  VectorType Delta_Time;         /*!< \brief Time step. */

  VectorOfMatrix Gradient;  /*!< \brief Gradient of the solution of the problem. */
  VectorOfMatrix Rmatrix;   /*!< \brief Geometry-based matrix for weighted least squares gradient calculations. */

  MatrixType Limiter;        /*!< \brief Limiter of the solution of the problem. */
  MatrixType Solution_Max;   /*!< \brief Max solution for limiter computation. */
  MatrixType Solution_Min;   /*!< \brief Min solution for limiter computation. */

  VectorType AuxVar;       /*!< \brief Auxiliar variable for gradient computation. */
  MatrixType Grad_AuxVar;  /*!< \brief Gradient of the auxiliar variable. */

  VectorType Max_Lambda_Inv;   /*!< \brief Maximun inviscid eingenvalue. */
  VectorType Max_Lambda_Visc;  /*!< \brief Maximun viscous eingenvalue. */
  VectorType Lambda;           /*!< \brief Value of the eingenvalue. */

  VectorType Sensor;               /*!< \brief Pressure sensor for high order central scheme and Roe dissipation. */
  MatrixType Undivided_Laplacian;  /*!< \brief Undivided laplacian of the solution. */

  MatrixType Res_TruncError;  /*!< \brief Truncation error for multigrid cycle. */
  MatrixType Residual_Old;    /*!< \brief Auxiliar structure for residual smoothing. */
  MatrixType Residual_Sum;    /*!< \brief Auxiliar structure for residual smoothing. */

  MatrixType Solution_Adj_Old;   /*!< \brief Solution of the problem in the previous AD-BGS iteration. */

  MatrixType Solution_BGS_k;     /*!< \brief Old solution container for BGS iterations. */

  su2matrix<int> Input_AdjIndices;    /*!< \brief Indices of Solution variables in the adjoint vector. */
  su2matrix<int> Output_AdjIndices;   /*!< \brief Indices of Solution variables in the adjoint vector after having been updated. */

  Idx_t nPoint = {0};  /*!< \brief Number of points in the domain. */
  Idx_t nDim = {0};      /*!< \brief Number of dimension of the problem. */
  Idx_t nVar = {0};        /*!< \brief Number of variables of the problem. */
  Idx_t nPrimVar = {0};      /*!< \brief Number of primitive variables. */
  Idx_t nPrimVarGrad = {0};    /*!< \brief Number of primitives for which a gradient is computed. */
  Idx_t nSecondaryVar = {0};     /*!< \brief Number of secondary variables. */
  Idx_t nSecondaryVarGrad = {0};   /*!< \brief Number of secondaries for which a gradient is computed. */
  
public:

  /*--- Disable default construction copy and assignment. ---*/
  CVariable() = delete;
  CVariable(const CVariable&) = delete;
  CVariable(CVariable&&) = delete;
  CVariable& operator= (const CVariable&) = delete;
  CVariable& operator= (CVariable&&) = delete;

  /*!
   * \overload
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CVariable(Idx_t npoint, Idx_t nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CVariable(Idx_t npoint, Idx_t ndim, Idx_t nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CVariable() = default;

  /*!
   * \brief Set the value of the solution, all variables.
   * \param[in] iPoint - Point index.
   * \param[in] solution - Solution of the problem.
   */
  inline void SetSolution(Idx_t iPoint, const su2double *solution) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++) Solution(iPoint,iVar) = solution[iVar];
  }

  /*!
   * \brief Set the value of the solution, one variable.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution for the index <i>iVar</i>.
   */
  inline void SetSolution(Idx_t iPoint, Idx_t iVar, su2double solution) { Solution(iPoint,iVar) = solution; }

  /*!
   * \brief Add the value of the solution vector to the previous solution (incremental approach).
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution for the index <i>iVar</i>.
   */
  inline void Add_DeltaSolution(Idx_t iPoint, Idx_t iVar, su2double solution) { Solution(iPoint,iVar) += solution; }

  /*!
   * \brief Set the value of the non-physical point.
   * \param[in] iPoint - Point index.
   * \param[in] value - identification of the non-physical point.
   */
  inline void SetNon_Physical(Idx_t iPoint, bool value) { Non_Physical(iPoint) = !value; }

  /*!
   * \brief Get the value of the non-physical point.
   * \param[in] iPoint - Point index.
   * \return Value of the Non-physical point.
   */
  inline su2double GetNon_Physical(Idx_t iPoint) const { return su2double(Non_Physical(iPoint)); }

  /*!
   * \brief Get the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution(Idx_t iPoint, Idx_t iVar) const { return Solution(iPoint,iVar); }

  /*!
   * \brief Get the old solution of the problem (Runge-Kutta method)
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Old(Idx_t iPoint, Idx_t iVar) const { return Solution_Old(iPoint,iVar); }

  /*!
   * \brief Get the old solution of the discrete adjoint problem (for multiphysics subiterations)
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Old_Adj(Idx_t iPoint, Idx_t iVar) const { return Solution_Adj_Old(iPoint,iVar); }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] iPoint - Point index.
   * \param[in] solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Old(Idx_t iPoint, const su2double *solution_old) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Solution_Old(iPoint,iVar) = solution_old[iVar];
  }

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] iPoint - Point index.
   * \param[in] solution_old - Value of the old solution for the index <i>iVar</i>.
   */
  inline void SetSolution_Old(Idx_t iPoint, Idx_t iVar, su2double solution_old) { Solution_Old(iPoint,iVar) = solution_old; }

  /*!
   * \brief Set old variables to the value of the current variables.
   */
  void Set_OldSolution();

  /*!
   * \brief Set variables to the value of the old variables.
   */
  void Set_Solution();

  /*!
   * \brief Set the variable solution at time n.
   */
  void Set_Solution_time_n();

  /*!
   * \brief Set the variable solution at time n-1.
   */
  void Set_Solution_time_n1();

  /*!
   * \brief Set the variable solution at time n.
   * \param[in] iPoint - Point index.
   */
  inline void Set_Solution_time_n(Idx_t iPoint, const su2double* val_sol) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++) Solution_time_n(iPoint,iVar) = val_sol[iVar];
  }

  /*!
   * \brief Set the variable solution at time n-1.
   * \param[in] iPoint - Point index.
   */
  inline void Set_Solution_time_n1(Idx_t iPoint, const su2double* val_sol) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Solution_time_n1(iPoint,iVar) = val_sol[iVar];
  }

  /*!
   * \brief Set the variable solution at time n.
   */
  inline void Set_Solution_time_n(Idx_t iPoint, Idx_t iVar, su2double val_sol) {
    Solution_time_n(iPoint,iVar) = val_sol;
  }

  /*!
   * \brief Set the variable solution at time n-1.
   */
  inline void Set_Solution_time_n1(Idx_t iPoint, Idx_t iVar, su2double val_sol) {
    Solution_time_n1(iPoint,iVar) = val_sol;
  }

  /*!
   * \brief Set to zero the velocity components of the solution.
   * \param[in] iPoint - Point index.
   */
  inline void SetVelSolutionZero(Idx_t iPoint) {
    for (Idx_t iDim = 0; iDim < nDim; iDim++) Solution(iPoint,iDim+1) = 0.0;
  }

  /*!
   * \brief Specify a vector to set the velocity components of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionVector(Idx_t iPoint, const su2double *val_vector) {
    for (Idx_t iDim = 0; iDim < nDim; iDim++) Solution(iPoint, iDim+1) = val_vector[iDim];
  }

  /*!
   * \brief Set to zero velocity components of the solution.
   */
  inline void SetVelSolutionOldZero(Idx_t iPoint) {
    for (Idx_t iDim = 0; iDim < nDim; iDim++) Solution_Old(iPoint, iDim+1) = 0.0;
  }

  /*!
   * \brief Add a value to the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Number of the variable.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline void AddSolution(Idx_t iPoint, Idx_t iVar, su2double solution) {
    Solution(iPoint, iVar) = Solution_Old(iPoint, iVar) + solution;
  }

  /*!
   * \brief Add a value to the solution.
   * \param[in] iPoint - Point index.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline void AddSolution(Idx_t iPoint, const su2double *solution) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++) Solution(iPoint, iVar) += solution[iVar];
  }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_New(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual su2double GetRoe_Dissipation(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetRoe_Dissipation(Idx_t iPoint, su2double val_dissipation) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetRoe_Dissipation_FD(Idx_t iPoint, su2double val_wall_dist) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_delta - A scalar measure of the grid size
   * \param[in] val_const_DES - The DES constant (C_DES)
   */
  inline virtual void SetRoe_Dissipation_NTS(Idx_t iPoint, su2double val_delta, su2double val_const_DES) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual su2double GetDES_LengthScale(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetDES_LengthScale(Idx_t iPoint, su2double val_des_lengthscale) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  virtual void SetSolution_New() {}

  /*!
   * \brief Set external contributions to zero.
   */
  void SetExternalZero();

  /*!
   * \brief Set old External to the value of the current variables.
   */
  void Set_OldExternal();

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Number of the variable.
   * \param[in] val_solution - Value that we want to add to the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Number of the variable.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline virtual void AddSolution_New(Idx_t iPoint, Idx_t iVar, su2double solution) {}

  /*!
   * \brief Add a value to the External vector.
   * \param[in] iPoint - Point index.
   * \param[in] val_sol - vector that has to be added component-wise
   */
  inline void Add_External(Idx_t iPoint, const su2double* val_sol) {
    for(Idx_t iVar = 0; iVar < nVar; iVar++) External(iPoint,iVar) += val_sol[iVar];
  }

  /*!
   * \brief Add a valaue to the old External vector.
   * \param[in] iPoint - Point index.
   * \param[in] val_solution - Value that we want to add to the solution.
   */
  inline void Add_ExternalOld(Idx_t iPoint, const su2double *val_sol) {
    for(Idx_t iVar = 0; iVar < nVar; iVar++) External_Old(iPoint,iVar) += val_sol[iVar];
  }

  /*!
   * \brief Add a value to the solution, clipping the values.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution change.
   * \param[in] lowerlimit - Lower value.
   * \param[in] upperlimit - Upper value.
   */
  inline void AddClippedSolution(Idx_t iPoint, Idx_t iVar, su2double solution,
                                 su2double lowerlimit, su2double upperlimit) {

    su2double val_new = Solution_Old(iPoint, iVar) + solution;
    Solution(iPoint,iVar) = min(max(val_new, lowerlimit), upperlimit);
  }

  /*!
   * \brief Update the variables using a conservative format.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution change.
   * \param[in] val_density - Value of the density.
   * \param[in] val_density_old - Value of the old density.
   * \param[in] lowerlimit - Lower value.
   * \param[in] upperlimit - Upper value.
   */
  inline void AddConservativeSolution(Idx_t iPoint, Idx_t iVar, su2double solution,
                                      su2double val_density, su2double val_density_old,
                                      su2double lowerlimit, su2double upperlimit) {

    su2double val_new = (Solution_Old(iPoint,iVar)*val_density_old + solution)/val_density;
    Solution(iPoint,iVar) = min(max(val_new, lowerlimit), upperlimit);
  }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline su2double *GetSolution(Idx_t iPoint) { return Solution[iPoint]; }

  /*!
   * \brief Get the old solution of the problem (Runge-Kutta method)
   * \return Pointer to the old solution vector.
   */
  inline su2double *GetSolution_Old(Idx_t iPoint) { return Solution_Old[iPoint]; }

  /*!
   * \brief Get the external contributions of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the External row for iPoint.
   */
  inline const su2double *Get_External(Idx_t iPoint) const { return External[iPoint]; }
  
  /*!
   * \brief Get the old external contributions of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the External_Old row for iPoint.
   */
  inline const su2double *Get_ExternalOld(Idx_t iPoint) const { return External_Old[iPoint]; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_time_n(Idx_t iPoint) { return Solution_time_n[iPoint]; }

  /*!
   * \brief Get the solution at time n-1.
   * \return Pointer to the solution (at time n-1) vector.
   */
  inline su2double *GetSolution_time_n1(Idx_t iPoint) { return Solution_time_n1[iPoint]; }

  /*!
   * \brief Set the value of the old residual.
   * \param[in] iPoint - Point index.
   * \param[in] val_residual_old - Pointer to the residual vector.
   */
  inline void SetResidual_Old(Idx_t iPoint, const su2double *val_residual_old) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Residual_Old(iPoint,iVar) = val_residual_old[iVar];
  }

  /*!
   * \brief Add a value to the summed residual vector.
   * \param[in] iPoint - Point index.
   * \param[in] val_residual - Pointer to the residual vector.
   */
  inline void AddResidual_Sum(Idx_t iPoint, const su2double *val_residual) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Residual_Sum(iPoint,iVar) += val_residual[iVar];
  }

  /*!
   * \brief Set summed residual vector to zero value.
   */
  void SetResidualSumZero();

  /*!
   * \brief Set the velocity of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetVel_ResTruncError_Zero(Idx_t iPoint, Idx_t iSpecies) {}

  /*!
   * \brief Get the value of the summed residual.
   * \param[in] iPoint - Point index.
   * \return Pointer to the summed residual.
   */
  inline su2double *GetResidual_Sum(Idx_t iPoint) { return Residual_Sum[iPoint]; }

  /*!
   * \brief Get the value of the old residual.
   * \param[in] iPoint - Point index.
   * \return Pointer to the old residual.
   */
  inline su2double *GetResidual_Old(Idx_t iPoint) { return Residual_Old[iPoint]; }

  /*!
   * \brief Get the value of the summed residual.
   * \param[in] iPoint - Point index.
   * \param[in] val_residual - Pointer to the summed residual.
   */
  inline void GetResidual_Sum(Idx_t iPoint, su2double *val_residual) const {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] = Residual_Sum(iPoint,iVar);
  }

  /*!
   * \brief Set auxiliar variables, we are looking for the gradient of that variable.
   * \param[in] iPoint - Point index.
   * \param[in] val_auxvar - Value of the auxiliar variable.
   */
  inline void SetAuxVar(Idx_t iPoint, su2double val_auxvar) { AuxVar(iPoint) = val_auxvar; }

  /*!
   * \brief Get the value of the auxiliary variable.
   * \param[in] iPoint - Point index.
   * \return Value of the auxiliary variable.
   */
  inline su2double GetAuxVar(Idx_t iPoint) const { return AuxVar(iPoint); }

  /*!
   * \brief Set the auxiliary variable gradient to zero value.
   */
  void SetAuxVarGradientZero();

  /*!
   * \brief Set the value of the auxiliary variable gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] val_gradient - Value of the gradient for the index <i>iDim</i>.
   */
  inline void SetAuxVarGradient(Idx_t iPoint, Idx_t iDim, su2double val_gradient) { Grad_AuxVar(iPoint,iDim) = val_gradient; }

  /*!
   * \brief Add a value to the auxiliary variable gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] val_value - Value of the gradient to be added for the index <i>iDim</i>.
   */
  inline void AddAuxVarGradient(Idx_t iPoint, Idx_t iDim, su2double val_value) { Grad_AuxVar(iPoint,iDim) += val_value;}

  /*!
   * \brief Subtract a value to the auxiliary variable gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] val_value - Value of the gradient to be subtracted for the index <i>iDim</i>.
   */
  inline void SubtractAuxVarGradient(Idx_t iPoint, Idx_t iDim, su2double val_value) { Grad_AuxVar(iPoint,iDim) -= val_value; }

  /*!
   * \brief Get the gradient of the auxiliary variable.
   * \param[in] iPoint - Point index.
   * \return Value of the gradient of the auxiliary variable.
   */
  inline su2double *GetAuxVarGradient(Idx_t iPoint) { return Grad_AuxVar[iPoint]; }

  /*!
   * \brief Get the gradient of the auxiliary variable.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the gradient of the auxiliary variable for the dimension <i>iDim</i>.
   */
  inline su2double GetAuxVarGradient(Idx_t iPoint, Idx_t iDim) const { return Grad_AuxVar(iPoint,iDim); }

  /*!
   * \brief Add a value to the truncation error.
   * \param[in] iPoint - Point index.
   * \param[in] val_truncation_error - Value that we want to add to the truncation error.
   */
  inline void AddRes_TruncError(Idx_t iPoint, const su2double *val_truncation_error) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Res_TruncError(iPoint, iVar) += val_truncation_error[iVar];
  }

  /*!
   * \brief Subtract a value to the truncation error.
   * \param[in] iPoint - Point index.
   * \param[in] val_truncation_error - Value that we want to subtract to the truncation error.
   */
  inline void SubtractRes_TruncError(Idx_t iPoint, const su2double *val_truncation_error) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Res_TruncError(iPoint, iVar) -= val_truncation_error[iVar];
  }

  /*!
   * \brief Set the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetRes_TruncErrorZero(Idx_t iPoint) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++) Res_TruncError(iPoint, iVar) = 0.0;
  }

  /*!
   * \brief Set the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVal_ResTruncError_Zero(Idx_t iPoint, Idx_t iVar) {Res_TruncError(iPoint, iVar) = 0.0;}

  /*!
   * \brief Set the velocity of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVel_ResTruncError_Zero(Idx_t iPoint) {
    for (Idx_t iDim = 0; iDim < nDim; iDim++) Res_TruncError(iPoint,iDim+1) = 0.0;
  }

  /*!
   * \brief Set the velocity of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetEnergy_ResTruncError_Zero(Idx_t iPoint) { Res_TruncError(iPoint,nDim+1) = 0.0;}

  /*!
   * \brief Get the truncation error.
   * \param[in] iPoint - Point index.
   * \return Pointer to the truncation error.
   */
  inline su2double *GetResTruncError(Idx_t iPoint) { return Res_TruncError[iPoint]; }

  /*!
   * \brief Get the truncation error.
   * \param[in] iPoint - Point index.
   * \param[in] val_trunc_error - Pointer to the truncation error.
   */
  inline void GetResTruncError(Idx_t iPoint, su2double *val_trunc_error) const {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      val_trunc_error[iVar] = Res_TruncError(iPoint, iVar);
  }

  /*!
   * \brief Set the gradient of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] gradient - Gradient of the solution.
   */
  inline void SetGradient(Idx_t iPoint, su2double **gradient) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      for (Idx_t iDim = 0; iDim < nDim; iDim++)
        Gradient(iPoint,iVar,iDim) = gradient[iVar][iDim];
  }

  /*!
   * \overload
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value of the gradient.
   */
  inline void SetGradient(Idx_t iPoint, Idx_t iVar, Idx_t iDim, su2double value) { Gradient(iPoint,iVar,iDim) = value; }

  /*!
   * \brief Set to zero the gradient of the solution.
   */
  void SetGradientZero();

  /*!
   * \brief Add <i>value</i> to the solution gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value to add to the solution gradient.
   */
  inline void AddGradient(Idx_t iPoint, Idx_t iVar, Idx_t iDim, su2double value) { Gradient(iPoint,iVar,iDim) += value; }

  /*!
   * \brief Subtract <i>value</i> to the solution gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value to subtract to the solution gradient.
   */
  inline void SubtractGradient(Idx_t iPoint, Idx_t iVar, Idx_t iDim, su2double value) { Gradient(iPoint,iVar,iDim) -= value; }

  /*!
   * \brief Get the value of the solution gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the gradient solution.
   */
  inline su2double **GetGradient(Idx_t iPoint) { return Gradient[iPoint]; }

  /*!
   * \brief Get the value of the solution gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the solution gradient.
   */
  inline su2double GetGradient(Idx_t iPoint, Idx_t iVar, Idx_t iDim) const { return Gradient(iPoint,iVar,iDim); }

  /*!
   * \brief Set the value of an entry in the Rmatrix for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] jDim - Index of the dimension.
   * \param[in] value - Value of the Rmatrix entry.
   */
  inline void SetRmatrix(Idx_t iPoint, Idx_t iDim, Idx_t jDim, su2double value) { Rmatrix(iPoint,iDim,jDim) = value; }

  /*!
   * \brief Set to zero the Rmatrix for least squares gradient calculations.
   */
  void SetRmatrixZero();

  /*!
   * \brief Add <i>value</i> to the Rmatrix for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] jDim - Index of the dimension.
   * \param[in] value - Value of the Rmatrix entry.
   */
  inline void AddRmatrix(Idx_t iPoint, Idx_t iDim, Idx_t jDim, su2double value) { Rmatrix(iPoint,iDim,jDim) += value; }

  /*!
   * \brief Get the value of the Rmatrix entry for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] jDim - Index of the dimension.
   * \return Value of the Rmatrix entry.
   */
  inline su2double GetRmatrix(Idx_t iPoint, Idx_t iDim, Idx_t jDim) const { return Rmatrix(iPoint,iDim,jDim); }

  /*!
   * \brief Get the value of the Rmatrix entry for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \return Value of the Rmatrix entry.
   */
  inline su2double **GetRmatrix(Idx_t iPoint) { return Rmatrix[iPoint]; }

  /*!
   * \brief Set the value of the limiter.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_limiter - Value of the limiter for the index <i>iVar</i>.
   */
  inline void SetLimiter(Idx_t iPoint, Idx_t iVar, su2double val_limiter) { Limiter(iPoint,iVar) = val_limiter; }

  /*!
   * \brief Set the value of the limiter.
   * \param[in] iPoint - Point index.
   * \param[in] val_species - Index of the species .
   * \param[in] iVar - Index of the variable.
   * \param[in] val_limiter - Value of the limiter for the index <i>iVar</i>.
   */
  inline virtual void SetLimiterPrimitive(Idx_t iPoint, Idx_t val_species, Idx_t iVar, su2double val_limiter) {}

  /*!
   * \brief Set the value of the limiter.
   * \param[in] iPoint - Point index.
   * \param[in] val_species - Index of the species .
   * \param[in] iVar - Index of the variable.
   */
  inline virtual su2double GetLimiterPrimitive(Idx_t iPoint, Idx_t val_species, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Set the value of the max solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the max solution for the index <i>iVar</i>.
   */
  inline void SetSolution_Max(Idx_t iPoint, Idx_t iVar, su2double solution) { Solution_Max(iPoint,iVar) = solution; }

  /*!
   * \brief Set the value of the min solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the min solution for the index <i>iVar</i>.
   */
  inline void SetSolution_Min(Idx_t iPoint, Idx_t iVar, su2double solution) { Solution_Min(iPoint,iVar) = solution; }

  /*!
   * \brief Get the value of the slope limiter.
   * \param[in] iPoint - Point index.
   * \return Pointer to the limiters vector.
   */
  inline su2double *GetLimiter(Idx_t iPoint) { return Limiter[iPoint]; }

  /*!
   * \brief Get the value of the slope limiter.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the limiter vector for the variable <i>iVar</i>.
   */
  inline su2double GetLimiter(Idx_t iPoint, Idx_t iVar) const { return Limiter(iPoint,iVar); }

  /*!
   * \brief Get the value of the min solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the min solution for the variable <i>iVar</i>.
   */
  inline su2double GetSolution_Max(Idx_t iPoint, Idx_t iVar) const { return Solution_Max(iPoint,iVar); }

  /*!
   * \brief Set the value of the preconditioner Beta.
   * \param[in] val_Beta - Value of the low Mach preconditioner variable Beta
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the min solution for the variable <i>iVar</i>.
   */
  inline su2double GetSolution_Min(Idx_t iPoint, Idx_t iVar) const { return Solution_Min(iPoint,iVar); }

  /*!
   * \brief Get the value of the wind gust
   * \param[in] iPoint - Point index.
   * \return Value of the wind gust
   */
  inline virtual su2double* GetWindGust(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief Set the value of the wind gust
   * \param[in] iPoint - Point index.
   * \param[in] val_WindGust - Value of the wind gust
   */
  inline virtual void SetWindGust(Idx_t iPoint, const su2double* val_WindGust) {}

  /*!
   * \brief Get the value of the derivatives of the wind gust
   * \return Value of the derivatives of the wind gust
   */
  inline virtual su2double* GetWindGustDer(Idx_t iPoint) { return nullptr;}

  /*!
   * \brief Set the value of the derivatives of the wind gust
   * \param[in] iPoint - Point index.
   * \param[in] val_WindGust - Value of the derivatives of the wind gust
   */
  inline virtual void SetWindGustDer(Idx_t iPoint, const su2double* val_WindGust) {}

  /*!
   * \brief Set the value of the time step.
   * \param[in] iPoint - Point index.
   * \param[in] val_delta_time - Value of the time step.
   */
  inline void SetDelta_Time(Idx_t iPoint, su2double val_delta_time) { Delta_Time(iPoint) = val_delta_time; }

  /*!
   * \brief Set the value of the time step.
   * \param[in] iPoint - Point index.
   * \param[in] val_delta_time - Value of the time step.
   * \param[in] iSpecies - Index of the Species.
   */
  inline virtual void SetDelta_Time(Idx_t iPoint, su2double val_delta_time, Idx_t iSpecies) {}

  /*!
   * \brief Get the value of the time step.
   * \param[in] iPoint - Point index.
   * \return Value of the time step.
   */
  inline su2double GetDelta_Time(Idx_t iPoint) const {return Delta_Time(iPoint); }

  /*!
   * \brief Get the value of the time step.
   * \param[in] iPoint - Point index.
   * \param[in] iSpecies - Index of the Species
   * \return Value of the time step.
   */
  inline virtual su2double GetDelta_Time(Idx_t iPoint, Idx_t iSpecies) { return 0.0; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline void SetMax_Lambda_Inv(Idx_t iPoint, su2double val_max_lambda) { Max_Lambda_Inv(iPoint) = val_max_lambda; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_species - Value of the species index to set the maximum eigenvalue.
   */
  inline virtual void SetMax_Lambda_Inv(Idx_t iPoint, su2double val_max_lambda, Idx_t val_species) {}

  /*!
   * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline void SetMax_Lambda_Visc(Idx_t iPoint, su2double val_max_lambda) { Max_Lambda_Visc(iPoint) = val_max_lambda; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_species - Index of the species to set the maximum eigenvalue of the viscous terms.
   */
  inline virtual void SetMax_Lambda_Visc(Idx_t iPoint, su2double val_max_lambda, Idx_t val_species) {}

  /*!
   * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline void AddMax_Lambda_Inv(Idx_t iPoint, su2double val_max_lambda) { Max_Lambda_Inv(iPoint) += val_max_lambda; }

  /*!
   * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline void AddMax_Lambda_Visc(Idx_t iPoint, su2double val_max_lambda) { Max_Lambda_Visc(iPoint) += val_max_lambda; }

  /*!
   * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline su2double GetMax_Lambda_Inv(Idx_t iPoint) const { return Max_Lambda_Inv(iPoint); }

  /*!
   * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline su2double GetMax_Lambda_Visc(Idx_t iPoint) const { return Max_Lambda_Visc(iPoint); }

  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline void SetLambda(Idx_t iPoint, su2double val_lambda) { Lambda(iPoint) = val_lambda; }

  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_lambda - Value of the spectral radius.
   * \param[in] val_iSpecies -Index of species
   */
  inline virtual void SetLambda(Idx_t iPoint, su2double val_lambda, Idx_t val_iSpecies) {}

  /*!
   * \brief Add the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline void AddLambda(Idx_t iPoint, su2double val_lambda) { Lambda(iPoint) += val_lambda; }

  /*!
   * \brief Add the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_iSpecies -Index of species
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline virtual void AddLambda(Idx_t iPoint, su2double val_lambda, Idx_t val_iSpecies) {}

  /*!
   * \brief Get the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \return Value of the spectral radius.
   */
  inline su2double GetLambda(Idx_t iPoint) const { return Lambda(iPoint); }

  /*!
   * \brief Get the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_iSpecies -Index of species
   * \return Value of the spectral radius.
   */
  inline virtual su2double GetLambda(Idx_t iPoint, Idx_t val_iSpecies) { return 0.0; }

  /*!
   * \brief Set pressure sensor.
   * \param[in] iPoint - Point index.
   * \param[in] val_sensor - Value of the pressure sensor.
   */
  inline void SetSensor(Idx_t iPoint, su2double val_sensor) { Sensor(iPoint) = val_sensor; }

  /*!
   * \brief Set pressure sensor.
   * \param[in] iPoint - Point index.
   * \param[in] val_sensor - Value of the pressure sensor.
   * \param[in] iSpecies - Index of the species.
   */
  inline virtual void SetSensor(Idx_t iPoint, su2double val_sensor, Idx_t iSpecies) {}

  /*!
   * \brief Get the pressure sensor.
   * \param[in] iPoint - Point index.
   * \return Value of the pressure sensor.
   */
  inline su2double GetSensor(Idx_t iPoint) const { return Sensor(iPoint); }

  /*!
   * \brief Get the pressure sensor.
   * \param[in] iPoint - Point index.
   * \param[in] iSpecies - index of species
   * \return Value of the pressure sensor.
   */
  inline virtual su2double GetSensor(Idx_t iPoint, Idx_t iSpecies) const { return 0.0; }

  /*!
   * \brief Set the value of the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_undivided_laplacian - Value of the undivided solution for the index <i>iVar</i>.
   */
  inline void SetUndivided_Laplacian(Idx_t iPoint, Idx_t iVar, su2double val_undivided_laplacian) {
    Undivided_Laplacian(iPoint,iVar) = val_undivided_laplacian;
  }

  /*!
   * \brief Add the value of the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void AddUnd_Lapl(Idx_t iPoint, const su2double *val_und_lapl) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Undivided_Laplacian(iPoint, iVar) += val_und_lapl[iVar];
  }

  /*!
   * \brief Subtract the value of the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void SubtractUnd_Lapl(Idx_t iPoint, const su2double *val_und_lapl) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      Undivided_Laplacian(iPoint, iVar) -= val_und_lapl[iVar];
  }

  /*!
   * \brief Subtract the value of the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Variable of the undivided laplacian.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void SubtractUnd_Lapl(Idx_t iPoint, Idx_t iVar, su2double val_und_lapl) {
    Undivided_Laplacian(iPoint, iVar) -= val_und_lapl;
  }

  /*!
   * \brief Set the undivided laplacian of the solution to zero.
   */
  void SetUnd_LaplZero();

  /*!
   * \brief Set a value to the undivided laplacian.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Variable of the undivided laplacian.
   * \param[in] val_und_lapl - Value of the undivided laplacian.
   */
  inline void SetUnd_Lapl(Idx_t iPoint, Idx_t iVar, su2double val_und_lapl) {
    Undivided_Laplacian(iPoint, iVar) = val_und_lapl;
  }

  /*!
   * \brief Get the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \return Pointer to the undivided laplacian vector.
   */
  inline su2double *GetUndivided_Laplacian(Idx_t iPoint) { return Undivided_Laplacian[iPoint]; }

  /*!
   * \brief Get the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Variable of the undivided laplacian.
   * \return Value of the undivided laplacian vector.
   */
  inline su2double GetUndivided_Laplacian(Idx_t iPoint, Idx_t iVar) const { return Undivided_Laplacian(iPoint, iVar); }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow density.
   */
  inline virtual su2double GetDensity(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Old value of the flow density.
   */
  inline virtual su2double GetDensity_Old(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow density.
   */
  inline virtual su2double GetDensity(Idx_t iPoint, Idx_t val_iSpecies) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline virtual su2double GetMassFraction(Idx_t iPoint, Idx_t val_Species) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow energy.
   */
  inline virtual su2double GetEnergy(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Pointer to the force projection vector.
   */
  inline virtual su2double *GetForceProj_Vector(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Pointer to the objective function source.
   */
  inline virtual su2double *GetObjFuncSource(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Pointer to the internal boundary vector.
   */
  inline virtual su2double *GetIntBoundary_Jump(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the eddy viscosity.
   */
  inline virtual su2double GetEddyViscosity(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow enthalpy.
   */
  inline virtual su2double GetEnthalpy(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow pressure.
   */
  inline virtual su2double GetPressure(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline virtual su2double GetProjVel(Idx_t iPoint, const su2double *val_vector) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Direction of projection.
   * \param[in] val_species - Index of the desired species.
   * \return Value of the projected velocity.
   */
  inline virtual su2double GetProjVel(Idx_t iPoint, su2double *val_vector, Idx_t val_species) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the sound speed.
   */
  inline virtual su2double GetSoundSpeed(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the beta for the incompressible flow.
   */
  inline virtual su2double GetBetaInc2(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the temperature.
   */
  inline virtual su2double GetTemperature(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the vibrational-electronic temperature.
   */
  inline virtual su2double GetTemperature_ve(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (trans.-rot.).
   * \param[in] iPoint - Point index.
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  inline virtual su2double GetRhoCv_tr(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (vib.-el.).
   * \param[in] iPoint - Point index.
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  inline virtual su2double GetRhoCv_ve(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline virtual su2double GetVelocity(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Norm 2 of the velocity vector.
   */
  inline virtual su2double GetVelocity2(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Norm 2 of the velocity vector of Fluid val_species.
   */
  inline virtual su2double GetVelocity2(Idx_t iPoint, Idx_t val_species) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return The laminar viscosity of the flow.
   */
  inline virtual su2double GetLaminarViscosity(Idx_t iPoint) const { return 0.0; }


  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return The laminar viscosity of the flow.
   */
  inline virtual su2double GetLaminarViscosity(Idx_t iPoint, Idx_t iSpecies) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the species diffusion coefficient.
   */
  inline virtual su2double* GetDiffusionCoeff(Idx_t iPoint) {return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the thermal conductivity (translational/rotational)
   */
  inline virtual su2double GetThermalConductivity(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the specific heat at constant P
   */
  inline virtual su2double GetSpecificHeatCp(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the specific heat at constant V
   */
  inline virtual su2double GetSpecificHeatCv(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the thermal conductivity (vibrational)
   */
  inline virtual su2double GetThermalConductivity_ve(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Sets separation intermittency
   */
  inline virtual void SetGammaSep(Idx_t iPoint, su2double gamma_sep) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Sets separation intermittency
   */
  inline virtual void SetGammaEff(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Returns intermittency
   */
  inline virtual su2double GetIntermittency(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the vorticity.
   */
  inline virtual su2double *GetVorticity(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the rate of strain magnitude.
   */
  inline virtual su2double GetStrainMag(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
   */
  inline virtual void SetForceProj_Vector(Idx_t iPoint, const su2double *val_ForceProj_Vector) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
   */
  inline virtual void SetObjFuncSource(Idx_t iPoint, const su2double *val_SetObjFuncSource) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump.
   */
  inline virtual void SetIntBoundary_Jump(Idx_t iPoint, const su2double *val_IntBoundary_Jump) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the gamma_BC of B-C transition model.
   */
  inline virtual su2double GetGammaBC(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetGammaBC(Idx_t iPoint, su2double val_gamma) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline virtual void SetEddyViscosity(Idx_t iPoint, su2double eddy_visc) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetEnthalpy(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, CFluidModel *FluidModel) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondaryVar(Idx_t iPoint, CFluidModel *FluidModel) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool Cons2PrimVar(CConfig *config, Idx_t iPoint, su2double *U, su2double *V, su2double *dPdU,
                                   su2double *dTdU, su2double *dTvedU) { return false; }
  /*!
   * \brief A virtual member.
   */
  inline virtual void Prim2ConsVar(CConfig *config, Idx_t iPoint, su2double *V, su2double *U) { }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, su2double SharpEdge_Distance, bool check, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, su2double Density_Inf, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(Idx_t iPoint, su2double Density_Inf, su2double Viscosity_Inf,
                                 su2double eddy_visc, su2double turb_ke, CConfig *config) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrimitive(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(Idx_t iPoint, Idx_t iVar, su2double val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(Idx_t iPoint, const su2double *val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetPrimitive(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetSecondary(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondary(Idx_t iPoint, Idx_t iVar, su2double val_secondary) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondary(Idx_t iPoint, const su2double *val_secondary) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdPdrho_e(Idx_t iPoint, su2double dPdrho_e) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdPde_rho(Idx_t iPoint, su2double dPde_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdTdrho_e(Idx_t iPoint, su2double dTdrho_e) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdTde_rho(Idx_t iPoint, su2double dTde_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Setdmudrho_T(Idx_t iPoint, su2double dmudrho_T) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdmudT_rho(Idx_t iPoint, su2double dmudT_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Setdktdrho_T(Idx_t iPoint, su2double dktdrho_T) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdktdT_rho(Idx_t iPoint, su2double dktdT_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetSecondary(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetDensity(Idx_t iPoint, su2double val_density) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetDensity(Idx_t iPoint) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPressure(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelocity(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetBetaInc2(Idx_t iPoint, su2double val_betainc2) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_phi - Value of the adjoint velocity.
   */
  inline virtual void SetPhi_Old(Idx_t iPoint, const su2double *val_phi) {}

  /*!
   * \brief A virtual member.
   * \param[in] Gamma - Ratio of Specific heats
   */
  inline virtual bool SetPressure(Idx_t iPoint, su2double Gamma) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config
   */
  inline virtual bool SetPressure(Idx_t iPoint, CConfig *config) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPressure(Idx_t iPoint, su2double Gamma, su2double turb_ke) { return false; }

  /*!
   * \brief Calculates vib.-el. energy per mass, \f$e^{vib-el}_s\f$, for input species (not including KE)
   */
  inline virtual su2double CalcEve(Idx_t iPoint, su2double *V, CConfig *config, Idx_t val_Species) { return 0.0; }

  /*!
   * \brief Calculates enthalpy per mass, \f$h_s\f$, for input species (not including KE)
   */
  inline virtual su2double CalcHs(Idx_t iPoint, su2double *V, CConfig *config, Idx_t val_Species) { return 0.0; }

  /*!
   * \brief Calculates enthalpy per mass, \f$Cv_s\f$, for input species (not including KE)
   */
  inline virtual su2double CalcCvve(Idx_t iPoint, su2double val_Tve, CConfig *config, Idx_t val_Species) { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dPdU
   */
  inline virtual void CalcdPdU(Idx_t iPoint, su2double *V, CConfig *config, su2double *dPdU) {}

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dTdU
   */
  inline virtual void CalcdTdU(Idx_t iPoint, su2double *V, CConfig *config, su2double *dTdU) {}

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dTdU
   */
  inline virtual void CalcdTvedU(Idx_t iPoint, su2double *V, CConfig *config, su2double *dTdU) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdPdU(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdTdU(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdTvedU(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] Gamma - Ratio of Specific heats
   */
  inline virtual void SetDeltaPressure(Idx_t iPoint, const su2double *val_velocity, su2double Gamma) {}

  /*!
   * \brief A virtual member.
   * \param[in] Gamma - Ratio of specific heats.
   */
  inline virtual bool SetSoundSpeed(Idx_t iPoint, su2double Gamma) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual bool SetSoundSpeed(Idx_t iPoint, CConfig *config) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetSoundSpeed(Idx_t iPoint) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] Gas_Constant - Value of the Gas Constant
   */
  inline virtual bool SetTemperature(Idx_t iPoint, su2double Gas_Constant) { return false; }

  /*!
   * \brief Sets the vibrational electronic temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline virtual bool SetTemperature_ve(Idx_t iPoint, su2double val_Tve) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual bool SetTemperature(Idx_t iPoint, CConfig *config) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual void SetPrimitive(Idx_t iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   * \param[in] Coord - Physical coordinates.
   */
  inline virtual void SetPrimitive(Idx_t iPoint, CConfig *config, su2double *Coord) {}

  /*!
   * \brief A virtual member.
   * \param[in] Temperature_Wall - Value of the Temperature at the wall
   */
  inline virtual void SetWallTemperature(Idx_t iPoint, su2double Temperature_Wall) {}

  /*!
   * \brief A virtual member.
   * \param[in] Temperature_Wall - Value of the Temperature at the wall
   */
  inline virtual void SetWallTemperature(Idx_t iPoint, su2double* Temperature_Wall) {}

  /*!
   * \brief Set the thermal coefficient.
   * \param[in] config - Configuration parameters.
   */
  inline virtual void SetThermalCoeff(Idx_t iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetStress_FEM(Idx_t iPoint, Idx_t iVar, su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void AddStress_FEM(Idx_t iPoint, Idx_t iVar, su2double val_stress) {}

  /*!
   * \brief A virtual member.

   */
  inline virtual su2double *GetStress_FEM(Idx_t iPoint) {return nullptr;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVonMises_Stress(Idx_t iPoint, su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetVonMises_Stress(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_SurfaceLoad_Res(Idx_t iPoint, const su2double *val_surfForce) {}

  /*!
   * \brief  A virtual member.
   */
  inline virtual void Set_SurfaceLoad_Res(Idx_t iPoint, Idx_t iVar, su2double val_surfForce) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_SurfaceLoad_Res(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_SurfaceLoad_Res(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   */
  virtual void Set_SurfaceLoad_Res_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_SurfaceLoad_Res_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_BodyForces_Res(Idx_t iPoint, const su2double *val_bodyForce) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_BodyForces_Res(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_BodyForces_Res(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_FlowTraction(Idx_t iPoint, const su2double *val_flowTraction) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_FlowTraction(Idx_t iPoint, const su2double *val_flowTraction) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_FlowTraction(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  virtual void Set_FlowTraction_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_FlowTraction_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_FlowTraction() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_isVertex(Idx_t iPoint, bool isVertex) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool Get_isVertex(Idx_t iPoint) const { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelocity2(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline virtual void SetVelocity_Old(Idx_t iPoint, const su2double *val_velocity) {}

  /*!
   * \brief A virtual member.
   * \param[in] laminarViscosity
   */
  inline virtual void SetLaminarViscosity(Idx_t iPoint, su2double laminarViscosity) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetLaminarViscosity(Idx_t iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] thermalConductivity
   */
  inline virtual void SetThermalConductivity(Idx_t iPoint, su2double thermalConductivity) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetThermalConductivity(Idx_t iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] Cp - Constant pressure specific heat.
   */
  inline virtual void SetSpecificHeatCp(Idx_t iPoint, su2double Cp) {}

  /*!
   * \brief A virtual member.
   * \param[in] Cv - Constant volume specific heat.
   */
  inline virtual void SetSpecificHeatCv(Idx_t iPoint, su2double Cv) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetVorticity_StrainMag() { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelSolutionDVector(Idx_t iPoint) {}

  /*!
   * \brief A virtual member.
   */
  virtual void SetGradient_PrimitiveZero() {}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] val_value - Value to add to the gradient of the primitive variables.
   */
  inline virtual void AddGradient_Primitive(Idx_t iPoint, Idx_t iVar, Idx_t iDim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
   */
  inline virtual void SubtractGradient_Primitive(Idx_t iPoint, Idx_t iVar, Idx_t iDim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double GetGradient_Primitive(Idx_t iPoint, Idx_t iVar, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double GetLimiter_Primitive(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline virtual void SetGradient_Primitive(Idx_t iPoint, Idx_t iVar, Idx_t iDim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_value - Value of the gradient.
   */
  inline virtual void SetLimiter_Primitive(Idx_t iPoint, Idx_t iVar, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double **GetGradient_Primitive(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double *GetLimiter_Primitive(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_density - Value of the density.
   * \param[in] val_dist - Value of the distance to the wall.
   */
  inline virtual void SetBlendingFunc(Idx_t iPoint, su2double val_viscosity, su2double val_dist, su2double val_density) {}

  /*!
   * \brief Get the first blending function of the SST model.
   */
  inline virtual su2double GetF1blending(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief Get the second blending function of the SST model.
   */
  inline virtual su2double GetF2blending(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline virtual su2double GetCrossDiff(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief Get the value of the eddy viscosity.
   * \return the value of the eddy viscosity.
   */
  inline virtual su2double GetmuT(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief Set the value of the eddy viscosity.
   * \param[in] val_muT
   */
  inline virtual void SetmuT(Idx_t iPoint, su2double val_muT) {}

  /*!
   * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iSpecies - Value of iSpecies to which the eigenvalue belongs
   */
  inline virtual void AddMax_Lambda_Inv(Idx_t iPoint, su2double val_max_lambda, Idx_t iSpecies) {}

  /*!
   * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iSpecies - Value of iSpecies to which the eigenvalue belongs
   */
  inline virtual void AddMax_Lambda_Visc(Idx_t iPoint, su2double val_max_lambda, Idx_t iSpecies) {}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_source - Value of the harmonic balance source.
   */
  inline virtual void SetHarmonicBalance_Source(Idx_t iPoint, Idx_t iVar, su2double val_source) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetHarmonicBalance_Source(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Set the Eddy Viscosity Sensitivity of the problem.
   * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
   * \param[in] numTotalVar - Number of variables.
   */
  inline virtual void SetEddyViscSens(Idx_t iPoint, const su2double *val_EddyViscSens, Idx_t numTotalVar) {}

  /*!
   * \brief Get the Eddy Viscosity Sensitivity of the problem.
   * \return Pointer to the Eddy Viscosity Sensitivity.
   */
  inline virtual su2double *GetEddyViscSens(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Set the direct solution for the adjoint solver.
   * \param[in] solution_direct - Value of the direct solution.
   */
  inline virtual void SetSolution_Direct(Idx_t iPoint, const su2double *solution_direct) {}

  /*!
   * \brief A virtual member. Get the direct solution for the adjoint solver.
   * \return Pointer to the direct solution vector.
   */
  inline virtual su2double *GetSolution_Direct(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Set the restart geometry (coordinate of the converged solution)
   * \param[in] val_coordinate_direct - Value of the restart coordinate.
   */
  inline virtual void SetGeometry_Direct(Idx_t iPoint, const su2double *val_coordinate_direct) {}

  /*!
   * \brief A virtual member. Get the restart geometry (coordinate of the converged solution).
   * \return Pointer to the restart coordinate vector.
   */
  inline virtual su2double *GetGeometry_Direct(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Get the restart geometry (coordinate of the converged solution).
   * \return Coordinate of the direct solver restart for .
   */
  inline virtual su2double GetGeometry_Direct(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetSolution_Geometry(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] solution - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Geometry(Idx_t iPoint, const su2double *solution_geometry) {}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] solution - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Geometry(Idx_t iPoint, Idx_t iVar, su2double solution_geometry) {}

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetGeometry_CrossTerm_Derivative(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] solution - Solution of the problem (acceleration).
   */
  inline virtual void SetGeometry_CrossTerm_Derivative(Idx_t iPoint, Idx_t iDim, su2double der) {}

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetGeometry_CrossTerm_Derivative_Flow(Idx_t iPoint, Idx_t iVar) const { return 0.0;}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] solution - Solution of the problem (acceleration).
   */
  inline virtual void SetGeometry_CrossTerm_Derivative_Flow(Idx_t iPoint, Idx_t iDim, su2double der) {}

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_OldSolution_Geometry() {}

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_OldSolution_Geometry(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_BGSSolution(Idx_t iPoint, Idx_t iDim, su2double solution) {}

  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  virtual void Set_BGSSolution_k();
  
  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  inline void Set_BGSSolution_k(Idx_t iPoint, Idx_t iVar, su2double val_var) {
    Solution_BGS_k(iPoint,iVar) = val_var;
  }

  /*!
   * \brief Get the value of the solution in the previous BGS subiteration.
   * \param[out] val_solution - solution in the previous BGS subiteration.
   */
  inline virtual su2double Get_BGSSolution_k(Idx_t iPoint, Idx_t iVar) const {
    return Solution_BGS_k(iPoint,iVar);
  }

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_BGSSolution(Idx_t iPoint, Idx_t iDim) const {return 0.0;}

  /*!
   * \brief  A virtual member. Set the contribution of crossed terms into the derivative.
   */
  inline virtual void SetCross_Term_Derivative(Idx_t iPoint, Idx_t iVar, su2double der) {}

  /*!
   * \brief  A virtual member. Get the contribution of crossed terms into the derivative.
   * \return The contribution of crossed terms into the derivative.
   */
  inline virtual su2double GetCross_Term_Derivative(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the direct velocity solution for the adjoint solver.
   * \param[in] solution_direct - Value of the direct velocity solution.
   */
  inline virtual void SetSolution_Vel_Direct(Idx_t iPoint, const su2double *sol) {}

  /*!
   * \brief A virtual member. Set the direct acceleration solution for the adjoint solver.
   * \param[in] solution_direct - Value of the direct acceleration solution.
   */
  inline virtual void SetSolution_Accel_Direct(Idx_t iPoint, const su2double *sol) {}

  /*!
   * \brief A virtual member. Get the direct velocity solution for the adjoint solver.
   * \return Pointer to the direct velocity solution vector.
   */
  inline virtual su2double* GetSolution_Vel_Direct(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Get the direct acceleraction solution for the adjoint solver.
   * \return Pointer to the direct acceleraction solution vector.
   */
  inline virtual su2double* GetSolution_Accel_Direct(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief Set the value of the velocity (Structural Analysis).
   * \param[in] solution - Solution of the problem (velocity).
   */
  inline virtual void SetSolution_Vel(Idx_t iPoint, const su2double *solution) {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_vel - Value of the solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Vel(Idx_t iPoint, Idx_t iVar, su2double solution_vel) {}

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   * \param[in] solution_vel_time_n - Value of the old solution.
   */
  inline virtual void SetSolution_Vel_time_n(Idx_t iPoint, const su2double *solution_vel_time_n) {}

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   */
  inline virtual void SetSolution_Vel_time_n() {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_vel_time_n - Value of the old solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Vel_time_n(Idx_t iPoint, Idx_t iVar, su2double solution_vel_time_n) {}

  /*!
   * \brief Get the solution at time n.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution_time_n(Idx_t iPoint, Idx_t iVar) const { return Solution_time_n(iPoint,iVar); }

  /*!
   * \brief Get the solution at time n-1.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution_time_n1(Idx_t iPoint, Idx_t iVar) const { return Solution_time_n1(iPoint,iVar); }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetSolution_Vel(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline virtual su2double *GetSolution_Vel(Idx_t iPoint) {return nullptr; }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Vel_time_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Vel_time_n(Idx_t iPoint) { return nullptr; }


  /*!
   * \brief Set the value of the acceleration (Structural Analysis).
   * \param[in] solution_accel - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Accel(Idx_t iPoint, const su2double *solution_accel) {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_accel - Value of the solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Accel(Idx_t iPoint, Idx_t iVar, su2double solution_accel) {}

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   * \param[in] solution_accel_time_n - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Accel_time_n(Idx_t iPoint, const su2double *solution_accel_time_n) {}

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   */
  inline virtual void SetSolution_Accel_time_n() {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_accel_time_n - Value of the old solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Accel_time_n(Idx_t iPoint, Idx_t iVar, su2double solution_accel_time_n) {}

  /*!
   * \brief Get the acceleration (Structural Analysis).
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetSolution_Accel(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline virtual su2double *GetSolution_Accel(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief Get the acceleration of the nodes (Structural Analysis) at time n.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Accel_time_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Accel_time_n(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_OldSolution_Vel() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_OldSolution_Accel() {}

  /*!
   * \brief  A virtual member. Set the value of the solution predictor.
   */
  inline virtual void SetSolution_Pred(Idx_t iPoint) {}

  /*!
   * \brief  A virtual member. Set the value of the old solution.
   * \param[in] solution_pred - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred(Idx_t iPoint, const su2double *solution_pred) {}

  /*!
   * \brief  A virtual member. Set the value of the solution predicted.
   * \param[in] solution_old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred(Idx_t iPoint, Idx_t iVar, su2double solution_pred) {}

  /*!
   * \brief  A virtual member. Get the value of the solution predictor.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Pred(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief  A virtual member. Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Pred(Idx_t iPoint) {return nullptr; }

  /*!
   * \brief  A virtual member. Set the value of the solution predictor.
   */
  inline virtual void SetSolution_Pred_Old(Idx_t iPoint) {}

  /*!
   * \brief  A virtual member. Set the value of the old solution.
   * \param[in] solution_pred_Old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred_Old(Idx_t iPoint, const su2double *solution_pred_Old) {}

  /*!
   * \brief  A virtual member. Set the value of the old solution predicted.
   * \param[in] solution_pred_old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred_Old(Idx_t iPoint, Idx_t iVar, su2double solution_pred_old) {}

  /*!
   * \brief  A virtual member. Get the value of the solution predictor.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Pred_Old(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief  A virtual member. Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Pred_Old(Idx_t iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetReference_Geometry(Idx_t iPoint, Idx_t iVar, su2double ref_geometry) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetReference_Geometry(Idx_t iPoint) {return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrestretch(Idx_t iPoint, Idx_t iVar, su2double val_prestretch) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetPrestretch(Idx_t iPoint) {return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrestretch(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetReference_Geometry(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the original coordinate iDim.
   */
  inline virtual su2double GetMesh_Coord(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the undeformed coordinates.
   * \return Pointer to the reference coordinates.
   */
  inline virtual const su2double *GetMesh_Coord(Idx_t iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member. Set the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \param[in] val_coord - Value of Mesh_Coord[nDim]
   */
  inline virtual void SetMesh_Coord(Idx_t iPoint, Idx_t iDim, su2double val_coord) { }

    /*!
   * \brief A virtual member. Get the value of the wall distance in reference coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the wall distance in reference coordinates.
   */
  inline virtual su2double GetWallDistance(Idx_t iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the value of the wall distance in reference coordinates.
   * \param[in] val_dist - Value of wall distance.
   */
  inline virtual void SetWallDistance(Idx_t iPoint, su2double val_dist) { }

  /*!
   * \brief A virtual member. Register the reference coordinates of the mesh.
   * \param[in] input - Defines whether we are registering the variable as input or as output.
   */
  inline virtual void Register_MeshCoord(bool input) { }

  /*!
   * \brief A virtual member. Recover the value of the adjoint of the mesh coordinates.
   */
  inline virtual void GetAdjoint_MeshCoord(Idx_t iPoint, su2double *adj_mesh) const { }

  /*!
   * \brief A virtual member. Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline virtual su2double GetBound_Disp(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the boundary displacement.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline virtual void SetBound_Disp(Idx_t iPoint, const su2double *val_BoundDisp) { }


  /*!
   * \brief A virtual member. Set the boundary displacement.
   * \param[in] iDim - Index of the dimension of interest.
   * \param[in] val_BoundDisp - Value of the boundary displacements.
   */
  inline virtual void SetBound_Disp(Idx_t iPoint, Idx_t iDim, const su2double val_BoundDisp) { }

  /*!
   * \brief A virtual member. Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline virtual const su2double* GetBoundDisp_Direct(Idx_t iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member. Set the solution for the boundary displacements.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline virtual void SetBoundDisp_Direct(Idx_t iPoint, const su2double *val_BoundDisp) { }

  /*!
   * \brief Set the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] val_sens - Pointer to the sensitivities of the boundary displacements.
   */
  inline virtual void SetBoundDisp_Sens(Idx_t iPoint, const su2double *val_sens) { }

  /*!
   * \brief A virtual member. Get the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord_Sens[nDim]
   * \return Value of the original Mesh_Coord_Sens iDim.
   */
  inline virtual su2double GetBoundDisp_Sens(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Register the boundary displacements of the mesh.
   * \param[in] input - Defines whether we are registering the variable as input or as output.
   */
  inline virtual void Register_BoundDisp(bool input) { }

  /*!
   * \brief A virtual member. Recover the value of the adjoint of the boundary displacements.
   */
  inline virtual void GetAdjoint_BoundDisp(Idx_t iPoint, su2double *adj_disp) const { }

   /*!
    * \brief A virtual member.
    */
  inline virtual void Register_femSolution_time_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void RegisterSolution_Vel(bool input) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void RegisterSolution_Vel_time_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void RegisterSolution_Accel(bool input) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void RegisterSolution_Accel_time_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Vel(Idx_t iPoint, const su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void RegisterFlowTraction() { }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double ExtractFlowTraction_Sensitivity(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Vel(Idx_t iPoint, su2double *adj_sol) const {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Vel_time_n(Idx_t iPoint, const su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Vel_time_n(Idx_t iPoint, su2double *adj_sol) const {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Accel(Idx_t iPoint, const su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Accel(Idx_t iPoint, su2double *adj_sol) const {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Accel_time_n(Idx_t iPoint, const su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Accel_time_n(Idx_t iPoint, su2double *adj_sol) const {}

  /*!
   * \brief Register the variables in the solution array as input/output variable.
   * \param[in] input - input or output variables.
   */
  void RegisterSolution(bool input);

  /*!
   * \brief Register the variables in the solution array as input/output variable.
   * \param[in] input - input or output variables.
   */
  void RegisterSolution_intIndexBased(bool input);

  /*!
   * \brief Saving the adjoint vector position with respect to the solution variables.
   * \param[in] input - input or output variables.
   */
  void SetAdjIndices(bool input);

  /*!
   * \brief Register the variables in the solution_time_n array as input/output variable.
   */
  void RegisterSolution_time_n();

  /*!
   * \brief Register the variables in the solution_time_n1 array as input/output variable.
   */
  void RegisterSolution_time_n1();

  /*!
   * \brief Set the adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution(Idx_t iPoint, const su2double *adj_sol) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Set the adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_intIndexBased(Idx_t iPoint, const su2double *adj_sol) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      AD::SetDerivative(Output_AdjIndices(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution.
   * \param[out] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution(Idx_t iPoint, su2double *adj_sol) const {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution(iPoint,iVar));
  }

  /*!
   * \brief Get the adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_intIndexBased(Idx_t iPoint, su2double *adj_sol) const {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = AD::GetDerivative(Input_AdjIndices(iPoint,iVar));
  }

  /*!
   * \brief Set the adjoint values of the solution at time n.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_time_n(Idx_t iPoint, const su2double *adj_sol) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution at time n.
   * \param[out] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_time_n(Idx_t iPoint, su2double *adj_sol) const {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n(iPoint,iVar));
  }

  /*!
   * \brief Set the adjoint values of the solution at time n-1.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_time_n1(Idx_t iPoint, const su2double *adj_sol) {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n1(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution at time n-1.
   * \param[out] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_time_n1(Idx_t iPoint, su2double *adj_sol) const {
    for (Idx_t iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n1(iPoint,iVar));
  }

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline virtual void SetSensitivity(Idx_t iPoint, Idx_t iDim, su2double val) {}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline virtual su2double GetSensitivity(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  inline virtual void SetDual_Time_Derivative(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual void SetDual_Time_Derivative_n(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual su2double GetDual_Time_Derivative(Idx_t iPoint, Idx_t iVar) const {return 0.0;}

  inline virtual su2double GetDual_Time_Derivative_n(Idx_t iPoint, Idx_t iVar) const {return 0.0;}

  inline virtual void SetTauWall(Idx_t iPoint, su2double val_tau_wall) {}

  inline virtual su2double GetTauWall(Idx_t iPoint) const { return 0.0; }

  inline virtual void SetVortex_Tilting(Idx_t iPoint, su2double **PrimGrad_Flow, su2double* Vorticity, su2double LaminarViscosity) {}

  inline virtual su2double GetVortex_Tilting(Idx_t iPoint) const { return 0.0; }

  inline virtual void SetDynamic_Derivative(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual void SetDynamic_Derivative_n(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual su2double GetDynamic_Derivative(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual su2double GetDynamic_Derivative_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual void SetDynamic_Derivative_Vel(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual void SetDynamic_Derivative_Vel_n(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual su2double GetDynamic_Derivative_Vel(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual su2double GetDynamic_Derivative_Vel_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual void SetDynamic_Derivative_Accel(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual void SetDynamic_Derivative_Accel_n(Idx_t iPoint, Idx_t iVar, su2double der) {}

  inline virtual su2double GetDynamic_Derivative_Accel(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual su2double GetDynamic_Derivative_Accel_n(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual su2double GetSolution_Old_Vel(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  inline virtual su2double GetSolution_Old_Accel(Idx_t iPoint, Idx_t iVar) const { return 0.0; }

  /*!
   * \brief Set the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  virtual void SetFlowTractionSensitivity(Idx_t iPoint, Idx_t iDim, su2double val) { }

  /*!
   * \brief Get the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  virtual su2double GetFlowTractionSensitivity(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

  /*!
   * \brief Set the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \param[in] val - value of the source term
   */
  virtual void SetSourceTerm_DispAdjoint(Idx_t iPoint, Idx_t iDim, su2double val) { }

  /*!
   * \brief Get the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \return value of the source term
   */
  virtual su2double GetSourceTerm_DispAdjoint(Idx_t iPoint, Idx_t iDim) const { return 0.0; }

};
