/*!
 * \file CVariable.hpp
 * \brief Declaration and inlines of the parent class for defining problem
          variables, function definitions in file <i>CVariable.cpp</i>.
          All variables are children of at least this class.
 * \author F. Palacios, T. Economon
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

#include "../../../Common/include/parallelization/mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/containers/container_decorators.hpp"

class CFluidModel;
class CNEMOGas;

/*!
 * \class CVariable
 * \ingroup Variable
 * \brief Main class for defining the variables.
 * \author F. Palacios
 */
class CVariable {
protected:
  using VectorType = C2DContainer<unsigned long, su2double, StorageType::ColumnMajor, 64, DynamicSize, 1>;
  using MatrixType = C2DContainer<unsigned long, su2double, StorageType::RowMajor,    64, DynamicSize, DynamicSize>;

  MatrixType Solution;       /*!< \brief Solution of the problem. */
  MatrixType Solution_Old;   /*!< \brief Old solution of the problem R-K. */

  MatrixType External;       /*!< \brief External (outer) contribution in discrete adjoint multizone problems. */

  su2vector<bool> Non_Physical;  /*!< \brief Non-physical points in the solution (force first order). */
  su2vector<unsigned short>
  Non_Physical_Counter;          /*!< \brief Number of consecutive iterations that a point has been treated first-order.
                                  After a specified number of successful reconstructions, the point can be returned to second-order. */

  VectorType UnderRelaxation;  /*!< \brief Value of the under-relaxation parameter local to the control volume. */
  VectorType LocalCFL;         /*!< \brief Value of the CFL number local to the control volume. */

  MatrixType Solution_time_n;    /*!< \brief Solution of the problem at time n for dual-time stepping technique. */
  MatrixType Solution_time_n1;   /*!< \brief Solution of the problem at time n-1 for dual-time stepping technique. */
  VectorType Delta_Time;         /*!< \brief Time step. */

  CVectorOfMatrix Gradient;  /*!< \brief Gradient of the solution of the problem. */
  C3DDoubleMatrix Rmatrix;   /*!< \brief Geometry-based matrix for weighted least squares gradient calculations. */

  MatrixType Limiter;        /*!< \brief Limiter of the solution of the problem. */
  MatrixType Solution_Max;   /*!< \brief Max solution for limiter computation. */
  MatrixType Solution_Min;   /*!< \brief Min solution for limiter computation. */

  MatrixType AuxVar;             /*!< \brief Auxiliary variable for gradient computation. */
  CVectorOfMatrix Grad_AuxVar;   /*!< \brief Gradient of the auxiliary variables of the problem. */

  VectorType Max_Lambda_Inv;   /*!< \brief Maximun inviscid eingenvalue. */
  VectorType Max_Lambda_Visc;  /*!< \brief Maximun viscous eingenvalue. */
  VectorType Lambda;           /*!< \brief Value of the eingenvalue. */

  VectorType Sensor;               /*!< \brief Pressure sensor for high order central scheme and Roe dissipation. */
  MatrixType Undivided_Laplacian;  /*!< \brief Undivided laplacian of the solution. */

  MatrixType Res_TruncError;  /*!< \brief Truncation error for multigrid cycle. */
  MatrixType Residual_Old;    /*!< \brief Auxiliary structure for residual smoothing. */
  MatrixType Residual_Sum;    /*!< \brief Auxiliary structure for residual smoothing. */

  MatrixType Solution_BGS_k;     /*!< \brief Old solution container for BGS iterations. */

  su2matrix<int> AD_InputIndex;    /*!< \brief Indices of Solution variables in the adjoint vector. */
  su2matrix<int> AD_OutputIndex;   /*!< \brief Indices of Solution variables in the adjoint vector after having been updated. */

  VectorType SolutionExtra; /*!< \brief Stores adjoint solution for extra solution variables.
                                        Currently only streamwise periodic pressure-drop for massflow prescribed flows. */
  VectorType ExternalExtra; /*!< \brief External storage for the adjoint value (i.e. for the OF mainly */

  VectorType SolutionExtra_BGS_k; /*!< \brief Intermediate storage, enables cross term extraction as that is also pushed to Solution. */

 protected:
  unsigned long nPoint = 0;  /*!< \brief Number of points in the domain. */
  unsigned long nDim = 0;      /*!< \brief Number of dimension of the problem. */
  unsigned long nVar = 0;        /*!< \brief Number of variables of the problem. */
  unsigned long nPrimVar = 0;      /*!< \brief Number of primitive variables. */
  unsigned long nPrimVarGrad = 0;    /*!< \brief Number of primitives for which a gradient is computed. */
  unsigned long nSecondaryVar = 0;     /*!< \brief Number of secondary variables. */
  unsigned long nSecondaryVarGrad = 0;   /*!< \brief Number of secondaries for which a gradient is computed. */
  unsigned long nAuxVar = 0; /*!< \brief Number of auxiliary variables. */

  /*--- Only allow default construction by derived classes. ---*/
  CVariable() = default;

  inline static void AssertOverride() {
    assert(false && "A base method of CVariable was used, but it should have been overridden by the derived class.");
  }

  void RegisterContainer(bool input, su2activematrix& variable, su2matrix<int>* ad_index = nullptr) {
    const auto nPoint = variable.rows();
    SU2_OMP_FOR_STAT(roundUpDiv(nPoint,omp_get_num_threads()))
    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      for(unsigned long iVar=0; iVar<variable.cols(); ++iVar) {

        if (input) AD::RegisterInput(variable(iPoint,iVar));
        else AD::RegisterOutput(variable(iPoint,iVar));

        if (ad_index) AD::SetIndex((*ad_index)(iPoint,iVar), variable(iPoint,iVar));
      }
    }
    END_SU2_OMP_FOR
  }

  void RegisterContainer(bool input, su2activematrix& variable, su2matrix<int>& ad_index) {
    RegisterContainer(input, variable, &ad_index);
  }

public:
  /*--- Disable copy and assignment. ---*/
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
  CVariable(unsigned long npoint, unsigned long nvar, const CConfig *config);

  /*!
   * \overload
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] adjoint - True if derived class is an adjoint variable.
   */
  CVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config, bool adjoint = false);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CVariable() = default;

  /*!
   * \brief Get the number of auxiliary variables.
   */
  inline unsigned long GetnAuxVar() const { return nAuxVar; }

  /*!
   * \brief Set the value of the solution, all variables.
   * \param[in] iPoint - Point index.
   * \param[in] solution - Solution of the problem.
   */
  inline void SetSolution(unsigned long iPoint, const su2double *solution) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution(iPoint,iVar) = solution[iVar];
  }

  /*!
   * \brief Set the value of the solution, one variable.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution for the index <i>iVar</i>.
   */
  inline void SetSolution(unsigned long iPoint, unsigned long iVar, su2double solution) { Solution(iPoint,iVar) = solution; }

  /*!
   * \brief Add the value of the solution vector to the previous solution (incremental approach).
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution for the index <i>iVar</i>.
   */
  inline void Add_DeltaSolution(unsigned long iPoint, unsigned long iVar, su2double solution) { Solution(iPoint,iVar) += solution; }

  /*!
   * \brief Set the value of the non-physical point.
   * \param[in] iPoint - Point index.
   * \param[in] value - identification of the non-physical point.
   */
  inline void SetNon_Physical(unsigned long iPoint, bool val_value) {
    if (val_value) {
      Non_Physical(iPoint) = val_value;
      Non_Physical_Counter(iPoint) = 0;
    } else {
      Non_Physical_Counter(iPoint)++;
      if (Non_Physical_Counter(iPoint) > 20) {
        Non_Physical(iPoint) = false;
      }
    }
  }

  /*!
   * \brief Get the value of the non-physical boolean at a point.
   * \param[in] iPoint - Point index.
   * \return Value of the Non-physical point.
   */
  inline bool GetNon_Physical(unsigned long iPoint) { return Non_Physical(iPoint); }

  /*!
   * \brief Get the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution(unsigned long iPoint, unsigned long iVar) const { return Solution(iPoint,iVar); }

  /*!
   * \brief Get the old solution of the problem (Runge-Kutta method)
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Old(unsigned long iPoint, unsigned long iVar) const { return Solution_Old(iPoint,iVar); }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] iPoint - Point index.
   * \param[in] solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Old(unsigned long iPoint, const su2double *solution_old) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution_Old(iPoint,iVar) = solution_old[iVar];
  }

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] iPoint - Point index.
   * \param[in] solution_old - Value of the old solution for the index <i>iVar</i>.
   */
  inline void SetSolution_Old(unsigned long iPoint, unsigned long iVar, su2double solution_old) { Solution_Old(iPoint,iVar) = solution_old; }

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
  inline void Set_Solution_time_n(unsigned long iPoint, const su2double* val_sol) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_time_n(iPoint,iVar) = val_sol[iVar];
  }

  /*!
   * \brief Set the variable solution at time n-1.
   * \param[in] iPoint - Point index.
   */
  inline void Set_Solution_time_n1(unsigned long iPoint, const su2double* val_sol) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution_time_n1(iPoint,iVar) = val_sol[iVar];
  }

  /*!
   * \brief Set the variable solution at time n.
   * \param[in] iPoint - Point index.
   */
  inline void Set_Solution_time_n(unsigned long iPoint, unsigned long iVar, su2double val_sol) {
    Solution_time_n(iPoint,iVar) = val_sol;
  }

  /*!
   * \brief Set the variable solution at time n-1.
   * \param[in] iPoint - Point index.
   */
  inline void Set_Solution_time_n1(unsigned long iPoint, unsigned long iVar, su2double val_sol) {
    Solution_time_n1(iPoint,iVar) = val_sol;
  }

  /*!
   * \brief Virtual Member. Specify a vector to set the velocity components of the solution.
   *        Multiplied by density for compressible cases.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline virtual void SetVelSolutionVector(unsigned long iPoint, const su2double *val_vector) { }

  /*!
   * \brief Add a value to the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Number of the variable.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline void AddSolution(unsigned long iPoint, unsigned long iVar, su2double solution) {
    Solution(iPoint, iVar) = Solution_Old(iPoint, iVar) + solution;
  }

  /*!
   * \brief Add a value to the solution.
   * \param[in] iPoint - Point index.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline void AddSolution(unsigned long iPoint, const su2double *solution) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution(iPoint, iVar) += solution[iVar];
  }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_New(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual su2double GetRoe_Dissipation(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetRoe_Dissipation(unsigned long iPoint, su2double val_dissipation) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetRoe_Dissipation_FD(unsigned long iPoint, su2double val_wall_dist) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_delta - A scalar measure of the grid size
   * \param[in] val_const_DES - The DES constant (C_DES)
   */
  inline virtual void SetRoe_Dissipation_NTS(unsigned long iPoint, su2double val_delta, su2double val_const_DES) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual su2double GetDES_LengthScale(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetDES_LengthScale(unsigned long iPoint, su2double val_des_lengthscale) {}

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
   * \brief Set Dual-time derivative contributions to the external.
   */
  inline virtual void Set_External_To_DualTimeDer() { SetExternalZero(); }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Number of the variable.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline virtual void AddSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution) {}

  /*!
   * \brief Add a value to the External vector.
   * \param[in] iPoint - Point index.
   * \param[in] val_sol - vector that has to be added component-wise
   */
  inline void Add_External(unsigned long iPoint, const su2double* val_sol) {
    for(unsigned long iVar = 0; iVar < nVar; iVar++) External(iPoint,iVar) += val_sol[iVar];
  }

  /*!
   * \brief Store the adjoint solution of the extra adjoint into the external container.
   */
  void Set_ExternalExtra_To_SolutionExtra() {
    assert(SolutionExtra.size() == ExternalExtra.size());
    for (auto iEntry = 0ul; iEntry < SolutionExtra.size(); iEntry++)
      ExternalExtra[iEntry] = SolutionExtra[iEntry];
  }

  /*!
   * \brief Add the external contribution to the solution for the extra adjoint solutions.
   */
  void Add_ExternalExtra_To_SolutionExtra() {
    assert(SolutionExtra.size() == ExternalExtra.size());
    for (auto iEntry = 0ul; iEntry < SolutionExtra.size(); iEntry++)
      SolutionExtra[iEntry] += ExternalExtra[iEntry];
  }

  /*!
   * \brief Return the extra adjoint solution.
   * \return Reference to extra adjoint solution.
   */
  inline VectorType& GetSolutionExtra() { return SolutionExtra; }
  inline const VectorType& GetSolutionExtra() const { return SolutionExtra; }

  /*!
   * \brief Update the variables using a conservative format.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] solution - Value of the solution change.
   * \param[in] lowerlimit - Lower value for Solution clipping.
   * \param[in] upperlimit - Upper value for Solution clipping.
   * \param[in] Sol2Conservative - Factor multiplied to Solution to get transported variable.
   * \param[in] Sol2Conservative_old - Factor multiplied to Solution to get transported variable, of the previous Iteration.
   */
  inline void AddClippedSolution(unsigned long iPoint, unsigned long iVar, su2double solution,
                                 su2double lowerlimit, su2double upperlimit,
                                 su2double Sol2Conservative = 1.0, su2double Sol2Conservative_old = 1.0) {

    su2double val_new = (Solution_Old(iPoint,iVar)*Sol2Conservative_old + solution)/Sol2Conservative;
    Solution(iPoint,iVar) = min(max(val_new, lowerlimit), upperlimit);
  }

  /*!
   * \brief Get the entire solution of the problem.
   * \return Reference to the solution matrix.
   */
  inline const MatrixType& GetSolution() const { return Solution; }
  inline MatrixType& GetSolution() { return Solution; }

  /*!
   * \brief Get the solution of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the solution vector.
   */
  inline su2double *GetSolution(unsigned long iPoint) { return Solution[iPoint]; }

  /*!
   * \brief Get the old solution of the problem (Runge-Kutta method)
   * \param[in] iPoint - Point index.
   * \return Pointer to the old solution vector.
   */
  inline su2double *GetSolution_Old(unsigned long iPoint) { return Solution_Old[iPoint]; }

  /*!
   * \brief Get the external contributions of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the External row for iPoint.
   */
  inline const su2double *Get_External(unsigned long iPoint) const { return External[iPoint]; }
  inline const MatrixType& Get_External() const { return External; }

  /*!
   * \brief Get the solution at time n.
   * \param[in] iPoint - Point index.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_time_n(unsigned long iPoint) { return Solution_time_n[iPoint]; }
  inline MatrixType& GetSolution_time_n() { return Solution_time_n; }

  /*!
   * \brief Get the solution at time n-1.
   * \param[in] iPoint - Point index.
   * \return Pointer to the solution (at time n-1) vector.
   */
  inline su2double *GetSolution_time_n1(unsigned long iPoint) { return Solution_time_n1[iPoint]; }
  inline MatrixType& GetSolution_time_n1() { return Solution_time_n1; }

  /*!
   * \brief Set the value of the old residual.
   * \param[in] iPoint - Point index.
   * \param[in] val_residual_old - Pointer to the residual vector.
   */
  inline void SetResidual_Old(unsigned long iPoint, const su2double *val_residual_old) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Residual_Old(iPoint,iVar) = val_residual_old[iVar];
  }

  /*!
   * \brief Add a value to the summed residual vector.
   * \param[in] iPoint - Point index.
   * \param[in] val_residual - Pointer to the residual vector.
   */
  inline void AddResidual_Sum(unsigned long iPoint, const su2double *val_residual) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Residual_Sum(iPoint,iVar) += val_residual[iVar];
  }

  /*!
   * \brief Set summed residual vector to zero value.
   */
  inline void SetResidualSumZero(unsigned long iPoint) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Residual_Sum(iPoint,iVar) = 0.0;
  }

  /*!
   * \brief Get the value of the summed residual.
   * \param[in] iPoint - Point index.
   * \return Pointer to the summed residual.
   */
  inline su2double *GetResidual_Sum(unsigned long iPoint) { return Residual_Sum[iPoint]; }

  /*!
   * \brief Get the value of the old residual.
   * \param[in] iPoint - Point index.
   * \return Pointer to the old residual.
   */
  inline su2double *GetResidual_Old(unsigned long iPoint) { return Residual_Old[iPoint]; }

  /*!
   * \brief Get the value of the summed residual.
   * \param[in] iPoint - Point index.
   * \param[in] val_residual - Pointer to the summed residual.
   */
  inline void GetResidual_Sum(unsigned long iPoint, su2double *val_residual) const {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] = Residual_Sum(iPoint,iVar);
  }

  /*!
   * \brief Set the value of the under-relaxation parameter for the current control volume (CV).
   * \param[in] iPoint - Point index.
   * \param[in] val_under_relaxation - the input value of the under-relaxation parameter for this CV.
   */
  inline void SetUnderRelaxation(unsigned long iPoint, su2double val_under_relaxation) { UnderRelaxation(iPoint) = val_under_relaxation; }

  /*!
   * \brief Get the value of the under-relaxation parameter for the current control volume (CV).
   * \param[in] iPoint - Point index.
   * \return Value of the under-relaxation parameter for this CV.
   */
  inline su2double GetUnderRelaxation(unsigned long iPoint) const { return UnderRelaxation(iPoint); }

  /*!
   * \brief Set the value of the local CFL number for the current control volume (CV).
   * \param[in] iPoint - Point index.
   * \param[in] val_cfl - the input value of the local CFL number for this CV.
   */
  inline void SetLocalCFL(unsigned long iPoint, su2double val_cfl) { LocalCFL(iPoint) = val_cfl; }

  /*!
   * \brief Get the value of the local CFL number for the current control volume (CV).
   * \param[in] iPoint - Point index.
   * \return Value of the local CFL number for this CV.
   */
  inline su2double GetLocalCFL(unsigned long iPoint) const { return LocalCFL(iPoint); }

  /*!
   * \brief Get the entire Aux matrix of the problem.
   * \return Reference to the aux var  matrix.
   */
  inline const MatrixType& GetAuxVar(void) const { return AuxVar; }

  /*!
   * \brief Get the Aux var value at Point i, variable j.
   */
  inline su2double GetAuxVar(unsigned long iPoint, unsigned long iVar = 0) const { return AuxVar(iPoint,iVar); }

  /*!
   * \brief Set auxiliary variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Varriable indexs
   * \param[in] val_auxvar - Value of the auxiliar variable.
   */
  inline void SetAuxVar(unsigned long iPoint, unsigned long iVar, const su2double auxvar) {
    AuxVar(iPoint,iVar) = auxvar;
  }

  /*!
   * \brief Set value of auxillary gradients.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value of the gradient.
   */
  inline void SetAuxVarGradient(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) {
    Grad_AuxVar(iPoint,iVar,iDim) = value;
  }

  /*!
   * \brief Get the gradient of the auxilary variables.
   * \return Reference to gradient.
   */
  inline CVectorOfMatrix& GetAuxVarGradient(void) { return Grad_AuxVar; }

  /*!
   * \brief Get the value of the auxilliary gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the solution gradient.
   */
  inline su2double GetAuxVarGradient(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const {
    return Grad_AuxVar(iPoint,iVar,iDim);
  }

  /*!
   * \brief Get the value of the auxilliary gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the solution gradient.
   */
  inline CMatrixView<su2double> GetAuxVarGradient(unsigned long iPoint) {
    return Grad_AuxVar[iPoint];
  }

  /*!
   * \brief Add a value to the truncation error.
   * \param[in] iPoint - Point index.
   * \param[in] val_truncation_error - Value that we want to add to the truncation error.
   */
  inline void AddRes_TruncError(unsigned long iPoint, const su2double *val_truncation_error) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Res_TruncError(iPoint, iVar) += val_truncation_error[iVar];
  }

  /*!
   * \brief Subtract a value to the truncation error.
   * \param[in] iPoint - Point index.
   * \param[in] val_truncation_error - Value that we want to subtract to the truncation error.
   */
  inline void SubtractRes_TruncError(unsigned long iPoint, const su2double *val_truncation_error) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Res_TruncError(iPoint, iVar) -= val_truncation_error[iVar];
  }

  /*!
   * \brief Set the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetRes_TruncErrorZero(unsigned long iPoint) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Res_TruncError(iPoint, iVar) = 0.0;
  }

  /*!
   * \brief Set the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVal_ResTruncError_Zero(unsigned long iPoint, unsigned long iVar) {Res_TruncError(iPoint, iVar) = 0.0;}

  /*!
   * \brief Set the momentum part of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetVel_ResTruncError_Zero(unsigned long iPoint) { }

  /*!
   * \brief Set the velocity of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetEnergy_ResTruncError_Zero(unsigned long iPoint) { Res_TruncError(iPoint,nDim+1) = 0.0;}

  /*!
   * \brief Get the truncation error.
   * \param[in] iPoint - Point index.
   * \return Pointer to the truncation error.
   */
  inline su2double *GetResTruncError(unsigned long iPoint) { return Res_TruncError[iPoint]; }

  /*!
   * \brief Get the truncation error.
   * \param[in] iPoint - Point index.
   * \param[in] val_trunc_error - Pointer to the truncation error.
   */
  inline void GetResTruncError(unsigned long iPoint, su2double *val_trunc_error) const {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      val_trunc_error[iVar] = Res_TruncError(iPoint, iVar);
  }

  /*!
   * \brief Set the gradient of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] gradient - Gradient of the solution.
   */
  inline void SetGradient(unsigned long iPoint, su2double **gradient) {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      for (unsigned long iDim = 0; iDim < nDim; iDim++)
        Gradient(iPoint,iVar,iDim) = gradient[iVar][iDim];
  }

  /*!
   * \brief Get the gradient of the entire solution.
   * \return Reference to gradient.
   */
  inline CVectorOfMatrix& GetGradient(void) { return Gradient; }

  /*!
   * \brief Get the value of the solution gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the gradient solution.
   */
  inline CMatrixView<su2double> GetGradient(unsigned long iPoint) { return Gradient[iPoint]; }

  /*!
   * \brief Get the value of the solution gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the solution gradient.
   */
  inline su2double GetGradient(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const { return Gradient(iPoint,iVar,iDim); }

  /*!
   * \brief Add <i>value</i> to the Rmatrix for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] jDim - Index of the dimension.
   * \param[in] value - Value of the Rmatrix entry.
   */
  inline void AddRmatrix(unsigned long iPoint, unsigned long iDim, unsigned long jDim, su2double value) { Rmatrix(iPoint,iDim,jDim) += value; }

  /*!
   * \brief Get the value of the Rmatrix entry for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] jDim - Index of the dimension.
   * \return Value of the Rmatrix entry.
   */
  inline su2double GetRmatrix(unsigned long iPoint, unsigned long iDim, unsigned long jDim) const { return Rmatrix(iPoint,iDim,jDim); }

  /*!
   * \brief Get the value Rmatrix for the entire domain.
   * \return Reference to the Rmatrix.
   */
  inline C3DDoubleMatrix& GetRmatrix(void) { return Rmatrix; }

  /*!
   * \brief Get the slope limiter.
   * \return Reference to the limiters vector.
   */
  inline MatrixType& GetLimiter(void) { return Limiter; }

  /*!
   * \brief Get the value of the slope limiter.
   * \param[in] iPoint - Point index.
   * \return Pointer to the limiters vector.
   */
  inline su2double *GetLimiter(unsigned long iPoint) { return Limiter[iPoint]; }

  /*!
   * \brief Get the value of the slope limiter.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the limiter vector for the variable <i>iVar</i>.
   */
  inline su2double GetLimiter(unsigned long iPoint, unsigned long iVar) const { return Limiter(iPoint,iVar); }

  /*!
   * \brief Get the min solution.
   * \return Value of the min solution for the domain.
   */
  inline MatrixType& GetSolution_Max() { return Solution_Max; }
  inline const MatrixType& GetSolution_Max() const { return Solution_Max; }

  /*!
   * \brief Get the min solution.
   * \return Value of the min solution for the domain.
   */
  inline MatrixType& GetSolution_Min() { return Solution_Min; }
  inline const MatrixType& GetSolution_Min() const { return Solution_Min; }

  /*!
   * \brief Set the value of the time step.
   * \param[in] iPoint - Point index.
   * \param[in] val_delta_time - Value of the time step.
   */
  inline void SetDelta_Time(unsigned long iPoint, su2double val_delta_time) { Delta_Time(iPoint) = val_delta_time; }

  /*!
   * \brief Get the value of the time step.
   * \param[in] iPoint - Point index.
   * \return Value of the time step.
   */
  inline su2double GetDelta_Time(unsigned long iPoint) const {return Delta_Time(iPoint); }

  /*!
   * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline void SetMax_Lambda_Inv(unsigned long iPoint, su2double val_max_lambda) { Max_Lambda_Inv(iPoint) = val_max_lambda; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline void SetMax_Lambda_Visc(unsigned long iPoint, su2double val_max_lambda) { Max_Lambda_Visc(iPoint) = val_max_lambda; }

  /*!
   * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline void AddMax_Lambda_Inv(unsigned long iPoint, su2double val_max_lambda) { Max_Lambda_Inv(iPoint) += val_max_lambda; }

  /*!
   * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline void AddMax_Lambda_Visc(unsigned long iPoint, su2double val_max_lambda) { Max_Lambda_Visc(iPoint) += val_max_lambda; }

  /*!
   * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iPoint - Point index.
   * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline su2double GetMax_Lambda_Inv(unsigned long iPoint) const { return Max_Lambda_Inv(iPoint); }

  /*!
   * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iPoint - Point index.
   * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline su2double GetMax_Lambda_Visc(unsigned long iPoint) const { return Max_Lambda_Visc(iPoint); }

  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline void SetLambda(unsigned long iPoint, su2double val_lambda) { Lambda(iPoint) = val_lambda; }

  /*!
   * \brief Add the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline void AddLambda(unsigned long iPoint, su2double val_lambda) { Lambda(iPoint) += val_lambda; }

  /*!
   * \brief Get the value of the spectral radius.
   * \param[in] iPoint - Point index.
   * \return Value of the spectral radius.
   */
  inline su2double GetLambda(unsigned long iPoint) const { return Lambda(iPoint); }
  inline const VectorType& GetLambda() const { return Lambda; }

  /*!
   * \brief Set pressure sensor.
   * \param[in] iPoint - Point index.
   * \param[in] val_sensor - Value of the pressure sensor.
   */
  inline void SetSensor(unsigned long iPoint, su2double val_sensor) { Sensor(iPoint) = val_sensor; }

  /*!
   * \brief Get the pressure sensor.
   * \param[in] iPoint - Point index.
   * \return Value of the pressure sensor.
   */
  inline su2double GetSensor(unsigned long iPoint) const { return Sensor(iPoint); }
  inline const VectorType& GetSensor() const { return Sensor; }

  /*!
   * \brief Increment the value of the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Variable of the undivided laplacian.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void AddUnd_Lapl(unsigned long iPoint, unsigned long iVar, su2double val_und_lapl) {
    Undivided_Laplacian(iPoint, iVar) += val_und_lapl;
  }

  /*!
   * \brief Set a value to the undivided laplacian.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Variable of the undivided laplacian.
   * \param[in] val_und_lapl - Value of the undivided laplacian.
   */
  inline void SetUnd_Lapl(unsigned long iPoint, unsigned long iVar, su2double val_und_lapl) {
    Undivided_Laplacian(iPoint, iVar) = val_und_lapl;
  }

  /*!
   * \brief Get the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \return Pointer to the undivided laplacian vector.
   */
  inline su2double *GetUndivided_Laplacian(unsigned long iPoint) { return Undivided_Laplacian[iPoint]; }

  /*!
   * \brief Get the undivided laplacian of the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Variable of the undivided laplacian.
   * \return Value of the undivided laplacian vector.
   */
  inline su2double GetUndivided_Laplacian(unsigned long iPoint, unsigned long iVar) const { return Undivided_Laplacian(iPoint, iVar); }
  inline const MatrixType& GetUndivided_Laplacian() const { return Undivided_Laplacian; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow density.
   */
  inline virtual su2double GetDensity(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow density.
   */
  inline virtual su2double GetDensity(unsigned long iPoint, unsigned long val_iSpecies) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline virtual su2double GetMassFraction(unsigned long iPoint, unsigned long val_Species) const { return 0.0; }

  /*!
   * \brief Get the species enthalpy.
   * \return Value of the species enthalpy.
   */
  inline virtual su2double* GetEnthalpys(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow energy.
   */
  inline virtual su2double GetEnergy(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Pointer to the force projection vector.
   */
  inline virtual su2double *GetForceProj_Vector(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Pointer to the objective function source.
   */
  inline virtual su2double *GetObjFuncSource(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the eddy viscosity.
   */
  inline virtual su2double GetEddyViscosity(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow enthalpy.
   */
  inline virtual su2double GetEnthalpy(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the flow pressure.
   */
  inline virtual su2double GetPressure(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline virtual su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the sound speed.
   */
  inline virtual su2double GetSoundSpeed(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the beta for the incompressible flow.
   */
  inline virtual su2double GetBetaInc2(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the temperature.
   */
  inline virtual su2double GetTemperature(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the vibrational-electronic temperature.
   */
  inline virtual su2double GetTemperature_ve(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (trans.-rot.).
   * \param[in] iPoint - Point index.
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  inline virtual su2double GetRhoCv_tr(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (vib.-el.).
   * \param[in] iPoint - Point index.
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  inline virtual su2double GetRhoCv_ve(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline virtual su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the velocity gradient.
   */
  inline virtual CMatrixView<const su2double> GetVelocityGradient(unsigned long iPoint) const {
    return CMatrixView<const su2double>();
  }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Norm 2 of the velocity vector.
   */
  inline virtual su2double GetVelocity2(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return The laminar viscosity of the flow.
   */
  inline virtual su2double GetLaminarViscosity(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the species diffusion coefficient.
   */
  inline virtual su2double* GetDiffusionCoeff(unsigned long iPoint) {return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the thermal conductivity (translational/rotational)
   */
  inline virtual su2double GetThermalConductivity(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the mass diffusivity.
   */
  inline virtual su2double GetDiffusivity(unsigned long iPoint, unsigned short val_ivar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the specific heat at constant P
   */
  inline virtual su2double GetSpecificHeatCp(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the specific heat at constant V
   */
  inline virtual su2double GetSpecificHeatCv(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the thermal conductivity (vibrational)
   */
  inline virtual su2double GetThermalConductivity_ve(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Sets separation intermittency
   */
  inline virtual void SetGammaSep(unsigned long iPoint, su2double gamma_sep) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Sets Effective intermittency
   */
  inline virtual void SetGammaEff(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the vorticity.
   */
  inline virtual su2double *GetVorticity(unsigned long iPoint) { return nullptr; }
  inline virtual const su2double *GetVorticity(unsigned long iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the rate of strain magnitude.
   */
  inline virtual su2double GetStrainMag(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
   */
  inline virtual void SetForceProj_Vector(unsigned long iPoint, const su2double *val_ForceProj_Vector) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
   */
  inline virtual void SetObjFuncSource(unsigned long iPoint, const su2double *val_SetObjFuncSource) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \return Value of the gamma_BC of B-C transition model.
   */
  inline virtual su2double GetGammaBC(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetGammaBC(unsigned long iPoint, su2double val_gamma) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline virtual void SetEddyViscosity(unsigned long iPoint, su2double eddy_visc) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual void SetEnthalpy(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) { return true; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] fluidmodel - fluid model.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, CNEMOGas *fluidmodel) {return false;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, su2double SharpEdge_Distance, bool check, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, su2double Density_Inf, CConfig *config) { return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(unsigned long iPoint, su2double Density_Inf, su2double Viscosity_Inf,
                                 su2double eddy_visc, su2double turb_ke, CConfig *config) {return true; }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline virtual const MatrixType& GetPrimitive() const { AssertOverride(); return Solution; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(unsigned long iPoint, const su2double *val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetPrimitive(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetSecondary(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondary(unsigned long iPoint, unsigned long iVar, su2double val_secondary) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondary(unsigned long iPoint, const su2double *val_secondary) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdPdrho_e(unsigned long iPoint, su2double dPdrho_e) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdPde_rho(unsigned long iPoint, su2double dPde_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdTdrho_e(unsigned long iPoint, su2double dTdrho_e) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdTde_rho(unsigned long iPoint, su2double dTde_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Setdmudrho_T(unsigned long iPoint, su2double dmudrho_T) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdmudT_rho(unsigned long iPoint, su2double dmudT_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Setdktdrho_T(unsigned long iPoint, su2double dktdrho_T) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdktdT_rho(unsigned long iPoint, su2double dktdT_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetSecondary(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetDensity(unsigned long iPoint, su2double val_density) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetDensity(unsigned long iPoint) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPressure(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelocity(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetBetaInc2(unsigned long iPoint, su2double val_betainc2) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_phi - Value of the adjoint velocity.
   */
  inline virtual void SetPhi_Old(unsigned long iPoint, const su2double *val_phi) {}

  /*!
   * \brief A virtual member.
   * \param[in] Gamma - Ratio of Specific heats
   */
  inline virtual bool SetPressure(unsigned long iPoint, su2double Gamma) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config
   */
  inline virtual bool SetPressure(unsigned long iPoint, CConfig *config) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPressure(unsigned long iPoint, su2double Gamma, su2double turb_ke) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdPdU(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdTdU(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdTvedU(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetSoundSpeed(unsigned long iPoint, su2double soundspeed2) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] Gas_Constant - Value of the Gas Constant
   */
  inline virtual bool SetTemperature(unsigned long iPoint, su2double Gas_Constant) { return false; }

  /*!
   * \brief Sets the vibrational electronic temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline virtual bool SetTemperature_ve(unsigned long iPoint, su2double val_Tve) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual bool SetTemperature(unsigned long iPoint, CConfig *config) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual void SetPrimitive(unsigned long iPoint, CConfig *config) {}

  /*!
   * \brief Set the thermal coefficient.
   * \param[in] config - Configuration parameters.
   */
  inline virtual void SetThermalCoeff(unsigned long iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetStress_FEM(unsigned long iPoint, unsigned long iVar, su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void AddStress_FEM(unsigned long iPoint, unsigned long iVar, su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual const su2double *GetStress_FEM(unsigned long iPoint) const {return nullptr;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVonMises_Stress(unsigned long iPoint, su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetVonMises_Stress(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_SurfaceLoad_Res(unsigned long iPoint, const su2double *val_surfForce) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_SurfaceLoad_Res(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_SurfaceLoad_Res() {}

  /*!
   * \brief A virtual member.
   */
  virtual void Set_SurfaceLoad_Res_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_SurfaceLoad_Res_n(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_BodyForces_Res(unsigned long iPoint, const su2double *val_bodyForce) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_BodyForces_Res(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_BodyForces_Res(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_FlowTraction(unsigned long iPoint, const su2double *val_flowTraction) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_FlowTraction(unsigned long iPoint, const su2double *val_flowTraction) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_FlowTraction(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  virtual void Set_FlowTraction_n() {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_FlowTraction_n(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_isVertex(unsigned long iPoint, bool isVertex) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool Get_isVertex(unsigned long iPoint) const { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelocity2(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline virtual void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) {}

  /*!
   * \brief A virtual member.
   * \param[in] laminarViscosity
   */
  inline virtual void SetLaminarViscosity(unsigned long iPoint, su2double laminarViscosity) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetLaminarViscosity(unsigned long iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] thermalConductivity
   */
  inline virtual void SetThermalConductivity(unsigned long iPoint, su2double thermalConductivity) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetThermalConductivity(unsigned long iPoint, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] Cp - Constant pressure specific heat.
   */
  inline virtual void SetSpecificHeatCp(unsigned long iPoint, su2double Cp) {}

  /*!
   * \brief A virtual member.
   * \param[in] Cv - Constant volume specific heat.
   */
  inline virtual void SetSpecificHeatCv(unsigned long iPoint, su2double Cv) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelSolutionDVector(unsigned long iPoint) {}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double GetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief Get the primitive variable gradients for all points.
   * \return Reference to primitive variable gradient.
   */
  inline virtual CVectorOfMatrix& GetGradient_Primitive() { AssertOverride(); return Gradient; }
  inline virtual const CVectorOfMatrix& GetGradient_Primitive() const { AssertOverride(); return Gradient; }

  /*!
   * \brief Get the primitive variables limiter.
   * \return Primitive variables limiter for the entire domain.
   */
  inline virtual MatrixType& GetLimiter_Primitive() { AssertOverride(); return Limiter; }
  inline virtual const MatrixType& GetLimiter_Primitive() const { AssertOverride(); return Limiter; }

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double GetLimiter_Primitive(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the primitive variables gradient.
   */
  inline virtual CMatrixView<su2double> GetGradient_Primitive(unsigned long iPoint, unsigned long iVar=0) { return nullptr; }

  /*!
   * \brief A virtual member.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double *GetLimiter_Primitive(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief Get the value of the primitive gradient for MUSCL reconstruction.
   * \return Value of the primitive gradient for MUSCL reconstruction.
   */
  inline virtual CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline virtual CVectorOfMatrix& GetGradient_Reconstruction() { AssertOverride(); return Gradient; }
  inline virtual const CVectorOfMatrix& GetGradient_Reconstruction() const { AssertOverride(); return Gradient; }

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_density - Value of the density.
   * \param[in] val_dist - Value of the distance to the wall.
   */
  inline virtual void SetBlendingFunc(unsigned long iPoint, su2double val_viscosity, su2double val_dist, su2double val_density, TURB_TRANS_MODEL trans_model) {}

  /*!
   * \brief Get the first blending function of the SST model.
   */
  inline virtual su2double GetF1blending(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Get the second blending function of the SST model.
   */
  inline virtual su2double GetF2blending(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline virtual su2double GetCrossDiff(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Get the value of the eddy viscosity.
   * \return the value of the eddy viscosity.
   */
  inline virtual su2double GetmuT(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Get the value of the intermittency.
   * \return the value of the intermittency.
   */
  inline virtual su2double GetIntermittency(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Set the intermittency.
   * \param[in] val_dist - Value of the  intermittency.
   */
  inline virtual void SetIntermittency(unsigned long iPoint, su2double val_Intermittency) {}

  /*!
   * \brief Get the value of the separation intermittency.
   * \return the value of the separation intermittency.
   */
  inline virtual su2double GetIntermittencySep(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Set the separation intermittency (gamma_sep).
   * \param[in] val_dist - Value of the separation intermittency (gamma_sep).
   */
  inline virtual void SetIntermittencySep(unsigned long iPoint, su2double val_Intermittency_sep) {}

  /*!
   * \brief Get the value of the effective intermittency.
   * \return the value of the effective intermittency.
   */
  inline virtual su2double GetIntermittencyEff(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Set the effective intermittency (gamma_eff).
   * \param[in] Value of the effective intermittency (gamma_eff).
   */
  inline virtual void SetIntermittencyEff(unsigned long iPoint, su2double val_Intermittency_eff) {}

  /*!
   * \brief Set the value of the eddy viscosity.
   * \param[in] val_muT
   */
  inline virtual void SetmuT(unsigned long iPoint, su2double val_muT) {}

  /*!
   * \brief Set the value of the turbulence index.
   * \param[in] val_turb_index - turbulence index
   */
  inline virtual void SetTurbIndex(unsigned long iPoint, su2double val_turb_index) {}

  /*!
   * \brief Get the value of the turbulence index.
   * \return val_turb_index - turbulence index
   */
  inline virtual su2double GetTurbIndex(unsigned long iPoint) const {return 0.0;}

  /*!
   * \brief A virtual member.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_source - Value of the harmonic balance source.
   */
  inline virtual void SetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar, su2double val_source) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief Set the Eddy Viscosity Sensitivity of the problem.
   * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
   * \param[in] numTotalVar - Number of variables.
   */
  inline virtual void SetEddyViscSens(unsigned long iPoint, const su2double *val_EddyViscSens, unsigned long numTotalVar) {}

  /*!
   * \brief Get the Eddy Viscosity Sensitivity of the problem.
   * \return Pointer to the Eddy Viscosity Sensitivity.
   */
  inline virtual su2double *GetEddyViscSens(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Set the direct solution for the adjoint solver.
   * \param[in] solution_direct - Value of the direct solution.
   */
  inline virtual void SetSolution_Direct(unsigned long iPoint, const su2double *solution_direct) {}

  /*!
   * \brief A virtual member. Get the direct solution for the adjoint solver.
   * \return Pointer to the direct solution vector.
   */
  inline virtual su2double *GetSolution_Direct(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Set the restart geometry (coordinate of the converged solution)
   * \param[in] val_coordinate_direct - Value of the restart coordinate.
   */
  inline virtual void SetGeometry_Direct(unsigned long iPoint, const su2double *val_coordinate_direct) {}

  /*!
   * \brief A virtual member. Get the restart geometry (coordinate of the converged solution).
   * \return Pointer to the restart coordinate vector.
   */
  inline virtual su2double *GetGeometry_Direct(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief A virtual member. Get the restart geometry (coordinate of the converged solution).
   * \return Coordinate of the direct solver restart for .
   */
  inline virtual su2double GetGeometry_Direct(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetSolution_Geometry(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] solution - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Geometry(unsigned long iPoint, const su2double *solution_geometry) {}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] solution - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Geometry(unsigned long iPoint, unsigned long iVar, su2double solution_geometry) {}

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_OldSolution_Geometry() {}

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_OldSolution_Geometry(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief Get BGS solution to compute the BGS residual (difference between BGS and BGS_k).
   * \note This is virtual because for some classes the result of a BGS iteration is not "Solution".
   *       If this method is overriden, the BGSSolution_k ones proabably have to be too.
   */
  inline virtual su2double Get_BGSSolution(unsigned long iPoint, unsigned long iVar) const {
    return Solution(iPoint, iVar);
  }

  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  virtual void Set_BGSSolution_k();

  /*!
   * \brief Restore the previous BGS subiteration to solution.
   */
  virtual void Restore_BGSSolution_k();

  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  inline virtual void Set_BGSSolution_k(unsigned long iPoint, unsigned long iVar, su2double val_var) {
    Solution_BGS_k(iPoint,iVar) = val_var;
  }

  /*!
   * \brief Get the value of the solution in the previous BGS subiteration.
   * \param[out] val_solution - solution in the previous BGS subiteration.
   */
  inline virtual su2double Get_BGSSolution_k(unsigned long iPoint, unsigned long iVar) const {
    return Solution_BGS_k(iPoint,iVar);
  }

  /*!
   * \brief Set the value of the velocity (Structural Analysis).
   * \param[in] solution - Solution of the problem (velocity).
   */
  inline virtual void SetSolution_Vel(unsigned long iPoint, const su2double *solution) {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_vel - Value of the solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Vel(unsigned long iPoint, unsigned long iVar, su2double solution_vel) {}

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   * \param[in] solution_vel_time_n - Value of the old solution.
   */
  inline virtual void SetSolution_Vel_time_n(unsigned long iPoint, const su2double *solution_vel_time_n) {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_vel_time_n - Value of the old solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Vel_time_n(unsigned long iPoint, unsigned long iVar, su2double solution_vel_time_n) {}

  /*!
   * \brief Get the solution at time n.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution_time_n(unsigned long iPoint, unsigned long iVar) const { return Solution_time_n(iPoint,iVar); }

  /*!
   * \brief Get the solution at time n-1.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution_time_n1(unsigned long iPoint, unsigned long iVar) const { return Solution_time_n1(iPoint,iVar); }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetSolution_Vel(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \return Pointer to the velocity vector at a point.
   */
  inline virtual su2double* GetSolution_Vel(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Vel_time_n(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \return Pointer to the velocity vector at a point.
   */
  inline virtual su2double* GetSolution_Vel_time_n(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief Set the value of the acceleration (Structural Analysis).
   * \param[in] solution_accel - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Accel(unsigned long iPoint, const su2double *solution_accel) {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_accel - Value of the solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Accel(unsigned long iPoint, unsigned long iVar, su2double solution_accel) {}

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   * \param[in] solution_accel_time_n - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Accel_time_n(unsigned long iPoint, const su2double *solution_accel_time_n) {}

  /*!
   * \overload
   * \param[in] iVar - Index of the variable.
   * \param[in] solution_accel_time_n - Value of the old solution for the index <i>iVar</i>.
   */
  inline virtual void SetSolution_Accel_time_n(unsigned long iPoint, unsigned long iVar, su2double solution_accel_time_n) {}

  /*!
   * \brief Get the acceleration (Structural Analysis).
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline virtual su2double GetSolution_Accel(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline virtual su2double *GetSolution_Accel(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief Get the acceleration of the nodes (Structural Analysis) at time n.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Accel_time_n(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Accel_time_n(unsigned long iPoint) { return nullptr; }

  /*!
   * \brief  A virtual member. Set the value of the velocity solution predictor.
   */
  inline virtual void SetSolution_Vel_Pred(unsigned long iPoint, const su2double *val_solution_pred) { }

  /*!
   * \brief  A virtual member. Set the value of the old solution.
   * \param[in] solution_pred - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred(unsigned long iPoint, const su2double *solution_pred) {}

  /*!
   * \brief  A virtual member. Get the velocity solution predictor.
   * \return Pointer to the velocity solution vector.
   */
  inline virtual const su2double *GetSolution_Vel_Pred(unsigned long iPoint) const {return nullptr; }

  /*!
   * \brief  A virtual member. Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual const su2double *GetSolution_Pred(unsigned long iPoint) const { return nullptr; }

  /*!
   * \brief  A virtual member. Set the value of the old solution.
   * \param[in] solution_pred_Old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred_Old(unsigned long iPoint, const su2double *solution_pred_Old) {}

  /*!
   * \brief  A virtual member. Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual const su2double *GetSolution_Pred_Old(unsigned long iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetReference_Geometry(unsigned long iPoint, unsigned long iVar, su2double ref_geometry) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual const su2double* GetReference_Geometry(unsigned long iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrestretch(unsigned long iPoint, unsigned long iVar, su2double val_prestretch) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual const su2double *GetPrestretch(unsigned long iPoint) const {return nullptr; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrestretch(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the original coordinate iDim.
   */
  inline virtual su2double GetMesh_Coord(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the undeformed coordinates.
   * \return Pointer to the reference coordinates.
   */
  inline virtual const su2double *GetMesh_Coord(unsigned long iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member. Get the undeformed coordinates for the entire domain.
   */
  inline virtual const MatrixType *GetMesh_Coord() const { return nullptr; }

  /*!
   * \brief A virtual member. Set the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \param[in] val_coord - Value of Mesh_Coord[nDim]
   */
  inline virtual void SetMesh_Coord(unsigned long iPoint, unsigned long iDim, su2double val_coord) { }

    /*!
   * \brief A virtual member. Get the value of the wall distance in reference coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the wall distance in reference coordinates.
   */
  inline virtual su2double GetWallDistance(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the value of the wall distance in reference coordinates.
   * \param[in] val_dist - Value of wall distance.
   */
  inline virtual void SetWallDistance(unsigned long iPoint, su2double val_dist) { }

  /*!
   * \brief A virtual member. Register the reference coordinates of the mesh.
   */
  inline virtual void Register_MeshCoord() { }

  /*!
   * \brief A virtual member. Recover the value of the adjoint of the mesh coordinates.
   */
  inline virtual void GetAdjoint_MeshCoord(unsigned long iPoint, su2double *adj_mesh) { }

  /*!
   * \brief A virtual member. Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline virtual su2double GetBound_Disp(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Get the value of the velocity imposed at the boundary.
   * \return Value of the boundary velocity.
   */
  inline virtual su2double GetBound_Vel(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Set the boundary displacement.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline virtual void SetBound_Disp(unsigned long iPoint, const su2double *val_BoundDisp) { }

  /*!
   * \brief A virtual member. Set the boundary velocity.
   * \param[in] val_BoundVel - Pointer to the boundary velocity.
   */
  inline virtual void SetBound_Vel(unsigned long iPoint, const su2double *val_BoundVel) { }

  /*!
   * \brief A virtual member. Set the boundary displacement.
   * \param[in] iDim - Index of the dimension of interest.
   * \param[in] val_BoundDisp - Value of the boundary displacements.
   */
  inline virtual void SetBound_Disp(unsigned long iPoint, unsigned long iDim, const su2double val_BoundDisp) { }

  /*!
   * \brief A virtual member. Set the boundary velocity.
   * \param[in] iDim - Index of the dimension of interest.
   * \param[in] val_BoundVel - Value of the boundary velocity.
   */
  inline virtual void SetBound_Vel(unsigned long iPoint, unsigned long iDim, const su2double val_BoundVel) { }

  /*!
   * \brief A virtual member. Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline virtual const su2double* GetBoundDisp_Direct(unsigned long iPoint) const { return nullptr; }

  /*!
   * \brief A virtual member. Set the solution for the boundary displacements.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline virtual void SetBoundDisp_Direct(unsigned long iPoint, const su2double *val_BoundDisp) { }

  /*!
   * \brief Set the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] val_sens - Pointer to the sensitivities of the boundary displacements.
   */
  inline virtual void SetBoundDisp_Sens(unsigned long iPoint, const su2double *val_sens) { }

  /*!
   * \brief A virtual member. Get the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord_Sens[nDim]
   * \return Value of the original Mesh_Coord_Sens iDim.
   */
  inline virtual su2double GetBoundDisp_Sens(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief A virtual member. Register the boundary displacements of the mesh.
   */
  inline virtual void Register_BoundDisp() { }

  /*!
   * \brief A virtual member. Recover the value of the adjoint of the boundary displacements.
   */
  inline virtual void GetAdjoint_BoundDisp(unsigned long iPoint, su2double *adj_disp) const { }

  /*!
   * \brief A virtual member.
   */
  inline virtual void RegisterFlowTraction(bool reset) { }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double ExtractFlowTractionSensitivity(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief Register the variables in the solution array as input/output variable.
   * \param[in] input - input or output variables.
   */
  void RegisterSolution(bool input);

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
  inline void SetAdjointSolution(unsigned long iPoint, const su2double *adj_sol) {
    for (unsigned long iVar = 0; iVar < AD_OutputIndex.cols(); iVar++)
      AD::SetDerivative(AD_OutputIndex(iPoint,iVar), SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution(unsigned long iPoint, su2double *adj_sol) const {
    for (unsigned long iVar = 0; iVar < AD_InputIndex.cols(); iVar++)
      adj_sol[iVar] = AD::GetDerivative(AD_InputIndex(iPoint,iVar));
  }

  inline void GetAdjointSolution_time_n(unsigned long iPoint, su2double *adj_sol) const {
    for (unsigned long iVar = 0; iVar < Solution_time_n.cols(); iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n(iPoint,iVar));
  }

  inline void GetAdjointSolution_time_n1(unsigned long iPoint, su2double *adj_sol) const {
    for (unsigned long iVar = 0; iVar < Solution_time_n1.cols(); iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n1(iPoint,iVar));
  }

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline virtual void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) {}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline virtual su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const { return 0.0; }
  inline virtual const MatrixType& GetSensitivity() const { AssertOverride(); return Solution; }

  inline virtual void SetTau_Wall(unsigned long iPoint, su2double tau_wall) {}

  inline virtual su2double GetTau_Wall(unsigned long iPoint) const { return 0.0; }

  inline virtual void SetVortex_Tilting(unsigned long iPoint, CMatrixView<const su2double> PrimGrad_Flow,
                                        const su2double* Vorticity, su2double LaminarViscosity) {}

  inline virtual su2double GetVortex_Tilting(unsigned long iPoint) const { return 0.0; }

   /*!
   * \brief A virtual member: Set the recovered pressure for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \param[in] val_pressure - pressure value.
   */
  inline virtual void SetStreamwise_Periodic_RecoveredPressure(unsigned long iPoint,su2double val_pressure) { }

  /*!
   * \brief A virtual member: Get the recovered pressure for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \return Recovered/Physical pressure for streamwise periodic flow.
   */
  inline virtual su2double GetStreamwise_Periodic_RecoveredPressure(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief A virtual member: Set the recovered temperature for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \param[in] val_temperature - temperature value.
   */
  inline virtual void SetStreamwise_Periodic_RecoveredTemperature(unsigned long iPoint, su2double val_temperature) { }

  /*!
   * \brief A virtual member: Get the recovered temperature for streamwise periodic flow.
   * \param[in] iPoint - Point index.
   * \return Recovered/Physical temperature for streamwise periodic flow.
   */
  inline virtual su2double GetStreamwise_Periodic_RecoveredTemperature(unsigned long iPoint) const { return 0.0; }

  /*!
   * \brief Virtual member: Set the Radiative source term at the node
   * \return value of the radiative source term
   */
  inline virtual const su2double *GetRadiative_SourceTerm(unsigned long iPoint) const { return nullptr;}

  /*!
   * \brief  Virtual member: Set the Radiative source term at the node
   * \param[in] val_RadSourceTerm - value of the radiative source term
   */
  inline virtual void SetRadiative_SourceTerm(unsigned long iPoint, unsigned long iVar, su2double val_RadSourceTerm) { }

  /*!
   * \brief Get whether a volumetric heat source is to be introduced in point iPoint
   * \return Bool, determines if this point introduces volumetric heat
   */
  inline virtual bool GetVol_HeatSource(unsigned long iPoint) const { return false; }

  /*!
   * \brief Set the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  virtual void SetFlowTractionSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) { }

  /*!
   * \brief Get the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  virtual su2double GetFlowTractionSensitivity(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief Set the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \param[in] val - value of the source term
   */
  virtual void SetSourceTerm_DispAdjoint(unsigned long iPoint, unsigned long iDim, su2double val) { }
  virtual void SetSourceTerm_VelAdjoint(unsigned long iPoint, unsigned long iDim, su2double val) { }

  /*!
   * \brief Get the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \return value of the source term
   */
  virtual su2double GetSourceTerm_DispAdjoint(unsigned long iPoint, unsigned long iDim) const { return 0.0; }
  virtual su2double GetSourceTerm_VelAdjoint(unsigned long iPoint, unsigned long iDim) const { return 0.0; }

  /*!
   * \brief Set fluid entropy
   * \param[in] iPoint - Node index
   * \param[in] entropy - fluid entropy value.
   */
  inline virtual void SetEntropy(unsigned long iPoint, su2double entropy) { };

  /*!
   * \brief Get fluid entropy
   * \param[in] iPoint - Node index
   * \return Entropy - Fluid entropy value
   */
  inline virtual su2double GetEntropy(unsigned long iPoint) const { return 0; }

  /*!
   * \brief Set dataset extrapolation instance
   * \param[in] iPoint - Node index
   * \param[in] extrapolation - Extrapolation instance (0 = within dataset, 1 = outside dataset)
   */
  inline virtual void SetDataExtrapolation(unsigned long iPoint, unsigned short extrapolation) { };

  /*!
   * \brief Get dataset extrapolation instance
   * \param[in] iPoint - Node index
   * \return extrapolation - Extrapolation instance (0 = within dataset, 1 = outside dataset)
   */
  inline virtual unsigned short GetDataExtrapolation(unsigned long iPoint) const { return 0; }

  /*!
   * \brief Set the number of iterations required by a Newton solver used by the fluid model.
   * \param[in] iPoint - Node index
   * \param[in] nIter - Number of iterations evaluated by the Newton solver
   */
  inline virtual void SetNewtonSolverIterations(unsigned long iPoint, unsigned long nIter) { }

  /*!
   * \brief Get the number of iterations required by a Newton solver used by the fluid model.
   * \param[in] iPoint - Node index
   * \return Number of iterations evaluated by the Newton solver
   */
  inline virtual unsigned long GetNewtonSolverIterations(unsigned long iPoint) const { return 0; }

  /*!
   * \brief LUT premixed flamelet: virtual functions for the speciesflameletvariable LUT
   */
  inline virtual void SetLookupScalar(unsigned long iPoint, su2double val_lookup_scalar, unsigned short val_ivar) { }

  inline virtual void SetScalarSource(unsigned long iPoint, unsigned short val_ivar, su2double val_source) { }

  inline virtual void SetTableMisses(unsigned long iPoint, unsigned short misses) { }

  inline virtual unsigned short GetTableMisses(unsigned long iPoint) const { return 0; }

  inline virtual const su2double *GetScalarSources(unsigned long iPoint) const { return nullptr; }
  inline virtual const su2double *GetScalarLookups(unsigned long iPoint) const { return nullptr; }
};
