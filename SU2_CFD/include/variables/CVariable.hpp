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


using namespace std;

/*!
 * \class CVariable
 * \brief Main class for defining the variables.
 * \author F. Palacios
 */
class CVariable {
protected:

  su2double *Solution,    /*!< \brief Solution of the problem. */
  *Solution_Old;      /*!< \brief Old solution of the problem R-K. */
  bool Non_Physical;      /*!< \brief Non-physical points in the solution (force first order). */
  su2double *Solution_time_n,  /*!< \brief Solution of the problem at time n for dual-time stepping technique. */
  *Solution_time_n1;      /*!< \brief Solution of the problem at time n-1 for dual-time stepping technique. */
  su2double **Gradient;    /*!< \brief Gradient of the solution of the problem. */
  su2double **Rmatrix;    /*!< \brief Geometry-based matrix for weighted least squares gradient calculations. */
  su2double *Limiter;        /*!< \brief Limiter of the solution of the problem. */
  su2double *Solution_Max;    /*!< \brief Max solution for limiter computation. */
  su2double *Solution_Min;    /*!< \brief Min solution for limiter computation. */
  su2double AuxVar;      /*!< \brief Auxiliar variable for gradient computation. */
  su2double *Grad_AuxVar;  /*!< \brief Gradient of the auxiliar variable. */
  su2double Delta_Time;  /*!< \brief Time step. */
  su2double Max_Lambda,  /*!< \brief Maximun eingenvalue. */
  Max_Lambda_Inv,    /*!< \brief Maximun inviscid eingenvalue. */
  Max_Lambda_Visc,  /*!< \brief Maximun viscous eingenvalue. */
  Lambda;        /*!< \brief Value of the eingenvalue. */
  su2double Sensor;  /*!< \brief Pressure sensor for high order central scheme and Roe dissipation. */
  su2double *Undivided_Laplacian;  /*!< \brief Undivided laplacian of the solution. */
  su2double *Res_TruncError,  /*!< \brief Truncation error for multigrid cycle. */
  *Residual_Old,    /*!< \brief Auxiliar structure for residual smoothing. */
  *Residual_Sum;    /*!< \brief Auxiliar structure for residual smoothing. */
  static unsigned short nDim;    /*!< \brief Number of dimension of the problem. */
  unsigned short nVar;    /*!< \brief Number of variables of the problem,
                           note that this variable cannnot be static, it is possible to
                           have different number of nVar in the same problem. */
  unsigned short nPrimVar, nPrimVarGrad;    /*!< \brief Number of variables of the problem,
                                             note that this variable cannnot be static, it is possible to
                                             have different number of nVar in the same problem. */
  unsigned short nSecondaryVar, nSecondaryVarGrad;    /*!< \brief Number of variables of the problem,
                                                       note that this variable cannnot be static, it is possible to
                                                       have different number of nVar in the same problem. */
  su2double *Solution_Adj_Old;    /*!< \brief Solution of the problem in the previous AD-BGS iteration. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CVariable(void);

  /*!
   * \overload
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CVariable(unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CVariable(void);

  /*!
   * \brief Set the value of the solution.
   * \param[in] val_solution - Solution of the problem.
   */
  inline void SetSolution(su2double *val_solution) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = val_solution[iVar];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
   */
  inline void SetSolution(unsigned short val_var, su2double val_solution) {Solution[val_var] = val_solution;}

  /*!
   * \brief Add the value of the solution vector to the previous solution (incremental approach).
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
   */
  inline void Add_DeltaSolution(unsigned short val_var, su2double val_solution) {Solution[val_var] += val_solution;}

  /*!
   * \brief Set the value of the non-physical point.
   * \param[in] val_value - identification of the non-physical point.
   */
  inline void SetNon_Physical(bool val_value) { Non_Physical = !val_value; }

  /*!
   * \brief Get the value of the non-physical point.
   * \return Value of the Non-physical point.
   */
  inline su2double GetNon_Physical(void) { return su2double(Non_Physical); }

  /*!
   * \brief Get the solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution(unsigned short val_var) {return Solution[val_var]; }

  /*!
   * \brief Get the old solution of the problem (Runge-Kutta method)
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Old(unsigned short val_var) {return Solution_Old[val_var]; }

  /*!
   * \brief Get the old solution of the discrete adjoint problem (for multiphysics subiterations=
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Old_Adj(unsigned short val_var) {return Solution_Adj_Old[val_var]; }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Old(su2double *val_solution_old) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_Old[iVar] = val_solution_old[iVar];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_old - Value of the old solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Old(unsigned short val_var, su2double val_solution_old) {Solution_Old[val_var] = val_solution_old; }

  /*!
   * \brief Set old variables to the value of the current variables.
   */
  inline void Set_OldSolution(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_Old[iVar] = Solution[iVar];
  }

  /*!
   * \brief Set variables to the value of the old variables.
   */
  inline void Set_Solution(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
     Solution[iVar] = Solution_Old[iVar];
  }

  /*!
   * \brief Set old discrete adjoint variables to the current value of the adjoint variables.
   */
  inline void Set_OldSolution_Adj(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_Adj_Old[iVar] = Solution[iVar];
  }

  /*!
   * \brief Set the variable solution at time n.
   */
  inline void Set_Solution_time_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_time_n[iVar] = Solution[iVar];
  }

  /*!
   * \brief Set the variable solution at time n-1.
   */
  inline void Set_Solution_time_n1(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_time_n1[iVar] = Solution_time_n[iVar];
  }

  /*!
   * \brief Set the variable solution at time n.
   */
  inline void Set_Solution_time_n(su2double* val_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_time_n[iVar] = val_sol[iVar];
  }

  /*!
   * \brief Set the variable solution at time n-1.
   */
  inline void Set_Solution_time_n1(su2double* val_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_time_n1[iVar] = val_sol[iVar];
  }

  /*!
   * \brief Set to zero the velocity components of the solution.
   */
  inline void SetVelSolutionZero(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = 0.0;
  }

  /*!
   * \brief Specify a vector to set the velocity components of the solution.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionVector(su2double *val_vector) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution[iDim+1] = val_vector[iDim];
  }

  /*!
   * \brief Set to zero velocity components of the solution.
   */
  inline void SetVelSolutionOldZero(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = 0.0;
  }

  /*!
   * \brief Specify a vector to set the velocity components of the old solution.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionOldVector(su2double *val_vector) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Old[iDim+1] = val_vector[iDim];
  }

  /*!
   * \brief Set to zero the solution.
   */
  inline void SetSolutionZero(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
  }

  /*!
   * \brief Set to zero a particular solution.
   */
  inline void SetSolutionZero(unsigned short val_var) {Solution[val_var] = 0.0;}

  /*!
   * \brief Add a value to the solution.
   * \param[in] val_var - Number of the variable.
   * \param[in] val_solution - Value that we want to add to the solution.
   */
  inline void AddSolution(unsigned short val_var, su2double val_solution) {
    Solution[val_var] = Solution_Old[val_var] + val_solution;
  }

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_New(unsigned short val_var) {return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetRoe_Dissipation(void) {return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetRoe_Dissipation(su2double val_dissipation) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetRoe_Dissipation_FD(su2double val_wall_dist) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_delta - A scalar measure of the grid size
   * \param[in] val_const_DES - The DES constant (C_DES)
   */
  inline virtual void SetRoe_Dissipation_NTS(su2double val_delta, su2double val_const_DES) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetDES_LengthScale(void) {return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetDES_LengthScale(su2double val_des_lengthscale) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSolution_New(void) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Number of the variable.
   * \param[in] val_solution - Value that we want to add to the solution.
   */
  inline virtual void AddSolution_New(unsigned short val_var, su2double val_solution) {}

  /*!
   * \brief Add a value to the solution, clipping the values.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the solution change.
   * \param[in] lowerlimit - Lower value.
   * \param[in] upperlimit - Upper value.
   */
  inline void AddClippedSolution(unsigned short val_var, su2double val_solution,
                                 su2double lowerlimit, su2double upperlimit) {

    su2double val_new = Solution_Old[val_var] + val_solution;
    Solution[val_var] = min(max(val_new, lowerlimit), upperlimit);
  }

  /*!
   * \brief Update the variables using a conservative format.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the solution change.
   * \param[in] val_density - Value of the density.
   * \param[in] val_density_old - Value of the old density.
   * \param[in] lowerlimit - Lower value.
   * \param[in] upperlimit - Upper value.
   */
  inline void AddConservativeSolution(unsigned short val_var, su2double val_solution,
                                      su2double val_density, su2double val_density_old,
                                      su2double lowerlimit, su2double upperlimit) {

    su2double val_new = (Solution_Old[val_var]*val_density_old + val_solution)/val_density;
    Solution[val_var] = min(max(val_new, lowerlimit), upperlimit);
  }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline su2double *GetSolution(void) {return Solution; }

  /*!
   * \brief Get the old solution of the problem (Runge-Kutta method)
   * \return Pointer to the old solution vector.
   */
  inline su2double *GetSolution_Old(void) {return Solution_Old; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_time_n(void) {return Solution_time_n; }

  /*!
   * \brief Get the solution at time n-1.
   * \return Pointer to the solution (at time n-1) vector.
   */
  inline su2double *GetSolution_time_n1(void) {return Solution_time_n1; }

  /*!
   * \brief Set the value of the old residual.
   * \param[in] val_residual_old - Pointer to the residual vector.
   */
  inline void SetResidual_Old(su2double *val_residual_old) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Residual_Old[iVar] = val_residual_old[iVar];
  }

  /*!
   * \brief Add a value to the summed residual vector.
   * \param[in] val_residual - Pointer to the residual vector.
   */
  inline void AddResidual_Sum(su2double *val_residual) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Residual_Sum[iVar] += val_residual[iVar];
  }

  /*!
   * \brief Set summed residual vector to zero value.
   */
  inline void SetResidualSumZero(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_Sum[iVar] = 0.0;
  }

  /*!
   * \brief Set the velocity of the truncation error to zero.
   */
  inline virtual void SetVel_ResTruncError_Zero(unsigned short iSpecies) {}

  /*!
   * \brief Get the value of the summed residual.
   * \return Pointer to the summed residual.
   */
  inline su2double *GetResidual_Sum(void) {return Residual_Sum; }

  /*!
   * \brief Get the value of the old residual.
   * \return Pointer to the old residual.
   */
  inline su2double *GetResidual_Old(void) {return Residual_Old; }

  /*!
   * \brief Get the value of the summed residual.
   * \param[in] val_residual - Pointer to the summed residual.
   */
  inline void GetResidual_Sum(su2double *val_residual) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] = Residual_Sum[iVar];
  }

  /*!
   * \brief Set auxiliar variables, we are looking for the gradient of that variable.
   * \param[in] val_auxvar - Value of the auxiliar variable.
   */
  inline void SetAuxVar(su2double val_auxvar) {AuxVar = val_auxvar; }

  /*!
   * \brief Get the value of the auxiliary variable.
   * \return Value of the auxiliary variable.
   */
  inline su2double GetAuxVar(void) {return AuxVar; }

  /*!
   * \brief Set the auxiliary variable gradient to zero value.
   */
  inline void SetAuxVarGradientZero(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Grad_AuxVar[iDim] = 0.0;
  }

  /*!
   * \brief Set the value of the auxiliary variable gradient.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_gradient - Value of the gradient for the index <i>val_dim</i>.
   */
  inline void SetAuxVarGradient(unsigned short val_dim, su2double val_gradient) {Grad_AuxVar[val_dim] = val_gradient;}

  /*!
   * \brief Add a value to the auxiliary variable gradient.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient to be added for the index <i>val_dim</i>.
   */
  inline void AddAuxVarGradient(unsigned short val_dim, su2double val_value) {Grad_AuxVar[val_dim] += val_value;}

  /*!
   * \brief Subtract a value to the auxiliary variable gradient.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient to be subtracted for the index <i>val_dim</i>.
   */
  inline void SubtractAuxVarGradient(unsigned short val_dim, su2double val_value) {Grad_AuxVar[val_dim] -= val_value; }

  /*!
   * \brief Get the gradient of the auxiliary variable.
   * \return Value of the gradient of the auxiliary variable.
   */
  inline su2double *GetAuxVarGradient(void) {return Grad_AuxVar; }

  /*!
   * \brief Get the gradient of the auxiliary variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the gradient of the auxiliary variable for the dimension <i>val_dim</i>.
   */
  inline su2double GetAuxVarGradient(unsigned short val_dim) {return Grad_AuxVar[val_dim]; }

  /*!
   * \brief Add a value to the truncation error.
   * \param[in] val_truncation_error - Value that we want to add to the truncation error.
   */
  inline void AddRes_TruncError(su2double *val_truncation_error) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Res_TruncError[iVar] += val_truncation_error[iVar];
  }

  /*!
   * \brief Subtract a value to the truncation error.
   * \param[in] val_truncation_error - Value that we want to subtract to the truncation error.
   */
  inline void SubtractRes_TruncError(su2double *val_truncation_error) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Res_TruncError[iVar] -= val_truncation_error[iVar];
  }

  /*!
   * \brief Set the truncation error to zero.
   */
  inline void SetRes_TruncErrorZero(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_TruncError[iVar] = 0.0;
  }

  /*!
   * \brief Set the truncation error to zero.
   */
  inline void SetVal_ResTruncError_Zero(unsigned short val_var) {Res_TruncError[val_var] = 0.0;}

  /*!
   * \brief Set the velocity of the truncation error to zero.
   */
  inline void SetVel_ResTruncError_Zero(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Res_TruncError[iDim+1] = 0.0;
  }

  /*!
   * \brief Set the velocity of the truncation error to zero.
   */
  inline void SetEnergy_ResTruncError_Zero(void) {Res_TruncError[nDim+1] = 0.0;}

  /*!
   * \brief Get the truncation error.
   * \return Pointer to the truncation error.
   */
  inline su2double *GetResTruncError(void) {return Res_TruncError; }

  /*!
   * \brief Get the truncation error.
   * \param[in] val_trunc_error - Pointer to the truncation error.
   */
  inline void GetResTruncError(su2double *val_trunc_error) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      val_trunc_error[iVar] = Res_TruncError[iVar];
  }

  /*!
   * \brief Set the gradient of the solution.
   * \param[in] val_gradient - Gradient of the solution.
   */
  inline void SetGradient(su2double **val_gradient) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Gradient[iVar][iDim] = val_gradient[iVar][iDim];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient[val_var][val_dim] = val_value; }

  /*!
   * \brief Set to zero the gradient of the solution.
   */
  inline void SetGradientZero(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Gradient[iVar][iDim] = 0.0;
  }

  /*!
   * \brief Add <i>val_value</i> to the solution gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to add to the solution gradient.
   */
  inline void AddGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient[val_var][val_dim] += val_value; }

  /*!
   * \brief Subtract <i>val_value</i> to the solution gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the solution gradient.
   */
  inline void SubtractGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient[val_var][val_dim] -= val_value; }

  /*!
   * \brief Get the value of the solution gradient.
   * \return Value of the gradient solution.
   */
  inline su2double **GetGradient(void) {return Gradient; }

  /*!
   * \brief Get the value of the solution gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the solution gradient.
   */
  inline su2double GetGradient(unsigned short val_var, unsigned short val_dim) {return Gradient[val_var][val_dim]; }

  /*!
   * \brief Set the value of an entry in the Rmatrix for least squares gradient calculations.
   * \param[in] val_iDim - Index of the dimension.
   * \param[in] val_jDim - Index of the dimension.
   * \param[in] val_value - Value of the Rmatrix entry.
   */
  inline void SetRmatrix(unsigned short val_iDim, unsigned short val_jDim, su2double val_value) {Rmatrix[val_iDim][val_jDim] = val_value; }

  /*!
   * \brief Set to zero the Rmatrix for least squares gradient calculations.
   */
  inline void SetRmatrixZero(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0; jDim < nDim; jDim++)
        Rmatrix[iDim][jDim] = 0.0;
  }

  /*!
   * \brief Add <i>val_value</i> to the Rmatrix for least squares gradient calculations.
   * \param[in] val_iDim - Index of the dimension.
   * \param[in] val_jDim - Index of the dimension.
   * \param[in] val_value - Value to add to the Rmatrix entry.
   */
  inline void AddRmatrix(unsigned short val_iDim, unsigned short val_jDim, su2double val_value) {Rmatrix[val_iDim][val_jDim] += val_value; }

  /*!
   * \brief Get the value of the Rmatrix entry for least squares gradient calculations.
   * \param[in] val_iDim - Index of the dimension.
   * \param[in] val_jDim - Index of the dimension.
   * \return Value of the Rmatrix entry.
   */
  inline su2double GetRmatrix(unsigned short val_iDim, unsigned short val_jDim) {return Rmatrix[val_iDim][val_jDim]; }

  /*!
   * \brief Set the value of the limiter.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_limiter - Value of the limiter for the index <i>val_var</i>.
   */
  inline void SetLimiter(unsigned short val_var, su2double val_limiter) {Limiter[val_var] = val_limiter; }

  /*!
   * \brief Set the value of the limiter.
   * \param[in] val_species - Index of the species .
   * \param[in] val_var - Index of the variable.
   * \param[in] val_limiter - Value of the limiter for the index <i>val_var</i>.
   */
  inline virtual void SetLimiterPrimitive(unsigned short val_species, unsigned short val_var, su2double val_limiter) {}

  /*!
   * \brief Set the value of the limiter.
   * \param[in] val_species - Index of the species .
   * \param[in] val_var - Index of the variable.
   */
  inline virtual su2double GetLimiterPrimitive(unsigned short val_species, unsigned short val_var) {return 0.0; }

  /*!
   * \brief Set the value of the max solution.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the max solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Max(unsigned short val_var, su2double val_solution) {Solution_Max[val_var] = val_solution; }

  /*!
   * \brief Set the value of the min solution.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the min solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Min(unsigned short val_var, su2double val_solution) {Solution_Min[val_var] = val_solution; }

  /*!
   * \brief Get the value of the slope limiter.
   * \return Pointer to the limiters vector.
   */
  inline su2double *GetLimiter(void) {return Limiter; }

  /*!
   * \brief Get the value of the slope limiter.
   * \param[in] val_var - Index of the variable.
   * \return Value of the limiter vector for the variable <i>val_var</i>.
   */
  inline su2double GetLimiter(unsigned short val_var) {return Limiter[val_var]; }

  /*!
   * \brief Get the value of the min solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the min solution for the variable <i>val_var</i>.
   */
  inline su2double GetSolution_Max(unsigned short val_var) {return Solution_Max[val_var]; }

  /*!
   * \brief Get the value of the min solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the min solution for the variable <i>val_var</i>.
   */
  inline su2double GetSolution_Min(unsigned short val_var) {return Solution_Min[val_var]; }

  /*!
   * \brief Get the value of the preconditioner Beta.
   * \return Value of the low Mach preconditioner variable Beta
   */
  inline virtual su2double GetPreconditioner_Beta() {return 0; }

  /*!
   * \brief Set the value of the preconditioner Beta.
   * \param[in] val_Beta - Value of the low Mach preconditioner variable Beta
   */
  inline virtual void SetPreconditioner_Beta(su2double val_Beta) {}

  /*!
   * \brief Get the value of the wind gust
   * \return Value of the wind gust
   */
  inline virtual su2double* GetWindGust() {return 0; }

  /*!
   * \brief Set the value of the wind gust
   * \param[in] val_WindGust - Value of the wind gust
   */
  inline virtual void SetWindGust(su2double* val_WindGust) {}

  /*!
   * \brief Get the value of the derivatives of the wind gust
   * \return Value of the derivatives of the wind gust
   */
  inline virtual su2double* GetWindGustDer() {return NULL;}

  /*!
   * \brief Set the value of the derivatives of the wind gust
   * \param[in] val_WindGust - Value of the derivatives of the wind gust
   */
  inline virtual void SetWindGustDer(su2double* val_WindGust) {}

  /*!
   * \brief Set the value of the time step.
   * \param[in] val_delta_time - Value of the time step.
   */
  inline void SetDelta_Time(su2double val_delta_time) {Delta_Time = val_delta_time; }

  /*!
   * \brief Set the value of the time step.
   * \param[in] val_delta_time - Value of the time step.
   * \param[in] iSpecies - Index of the Species .
   */
  inline virtual void SetDelta_Time(su2double val_delta_time, unsigned short iSpecies) {}

  /*!
   * \brief Get the value of the time step.
   * \return Value of the time step.
   */
  inline su2double GetDelta_Time(void) {return Delta_Time; }

  /*!
   * \brief Get the value of the time step.
   * \param[in] iSpecies - Index of the Species
   * \return Value of the time step.
   */
  inline virtual su2double GetDelta_Time(unsigned short iSpecies) {return 0;}

  /*!
   * \brief Set the value of the maximum eigenvalue.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue.
   */
  inline void SetMax_Lambda(su2double val_max_lambda) {Max_Lambda = val_max_lambda; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline void SetMax_Lambda_Inv(su2double val_max_lambda) {Max_Lambda_Inv = val_max_lambda; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_species - Value of the species index to set the maximum eigenvalue.
   */
  inline virtual void SetMax_Lambda_Inv(su2double val_max_lambda, unsigned short val_species) {}

  /*!
   * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline void SetMax_Lambda_Visc(su2double val_max_lambda) {Max_Lambda_Visc = val_max_lambda; }

  /*!
   * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_species - Index of the species to set the maximum eigenvalue of the viscous terms.
   */
  inline virtual void SetMax_Lambda_Visc(su2double val_max_lambda, unsigned short val_species) {}

  /*!
   * \brief Add a value to the maximum eigenvalue.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue.
   */
  inline void AddMax_Lambda(su2double val_max_lambda) {Max_Lambda += val_max_lambda; }

  /*!
   * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline void AddMax_Lambda_Inv(su2double val_max_lambda) {Max_Lambda_Inv += val_max_lambda; }

  /*!
   * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline void AddMax_Lambda_Visc(su2double val_max_lambda) {Max_Lambda_Visc += val_max_lambda; }

  /*!
   * \brief Get the value of the maximum eigenvalue.
   * \return the value of the maximum eigenvalue.
   */
  inline su2double GetMax_Lambda(void) {return Max_Lambda; }

  /*!
   * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
   */
  inline su2double GetMax_Lambda_Inv(void) {return Max_Lambda_Inv; }

  /*!
   * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
   * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
   */
  inline su2double GetMax_Lambda_Visc(void) {return Max_Lambda_Visc; }

  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline void SetLambda(su2double val_lambda) {Lambda = val_lambda; }

  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] val_lambda - Value of the spectral radius.
   * \param[in] val_iSpecies -Index of species
   */
  inline virtual void SetLambda(su2double val_lambda, unsigned short val_iSpecies) {}

  /*!
   * \brief Add the value of the spectral radius.
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline void AddLambda(su2double val_lambda) {Lambda += val_lambda; }

  /*!
   * \brief Add the value of the spectral radius.
   * \param[in] val_iSpecies -Index of species
   * \param[in] val_lambda - Value of the spectral radius.
   */
  inline virtual void AddLambda(su2double val_lambda, unsigned short val_iSpecies) {}

  /*!
   * \brief Get the value of the spectral radius.
   * \return Value of the spectral radius.
   */
  inline su2double GetLambda(void) {return Lambda; }

  /*!
   * \brief Get the value of the spectral radius.
   * \param[in] val_iSpecies -Index of species
   * \return Value of the spectral radius.
   */
  inline virtual su2double GetLambda(unsigned short val_iSpecies) {return 0.0;}

  /*!
   * \brief Set pressure sensor.
   * \param[in] val_sensor - Value of the pressure sensor.
   */
  inline void SetSensor(su2double val_sensor) {Sensor = val_sensor; }

  /*!
   * \brief Set pressure sensor.
   * \param[in] val_sensor - Value of the pressure sensor.
   * \param[in] iSpecies - Index of the species.
   */
  inline virtual void SetSensor(su2double val_sensor, unsigned short iSpecies) {}

  /*!
   * \brief Get the pressure sensor.
   * \return Value of the pressure sensor.
   */
  inline su2double GetSensor(void) {return Sensor; }

  /*!
   * \brief Get the pressure sensor.
   * \param[in] iSpecies - index of species
   * \return Value of the pressure sensor.
   */
  inline virtual su2double GetSensor(unsigned short iSpecies) {return 0;}

  /*!
   * \brief Set the value of the undivided laplacian of the solution.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_undivided_laplacian - Value of the undivided solution for the index <i>val_var</i>.
   */
  inline void SetUndivided_Laplacian(unsigned short val_var, su2double val_undivided_laplacian) {
    Undivided_Laplacian[val_var] = val_undivided_laplacian;
  }

  /*!
   * \brief Add the value of the undivided laplacian of the solution.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void AddUnd_Lapl(su2double *val_und_lapl) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Undivided_Laplacian[iVar] += val_und_lapl[iVar];
  }

  /*!
   * \brief Subtract the value of the undivided laplacian of the solution.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void SubtractUnd_Lapl(su2double *val_und_lapl) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Undivided_Laplacian[iVar] -= val_und_lapl[iVar];
  }

  /*!
   * \brief Subtract the value of the undivided laplacian of the solution.
   * \param[in] val_var - Variable of the undivided laplacian.
   * \param[in] val_und_lapl - Value of the undivided solution.
   */
  inline void SubtractUnd_Lapl(unsigned short val_var, su2double val_und_lapl) {
    Undivided_Laplacian[val_var] -= val_und_lapl;
  }

  /*!
   * \brief Set the undivided laplacian of the solution to zero.
   */
  inline void SetUnd_LaplZero(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Undivided_Laplacian[iVar] = 0.0;
  }

  /*!
   * \brief Set a value to the undivided laplacian.
   * \param[in] val_var - Variable of the undivided laplacian.
   * \param[in] val_und_lapl - Value of the undivided laplacian.
   */
  inline void SetUnd_Lapl(unsigned short val_var, su2double val_und_lapl) {
    Undivided_Laplacian[val_var] = val_und_lapl;
  }

  /*!
   * \brief Get the undivided laplacian of the solution.
   * \return Pointer to the undivided laplacian vector.
   */
  inline su2double *GetUndivided_Laplacian(void) {return Undivided_Laplacian; }

  /*!
   * \brief Get the undivided laplacian of the solution.
   * \param[in] val_var - Variable of the undivided laplacian.
   * \return Value of the undivided laplacian vector.
   */
  inline su2double GetUndivided_Laplacian(unsigned short val_var) {return Undivided_Laplacian[val_var]; }

  /*!
   * \brief A virtual member.
   * \return Value of the flow density.
   */
  inline virtual su2double GetDensity(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Old value of the flow density.
   */
  inline virtual su2double GetDensity_Old(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the flow density.
   */
  inline virtual su2double GetDensity(unsigned short val_iSpecies) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline virtual su2double GetMassFraction(unsigned short val_Species) {return 0.0;}

  /*!
   * \brief A virtual member.
   * \return Value of the flow energy.
   */
  inline virtual su2double GetEnergy(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Pointer to the force projection vector.
   */
  inline virtual su2double *GetForceProj_Vector(void) {return NULL; }

  /*!
   * \brief A virtual member.
   * \return Pointer to the objective function source.
   */
  inline virtual su2double *GetObjFuncSource(void) {return NULL; }

  /*!
   * \brief A virtual member.
   * \return Pointer to the internal boundary vector.
   */
  inline virtual su2double *GetIntBoundary_Jump(void) {return NULL; }

  /*!
   * \brief A virtual member.
   * \return Value of the eddy viscosity.
   */
  inline virtual su2double GetEddyViscosity(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the flow enthalpy.
   */
  inline virtual su2double GetEnthalpy(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the flow pressure.
   */
  inline virtual su2double GetPressure(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline virtual su2double GetProjVel(su2double *val_vector) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_vector - Direction of projection.
   * \param[in] val_species - Index of the desired species.
   * \return Value of the projected velocity.
   */
  inline virtual su2double GetProjVel(su2double *val_vector, unsigned short val_species) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sound speed.
   */
  inline virtual su2double GetSoundSpeed(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the beta for the incompressible flow.
   */
  inline virtual su2double GetBetaInc2(void) { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the temperature.
   */
  inline virtual su2double GetTemperature(void) {return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the vibrational-electronic temperature.
   */
  inline virtual su2double GetTemperature_ve(void) {return 0; }

  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (trans.-rot.).
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  inline virtual su2double GetRhoCv_tr(void) {return 0; }

  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (vib.-el.).
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  inline virtual su2double GetRhoCv_ve(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>val_dim</i>.
   */
  inline virtual su2double GetVelocity(unsigned short val_dim) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Norm 2 of the velocity vector.
   */
  inline virtual su2double GetVelocity2(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Norm 2 of the velocity vector of Fluid val_species.
   */
  inline virtual su2double GetVelocity2(unsigned short val_species) {return 0;}

  /*!
   * \brief A virtual member.
   * \return The laminar viscosity of the flow.
   */
  inline virtual su2double GetLaminarViscosity(void) {return 0; }


  /*!
   * \brief A virtual member.
   * \return The laminar viscosity of the flow.
   */
  inline virtual su2double GetLaminarViscosity(unsigned short iSpecies) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the species diffusion coefficient.
   */
  inline virtual su2double* GetDiffusionCoeff(void) {return NULL; }

  /*!
   * \brief A virtual member.
   * \return Value of the thermal conductivity (translational/rotational)
   */
  inline virtual su2double GetThermalConductivity(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the specific heat at constant P
   */
  inline virtual su2double GetSpecificHeatCp(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the specific heat at constant V
   */
  inline virtual su2double GetSpecificHeatCv(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the thermal conductivity (vibrational)
   */
  inline virtual su2double GetThermalConductivity_ve(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Sets separation intermittency
   */
  inline virtual void SetGammaSep(su2double gamma_sep) {}

  /*!
   * \brief A virtual member.
   * \return Sets separation intermittency
   */
  inline virtual void SetGammaEff(void) {}

  /*!
   * \brief A virtual member.
   * \return Returns intermittency
   */
  inline virtual su2double GetIntermittency() { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the vorticity.
   */
  inline virtual su2double *GetVorticity(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the rate of strain magnitude.
   */
  inline virtual su2double GetStrainMag(void) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
   */
  inline virtual void SetForceProj_Vector(su2double *val_ForceProj_Vector) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
   */
  inline virtual void SetObjFuncSource(su2double *val_SetObjFuncSource) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump.
   */
  inline virtual void SetIntBoundary_Jump(su2double *val_IntBoundary_Jump) {}

  /*!
   * \brief A virtual member.
   * \return Value of the gamma_BC of B-C transition model.
   */
  inline virtual su2double GetGammaBC(void) {return 0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetGammaBC(su2double val_gamma) {}

  /*!
   * \brief A virtual member.
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline virtual void SetEddyViscosity(su2double eddy_visc) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetEnthalpy(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(CConfig *config) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(CFluidModel *FluidModel) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondaryVar(CFluidModel *FluidModel) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool Cons2PrimVar(CConfig *config, su2double *U, su2double *V, su2double *dPdU,
                                   su2double *dTdU, su2double *dTvedU) { return false; }
  /*!
   * \brief A virtual member.
   */
  inline virtual void Prim2ConsVar(CConfig *config, su2double *V, su2double *U) {return; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(su2double SharpEdge_Distance, bool check, CConfig *config) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(su2double eddy_visc, su2double turb_ke, CConfig *config) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(su2double Density_Inf, CConfig *config) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPrimVar(su2double Density_Inf, su2double Viscosity_Inf, su2double eddy_visc, su2double turb_ke, CConfig *config) {return true; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrimitive(unsigned short val_var) {return 0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(unsigned short val_var, su2double val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(su2double *val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetPrimitive(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetSecondary(unsigned short val_var) {return 0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondary(unsigned short val_var, su2double val_secondary) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetSecondary(su2double *val_secondary) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdPdrho_e(su2double dPdrho_e) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdPde_rho(su2double dPde_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdTdrho_e(su2double dTdrho_e) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdTde_rho(su2double dTde_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Setdmudrho_T(su2double dmudrho_T) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdmudT_rho(su2double dmudT_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Setdktdrho_T(su2double dktdrho_T) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetdktdT_rho(su2double dktdT_rho) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetSecondary(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetDensity(su2double val_density) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetDensity(void) { return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPressure(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelocity(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetBetaInc2(su2double val_betainc2) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_phi - Value of the adjoint velocity.
   */
  inline virtual void SetPhi_Old(su2double *val_phi) {}

  /*!
   * \brief A virtual member.
   * \param[in] Gamma - Ratio of Specific heats
   */
  inline virtual bool SetPressure(su2double Gamma) {return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config
   */
  inline virtual bool SetPressure(CConfig *config) {return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetPressure(su2double Gamma, su2double turb_ke) {return false; }

  /*!
   * \brief Calculates vib.-el. energy per mass, \f$e^{vib-el}_s\f$, for input species (not including KE)
   */
  inline virtual su2double CalcEve(su2double *V, CConfig *config, unsigned short val_Species) {return 0; }

  /*!
   * \brief Calculates enthalpy per mass, \f$h_s\f$, for input species (not including KE)
   */
  inline virtual su2double CalcHs(su2double *V, CConfig *config, unsigned short val_Species) {return 0; }

  /*!
   * \brief Calculates enthalpy per mass, \f$Cv_s\f$, for input species (not including KE)
   */
  inline virtual su2double CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dPdU
   */
  inline virtual void CalcdPdU(su2double *V, CConfig *config, su2double *dPdU) {}

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dTdU
   */
  inline virtual void CalcdTdU(su2double *V, CConfig *config, su2double *dTdU) {}

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dTdU
   */
  inline virtual void CalcdTvedU(su2double *V, CConfig *config, su2double *dTdU) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdPdU(void) { return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdTdU(void) { return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetdTvedU(void) { return NULL; }

  /*!
   * \brief A virtual member.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] Gamma - Ratio of Specific heats
   */
  inline virtual void SetDeltaPressure(su2double *val_velocity, su2double Gamma) {}

  /*!
   * \brief A virtual member.
   * \param[in] Gamma - Ratio of specific heats.
   */
  inline virtual bool SetSoundSpeed(su2double Gamma) {return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual bool SetSoundSpeed(CConfig *config) {return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetSoundSpeed(void) { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] Gas_Constant - Value of the Gas Constant
   */
  inline virtual bool SetTemperature(su2double Gas_Constant) {return false; }

  /*!
   * \brief Sets the vibrational electronic temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline virtual bool SetTemperature_ve(su2double val_Tve) {return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual bool SetTemperature(CConfig *config) {return false; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  inline virtual void SetPrimitive(CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   * \param[in] Coord - Physical coordinates.
   */
  inline virtual void SetPrimitive(CConfig *config, su2double *Coord) {}

  /*!
   * \brief A virtual member.
   * \param[in] Temperature_Wall - Value of the Temperature at the wall
   */
  inline virtual void SetWallTemperature(su2double Temperature_Wall) {}

  /*!
   * \brief A virtual member.
   * \param[in] Temperature_Wall - Value of the Temperature at the wall
   */
  inline virtual void SetWallTemperature(su2double* Temperature_Wall) {}

  /*!
   * \brief Set the thermal coefficient.
   * \param[in] config - Configuration parameters.
   */
  inline virtual void SetThermalCoeff(CConfig *config) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetStress_FEM(unsigned short iVar, su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void AddStress_FEM(unsigned short iVar, su2double val_stress) {}

  /*!
   * \brief A virtual member.

   */
  inline virtual su2double *GetStress_FEM(void) {return NULL;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVonMises_Stress(su2double val_stress) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetVonMises_Stress(void) {return 0.0;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_SurfaceLoad_Res(su2double *val_surfForce) {}

  /*!
   * \brief  A virtual member.
   */
  inline virtual void Set_SurfaceLoad_Res(unsigned short iVar, su2double val_surfForce) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_SurfaceLoad_Res(unsigned short iVar) {return 0.0;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_SurfaceLoad_Res(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_SurfaceLoad_Res_n(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_SurfaceLoad_Res_n(unsigned short iVar) {return 0.0;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_BodyForces_Res(su2double *val_bodyForce) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_BodyForces_Res(unsigned short iVar) {return 0.0;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_BodyForces_Res(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_FlowTraction(su2double *val_flowTraction) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Add_FlowTraction(su2double *val_flowTraction) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_FlowTraction(unsigned short iVar) {return 0.0;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_FlowTraction_n(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double Get_FlowTraction_n(unsigned short iVar) {return 0.0;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Clear_FlowTraction(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool Get_isVertex(void) {return false;}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelocity2(void) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline virtual void SetVelocity_Old(su2double *val_velocity) {}

  /*!
   * \brief A virtual member.
   * \param[in] laminarViscosity
   */
  inline virtual void SetLaminarViscosity(su2double laminarViscosity) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetLaminarViscosity(CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] thermalConductivity
   */
  inline virtual void SetThermalConductivity(su2double thermalConductivity) {}

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetThermalConductivity(CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] Cp - Constant pressure specific heat.
   */
  inline virtual void SetSpecificHeatCp(su2double Cp) {}

  /*!
   * \brief A virtual member.
   * \param[in] Cv - Constant volume specific heat.
   */
  inline virtual void SetSpecificHeatCv(su2double Cv) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetVorticity(void) {return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual bool SetStrainMag(void) {return false; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelSolutionOldDVector(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetVelSolutionDVector(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetGradient_PrimitiveZero(unsigned short val_primvar) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to add to the gradient of the primitive variables.
   */
  inline virtual void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
   */
  inline virtual void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double GetLimiter_Primitive(unsigned short val_var) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline virtual void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_value - Value of the gradient.
   */
  inline virtual void SetLimiter_Primitive(unsigned short val_var, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double **GetGradient_Primitive(void) {return NULL; }

  /*!
   * \brief A virtual member.
   * \return Value of the primitive variables gradient.
   */
  inline virtual su2double *GetLimiter_Primitive(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetGradient_SecondaryZero(unsigned short val_secondaryvar) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to add to the gradient of the Secondary variables.
   */
  inline virtual void AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the gradient of the Secondary variables.
   */
  inline virtual void SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the Secondary variables gradient.
   */
  inline virtual su2double GetGradient_Secondary(unsigned short val_var, unsigned short val_dim) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \return Value of the Secondary variables gradient.
   */
  inline virtual su2double GetLimiter_Secondary(unsigned short val_var) {return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline virtual void SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_value - Value of the gradient.
   */
  inline virtual void SetLimiter_Secondary(unsigned short val_var, su2double val_value) {}

  /*!
   * \brief A virtual member.
   * \return Value of the Secondary variables gradient.
   */
  inline virtual su2double **GetGradient_Secondary(void) {return NULL; }

  /*!
   * \brief A virtual member.
   * \return Value of the Secondary variables gradient.
   */
  inline virtual su2double *GetLimiter_Secondary(void) {return NULL; }

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_density - Value of the density.
   * \param[in] val_dist - Value of the distance to the wall.
   */
  inline virtual void SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) {}

  /*!
   * \brief Get the first blending function of the SST model.
   */
  inline virtual su2double GetF1blending(void) {return 0; }

  /*!
   * \brief Get the second blending function of the SST model.
   */
  inline virtual su2double GetF2blending(void) {return 0; }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline virtual su2double GetCrossDiff(void) { return 0.0; }

  /*!
   * \brief Get the value of the eddy viscosity.
   * \return the value of the eddy viscosity.
   */
  inline virtual su2double GetmuT(void) { return 0.0; }

  /*!
   * \brief Set the value of the eddy viscosity.
   * \param[in] val_muT
   */
  inline virtual void SetmuT(su2double val_muT) {}

  /*!
   * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
   * \param[in] iSpecies - Value of iSpecies to which the eigenvalue belongs
   */
  inline virtual void AddMax_Lambda_Inv(su2double val_max_lambda, unsigned short iSpecies) {}

  /*!
   * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
   * \param[in] iSpecies - Value of iSpecies to which the eigenvalue belongs
   */
  inline virtual void AddMax_Lambda_Visc(su2double val_max_lambda, unsigned short iSpecies) {}

  /*!
   * \brief A virtual member.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_source - Value of the harmonic balance source.
   */
  inline virtual void SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetHarmonicBalance_Source(unsigned short val_var) {return 0; }

  /*!
   * \brief Set the Eddy Viscosity Sensitivity of the problem.
   * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
   * \param[in] numTotalVar - Number of variables.
   */
  inline virtual void SetEddyViscSens(su2double *val_EddyViscSens, unsigned short numTotalVar) {}

  /*!
   * \brief Get the Eddy Viscosity Sensitivity of the problem.
   * \return Pointer to the Eddy Viscosity Sensitivity.
   */
  inline virtual su2double *GetEddyViscSens(void) {return NULL; }

  /*!
   * \brief A virtual member. Set the direct solution for the adjoint solver.
   * \param[in] val_solution_direct - Value of the direct solution.
   */
  inline virtual void SetSolution_Direct(su2double *val_solution_direct) {}

  /*!
   * \brief A virtual member. Get the direct solution for the adjoint solver.
   * \return Pointer to the direct solution vector.
   */
  inline virtual su2double *GetSolution_Direct(void) { return NULL; }

  /*!
   * \brief A virtual member. Set the restart geometry (coordinate of the converged solution)
   * \param[in] val_coordinate_direct - Value of the restart coordinate.
   */
  inline virtual void SetGeometry_Direct(su2double *val_coordinate_direct) {}

  /*!
   * \brief A virtual member. Get the restart geometry (coordinate of the converged solution).
   * \return Pointer to the restart coordinate vector.
   */
  inline virtual su2double *GetGeometry_Direct(void) { return NULL; }

  /*!
   * \brief A virtual member. Get the restart geometry (coordinate of the converged solution).
   * \return Coordinate of the direct solver restart for .
   */
  inline virtual su2double GetGeometry_Direct(unsigned short val_dim) {return 0.0; }

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline virtual su2double GetSolution_Geometry(unsigned short val_var) {return 0.0;}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Geometry(su2double *val_solution_geometry) {}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Geometry(unsigned short val_var, su2double val_solution_geometry) {}

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline virtual su2double GetGeometry_CrossTerm_Derivative(unsigned short val_var) {return 0.0;}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline virtual void SetGeometry_CrossTerm_Derivative(unsigned short iDim, su2double der) {}

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline virtual su2double GetGeometry_CrossTerm_Derivative_Flow(unsigned short val_var) {return 0.0;}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline virtual void SetGeometry_CrossTerm_Derivative_Flow(unsigned short iDim, su2double der) {}

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_OldSolution_Geometry(void) {}

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_OldSolution_Geometry(unsigned short iDim) {return 0.0;}

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_BGSSolution(unsigned short iDim, su2double val_solution) {}

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_BGSSolution_k(void) {}

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_BGSSolution(unsigned short iDim) {return 0.0;}

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_BGSSolution_k(unsigned short iDim) {return 0.0;}

  /*!
   * \brief A virtual member. Set the value of the old geometry solution (adjoint).
   */
  inline virtual void Set_BGSSolution_Geometry(void) {}

  /*!
   * \brief A virtual member. Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline virtual su2double Get_BGSSolution_Geometry(unsigned short iDim) {return 0.0;}

  /*!
   * \brief  A virtual member. Set the contribution of crossed terms into the derivative.
   */
  inline virtual void SetCross_Term_Derivative(unsigned short iVar, su2double der) {}

  /*!
   * \brief  A virtual member. Get the contribution of crossed terms into the derivative.
   * \return The contribution of crossed terms into the derivative.
   */
  inline virtual su2double GetCross_Term_Derivative(unsigned short iVar) {return 0.0; }

  /*!
   * \brief A virtual member. Set the direct velocity solution for the adjoint solver.
   * \param[in] val_solution_direct - Value of the direct velocity solution.
   */
  inline virtual void SetSolution_Vel_Direct(su2double *sol) {}

  /*!
   * \brief A virtual member. Set the direct acceleration solution for the adjoint solver.
   * \param[in] val_solution_direct - Value of the direct acceleration solution.
   */
  inline virtual void SetSolution_Accel_Direct(su2double *sol) {}

  /*!
   * \brief A virtual member. Get the direct velocity solution for the adjoint solver.
   * \return Pointer to the direct velocity solution vector.
   */
  inline virtual su2double* GetSolution_Vel_Direct() {return NULL; }

  /*!
   * \brief A virtual member. Get the direct acceleraction solution for the adjoint solver.
   * \return Pointer to the direct acceleraction solution vector.
   */
  inline virtual su2double* GetSolution_Accel_Direct() {return NULL; }

  /*!
   * \brief Set the value of the old solution.
   */
  inline virtual void SetSolution_time_n(void) {}

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_time_n - Pointer to the residual vector.
   */
  inline virtual void SetSolution_time_n(unsigned short val_var, su2double val_solution) {}

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_time_n(su2double *val_solution_time_n) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_time_n[iVar] = val_solution_time_n[iVar];
  }

  /*!
   * \brief Set the value of the velocity (Structural Analysis).
   * \param[in] val_solution - Solution of the problem (velocity).
   */
  inline virtual void SetSolution_Vel(su2double *val_solution) {}

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_vel - Value of the solution for the index <i>val_var</i>.
   */
  inline virtual void SetSolution_Vel(unsigned short val_var, su2double val_solution_vel) {}

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   * \param[in] val_solution_vel_time_n - Value of the old solution.
   */
  inline virtual void SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) {}

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   */
  inline virtual void SetSolution_Vel_time_n(void) {}

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_vel_time_n - Value of the old solution for the index <i>val_var</i>.
   */
  inline virtual void SetSolution_Vel_time_n(unsigned short val_var, su2double val_solution_vel_time_n) {}

  /*!
   * \brief Get the solution at time n.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution_time_n(unsigned short val_var) {return Solution_time_n[val_var]; }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline virtual su2double GetSolution_Vel(unsigned short val_var) {return 0; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline virtual su2double *GetSolution_Vel(void) {return NULL; }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Vel_time_n(unsigned short val_var) {return 0; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Vel_time_n(void) {return NULL; }


  /*!
   * \brief Set the value of the acceleration (Structural Analysis).
   * \param[in] val_solution_accel - Solution of the problem (acceleration).
   */
  inline virtual void SetSolution_Accel(su2double *val_solution_accel) {}

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_accel - Value of the solution for the index <i>val_var</i>.
   */
  inline virtual void SetSolution_Accel(unsigned short val_var, su2double val_solution_accel) {}

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   * \param[in] val_solution_accel_time_n - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) {}

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   */
  inline virtual void SetSolution_Accel_time_n(void) {}

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_accel_time_n - Value of the old solution for the index <i>val_var</i>.
   */
  inline virtual void SetSolution_Accel_time_n(unsigned short val_var, su2double val_solution_accel_time_n) {}

  /*!
   * \brief Get the acceleration (Structural Analysis).
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline virtual su2double GetSolution_Accel(unsigned short val_var) {return 0; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline virtual su2double *GetSolution_Accel(void) {return NULL; }

  /*!
   * \brief Get the acceleration of the nodes (Structural Analysis) at time n.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Accel_time_n(unsigned short val_var) {return 0; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Accel_time_n(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_OldSolution_Vel(void) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void Set_OldSolution_Accel(void) {}

  /*!
   * \brief  A virtual member. Set the value of the solution predictor.
   */
  inline virtual void SetSolution_Pred(void) {}

  /*!
   * \brief  A virtual member. Set the value of the old solution.
   * \param[in] val_solution_pred - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred(su2double *val_solution_pred) {}

  /*!
   * \brief  A virtual member. Set the value of the solution predicted.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred(unsigned short val_var, su2double val_solution_pred) {}

  /*!
   * \brief  A virtual member. Get the value of the solution predictor.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Pred(unsigned short val_var) {return 0.0; }

  /*!
   * \brief  A virtual member. Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Pred(void) {return NULL; }

  /*!
   * \brief  A virtual member. Set the value of the solution predictor.
   */
  inline virtual void SetSolution_Pred_Old(void) {}

  /*!
   * \brief  A virtual member. Set the value of the old solution.
   * \param[in] val_solution_pred_Old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred_Old(su2double *val_solution_pred_Old) {}

  /*!
   * \brief  A virtual member. Set the value of the old solution predicted.
   * \param[in] val_solution_pred_old - Pointer to the residual vector.
   */
  inline virtual void SetSolution_Pred_Old(unsigned short val_var, su2double val_solution_pred_old) {}

  /*!
   * \brief  A virtual member. Get the value of the solution predictor.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline virtual su2double GetSolution_Pred_Old(unsigned short val_var) {return 0.0; }

  /*!
   * \brief  A virtual member. Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline virtual su2double *GetSolution_Pred_Old(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetReference_Geometry(unsigned short iVar, su2double ref_geometry) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetReference_Geometry(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrestretch(unsigned short iVar, su2double val_prestretch) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetPrestretch(void) {return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrestretch(unsigned short iVar) {return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetReference_Geometry(unsigned short iVar) {return 0.0; }

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
  inline virtual void SetAdjointSolution_Vel(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Vel(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Vel_time_n(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Vel_time_n(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Accel(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Accel(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetAdjointSolution_Accel_time_n(su2double *adj_sol) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void GetAdjointSolution_Accel_time_n(su2double *adj_sol) {}

  /*!
   * \brief Register the variables in the solution array as input/output variable.
   * \param[in] input - input or output variables.
   */
  inline void RegisterSolution(bool input) {
    if (input) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        AD::RegisterInput(Solution[iVar]);
    }
    else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
        AD::RegisterOutput(Solution[iVar]);}
  }

  /*!
   * \brief Register the variables in the solution_time_n array as input/output variable.
   */
  inline void RegisterSolution_time_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution_time_n[iVar]);
  }

  /*!
   * \brief Register the variables in the solution_time_n1 array as input/output variable.
   */
  inline void RegisterSolution_time_n1(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution_time_n1[iVar]);
  }

  /*!
   * \brief Set the adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution.
   * \param[out] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution[iVar]);
  }

  /*!
   * \brief Set the adjoint values of the solution at time n.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_time_n(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution at time n.
   * \param[out] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_time_n(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n[iVar]);
  }

  /*!
   * \brief Set the adjoint values of the solution at time n-1.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_time_n1(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n1[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the adjoint values of the solution at time n-1.
   * \param[out] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_time_n1(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n1[iVar]);
  }

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline virtual void SetSensitivity(unsigned short iDim, su2double val) {}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline virtual su2double GetSensitivity(unsigned short iDim) {return 0.0; }

  inline virtual void SetDual_Time_Derivative(unsigned short iVar, su2double der) {}

  inline virtual void SetDual_Time_Derivative_n(unsigned short iVar, su2double der) {}

  inline virtual su2double GetDual_Time_Derivative(unsigned short iVar) {return 0.0;}

  inline virtual su2double GetDual_Time_Derivative_n(unsigned short iVar) {return 0.0;}

  inline virtual void SetTauWall(su2double val_tau_wall) {}

  inline virtual su2double GetTauWall() {return 0.0; }

  inline virtual void SetVortex_Tilting(su2double **PrimGrad_Flow, su2double* Vorticity, su2double LaminarViscosity) {}

  inline virtual su2double GetVortex_Tilting() {return 0.0; }

  inline virtual void SetDynamic_Derivative(unsigned short iVar, su2double der) {}

  inline virtual void SetDynamic_Derivative_n(unsigned short iVar, su2double der) {}

  inline virtual su2double GetDynamic_Derivative(unsigned short iVar) {return 0.0; }

  inline virtual su2double GetDynamic_Derivative_n(unsigned short iVar) {return 0.0; }

  inline virtual void SetDynamic_Derivative_Vel(unsigned short iVar, su2double der) {}

  inline virtual void SetDynamic_Derivative_Vel_n(unsigned short iVar, su2double der) {}

  inline virtual su2double GetDynamic_Derivative_Vel(unsigned short iVar) {return 0.0; }

  inline virtual su2double GetDynamic_Derivative_Vel_n(unsigned short iVar) {return 0.0; }

  inline virtual void SetDynamic_Derivative_Accel(unsigned short iVar, su2double der) {}

  inline virtual void SetDynamic_Derivative_Accel_n(unsigned short iVar, su2double der) {}

  inline virtual su2double GetDynamic_Derivative_Accel(unsigned short iVar) {return 0.0; }

  inline virtual su2double GetDynamic_Derivative_Accel_n(unsigned short iVar) {return 0.0; }

  inline virtual su2double GetSolution_Old_Vel(unsigned short iVar) {return 0.0; }

  inline virtual su2double GetSolution_Old_Accel(unsigned short iVar) {return 0.0; }

};
