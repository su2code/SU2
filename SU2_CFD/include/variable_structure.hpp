/*!
 * \file variable_structure.hpp
 * \brief Headers of the main subroutines for storing all the variables for 
 *        each kind of governing equation (direct, adjoint and linearized).
 *        The subroutines and functions are in the <i>variable_structure.cpp</i> file.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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

#include <cmath>
#include <iostream>

#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \class CVariable
 * \brief Main class for defining the variables.
 * \author F. Palacios.
 * \version 1.1.
 */
class CVariable {
protected:

	double *Solution,		/*!< \brief Solution of the problem. */
	*Solution_Old;			/*!< \brief Old solution of the problem R-K. */
	double **Gradient;		/*!< \brief Gradient of the solution of the problem. */ 
	double *Limiter;		/*!< \brief Limiter of the solution of the problem. */ 
	double *Solution_time_n,	/*!< \brief Solution of the problem at time n for dual-time stepping technique. */
	*Solution_time_n1;			/*!< \brief Solution of the problem at time n-1 for dual-time stepping technique. */
	double AuxVar;			/*!< \brief Auxiliar variable for gradient computation. */
	double *Grad_AuxVar;	/*!< \brief Gradient of the auxiliar variable. */
	double Delta_Time;	/*!< \brief Time step. */
	double Max_Lambda,	/*!< \brief Maximun eingenvalue. */
	Max_Lambda_Inv,		/*!< \brief Maximun inviscid eingenvalue. */
	Max_Lambda_Visc,	/*!< \brief Maximun viscous eingenvalue. */
	Lambda;				/*!< \brief Value of the eingenvalue. */
	double Sensor;	/*!< \brief Pressure sensor for high order central scheme. */
	double *Undivided_Laplacian;	/*!< \brief Undivided laplacian of the solution. */
	double *Residual,	/*!< \brief Residual of the problem. */
	*Res_Visc,			/*!< \brief Viscous Residual of the problem. */
	*Res_Sour,			/*!< \brief Source Residual of the problem. */
	**Res_Visc_RK,		/*!< \brief Viscous Residual for RK steps. */
	*Res_Conv,			/*!< \brief Convective Residual of the problem. */
	*Residual_Chemistry,	/*!< \brief Residual of the chemistry source terms. */
	*Residual_ElecForce,	/*!< \brief Residual of the electrostatic force source terms. */
	*Residual_MomentumExch, /*!< \brief Residual of the collisional momentum exchange source terms. */
	*Residual_EnergyExch,	/*! < \brief Residual of the collisional energy exchange source terms. */
	*Residual_Old,		/*!< \brief Auxiliar structure for residual smoothing. */
	*Residual_Sum;		/*!< \brief Auxiliar structure for residual smoothing. */
	double *TruncationError;	/*!< \brief Truncation error for multigrid cycle. */

	unsigned short nDim,	/*!< \brief Number of dimension of the problem. */
	nVar;					/*!< \brief Number of variables of the problem. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CVariable(void);

	/*!
	 * \overload 
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	CVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */
	virtual ~CVariable(void);

	/*!
	 * \brief Set the value of the solution.
	 * \param[in] val_solution - Solution of the problem.
	 */
	void SetSolution(double *val_solution);

	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
	 */
	void SetSolution(unsigned short val_var, double val_solution);

	/*!
	 * \brief Get the solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the solution for the index <i>val_var</i>.
	 */
	double GetSolution(unsigned short val_var);

	/*!
	 * \brief Set the value of the old solution.
	 * \param[in] val_solution_old - Pointer to the residual vector.
	 */
	void SetSolution_Old(double *val_solution_old);

	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution_old - Value of the old solution for the index <i>val_var</i>.
	 */	
	void SetSolution_Old(unsigned short val_var, double val_solution_old);

	/*!
	 * \brief Set old variables to the value of the current variables.
	 */
	void Set_OldSolution(void);

	/*!
	 * \brief Set variables to the value of the old variables.
	 */
	void Set_Solution(void);	

	/*!
	 * \brief Set the variable solution at time n.
	 */	
	void Set_Solution_time_n(void);

	/*!
	 * \brief Set the variable solution at time n-1.
	 */	
	void Set_Solution_time_n1(void);

	/*!
	 * \brief Set to zero velocity components of the solution.
	 */
	void SetVelSolutionZero(void);

	/*!
	 * \brief Set to zero velocity components of the solution.
	 */
	void SetVelSolutionOldZero(void);

	/*!
	 * \brief Set to zero the solution.
	 */	
	void SetSolutionZero(void);

	/*!
	 * \brief Add a value to the solution.
	 * \param[in] val_var - Number of the variable.
	 * \param[in] val_solution - Value that we want to add to the solution.
	 */
	void AddSolution(unsigned short val_var, double val_solution);

	/*!
	 * \brief Get the solution of the problem.
	 * \return Pointer to the solution vector.
	 */
	double *GetSolution(void);

	/*!
	 * \brief Get the old solution of the problem (Runge-Kutta method)
	 * \return Pointer to the old solution vector.
	 */
	double *GetSolution_Old(void);

	/*!
	 * \brief Get the solution at time n.
	 * \return Pointer to the solution (at time n) vector.
	 */	
	double *GetSolution_time_n(void);

	/*!
	 * \brief Get the solution at time n-1.
	 * \return Pointer to the solution (at time n-1) vector.
	 */	
	double *GetSolution_time_n1(void);

	/*!
	 * \brief Set the value of the residual.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual for the index <i>val_var</i>.
	 */
	void SetResidual(unsigned short val_var, double val_residual);

	/*!
	 * \overload
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SetResidual(double *val_residual);

	/*!
	 * \brief Set the value of the convective residual.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SetRes_Conv(double *val_residual);

	/*!
	 * \brief Set the value of the viscous residual.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SetRes_Visc(double *val_residual);
	
	/*!
	 * \brief Set the value of the source residual.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SetRes_Sour(double *val_residual);

	/*!
	 * \brief Set the value of the viscous residual in the Runge-Kutta iteration.
	 * \param[in] val_residual - Pointer to the residual vector.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void SetRes_Visc_RK(double *val_residual, unsigned short iRKStep);

	/*!
	 * \brief Set the value of the old residual.
	 * \param[in] val_residual_old - Pointer to the residual vector.
	 */
	void SetResidual_Old(double *val_residual_old);

	/*!
	 * \brief Add a value to the residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void AddResidual(double *val_residual);

	/*!
	 * \brief Add a value to the residual vector.
	 * \param[in] val_residual - residual value.
	 */
	void AddResidual(double val_residual);

	/*!
	 * \brief Add a value to the convective residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void AddRes_Conv(double *val_residual);

	/*!
	 * \brief Add a value to the viscous residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */	
	void AddRes_Visc(double *val_residual);
	
	/*!
	 * \brief Add a value to the source residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */	
	void AddRes_Sour(double *val_residual);

	/*!
	 * \brief Add a value to the summed residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void AddResidual_Sum(double *val_residual);

	/*!
	 * \brief Set summed residual vector to zero value.
	 */
	void SetResidualSumZero(void);

	/*!
	 * \brief Subtract a value to the residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SubtractResidual(double *val_residual);

	/*!
	 * \brief Subtract a value to the viscous vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SubtractRes_Visc(double *val_residual);
	
	/*!
	 * \brief Subtract a value to the source vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void SubtractRes_Sour(double *val_residual);

	/*!
	 * \brief Subtract a value to the convective residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */	
	void SubtractRes_Conv(double *val_residual);

	/*!
	 * \brief Set the residual to zero value.
	 */
	void SetResidualZero(void);

	/*!
	 * \brief Set the residual of the velocity to zero value.
	 */
	void SetVelResidualZero(void);

	/*!
	 * \brief Set the viscous velocity residual to zero value.
	 */		
	void SetVel_ResVisc_Zero(void);
	
	/*!
	 * \brief Set the viscous velocity residual to zero value.
	 */		
	void SetVel_ResSour_Zero(void);

	/*!
	 * \brief Set the viscous velocity residual to zero value.
	 */
	void SetVel_ResVisc_Zero(unsigned short iSpecies);
	
	/*!
	 * \brief Set the viscous velocity residual to zero value.
	 */
	void SetVel_ResSour_Zero(unsigned short iSpecies);

	/*!
	 * \brief Set the convective velocity residual to zero value.
	 */		
	void SetVel_ResConv_Zero(void);

	/*!
	 * \brief Set the convective velocity residual to zero value.
	 */
	void SetVel_ResConv_Zero(unsigned short iSpecies);

	/*!
	 * \brief Set the viscous residual to zero value.
	 */	
	void Set_ResVisc_Zero(void);

	/*!
	 * \brief Set the convective residual to zero value.
	 */	
	void Set_ResConv_Zero(void);
	
	/*!
	 * \brief Set the convective residual to zero value.
	 */	
	void Set_ResSour_Zero(void);

	/*!
	 * \brief Get the residual value.
	 * \return Pointer to the residual vector.
	 */
	double *GetResidual(void);

	/*!
	 * \brief Get the convective residual value.
	 * \return Pointer to the convective residual vector.
	 */	
	double *GetResConv(void);
	
	/*!
	 * \brief Get the viscous residual value.
	 * \return Pointer to the viscous residual vector.
	 */	
	double *GetResVisc(void);
	
	/*!
	 * \brief Get the source residual value.
	 * \return Pointer to the viscous residual vector.
	 */	
	double *GetResSour(void);

	/*!
	 * \brief Get the viscous R-K residual value.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \return Pointer to the viscous R-K  residual vector.
	 */
	double GetRes_Visc_RK(unsigned short val_var, unsigned short iRKStep);

	/*!
	 * \brief Get the residual value.
	 * \param[out] val_residual - Pointer to the residual vector.
	 */	
	void GetResidual(double *val_residual);

	/*!
	 * \brief Get the residual value.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the residual vector for the variable <i>val_var</i>.
	 */	
	double GetResidual(unsigned short val_var);

	/*!
	 * \brief Get the value of the summed residual.
	 * \return Pointer to the summed residual.
	 */	
	double *GetResidual_Sum(void);

	/*!
	 * \brief Get the value of the old residual.
	 * \return Pointer to the old residual.
	 */	
	double *GetResidual_Old(void);

	/*!
	 * \brief Get the value of the summed residual.
	 * \param[out] val_residual - Pointer to the summed residual.
	 */	
	void GetResidual_Sum(double *val_residual);

	/*!
	 * \brief Get the projected residual.
	 * \param[in] val_vector - Direction of the projection.
	 * \return Value of the projected residual.
	 */
	double GetProjRes(double *val_vector);

	/*!
	 * \brief Set auxiliar variables, we are looking for the gradient of that variable.
	 * \param[in] val_auxvar - Value of the auxiliar variable.
	 */
	void SetAuxVar(double val_auxvar);

	/*!
	 * \brief Get the value of the auxiliary variable.
	 * \return Value of the auxiliary variable.
	 */
	double GetAuxVar(void);

	/*!
	 * \brief Set the auxiliary variable gradient to zero value.
	 */
	void SetAuxVarGradientZero(void);

	/*!
	 * \brief Set the value of the auxiliary variable gradient.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_gradient - Value of the gradient for the index <i>val_dim</i>.
	 */
	void SetAuxVarGradient(unsigned short val_dim, double val_gradient);

	/*!
	 * \brief Add a value to the auxiliary variable gradient.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient to be added for the index <i>val_dim</i>.
	 */		
	void AddAuxVarGradient(unsigned short val_dim, double val_value);

	/*!
	 * \brief Subtract a value to the auxiliary variable gradient.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient to be subtracted for the index <i>val_dim</i>.
	 */		
	void SubtractAuxVarGradient(unsigned short val_dim, double val_value);

	/*!
	 * \brief Get the gradient of the auxiliary variable.
	 * \return Value of the gradient of the auxiliary variable.
	 */		
	double *GetAuxVarGradient(void);

	/*!
	 * \brief Get the gradient of the auxiliary variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the gradient of the auxiliary variable for the dimension <i>val_dim</i>.
	 */		
	double GetAuxVarGradient(unsigned short val_dim);	

	/*!
	 * \brief Add a value to the truncation error.
	 * \param[in] val_truncation_error - Value that we want to add to the truncation error.
	 */		
	void AddTruncationError(double *val_truncation_error);

	/*!
	 * \brief Subtract a value to the truncation error.
	 * \param[in] val_truncation_error - Value that we want to subtract to the truncation error.
	 */		
	void SubtractTruncationError(double *val_truncation_error);

	/*!
	 * \brief Set the truncation error to zero.
	 */		
	void SetTruncationErrorZero(void);

	/*!
	 * \brief Set the velocity of the truncation error to zero.
	 */		
	void SetVelTruncationErrorZero(void);

	/*!
	 * \brief Set the velocity of the truncation error to zero.
	 */
	void SetVelTruncationErrorZero(unsigned short iSpecies);

	/*!
	 * \brief Get the truncation error.
	 * \return Pointer to the truncation error.
	 */	
	double *GetTruncationError(void);

	/*!
	 * \brief Get the truncation error.
	 * \param[out] val_trunc_error - Pointer to the truncation error.
	 */	
	void GetTruncationError(double *val_trunc_error);

	/*!
	 * \brief Set the gradient of the solution.
	 * \param[in] val_gradient - Gradient of the solution.
	 */
	void SetGradient(double **val_gradient);

	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Set to zero the gradient of the solution.
	 */
	void SetGradientZero(void);

	/*!
	 * \brief Add <i>val_value</i> to the solution gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the solution gradient.
	 */
	void AddGradient(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Subtract <i>val_value</i> to the solution gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the solution gradient.
	 */
	void SubtractGradient(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Get the value of the solution gradient.
	 * \return Value of the gradient solution.
	 */
	double **GetGradient(void);

	/*!
	 * \brief Get the value of the solution gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the solution gradient.
	 */
	double GetGradient(unsigned short val_var, unsigned short val_dim);

	/*!
	 * \brief Set the value of the limiter.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_limiter - Value of the limiter for the index <i>val_var</i>.
	 */
	void SetLimiter(unsigned short val_var, double val_limiter);

	/*!
	 * \brief Get the value of the slope limiter.
	 * \return Pointer to the limiters vector.
	 */
	double *GetLimiter(void);

	/*!
	 * \brief Get the value of the slope limiter.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the limiter vector for the variable <i>val_var</i>.
	 */
	double GetLimiter(unsigned short val_var);

	/*!
	 * \brief Set the value of the time step.
	 * \param[in] val_delta_time - Value of the time step.
	 */
	void SetDelta_Time(double val_delta_time);

	/*!
	 * \brief Set the value of the time step.
	 * \param[in] val_delta_time - Value of the time step.
	 * \param[in] iSpecies - Index of the Species .
	 */
	virtual void SetDelta_Time(double val_delta_time, unsigned short iSpecies);

	/*!
	 * \brief Get the value of the time step.
	 * \return Value of the time step.
	 */
	double GetDelta_Time(void);

	/*!
	 * \brief Get the value of the time step.
	 * \param[in] iSpecies - Index of the Species
	 * \return Value of the time step.
	 */
	virtual double GetDelta_Time(unsigned short iSpecies);

	/*!
	 * \brief Set the value of the maximum eigenvalue.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue.
	 */
	void SetMax_Lambda(double val_max_lambda);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	void SetMax_Lambda_Inv(double val_max_lambda);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_Species - Value of the species index to set the maximum eigenvalue.
	 */
	virtual void SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_Species);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	void SetMax_Lambda_Visc(double val_max_lambda);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_Species - Index of the species to set the maximum eigenvalue of the viscous terms.
	 */
	virtual void SetMax_Lambda_Visc(double val_max_lambda, unsigned short val_Species);

	/*!
	 * \brief Add a value to the maximum eigenvalue.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue.
	 */
	void AddMax_Lambda(double val_max_lambda);

	/*!
	 * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	void AddMax_Lambda_Inv(double val_max_lambda);

	/*!
	 * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	void AddMax_Lambda_Visc(double val_max_lambda);

	/*!
	 * \brief Get the value of the maximum eigenvalue.
	 * \return the value of the maximum eigenvalue.
	 */
	double GetMax_Lambda(void);

	/*!
	 * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */	
	double GetMax_Lambda_Inv(void);

	/*!
	 * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	double GetMax_Lambda_Visc(void);

	/*!
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda - Value of the spectral radius.
	 */
	void SetLambda(double val_lambda);

	/*!
	 * \brief Add the value of the spectral radius.
	 * \param[in] val_lambda - Value of the spectral radius.
	 */
	void AddLambda(double val_lambda);

	/*!
	 * \brief Get the value of the spectral radius.
	 * \return Value of the spectral radius.
	 */
	double GetLambda(void);

	/*!
	 * \brief Set pressure sensor.
	 * \param[in] val_sensor - Value of the pressure sensor.
	 */
	void SetSensor(double val_sensor);

	/*!
	 * \brief Get the pressure sensor.
	 * \return Value of the pressure sensor.
	 */	
	double GetSensor(void);

	/*!
	 * \brief Set the value of the undivided laplacian of the solution.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_undivided_laplacian - Value of the undivided solution for the index <i>val_var</i>.
	 */
	void SetUndivided_Laplacian(unsigned short val_var, double val_undivided_laplacian);

	/*!
	 * \brief Add the value of the undivided laplacian of the solution.
	 * \param[in] val_und_lapl - Value of the undivided solution.
	 */	
	void AddUnd_Lapl(double *val_und_lapl);

	/*!
	 * \brief Subtract the value of the undivided laplacian of the solution.
	 * \param[in] val_und_lapl - Value of the undivided solution.
	 */		
	void SubtractUnd_Lapl(double *val_und_lapl);
	
	/*!
	 * \brief Subtract the value of the undivided laplacian of the solution.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \param[in] val_und_lapl - Value of the undivided solution.
	 */		
	void SubtractUnd_Lapl(unsigned short val_var, double val_und_lapl);

	/*!
	 * \brief Set the undivided laplacian of the solution to zero.
	 */			
	void SetUnd_LaplZero(void);

	/*!
	 * \brief Set a value to the undivided laplacian.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \param[in] val_und_lapl - Value of the undivided laplacian.
	 */	
	void SetUnd_Lapl(unsigned short val_var, double val_und_lapl);

	/*!
	 * \brief Get the undivided laplacian of the solution.
	 * \return Pointer to the undivided laplacian vector.
	 */
	double *GetUnd_Lapl(void);

	/*!
	 * \brief Get the undivided laplacian of the solution.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \return Value of the undivided laplacian vector.
	 */
	double GetUnd_Lapl(unsigned short val_var);

	/*!
	 * \brief A virtual member.
	 * \return Value of the flow density.
	 */
	virtual double GetDensity(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the flow density.
	 */
	virtual double GetDensity(unsigned short val_iSpecies);

	/*!
	 * \brief A virtual member.
	 * \return Value of the flow energy.
	 */	
	virtual double GetEnergy(void);

	/*!
	 * \brief A virtual member.
	 * \return Pointer to the force projection vector.
	 */
	virtual double *GetForceProj_Vector(void);

	/*!
	 * \brief A virtual member.
	 * \return Pointer to the internal boundary vector.
	 */
	virtual double *GetIntBoundary_Jump(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the eddy viscosity.
	 */		
	virtual double GetEddyViscosity(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the flow enthalpy.
	 */		
	virtual double GetEnthalpy(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the flow pressure.
	 */		
	virtual double GetPressure(void);
	
	/*!
	 * \brief Get the species energy of formation.
	 * \return Value of the energy of formation.
	 * \param[in] val_Species - Index of the species.
	 */
	virtual double GetEnthalpyFormation(unsigned short val_Species);

	/*!
	 * \brief A virtual member.
	 * \return Value of the linearized pressure.
	 */		
	virtual double GetDeltaPressure(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_vector - Direction of projection.
	 * \return Value of the projected velocity.
	 */		
	virtual double GetProjVel(double *val_vector);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_vector - Direction of projection.
	 * \param[in] val_Species - Index of the desired species.
	 * \param[in] val_nDiatomics - Number of diatomic species.
	 * \return Value of the projected velocity.
	 */		
	virtual double GetProjVel(double *val_vector, unsigned short val_Species, unsigned short val_nDiatomics);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_vector - Direction of projection.
	 * \return Value of the projected velocity for the incompressible flow.
	 */		
	virtual double GetProjVelInc(double *val_vector);

	/*! 
	 * \ Overloaded for multi species equations
	 * \param[in] val_vector - Direction of projection.
	 * \param[in] val_Fluid - Index of the fluid to be projected.
	 * \return Value of the projected velocity.
	 */
	virtual double GetProjVel(double *val_vector, unsigned short val_Fluid);	

	/*!
	 * \brief A virtual member.
	 * \return Value of the sound speed.
	 */		
	virtual double GetSoundSpeed(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the sound speed for the incompressible flow.
	 */		
	virtual double GetSoundSpeedInc(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density for the incompressible flow.
	 */		
	virtual double GetDensityInc(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the beta for the incompressible flow.
	 */		
	virtual double GetBetaInc2(void);

	/*! Overloaded for plasma equations
	 * \brief A virtual member.
	 * \return Value of the sound speed of Fluid val_Fluid
	 */
	virtual double GetSoundSpeed(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return Value of the temperature.
	 */		
	virtual double GetTemperature(void);

	/*! Overloaded for plasma equations
	 * \brief A virtual member.
	 * \return Value of the temperature.
	 */
	virtual double GetTemperature(unsigned short val_iSpecies);
	
	/*!
	 * \brief Overloaded for plasma equations
	 * \param[in] iSpecies - Index of species for desired TR temperature.
	 */
	virtual double GetTemperature_TR(unsigned short iSpecies);

	/*
	 * \brief A virtual member.
	 * \return Value of the thermal coefficient.
	 */
	virtual double GetThermalCoeff(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the Theta variable of the adjoint problem.
	 */		
	virtual double GetTheta(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the velocity for the dimension <i>val_dim</i>.
	 */		
	virtual double GetVelocity(unsigned short val_dim);

	/*!
	 * \brief A virtual member.
	 * \return Norm 2 of the velocity vector.
	 */		
	virtual double GetVelocity2(void);

	/*!
	 * \brief A virtual member.
	 * \return Pressure of Fluid val_Fluid
	 */	
	virtual double GetPressure(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return Norm 2 of the velocity vector of Fluid val_Fluid.
	 */	
	virtual double GetVelocity2(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return val_dim component of velocity vector of Fluid val_Fluid.
	 */	
	virtual double GetVelocity(unsigned short val_dim,unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return The laminar viscosity of the flow.
	 */		
	virtual double GetLaminarViscosity(void);

	/*!
	 * \brief A virtual member.
	 * \return The laminar viscosity of the incompressible flow.
	 */		
	virtual double GetLaminarViscosityInc(void);

	/*!
	 * \brief A virtual member.
	 * \return The laminar viscosity of the flow.
	 */
	virtual double GetLaminarViscosity(unsigned short iSpecies);

	/*!
	 * \brief A virtual member.
	 * \return The Eddy viscosity of the flow.
	 */
	virtual double GetEddyViscosity(unsigned short iSpecies);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */		
	virtual double GetVorticity(unsigned short val_dim);

	/*!
	 * \brief A virtual member.
	 * \return Value of the rate of strain magnitude.
	 */
	virtual double GetStrainMag(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
	 */		
	virtual void SetForceProj_Vector(double *val_ForceProj_Vector);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump.
	 */		
	virtual void SetIntBoundary_Jump(double *val_IntBoundary_Jump);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_Kind_Turb_Model - Kind of turbulence model.
	 * \param[in] Turb_Solution - Solution of the turbulence model.
	 */		
	virtual void SetEddyViscosity(unsigned short val_Kind_Turb_Model, CVariable *Turb_Solution);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 */		
	virtual void SetEddyViscosity(double val_eddy_viscosity);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetEnthalpy(void);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetDensityInc(double val_densityinc);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetBetaInc2(double val_betainc2);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */		
	virtual void SetPhi_Old(double *val_phi);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Ratio of Specific heats
	 * \param[in] Coord - Coordinates of the point where the pressure is computed
	 */
	virtual void SetPressure(double Gamma, double *Coord);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Ratio of Specific heats
	 */
	virtual void SetPressure(double Gamma);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetPressure(double GammaMonatomic, double GammaDiatomic, double *Coord);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetPressure(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] Gamma - Ratio of Specific heats
	 */		
	virtual void SetDeltaPressure(double *val_velocity, double Gamma);

	/*!
	 * \brief A virtual member.
	 * \param[in] GammaMonatomic - Ratio of Specific heats for Monatomic.
	 * \param[in] GammaDiatomic - Ratio of Specific heats for Diatomic.
	 */		
	virtual void SetSoundSpeed(double GammaMonatomic, double GammaDiatomic);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetSoundSpeed(double Gamma);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetSoundSpeed(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gas_Constant - Value of the Gas Constant
	 */		
	virtual void SetTemperature(double Gas_Constant);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gas_Constant - Value of the Gas Constant
	 */
	virtual void SetTemperature(double *Gas_Constant);

	/*!
	 * \brief A virtual member.
	 * \param[in] Molar_Mass - Value of the molar mass (kg/kmol) for each species.
	 * \param[in] GammaMonatomic - Value of the ratio of specific heats for monatomic species.
	 * \param[in] GammaDiatomic - Value of the ratio of specific heats for diatomic species.
	 */
	virtual void SetTemperature_TR(double *Molar_Mass, double GammaMonatomic, double GammaDiatomic);
	
	/*!
	 * \brief A virtual member.
	 * \param[in] Temperature_Wall - Value of the Temperature at the wall
	 */
	virtual void SetWallTemperature(double Temperature_Wall);

	/*!
	 * \brief A virtual member.
	 * \param[in] Temperature_Wall - Value of the Temperature at the wall
	 */
	virtual void SetWallTemperature(double* Temperature_Wall);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Value of the Gamma
	 * \param[in] Gas_Constant - Value of the Gas Constant
	 */		
	virtual void SetThermalCoeff(double Gamma, double Gas_Constant);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Value of the Gamma
	 * \param[in] Gas_Constant - Value of the Gas Constant of the various species
	 */
	virtual void SetThermalCoeff(double Gamma, double *Gas_Constant);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_enthalpy - Value of the enthalpy.
	 */		
	virtual void SetTheta(double val_density, double *val_velocity, double val_enthalpy);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */		
	virtual void SetVelocity(double *val_velocity);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVelocityInc2(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetBetaInc2(CConfig *config);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetPressureInc(void);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetPressureInc(double val_pressure);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetSoundSpeedInc(double *val_normal);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetVelocity2(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */		
	virtual void SetVelocity_Old(double *val_velocity);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] iSpecies - Index of the species to set the velocity.
	 */
	virtual void SetVelocity_Old(double *val_velocity, unsigned short iSpecies);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	virtual void SetLaminarViscosity(CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 */		
	virtual void SetLaminarViscosity(double val_laminar_viscosity);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_laminar_viscosity_inc - Value of the laminar viscosity (incompressible flows).
	 */		
	virtual void SetLaminarViscosityInc(double val_laminar_viscosity_inc);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetVorticity(void);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetStrainMag(void);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVelSolutionOldDVector(void);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVelSolutionDVector(void);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetGradient_PrimitiveZero(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the primitive variables.
	 */
	virtual void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
	 */
	virtual void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief A virtual member.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double **GetGradient_Primitive(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of source term contribution.
	 */		
	virtual double GetSource(unsigned short val_var);

	/*!
	 * \brief A virtual member.
	 * \return Sets the value source terms contribution.
	 */	
	virtual void SetSource(double * res);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_ivar - Index for i variable number.
	 * \param[in] val_jvar - Index for j variable number.
	 * \return Value of source term contribution.
	 */		
	virtual double GetSource_Jacobian(unsigned short val_ivar,unsigned short val_jvar);

	/*!
	 * \brief A virtual member.
	 * \return Sets the value source terms contribution.
	 */	
	virtual void SetSourceJacobian(double ** src_jacobian);

	/*!
	 * \brief Set the blending function for the blending of k-w and k-eps.
	 * \param[in] val_viscosity - Value of the vicosity.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_dist - Value of the distance to the wall.
	 */
	virtual void SetBlendingFunc(double val_viscosity, double val_dist, double val_density);

	/*!
	 * \brief Get the first blending function of the SST model.
	 */
	virtual double GetF1blending(void);

	/*!
	 * \brief Get the second blending function of the SST model.
	 */
	virtual double GetF2blending(void);

	/*!
	 * \brief Get the value of the cross diffusion of tke and omega.
	 */
	virtual double GetCrossDiff(void){ return 0.0; };

	/*!
	 * \brief Get the value of the eddy viscosity.
	 * \return the value of the eddy viscosity.
	 */
	virtual double GetmuT(void);

	/*!
	 * \brief Set the value of the eddy viscosity.
	 * \param[in] val_muT
	 */
	virtual void SetmuT(double val_muT);

	/*!
	 * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] iFluids - ID of the fluid in multispecies solver
	 * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	virtual double GetMax_Lambda_Inv(unsigned short iFluids);

	/*!
	 * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] iFluids - ID of the fluid in multispecies solver
	 * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	virtual double GetMax_Lambda_Visc(unsigned short iFluids);

	/*!
	 * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] iSpecies - Value of iSpecies to which the eigenvalue belongs
	 */
	virtual void AddMax_Lambda_Inv(double val_max_lambda, unsigned short iSpecies);

	/*!
	 * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] iSpecies - Value of iSpecies to which the eigenvalue belongs
	 */
	virtual void AddMax_Lambda_Visc(double val_max_lambda, unsigned short iSpecies);
	
	/*!
	 * \brief A virtual member.
	 * \param[in] val_difflevelset - Value of the diff level set (value-target).
	 */	
	virtual void SetDiffLevelSet(double val_difflevelset);
	
	/*!
	 * \brief A virtual member.
	 */		
	virtual double GetDiffLevelSet(void);
};

/*! 
 * \class CPotentialVariable
 * \brief Main class for defining the variables of the potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios.
 * \version 1.1.
 */
class CPotentialVariable : public CVariable {
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CPotentialVariable(void);

	/*!
	 * \overload
	 * \param[in] val_potential - Value of the potential solution (initialization value).		 
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CPotentialVariable(double val_potential, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CPotentialVariable(void);
};

/*! 
 * \class CWaveVariable
 * \brief Main class for defining the variables of the wave equation solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios.
 * \version 1.1.
 */
class CWaveVariable : public CVariable {
public:
	
	/*!
	 * \brief Constructor of the class. 
	 */
	CWaveVariable(void);
	
	/*!
	 * \overload
	 * \param[in] val_potential - Value of the potential solution (initialization value).		 
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CWaveVariable(double val_wave, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
	
	/*!
	 * \brief Destructor of the class. 
	 */	
	~CWaveVariable(void);
};

/*! 
 * \class CEulerVariable
 * \brief Main class for defining the variables of the Euler's solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CEulerVariable : public CVariable {
protected:
	double Velocity2;	/*!< \brief Square of the velocity vector. */
	double BetaInc2;	/*!< \brief Square of the artificial compresibility factor. */
	double SoundSpeedInc;	/*!< \brief Speed of Sound in the fluid (incompressible flow). */
	double Pressure;	/*!< \brief Pressure of the fluid. */
	double SoundSpeed;	/*!< \brief Speed of Sound in the fluid. */
	double Enthalpy;	/*!< \brief Enthalpy of the fluid. */
	double Temperature; /*!< \brief Temperature of the fluid. */
	double DensityInc; /*!< \brief Density of the fluid (used in free surface incompressible flows). */
	double Kappa;		/*!< \brief Thermal conductivity of the fluid. */
	double **Gradient_Primitive;	/*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */ 
	bool incompressible;
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CEulerVariable(void);

	/*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_energy - Value of the flow energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.		 
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_ndim, 
			unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */		
	virtual ~CEulerVariable(void);

	/*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	void SetGradient_PrimitiveZero(void);

	/*!
	 * \brief Add <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the primitive variables.
	 */
	void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
	 */
	void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim);

	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double **GetGradient_Primitive(void);

	/*!
	 * \brief Set the value of the velocity*velocity for the incompressible solver.
	 */
	void SetVelocityInc2(void);

	/*!
	 * \brief Set the value of the beta*beta for the incompressible solver.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetBetaInc2(CConfig *config);

	/*!
	 * \brief Set the value of the pressure for the incompressible solver.
	 */
	void SetPressureInc(void);

	/*!
	 * \brief Set the value of the pressure for the incompressible solver.
	 */
	void SetPressureInc(double val_pressure);

	/*!
	 * \brief Set the value of the sound speed for the incompressible solver.
	 * \param[in] val_normal - Value of the normal vector with modulus.
	 */
	void SetSoundSpeedInc(double *val_normal);

	/*!
	 * \brief Set the value of the velocity*velocity.
	 */
	void SetVelocity2(void);

	/*!
	 * \brief Set the value of the pressure.
	 */
	void SetPressure(double Gamma, double *Coord);

	/*!
	 * \brief Set the value of the speed of the sound.
	 * \param[in] Gamma - Value of Gamma.
	 */
	void SetSoundSpeed(double Gamma);

	/*!
	 * \brief Set the value of the enthalpy.
	 */
	void SetEnthalpy(void);

	/*!
	 * \brief Set the value of the density for the incompressible flows.
	 */
	void SetDensityInc(double val_densityinc);

	/*!
	 * \brief Set the value of the beta coeffient for incompressible flows.
	 */
	void SetBetaInc2(double val_betainc2);

	/*!
	 * \brief Set the value of the temperature.
	 * \param[in] Gas_Constant - Value of Gas Constant
	 */
	void SetTemperature(double Gas_Constant);

	/*!
	 * \brief Get the norm 2 of the velocity.
	 * \return Norm 2 of the velocity vector.
	 */
	double GetVelocity2(void);

	/*!
	 * \brief Get the flow pressure.
	 * \return Value of the flow pressure.
	 */
	double GetPressure(void);

	/*!
	 * \brief Get the speed of the sound.
	 * \return Value of speed of the sound.
	 */
	double GetSoundSpeed(void);

	/*!
	 * \brief Get the speed of the sound for the incompressible flow
	 * \return Value of speed of the sound.
	 */
	double GetSoundSpeedInc(void);

	/*!
	 * \brief Get the value of density for the incompressible flow
	 * \return Value of beta squared.
	 */
	double GetDensityInc(void);

	/*!
	 * \brief Get the value of beta squared for the incompressible flow
	 * \return Value of beta squared.
	 */
	double GetBetaInc2(void);

	/*!
	 * \brief Get the enthalpy of the flow.
	 * \return Value of the enthalpy of the flow.
	 */
	double GetEnthalpy(void);

	/*!
	 * \brief Get the density of the flow.
	 * \return Value of the density of the flow.
	 */
	double GetDensity(void);

	/*!
	 * \brief Get the energy of the flow.
	 * \return Value of the energy of the flow.
	 */
	double GetEnergy(void);

	/*!
	 * \brief Get the temperature of the flow.
	 * \return Value of the temperature of the flow.
	 */
	double GetTemperature(void);

	/*!
	 * \brief Get the thermal coefficient of the flow.
	 * \return Value of the thermal coefficient of the flow.
	 */	
	double GetThermalCoeff(void);

	/*!
	 * \brief Get the velocity of the flow.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the velocity for the dimension <i>val_dim</i>.
	 */
	double GetVelocity(unsigned short val_dim);

	/*!
	 * \brief Get the projected velocity in a unitary vector direction (compressible solver).
	 * \param[in] val_vector - Direction of projection.
	 * \return Value of the projected velocity.
	 */
	double GetProjVel(double *val_vector);

	/*!
	 * \brief Get the projected velocity in a unitary vector direction (incompressible solver).
	 * \param[in] val_vector - Direction of projection.
	 * \return Value of the projected velocity.
	 */
	double GetProjVelInc(double *val_vector);

	/*!
	 * \brief Set the velocity vector from the solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */	
	void SetVelocity(double *val_velocity);

	/*!
	 * \brief Set the velocity vector from the old solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */		
	void SetVelocity_Old(double *val_velocity);
};

/*! 
 * \class CNSVariable
 * \brief Main class for defining the variables of the Navier-Stokes' solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CNSVariable : public CEulerVariable {
private:
	double Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
	double Prandtl_Turb;      /*!< \brief Turbulent Prandtl number. */
	double Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
	double Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
	double Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
	double LaminarViscosity;	/*!< \brief Viscosity of the fluid. */
	double LaminarViscosityInc;	/*!< \brief Viscosity of the fluid (incompressible flows). */
	double EddyViscosity;		/*!< \brief Eddy viscosity of the fluid. */
	double Vorticity[3];		/*!< \brief Vorticity of the fluid. */
	double StrainMag;           /*!< \brief Magnitude of rate of strain tensor. */
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CNSVariable(void);

	/*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_energy - Value of the flow energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.		 
	 * \param[in] config - Definition of the particular problem.
	 */
	CNSVariable(double val_density, double *val_velocity, 
			double val_energy, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	
	 */
	CNSVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CNSVariable(void);

	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] config - Definition of the particular problem.	
	 */
	void SetLaminarViscosity(CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity);

	/*!
	 * \overload
	 * \param[in] val_laminar_viscosity_inc - Value of the laminar viscosity (incompressible flows).
	 */
	void SetLaminarViscosityInc(double val_laminar_viscosity_inc);

	/*!
	 * \brief Set the vorticity value.
	 */
	void SetVorticity(void);

	/*!
	 * \brief Set the rate of strain magnitude.
	 */
	void SetStrainMag(void);

	/*!
	 * \brief Set the thermal coefficient.
	 * \param[in] Gamma - Ratio of specific heats
	 * \param[in] Gas_Constant - The Gas Constant
	 */
	void SetThermalCoeff(double Gamma, double Gas_Constant);

	/*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_Kind_Turb_Model - Kind of turbulence model.
	 * \param[in] Turb_Solution - Solution of the turbulence model.
	 */
	void SetEddyViscosity(unsigned short val_Kind_Turb_Model, CVariable *Turb_Solution);

	/*!
	 * \overload
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 */
	void SetEddyViscosity(double val_eddy_viscosity);

	/*!
	 * \brief Get the laminar viscosity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetLaminarViscosity(void);

	/*!
	 * \brief Get the laminar viscosity of the incompressible flow.
	 * \return Value of the laminar viscosity of the incompressible flow.
	 */
	double GetLaminarViscosityInc(void);

	/*!
	 * \brief Get the eddy viscosity of the flow.
	 * \return The eddy viscosity of the flow.
	 */
	double GetEddyViscosity(void);

	/*!
	 * \brief Set the temperature at the wall
	 */
	void SetWallTemperature(double temperature_wall);

	/*!
	 * \brief Get the value of the vorticity.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */	
	double GetVorticity(unsigned short val_dim);

	/*!
	 * \brief Get the value of the magnitude of rate of strain.
	 * \return Value of the rate of strain magnitude.
	 */
	double GetStrainMag(void);
};

/*! 
 * \class CTurbVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.1.
 */
class CTurbVariable : public CVariable {
protected:
	double muT;                /*!< \brief Eddy viscosity. */
public:
	/*!
	 * \brief Constructor of the class. 
	 */	
	CTurbVariable(void);

	/*!
	 * \overload
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CTurbVariable(void);

	/*!
	 * \brief Get the value of the eddy viscosity.
	 * \return the value of the eddy viscosity.
	 */
	double GetmuT();

	/*!
	 * \brief Set the value of the eddy viscosity.
	 * \param[in] val_muT - Value of the eddy viscosity.
	 */
	void SetmuT(double val_muT);
};

/*!
 * \class CTurbSAVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.1.
 */

class CTurbSAVariable : public CTurbVariable {
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSAVariable(void);

	/*!
	 * \overload
	 * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
	 * \param[in] val_muT  - The eddy viscosity
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */
	~CTurbSAVariable(void);
};

/*! 
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.1.
 */

class CTurbSSTVariable : public CTurbVariable {
protected:
	double F1;		/*!< \brief Menter blending function for blending of k-w and k-eps. */
	double F2;		/*!< \brief Menter blending function for stress limiter. */
	double CDkw;    /*!< \brief Cross-diffusion. */
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSSTVariable(void);

	/*!
	 * \overload
	 * \param[in] val_rho_kine - Turbulent variable value (initialization value).
	 * \param[in] val_rho_omega - Turbulent variable value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSSTVariable(double val_rho_kine, double val_rho_omega, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSSTVariable(void);

	/*!
	 * \brief Set the blending function for the blending of k-w and k-eps.
	 * \param[in] val_viscosity - Value of the vicosity.
	 * \param[in] val_dist - Value of the distance to the wall.
	 * \param[in] val_density - Value of the density.
	 */
	void SetBlendingFunc(double val_viscosity, double val_dist, double val_density);

	/*!
	 * \brief Get the first blending function.
	 */
	double GetF1blending(void);

	/*!
	 * \brief Get the second blending function.
	 */
	double GetF2blending(void);

	/*!
	 * \brief Get the value of the cross diffusion of tke and omega.
	 */
	double GetCrossDiff(void);
};

/*!
 * \class CAdjPotentialVariable
 * \brief Main class for defining the variables of the adjoint potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios.
 * \version 1.1.
 */
class CAdjPotentialVariable : public CVariable {
private:
	double Psi;			/*!< \brief Value of the adjoint variable. */
	double *ForceProj_Vector;	/*!< \brief Vector d. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CAdjPotentialVariable(void);

	/*!
	 * \overload
	 * \param[in] val_psi - Potential adjoint variable value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAdjPotentialVariable(double val_psi, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */

	~CAdjPotentialVariable(void);
};

/*! 
 * \class CAdjEulerVariable
 * \brief Main class for defining the variables of the adjoint Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CAdjEulerVariable : public CVariable {
protected:
	double *Psi;		/*!< \brief Vector of the adjoint variables. */
	double *ForceProj_Vector;	/*!< \brief Vector d. */
	double *IntBoundary_Jump;	/*!< \brief Interior boundary jump vector. */
	double Theta;		/*!< \brief Theta variable. */
	bool incompressible;
public:

	/*!
	 * \brief Constructor of the class. 
	 */		
	CAdjEulerVariable(void);

	/*!
	 * \overload
	 * \param[in] val_psirho - Value of the adjoint density (initialization value).
	 * \param[in] val_phi - Value of the adjoint velocity (initialization value).
	 * \param[in] val_psie - Value of the adjoint energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjEulerVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	virtual ~CAdjEulerVariable(void);

	/*!
	 * \brief Set the value of the adjoint velocity.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */	
	void SetPhi_Old(double *val_phi);

	/*!
	 * \brief Set the value of theta.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_enthalpy - Value of the enthalpy.
	 */		
	void SetTheta(double val_density, double *val_velocity, double val_enthalpy);

	/*!
	 * \brief Get the value of theta.
	 */		
	double GetTheta(void);

	/*!
	 * \brief Set the value of the force projection vector.
	 * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
	 */		
	void SetForceProj_Vector(double *val_ForceProj_Vector);

	/*!
	 * \brief Set the value of the interior boundary jump vector vector.
	 * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump vector.
	 */		
	void SetIntBoundary_Jump(double *val_IntBoundary_Jump);

	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */		
	double *GetForceProj_Vector(void);

	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */		
	double *GetIntBoundary_Jump(void);
};

/*! 
 * \class CAdjNSVariable
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CAdjNSVariable : public CAdjEulerVariable {	
public:

	/*!
	 * \brief Constructor of the class. 
	 */	
	CAdjNSVariable(void);

	/*!
	 * \overload
	 * \param[in] val_psirho - Value of the adjoint density (initialization value).
	 * \param[in] val_phi - Value of the adjoint velocity (initialization value).
	 * \param[in] val_psie - Value of the adjoint energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAdjNSVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	CAdjNSVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CAdjNSVariable(void);

	/*!
	 * \brief Set the value of the adjoint velocity.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */	
	void SetPhi_Old(double *val_phi);

	/*!
	 * \brief Set the value of theta.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_enthalpy - Value of the enthalpy.
	 */	
	void SetTheta(double val_density, double *val_velocity, double val_enthalpy);

	/*!
	 * \brief Get the value of theta.
	 */		
	double GetTheta(void);

	/*!
	 * \brief Set the value of the force projection vector.
	 * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
	 */
	void SetForceProj_Vector(double *val_ForceProj_Vector);

	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */
	double *GetForceProj_Vector(void);

	/*!
	 * \brief Set the value of the force projection vector on the solution vector.
	 */
	void SetVelSolutionOldDVector(void);

	/*!
	 * \brief Set the value of the force projection vector on the old solution vector.
	 */
	void SetVelSolutionDVector(void);
};

/*! 
 * \class CAdjTurbVariable
 * \brief Main class for defining the variables of the adjoint turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.1.
 */
class CAdjTurbVariable : public CVariable {
public:

	/*!
	 * \brief Constructor of the class. 
	 */		
	CAdjTurbVariable(void);

	/*!
	 * \overload
	 * \param[in] val_psinu_inf - Value of the adjoint turbulence variable at the infinity (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjTurbVariable(double val_psinu_inf, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */		
	~CAdjTurbVariable(void);
};

/*! 
 * \class CLinPotentialVariable
 * \brief Main class for defining the variables of the linearized potential equation.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios.
 * \version 1.1.
 */
class CLinPotentialVariable : public CVariable {
public:	
};

/*! 
 * \class CLinEulerVariable
 * \brief Main class for defining the variables of the linearized Euler's equations.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CLinEulerVariable : public CVariable {
private:
	double *DeltaU;			/*!< \brief Vector of the linearized variables. */
	double *ForceProj_Vector;		/*!< \brief Vector d. */
	double DeltaPressure;	/*!< \brief Linearized pressure variable. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */		
	CLinEulerVariable(void);

	/*!
	 * \overload
	 * \param[in] val_deltarho - Value of the linearized density (initialization value).
	 * \param[in] val_deltavel - Value of the linearized velocity (initialization value).
	 * \param[in] val_deltae - Value of the linearized energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.	
	 * \param[in] config - Definition of the particular problem.
	 */		
	CLinEulerVariable(double val_deltarho, double *val_deltavel, double val_deltae, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the linearized value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	CLinEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CLinEulerVariable(void);

	/*!
	 * \brief Set the value of the linearized velocity.
	 * \param[in] val_deltavel - Value of the linearized velocity.
	 */	
	void SetDeltaVel_Old(double *val_deltavel);

	/*!
	 * \brief Set the value of the force projection vector.
	 * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
	 */		
	void SetForceProj_Vector(double *val_ForceProj_Vector);

	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */		
	double *GetForceProj_Vector(void);

	/*!
	 * \brief Set the value of the linearized pressure.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] Gamma - The ratio of specific heats.
	 */		
	void SetDeltaPressure(double *val_velocity, double Gamma);

	/*!
	 * \brief Get the value of the linearized pressure.
	 * \return Value of the linearized pressure.
	 */		
	double GetDeltaPressure(void);
};

/*! 
 * \class CLinNSVariable
 * \brief Main class for defining the variables of the linearized Navier-Stokes' equations.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CLinNSVariable : public CLinEulerVariable {
public:
};

/*!
 * \class CPlasmaVariable
 * \brief Main class for defining the variables of the Plasma solver.
 * \version 1.1.
 */
class CPlasmaVariable : public CVariable {
protected:
	double *Pressure;				/*!< \brief  partial Pressure of species in flow. */
	double *Temperature;			/*!< \brief temperature of species in flow. */
	double *Temperature_tr; /*!< \brief translation-rotation temperature of species in flow. */

	double **Gradient_Primitive;	/*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */
	double *Velocity2;
	double **Velocity;
	double *SoundSpeed;
	double *Enthalpy;
	double *Enthalpy_Formation;
	double *Source;
	double *Max_Lambda_Inv_MultiSpecies; /*!< \brief Vector of Maximun inviscid eingenvalues for different fluids */
	double *Max_Lambda_Visc_MultiSpecies; /*!< \brief Vector of Maximun viscous eingenvalues for different fluids */

	double **Source_Jacobian;
	unsigned short nSpecies,	/*!< \brief Number of species in plasma and number of fluids to model plasma with. */
	nFluids;
	unsigned short nMonatomics,	/*!< \brief Number of species in plasma and number of fluids to model plasma with. */
	nDiatomics;
	double *LaminarViscosity_MultiSpecies;	/*!< \brief Viscosity of the fluid. */
	double *EddyViscosity_MultiSpecies;	/*!< \brief Viscosity of the fluid. */
	double Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
	double *Kappa;		/*!< \brief Thermal conductivity of the fluid. */
	double *Species_Delta_Time;


public:

	/*!
	 * \brief Constructor of the class.
	 */
	CPlasmaVariable(void);

	/*!
	 * \overload
	 * \param[in] *val_density - Vector of value of densities for all species (initialization value).
	 * \param[in] **val_velocity - Matrix of velocities for all species  [Row: Species, Col: Dim](initialization value).
	 * \param[in] *val_energy -Vector of value of Energies for all species
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Total number of total chemical constitutents in the multi-species flow.
	 * \param[in] val_nFluids - Total number of fluids in the multi-species flow.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPlasmaVariable(double *val_density, double **val_velocity, double *val_energy, unsigned short val_ndim,
			unsigned short val_nvar, unsigned short val_nSpecies, unsigned short  val_nFluids, unsigned short val_nMonatomics,
			unsigned short val_nDiatomics, CConfig *config);

	/*!
	 * \overload
	 * \param[in] *val_density - Vector of value of densities for all species (initialization value).
	 * \param[in] **val_velocity - Matrix of velocities for all species  [Row: Species, Col: Dim](initialization value).
	 * \param[in] *val_energy -Vector of value of total Energies for all species.
	 * \param[in] *val_energy_vib -Vector of values of vibrational energies for all species.
	 * \param[in] *val_enthalpy_formation -Vector of values of chemical enthalpies of formation for all species.
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Total number of total chemical constitutents in the multi-species flow.
	 * \param[in] val_nMonatomics - Number of monatomic species in the multi-species flow.
	 * \param[in] val_nDiatomics - Number of diatomic species in the multi-species flow.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPlasmaVariable(double *val_density, double **val_velocity, double *val_energy, double *val_energy_vib,
			double *val_enthalpy_formation, unsigned short val_ndim,
			unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nMonatomics,
			unsigned short val_nDiatomics, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Total number of total chemical constitutents in the multi-species flow.
	 * \param[in] val_nFluids - Total number of fluids in the multi-species flow.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPlasmaVariable(double *val_solution, unsigned short val_ndim,
				unsigned short val_nvar, unsigned short val_nSpecies, unsigned short  val_nFluids, unsigned short val_nMonatomics, unsigned short val_nDiatomics,CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] val_nSpecies - Total number of species in the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, unsigned short val_nSpecies, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CPlasmaVariable(void);

	/*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	void SetGradient_PrimitiveZero(void);

	/*!
	 * \param[in] val_delta_time - Value of the time step
	 * \param[in] iSpecies - Index of the Species.
	 * \brief Time step of species iSpecies.
	 */
	void SetDelta_Time(double val_delta_time, unsigned short iSpecies);

	/*!
	 * \param[in] iSpecies - Index of the Species.
	 * \brief Time step of species iSpecies.
	 */
	double GetDelta_Time(unsigned short iSpecies);

	/*!
	 * \brief Add <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the primitive variables.
	 */
	void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
	 */
	void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim);

	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double **GetGradient_Primitive(void);

	/*!
	 * \brief Set the value of the pressure.
	 * \param[in] Gamma - Ratio of specific heats
	 */
	void SetPressure(double Gamma);

	/*!
	 * \brief Set the value of the pressure.
	 */
	void SetPressure(double GammaMonatomic, double GammaDiatomic, double *Coord);

	/*!
	 * \brief Set the value of the temperature.
	 * \param[in] Gas_Constant - Vector of Gas_Constants
	 */
	void SetTemperature(double * Gas_Constant);
	
	/*!
	 * \brief Set the value of the translational-rotational temperature
	 * \param[in] Molar_Mass - Value of the molar mass (kg/kmol) for each species
	 * \param[in] GammaMonatomic - Value of the ratio of specific heats for monatomic species.
	 * \param[in] GammaDiatomic - Value of the ratio of specific heats for diatomic species.
	 */
	void SetTemperature_TR(double *Molar_Mass, double GammaMonatomic, double GammaDiatomic);

	/*!
	 * \brief Get the value of the temperature.
	 * \param[in] iSpecies - Species index of desired temperature.
	 */
	double GetTemperature(unsigned short iSpecies);
	
	/*!
	 * \brief Get the value of the trans-rot temperature.
	 * \param[in] iSpecies - Species index of desired temperature.
	 */
	double GetTemperature_TR(unsigned short iSpecies);	 

	/*!
	 * \brief Set the value of the enthalpy.
	 */
	void SetEnthalpy(void);

	/*!
	 * \brief Get the flow pressure.
	 * \return Value of the flow pressure.
	 * \param[in] val_Fluid - Index of the fluid.
	 */
	double GetPressure(unsigned short val_Fluid);

	/*!
	 * \brief Get the species energy of formation.
	 * \return Value of the energy of formation.
	 * \param[in] val_Species - Index of the species.
	 */
	double GetEnthalpyFormation(unsigned short val_Species);
	
	/*!
	 * \brief Get the velocity of the flow.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_Species - Index of the species.
	 * \return Value of the velocity for the dimension <i>val_dim</i>.
	 */
	double GetVelocity(unsigned short val_dim, unsigned short val_Species);

	/*!
	 * \brief Get the square of velocity of that particular flow.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the velocity for the dimension <i>val_dim</i>.
	 */
	double GetVelocity2(unsigned short val_Fluid);

	/*! overload
	 * \brief Get the projected velocity in a unitary vector direction.
	 * \param[in] val_vector - Direction of projection.
	 * \param[in] val_var - Variable index desired.
	 * \return Value of the projected velocity.
	 */
	double GetProjVel(double *val_vector, unsigned short val_var);

	/*! overload
	 * \brief Get the projected velocity in a unitary vector direction.
	 * \param[in] val_vector - Direction of projection.
	 * \param[in] val_Species - Index of the desired species.
	 * \param[in] val_nDiatomics - Number of diatomic species in the multi-species flow.
	 * \return Value of the projected velocity.
	 */
	double GetProjVel(double *val_vector, unsigned short val_Species, unsigned short val_nDiatomics);

	/*!
	 * \brief Get the speed of the sound.
	 * \return Value of speed of the sound pf val_Fluid
	 */
	double GetSoundSpeed(unsigned short val_Species);

	/*!
	 * \brief Get the Density of species
	 * \return Value of speed of the sound pf val_Fluid
	 */
	double GetDensity(unsigned short val_Species);

	/*!
	 * \brief Get the square of velocity of fluid val_var
	 * \return Value of speed of the sound of val_var fluid
	 */

	void SetVelocity2(void);

	/*!
	 * \brief Calculate the speed of sound of all the fluids
	 * \param[in] Gamma - Ratio of specific heats
	 * \return Value of speed of the sound of all the fluids
	 */
	void SetSoundSpeed(double Gamma);

	/*!
	 * \brief Calculate the speed of sound of all the fluids
	 * \return Value of speed of the sound of all the fluids
	 */
	void SetSoundSpeed(double GammaMonatomic, double GammaDiatomic);

	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetLaminarViscosity(CConfig *config);

	/*!
	 * \brief Get the contribution of source term to plasma equations
	 * \param[in] val_var - Index of the variable.
	 * \return Value of source added to plasma equation val_var
	 */
	double GetSource(unsigned short val_var);

	/*!
	 * \brief Set the contribution of source term to plasma equations
	 */
	void SetSource(double * res);

	/*!
	 * \brief Get the contribution of source term to plasma equations
	 * \param[in] val_ivar - Index for i variable number.
	 * \param[in] val_jvar - Index for j variable number.
	 * \return Value of source added to plasma equation val_var
	 */
	double GetSource_Jacobian(unsigned short val_ivar,unsigned short val_jvar);

	/*!
	 * \brief Set the contribution of source term to plasma equations
	 */
	void SetSourceJacobian(double ** src_jacobian);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_iSpecies - Index of the species to set the value of max lambda.
	 */
	void SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_iSpecies);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_iSpecies - Index of the species to set the value of max lambda.
	 */
	void SetMax_Lambda_Visc(double val_max_lambda, unsigned short val_iSpecies);

	/*!
	 * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] iSpecies - ID of the species in multispecies solver
	 * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	double GetMax_Lambda_Inv(unsigned short iSpecies);

	/*!
	 * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] iSpecies - ID of the species in multispecies solver
	 * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	double GetMax_Lambda_Visc(unsigned short iSpecies);

	/*!
	 * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] iSpecies - Value of fluid to which the eigenvalue belongs
	 */
	void AddMax_Lambda_Inv(double val_max_lambda, unsigned short iSpecies);

	/*!
	 * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] iSpecies - Value of fluid to which the eigenvalue belongs
	 */
	void AddMax_Lambda_Visc(double val_max_lambda, unsigned short iSpecies);

	/*!
	 * \brief Get the laminar viscosity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetLaminarViscosity(unsigned short iSpecies);

	/*!
	 * \brief Set the thermal coefficient.
	 * \param[in] Gamma - Ratio of specific heats
	 * \param[in] Gas_Constant - The Gas Constant of various species
	 */
	void SetThermalCoeff(double Gamma, double *Gas_Constant);

	/*!
	 * \brief Get the laminar viscosity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetEddyViscosity(unsigned short iSpecies);

	/*!
	 * \brief Set the velocity vector from the old solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] iSpecies - Index of the species to set the old velocity.
	 */
	void SetVelocity_Old(double *val_velocity, unsigned short iSpecies);

	/*!
	 * \brief A virtual member.
	 * \param[in] Temperature_Wall - Value of the Temperature at the wall
	 */
	void SetWallTemperature(double* Temperature_Wall);

};


/*! 
 * \class CLevelSetVariable
 * \brief Main class for defining the variables of the Level Set.
 * \ingroup LevelSet_Model
 * \author F. Palacios.
 * \version 1.1.
 */
class CLevelSetVariable : public CVariable {
protected:
	double DiffLevelSet;		/*!< \brief  Diff LevelSet Distribution. */

public:
	/*!
	 * \brief Constructor of the class. 
	 */	
	CLevelSetVariable(void);

	/*!
	 * \overload
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CLevelSetVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_levelset - Level set variable value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CLevelSetVariable(double val_levelset, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
	
	/*!
	 * \brief Set the value of the diff level set (value-target).
	 * \param[in] val_difflevelset - Value of the diff level set (value-target).
	 */	
	void SetDiffLevelSet(double val_difflevelset);

	/*!
	 * \brief Get the value of theta.
	 */		
	double GetDiffLevelSet(void);
	
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CLevelSetVariable(void);

};

/*! 
 * \class CAdjLevelSetVariable
 * \brief Main class for defining the variables of the Level Set.
 * \ingroup LevelSet_Model
 * \author F. Palacios.
 * \version 1.1.
 */
class CAdjLevelSetVariable : public CVariable {
public:
	/*!
	 * \brief Constructor of the class. 
	 */	
	CAdjLevelSetVariable(void);
	
	/*!
	 * \overload
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjLevelSetVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
	
	/*!
	 * \overload
	 * \param[in] val_levelset - Level set variable value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAdjLevelSetVariable(double val_levelset, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
	
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjLevelSetVariable(void);
	
};

/*! 
 * \class CAdjPlasmaVariable
 * \brief Main class for defining the variables of the adjoint Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.1.
 */
class CAdjPlasmaVariable : public CVariable {
protected:
	double *Psi;		/*!< \brief Vector of the adjoint variables. */
	double *ForceProj_Vector;	/*!< \brief Vector d. */
	double *IntBoundary_Jump;	/*!< \brief Interior boundary jump vector. */
	double Theta;		/*!< \brief Theta variable. */
	bool incompressible;
public:
	
	/*!
	 * \brief Constructor of the class. 
	 */		
	CAdjPlasmaVariable(void);
	
	/*!
	 * \overload
	 * \param[in] val_psirho - Value of the adjoint density (initialization value).
	 * \param[in] val_phi - Value of the adjoint velocity (initialization value).
	 * \param[in] val_psie - Value of the adjoint energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjPlasmaVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
	
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
	
	/*!
	 * \brief Destructor of the class. 
	 */	
	virtual ~CAdjPlasmaVariable(void);
	
	/*!
	 * \brief Set the value of the adjoint velocity.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */	
	void SetPhi_Old(double *val_phi);
	
	/*!
	 * \brief Set the value of theta.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_enthalpy - Value of the enthalpy.
	 */		
	void SetTheta(double val_density, double *val_velocity, double val_enthalpy);
	
	/*!
	 * \brief Get the value of theta.
	 */		
	double GetTheta(void);
	
	/*!
	 * \brief Set the value of the force projection vector.
	 * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
	 */		
	void SetForceProj_Vector(double *val_ForceProj_Vector);
	
	/*!
	 * \brief Set the value of the interior boundary jump vector vector.
	 * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump vector.
	 */		
	void SetIntBoundary_Jump(double *val_IntBoundary_Jump);
	
	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */		
	double *GetForceProj_Vector(void);
	
	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */		
	double *GetIntBoundary_Jump(void);
};

/*! 
 * \class CTemplateVariable
 * \brief Main class for defining the variables of the potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios.
 * \version 1.1.
 */
class CTemplateVariable : public CVariable {
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CTemplateVariable(void);

	/*!
	 * \overload
	 * \param[in] val_potential - Value of the potential solution (initialization value).		 
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CTemplateVariable(double val_potential, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CTemplateVariable(void);
};

#include "variable_structure.inl"
