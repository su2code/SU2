/*!
 * \file variable_structure.hpp
 * \brief Headers of the main subroutines for storing all the variables for 
 *        each kind of governing equation (direct, adjoint and linearized).
 *        The subroutines and functions are in the <i>variable_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 3.2.8.1 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../Common/include/config_structure.hpp"
#include "fluid_model.hpp"


using namespace std;

/*! 
 * \class CVariable
 * \brief Main class for defining the variables.
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CVariable {
protected:

	double *Solution,		/*!< \brief Solution of the problem. */
	*Solution_Old;			/*!< \brief Old solution of the problem R-K. */
  bool Non_Physical;			/*!< \brief Non-physical points in the solution (force first order). */
	double *Solution_time_n,	/*!< \brief Solution of the problem at time n for dual-time stepping technique. */
	*Solution_time_n1;			/*!< \brief Solution of the problem at time n-1 for dual-time stepping technique. */
	double **Gradient;		/*!< \brief Gradient of the solution of the problem. */ 
	double *Limiter;				/*!< \brief Limiter of the solution of the problem. */
	double *Solution_Max;		/*!< \brief Max solution for limiter computation. */
	double *Solution_Min;		/*!< \brief Min solution for limiter computation. */
	double AuxVar;			/*!< \brief Auxiliar variable for gradient computation. */
	double *Grad_AuxVar;	/*!< \brief Gradient of the auxiliar variable. */
	double Delta_Time;	/*!< \brief Time step. */
	double Max_Lambda,	/*!< \brief Maximun eingenvalue. */
	Max_Lambda_Inv,		/*!< \brief Maximun inviscid eingenvalue. */
	Max_Lambda_Visc,	/*!< \brief Maximun viscous eingenvalue. */
	Lambda;				/*!< \brief Value of the eingenvalue. */
	double Sensor;	/*!< \brief Pressure sensor for high order central scheme. */
	double *Undivided_Laplacian;	/*!< \brief Undivided laplacian of the solution. */
	double *Res_TruncError,	/*!< \brief Truncation error for multigrid cycle. */
	*Residual_Old,		/*!< \brief Auxiliar structure for residual smoothing. */
	*Residual_Sum;		/*!< \brief Auxiliar structure for residual smoothing. */
	static unsigned short nDim;		/*!< \brief Number of dimension of the problem. */
	unsigned short nVar;		/*!< \brief Number of variables of the problem, 
													 note that this variable cannnot be static, it is possible to 
													 have different number of nVar in the same problem. */
  unsigned short nPrimVar, nPrimVarGrad;		/*!< \brief Number of variables of the problem,
                                             note that this variable cannnot be static, it is possible to
                                             have different number of nVar in the same problem. */
  unsigned short nSecondaryVar, nSecondaryVarGrad;		/*!< \brief Number of variables of the problem,
                                             note that this variable cannnot be static, it is possible to
                                             have different number of nVar in the same problem. */
  
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
	void SetSolution(double *val_solution);

	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
	 */
	void SetSolution(unsigned short val_var, double val_solution);

  /*!
   * \brief Set the value of the non-physical point.
   * \param[in] val_value - identification of the non-physical point.
   */
  void SetNon_Physical(bool val_value);
  
  /*!
   * \brief Get the value of the non-physical point.
   * \return Value of the Non-physical point.
   */
  double GetNon_Physical(void);
  
	/*!
	 * \brief Get the solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the solution for the index <i>val_var</i>.
	 */
	double GetSolution(unsigned short val_var);

	/*!
	 * \brief Get the old solution of the problem (Runge-Kutta method)
	 * \param[in] val_var - Index of the variable.
	 * \return Pointer to the old solution vector.
	 */
	double GetSolution_Old(unsigned short val_var);

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
	 * \brief Set to zero the velocity components of the solution.
	 */
	void SetVelSolutionZero(void);

  /*!
	 * \brief Specify a vector to set the velocity components of the solution.
   * \param[in] val_vector - Pointer to the vector.
	 */
	void SetVelSolutionVector(double *val_vector);
  
	/*!
	 * \brief Set to zero velocity components of the solution.
	 */
	void SetVelSolutionOldZero(void);

  /*!
	 * \brief Specify a vector to set the velocity components of the old solution.
   * \param[in] val_vector - Pointer to the vector.
	 */
	void SetVelSolutionOldVector(double *val_vector);
  
	/*!
	 * \brief Set to zero the solution.
	 */	
	void SetSolutionZero(void);
  
  /*!
	 * \brief Set to zero a particular solution.
	 */
  void SetSolutionZero(unsigned short val_var);

	/*!
	 * \brief Add a value to the solution.
	 * \param[in] val_var - Number of the variable.
	 * \param[in] val_solution - Value that we want to add to the solution.
	 */
	void AddSolution(unsigned short val_var, double val_solution);

  /*!
	 * \brief Add a value to the solution, clipping the values.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution change.
   * \param[in] lowerlimit - Lower value.
   * \param[in] upperlimit - Upper value.
	 */
	void AddClippedSolution(unsigned short val_var, double val_solution,
                          double lowerlimit, double upperlimit);
  
	/*!
	 * \brief Update the variables using a conservative format.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution change.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_density_old - Value of the old density.
	 */
	void AddConservativeSolution(unsigned short val_var, double val_solution,
			double val_density, double val_density_old, double lowerlimit,
			double upperlimit);

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
	 * \brief Set the value of the old residual.
	 * \param[in] val_residual_old - Pointer to the residual vector.
	 */
	void SetResidual_Old(double *val_residual_old);

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
	 * \brief Set the velocity of the truncation error to zero.
	 */
	virtual void SetVel_ResTruncError_Zero(unsigned short iSpecies);

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
	 * \param[in] val_residual - Pointer to the summed residual.
	 */	
	void GetResidual_Sum(double *val_residual);

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
	void AddRes_TruncError(double *val_truncation_error);

	/*!
	 * \brief Subtract a value to the truncation error.
	 * \param[in] val_truncation_error - Value that we want to subtract to the truncation error.
	 */		
	void SubtractRes_TruncError(double *val_truncation_error);

	/*!
	 * \brief Set the truncation error to zero.
	 */		
	void SetRes_TruncErrorZero(void);
  
  /*!
	 * \brief Set the truncation error to zero.
	 */
	void SetVal_ResTruncError_Zero(unsigned short val_var);

	/*!
	 * \brief Set the velocity of the truncation error to zero.
	 */		
	void SetVel_ResTruncError_Zero(void);
  
  /*!
	 * \brief Set the velocity of the truncation error to zero.
	 */
	void SetEnergy_ResTruncError_Zero(void);

	/*!
	 * \brief Get the truncation error.
	 * \return Pointer to the truncation error.
	 */	
	double *GetResTruncError(void);

	/*!
	 * \brief Get the truncation error.
	 * \param[in] val_trunc_error - Pointer to the truncation error.
	 */	
	void GetResTruncError(double *val_trunc_error);

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
	 * \brief Set the value of the limiter.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_limiter - Value of the limiter for the index <i>val_var</i>.
	 */
	virtual void SetLimiterPrimitive(unsigned short val_species, unsigned short val_var, double val_limiter);
  
  /*!
	 * \brief Set the value of the limiter.
   * \param[in] val_species - Value of the limiter for the index <i>val_var</i>.
	 * \param[in] val_var - Index of the variable.
	 */
  virtual double GetLimiterPrimitive(unsigned short val_species, unsigned short val_var);
	
	/*!
	 * \brief Set the value of the max solution.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the max solution for the index <i>val_var</i>.
	 */
	void SetSolution_Max(unsigned short val_var, double val_solution);
	
	/*!
	 * \brief Set the value of the min solution.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the min solution for the index <i>val_var</i>.
	 */
	void SetSolution_Min(unsigned short val_var, double val_solution);

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
	 * \brief Get the value of the min solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the min solution for the variable <i>val_var</i>.
	 */
	double GetSolution_Max(unsigned short val_var);
	
	/*!
	 * \brief Get the value of the min solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the min solution for the variable <i>val_var</i>.
	 */
	double GetSolution_Min(unsigned short val_var);

	/*!
	 * \brief Get the value of the preconditioner Beta.
	 * \return Value of the low Mach preconditioner variable Beta
	 */
	virtual double GetPreconditioner_Beta();

	/*!
	 * \brief Set the value of the preconditioner Beta.
	 * \param[in] val_Beta - Value of the low Mach preconditioner variable Beta
	 */
	virtual void SetPreconditioner_Beta(double val_Beta);

       /*!
	 * \brief Get the value of the wind gust
	 * \return Value of the wind gust
	 */
	virtual double* GetWindGust();
    
	/*!
	 * \brief Set the value of the wind gust
	 * \param[in] val_WindGust - Value of the wind gust
	 */
	virtual void SetWindGust(double* val_WindGust);
    
    /*!
	 * \brief Get the value of the derivatives of the wind gust
	 * \return Value of the derivatives of the wind gust
	 */
	virtual double* GetWindGustDer();
    
	/*!
	 * \brief Set the value of the derivatives of the wind gust
	 * \param[in] val_WindGust - Value of the derivatives of the wind gust
	 */
	virtual void SetWindGustDer(double* val_WindGust);
    
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
	 * \param[in] val_species - Value of the species index to set the maximum eigenvalue.
	 */
	virtual void SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_species);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	void SetMax_Lambda_Visc(double val_max_lambda);

	/*!
	 * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_species - Index of the species to set the maximum eigenvalue of the viscous terms.
	 */
	virtual void SetMax_Lambda_Visc(double val_max_lambda, unsigned short val_species);

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
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda - Value of the spectral radius.
	 * \param[in] val_iSpecies -Index of species
	 */
	virtual void SetLambda(double val_lambda, unsigned short val_iSpecies);

	/*!
	 * \brief Add the value of the spectral radius.
	 * \param[in] val_lambda - Value of the spectral radius.
	 */
	void AddLambda(double val_lambda);

	/*!
	 * \brief Add the value of the spectral radius.
	 * \param[in] val_iSpecies -Index of species
	 * \param[in] val_lambda - Value of the spectral radius.
	 */
	virtual void AddLambda(double val_lambda, unsigned short val_iSpecies);

	/*!
	 * \brief Get the value of the spectral radius.
	 * \return Value of the spectral radius.
	 */
	double GetLambda(void);

	/*!
	 * \brief Get the value of the spectral radius.
	 * \param[in] val_iSpecies -Index of species
	 * \return Value of the spectral radius.
	 */
	virtual double GetLambda(unsigned short val_iSpecies);

	/*!
	 * \brief Set pressure sensor.
	 * \param[in] val_sensor - Value of the pressure sensor.
	 */
	void SetSensor(double val_sensor);

	/*!
	 * \brief Set pressure sensor.
	 * \param[in] val_sensor - Value of the pressure sensor.
	 * \param[in] val_sensor - Index of the Species.
	 */
	virtual void SetSensor(double val_sensor, unsigned short iSpecies);

	/*!
	 * \brief Get the pressure sensor.
	 * \return Value of the pressure sensor.
	 */	
	double GetSensor(void);

	/*!
	 * \brief Get the pressure sensor.
	 * \param[in] iSpecies - index of species
	 * \return Value of the pressure sensor.
	 */
	virtual double GetSensor(unsigned short iSpecies);

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
	double *GetUndivided_Laplacian(void);

	/*!
	 * \brief Get the undivided laplacian of the solution.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \return Value of the undivided laplacian vector.
	 */
	double GetUndivided_Laplacian(unsigned short val_var);

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
   * \param[in] val_Species - Index of species s.
	 * \return Value of the mass fraction of species s.
	 */
	virtual double GetMassFraction(unsigned short val_Species);

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
	 * \return Pointer to the objective function source.
	 */
	virtual double *GetObjFuncSource(void);

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
	 * \return Value of the eddy viscosity.
	 */
	virtual double GetEddyViscosityInc(void);

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
	 * \brief A virtual member.
	 * \return Value of the flow pressure.
	 */
	virtual double GetPressureInc(void);

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
	 * \param[in] val_species - Index of the desired species.
	 * \return Value of the projected velocity.
	 */		
	virtual double GetProjVel(double *val_vector, unsigned short val_species);

	/*!
	 * \brief A virtual member.
	 * \return Value of the sound speed.
	 */		
	virtual double GetSoundSpeed(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density for the incompressible flow.
	 */		
	virtual double GetDensityInc(void);

  /*!
	 * \brief A virtual member.
	 * \return Value of the levelset for the freesurface flows.
	 */
	virtual double GetLevelSet(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the distance for the freesurface flows.
	 */
	virtual double GetDistance(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the beta for the incompressible flow.
	 */		
	virtual double GetBetaInc2(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the temperature.
	 */		
	virtual double GetTemperature(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the vibrational-electronic temperature.
	 */
	virtual double GetTemperature_ve(void);
  
  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (trans.-rot.).
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  virtual double GetRhoCv_tr(void);
  
  /*!
   * \brief A virtual member -- Get the mixture specific heat at constant volume (vib.-el.).
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  virtual double GetRhoCv_ve(void);
  
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
	 * \return Norm 2 of the velocity vector of Fluid val_species.
	 */	
	virtual double GetVelocity2(unsigned short val_species);

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
   * \return Value of the species diffusion coefficient.
   */
  virtual double* GetDiffusionCoeff(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the thermal conductivity (translational/rotational)
   */
  virtual double GetThermalConductivity(void);

  /*!
   * \brief A virtual member.
   * \return Value of the specific heat at constant P
   */
   virtual double GetSpecificHeatCp(void);

  /*!
   * \brief A virtual member.
   * \return Value of the thermal conductivity (vibrational)
   */
  virtual double GetThermalConductivity_ve(void);

	/*!
	 * \brief A virtual member.
	 * \return Sets separation intermittency
	 */
	virtual void SetGammaSep(double gamma_sep);

	/*!
	 * \brief A virtual member.
	 * \return Sets separation intermittency
	 */
	virtual void SetGammaEff(void);

	/*!
	 * \brief A virtual member.
	 * \return Returns intermittency
	 */
	virtual double GetIntermittency();

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */		
	virtual double *GetVorticity(void);

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
	 * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
	 */
	virtual void SetObjFuncSource(double *val_SetObjFuncSource);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump.
	 */		
	virtual void SetIntBoundary_Jump(double *val_IntBoundary_Jump);

	/*!
	 * \brief A virtual member.
	 * \param[in] eddy_visc - Value of the eddy viscosity.
	 */		
	virtual void SetEddyViscosity(double eddy_visc);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] eddy_visc - Value of the eddy viscosity.
	 */
	virtual void SetEddyViscosityInc(double eddy_visc);

	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetEnthalpy(void);
	
	/*!
	 * \brief A virtual member.
	 */		
	virtual bool SetPrimVar_Compressible(CConfig *config);

	/*!
	 * \brief A virtual member.
	 */
	 virtual bool SetPrimVar_Compressible(CFluidModel *FluidModel);

	/*!
     * \brief A virtual member.
	 */
     virtual void SetSecondaryVar_Compressible(CFluidModel *FluidModel);

  /*!
	 * \brief A virtual member.
	 */
  virtual bool Cons2PrimVar(CConfig *config, double *U, double *V,
                            double *dPdU, double *dTdU,
                            double *dTvedU);
  /*!
	 * \brief A virtual member.
	 */
  virtual void Prim2ConsVar(CConfig *config, double *V, double *U);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_Compressible(double SharpEdge_Distance, bool check, CConfig *config);
	
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_Incompressible(double SharpEdge_Distance, bool check, CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_FreeSurface(double SharpEdge_Distance, bool check, CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 */		
	virtual bool SetPrimVar_Compressible(double eddy_visc, double turb_ke, CConfig *config);
	
	/*!
	 * \brief A virtual member.
	 */		
	virtual bool SetPrimVar_Compressible(double eddy_visc, double turb_ke, CFluidModel *FluidModel);

	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_Incompressible(double Density_Inf, CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_FreeSurface(CConfig *config);
	
	/*!
	 * \brief A virtual member.
	 */		
	virtual bool SetPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf, double eddy_visc, double turb_ke, CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_FreeSurface(double eddy_visc, double turb_ke, CConfig *config);
	
	/*!
	 * \brief A virtual member.
	 */
	virtual double GetPrimitive(unsigned short val_var);
  
  /*!
	 * \brief A virtual member.
	 */
  virtual void SetPrimitive(unsigned short val_var, double val_prim);
  
  /*!
	 * \brief A virtual member.
	 */
  virtual void SetPrimitive(double *val_prim);

	/*!
	 * \brief A virtual member.
	 */
	virtual double *GetPrimitive(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual double GetSecondary(unsigned short val_var);
  
  /*!
	 * \brief A virtual member.
	 */
  virtual void SetSecondary(unsigned short val_var, double val_secondary);
  
    /*!
	 * \brief A virtual member.
	 */
    virtual void SetSecondary(double *val_secondary);

    /*!
 	 * \brief A virtual member.
 	 */
     virtual void SetdPdrho_e(double dPdrho_e);

     /*!
   	 * \brief A virtual member.
   	 */
     virtual void SetdPde_rho(double dPde_rho);

    /*!
  	 * \brief A virtual member.
  	 */
     virtual void SetdTdrho_e(double dTdrho_e);

    /*!
     * \brief A virtual member.
     */
     virtual void SetdTde_rho(double dTde_rho);

     /*!
   	 * \brief A virtual member.
   	 */
      virtual void Setdmudrho_T(double dmudrho_T);

     /*!
      * \brief A virtual member.
      */
      virtual void SetdmudT_rho(double dmudT_rho);

    /*!
     * \brief A virtual member.
     */
     virtual void Setdktdrho_T(double dktdrho_T);

    /*!
     * \brief A virtual member.
     */
     virtual void SetdktdT_rho(double dktdT_rho);

	/*!
	 * \brief A virtual member.
	 */
	 virtual double *GetSecondary(void);
	
	/*!
	 * \brief A virtual member.
	 */		
	virtual void SetDensityInc(double val_density);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual void SetPressureInc(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual void SetVelocityInc(void);

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
	 */
	virtual bool SetPressure(double Gamma);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] config
	 */
	virtual bool SetPressure(CConfig *config);

	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetPressure(double Gamma, double turb_ke);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetPressure(void);
  
  /*!
   * \brief Calculates vib.-el. energy per mass, \f$e^{vib-el}_s\f$, for input species (not including KE)
   */
  virtual double CalcEve(double *V, CConfig *config, unsigned short val_Species);
  
  /*!
   * \brief Calculates enthalpy per mass, \f$h_s\f$, for input species (not including KE)
   */
  virtual double CalcHs(double *V, CConfig *config, unsigned short val_Species);
  
  /*!
   * \brief Calculates enthalpy per mass, \f$Cv_s\f$, for input species (not including KE)
   */
  virtual double CalcCvve(double val_Tve, CConfig *config, unsigned short val_Species);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration settings
	 */
	virtual void CalcdPdU(double *V, CConfig *config, double *dPdU);
  
  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  virtual void CalcdTdU(double *V, CConfig *config, double *dTdU);
  
  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  virtual void CalcdTvedU(double *V, CConfig *config, double *dTdU);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual double *GetdPdU(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual double *GetdTdU(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual double *GetdTvedU(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetDensity(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] Gamma - Ratio of Specific heats
	 */		
	virtual void SetDeltaPressure(double *val_velocity, double Gamma);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Ratio of specific heats.
	 */		
	virtual bool SetSoundSpeed(double Gamma);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */		
	virtual bool SetSoundSpeed(CConfig *config);

	/*!
	 * \brief A virtual member.
	 */		
	virtual bool SetSoundSpeed(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] Gas_Constant - Value of the Gas Constant
	 */		
	virtual bool SetTemperature(double Gas_Constant);
  
  /*!
	 * \brief Sets the vibrational electronic temperature of the flow.
	 * \return Value of the temperature of the flow.
	 */
  virtual bool SetTemperature_ve(double val_Tve);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	virtual bool SetTemperature(CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */	
	virtual void SetPrimitive(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	virtual void SetPrimitive(CConfig *config, double *Coord);
	
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
	 * \brief Set the thermal coefficient.
	 * \param[in] config - Configuration parameters.
	 */
	virtual void SetThermalCoeff(CConfig *config);

	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVelocity(void);
  
	/*!
	 * \brief A virtual member.
	 */
  virtual void SetStress(unsigned short iVar, unsigned short jVar, double val_stress);
  
	/*!
	 * \brief A virtual member.
   
	 */
  virtual double **GetStress(void);
  
	/*!
	 * \brief A virtual member.
	 */
  virtual void SetVonMises_Stress(double val_stress);
  
	/*!
	 * \brief A virtual member.
   
	 */
  virtual double GetVonMises_Stress(void);
  
  /*!
	 * \brief A virtual member.
	 */
  virtual void SetFlow_Pressure(double val_pressure);
  
	/*!
	 * \brief A virtual member.
   
	 */
  virtual double GetFlow_Pressure(void);

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
	 */
	virtual void SetVelocityInc_Old(double *val_velocity);

	/*!
	 * \brief A virtual member.
	 * \param[in] laminarViscosity
	 */	
	virtual void SetLaminarViscosity(double laminarViscosity);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetLaminarViscosity(CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_laminar_viscosity_inc - Value of the laminar viscosity (incompressible flows).
	 */		
	virtual void SetLaminarViscosityInc(double val_laminar_viscosity_inc);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetThermalConductivity(double thermalConductivity);
  
  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetThermalConductivity(CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetSpecificHeatCp(double Cp);

	/*!
	 * \brief A virtual member.
	 */		
	virtual bool SetVorticity(bool val_limiter);

	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetStrainMag(bool val_limiter);

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
	virtual void SetGradient_PrimitiveZero(unsigned short val_primvar);

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
	 * \return Value of the primitive variables gradient.
	 */
	virtual double GetLimiter_Primitive(unsigned short val_var);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

  /*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetLimiter_Primitive(unsigned short val_var, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double **GetGradient_Primitive(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double *GetLimiter_Primitive(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual void SetGradient_SecondaryZero(unsigned short val_secondaryvar);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the Secondary variables.
	 */
	virtual void AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the Secondary variables.
	 */
	virtual void SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the Secondary variables gradient.
	 */
	virtual double GetGradient_Secondary(unsigned short val_var, unsigned short val_dim);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the Secondary variables gradient.
	 */
	virtual double GetLimiter_Secondary(unsigned short val_var);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, double val_value);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetLimiter_Secondary(unsigned short val_var, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the Secondary variables gradient.
	 */
	virtual double **GetGradient_Secondary(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the Secondary variables gradient.
	 */
	virtual double *GetLimiter_Secondary(void);
  
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
	virtual double GetCrossDiff(void) { return 0.0; };

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

	/*!
	 * \brief A virtual member.
	 * \param[in] val_source - Value of the time spectral source.
	 */
	virtual void SetTimeSpectral_Source(unsigned short val_var, double val_source);

	/*!
	 * \brief A virtual member.
	 */
	virtual double GetTimeSpectral_Source(unsigned short val_var);

	/*!
	 * \brief Set the Eddy Viscosity Sensitivity of the problem.
	 * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
	 */
	virtual void SetEddyViscSens(double *val_EddyViscSens, unsigned short numTotalVar);

	/*!
	 * \brief Get the Eddy Viscosity Sensitivity of the problem.
	 * \return Pointer to the Eddy Viscosity Sensitivity.
	 */
	virtual double *GetEddyViscSens(void);
  
  /*!
   * \brief A virtual member.  Retrieves index of species densities in the TNE2 solver.
   */
  virtual unsigned short GetRhosIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 *  iRho_s
	 */
  virtual unsigned short GetRhoIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 *  iRho_s
	 */
  virtual unsigned short GetPIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  virtual unsigned short GetTIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  virtual unsigned short GetTveIndex(void);
  
  /*!
	 * \brief Retrieves the value of the velocity index in the primitive variable vector.
	 * iRho*u
	 */
  virtual unsigned short GetVelIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  virtual unsigned short GetHIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  virtual unsigned short GetAIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  virtual unsigned short GetRhoCvtrIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  virtual unsigned short GetRhoCvveIndex(void);
  
  /*!
	 * \brief A virtual member. Set the direct solution for the adjoint solver.
	 * \param[in] val_solution_direct - Value of the direct solution.
	 */
	virtual void SetSolution_Direct(double *val_solution_direct);
  
	/*!
	 * \brief A virtual member. Get the direct solution for the adjoint solver.
	 * \return Pointer to the direct solution vector.
	 */
	virtual double *GetSolution_Direct(void);
  
};

/*!
 * \class CBaselineVariable
 * \brief Main class for defining the variables of a baseline solution from a restart file (for output).
 * \author F. Palacios, T. Economon.
 * \version 3.2.8.1 "eagle"
 */
class CBaselineVariable : public CVariable {
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CBaselineVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CBaselineVariable(double *val_solution, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CBaselineVariable(void);
  
};

/*! 
 * \class CPotentialVariable
 * \brief Main class for defining the variables of the potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CPotentialVariable : public CVariable {
	double *Charge_Density;
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CPotentialVariable(void);

	/*!
	 * \overload
	 * \param[in] val_potential - Value of the potential solution (initialization value).		 
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CPotentialVariable(double val_potential, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
	 * \brief Destructor of the class.
	 */
	~CPotentialVariable(void);
  
	/*!
	 * \brief A virtual member.
	 */
	double* GetChargeDensity();

	/*!
	 * \brief A virtual member.
	 * \param[in] positive_charge - Mass density of positive charge.
	 * \param[in] negative_charge - Mass density of negative charge.
	 */
	void SetChargeDensity(double positive_charge, double negative_charge);

};

/*! 
 * \class CWaveVariable
 * \brief Main class for defining the variables of the wave equation solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CWaveVariable : public CVariable {
protected:
	double *Solution_Direct;  /*!< \brief Direct solution container for use in the adjoint wave solver. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CWaveVariable(void);

	/*!
	 * \overload
	 * \param[in] val_wave - Values of the wave solution (initialization value).		 
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CWaveVariable(double *val_wave, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CWaveVariable(void);

	/*!
	 * \brief Set the direct solution for the adjoint solver.
	 * \param[in] val_solution_direct - Value of the direct solution.
	 */
	void SetSolution_Direct(double *val_solution_direct);

	/*!
	 * \brief Get the direct solution for the adjoint solver.
	 * \return Pointer to the direct solution vector.
	 */
	double *GetSolution_Direct(void);

};

/*! 
 * \class CHeatVariable
 * \brief Main class for defining the variables of the Heat equation solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CHeatVariable : public CVariable {
protected:
	double *Solution_Direct;  /*!< \brief Direct solution container for use in the adjoint Heat solver. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CHeatVariable(void);

	/*!
	 * \overload
	 * \param[in] val_Heat - Values of the Heat solution (initialization value).		 
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CHeatVariable(double *val_Heat, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CHeatVariable(void);

	/*!
	 * \brief Set the direct solution for the adjoint solver.
	 * \param[in] val_solution_direct - Value of the direct solution.
	 */
	void SetSolution_Direct(double *val_solution_direct);

	/*!
	 * \brief Get the direct solution for the adjoint solver.
	 * \return Pointer to the direct solution vector.
	 */
	double *GetSolution_Direct(void);

};

/*! 
 * \class CFEAVariable
 * \brief Main class for defining the variables of the FEA equation solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CFEAVariable : public CVariable {
protected:
	double Flow_Pressure;	/*!< \brief Pressure of the fluid. */
  double **Stress;  /*!< \brief Stress tensor. */
  double VonMises_Stress; /*!< \brief Von Mises stress. */
  
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CFEAVariable(void);

	/*!
	 * \overload
	 * \param[in] val_fea - Values of the fea solution (initialization value).		 
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CFEAVariable(double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CFEAVariable(void);
  
  /*!
	 * \brief Set the value of the stress.
   * \param[in] iVar - i index.
	 * \param[in] jVar - j index.
	 * \param[in] val_stress - Value of the stress.
	 */
  void SetStress(unsigned short iVar, unsigned short jVar, double val_stress);
  
  /*!
	 * \brief Get the value of the stress.
   * \return Value of the stress.
	 */
  double **GetStress(void);

  /*!
	 * \brief Set the value of the Von Mises stress.
	 * \param[in] val_stress - Value of the Von Mises stress.
	 */
  void SetVonMises_Stress(double val_stress);
  
  /*!
	 * \brief Get the value of the Von Mises stress.
   * \return Value of the Von Mises stress.
	 */
  double GetVonMises_Stress(void);
  
  /*!
	 * \brief Set the value of the Von Mises stress.
	 * \param[in] val_stress - Value of the Von Mises stress.
	 */
  void SetFlow_Pressure(double val_pressure);
  
  /*!
	 * \brief Get the value of the Von Mises stress.
   * \return Value of the Von Mises stress.
	 */
  double GetFlow_Pressure(void);

};

/*! 
 * \class CEulerVariable
 * \brief Main class for defining the variables of the Euler's solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CEulerVariable : public CVariable {
protected:
	double Velocity2;			/*!< \brief Square of the velocity vector. */
	double *TS_Source;		/*!< \brief Time spectral source term. */
	double Precond_Beta;	/*!< \brief Low Mach number preconditioner value, Beta. */
  double *WindGust;           /*! < \brief Wind gust value */
  double *WindGustDer;        /*! < \brief Wind gust derivatives value */

	/*--- Primitive variable definition ---*/
  
	double *Primitive;	/*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,c) in compressible flows. */
	double **Gradient_Primitive;	/*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */ 
  double *Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T,vx,vy,vz,P,rho). */ 

  /*--- Secondary variable definition ---*/
  
	double *Secondary;	/*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,c) in compressible flows. */
	double **Gradient_Secondary;	/*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */
  double *Limiter_Secondary;    /*!< \brief Limiter of the primitive variables (T,vx,vy,vz,P,rho). */

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
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.		 
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_nDim, 
			unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CEulerVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */		
	virtual ~CEulerVariable(void);

	/*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	void SetGradient_PrimitiveZero(unsigned short val_primvar);

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
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the primitive variables gradient.
	 */
	double GetLimiter_Primitive(unsigned short val_var);

	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);

  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetLimiter_Primitive(unsigned short val_var, double val_value);
  
	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double **GetGradient_Primitive(void);
  
  /*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double *GetLimiter_Primitive(void);

  /*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	void SetGradient_SecondaryZero(unsigned short val_secondaryvar);
  
	/*!
	 * \brief Add <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the primitive variables.
	 */
	void AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
	 */
	void SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	double GetGradient_Secondary(unsigned short val_var, unsigned short val_dim);
  
  /*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	double GetLimiter_Secondary(unsigned short val_var);
  
	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, double val_value);
  
  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetLimiter_Secondary(unsigned short val_var, double val_value);
  
	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double **GetGradient_Secondary(void);
  
  /*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double *GetLimiter_Secondary(void);

    /*!
 	 * \brief A virtual member.
 	 */
    void SetdPdrho_e(double dPdrho_e);

     /*!
   	 * \brief A virtual member.
   	 */
    void SetdPde_rho(double dPde_rho);
  
	/*!
	 * \brief Set the value of the pressure.
	 */
	bool SetPressure(double Gamma);

	/*!
	 * \brief Set the value of the speed of the sound.
	 * \param[in] Gamma - Value of Gamma.
	 */
	bool SetSoundSpeed(double Gamma);

	/*!
	 * \brief Set the value of the enthalpy.
	 */
	void SetEnthalpy(void);
	
//	/*!
//	 * \brief Set all the primitive variables for compressible flows.
//	 */
//	bool SetPrimVar_Compressible(CConfig *config);

	/*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_Compressible(CFluidModel *FluidModel);

	/*!
	 * \brief A virtual member.
	 */
	void SetSecondaryVar_Compressible(CFluidModel *FluidModel);

	/*!
	 * \brief Set all the primitive variables for incompressible flows.
	 */
	bool SetPrimVar_Incompressible(double Density_Inf, CConfig *config);
  
  /*!
	 * \brief Set all the primitive variables for incompressible flows.
	 */
	bool SetPrimVar_FreeSurface(CConfig *config);
	
	/*!
	 * \brief Get the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the primitive variable for the index <i>val_var</i>.
	 */
	double GetPrimitive(unsigned short val_var);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetPrimitive(unsigned short val_var, double val_prim);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_prim - Primitive variables.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetPrimitive(double *val_prim);

	/*!
	 * \brief Get the primitive variables of the problem.
	 * \return Pointer to the primitive variable vector.
	 */
	double *GetPrimitive(void);
  
  /*!
	 * \brief Get the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the primitive variable for the index <i>val_var</i>.
	 */
	double GetSecondary(unsigned short val_var);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetSecondary(unsigned short val_var, double val_secondary);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_prim - Primitive variables.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetSecondary(double *val_secondary);
  
	/*!
	 * \brief Get the primitive variables of the problem.
	 * \return Pointer to the primitive variable vector.
	 */
	double *GetSecondary(void);
  
	/*!
	 * \brief Set the value of the density for the incompressible flows.
	 */
	void SetDensityInc(double val_density);
  
  /*!
	 * \brief Set the value of the density for the incompressible flows.
	 */
	bool SetDensity(void);
  
  /*!
	 * \brief Set the value of the density for the incompressible flows.
	 */
	void SetPressureInc(void);
  
  /*!
	 * \brief Set the value of the density for the incompressible flows.
	 */
	void SetVelocityInc(void);

	/*!
	 * \brief Set the value of the beta coeffient for incompressible flows.
	 */
	void SetBetaInc2(double val_betainc2);

	/*!
	 * \brief Set the value of the temperature.
	 * \param[in] Gas_Constant - Value of Gas Constant
	 */
	bool SetTemperature(double Gas_Constant);

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
	 * \brief Get the flow pressure.
	 * \return Value of the flow pressure.
	 */
	double GetPressureInc(void);
  
	/*!
	 * \brief Get the speed of the sound.
	 * \return Value of speed of the sound.
	 */
	double GetSoundSpeed(void);

	/*!
	 * \brief Get the value of density for the incompressible flow
	 * \return Value of beta squared.
	 */
	double GetDensityInc(void);

  /*!
	 * \brief Get the value of levelset for the freesurface flows
	 * \return Value of beta squared.
	 */
	double GetLevelSet(void);
  
  /*!
	 * \brief Get the value of distance for the freesurface flows
	 * \return Value of beta squared.
	 */
	double GetDistance(void);
  
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
	 * \brief Set the velocity vector from the solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */	
	void SetVelocity(void);

	/*!
	 * \brief Set the velocity vector from the old solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */		
	void SetVelocity_Old(double *val_velocity);
  
  /*!
	 * \brief Set the velocity vector from the old solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */
	void SetVelocityInc_Old(double *val_velocity);

	/*!
	 * \brief Set the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the time spectral source term. for the index <i>val_var</i>.
	 */
	void SetTimeSpectral_Source(unsigned short val_var, double val_source);

	/*!
	 * \brief Get the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the time spectral source term for the index <i>val_var</i>.
	 */
	double GetTimeSpectral_Source(unsigned short val_var);

	/*!
	 * \brief Get the value of the preconditioner Beta.
	 * \return Value of the low Mach preconditioner variable Beta
	 */
	double GetPreconditioner_Beta();

	/*!
	 * \brief Set the value of the preconditioner Beta.
	 * \param[in] Value of the low Mach preconditioner variable Beta
	 */
	void SetPreconditioner_Beta(double val_Beta);

	/*!
	 * \brief Set the value of the magnetic field
	 * \param[in] Value of the magnetic field
	 */
	void SetMagneticField(double* val_B);
    
    /*!
	 * \brief Get the value of the wind gust
	 * \return Value of the wind gust
	 */
	double* GetWindGust();
    
	/*!
	 * \brief Set the value of the wind gust
	 * \param[in] Value of the wind gust
	 */
	void SetWindGust(double* val_WindGust);
    
    /*!
	 * \brief Get the value of the derivatives of the wind gust
	 * \return Value of the derivatives of the wind gust
	 */
	double* GetWindGustDer();
    
	/*!
	 * \brief Set the value of the derivatives of the wind gust
	 * \param[in] Value of the derivatives of the wind gust
	 */
	void SetWindGustDer(double* val_WindGust);
};

/*! 
 * \class CNSVariable
 * \brief Main class for defining the variables of the Navier-Stokes' solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CNSVariable : public CEulerVariable {
private:
	double Prandtl_Lam;     /*!< \brief Laminar Prandtl number. */
	double Prandtl_Turb;    /*!< \brief Turbulent Prandtl number. */
	double Temperature_Ref; /*!< \brief Reference temperature of the fluid. */
	double Viscosity_Ref;   /*!< \brief Reference viscosity of the fluid. */
	double Viscosity_Inf;   /*!< \brief Viscosity of the fluid at the infinity. */
	double Vorticity[3];    /*!< \brief Vorticity of the fluid. */
	double StrainMag;       /*!< \brief Magnitude of rate of strain tensor. */
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
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.		 
	 * \param[in] config - Definition of the particular problem.
	 */
	CNSVariable(double val_density, double *val_velocity, 
			double val_energy, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	
	 */
	CNSVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CNSVariable(void);

	/*!
	 * \brief Set the laminar viscosity.
	 */
	void SetLaminarViscosity(double laminarViscosity);

	/*!
	 * \overload
	 * \param[in] val_laminar_viscosity_inc - Value of the laminar viscosity (incompressible flows).
	 */
	void SetLaminarViscosityInc(double val_laminar_viscosity_inc);

	/*!
	 * \brief Set the laminar viscosity.
	 */
	void SetThermalConductivity(double thermalConductivity);

	/*!
	 * \brief Set the specific heat Cp.
	 */
	void SetSpecificHeatCp(double Cp);

	/*!
	 * \brief Set the vorticity value.
	 */
	bool SetVorticity(bool val_limiter);

	/*!
	 * \brief Set the rate of strain magnitude.
	 */
	bool SetStrainMag(bool val_limiter);

	/*!
	 * \overload
	 * \param[in] eddy_visc - Value of the eddy viscosity.
	 */
	void SetEddyViscosity(double eddy_visc);
  
  /*!
	 * \overload
	 * \param[in] eddy_visc - Value of the eddy viscosity.
	 */
	void SetEddyViscosityInc(double eddy_visc);

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
	 * \brief Get the thermal conductivity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetThermalConductivity(void);

	/*!
	 * \brief Get the eddy viscosity of the flow.
	 * \return The eddy viscosity of the flow.
	 */
	double GetEddyViscosity(void);

	/*!
	 * \brief Get the specific heat at constant P of the flow.
	 * \return Value of the specific heat at constant P  of the flow.
	 */
	double GetSpecificHeatCp(void);

    /*!
	 * \brief Get the eddy viscosity of the flow.
	 * \return The eddy viscosity of the flow.
	 */
	double GetEddyViscosityInc(void);

	/*!
	 * \brief Set the temperature at the wall
	 */
	void SetWallTemperature(double temperature_wall);

	/*!
	 * \brief Get the value of the vorticity.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */	
	double *GetVorticity(void);

	/*!
	 * \brief Get the value of the magnitude of rate of strain.
	 * \return Value of the rate of strain magnitude.
	 */
	double GetStrainMag(void);
  
  /*!
   * \brief Set the derivative of temperature with respect to density (at constant internal energy).
   */
  void SetdTdrho_e(double dTdrho_e);
  
  /*!
   * \brief Set the derivative of temperature with respect to internal energy (at constant density).
   */
  void SetdTde_rho(double dTde_rho);
  
  /*!
   * \brief Set the derivative of laminar viscosity with respect to density (at constant temperature).
   */
  void Setdmudrho_T(double dmudrho_T);
  
  /*!
   * \brief Set the derivative of laminar viscosity with respect to temperature (at constant density).
   */
  void SetdmudT_rho(double dmudT_rho);
  
  /*!
   * \brief Set the derivative of thermal conductivity with respect to density (at constant temperature).
   */
  void Setdktdrho_T(double dktdrho_T);
  
  /*!
   * \brief Set the derivative of thermal conductivity with respect to temperature (at constant density).
   */
  void SetdktdT_rho(double dktdT_rho);
  
  /*!
   * \brief Set all the primitive variables for compressible flows
   */
  bool SetPrimVar_Compressible(double eddy_visc, double turb_ke, CFluidModel *FluidModel);

	/*!
	 * \brief Set all the secondary variables (partial derivatives) for compressible flows
	 */
	void SetSecondaryVar_Compressible(CFluidModel *FluidModel);

	/*!
	 * \brief Set all the primitive variables for incompressible flows
	 */
	bool SetPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf, double eddy_visc, double turb_ke, CConfig *config);
  
  /*!
	 * \brief Set all the primitive variables for incompressible flows
	 */
	bool SetPrimVar_FreeSurface(double eddy_visc, double turb_ke, CConfig *config);
};

/*! 
 * \class CTurbVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 3.2.8.1 "eagle"
 */
class CTurbVariable : public CVariable {
protected:
	double muT;                /*!< \brief Eddy viscosity. */
	double *TS_Source; 	       /*!< \brief Time spectral source term. */

public:
	/*!
	 * \brief Constructor of the class. 
	 */	
	CTurbVariable(void);

	/*!
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

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
 * \version 3.2.8.1 "eagle"
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
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */
	~CTurbSAVariable(void);

	/*!
	 * \brief Set the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_source - Value of the time spectral source term. for the index <i>val_var</i>.
	 */
	void SetTimeSpectral_Source(unsigned short val_var, double val_source);

	/*!
	 * \brief Get the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the time spectral source term for the index <i>val_var</i>.
	 */
	double GetTimeSpectral_Source(unsigned short val_var);

};


/*!
 * \class CTurbMLVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 3.2.8.1 "eagle"
 */

class CTurbMLVariable : public CTurbVariable {
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbMLVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
	 * \param[in] val_muT  - The eddy viscosity
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbMLVariable(double val_nu_tilde, double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbMLVariable(void);
  
	/*!
	 * \brief Set the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_source - Value of the time spectral source term. for the index <i>val_var</i>.
	 */
	void SetTimeSpectral_Source(unsigned short val_var, double val_source);
  
	/*!
	 * \brief Get the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the time spectral source term for the index <i>val_var</i>.
	 */
	double GetTimeSpectral_Source(unsigned short val_var);
  
};

/*!
 * \class CTransLMVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 3.2.8.1 "eagle"
 */

class CTransLMVariable : public CTurbVariable {
protected:
  double gamma_sep;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CTransLMVariable(void);

	/*!
	 * \overload
	 * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
	 * \param[in] val_REth
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CTransLMVariable(double val_nu_tilde, double val_intermittency, double val_REth, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */
	~CTransLMVariable(void);

  /*!
	 * \brief ________________.
	 */
  double GetIntermittency(void);
  
  /*!
	 * \brief ________________.
	 * \param[in] gamma_sep_in
	 */
  void SetGammaSep(double gamma_sep_in);
  
  /*!
	 * \brief ________________.
	 */
  void SetGammaEff(void);
  
};

/*! 
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 3.2.8.1 "eagle"
 */

class CTurbSSTVariable : public CTurbVariable {
protected:
	double sigma_om2,
	beta_star;
	double F1,		/*!< \brief Menter blending function for blending of k-w and k-eps. */
	F2,		        /*!< \brief Menter blending function for stress limiter. */
	CDkw;           /*!< \brief Cross-diffusion. */
  
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSSTVariable(void);

	/*!
	 * \overload
	 * \param[in] val_rho_kine - Turbulent variable value (initialization value).
	 * \param[in] val_rho_omega - Turbulent variable value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSSTVariable(double val_rho_kine, double val_rho_omega, double val_muT, unsigned short val_nDim, unsigned short val_nvar,
			double *constants, CConfig *config);

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
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
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
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAdjPotentialVariable(double val_psi, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */

	~CAdjPotentialVariable(void);
};

/*! 
 * \class CAdjEulerVariable
 * \brief Main class for defining the variables of the adjoint Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CAdjEulerVariable : public CVariable {
protected:
	double *Psi;		/*!< \brief Vector of the adjoint variables. */
	double *ForceProj_Vector;	/*!< \brief Vector d. */
	double *ObjFuncSource;    /*!< \brief Vector containing objective function sensitivity for discrete adjoint. */
	double *IntBoundary_Jump;	/*!< \brief Interior boundary jump vector. */
	double *TS_Source;		/*!< \brief Time spectral source term. */
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
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjEulerVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjEulerVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	virtual ~CAdjEulerVariable(void);

  /*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_Compressible(double SharpEdge_Distance, bool check, CConfig *config);
  
  /*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_Incompressible(double SharpEdge_Distance, bool check, CConfig *config);
  
  /*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_FreeSurface(double SharpEdge_Distance, bool check, CConfig *config);
  
	/*!
	 * \brief Set the value of the adjoint velocity.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */	
	void SetPhi_Old(double *val_phi);

	/*!
	 * \brief Set the value of the force projection vector.
	 * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
	 */		
	void SetForceProj_Vector(double *val_ForceProj_Vector);

	/*!
	 * \brief Set the value of the objective function source.
	 * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
	 */
	void SetObjFuncSource(double *val_SetObjFuncSource);

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
	 * \brief Get the value of the objective function source.
	 * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
	 */
	double *GetObjFuncSource(void);

	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */		
	double *GetIntBoundary_Jump(void);

	/*!
	 * \brief Set the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the time spectral source term. for the index <i>val_var</i>.
	 */
	void SetTimeSpectral_Source(unsigned short val_var, double val_source);

	/*!
	 * \brief Get the time spectral source term.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the time spectral source term for the index <i>val_var</i>.
	 */
	double GetTimeSpectral_Source(unsigned short val_var);
};

/*! 
 * \class CAdjNSVariable
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CAdjNSVariable : public CAdjEulerVariable {	
private:
  
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
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAdjNSVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	CAdjNSVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

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
 * \version 3.2.8.1 "eagle"
 */
class CAdjTurbVariable : public CVariable {
protected:
	double *dmuT_dUTvar;       /*!< \brief Sensitivity of eddy viscosity to mean flow and turbulence vars. */
	double **dRTstar_dUTvar; 	/*!< \brief Sensitivity of modified turbulence residual (no boundary flux)
	 	 	 	 	 	 	 	 to mean flow and turbulence vars. */
	double **dFT_dUTvar; 	/*!< \brief Sensitivity of boundary flux
		 	 	 	 	 	 	 	 to mean flow and turbulence vars. */
	double *EddyViscSens;    /*!< \brief Eddy Viscosity Sensitivity. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */		
	CAdjTurbVariable(void);

	/*!
	 * \overload
	 * \param[in] val_psinu_inf - Value of the adjoint turbulence variable at the infinity (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CAdjTurbVariable(double val_psinu_inf, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */		
	~CAdjTurbVariable(void);

	/*!
	 * \brief Set the Eddy Viscosity Sensitivity of the problem.
	 * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
	 */
	void SetEddyViscSens(double *val_EddyViscSens, unsigned short numTotalVar);

	/*!
	 * \brief Get the Eddy Viscosity Sensitivity of the problem.
	 * \return Pointer to the Eddy Viscosity Sensitivity.
	 */
	double *GetEddyViscSens(void);
};

/*! 
 * \class CLinPotentialVariable
 * \brief Main class for defining the variables of the linearized potential equation.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CLinPotentialVariable : public CVariable {
public:	
};

/*! 
 * \class CLinEulerVariable
 * \brief Main class for defining the variables of the linearized Euler's equations.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
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
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.	
	 * \param[in] config - Definition of the particular problem.
	 */		
	CLinEulerVariable(double val_deltarho, double *val_deltavel, double val_deltae, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the linearized value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	CLinEulerVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

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
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CLinNSVariable : public CLinEulerVariable {
public:
};

/*! 
 * \class CAdjLevelSetVariable
 * \brief Main class for defining the variables of the Level Set.
 * \ingroup LevelSet_Model
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 */
class CAdjLevelSetVariable : public CVariable {
public:
	/*!
	 * \brief Constructor of the class. 
	 */	
	CAdjLevelSetVariable(void);

	/*!
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjLevelSetVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \overload
	 * \param[in] val_levelset - Level set variable value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CAdjLevelSetVariable(double val_levelset, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjLevelSetVariable(void);

};

/*!
 * \class CTNE2EulerVariable
 * \brief Main class for defining the variables of the TNE2 Euler's solver.
 * \ingroup Euler_Equations
 * \author S. R. Copeland, F. Palacios
 * \version 2.0.6
 */
class CTNE2EulerVariable : public CVariable {
protected:
  bool ionization;       /*!< \brief Presence of charged species in gas mixture. */
  unsigned short nSpecies;  /*!< \brief Number of species in the gas mixture. */
	double Velocity2;			/*!< \brief Square of the velocity vector. */
	double Precond_Beta;	/*!< \brief Low Mach number preconditioner value, Beta. */
  
	/*--- Primitive variable definition ---*/
	double *Primitive;	/*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,c) in compressible flows. */
	double **Gradient_Primitive;	/*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */
  double *Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T,vx,vy,vz,P,rho). */
  double *dPdU;                 /*!< \brief Partial derivative of pressure w.r.t. conserved variables. */
  double *dTdU;  /*!< \brief Partial derivative of temperature w.r.t. conserved variables. */
  double *dTvedU; /*!< \brief Partial derivative of vib.-el. temperature w.r.t. conserved variables. */
  
  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;

  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CTNE2EulerVariable(void);

  /*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of conserved variables.
	 */
  CTNE2EulerVariable(unsigned short val_nDim, unsigned short val_nVar,
                     unsigned short val_nPrimVar,
                     unsigned short val_nPrimVarGrad,
                     CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_pressure
   * \param[in] val_massfrac
	 * \param[in] val_mach
	 * \param[in] val_temperature
	 * \param[in] val_temperature_ve
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nvarprim - Number of primitive variables.
   * \param[in] val_nvarprimgrad - Number of primitive gradient variables.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTNE2EulerVariable(double val_pressure, double *val_massfrac,
                     double *val_mach, double val_temperature,
                     double val_temperature_ve, unsigned short val_nDim,
                     unsigned short val_nvar, unsigned short val_nvarprim,
                     unsigned short val_nvarprimgrad, CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTNE2EulerVariable(double *val_solution, unsigned short val_nDim,
                     unsigned short val_nvar, unsigned short val_nvarprim,
                     unsigned short val_nvarprimgrad, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CTNE2EulerVariable(void);
  
	/*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	void SetGradient_PrimitiveZero(unsigned short val_primvar);
  
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
	 * \brief Set the value of the velocity*velocity.
	 */
	void SetVelocity2(void);
  
  /*!
	 * \brief Set the value of the mixture density.
	 */
	bool SetDensity(void);
  
	/*!
	 * \brief Set the value of the pressure.  Requires T&Tve calculation.
   * \param[in] config
	 */
	bool SetPressure(CConfig *config);
  
	/*!
	 * \brief Set the value of the speed of the sound.
	 * \param[in] config
	 */
	bool SetSoundSpeed(CConfig *config);
  
	/*!
	 * \brief Set the value of the enthalpy.
	 */
	void SetEnthalpy(void);
  
  /*!
   * \brief Sets gas mixture quantities (\f$\rho C^{trans-rot}_v\f$ & \f$\rho C^{vib-el}_v\f$)
   * \param[in] config
   */
  void SetGasProperties(CConfig *config);

  /*!
   * \brief Calculates vib.-el. energy per mass, \f$e^{vib-el}_s\f$, for input species (not including KE)
   */
  double CalcEve(double *V, CConfig *config, unsigned short val_Species);
  
  /*!
   * \brief Calculates enthalpy per mass, \f$h^{vib-el}_s\f$, for input species (not including KE)
   */
  double CalcHs(double *V, CConfig *config, unsigned short val_Species);
  
  /*!
   * \brief Calculates enthalpy per mass, \f$C^{vib-el}_{v_s}\f$, for input species (not including KE)
   */
  double CalcCvve(double val_Tve, CConfig *config, unsigned short val_Species);

  /*!
   * \brief Calculates partial derivative of pressure w.r.t. conserved variables \f$\frac{\partial P}{\partial U}\f$
   * \param[in] config - Configuration settings
   * \param[in] dPdU - Passed-by-reference array to assign the derivatives
   */
  void CalcdPdU(double *V, CConfig *config, double *dPdU);

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void CalcdTdU(double *V, CConfig *config, double *dTdU);
  
  /*!
   * \brief Set partial derivative of vib.-el. temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void CalcdTvedU(double *V, CConfig *config, double *dTvedU);
  
  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void SetdTdrhos(CConfig *config);
  
  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void SetdTvedrhos(CConfig *config);
  
  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  double *GetdPdU(void);
  
  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial T}{\partial \rho_s}\f$
   */
  double *GetdTdU(void);
  
  /*!
   * \brief Set partial derivative of vib.-el. temperature w.r.t. density \f$\frac{\partial T^{V-E}}{\partial \rho_s}\f$
   */
  double *GetdTvedU(void);
  
	/*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_Compressible(CConfig *config);
  
  /*!
	 * \brief Set all the conserved variables.
	 */
	bool Cons2PrimVar(CConfig *config, double *U, double *V, double *dPdU,
                    double *dTdU, double *dTvedU);
  
  /*!
	 * \brief Set all the conserved variables.
	 */
	void Prim2ConsVar(CConfig *config, double *V, double *U);
	
	/*!
	 * \brief Get the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the primitive variable for the index <i>val_var</i>.
	 */
	double GetPrimitive(unsigned short val_var);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetPrimitive(unsigned short val_var, double val_prim);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_prim - Primitive variables.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetPrimitive(double *val_prim);
  
	/*!
	 * \brief Get the primitive variables of the problem.
	 * \return Pointer to the primitive variable vector.
	 */
	double *GetPrimitive(void);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	bool SetTemperature(CConfig *config);
  
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
	 * \brief Get the mass fraction \f$\rho_s / \rho \f$ of species s.
   * \param[in] val_Species - Index of species s.
	 * \return Value of the mass fraction of species s.
	 */
	double GetMassFraction(unsigned short val_Species);
  
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
	 * \brief Sets the temperature of the flow.
	 * \return Value of the temperature of the flow.
	 */
	bool SetTemperature(double val_T);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the vibrational-electronic temperature.
	 */
	double GetTemperature_ve(void);

  /*!
	 * \brief Sets the vibrational electronic temperature of the flow.
	 * \return Value of the temperature of the flow.
	 */
  bool SetTemperature_ve(double val_Tve);
  
  /*!
   * \brief Get the mixture specific heat at constant volume (trans.-rot.).
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  double GetRhoCv_tr(void);
  
  /*!
   * \brief Get the mixture specific heat at constant volume (vib.-el.).
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  double GetRhoCv_ve(void);
  
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
	 * \brief Set the velocity vector from the old solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */
	void SetVelocity_Old(double *val_velocity);
  
  /*!
	 * \brief Get the value of the limiter.
	 */
  double *GetLimiter_Primitive(void);
  
  /*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the primitive variables gradient.
	 */
	double GetLimiter_Primitive(unsigned short val_var);
  
  /*!
	 * \brief Set the value of the limiter.
	 */
  void SetLimiter_Primitive(unsigned short val_var, double val_value);
  
  /*!
	 * \brief Set the value of the limiter.
	 */
  void SetLimiter(unsigned short val_var, double val_value);
  
	/*!
	 * \brief Get the value of the preconditioner Beta.
	 * \return Value of the low Mach preconditioner variable Beta
	 */
	double GetPreconditioner_Beta();
  
	/*!
	 * \brief Set the value of the preconditioner Beta.
	 * \param[in] val_Beta Value of the low Mach preconditioner variable Beta
	 */
	void SetPreconditioner_Beta(double val_Beta);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetRhosIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetRhoIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetPIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetTIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetTveIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho*u
	 */
  unsigned short GetVelIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetHIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetAIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetRhoCvtrIndex(void);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * iRho_s
	 */
  unsigned short GetRhoCvveIndex(void);
  
};

/*!
 * \class CTNE2NSVariable
 * \brief Main class for defining the variables of the TNE2 Navier-Stokes' solver.
 * \ingroup Navier_Stokes_Equations
 * \author S. R. Copeland, F. Palacios
 * \version 2.0.6
 */
class CTNE2NSVariable : public CTNE2EulerVariable {
private:
	double Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
	double Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
	double Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
	double Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
  double *DiffusionCoeff;    /*!< \brief Diffusion coefficient of the mixture. */
	double LaminarViscosity;	/*!< \brief Viscosity of the fluid. */
  double ThermalCond;       /*!< \brief T-R thermal conductivity of the gas mixture. */
  double ThermalCond_ve;    /*!< \brief V-E thermal conductivity of the gas mixture. */
	double Vorticity[3];		/*!< \brief Vorticity of the fluid. */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CTNE2NSVariable(void);
  
  
  /*!
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nprimvar - Number of primitive variables.
   * \param[in] val_nprimvargrad - Number of primitive gradient variables.
	 * \param[in] config - Definition of the particular problem.
	 */
  CTNE2NSVariable(unsigned short val_nDim, unsigned short val_nvar,
                  unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                  CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
   * \param[in] val_massfrac -
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature -
	 * \param[in] val_temperature_ve -
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nvarprim - Number of primitive variables.
   * \param[in] val_nvarprimgrad - Number of primitive gradient variables.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTNE2NSVariable(double val_density, double *val_massfrac, double *val_velocity,
                  double val_temperature, double val_temperature_ve, unsigned short val_nDim,
                  unsigned short val_nvar, unsigned short val_nvarprim,
                  unsigned short val_nvarprimgrad, CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nvarprim - Number of primitive variables.
   * \param[in] val_nvarprimgrad - Number of primitive gradient variables.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTNE2NSVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar,
                  unsigned short val_nvarprim, unsigned short val_nvarprimgrad,
                  CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CTNE2NSVariable(void);

  /*!
	 * \brief Set the laminar viscosity.
	 */
	void SetDiffusionCoeff(CConfig *config);
  
	/*!
	 * \brief Set the laminar viscosity.
	 */
	void SetLaminarViscosity(CConfig *config);
  
  /*!
	 * \brief Get the laminar viscosity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	void SetThermalConductivity(CConfig *config);
  
	/*!
	 * \brief Set the vorticity value.
	 */
	bool SetVorticity(bool val_limiter);
  
  /*!
	 * \brief Get the species diffusion coefficient.
	 * \return Value of the species diffusion coefficient.
	 */
  double* GetDiffusionCoeff(void);
  
	/*!
	 * \brief Get the laminar viscosity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetLaminarViscosity(void);
  
  /*!
	 * \brief Get the thermal conductivity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetThermalConductivity(void);

  /*!
	 * \brief Get the vib-el. thermal conductivity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetThermalConductivity_ve(void);
  
	/*!
	 * \brief Set the temperature at the wall
	 */
	void SetWallTemperature(double temperature_wall);
  
	/*!
	 * \brief Get the value of the vorticity.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */
	double *GetVorticity(void);
	
	/*!
	 * \brief Set all the primitive variables for compressible flows
	 */
	bool SetPrimVar_Compressible(CConfig *config);
  
};


/*!
 * \class CAdjTNE2EulerVariable
 * \brief Main class for defining the variables of the adjoint Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 2.0.6
 */
class CAdjTNE2EulerVariable : public CVariable {
protected:
  unsigned short nSpecies;
	double *Psi;		/*!< \brief Vector of the adjoint variables. */
	double *ForceProj_Vector;	/*!< \brief Vector d. */
	double *ObjFuncSource;    /*!< \brief Vector containing objective function sensitivity for discrete adjoint. */
	double *IntBoundary_Jump;	/*!< \brief Interior boundary jump vector. */
	double *TS_Source;		/*!< \brief Time spectral source term. */
	double Theta;		/*!< \brief Theta variable. */
	bool incompressible;
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjTNE2EulerVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_psirho - Value of the adjoint density (initialization value).
	 * \param[in] val_phi - Value of the adjoint velocity (initialization value).
	 * \param[in] val_psie - Value of the adjoint energy (initialization value).
   * \param[in] val_psieve - Value of the adjoint vibrational energy (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjTNE2EulerVariable(double *val_psirho, double *val_phi,
                        double val_psie, double val_psieve,
                        unsigned short val_nDim, unsigned short val_nvar,
                        CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjTNE2EulerVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjTNE2EulerVariable(void);
  
  /*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_Compressible(double SharpEdge_Distance,
                               bool check,
                               CConfig *config);
  
	/*!
	 * \brief Set the value of the adjoint velocity.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */
	void SetPhi_Old(double *val_phi);
  
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
	 * \brief Set the value of the objective function source.
	 * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
	 */
	void SetObjFuncSource(double *val_SetObjFuncSource);
  
	/*!
	 * \brief Get the value of the force projection vector.
	 * \return Pointer to the force projection vector.
	 */
	double *GetForceProj_Vector(void);
  
	/*!
	 * \brief Get the value of the objective function source.
	 * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
	 */
	double *GetObjFuncSource(void);
  
};

/*!
 * \class CAdjNSVariable
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author S. R. Copeland, F. Palacios
 * \version 2.0.6
 */
class CAdjTNE2NSVariable : public CAdjTNE2EulerVariable {
private:
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjTNE2NSVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_psirho - Value of the adjoint density (initialization value).
	 * \param[in] val_phi - Value of the adjoint velocity (initialization value).
	 * \param[in] val_psie - Value of the adjoint energy (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjTNE2NSVariable(double *val_psirho, double *val_phi,
                     double val_psie, double val_psieve,
                     unsigned short val_nDim, unsigned short val_nvar,
                     CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the adjoint value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjTNE2NSVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAdjTNE2NSVariable(void);
  
	/*!
	 * \brief Set the value of the adjoint velocity.
	 * \param[in] val_phi - Value of the adjoint velocity.
	 */
	void SetPhi_Old(double *val_phi);
  
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
 * \class CTemplateVariable
 * \brief Main class for defining the variables of the potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
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
	 * \param[in] val_nDim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */	
	CTemplateVariable(double val_potential, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */	
	~CTemplateVariable(void);
};

#include "variable_structure.inl"
