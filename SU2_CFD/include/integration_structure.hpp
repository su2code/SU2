/*!
 * \file integration_structure.hpp
 * \brief Headers of the main subroutines for space and time integration. 
 *        The subroutines and functions are in the <i>integration_structure.cpp</i>, 
 *        <i>integration_time.cpp</i>, and <i>integration_notime.cpp</i> files.
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

#ifndef NO_MPI
#include <mpi.h>
#endif
#include <cmath>
#include <iostream>

#include "solver_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \class CIntegration
 * \brief Main class for doing the space integration, time integration, and monitoring 
 *        of a system of Partial Differential Equations (PDE).
 * \author F. Palacios.
 * \version 2.0.6
 */
class CIntegration {
protected:
	double Cauchy_Value,	/*!< \brief Summed value of the convergence indicator. */
	Cauchy_Func;			/*!< \brief Current value of the convergence indicator at one iteration. */
	unsigned short Cauchy_Counter;	/*!< \brief Number of elements of the Cauchy serial. */
	double *Cauchy_Serie;			/*!< \brief Complete Cauchy serial. */
	double Old_Func,	/*!< \brief Old value of the objective function (the function which is monitored). */
	New_Func;			/*!< \brief Current value of the objective function (the function which is monitored). */
	bool Convergence,		/*!< \brief To indicate if the flow solver (direct, adjoint, or linearized) has converged or not. */
	Convergence_OneShot,	/*!< \brief To indicate if the one-shot method has converged. */
	Convergence_FullMG;		/*!< \brief To indicate if the Full Multigrid has converged and it is necessary to add a new level. */
	double InitResidual;	/*!< \brief Initial value of the residual to evaluate the convergence level. */

public:
	
	/*! 
	 * \brief Constructor of the class. 
	 */
	CIntegration(CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CIntegration(void);
	
	/*! 
	 * \brief Do the space integration of the numerical system.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Space_Integration(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config, 
						   unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);

	/*! 
	 * \brief Do the time integration (explicit or implicit) of the numerical system.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	void Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
						  unsigned short iRKStep, unsigned short RunTime_EqSystem, unsigned long Iteration);
	
	/*! 
	 * \brief Initialize the adjoint solution using the primal problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	void Adjoint_Setup(CGeometry ***geometry, CSolver ****solver_container, CConfig **config,
                     unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);
	
	/*!
	 * \brief Numerical method for solving a non-time depending equation, like the potential equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution - Solution of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solving_Linear_System(CGeometry *geometry, CSolver *solver, CSolver **solver_container, 
							   CConfig *config, unsigned short iMesh);

	/*! 
	 * \brief Do the convergence analisys to determine if the code must stop the execution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] Iteration - Current iteration.
	 * \param[in] monitor - Objective function that is use to study its convergence.
	 */
	void Convergence_Monitoring(CGeometry *geometry, CConfig *config, 
								unsigned long Iteration, double monitor);
	
	/*! 
	 * \brief Get the value of the convergence.
	 * \return Level of convergence of the solution.
	 */
	double GetCauchy_Value(void);
	
	/*! 
	 * \brief Get the indicator of the convergence for the direct, adjoint and linearized problem.
	 * \return <code>TRUE</code> means that the convergence criteria is satisfied; 
	 *         otherwise <code>FALSE</code>.
	 */
	bool GetConvergence(void);
	
	/*! 
	 * \brief Set the indicator of the convergence.
	 * \param[in] value - <code>TRUE</code> means that the convergence criteria is satisfied; 
	 *            otherwise <code>FALSE</code>.
	 */
	void SetConvergence(bool value);
	
	/*! 
	 * \brief Get the indicator of the convergence for the one-shot problem.
	 * \return <code>TRUE</code> means that the convergence criteria is satisfied; 
	 *         otherwise <code>FALSE</code>.
	 */
	bool GetConvergence_OneShot(void);
	
	/*! 
	 * \brief Get the indicator of the convergence for the full multigrid problem.
	 * \return <code>TRUE</code> means that the convergence criteria is satisfied; 
	 *         otherwise <code>FALSE</code>.
	 */
	bool GetConvergence_FullMG(void);
	
	/*! 
	 * \brief Save the solution, and volume at different time steps. 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution - Flow solution.
   * \param[in] config - Definition of the particular problem.
	 */
	void SetDualTime_Solver(CGeometry *geometry, CSolver *solver, CConfig *config);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	virtual void MultiGrid_Iteration(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
								  CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] mu - Variable for controlling the kind of multigrid algorithm.	 
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	virtual void MultiGrid_Cycle(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
							   CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
							   unsigned long Iteration, unsigned short iZone);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 * \param[in] monitor - value of the non-dimensional parameters for monitoring the convergence.
	 */
	virtual void NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container, 
																				 CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, unsigned long Iteration, 
																				 double *monitor);
	
	/*! 
	 * \brief A virtual member.
	 * \param[out] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[out] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver **sol_fine, CSolver **sol_coarse, 
										 CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
										CGeometry *geo_coarse, CConfig *config);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] InclSharedDomain - Include the shared domain in the interpolation.
	 */
	virtual void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver **sol_fine, CSolver **sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver **sol_fine, CSolver **sol_coarse, 
										CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] flow - Flow solution.
	 */
	virtual void SetResidual_Term(CGeometry *geometry, CSolver *flow);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, 
								 CConfig *config);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	virtual void SingleGrid_Iteration(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
								  CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);
	
	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void SetPotential_Solver(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
                                   CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] solution - Container vector with all the solutions on the finest grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Smooth_Solution(unsigned short RunTime_EqSystem, CSolver **solver, CGeometry *geometry,
                       unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config);

};

/*!
 * \class CMultiGridIntegration
 * \brief Class for doing the numerical integration using a multigrid method.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CMultiGridIntegration : public CIntegration {
protected:
    
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] config - Definition of the particular problem.
	 */
	CMultiGridIntegration(CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CMultiGridIntegration(void);
	
	/*! 
	 * \brief This subroutine calls the MultiGrid_Cycle and also prepare the multigrid levels and the monitoring.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	void MultiGrid_Iteration(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
							 CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);
	
	/*! 
	 * \brief Perform a Full-Approximation Storage (FAS) Multigrid. 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] mu - Variable for controlling the kind of multigrid algorithm.	 
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	void MultiGrid_Cycle(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
                     CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
                     unsigned long Iteration, unsigned short iZone);
	
	/*! 
	 * \brief Compute the non-dimensional parameters.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 * \param[in] monitor - value of the non-dimensional parameters for monitoring the convergence.
	 */
	void NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container, 
																 CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, unsigned long Iteration, 
																 double *monitor);

	/*! 
	 * \brief Compute the fine solution from a coarse solution. 
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[out] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, 
								 CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
	
	/*! 
	 * \brief Compute the fine grid correction from the coarse solution. 
	 * \param[out] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	void GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
																 CGeometry *geo_coarse, CConfig *config);
	
	/*! 
	 * \brief Do an implicit smoothing of the prolongated correction. 
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] solution - Container vector with all the solutions on the finest grid.	 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.		 
	 * \param[in] config - Definition of the particular problem.
	 */
	void SmoothProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
																		 unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config);
  
  /*!
	 * \brief Do an implicit smoothing of the solution.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] solution - Container vector with all the solutions on the finest grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Smooth_Solution(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                    unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config);

	/*!
	 * \brief Set the value of the corrected fine grid solution.
	 * \param[out] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config);

	/*! 
	 * \brief Compute truncation error in the coarse grid using the fine grid information. 
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, 
								CGeometry *geo_coarse, CConfig *config);

	/*! 
	 * \brief Restrict solution from fine grid to a coarse grid. 
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] InclSharedDomain - Include the shared domain in the interpolation.
	 */
	void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
	
	/*! 
	 * \brief Compute the gradient in coarse grid using the fine grid information. 
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver **sol_fine, CSolver **sol_coarse, 
								CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

	/*! 
	 * \brief Add the truncation error to the residual. 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] flow - Flow solution.
	 */
	void SetResidual_Term(CGeometry *geometry, CSolver *flow);
	
	/*! 
	 * \brief Compute the forcing term. 
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, 
						 CConfig *config);
};

/*! 
 * \class CSingleGridIntegration
 * \brief Class for doing the numerical integration of the turbulence model.
 * \author A. Bueno.
 * \version 2.0.6
 */
class CSingleGridIntegration : public CIntegration {
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSingleGridIntegration(CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSingleGridIntegration(void);
	
	/*! 
	 * \brief Do the numerical integration (implicit) of the turbulence solver. 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] Iteration - Current iteration.
	 */
	void SingleGrid_Iteration(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
							 CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);
  
  /*!
	 * \brief Restrict solution from fine grid to a coarse grid.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] sol_fine - Pointer to the solution on the fine grid.
	 * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
	 * \param[in] geo_fine - Geometrical definition of the fine grid.
	 * \param[in] geo_coarse - Geometrical definition of the coarse grid.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] InclSharedDomain - Include the shared domain in the interpolation.
	 */
	void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver **sol_fine, CSolver **sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);
  
};

/*! 
 * \class CPotentialIntegration
 * \brief Class for doing the numerical integration of the potential equation.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CPotentialIntegration : public CIntegration {
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPotentialIntegration(CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CPotentialIntegration(void);
	
	/*! 
	 * \brief Do the numerical integration (Galerkin explicit) of the potential equation. 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetPotential_Solver(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
                           CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone);
};

#include "integration_structure.inl"
