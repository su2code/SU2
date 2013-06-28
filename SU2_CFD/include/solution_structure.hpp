/*!
 * \file solution_structure.hpp
 * \brief Headers of the main subroutines for solving the partial differential equations. 
 *        The subroutines and functions are in the <i>solution_structure.cpp</i>, 
 *        <i>solution_direct.cpp</i>, <i>solution_adjoint.cpp</i>, and 
 *        <i>solution_linearized.cpp</i> files.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#ifndef NO_MPI
#include <mpi.h>
#endif

#include "numerics_structure.hpp"
#include "variable_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/sparse_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"

using namespace std;

/*! 
 * \class CSolution
 * \brief Main class for defining Computational 
 *		  Fluid Dynamics subroutines.
 * \author F. Palacios.
 * \version 1.0.
 */
class CSolution {
protected:
	unsigned short nVar,	/*!< \brief Number of variables of the problem. */
	nDim;					/*!< \brief Number of dimensions of the problem. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	unsigned short nSpecies, nFluids;/*!< \brief Number of Species present in plasma and number of fluids to model the plasma as. */
	unsigned short nMonatomics, nDiatomics;
	unsigned long nPoint;	/*!< \brief Number of points of the computational grid. */
	double *Residual_Max,	/*!< \brief Vector with the maximal residual for each variable. */
	*Residual,				/*!< \brief Auxiliary vector for storing the residual. */
	*Residual_i,			/*!< \brief Auxiliary vector for storing the residual at point i. */
	*Residual_j;			/*!< \brief Auxiliary vector for storing the residual at point j. */
	double *Solution,	/*!< \brief Auxiliary vector for storing the solution. */ 
	*Solution_i,		/*!< \brief Auxiliary vector for storing the solution at point i. */ 
	*Solution_j;		/*!< \brief Auxiliary vector for storing the solution at point j. */ 
	double *Vector,		/*!< \brief Auxiliary vector to do the reconstruction of the variables. */ 
	*Vector_i,			/*!< \brief Auxiliary vector to do the reconstruction of the variables at point i. */ 
	*Vector_j;			/*!< \brief Auxiliary vector to do the reconstruction of the variables at point j. */ 
	double *Res_Conv_i,		/*!< \brief Auxiliary vector for storing the convective residual at point i (adjoint). */
	*Res_Visc_i,			/*!< \brief Auxiliary vector for storing the viscous residual at point i (adjoint). */
	*Res_Conv_j,			/*!< \brief Auxiliary vector for storing the convective residual at point j (adjoint). */
	*Res_Visc_j,			/*!< \brief Auxiliary vector for storing the viscous residual at point j (adjoint). */
	*Res_Conv,				/*!< \brief Auxiliary vector for storing the convective residual (direct). */
	*Res_Visc;				/*!< \brief Auxiliary vector for storing the viscous residual (direct). */
	CSparseMatrix Jacobian; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
	double **Jacobian_i,	/*!< \brief Auxiliary matrices for storing point to point Jacobians at point i. */
	**Jacobian_j;			/*!< \brief Auxiliary matrices for storing point to point Jacobians at point j. */
	double **Jacobian_ii,	/*!< \brief Auxiliary matrices for storing point to point Jacobians (adjoint problems). */
	**Jacobian_ij,			/*!< \brief Auxiliary matrices for storing point to point Jacobians (adjoint problems). */
	**Jacobian_ji,			/*!< \brief Auxiliary matrices for storing point to point Jacobians (adjoint problems). */
	**Jacobian_jj;			/*!< \brief Auxiliary matrices for storing point to point Jacobians (adjoint problems). */
	double *xsol; /*!< \brief vector to store iterative solution of implicit linear system. */
	double *xres; /*!< \brief vector to store iterative residual of implicit linear system. */
	double *rhs;  /*!< \brief right hand side of implicit linear system. */
	double **Smatrix,	/*!< \brief Auxiliary structure for computing gradients by least-squares */
	**cvector;			/*!< \brief Auxiliary structure for computing gradients by least-squares */
	double **tau;		/*!< \brief Matrix to store viscous stresses for forces computation. */
	CSparseMatrix StiffMatrix; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	double **StiffMatrix_Elem,	/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	**StiffMatrix_Node;			/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	double *Source_Vector;		/*!< \brief Auxiliary vector for storing element source vector. */
	double Max_Delta_Time,	/*!< \brief Maximum value of the delta time for all the control volumes. */
	Min_Delta_Time,					/*!< \brief Minimum value of the delta time for all the control volumes. */
	Sum_Delta_Time;					/*!< \brief Sum of the minimum value of the delta time for all the control volumes. */

public:
	CVariable** node;	/*!< \brief Vector which the define the variables for each problem. */

	/*! 
	 * \brief Constructor of the class. 
	 */
	CSolution(void);

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CSolution(void);

	/*!
	 * \brief Get the value of the maximum delta time.
	 * \return Value of the maximum delta time.
	 */
	double GetMax_Delta_Time(void);

	/*!
	 * \brief Get the value of the minimum delta time.
	 * \return Value of the minimum delta time.
	 */
	double GetMin_Delta_Time(void);

	/*!
	 * \brief Get the value of the summatory delta time.
	 * \return summatory delta time.
	 */
	double GetSum_Delta_Time(void);

	/*! 
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnVar(void);

	/*! 
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnSpecies(void);

	/*! 
	 * \brief Get the number of fluids to model the problem with.
	 */
	unsigned short GetnFluids(void);

	/*! 
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnMonatomics(void);

	/*! 
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnDiatomics(void);

	/*! 
	 * \brief Compute the combination between convective and viscous residual in a RK-Stage.
	 * \attention It is not necessary to store all the viscous residual, just storing the k-1 is enougth
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */	
	void SetResidual_RKCoeff(CGeometry *geometry, CConfig *config, unsigned short iRKStep);

	/*! 
	 * \brief Set the total residual (convective + viscous).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */	
	void SetResidual_Total(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep, unsigned short iMesh);

	/*! 
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */	
	void SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep, unsigned short iMesh);

	/*! 
	 * \brief Set the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void SetRes_Max(unsigned short val_var, double val_residual);

	/*! 
	 * \brief Adds the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void AddRes_Max(unsigned short val_var, double val_residual);

	/*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the bigest residual for the variable in the position <i>val_var</i>.
	 */
	double GetRes_Max(unsigned short val_var);

	/*!
	 * \brief Compute a Jacobi implicit smoothing of the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.		 
	 */
	void SetResidual_Smoothing(CGeometry *geometry, unsigned short val_nSmooth, double val_smooth_coeff);

	/*! 
	 * \brief Compute a Jacobi implicit smoothing of the solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] val_nSmooth - Number of smoothing iterations.
	 * \param[in] val_smooth_coeff - Relaxation factor.		 
	 */
	void SetSolution_Smoothing(CGeometry *geometry, unsigned short val_nSmooth, double val_smooth_coeff);

	/*!
	 * \brief Set Value of the residual if there is a grid movement.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetGrid_Movement_Residual(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the Green-Gauss gradient of the auxiliary variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetAuxVar_Gradient_GG(CGeometry *geometry);

	/*!
	 * \brief Compute the Least Squares gradient of the auxiliary variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetAuxVar_Gradient_LS(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the Least Squares gradient of an auxiliar variable on the profile surface.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetAuxVar_Surface_Gradient(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the Green-Gauss gradient of the solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetSolution_Gradient_GG(CGeometry *geometry);

	/*!
	 * \brief Compute the Least Squares gradient of the solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the Least Squares gradient of the solution on the profile surface.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSurface_Gradient(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the Venkatakrishnan slope limiter.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSolution_Limiter(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the pressure laplacian using in a incompressible solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] PressureLaplacian - Pressure laplacian.
	 */
	void SetPressureLaplacian(CGeometry *geometry, double *PressureLaplacian);

	/*!
	 * \brief Set the old solution variables to the current solution value for Runge-Kutta iteration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void Set_OldSolution(CGeometry *geometry);

	/*!
	 * \brief Creates a copy of the solution in <i>OldSolution</i> to do a residual computation starting 
	 *        with a interpolated solution, not the original one (this is useful in the multigrid cycle).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] index - If greater that 0 means copy from <i>Solution</i> to <i>OldSolution</i>; 
	 *            Otherwise means the reverse.
	 */
	void Set_MultiSolution(CGeometry *geometry, short index);

	/*!
	 * \brief Initializes the structure of the whole problem Jacobian.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void InitializeJacobianStructure(CGeometry *geometry, CConfig *config);
	/*!
	 * \brief Initializes the structure of the whole problem Stiffness Matrix.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void InitializeStiffMatrixStructure(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetSpectral_Radius(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	virtual void SetPress_Switch(CGeometry *geometry);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Dirichlet(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	virtual void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	virtual void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Dielectric(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Electrode(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Inviscid_Forces(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Inviscid_DeltaForces(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Viscous_Forces(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Viscous_DeltaForces(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void SourceConserv_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Charge_Dist_SourceTerm(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \param[in] val_sensitivity - Value of the sensitivity coefficient.
	 */
	virtual void SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetForceProj_Vector(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*! 
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	virtual void SetIntBoundary_Jump(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CDrag - Value of the total drag coefficient.
	 */
	virtual void SetTotal_CDrag(double val_Total_CDrag);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CCharge - Value of the total charge coefficient.
	 */
	virtual void SetTotal_CCharge(double val_Total_CCharge);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CLift - Value of the total lift coefficient.
	 */
	virtual void SetTotal_CLift(double val_Total_CLift);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetDistance(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void FEMEikonalSolver(CGeometry *geometry, CConfig *config);	

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Inviscid_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Viscous_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCLift_Inv(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCLift_Visc(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCDrag_Inv(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCSideForce_Inv(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the pressure coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCPress_Inv(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCEff_Inv(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCDrag_Visc(unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CLift(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the sideforce coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CSideForce(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the pressure coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CPress(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the efficiency coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CEff(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CEquivArea(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the Near-Field Pressure coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CNearFieldPress(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	virtual void SetTotal_CEquivArea(double val_cequivarea);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
	 */
	virtual void SetTotal_CNearFieldPress(double val_cnearfieldpress);

	/*!
	 * \brief A virtual member.
	 * \return Value of the linearized lift coefficient (inviscid contribution).
	 */
	virtual double GetTotal_CDeltaLift(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the charge coefficient.
	 */
	virtual double GetTotal_CCharge(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CDrag(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the moment x coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CMx(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the moment y coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CMy(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the moment y coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CMz(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the linearized drag coefficient (inviscid contribution).
	 */
	virtual double GetTotal_CDeltaDrag(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (inviscid contribution).
	 */
	virtual double GetAllBound_CLift_Inv(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual double GetAllBound_CDrag_Inv(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual double GetAllBound_CSideForce_Inv(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual double GetAllBound_CPress_Inv(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual double GetAllBound_CEff_Inv(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	virtual double GetAllBound_CLift_Visc(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (viscous contribution).
	 */
	virtual double GetAllBound_CDrag_Visc(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual double GetCPressure(unsigned short val_marker, unsigned short val_vertex);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the skin friction coefficient.
	 */
	virtual double GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex);

	/*!
	 * \brief A virtual member.
	 * \return Value of the adjoint density at the infinity.
	 */
	virtual double GetPsiRho_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the adjoint energy at the infinity.
	 */
	virtual double GetPsiE_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the adjoint velocity vector at the infinity.
	 */
	virtual double GetPhi_Inf(unsigned short val_dim);

	/*!
	 * \brief A virtual member.
	 * \return Value of the geometrical sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CSens_Geo(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the Mach sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CSens_Mach(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the angle of attack sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CSens_AoA(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density at the infinity.
	 */
	virtual double GetDensity_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable for the density.
	 * \return Value of the density at the infinity.
	 */
	virtual double GetDensity_Inf(unsigned short val_var);

	/*!
	 * \brief A virtual member.
	 * \return Value of the velocity at the infinity.
	 */
	virtual double GetModVelocity_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density x energy at the infinity.
	 */
	virtual double GetDensity_Energy_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable for the energy.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	virtual double GetDensity_Energy_Inf(unsigned short val_var);

	/*!
	 * \brief A virtual member.
	 * \return Value of the pressure at the infinity.
	 */
	virtual double GetPressure_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the density x velocity at the infinity.
	 */
	virtual double GetDensity_Velocity_Inf(unsigned short val_dim);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \param[in] val_var - Index of the variable for the velocity.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	virtual double GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	virtual double GetVelocity_Inf(unsigned short val_dim);

	/*!
	 * \brief A virtual member.
	 * \return Value of the viscosity at the infinity.
	 */
	virtual double GetViscosity_Inf(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the sensitivity coefficient.
	 */
	virtual double GetCSensitivity(unsigned short val_marker, unsigned short val_vertex);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density at the inlet.
	 */
	virtual double GetDensity_Inlet(void);

	/*!
	 * \overload
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density at the outlet.
	 */
	virtual double GetDensity_Inlet(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density x energy at the inlet.
	 */
	virtual double GetDensity_Energy_Inlet(void);

	/*!
	 * \overload
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density x energy at the outlet.
	 */
	virtual double GetDensity_Energy_Inlet(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the density x velocity at the inlet.
	 */
	virtual double GetDensity_Velocity_Inlet(unsigned short val_dim);

	/*!
	 * \overload
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density x velocity at the outlet.
	 */
	virtual double GetDensity_Velocity_Inlet(unsigned short val_dim,unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density at the outlet.
	 */
	virtual double GetDensity_Outlet(void);

	/*!
	 * \overload
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density at the outlet.
	 */
	virtual double GetDensity_Outlet(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \return Value of the density x energy at the outlet.
	 */
	virtual double GetDensity_Energy_Outlet(void);

	/*!
	 * \overload
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density x energy at the outlet.
	 */
	virtual double GetDensity_Energy_Outlet(unsigned short val_Fluid);

	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the density x velocity at the outlet.
	 */
	virtual double GetDensity_Velocity_Outlet(unsigned short val_dim);

	/*!
	 * \overload
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density x velocity at the outlet.
	 */
	virtual double GetDensity_Velocity_Outlet(unsigned short val_dim,unsigned short val_Fluid);
};

/*! 
 * \class CEulerSolution
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.0.
 */
class CEulerSolution : public CSolution {
protected:
	double Mach_Inf,	/*!< \brief Mach number at the infinity. */
	Mach_Inlet,		/*!< \brief Mach number at the inlet. */
	Mach_Outlet,		/*!< \brief Mach number at the outlet. */
	Density_Inf,	/*!< \brief Density at the infinity. */
	Density_Inlet,		/*!< \brief Density at the inlet. */
	Density_Outlet,		/*!< \brief Density at the outlet. */
	Density_Back,		/*!< \brief Density at infinity behind the Shock. */
	Energy_Inf,			/*!< \brief Energy at the infinity. */
	Energy_Inlet,		/*!< \brief Energy at the inlet. */
	Energy_Outlet,		/*!< \brief Energy at the outlet. */
	Energy_Back,		/*!< \brief Energy at the infinity behind the Shock. */
	Pressure_Inf,		/*!< \brief Pressure at the infinity. */
	Pressure_Inlet,		/*!< \brief Pressure at the inlet. */
	Pressure_Outlet,	/*!< \brief Pressure at the outlet. */
	Pressure_Back,		/*!< \brief Pressure at the infinity behind the Shock. */
	*Velocity_Inf,		/*!< \brief Flow Velocity vector at the infinity. */
	*Velocity_Inlet,	/*!< \brief Flow Velocity vector at the inlet. */
	*Velocity_Outlet,	/*!< \brief Flow Velocity vector at the outlet. */
	*Velocity_Back;		/*!< \brief Flow Velocity vector at the infinity behind the Shock. */
	double *CDrag_Inv,	/*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
	*CLift_Inv,			/*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
	*CSideForce_Inv,		/*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
	*CPress_Inv,		/*!< \brief Pressure coefficient (inviscid contribution) for each boundary. */
	*CMx_Inv,			/*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
	*CMy_Inv,			/*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
	*CMz_Inv,			/*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
	*CEff_Inv,				/*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
	*CEquivArea_Inv,				/*!< \brief Equivalent area (inviscid contribution) for each boundary. */
	*CNearFieldPress_Inv,				/*!< \brief Near field pressure (inviscid contribution) for each boundary. */
	**CPressure,		/*!< \brief Pressure coefficient for each boundary and vertex. */
	*ForceInviscid,		/*!< \brief Inviscid force for each boundary. */
	*MomentInviscid,	/*!< \brief Inviscid moment for each boundary. */
	InverseDesign;	/*!< \brief Inverse design functional for each boundary. */
	double AllBound_CDrag_Inv,	/*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CLift_Inv,			/*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CSideForce_Inv,			/*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CPress_Inv,			/*!< \brief Total press coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMx_Inv,			/*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMy_Inv,			/*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMz_Inv,			/*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEff_Inv,			/*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEquivArea_Inv,			/*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CNearFieldPress_Inv;			/*!< \brief Near-Field press coefficient (inviscid contribution) for all the boundaries. */
	double Total_CDrag, /*!< \brief Total drag coefficient for all the boundaries. */
	Total_CLift,		/*!< \brief Total lift coefficient for all the boundaries. */
	Total_CSideForce,		/*!< \brief Total sideforce coefficient for all the boundaries. */	
	Total_CPress,		/*!< \brief Total press coefficient for all the boundaries. */	
	Total_CMx,			/*!< \brief Total x moment coefficient for all the boundaries. */
	Total_CMy,			/*!< \brief Total y moment coefficient for all the boundaries. */
	Total_CMz,			/*!< \brief Total z moment coefficient for all the boundaries. */
	Total_CEff,			/*!< \brief Total efficiency coefficient for all the boundaries. */
	Total_CEquivArea,			/*!< \brief Total Equivalent Area coefficient for all the boundaries. */
	Total_CNearFieldPress;			/*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
	double *p1_Und_Lapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */ 
	*p2_Und_Lapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */ 
	double *PrimVar_i,	/*!< \brief Auxiliary vector for storing the solution at point i. */
	*PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	unsigned short nPoint, /*!< \brief Number of points of the mesh. */
	nMarker;				/*!< \brief Total number of markers using the grid information. */
	bool space_centered,  /*!< \brief True if space centered scheeme used. */
	euler_implicit,			/*!< \brief True if euler implicit scheme used. */
	least_squares;        /*!< \brief True if computing gradients by least squares. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CEulerSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CEulerSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CEulerSolution(void);

	/*!
	 * \brief Compute the density at the inlet.
	 * \return Value of the density at the infinity.
	 */
	double GetDensity_Inlet(void);

	/*!
	 * \brief Compute the density multiply by energy at the inlet.
	 * \return Value of the density multiply by energy at the inlet.
	 */
	double GetDensity_Energy_Inlet(void);

	/*!
	 * \brief Compute the density multiply by velocity at the inlet.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at the inlet.
	 */
	double GetDensity_Velocity_Inlet(unsigned short val_dim);

	/*!
	 * \brief Compute the density at the inlet.
	 * \return Value of the density at the infinity.
	 */
	double GetDensity_Outlet(void);

	/*!
	 * \brief Compute the density multiply by energy at the inlet.
	 * \return Value of the density multiply by energy at the inlet.
	 */
	double GetDensity_Energy_Outlet(void);

	/*!
	 * \brief Compute the density multiply by velocity at the inlet.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at the inlet.
	 */
	double GetDensity_Velocity_Outlet(unsigned short val_dim);

	/*!
	 * \brief Compute the density at the infinity.
	 * \return Value of the density at the infinity.
	 */
	double GetDensity_Inf(void);

	/*!
	 * \brief Compute the density at infinity  behind the Shock.
	 * \return Value of the density at infinity behind the Shock .
	 */
	double GetDensity_Back(void);

	/*!
	 * \brief Compute 2-norm of the velocity at the infinity.
	 * \return Value of the 2-norm of the velocity at the infinity.
	 */
	double GetModVelocity_Inf(void);

	/*!
	 * \brief Compute the density multiply by energy at the infinity.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	double GetDensity_Energy_Inf(void);

	/*!
	 * \brief Compute the density multiply by energy at infinity behind the Shock.
	 * \return Value of the density multiply by  energy at infinity behind the Shock.
	 */
	double GetDensity_Energy_Back(void);

	/*!
	 * \brief Compute the pressure at the infinity.
	 * \return Value of the pressure at the infinity.
	 */
	double GetPressure_Inf(void);

	/*!
	 * \brief Compute the pressure at infinity behind the Shock.
	 * \return Value of the pressure at infinity behind the Shock.
	 */
	double GetPressure_Back(void);

	/*!
	 * \brief Compute the density multiply by velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	double GetDensity_Velocity_Inf(unsigned short val_dim);

	/*!
	 * \brief Compute the density multiply by velocity at infinity  behind the Shock.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at infinity  behind the Shock.
	 */
	double GetDensity_Velocity_Back(unsigned short val_dim);

	/*!
	 * \brief Get the velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	double GetVelocity_Inf(unsigned short val_dim);

	/*!
	 * \brief Get the velocity at infinity behind the Shock.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the Shock.
	 */
	double GetVelocity_Back(unsigned short val_dim);

	/*!
	 * \brief Compute the time step for solving the Euler equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the spatial integration using a centred scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spectral radius.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSpectral_Radius(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute a pressure sensor switch.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetPress_Switch(CGeometry *geometry);

	/*!
	 * \brief Compute the gradient of the primitive variables using Green-Gauss method, 
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method, 
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the undivided laplacian for the solution, except the energy equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Do the send-receive communications in the MPI parallelization.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Set the boundary contition of the interprocessor boundaries.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the symmetry boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the interface boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the dirichlet boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Dirichlet(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.

	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Compute the pressure forces and all the adimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_Forces(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Provide the non dimensional lift coefficient (inviscid contribution).
	 * \param val_marker Surface where the coefficient is going to be computed.
	 * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCLift_Inv(unsigned short val_marker);

	/*!
	 * \brief Provide the non dimensional drag coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCDrag_Inv(unsigned short val_marker);

	/*!
	 * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCSideForce_Inv(unsigned short val_marker);

	/*!
	 * \brief Provide the non dimensional pressure coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the pressure coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCPress_Inv(unsigned short val_marker);

	/*!
	 * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCEff_Inv(unsigned short val_marker);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
	 * \return Value of the lift coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CLift(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
	 * \return Value of the sideforce coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CSideForce(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional pressure coefficient.
	 * \return Value of the pressure coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CPress(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
	 * \return Value of the efficiency coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CEff(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CEquivArea(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Near-Field pressure coefficient.
	 * \return Value of the NearField pressure coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CNearFieldPress(void);

	/*!
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_CEquivArea(double val_cequivarea);

	/*!
	 * \brief Set the value of the Near-Field pressure oefficient.
	 * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
	 */
	void SetTotal_CNearFieldPress(double val_cnearfieldpress);

	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
	 * \param[in] val_Total_CLift - Value of the total lift coefficient.
	 */
	void SetTotal_CLift(double val_Total_CLift);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CDrag(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
	 * \return Value of the moment x coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CMx(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
	 * \return Value of the moment y coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CMy(void);

	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
	 * \return Value of the moment z coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CMz(void);

	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
	 * \param[in] val_Total_CDrag - Value of the total drag coefficient.
	 */
	void SetTotal_CDrag(double val_Total_CDrag);

	/*!
	 * \brief Get the inviscid contribution to the lift coefficient.
	 * \return Value of the lift coefficient (inviscid contribution).
	 */
	double GetAllBound_CLift_Inv(void);

	/*!
	 * \brief Get the inviscid contribution to the drag coefficient.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	double GetAllBound_CDrag_Inv(void);

	/*!
	 * \brief Get the inviscid contribution to the sideforce coefficient.
	 * \return Value of the sideforce coefficient (inviscid contribution).
	 */
	double GetAllBound_CSideForce_Inv(void);

	/*!
	 * \brief Get the inviscid contribution to the pressure coefficient.
	 * \return Value of the pressure coefficient (inviscid contribution).
	 */
	double GetAllBound_CPress_Inv(void);

	/*!
	 * \brief Get the inviscid contribution to the efficiency coefficient.
	 * \return Value of the efficiency coefficient (inviscid contribution).
	 */
	double GetAllBound_CEff_Inv(void);

	/*!
	 * \brief Provide the Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	double GetCPressure(unsigned short val_marker, unsigned short val_vertex);

};

/*! 
 * \class CNSSolution
 * \brief Main class for defining the Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 1.0.
 */
class CNSSolution : public CEulerSolution {
private:
	double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */
	double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
	Prandtl_Turb;         /*!< \brief Turbulent Prandtl number. */
	double *CDrag_Visc,	/*!< \brief Drag coefficient (viscous contribution) for each boundary. */
	*CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for each boundary. */
	*CMx_Visc,			/*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
	*CMy_Visc,			/*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
	*CMz_Visc,			/*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
	*CEff_Visc,			/*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
	**CSkinFriction;	/*!< \brief Skin friction coefficient for each boundary and vertex. */
	double *ForceViscous,	/*!< \brief Viscous force for each boundary. */
	*MomentViscous;			/*!< \brief Inviscid moment for each boundary. */
	double AllBound_CDrag_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
	AllBound_CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
	AllBound_CMx_Visc,			/*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMy_Visc,			/*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMz_Visc,			/*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEff_Visc;			/*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
	double *PrimVar,	/*!< \brief Auxiliary vector for storing the primitive solution. */ 
	*PrimVar_i,			/*!< \brief Auxiliary vector for storing the primitive solution at point i. */ 
	*PrimVar_j;			/*!< \brief Auxiliary vector for storing the primitive solution at point j. */ 

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CNSSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNSSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CNSSolution(void);

	/*!
	 * \brief Compute the viscosity at the infinity.
	 * \return Value of the viscosity at the infinity.
	 */
	double GetViscosity_Inf(void);

	/*!
	 * \brief Compute the time step for solving the Navier-Stokes equations with turbulence model.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Compute the viscous forces and all the addimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Viscous_Forces(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Get the non dimensional lift coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	double GetCLift_Visc(unsigned short val_marker);

	/*!
	 * \brief Get the non dimensional drag coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	double GetCDrag_Visc(unsigned short val_marker);

	/*!
	 * \brief Get the total non dimensional lift coefficient (viscous contribution).
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	double GetAllBound_CLift_Visc(void);

	/*!
	 * \brief Get the total non dimensional drag coefficient (viscous contribution).
	 * \return Value of the drag coefficient (viscous contribution).
	 */
	double GetAllBound_CDrag_Visc(void);

	/*!
	 * \brief Compute the viscous residuals.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the skin friction coefficient.
	 */
	double GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex);
};

/*! 
 * \class CTurbSolution
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.0.
 */
class CTurbSolution : public CSolution {
protected:
	double *FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j;		/*!< \brief Store the flow solution at point j. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CTurbSolution(void);

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CTurbSolution(void);

	/*! 
	 * \brief Constructor of the class. 
	 */
	CTurbSolution(CConfig *config);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the interprocessor boundary condition (dirichlet).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

};

/*!
 * \class CTurbSASolution
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.0.
 */

class CTurbSASolution: public CTurbSolution {
private:
	double nu_tilde_Inf;
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSASolution(void);

	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSASolution(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSASolution(void);

	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short iMesh);

	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);
};

/*! 
 * \class CTurbSSTSolution
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Campos.
 * \version 1.0.
 */

class CTurbSSTSolution: public CTurbSolution {
private:
	double kine_Inf, omega_Inf;
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSSTSolution(void);

	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSSTSolution(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSSTSolution(void);

	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short iMesh);

	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);
};

/*!
 * \class CAdjEulerSolution
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.0.
 */
class CAdjEulerSolution : public CSolution {
protected:
	double PsiRho_Inf,	/*!< \brief PsiRho variable at the infinity. */
	PsiE_Inf,			/*!< \brief PsiE variable at the infinity. */
	*Phi_Inf;			/*!< \brief Phi vector at the infinity. */
	double *CSens_Mach, /*!< \brief Mach sensitivity coefficient for each boundary. */ 
	*CSens_AoA,			/*!< \brief Angle of attack sensitivity coefficient for each boundary. */ 
	*CSens_Geo,			/*!< \brief Shape sensitivity coefficient for each boundary. */
	**CSensitivity;		/*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
	double Total_CSens_Mach;	/*!< \brief Total mach sensitivity coefficient for all the boundaries. */
	double Total_CSens_AoA;		/*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
	double Total_CSens_Geo;		/*!< \brief Total shape sensitivity coefficient for all the boundaries. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CAdjEulerSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjEulerSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CAdjEulerSolution(void);

	/*! 
	 * \brief Created the force projection vector for adjoint boundary conditions.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	void SetForceProj_Vector(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*! 
	 * \brief Compute the jump for the interior boundary problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	void SetIntBoundary_Jump(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Compute adjoint density at the infinity.
	 * \return Value of the adjoint density at the infinity.
	 */
	double GetPsiRho_Inf(void);

	/*!
	 * \brief Compute the adjoint energy at the infinity.
	 * \return Value of the adjoint energy at the infinity.
	 */
	double GetPsiE_Inf(void);

	/*!
	 * \brief Compute Phi (adjoint velocity) at the infinity.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the adjoint velocity vector at the infinity.
	 */
	double GetPhi_Inf(unsigned short val_dim);

	/*!
	 * \brief Compute the spatial integration using a centred scheme for the adjoint equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the undivided laplacian for the adjoint solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Impose via the residual the adjoint Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose via the residual the interface adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose via the residual the near-field adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the interprocessor boundary condition (dirichlet).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose via the residual the adjoint symmetry boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using a Runge-Kutta strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Update the solution using a explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Initialize the residual vectors.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep);

	/*!
	 * \brief Compute the inviscid sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Get the shape sensitivity coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the sensitivity coefficient.
	 */
	double GetCSensitivity(unsigned short val_marker, unsigned short val_vertex);

	/*!
	 * \brief Set the shape sensitivity coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \param[in] val_sensitivity - Value of the sensitivity coefficient.
	 */
	void SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity);

	/*!
	 * \brief Provide the total shape sensitivity coefficient.
	 * \return Value of the geometrical sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_CSens_Geo(void);

	/*!
	 * \brief Set the total Mach number sensitivity coefficient.
	 * \return Value of the Mach sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_CSens_Mach(void);

	/*!
	 * \brief Set the total angle of attack sensitivity coefficient.
	 * \return Value of the angle of attack sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_CSens_AoA(void);
};

/*! 
 * \class CAdjNSSolution
 * \brief Main class for defining the Navier-Stokes' adjoint flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 1.0.
 */
class CAdjNSSolution : public CAdjEulerSolution {
public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CAdjNSSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjNSSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAdjNSSolution(void);

	/*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep);

	/*!
	 * \brief Compute the viscous sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Viscous_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Compute the viscous residuals for the adjoint equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Conservative source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourceConserv_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

};

/*! 
 * \class CAdjTurbSolution
 * \brief Main class for defining the adjoint turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.0.
 */
class CAdjTurbSolution : public CSolution {
private:
	double PsiNu_Inf,	/*!< \brief PsiNu variable at the infinity. */
	*FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j;	/*!< \brief Store the flow solution at point j. */

public:

	/*! 
	 * \brief Default constructor of the class. 
	 */
	CAdjTurbSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjTurbSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Default destructor of the class. 
	 */
	virtual ~CAdjTurbSolution(void);

	/*! 
	 * \brief Impose the Navier-Stokes turbulent adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*! 
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*! 
	 * \brief Initializate the residual vectors.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*! 
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh);

	/*! 
	 * \brief Compute the viscous residuals for the turbulent adjoint equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh, unsigned short iRKStep);

	/*! 
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*! 
	 * \brief Conservative source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourceConserv_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*! 
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);
};

/*! 
 * \class CLinEulerSolution
 * \brief Main class for defining the linearized Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.0.
 */
class CLinEulerSolution : public CSolution {
private:
	double DeltaRho_Inf,	/*!< \brief Linearized density variable at the infinity. */
	DeltaE_Inf,				/*!< \brief Linearized energy at the infinity. */
	*DeltaVel_Inf;			/*!< \brief Linearized velocity vector at the infinity. */
	double *p1_Und_Lapl,	/*!< \brief Undivided Laplacians for centred scheme. */ 
	*p2_Und_Lapl;			/*!< \brief Undivided Laplacians for centred scheme. */ 
	double *CDeltaDrag_Inv, /*!< \brief Linearized drag coefficient (inviscid contribution) for each boundary. */
	*CDeltaLift_Inv,		/*!< \brief Linearized lift coefficient (inviscid contribution) for each boundary. */
	*DeltaForceInviscid;	/*!< \brief Linearized inviscid force for each boundary. */
	double AllBound_CDeltaDrag_Inv, /*!< \brief Total linearized drag coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CDeltaLift_Inv;		/*!< \brief Total linearized lift coefficient (inviscid contribution) for all the boundaries. */
	double Total_CDeltaDrag,	/*!< \brief Total linearized drag coefficient for all the boundaries. */
	Total_CDeltaLift;			/*!< \brief Total linearized lift coefficient for all the boundaries. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CLinEulerSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CLinEulerSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CLinEulerSolution(void);

	/*!
	 * \brief Compute the spatial integration using a centred scheme for the linearized equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Compute the undivided laplacian for the linearized solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Impose via the residual the linearized Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the interprocessor boundary condition (dirichlet).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Restart residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep);

	/*!
	 * \brief Compute the linearization of the pressure forces and all the adimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_DeltaForces(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Get the total non dimensional drag coefficient.
	 * \return Value of the linearized drag coefficient (inviscid contribution).
	 */
	double GetTotal_CDeltaDrag(void);

	/*!
	 * \brief Get the total non dimensional lift coefficient.
	 * \return Value of the linearized lift coefficient (inviscid contribution).
	 */
	double GetTotal_CDeltaLift(void);
};

/*! 
 * \class CEikonalSolution
 * \brief Main class for defining the Eikonal solver.
 * \author A. Bueno.
 * \version 1.0.
 */
class CEikonalSolution : public CSolution {
public:

	/*! 
	 * \brief Default constructor of the class. 
	 */
	CEikonalSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CEikonalSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Default destructor of the class. 
	 */
	~CEikonalSolution(void);

	/*!
	 * \brief Compute the distance a brute force strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetDistance(CGeometry *geometry, CConfig *config);

};

/*! 
 * \class CPotentialSolution
 * \brief Main class for defining the potential model solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios.
 * \version 1.0.
 */
class CPotentialSolution : public CSolution {
private:
	double Density_Inf,		/*!< \brief Density at the infinity. */
	Energy_Inf,				/*!< \brief Energy at the infinity. */
	Pressure_Inf,			/*!< \brief Pressure at the infinity. */
	*Velocity_Inf;			/*!<\brief  Flow Velocity vector at the infinity. */
	double *FlowSolution_i,		/*!<\brief  Flow Solution at point i. */
	*FlowSolution_j;			/*!<\brief  Flow Solution at point j. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CPotentialSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPotentialSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CPotentialSolution(void);

	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the density x velocity at the infinity.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the density x velocity at the infinity.
	 */
	double GetDensity_Velocity_Inf(unsigned short val_dim);

	/*!
	 * \brief Compute the density x velocity at the infinity.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	double GetVelocity_Inf(unsigned short val_dim);

	/*!
	 * \brief Impose the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the far-field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the symmetry plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep);
};

/*! \class CElectricSolution
 *  \brief Main class for defining the electric potential solver.
 *  \author F. Palacios.
 *  \version 1.0.
 *  \date May 3, 2010.
 */
class CElectricSolution : public CSolution {
private:
	double Total_CCharge;	/*!< \brief Total charge coefficient for all the domain. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CElectricSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CElectricSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CElectricSolution(void);

	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose a Neumann boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);
	/*!
	 * \brief Compute the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep);

	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Get the value of the charge coefficient.
	 * \return Value of the charge coefficient.
	 */
	double GetTotal_CCharge(void);

	/*!
	 * \brief Set the value of the charge coefficient.
	 * \param[in] val_Total_CCharge - Value of the total charge coefficient.
	 */
	void SetTotal_CCharge(double val_Total_CCharge);
};

/*! 
 * \class CAdjElectricSolution
 * \brief Main class for defining the adjoint electric potential solver.
 * \author F. Palacios.
 * \version 1.0.
 */
class CAdjElectricSolution : public CSolution {
private:
	double Total_CCharge;	/*!< \brief Total charge coefficient for all the domain. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CAdjElectricSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjElectricSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CAdjElectricSolution(void);

	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */	
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);
};

/*! 
 * \class CLinElectricSolution
 * \brief Main class for defining the linearized electric potential solver.
 * \author F. Palacios.
 * \version 1.0.
 */
class CLinElectricSolution : public CSolution {
private:
	double Total_CCharge;	/*!< \brief Total charge coefficient for all the domain. */

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CLinElectricSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CLinElectricSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CLinElectricSolution(void);

	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */	
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);
};


/*! 
 * \class CPlasmaMonatomicSolution
 * \brief Main class for defining the plasma solver.
 * \author ADL Stanford.
 * \version 1.0.
 */
class CPlasmaMonatomicSolution : public CSolution {
protected:
	double *Density_Inf;		/*!< \brief Density of all species at the infinity. */
	double *Density_Inlet;		/*!< \brief Density of all species at the inlet. */
	double *Density_Outlet;		/*!< \brief Density of all species at the outlet. */

	double *Mach_Inf;			/*!< \brief Mach number of all fluids at the infinity. */
	double *Mach_Inlet;			/*!< \brief Mach number of all fluids at the inlet. */
	double *Mach_Outlet;		/*!< \brief Mach number of all fluids at the outlet. */

	double *Energy_Inf;			/*!< \brief Energy of all fluids at the infinity. */
	double *Energy_Inlet;		/*!< \brief Energy of all fluids at the inlet. */
	double *Energy_Outlet;		/*!< \brief Energy of all fluids at the outlet. */

	double *Pressure_Inf;		/*!< \brief Pressure of all fluids at the infinity. */
	double *Pressure_Inlet;		/*!< \brief Pressure of all fluids at the inlet. */
	double *Pressure_Outlet;	/*!< \brief Pressure of all fluids at the outlet. */

	double **Velocity_Inf;		/*!< \brief Flow Velocity vector of all species at the infinity. */
	double **Velocity_Inlet;	/*!< \brief Flow Velocity vector of all species at the inlet. */
	double **Velocity_Outlet;	/*!< \brief Flow Velocity vector of all species at the outlet. */

	double *p1_Und_Lapl;		/*!< \brief Auxiliary variable for the undivided Laplacians. */ 
	double *p2_Und_Lapl;		/*!< \brief Auxiliary variable for the undivided Laplacians. */ 
	double *PrimVar_i;			/*!< \brief Auxiliary vector for storing the solution at point i. */
	double *PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	//double ***CPressure;


public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CPlasmaMonatomicSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPlasmaMonatomicSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CPlasmaMonatomicSolution(void);

	/*!
	 * \brief Compute the density at the infinity.
	 * \param[in] val_var - Index of the variable for the density.
	 * \return Value of the density at the infinity.
	 */
	double GetDensity_Inf(unsigned short val_var);

	/*!
	 * \brief Compute the density multiply by energy at the infinity.
	 * \param[in] val_var - Index of the variable for the energy.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	double GetDensity_Energy_Inf(unsigned short val_var);

	/*!
	 * \brief Compute the density multiply by velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \param[in] val_var - Index of the variable for the velocity.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	double GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var);

	/*!
	 * \brief Get the density at the outlet for fluid val_Fluid.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density at the outlet.
	 */
	double GetDensity_Outlet(unsigned short val_Fluid);

	/*!
	 * \brief Get the density at the inlet for fluid val_Fluid.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density at the outlet.
	 */
	double GetDensity_Inlet(unsigned short val_Fluid);

	/*!
	 * \brief Get the density*Energy at the outlet for fluid val_Fluid.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density*Energy at the outlet.
	 */
	double GetDensity_Energy_Outlet(unsigned short val_Fluid);

	/*!
	 * \brief Get the density*Energy at the inlet for fluid val_Fluid.
	 * \param[in] val_Fluid - Index of the fluid.
	 * \return Value of the density*Energy at the outlet.
	 */
	double GetDensity_Energy_Inlet(unsigned short val_Fluid);

	/*!
	 * \brief Get the density*Velocity at the outlet for fluid val_Fluid in val_dim direction
	 * \param[in] val_Fluid - Index of the fluid.
	 * \param[in] val_dim - Index of the direction.
	 * \return Value of the density*Velocity at the outlet.
	 */
	double GetDensity_Velocity_Outlet(unsigned short val_dim,unsigned short val_Fluid);

	/*!
	 * \brief Get the density*Velocity at the inlet for fluid val_Fluid in val_dim direction
	 * \param[in] val_Fluid - Index of the fluid.
	 * \param[in] val_dim - Index of the direction.
	 * \return Value of the density*Velocity at the outlet.
	 */
	double GetDensity_Velocity_Inlet(unsigned short val_dim,unsigned short val_Fluid);

	/*!
	 * \brief Compute the time step for solving the Euler equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the spatial integration using a centred scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spectral radius.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSpectral_Radius(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute a pressure sensor switch.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetPress_Switch(CGeometry *geometry);

	/*!
	 * \brief Compute the gradient of the primitive variables using Green-Gauss method, 
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method, 
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the undivided laplacian for the solution, except the energy equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the symmetry boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Do the send-receive communications in the MPI parallelization.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Set the boundary contition of the interprocessor boundaries.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);


	/*!
	 * \brief Impose the subsonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Electrode(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the subsonic outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.

	 */
	void BC_Dielectric(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.

	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker) ;
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.

	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker) ;

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Compute the pressure forces and all the adimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Electrostatic_Forces(CGeometry *geometry, CConfig *config);

};

/*! 
 * \class CLevelSetSolution
 * \brief Main class for defining the level set solver.
 * \ingroup LevelSet_Model
 * \author F. Palacios.
 * \version 1.0.
 */
class CLevelSetSolution : public CSolution {
protected:
	double *FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j;		/*!< \brief Store the flow solution at point j. */
	double levelset_Inf;

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CLevelSetSolution(void);

	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CLevelSetSolution(CGeometry *geometry, CConfig *config);	

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CLevelSetSolution(void);

	/*! 
	 * \brief Constructor of the class. 
	 */
	CLevelSetSolution(CConfig *config);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the interprocessor boundary condition (dirichlet).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short iMesh);

	/*!
	 * \brief Impose the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

};


/*! 
 * \class CCombustionSolution
 * \brief Main class for defining the combustion solver.
 * \version 1.0.
 */
class CCombustionSolution : public CSolution {
protected:

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CCombustionSolution(void);

	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCombustionSolution(CGeometry *geometry, CConfig *config);	

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CCombustionSolution(void);

	/*! 
	 * \brief Constructor of the class. 
	 */
	CCombustionSolution(CConfig *config);

	/*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the interprocessor boundary condition (dirichlet).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short iMesh);

	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Impose the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config,
			unsigned short val_marker);

	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
			unsigned short val_marker);

};

/*! 
 * \class CTemplateSolution
 * \brief Main class for defining the template model solver.
 * \ingroup Template_Flow_Equation
 * \author F. Palacios.
 * \version 1.0.
 */
class CTemplateSolution : public CSolution {
private:

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CTemplateSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTemplateSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CTemplateSolution(void);

	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep);

	/*!
	 * \brief Compute the time step for solving the Euler equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the spatial integration using a centred scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Impose the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the far-field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the symmetry plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Impose the dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

};

/*! 
 * \class CPlasmaDiatomicSolution
 * \brief Main class for defining the plasma solver.
 * \author ADL Stanford.
 * \version 1.0.
 */
class CPlasmaDiatomicSolution : public CSolution {
protected:
	double *Gamma;					/*! \brief Ratio of specific heats for each species. */
	double *Gas_Constant;		/*! \brief Gas constant for each species. */
	double *Molecular_Mass;	/*! \brief Molecular mass of each species. */
	double *Energy_Formation;		/*! \brief Chemical energy of formation for all species. */
	double *Enthalpy_Formation;		/*! \brief Chemical enthalpy of formation for all species. */

	double *Density_Inf;		/*!< \brief Density of all species at the infinity. */
	double *Mach_Inf;			/*!< \brief Mach number of all fluids at the infinity. */
	double *Energy_Inf;			/*!< \brief Energy of all fluids at the infinity. */
	double *Energy_vib_Inf;			/*!< \brief Energy of all fluids at the infinity. */
	double *Pressure_Inf;		/*!< \brief Pressure of all fluids at the infinity. */
	double **Velocity_Inf;		/*!< \brief Flow Velocity vector of all species at the infinity. */

	double *PrimVar_i;			/*!< \brief Auxiliary vector for storing the solution at point i. */
	double *PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */	

public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CPlasmaDiatomicSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPlasmaDiatomicSolution(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CPlasmaDiatomicSolution(void);

	/*!
	 * \brief Compute the density at the infinity.
	 * \param[in] val_var - Index of the variable for the density.
	 * \return Value of the density at the infinity.
	 */
	double GetDensity_Inf(unsigned short val_var);

	/*!
	 * \brief Compute the density multiply by energy at the infinity.
	 * \param[in] val_var - Index of the variable for the energy.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	double GetDensity_Energy_Inf(unsigned short val_var);

	/*!
	 * \brief Compute the density multiply by energy at the infinity.
	 * \param[in] val_var - Index of the variable for the energy.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	double GetDensity_Energy_vib_Inf(unsigned short val_var);

	/*!
	 * \brief Compute the density multiply by velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \param[in] val_var - Index of the variable for the velocity.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	double GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var);

	/*!
	 * \brief Compute the time step for solving the Euler equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iMesh);

	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
			CConfig *config, unsigned short iMesh);

	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Compute the spectral radius.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSpectral_Radius(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the gradient of the primitive variables using Green-Gauss method, 
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method, 
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the undivided laplacian for the solution, except the energy equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */	
	void BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */	
	void BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker);

	/*!
	 * \brief Do the send-receive communications in the MPI parallelization.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Set the boundary contition of the interprocessor boundaries.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in] val_mesh - Index of the mesh in multigrid computations.
	 */
	void BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short val_marker, unsigned short val_mesh);

	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
			CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
			unsigned short iRKStep);

	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config);

	/*!
	 * \brief Compute the pressure forces and all the adimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Electrostatic_Forces(CGeometry *geometry, CConfig *config);

};

#include "solution_structure.inl"
