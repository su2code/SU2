/*!
 * \file solver_structure.hpp
 * \brief Headers of the main subroutines for solving partial differential equations.
 *        The subroutines and functions are in the <i>solver_structure.cpp</i>,
 *        <i>solution_direct.cpp</i>, <i>solution_adjoint.cpp</i>, and
 *        <i>solution_linearized.cpp</i> files.
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
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "numerics_structure.hpp"
#include "variable_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"

using namespace std;

/*!
 * \class CSolver
 * \brief Main class for defining the PDE solution, it requires
 * a child class for each particular solver (Euler, Navier-Stokes, Plasma, etc.)
 * \author F. Palacios.
 * \version 2.0.6
 */
class CSolver {
protected:
	unsigned short IterLinSolver;	/*!< \brief Linear solver iterations. */
	unsigned short nVar,					/*!< \brief Number of variables of the problem. */
  nPrimVar,                     /*!< \brief Number of primitive variables of the problem. */
  nPrimVarGrad,                    /*!< \brief Number of primitive variables of the problem in the gradient computation. */
	nDim;													/*!< \brief Number of dimensions of the problem. */
	unsigned long nPoint;					/*!< \brief Number of points of the computational grid. */
  unsigned long nPointDomain; 	/*!< \brief Number of points of the computational grid. */
	double Max_Delta_Time,	/*!< \brief Maximum value of the delta time for all the control volumes. */
	Min_Delta_Time;					/*!< \brief Minimum value of the delta time for all the control volumes. */
	double *Residual_RMS,	/*!< \brief Vector with the mean residual for each variable. */
  *Residual_Max,        /*!< \brief Vector with the maximal residual for each variable. */
	*Residual,						/*!< \brief Auxiliary nVar vector. */
	*Residual_i,					/*!< \brief Auxiliary nVar vector for storing the residual at point i. */
	*Residual_j;					/*!< \brief Auxiliary nVar vector for storing the residual at point j. */
  unsigned long *Point_Max; /*!< \brief Vector with the maximal residual for each variable. */
	double *Solution,		/*!< \brief Auxiliary nVar vector. */
	*Solution_i,				/*!< \brief Auxiliary nVar vector for storing the solution at point i. */
	*Solution_j;				/*!< \brief Auxiliary nVar vector for storing the solution at point j. */
	double *Vector,	/*!< \brief Auxiliary nDim vector. */
	*Vector_i,			/*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point i. */
	*Vector_j;			/*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point j. */
	double *Res_Conv,	/*!< \brief Auxiliary nVar vector for storing the convective residual. */
	*Res_Visc,				/*!< \brief Auxiliary nVar vector for storing the viscous residual. */
	*Res_Sour,				/*!< \brief Auxiliary nVar vector for storing the viscous residual. */
	*Res_Conv_i,		  /*!< \brief Auxiliary vector for storing the convective residual at point i. */
	*Res_Visc_i,			/*!< \brief Auxiliary vector for storing the viscous residual at point i. */
	*Res_Conv_j,			/*!< \brief Auxiliary vector for storing the convective residual at point j. */
	*Res_Visc_j;			/*!< \brief Auxiliary vector for storing the viscous residual at point j. */
	double **Jacobian_i,	/*!< \brief Auxiliary matrices for storing point to point Jacobians at point i. */
	**Jacobian_j;			    /*!< \brief Auxiliary matrices for storing point to point Jacobians at point j. */
	double **Jacobian_ii,	/*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_ij,			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_ji,			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_jj;			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  
	double **Smatrix,	/*!< \brief Auxiliary structure for computing gradients by least-squares */
	**cvector;			 /*!< \brief Auxiliary structure for computing gradients by least-squares */

public:
  
  CSysVector LinSysSol;		/*!< \brief vector to store iterative solution of implicit linear system. */
  CSysVector LinSysRes;		/*!< \brief vector to store iterative residual of implicit linear system. */
	CSysMatrix Jacobian; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
  
	CSysMatrix StiffMatrix; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations, and grid movement. */

	CVariable** node;	/*!< \brief Vector which the define the variables for each problem. */
  
	/*!
	 * \brief Constructor of the class.
	 */
	CSolver(void);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CSolver(void);
     
  /*!
	 * \brief Set number of linear solver iterations.
	 * \param[in] val_iterlinsolver - Number of linear iterations.
	 */
	void SetIterLinSolver(unsigned short val_iterlinsolver);
    
	/*!
	 * \brief Set number of linear solver iterations.
	 * \param[in] val_iterlinsolver - Number of linear iterations.
	 */
	virtual void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Set the value of the max residual and RMS residual.
	 * \param[in] val_iterlinsolver - Number of linear iterations.
	 */
	void SetResidual_RMS(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Set number of linear solver iterations.
	 * \param[in] val_iterlinsolver - Number of linear iterations.
	 */
	virtual void Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Get number of linear solver iterations.
	 * \return Number of linear solver iterations.
	 */
	unsigned short GetIterLinSolver(void);
    
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
	 * \brief Get the value of the maximum delta time.
	 * \return Value of the maximum delta time.
	 */
	virtual double GetMax_Delta_Time(unsigned short val_Species);
    
	/*!
	 * \brief Get the value of the minimum delta time.
	 * \return Value of the minimum delta time.
	 */
	virtual double GetMin_Delta_Time(unsigned short val_Species);
    
	/*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnVar(void);
    
    /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnPrimVar(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	virtual void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Set the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void SetRes_RMS(unsigned short val_var, double val_residual);
    
	/*!
	 * \brief Adds the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void AddRes_RMS(unsigned short val_var, double val_residual);
    
	/*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
	 */
	double GetRes_RMS(unsigned short val_var);
    
    /*!
	 * \brief Set the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void SetRes_Max(unsigned short val_var, double val_residual, unsigned long val_point);
    
	/*!
	 * \brief Adds the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void AddRes_Max(unsigned short val_var, double val_residual, unsigned long val_point);
    
	/*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
	 */
	double GetRes_Max(unsigned short val_var);
    
    /*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
	 */
	unsigned long GetPoint_Max(unsigned short val_var);
    
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
	void SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the Least Squares gradient of the solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief MPI gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    virtual void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Compute the Least Squares gradient of the grid velocity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetRotVel_Gradient(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the Least Squares gradient of the solution on the profile surface.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSurface_Gradient(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute slope limiter.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSolution_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPrimVar_Limiter(CGeometry *geometry, CConfig *config);
    
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
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] Iteration - Index of the current iteration.
	 */
	virtual void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iMesh, unsigned long Iteration);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	virtual void Preprocessing(CGeometry *geometry, CSolver **solver_container,
                               CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetDissipation_Switch(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                               unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                 unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_FlowLoad(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                    unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Definition of hte solver settings.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Dielectric(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                              CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_FWH(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                        unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Observer(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iRKStep);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config,
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
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	virtual void Wave_Strength(CGeometry *geometry, CConfig *config);
    
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
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPrimVar_Limiter_MPI(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] iPoint - Index of the grid point.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPreconditioner(CConfig *config, unsigned short iPoint);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                  CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] StiffMatrix_Elem - Stiffness matrix of an element
	 */
	virtual void AddStiffMatrix(double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3 );
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                 CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void SourceConserv_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Charge_Dist_SourceTerm(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
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
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	 * \param[in] val_Total_CT - Value of the total thrust coefficient.
	 */
	virtual void SetTotal_CT(double val_Total_CT);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CQ - Value of the total torque coefficient.
	 */
	virtual void SetTotal_CQ(double val_Total_CQ);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_Total_Q - Value of the total heat load.
	 */
	virtual void SetTotal_Q(double val_Total_Q);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_Total_MaxQ - Value of the total heat load.
	 */
	virtual void SetTotal_MaxQ(double val_Total_MaxQ);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetDistance(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
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
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	virtual double GetFanFace_MassFlow(unsigned short val_marker);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	virtual double GetExhaust_MassFlow(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the fan face pressure on the surface <i>val_marker</i>.
	 */
	virtual double GetFanFace_Pressure(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the fan face mach on the surface <i>val_marker</i>.
	 */
	virtual double GetFanFace_Mach(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCSideForce_Inv(unsigned short val_marker);
    
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
	 * \return Value of the efficiency coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CEff(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the thrust coefficient (force in the -x direction, inviscid + viscous contribution).
	 */
	virtual double GetTotal_CT(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the torque coefficient (moment in the -x direction, inviscid + viscous contribution).
	 */
	virtual double GetTotal_CQ(void);
    
    /*!
	 * \brief A virtual member.
	 * \return Value of the heat load (integrated heat flux).
	 */
	virtual double GetTotal_Q(void);
    
    /*!
	 * \brief A virtual member.
	 * \return Value of the heat load (integrated heat flux).
	 */
	virtual double GetTotal_MaxQ(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual double Get_PressureDrag(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual double Get_ViscDrag(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual double Get_MagnetDrag(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the rotor Figure of Merit (FM) (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CMerit(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CEquivArea(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Free Surface coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CFreeSurface(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the FEA coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CFEA(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Near-Field Pressure coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CNearFieldOF(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	virtual void SetTotal_CEquivArea(double val_cequivarea);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cfreesurface - Value of the Free Surface coefficient.
	 */
	virtual void SetTotal_CFreeSurface(double val_cfreesurface);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cfea - Value of the FEA coefficient.
	 */
	virtual void SetTotal_CFEA(double val_cfea);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
	 */
	virtual void SetTotal_CNearFieldOF(double val_cnearfieldpress);
    
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
	 * \return Value of the force x coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CFx(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CFy(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	virtual double GetTotal_CFz(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the wave strength.
	 */
	virtual double GetTotal_CWave(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the wave strength.
	 */
	virtual double GetTotal_CHeat(void);
    
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
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual double GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_iSpecies - index of the chemical species
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual double GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_iSpecies, unsigned short val_vertex);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_iSpecies - index of the chemical species
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual double GetViscForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_iSpecies - index of the chemical species
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual double GetPressureForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex);
    
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the y plus.
	 */
	virtual double GetYPlus(unsigned short val_marker, unsigned short val_vertex);
	
    
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
	virtual double GetTotal_Sens_Geo(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Mach sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_Sens_Mach(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the angle of attack sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_Sens_AoA(void);
    
	/*!
	 * \brief Set the total farfield pressure sensitivity coefficient.
	 * \return Value of the farfield pressure sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_Sens_Press(void);
    
	/*!
	 * \brief Set the total farfield temperature sensitivity coefficient.
	 * \return Value of the farfield temperature sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual double GetTotal_Sens_Temp(void);
    
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
	 * \brief Compute the density multiply by energy at the infinity.
	 * \param[in] val_var - Index of the variable for the energy.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	virtual double GetDensity_Energy_vib_Inf(unsigned short val_var);
    
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
	virtual double GetDensity_Velocity_Outlet(unsigned short val_dim, unsigned short val_Fluid);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetLevelSet_Distance(CGeometry *geometry, CConfig *config, bool Initialization, bool WriteLevelSet);
    
	/*!
	 * \brief A virtual member.
	 * \return A pointer to an array containing a set of constants
	 */
	virtual double* GetConstants();
    
	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] fea_config - Definition of the particular problem.
	 */
	virtual void SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] solver1_geometry - Geometrical definition of the problem.
	 * \param[in] solver1_solution - Container vector with all the solutions.
	 * \param[in] solver1_config - Definition of the particular problem.
	 * \param[in] solver2_geometry - Geometrical definition of the problem.
	 * \param[in] solver2_solution - Container vector with all the solutions.
	 * \param[in] solver2_config - Definition of the particular problem.
	 */
	virtual void Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config,CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] fea_config - Definition of the particular problem.
	 */
	virtual void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] flow_geometry - Geometrical definition of the problem.
	 * \param[in] flow_grid_movement - Geometrical definition of the problem.
	 * \param[in] flow_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement, CConfig *flow_config,
                                      CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] wave_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] wave_config - Definition of the particular problem.
	 */
	virtual void SetNoise_Source(CSolver ***flow_solution, CGeometry **wave_geometry, CConfig *wave_config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] wave_solution - Solution container from the wave problem.
	 * \param[in] flow_solution - Solution container from the flow problem.
	 * \param[in] flow_geometry - Geometrical definition for the flow problem.
	 * \param[in] flow_config - Definition of the particular problem.
	 */
	virtual void SetAeroacoustic_Coupling(CSolver ***wave_solution, CSolver ***flow_solution, CNumerics *numerics,
                                          CGeometry **flow_geometry, CConfig *flow_config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iZone - Current zone in the mesh.
	 */
	virtual void GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone);
    
	/*!
	 * \brief Gauss method for solving a linear system.
	 * \param[in] A - Matrix Ax = b.
	 * \param[in] rhs - Right hand side.
	 * \param[in] nVar - Number of variables.
	 */
	void Gauss_Elimination(double** A, double* rhs, unsigned long nVar);
    
    /*!
	 * \brief Get the number of Species present in the flow.
	 */
	virtual unsigned short GetnSpecies(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	virtual unsigned short GetnMonatomics(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	virtual unsigned short GetnDiatomics(void);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution - Container vector with all the solutions.
	 */
	virtual void GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh);

};

/*!
 * \class CBaselineSolver
 * \brief Main class for defining a baseline solution from a restart file (for output).
 * \author F. Palacios, T. Economon.
 * \version 2.0.6
 */
class CBaselineSolver : public CSolver {
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CBaselineSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Load a solution from a restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iZone - Current zone in the mesh.
	 */
	void GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CBaselineSolver(void);
    
};

/*!
 * \class CEulerSolver
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 2.0.6
 */
class CEulerSolver : public CSolver {
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
	*CMx_Inv,			/*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
	*CMy_Inv,			/*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
	*CMz_Inv,			/*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
	*CFx_Inv,			/*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
	*CFy_Inv,			/*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
	*CFz_Inv,			/*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
	*CEff_Inv,				/*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
	*CMerit_Inv,				/*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
	*CT_Inv,			/*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
	*CQ_Inv,			/*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
	*CEquivArea_Inv,				/*!< \brief Equivalent area (inviscid contribution) for each boundary. */
	*CNearFieldOF_Inv,				/*!< \brief Near field pressure (inviscid contribution) for each boundary. */
	**CPressure,		/*!< \brief Pressure coefficient for each boundary and vertex. */
	**CHeatTransfer,		/*!< \brief Heat transfer coefficient for each boundary and vertex. */
	**YPlus,		/*!< \brief Yplus for each boundary and vertex. */
	*ForceInviscid,		/*!< \brief Inviscid force for each boundary. */
	*MomentInviscid,	/*!< \brief Inviscid moment for each boundary. */
	*FanFace_MassFlow,	/*!< \brief Mass flow rate for each boundary. */
	*Exhaust_MassFlow,	/*!< \brief Mass flow rate for each boundary. */
	*FanFace_Pressure,	/*!< \brief Fan face pressure for each boundary. */
	*FanFace_Mach,	/*!< \brief Fan face mach number for each boundary. */
	*FanFace_Area,	/*!< \brief Boundary total area. */
    *Exhaust_Area,	/*!< \brief Boundary total area. */
    FanFace_MassFlow_Total,	/*!< \brief Mass flow rate for each boundary. */
    Exhaust_MassFlow_Total,	/*!< \brief Mass flow rate for each boundary. */
	FanFace_Pressure_Total,	/*!< \brief Fan face pressure for each boundary. */
	FanFace_Mach_Total,	/*!< \brief Fan face mach number for each boundary. */
	InverseDesign;	/*!< \brief Inverse design functional for each boundary. */
	double AllBound_CDrag_Inv,	/*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CLift_Inv,			/*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CSideForce_Inv,			/*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMx_Inv,			/*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMy_Inv,			/*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMz_Inv,			/*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFx_Inv,			/*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFy_Inv,			/*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFz_Inv,			/*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEff_Inv,			/*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMerit_Inv,			/*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
	AllBound_CT_Inv,			/*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CQ_Inv,			/*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEquivArea_Inv,			/*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CNearFieldOF_Inv;			/*!< \brief Near-Field press coefficient (inviscid contribution) for all the boundaries. */
	double Total_CDrag, /*!< \brief Total drag coefficient for all the boundaries. */
	Total_CLift,		/*!< \brief Total lift coefficient for all the boundaries. */
	Total_CSideForce,		/*!< \brief Total sideforce coefficient for all the boundaries. */
	Total_CMx,			/*!< \brief Total x moment coefficient for all the boundaries. */
	Total_CMy,			/*!< \brief Total y moment coefficient for all the boundaries. */
	Total_CMz,			/*!< \brief Total z moment coefficient for all the boundaries. */
	Total_CFx,			/*!< \brief Total x force coefficient for all the boundaries. */
	Total_CFy,			/*!< \brief Total y force coefficient for all the boundaries. */
	Total_CFz,			/*!< \brief Total z force coefficient for all the boundaries. */
	Total_CEff,			/*!< \brief Total efficiency coefficient for all the boundaries. */
	Total_CMerit,			/*!< \brief Total rotor Figure of Merit for all the boundaries. */
	Total_CT,		/*!< \brief Total thrust coefficient for all the boundaries. */
	Total_CQ,		/*!< \brief Total torque coefficient for all the boundaries. */
    Total_Q,    /*!< \brief Total heat load for all the boundaries. */
    Total_Maxq, /*!< \brief Maximum heat flux on all boundaries. */
	Total_CEquivArea,			/*!< \brief Total Equivalent Area coefficient for all the boundaries. */
	Total_CNearFieldOF;			/*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
	double *p1_Und_Lapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*p2_Und_Lapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	double *PrimVar_i,	/*!< \brief Auxiliary vector for storing the solution at point i. */
	*PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	double **Precon_Mat_inv; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
	bool space_centered,  /*!< \brief True if space centered scheeme used. */
	euler_implicit,			/*!< \brief True if euler implicit scheme used. */
	roe_turkel,         /*!< \brief True if computing preconditioning matrix for roe-turkel method. */
	least_squares;        /*!< \brief True if computing gradients by least squares. */
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CEulerSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CEulerSolver(void);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    virtual void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
    
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
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] Iteration - Value of the current iteration.
	 */
	void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh, unsigned long Iteration);
    
	/*!
	 * \brief Compute the spatial integration using a centered scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute a pressure sensor switch.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetDissipation_Switch(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Parallelization of SetDissipation_Switch.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config);
    
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
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the limiter of the primitive variables.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
	 * \param[in] iPoint - Index of the grid point
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPreconditioner(CConfig *config, unsigned short iPoint);
    
	/*!
	 * \brief Compute the undivided laplacian for the solution, except the energy equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Parallelization of Undivided Laplacian.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Compute the max eigenvalue.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Parallelization of the Max eigenvalue.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the interface boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the dirichlet boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose a subsonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose a supersonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                   CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
     
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                   CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the nacelle inflow boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the ancelle exhaust boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep);
    
    /*!
	 * \brief Compute the Fan face Mach number.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution - Container vector with all the solutions.
	 */
	void GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	double GetFanFace_MassFlow(unsigned short val_marker);
    
    /*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	double GetExhaust_MassFlow(unsigned short val_marker);
    
	/*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the fan face pressure on the surface <i>val_marker</i>.
	 */
	double GetFanFace_Pressure(unsigned short val_marker);
    
	/*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the fan face mach on the surface <i>val_marker</i>.
	 */
	double GetFanFace_Mach(unsigned short val_marker);
    
	/*!
	 * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCSideForce_Inv(unsigned short val_marker);
    
	/*!
	 * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCEff_Inv(unsigned short val_marker);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
	 * \return Value of the sideforce coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CSideForce(void);
    
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
	double GetTotal_CNearFieldOF(void);
    
	/*!
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_CEquivArea(double val_cequivarea);
    
	/*!
	 * \brief Set the value of the Near-Field pressure oefficient.
	 * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
	 */
	void SetTotal_CNearFieldOF(double val_cnearfieldpress);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
	 * \param[in] val_Total_CLift - Value of the total lift coefficient.
	 */
	void SetTotal_CLift(double val_Total_CLift);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
	 * \return Value of the lift coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CLift(void);
    
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
	 * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
	 * \return Value of the force x coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFx(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFy(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
	 * \return Value of the force z coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFz(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional thrust coefficient.
	 * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CT(void);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional thrust coefficient.
	 * \param[in] val_Total_CT - Value of the total thrust coefficient.
	 */
	void SetTotal_CT(double val_Total_CT);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional torque coefficient.
	 * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CQ(void);
    
    /*!
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	double GetTotal_Q(void);
    
    /*!
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	double GetTotal_MaxQ(void);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional torque coefficient.
	 * \param[in] val_Total_CQ - Value of the total torque coefficient.
	 */
	void SetTotal_CQ(double val_Total_CQ);
    
    /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Q - Value of the heat load.
	 */
	void SetTotal_Q(double val_Total_Q);
    
    /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Q - Value of the heat load.
	 */
	void SetTotal_MaxQ(double val_Total_MaxQ);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional rotor Figure of Merit.
	 * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CMerit(void);
    
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
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] flow_geometry - Geometrical definition of the problem.
	 * \param[in] flow_grid_movement - Geometrical definition of the problem.
	 * \param[in] flow_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	void SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement, CConfig *flow_config, CConfig *fea_config,
                              CGeometry **fea_geometry, CSolver ***fea_solution);
    
	/*!
	 * \brief Load a direct flow solution for use with the adjoint solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iZone - Current zone in the mesh.
	 */
	void GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone);
    
	/*!
	 * \brief Set the initial condition for the Euler Equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] ExtIter - External iteration.
	 */
	void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
    
};

/*!
 * \class CNSSolver
 * \brief Main class for defining the Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 2.0.6
 */
class CNSSolver : public CEulerSolver {
private:
	double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */
	double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
	Prandtl_Turb;         /*!< \brief Turbulent Prandtl number. */
	double *CDrag_Visc,	/*!< \brief Drag coefficient (viscous contribution) for each boundary. */
	*CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for each boundary. */
	*CMx_Visc,			/*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
	*CMy_Visc,			/*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
	*CMz_Visc,			/*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
	*CFx_Visc,			/*!< \brief Force x coefficient (viscous contribution) for each boundary. */
	*CFy_Visc,			/*!< \brief Force y coefficient (viscous contribution) for each boundary. */
	*CFz_Visc,			/*!< \brief Force z coefficient (viscous contribution) for each boundary. */
	*CEff_Visc,			/*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
	*CMerit_Visc,			/*!< \brief Rotor Figure of Merit (Viscous contribution) for each boundary. */
	*CT_Visc,		/*!< \brief Thrust coefficient (viscous contribution) for each boundary. */
	*CQ_Visc,		/*!< \brief Torque coefficient (viscous contribution) for each boundary. */
    *Q_Visc,		/*!< \brief Heat load (viscous contribution) for each boundary. */
    *Maxq_Visc, /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
    
	**CSkinFriction;	/*!< \brief Skin friction coefficient for each boundary and vertex. */
	double *ForceViscous,	/*!< \brief Viscous force for each boundary. */
	*MomentViscous;			/*!< \brief Inviscid moment for each boundary. */
	double AllBound_CDrag_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
	AllBound_CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
	AllBound_CMx_Visc,			/*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMy_Visc,			/*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMz_Visc,			/*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEff_Visc,			/*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
	AllBound_CFx_Visc,			/*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFy_Visc,			/*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFz_Visc,			/*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMerit_Visc,			/*!< \brief Rotor Figure of Merit coefficient (Viscous contribution) for all the boundaries. */
	AllBound_CT_Visc,		/*!< \brief Thrust coefficient (viscous contribution) for all the boundaries. */
	AllBound_CQ_Visc,		/*!< \brief Torque coefficient (viscous contribution) for all the boundaries. */
    AllBound_Q_Visc,		/*!< \brief Heat load (viscous contribution) for all the boundaries. */
    AllBound_Maxq_Visc; /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CNSSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CNSSolver(void);
    
	/*!
	 * \brief Compute the viscosity at the infinity.
	 * \return Value of the viscosity at the infinity.
	 */
	double GetViscosity_Inf(void);
    
	/*!
	 * \brief Compute the time step for solving the Navier-Stokes equations with turbulence model.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] Iteration - Index of the current iteration.
	 */
	void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh, unsigned long Iteration);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
    /*!
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
    /*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                            unsigned short val_marker);
    
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
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                          CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the skin friction coefficient.
	 */
	double GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex);
    
	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex);
	
	/*!
	 * \brief Get the y plus.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the y plus.
	 */
	double GetYPlus(unsigned short val_marker, unsigned short val_vertex);
};

/*!
 * \class CTurbSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 2.0.6
 */
class CTurbSolver : public CSolver {
protected:
	double *FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j,        /*!< \brief Store the flow solution at point j. */
	*lowerlimit,            /*!< \brief contains lower limits for turbulence variables. */
	*upperlimit;            /*!< \brief contains upper limits for turbulence variables. */
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSolver(void);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CTurbSolver(void);
    
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSolver(CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
    /*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Alternative Green Gauss calculation of gradient of solution. Required for discrete adjoint
	 * \param[in] U_i - Solution at i.
	 * \param[in] U_js - Solutions at js (neighbours of i).
	 * \param[in] nNeigh - Number of neighbours of i.
	 * \param[in] Normals - Normals from i to neighbours.
	 * \param[out] grad_U_i - gradient at i.
	 * \param[in] config - Definition of the particular problem.
	 */
    
	void CalcGradient_GG(double *val_U_i, double **val_U_js, unsigned short nNeigh, unsigned short numVar,
                         double **Normals, double **grad_U_i, CConfig *config, CGeometry *geometry, unsigned long iPoint);
    
	/*!
	 * \brief Alternative Least Squares calculation of gradient of solution. Required for discrete adjoint
	 * \param[in] U_i - Solution at i.
	 * \param[in] U_js - Solutions at js (neighbours of i).
	 * \param[in] nNeigh - Number of neighbours of i.
	 * \param[in] coords_i - coords of i.
	 * \param[in] coords_js - coords of neighbours.
	 * \param[out] grad_U_i - gradient at i.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CalcGradient_LS(double *U_i, double **U_js, unsigned long nNeigh,
                         double *coords_i, double **coords_js, double **grad_U_i, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 */
	void CalcPrimVar_Compressible(double *val_Vars, double Gamma, double Gas_Constant, unsigned short numVar, double turb_ke,
                                  double* Primitive, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 */
    //	void CalcPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf,
    //			double ArtComp_Factor, double turb_ke, bool freesurface,
    //			double* LaminarViscosityInc, double* Primitive);
    
	/*!
	 * \brief Calculate the laminar viscosity (used for AD).
	 * \param[in] val_U_i - Value of the flow variables at point i.
	 * \param[out] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CalcLaminarViscosity(double *val_U_i, double *val_laminar_viscosity_i, CConfig *config);
    
	/*!
	 * \brief Calculate the eddy viscosity (used for AD).
	 */
	virtual void CalcEddyViscosity(double *val_FlowVars, double val_laminar_viscosity,
                                   double *val_TurbVar, double *val_eddy_viscosity);
    
    //	/*!
    //	 * \brief Alternative MUSCL reconstruction of solution. Required for discrete adjoint
    //	 * \param[in] val_grad_U_i - Solution gradient at i.
    //	 * \param[in] val_grad_U_j - Solution gradient at j.
    //	 * \param[in] Vector_i - 1/2 Vector i -> j
    //	 * \param[in] Vector_j - 1/2 Vector j -> i
    //	 * \param[in] config - Definition of the particular problem.
    //	 */
    //	virtual void MUSCL_Reconstruction(double **val_grad_U_i, double **val_grad_U_j,
    //			double *Vector_i, double *Vector_j, CConfig *config);
    
    
};

/*!
 * \class CTurbSASolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 2.0.6
 */

class CTurbSASolver: public CTurbSolver {
private:
	double nu_tilde_Inf;
	
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSASolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSASolver(void);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
    
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                          CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
    /*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
    /*!
	 * \brief Impose the nacelle inflow boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the ancelle exhaust boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
    /*!
	 * \brief Impose the interface boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
	 * \param[in] geometry - Geometric definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Calculate the laminar viscosity (used for AD).
	 * \param[in] val_U_i - Value of the flow variables at point i.
	 * \param[out] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CalcEddyViscosity(double *val_FlowVars, double val_laminar_viscosity,
                           double *val_TurbVar, double *val_eddy_viscosity);
    
    
    
};

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 2.0.6
 */

class CTransLMSolver: public CTurbSolver {
private:
	double Intermittency_Inf, REth_Inf;
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTransLMSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CTransLMSolver(void);
    
	/*!
	 * \brief Correlation function to relate turbulence intensity to transition onset reynolds number
	 * \param[in]  tu - turbulence intensity
	 * \param[out] REth - corresponding transition onset reynolds number
	 */
	double REthCorrelation(double tu);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
    
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                          CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	// Another set of matrix structures for the Lm equations
	CSysMatrix JacobianItmc; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
	double *LinSysSolItmc;		/*!< \brief vector to store iterative solution of implicit linear system. */
	double *LinSysResItmc;		/*!< \brief vector to store iterative residual of implicit linear system. */
	double *rhsItmc;		/*!< \brief right hand side of implicit linear system. */
	CSysMatrix JacobianReth; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
	double *LinSysSolReth;		/*!< \brief vector to store iterative solution of implicit linear system. */
	double *LinSysResReth;		/*!< \brief vector to store iterative residual of implicit linear system. */
	double *rhsReth;		/*!< \brief right hand side of implicit linear system. */
};

/*!
 * \class CTurbSSTSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Campos, F. Palacios, T. Economon
 * \version 2.0.6
 */

class CTurbSSTSolver: public CTurbSolver {
private:
	double *constants,  /*!< \brief Constants for the model. */
	kine_Inf,           /*!< \brief Free-stream turbulent kinetic energy. */
	omega_Inf;          /*!< \brief Free-stream specific dissipation. */
    
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSSTSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSSTSolver(void);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Computes the eddy viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                         CNumerics *numerics, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                          CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
    /*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
	 * \param[in] geometry - Geometric definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Get the constants for the SST model.
	 * \return A pointer to an array containing a set of constants
	 */
	double* GetConstants();
    
};

/*!
 * \class CAdjEulerSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 2.0.6
 */
class CAdjEulerSolver : public CSolver {
protected:
	double PsiRho_Inf,	/*!< \brief PsiRho variable at the infinity. */
	PsiE_Inf,			/*!< \brief PsiE variable at the infinity. */
	*Phi_Inf;			/*!< \brief Phi vector at the infinity. */
	double *Sens_Mach, /*!< \brief Mach sensitivity coefficient for each boundary. */
	*Sens_AoA,			/*!< \brief Angle of attack sensitivity coefficient for each boundary. */
	*Sens_Geo,			/*!< \brief Shape sensitivity coefficient for each boundary. */
	*Sens_Press,			/*!< \brief Pressure sensitivity coefficient for each boundary. */
	*Sens_Temp,			/*!< \brief Temperature sensitivity coefficient for each boundary. */
	**CSensitivity;		/*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
	double Total_Sens_Mach;	/*!< \brief Total mach sensitivity coefficient for all the boundaries. */
	double Total_Sens_AoA;		/*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
	double Total_Sens_Geo;		/*!< \brief Total shape sensitivity coefficient for all the boundaries. */
	double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
	double Total_Sens_Temp;    /*!< \brief Total farfield sensitivity to temperature. */
	double *p1_Und_Lapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*p2_Und_Lapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	bool space_centered;  /*!< \brief True if space centered scheeme used. */
    double **Jacobian_Axisymmetric; /*!< \brief Storage for axisymmetric Jacobian. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjEulerSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CAdjEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjEulerSolver(void);
    
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Compute the Fan face Mach number.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solution - Container vector with all the solutions.
	 */
	void GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Created the force projection vector for adjoint boundary conditions.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Compute the jump for the interior boundary problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	 * \brief Compute the spatial integration using a centered scheme for the adjoint equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                           unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term from objective function integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	//void SourceObjFunc_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
	//	CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Compute the undivided laplacian for the adjoint solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Parallelization of Undivided Laplacian.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the sensor for higher order dissipation control in rotating problems.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetDissipation_Switch(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the adjoint Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the interface adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                               unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the near-field adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                               unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the adjoint symmetry boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Impose the nacelle inflow adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the nacelle exhaust boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose a dirchlet boundary condition to set the couple aeroacoustic solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_FWH(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                unsigned short val_marker);
    
	/*!
	 * \brief Set the coupling terms for the adjoint aeroacoustic problem.
	 * \param[in] wave_solution - Solution container from the wave problem.
	 * \param[in] flow_solution - Solution container from the flow problem.
	 * \param[in] flow_geometry - Geometrical definition for the flow problem.
	 * \param[in] flow_config - Definition of the particular problem.
	 */
	void SetAeroacoustic_Coupling(CSolver ***wave_solution,  CSolver ***flow_solution, CNumerics *numerics,
                                  CGeometry **flow_geometry, CConfig *flow_config);
    
	/*!
	 * \brief Update the solution using a Runge-Kutta strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep);
    
	/*!
	 * \brief Update the solution using a explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Initialize the residual vectors.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the inviscid sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Smooth the inviscid sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
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
	double GetTotal_Sens_Geo(void);
    
	/*!
	 * \brief Set the total Mach number sensitivity coefficient.
	 * \return Value of the Mach sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_Mach(void);
    
	/*!
	 * \brief Set the total angle of attack sensitivity coefficient.
	 * \return Value of the angle of attack sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_AoA(void);
    
	/*!
	 * \brief Set the total farfield pressure sensitivity coefficient.
	 * \return Value of the farfield pressure sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_Press(void);
    
	/*!
	 * \brief Set the total farfield temperature sensitivity coefficient.
	 * \return Value of the farfield temperature sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_Temp(void);
    
    
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Get the value of kappapsi_Volume
	 * \return kappapsi_Volume
	 */
	double GetKappaPsiVolume(void);
    
    /*!
	 * \brief Set the initial condition for the Euler Equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] ExtIter - External iteration.
	 */
	void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
    
};

/*!
 * \class CAdjNSSolver
 * \brief Main class for defining the Navier-Stokes' adjoint flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 2.0.6
 */
class CAdjNSSolver : public CAdjEulerSolver {
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjNSSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CAdjNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CAdjNSSolver(void);
    
	/*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition (heat flux).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
    /*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition (heat flux).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the viscous sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Compute the viscous residuals for the adjoint equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                          CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
};

/*!
 * \class CAdjTurbSolver
 * \brief Main class for defining the adjoint turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 2.0.6
 */
class CAdjTurbSolver : public CSolver {
private:
	double PsiNu_Inf,	/*!< \brief PsiNu variable at the infinity. */
	*FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j;	/*!< \brief Store the flow solution at point j. */
    
	double *Jacobian_mui,
	*Jacobian_muj;
    
	double ***Jacobian_gradi,
	***Jacobian_gradj;
    
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
    double **Jacobian_k;			    /*!< \brief Auxiliary matrices for storing point to point Jacobians at point k. */
	double  **Jacobian_ik,			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_jk;			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
    
public:
    
	/*!
	 * \brief Default constructor of the class.
	 */
	CAdjTurbSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjTurbSolver(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Default destructor of the class.
	 */
	virtual ~CAdjTurbSolver(void);
    
	/*!
	 * \brief Impose the Navier-Stokes turbulent adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Initializate the residual vectors.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Compute the viscous residuals for the turbulent adjoint equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                          unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Conservative source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SourceConserv_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
};

/*!
 * \class CLinEulerSolver
 * \brief Main class for defining the linearized Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 2.0.6
 */
class CLinEulerSolver : public CSolver {
private:
	double DeltaRho_Inf,	/*!< \brief Linearized density variable at the infinity. */
	DeltaE_Inf,				/*!< \brief Linearized energy at the infinity. */
	*DeltaVel_Inf;			/*!< \brief Linearized velocity vector at the infinity. */
	double *p1_Und_Lapl,	/*!< \brief Undivided Laplacians for centered scheme. */
	*p2_Und_Lapl;			/*!< \brief Undivided Laplacians for centered scheme. */
	double *CDeltaDrag_Inv, /*!< \brief Linearized drag coefficient (inviscid contribution) for each boundary. */
	*CDeltaLift_Inv,		/*!< \brief Linearized lift coefficient (inviscid contribution) for each boundary. */
	*DeltaForceInviscid;	/*!< \brief Linearized inviscid force for each boundary. */
	double AllBound_CDeltaDrag_Inv, /*!< \brief Total linearized drag coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CDeltaLift_Inv;		/*!< \brief Total linearized lift coefficient (inviscid contribution) for all the boundaries. */
	double Total_CDeltaDrag,	/*!< \brief Total linearized drag coefficient for all the boundaries. */
	Total_CDeltaLift;			/*!< \brief Total linearized lift coefficient for all the boundaries. */
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CLinEulerSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CLinEulerSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CLinEulerSolver(void);
    
	/*!
	 * \brief Compute the spatial integration using a centered scheme for the linearized equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                           unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Compute the undivided laplacian for the linearized solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the linearized Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep);
    
	/*!
	 * \brief Restart residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the linearization of the pressure forces and all the adimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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

/*! \class CElectricSolver
 *  \brief Main class for defining the electric potential solver.
 *  \author F. Palacios.
 *  \version 2.0.6
 *  \date May 3, 2010.
 */
class CElectricSolver : public CSolver {
private:
	double Total_CCharge;			/*!< \brief Total charge coefficient for all the domain. */
	double *Source_Vector;		/*!< \brief Auxiliary vector for storing element source vector. */
    
    double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
    double **StiffMatrix_Elem,			/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	**StiffMatrix_Node;							/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
    
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CElectricSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CElectricSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] solver1_geometry - Geometrical definition of the problem.
	 * \param[in] solver1_solution - Container vector with all the solutions.
	 * \param[in] solver1_config - Definition of the particular problem.
	 * \param[in] solver2_geometry - Geometrical definition of the problem.
	 * \param[in] solver2_solution - Container vector with all the solutions.
	 * \param[in] solver2_config - Definition of the particular problem.
	 */
	void Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config, CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CElectricSolver(void);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] StiffMatrix_Elem - Element stiffness matrix
	 */
	void AddStiffMatrix(double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3);
    
	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh);
	/*!
	 * \brief Compute the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                          unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Definition of hte solver settings.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                   CConfig *config, unsigned short val_marker) ;
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker) ;
    
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
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

/*! \class CWaveSolver
 *  \brief Main class for defining the wave solver.
 *  \author F. Palacios.
 *  \version 2.0.6
 *  \date May 3, 2010.
 */
class CWaveSolver : public CSolver {
private:
	double *CWave;	/*!< \brief Wave strength for each boundary. */
	double AllBound_CWave;	/*!< \brief Total wave strength for all the boundaries. */
	double Total_CWave; /*!< \brief Total wave strength for all the boundaries. */
    
    CSysMatrix StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	CSysMatrix StiffMatrixTime;	/*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
    
    double **StiffMatrix_Elem,			/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	**StiffMatrix_Node;							/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CWaveSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CWaveSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CWaveSolver(void);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose a Neumann BC for the adjoint aeroacoustic problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Observer(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                     unsigned short val_marker) ;
    
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Set the noise sources from the flow problem for aeroacoustic computations.
	 * \param[in] wave_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] wave_config - Definition of the particular problem.
	 */
	void SetNoise_Source(CSolver ***flow_solution, CGeometry **wave_geometry, CConfig *wave_config);
    
	/*!
	 * \brief Compute the total wave strength coefficient.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Wave_Strength(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Load the direct solution from file for the adjoint problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void GetRestart(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Build stiffness matrix in space.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetSpace_Matrix(CGeometry *geometry,
                         CConfig   *config);
    
	/*!
	 * \brief Build stiffness matrix & Jacobian in time.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Matrix(CGeometry *geometry,
                        CConfig   *config);
	/*!
	 * \brief Provide the total wave strength.
	 * \return Value of the wave strength.
	 */
	double GetTotal_CWave(void);
    
};

/*! \class CHeatSolver
 *  \brief Main class for defining the heat solver.
 *  \author F. Palacios.
 *  \version 2.0.6
 *  \date May 3, 2010.
 */
class CHeatSolver : public CSolver {
private:
	double *CHeat;	/*!< \brief Heat strength for each boundary. */
	double AllBound_CHeat;	/*!< \brief Total Heat strength for all the boundaries. */
	double Total_CHeat; /*!< \brief Total Heat strength for all the boundaries. */
    
    CSysMatrix StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	CSysMatrix StiffMatrixTime;	/*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
    
    double **StiffMatrix_Elem,			/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	**StiffMatrix_Node;							/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CHeatSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CHeatSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CHeatSolver(void);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Load the direct solution from file for the adjoint problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void GetRestart(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Build stiffness matrix & Jacobian in time.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void SetTime_Matrix(CGeometry *geometry,
                        CConfig   *config);
	/*!
	 * \brief Provide the total heat strength.
	 * \return Value of the heat strength.
	 */
	double GetTotal_CHeat(void);
    
};

/*! \class CFEASolver
 *  \brief Main class for defining the FEA solver.
 *  \author F. Palacios.
 *  \version 2.0.6
 *  \date May 3, 2010.
 */
class CFEASolver : public CSolver {
private:
	double  Total_CFEA;			/*!< \brief Total FEA coefficient for all the boundaries. */
    
    CSysMatrix StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	CSysMatrix StiffMatrixTime;	/*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
    
    double **StiffMatrix_Elem,			/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	**StiffMatrix_Node;							/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CFEASolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CFEASolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CFEASolver(void);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Impose a displacement (constraint) boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);
    
	/*!
	 * \brief Impose a load boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_FlowLoad(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                     unsigned short val_marker);
    
	/*!
	 * \brief Impose a load boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                 unsigned short val_marker);
    
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Set the the pressure load in the FEA solver.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] fea_config - Definition of the particular problem.
	 */
	void SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config);
    
	/*!
	 * \brief Set the initial condition for the FEA Equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] ExtIter - External iteration.
	 */
	void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional FEA coefficient.
	 * \return Value of the FEA coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFEA(void);
    
	/*!
	 * \brief Set the value of the FEA coefficient.
	 * \param[in] val_cfea - Value of the FEA coefficient.
	 */
	void SetTotal_CFEA(double val_cfea);
    
};

/*!
 * \class CLevelSetSolver
 * \brief Main class for defining the level set solver.
 * \ingroup LevelSet_Model
 * \author F. Palacios.
 * \version 2.0.6
 */
class CLevelSetSolver : public CSolver {
protected:
	double *FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j,		/*!< \brief Store the flow solution at point j. */
	Total_CFreeSurface;			/*!< \brief Total Free Surface coefficient for all the boundaries. */
    
    double **Jacobian_MeanFlow_i,	/*!< \brief Auxiliary matrices for storing point to point Jacobians of the mean flow at point i. */
	**Jacobian_MeanFlow_j;			    /*!< \brief Auxiliary matrices for storing point to point Jacobians of the mean flow at point j. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CLevelSetSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CLevelSetSolver(void);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Free Surface coefficient.
	 * \return Value of the Free Surface coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFreeSurface(void);
    
	/*!
	 * \brief Set the value of the Free Surface coefficient.
	 * \param[in] val_cfreesurface - Value of the Free Surface coefficient.
	 */
	void SetTotal_CFreeSurface(double val_cfreesurface);
    
	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Recompute distance to the level set 0.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetLevelSet_Distance(CGeometry *geometry, CConfig *config, bool Initialization, bool WriteLevelSet);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
};

/*!
 * \class CAdjLevelSetSolver
 * \brief Main class for defining the level set solver.
 * \ingroup LevelSet_Model
 * \author F. Palacios.
 * \version 2.0.6
 */
class CAdjLevelSetSolver : public CSolver {
protected:
	double *FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j,		/*!< \brief Store the flow solution at point j. */
	Total_CFreeSurface;			/*!< \brief Total Free Surface coefficient for all the boundaries. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjLevelSetSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjLevelSetSolver(void);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
};

/*!
 * \class CTemplateSolver
 * \brief Main class for defining the template model solver.
 * \ingroup Template_Flow_Equation
 * \author F. Palacios.
 * \version 2.0.6
 */
class CTemplateSolver : public CSolver {
private:
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CTemplateSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTemplateSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CTemplateSolver(void);
    
	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the time step for solving the Euler equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] Iteration - Index of the current iteration.
	 */
	void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh, unsigned long Iteration);
    
	/*!
	 * \brief Compute the spatial integration using a centered scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using a linear solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                   CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using a Runge-Kutta scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep);
    
	/*!
	 * \brief Update the solution using the explicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
};

/*!
 * \class CPlasmaSolver
 * \brief Main class for defining the plasma solver.
 * \author ADL Stanford.
 * \version 2.0.6
 */
class CPlasmaSolver : public CSolver {
protected:
	bool restart,
	restart_from_Euler,
	implicit,
	centered_scheme,
	axisymmetric,
	weighted_LS;
    
	unsigned short nMarker; /*!< \brief Total number of markers using the grid information. */
	double Gas_Constant;		/*! \brief Gas constant for each species. */
	double Enthalpy_formation;		/*! \brief Chemical enthalpy of formation for all species [?? From Candler]. */
	double Molar_Mass; /*!< \brief Molar mass of each chemical constituent [kg/kmol].*/
	double *Gas_Composition; /*!< \brief Initial mass fraction of each gas constituent. */
	double **CSkinFriction;	/*!< \brief Skin friction coefficient for each boundary and vertex. */
	double ***CHeatTransfer;	/*!< \brief Heat Transfer coefficient for each boundary and vertex. */
	double **CPressure;		/*!< \brief Pressure coefficient for each boundary and vertex. */
    double ****CViscForce;  /*!< \brief Viscous Forces for each boundary and vertex. */
    double ****CPressForce; /*!< \brief Pressure Forces for each boundary and vertex. */
    
    
	double *Density_Inf;		/*!< \brief Density of all species at the infinity. */
	double *Temperature_Inf;		/*!< \brief Temperature of all species at the infinity. */
	double *Mach_Inf;			/*!< \brief Mach number of all fluids at the infinity. */
	double *Energy_Inf;			/*!< \brief Energy of all fluids at the infinity. */
	double *Energy_vib_Inf;			/*!< \brief Energy of all fluids at the infinity. */
	double Energy_el_Inf;		/*!< \brief Energy of electrons at the infinity. */
	double *Pressure_Inf;		/*!< \brief Pressure of all fluids at the infinity. */
	double **Velocity_Inf;		/*!< \brief Flow Velocity vector of all species at the infinity. */
	double *Velocity_Inf_mean;	/*!< \brief Flow Velocity vector of mass averaged fluid at the infinity. */
    
	double *Density_Inlet;		/*!< \brief Density of all species at the inlet. */
	double *Density_Outlet;		/*!< \brief Density of all species at the outlet. */
	double *Mach_Inlet;			/*!< \brief Mach number of all fluids at the inlet. */
	double *Mach_Outlet;		/*!< \brief Mach number of all fluids at the outlet. */
	double *Energy_Inlet;		/*!< \brief Energy of all fluids at the inlet. */
	double *Energy_vib_Inlet; /*!< \brief Vibrational energy of all fluids at the inlet. */
	double *Energy_Outlet;		/*!< \brief Energy of all fluids at the outlet. */
	double *Energy_vib_Outlet; /*!< \brief Vibrational energy of all fluids at the outlet. */
	double **Velocity_Inlet;	/*!< \brief Flow Velocity vector of all species at the inlet. */
	double **Velocity_Outlet;	/*!< \brief Flow Velocity vector of all species at the outlet. */
	double *Species_Delta;
	double *Mag_Force;
    
	double **p1_Und_Lapl;		/*!< \brief Auxiliary variable for the undivided Laplacians. */
	double **p2_Und_Lapl;		/*!< \brief Auxiliary variable for the undivided Laplacians. */
	double *PrimVar_i;			/*!< \brief Auxiliary vector for storing the solution at point i. */
	double *PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	double **Precon_Mat_inv; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
    double ***PrimGrad_i, ***PrimGrad_j, **PrimGradLimiter_i, **PrimGradLimiter_j;
    double ***PrimVar_max, ***PrimVar_min;
    
    double *Min_Delta_Time; /*!< \brief Minimum time step [nSpecies]. */
    double *Max_Delta_Time; /*!< \brief Maximum time step [nSpecies]. */
    
    
	double *Residual_Chemistry,	/*!< \brief Auxiliary vector for storing the residual from the chemistry source terms. */
	*Residual_ElecForce,	/*!< \brief Auxiliary vector for storing the electrostatic force source terms. */
	*Residual_MomentumExch,	/*!< \brief Auxiliary vector for storing the momentum exhange source terms. */
	*Residual_EnergyExch,	/*!< \brief Auxiliary vector for storing the energy exchange source terms. */
	*Residual_Axisymmetric; /*!< \brief Auxiliary vector for storing the axisymmetric source terms. */
    
	double **Jacobian_Chemistry,	/*!< \brief Auxiliary matrix for storing point to point jacobians for chemistry terms at pt i. */
	**Jacobian_ElecForce,					/*!< \brief Auxiliary matrix for storing point to point jacobians for electrostatic force terms at pt i. */
	**Jacobian_MomentumExch,			/*!< \brief Auxiliary matrix for storing point to point jacobians for momentum exchange terms at pt i. */
	**Jacobian_EnergyExch,				/*!< \brief Auxiliary matrix for storing point to point jacobians for energy exchange terms at pt i. */
	**Jacobian_Axisymmetric;			/*!< \brief Auxiliary matrix for storing axismmetric jacobians. */
    
	double *CDrag_Inv,	/*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
	*CLift_Inv,			/*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
	*CSideForce_Inv,		/*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
	*CMx_Inv,			/*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
	*CMy_Inv,			/*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
	*CMz_Inv,			/*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
	*CFx_Inv,			/*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
	*CFy_Inv,			/*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
	*CFz_Inv,			/*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
	*CEff_Inv,				/*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
	*CMerit_Inv,				/*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
	*CT_Inv,			/*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
	*CQ_Inv,			/*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
	*CEquivArea_Inv,				/*!< \brief Equivalent area (inviscid contribution) for each boundary. */
	*CNearFieldOF_Inv,				/*!< \brief Near field pressure (inviscid contribution) for each boundary. */
	*ForceInviscid,		/*!< \brief Inviscid force for each boundary. */
	*MomentInviscid;	/*!< \brief Inviscid moment for each boundary. */
	double *ForceViscous,	/*!< \brief Viscous force for each boundary. */
	*MomentViscous;			/*!< \brief Inviscid moment for each boundary. */
    
	double Total_CDrag, /*!< \brief Total drag coefficient for all the boundaries. */
	Total_CLift,		/*!< \brief Total lift coefficient for all the boundaries. */
	Total_CSideForce,		/*!< \brief Total sideforce coefficient for all the boundaries. */
	Total_CMx,			/*!< \brief Total x moment coefficient for all the boundaries. */
	Total_CMy,			/*!< \brief Total y moment coefficient for all the boundaries. */
	Total_CMz,			/*!< \brief Total z moment coefficient for all the boundaries. */
	Total_CFx,			/*!< \brief Total x force coefficient for all the boundaries. */
	Total_CFy,			/*!< \brief Total y force coefficient for all the boundaries. */
	Total_CFz,			/*!< \brief Total z force coefficient for all the boundaries. */
	Total_CEff,			/*!< \brief Total efficiency coefficient for all the boundaries. */
	Total_CMerit,			/*!< \brief Total rotor Figure of Merit for all the boundaries. */
	Total_CT,		/*!< \brief Total thrust coefficient for all the boundaries. */
	Total_CQ,		/*!< \brief Total torque coefficient for all the boundaries. */
    Total_Q,    /*!< \brief Total heat load for all the boundaries. */
    Total_Maxq, /*!< \brief Maximum heat flux on all boundaries. */
	Total_CEquivArea,			/*!< \brief Total Equivalent Area coefficient for all the boundaries. */
	Total_CNearFieldOF;			/*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
    
    double PressureDrag, /*!< \brief Total Pressure Drag coefficient for all the boundaries. */
    ViscDrag,           /*!< \brief Total Visc Drag coefficient for all the boundaries. */
    MagnetDrag;         /*!< \brief Total Magnet Drag coefficient for all the boundaries. */
    
	double AllBound_CDrag_Inv,	/*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CLift_Inv,			/*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CSideForce_Inv,			/*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMx_Inv,			/*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMy_Inv,			/*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMz_Inv,			/*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFx_Inv,			/*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFy_Inv,			/*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFz_Inv,			/*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEff_Inv,			/*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMerit_Inv,			/*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
	AllBound_CT_Inv,			/*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CQ_Inv,			/*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEquivArea_Inv,			/*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CNearFieldOF_Inv;			/*!< \brief Near-Field press coefficient (inviscid contribution) for all the boundaries. */
    
	double AllBound_CDrag_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
	AllBound_CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
	AllBound_CMx_Visc,			/*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMy_Visc,			/*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMz_Visc,			/*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CEff_Visc,			/*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
	AllBound_CFx_Visc,			/*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFy_Visc,			/*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CFz_Visc,			/*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
	AllBound_CMerit_Visc,			/*!< \brief Rotor Figure of Merit coefficient (Viscous contribution) for all the boundaries. */
	AllBound_CT_Visc,		/*!< \brief Thrust coefficient (viscous contribution) for all the boundaries. */
	AllBound_Q_Visc;		/*!< \brief Heat load (integrated heat flux) for all the boundaries. */
    
	double *CDrag_Visc,	/*!< \brief Drag coefficient (viscous contribution) for each boundary. */
	*CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for each boundary. */
	*CMx_Visc,			/*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
	*CMy_Visc,			/*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
	*CMz_Visc,			/*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
	*CFx_Visc,			/*!< \brief Force x coefficient (viscous contribution) for each boundary. */
	*CFy_Visc,			/*!< \brief Force y coefficient (viscous contribution) for each boundary. */
	*CFz_Visc,			/*!< \brief Force z coefficient (viscous contribution) for each boundary. */
	*CEff_Visc,			/*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
	*CMerit_Visc,			/*!< \brief Rotor Figure of Merit (Viscous contribution) for each boundary. */
	*CT_Visc,		/*!< \brief Thrust coefficient (viscous contribution) for each boundary. */
	*CQ_Visc,		/*!< \brief Torque coefficient (viscous contribution) for each boundary. */
    *Q_Visc, /*!< \brief Heat load (viscous contribution) for each boundary. */
    *Maxq_Visc; /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
    
    
	double Prandtl_Lam;   	/*!< \brief Laminar Prandtl number. */
	bool roe_turkel;         /*!< \brief True if computing preconditioning matrix for roe-turkel method. */
	bool magnet;         /*!< \brief True if including magnetic field in simulation */
    
    double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
    
	unsigned short nSpecies,				/*! \brief Number of species in the plasma simulation. */
	nMonatomics,										/*! \brief Number of monatomic species in the flow. */
	nDiatomics;											/*! \brief Number of diatomic species in the flow. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CPlasmaSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPlasmaSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CPlasmaSolver(void);
    
    void SetVel_Residual_Zero(unsigned long val_ipoint, unsigned short iSpecies);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnSpecies(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnMonatomics(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnDiatomics(void);
    
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
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the skin friction coefficient.
	 */
	double GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex);
    
	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_iSpecies, unsigned short val_vertex);
    
	/*!
	 * \brief Compute Viscous Forces
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_iSpecies - index of the chemical species
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetViscForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex);
    
	/*!
	 * \brief Compute Pressure Forces
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_iSpecies - index of the chemical species
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetPressureForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex);
    
	/*!
	 * \brief Compute the time step for solving the Euler equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] Iteration - Index of the current iteration.
	 */
	void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh, unsigned long Iteration);
    
    /*!
     * \brief Compute the time step for solving the Euler equations.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
     * \param[in] Iteration - Index of the current iteration.
     */
    void SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the spatial integration using a centered scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
    
	/*!
	 * \brief Compute the viscous residuals.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                          CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Get the value of the maximum delta time.
     * \param[in] val_Species - Value of the species
	 * \return Value of the maximum delta time for val_Species.
	 */
	double GetMax_Delta_Time(unsigned short val_Species);
    
    /*!
	 * \brief Get the value of the maximum delta time.
     * \param[in] val_Species - Value of the species
	 * \return Value of the minimum delta time for val_Species.
	 */
    double GetMin_Delta_Time(unsigned short val_Species);
    
	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute a pressure sensor switch.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetDissipation_Switch(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the gradient of the primitive variables using Green-Gauss method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Compute the gradient of the primitive variables using Green-Gauss method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
    void SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config, unsigned long iVertex, unsigned short val_marker, double *val_PrimVar_i);
    
	/*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config, unsigned long iPoint, double *val_PrimVar_i);
    
    /*!
     * \brief Send the gradient of the primitive variables to the other processors
     */
    void Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimVar_Limiter_MPI(CGeometry *geometry, CConfig *config);
    
    void SetPrimVar_Limiter(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
	 * \param[in] iPoint - Index of the grid point
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPreconditioner(CConfig *config, unsigned short iPoint);
    
	/*!
	 * \brief Compute the undivided laplacian for the solution, except the energy equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Definition of hte solver settings.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Compute the inviscid forces and all the addimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_Forces(CGeometry *geometry, CConfig *config);
    
    
	/*!
	 * \brief Compute the viscous forces and all the addimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Viscous_Forces(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose an isothermal wall boundary condition (no-slip).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief Impose a constant heat flux wall boundary condition (no-slip).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                    CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the subsonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                      CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the subsonic outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
     
	 */
	void BC_Dielectric(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                       CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                   CConfig *config, unsigned short val_marker) ;
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker) ;
    
	/*!
	 * \brief Update the solution using an implicit Euler scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Compute the pressure forces and all the adimensional coefficients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Electrostatic_Forces(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
	 * \return Value of the lift coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CLift(void);
    
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
	 * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
	 * \return Value of the force x coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFx(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFy(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
	 * \return Value of the force z coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CFz(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
	 * \return Value of the sideforce coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CSideForce(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
	 * \return Value of the efficiency coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CEff(void);
    
    /*!
	 * \brief Provide the total integrated heat flux.
	 * \return Value of the integrated surface heat flux.
	 */
	double GetTotal_Q(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	double Get_PressureDrag(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	double Get_ViscDrag(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	double Get_MagnetDrag(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] solver1_geometry - Geometrical definition of the problem.
	 * \param[in] solver1_solution - Container vector with all the solutions.
	 * \param[in] solver1_config - Definition of the particular problem.
	 * \param[in] solver2_geometry - Geometrical definition of the problem.
	 * \param[in] solver2_solution - Container vector with all the solutions.
	 * \param[in] solver2_config - Definition of the particular problem.
	 */
	void Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config, CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config);
    
};

/*!
 * \class CAdjPlasmaSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 2.0.6
 */
class CAdjPlasmaSolver : public CSolver {
protected:
	double *Residual_Chemistry,	/*!< \brief Auxiliary vector for storing the residual from the chemistry source terms. */
	*Residual_ElecForce,	/*!< \brief Auxiliary vector for storing the electrostatic force source terms. */
	*Residual_MomentumExch,	/*!< \brief Auxiliary vector for storing the momentum exhange source terms. */
	*Residual_EnergyExch,	/*!< \brief Auxiliary vector for storing the energy exchange source terms. */
	*Residual_Axisymmetric; /*!< \brief Auxiliary vector for storing the axisymmetric source terms. */
    
	double **Jacobian_Chemistry,	/*!< \brief Auxiliary matrix for storing point to point jacobians for chemistry terms at pt i. */
	**Jacobian_ElecForce,					/*!< \brief Auxiliary matrix for storing point to point jacobians for electrostatic force terms at pt i. */
	**Jacobian_MomentumExch,			/*!< \brief Auxiliary matrix for storing point to point jacobians for momentum exchange terms at pt i. */
	**Jacobian_EnergyExch,				/*!< \brief Auxiliary matrix for storing point to point jacobians for energy exchange terms at pt i. */
	**Jacobian_Axisymmetric;			/*!< \brief Auxiliary matrix for storing axismmetric jacobians. */
    
	double PsiRho_Inf,	/*!< \brief PsiRho variable at the infinity. */
	PsiE_Inf,			/*!< \brief PsiE variable at the infinity. */
	PsiEvib_Inf,    /*1< \brief PsiEvib variable at the infinity. */
	*Phi_Inf;			/*!< \brief Phi vector at the infinity. */
	double *Sens_Mach, /*!< \brief Mach sensitivity coefficient for each boundary. */
	*Sens_AoA,			/*!< \brief Angle of attack sensitivity coefficient for each boundary. */
	*Sens_Geo,			/*!< \brief Shape sensitivity coefficient for each boundary. */
	*Sens_Press,			/*!< \brief Pressure sensitivity coefficient for each boundary. */
	*Sens_Temp,			/*!< \brief Temperature sensitivity coefficient for each boundary. */
	**CSensitivity;		/*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
	double Total_Sens_Mach;	/*!< \brief Total mach sensitivity coefficient for all the boundaries. */
	double Total_Sens_AoA;		/*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
	double Total_Sens_Geo;		/*!< \brief Total shape sensitivity coefficient for all the boundaries. */
	double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
	double Total_Sens_Temp;    /*!< \brief Total farfield sensitivity to temperature. */
	double *p1_Und_Lapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*p2_Und_Lapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	bool space_centered;  /*!< \brief True if space centered scheeme used. */
    
	unsigned short nSpecies,				/*! \brief Number of species in the plasma simulation. */
	nMonatomics,										/*! \brief Number of monatomic species in the flow. */
	nDiatomics;											/*! \brief Number of diatomic species in the flow. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjPlasmaSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAdjPlasmaSolver(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjPlasmaSolver(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnSpecies(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnMonatomics(void);
    
	/*!
	 * \brief Get the number of Species present in the flow.
	 */
	unsigned short GetnDiatomics(void);
    
    /*!
     * \overload
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the problem.
     */
    void Set_MPI_Solution(CGeometry *geometry, CConfig *config);
    
    /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Created the force projection vector for adjoint boundary conditions.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
                         unsigned short iMesh);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the adjoint Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the adjoint symmetry boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
                      unsigned short val_marker);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
	/*!
	 * \brief Initialize the residual vectors.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem);
    
	/*!
	 * \brief Compute the inviscid sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
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
	double GetTotal_Sens_Geo(void);
    
	/*!
	 * \brief Set the total Mach number sensitivity coefficient.
	 * \return Value of the Mach sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_Mach(void);
    
	/*!
	 * \brief Set the total angle of attack sensitivity coefficient.
	 * \return Value of the angle of attack sensitivity coefficient 
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_AoA(void);
    
	/*!
	 * \brief Set the total farfield pressure sensitivity coefficient.
	 * \return Value of the farfield pressure sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_Press(void);
    
	/*!
	 * \brief Set the total farfield temperature sensitivity coefficient.
	 * \return Value of the farfield temperature sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	double GetTotal_Sens_Temp(void);
};


#include "solver_structure.inl"
