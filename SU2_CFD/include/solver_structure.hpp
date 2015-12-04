/*!
 * \file solver_structure.hpp
 * \brief Headers of the main subroutines for solving partial differential equations.
 *        The subroutines and functions are in the <i>solver_structure.cpp</i>,
 *        <i>solution_direct.cpp</i>, <i>solution_adjoint.cpp</i>, and
 *        <i>solution_linearized.cpp</i> files.
 * \author F. Palacios, T. Economon
 * \version 4.0.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
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

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>

#include "fluid_model.hpp"
#include "numerics_structure.hpp"
#include "variable_structure.hpp"
#include "../../Common/include/gauss_structure.hpp"
#include "../../Common/include/element_structure.hpp"
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
 * a child class for each particular solver (Euler, Navier-Stokes, etc.)
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
 */
class CSolver {
protected:
	unsigned short IterLinSolver;	/*!< \brief Linear solver iterations. */
	unsigned short nVar,					/*!< \brief Number of variables of the problem. */
  nPrimVar,                     /*!< \brief Number of primitive variables of the problem. */
  nPrimVarGrad,                 /*!< \brief Number of primitive variables of the problem in the gradient computation. */
  nSecondaryVar,                     /*!< \brief Number of primitive variables of the problem. */
  nSecondaryVarGrad,                 /*!< \brief Number of primitive variables of the problem in the gradient computation. */
	nDim;													/*!< \brief Number of dimensions of the problem. */
	unsigned long nPoint;					/*!< \brief Number of points of the computational grid. */
  unsigned long nPointDomain; 	/*!< \brief Number of points of the computational grid. */
	su2double Max_Delta_Time,	/*!< \brief Maximum value of the delta time for all the control volumes. */
	Min_Delta_Time;					/*!< \brief Minimum value of the delta time for all the control volumes. */
	su2double *Residual_RMS,	/*!< \brief Vector with the mean residual for each variable. */
  *Residual_Max,        /*!< \brief Vector with the maximal residual for each variable. */
	*Residual,						/*!< \brief Auxiliary nVar vector. */
	*Residual_i,					/*!< \brief Auxiliary nVar vector for storing the residual at point i. */
	*Residual_j;					/*!< \brief Auxiliary nVar vector for storing the residual at point j. */
  unsigned long *Point_Max; /*!< \brief Vector with the maximal residual for each variable. */
  su2double **Point_Max_Coord; /*!< \brief Vector with pointers to the coords of the maximal residual for each variable. */
	su2double *Solution,		/*!< \brief Auxiliary nVar vector. */
	*Solution_i,				/*!< \brief Auxiliary nVar vector for storing the solution at point i. */
	*Solution_j;				/*!< \brief Auxiliary nVar vector for storing the solution at point j. */
	su2double *Vector,	/*!< \brief Auxiliary nDim vector. */
	*Vector_i,			/*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point i. */
	*Vector_j;			/*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point j. */
	su2double *Res_Conv,	/*!< \brief Auxiliary nVar vector for storing the convective residual. */
	*Res_Visc,				/*!< \brief Auxiliary nVar vector for storing the viscous residual. */
	*Res_Sour,				/*!< \brief Auxiliary nVar vector for storing the viscous residual. */
	*Res_Conv_i,		  /*!< \brief Auxiliary vector for storing the convective residual at point i. */
	*Res_Visc_i,			/*!< \brief Auxiliary vector for storing the viscous residual at point i. */
	*Res_Conv_j,			/*!< \brief Auxiliary vector for storing the convective residual at point j. */
	*Res_Visc_j;			/*!< \brief Auxiliary vector for storing the viscous residual at point j. */
	su2double **Jacobian_i,	/*!< \brief Auxiliary matrices for storing point to point Jacobians at point i. */
	**Jacobian_j;			    /*!< \brief Auxiliary matrices for storing point to point Jacobians at point j. */
	su2double **Jacobian_ii,	/*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_ij,			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_ji,			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
	**Jacobian_jj;			  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  
	su2double **Smatrix,	/*!< \brief Auxiliary structure for computing gradients by least-squares */
	**cvector;			 /*!< \brief Auxiliary structure for computing gradients by least-squares */

    unsigned short nOutputVariables;  /*!< \brief Number of variables to write. */

public:
  
  CSysVector LinSysSol;		/*!< \brief vector to store iterative solution of implicit linear system. */
  CSysVector LinSysRes;		/*!< \brief vector to store iterative residual of implicit linear system. */
  CSysVector LinSysAux;		/*!< \brief vector to store iterative residual of implicit linear system. */
  CSysMatrix Jacobian; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
  
  CSysMatrix StiffMatrix; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations, and grid movement. */

  CSysVector OutputVariables;		/*!< \brief vector to store the extra variables to be written. */
  string* OutputHeadingNames; /*< \brief vector of strings to store the headings for the exra variables */
  
  CVariable** node;	/*!< \brief Vector which the define the variables for each problem. */
  CVariable* node_infty; /*!< \brief CVariable storing the free stream conditions. */
  
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
	 * \brief Set number of linear solver iterations.
	 * \param[in] val_iterlinsolver - Number of linear iterations.
	 */
	virtual void Set_MPI_Primitive(CGeometry *geometry, CConfig *config);
  
//  /*!
//	 * \brief Set number of linear solver iterations.
//	 * \param[in] val_iterlinsolver - Number of linear iterations.
//	 */
//	virtual void Set_MPI_Secondary(CGeometry *geometry, CConfig *config);

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
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config);
 
//  /*!
//	 * \brief Impose the send-receive boundary condition.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//  virtual void Set_MPI_Secondary_Limiter(CGeometry *geometry, CConfig *config);

  /*!
	 * \brief Set the fluid solver nondimensionalization.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Compute the pressure at the infinity.
	 * \return Value of the pressure at the infinity.
	 */
  virtual CFluidModel* GetFluidModel(void);

  	/*!
	 * \brief Get number of linear solver iterations.
	 * \return Number of linear solver iterations.
	 */
	unsigned short GetIterLinSolver(void);
    
	/*!
	 * \brief Get the value of the maximum delta time.
	 * \return Value of the maximum delta time.
	 */
	su2double GetMax_Delta_Time(void);
    
	/*!
	 * \brief Get the value of the minimum delta time.
	 * \return Value of the minimum delta time.
	 */
	su2double GetMin_Delta_Time(void);
    
    /*!
	 * \brief Get the value of the maximum delta time.
	 * \return Value of the maximum delta time.
	 */
	virtual su2double GetMax_Delta_Time(unsigned short val_Species);
    
	/*!
	 * \brief Get the value of the minimum delta time.
	 * \return Value of the minimum delta time.
	 */
	virtual su2double GetMin_Delta_Time(unsigned short val_Species);
    
	/*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnVar(void);
  
  /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnPrimVar(void);
  
  /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnPrimVarGrad(void);
  
  /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnSecondaryVar(void);
  
  /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnSecondaryVarGrad(void);
  
  /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnOutputVariables(void);
    
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
	void SetRes_RMS(unsigned short val_var, su2double val_residual);
    
	/*!
	 * \brief Adds the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void AddRes_RMS(unsigned short val_var, su2double val_residual);
    
	/*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
	 */
	su2double GetRes_RMS(unsigned short val_var);
    
    /*!
	 * \brief Set the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
	 */
	void SetRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point);
    
	/*!
	 * \brief Adds the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   * \param[in] val_point - Value of the point index for the max residual.
   * \param[in] val_coord - Location (x, y, z) of the max residual point.
	 */
	void AddRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point, su2double* val_coord);
    
	/*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
	 */
	su2double GetRes_Max(unsigned short val_var);
    
    /*!
	 * \brief Get the maximal residual, this is useful for the convergence history.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
	 */
	unsigned long GetPoint_Max(unsigned short val_var);
  
  /*!
   * \brief Get the location of the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the location (x, y, z) of the biggest residual for the variable <i>val_var</i>.
   */
  su2double* GetPoint_Max_Coord(unsigned short val_var);
  
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
	void SetGridVel_Gradient(CGeometry *geometry, CConfig *config);
    
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
	virtual void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config);
  
//	/*!
//	 * \brief A virtual member.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	virtual void SetSecondary_Limiter(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Compute the pressure laplacian using in a incompressible solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] PressureLaplacian - Pressure laplacian.
	 */
	void SetPressureLaplacian(CGeometry *geometry, su2double *PressureLaplacian);
    
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
	 * \brief A virtual member, overloaded.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 *
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics,
                                unsigned short iMesh);
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
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
	virtual void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief A virtual member overloaded.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Container vector of the numerics of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	virtual void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output);

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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                               unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */


	virtual void BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                 unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */


	virtual void BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                 unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */


	virtual void BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                 unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	virtual void BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);

  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */

	virtual void BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
						 unsigned short val_marker);

  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */


	virtual void BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                              unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_ActDisk_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Isothermal_Wall(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *conv_numerics,
                                  CNumerics *visc_numerics,
                                  CConfig *config,
                                  unsigned short val_marker);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker);
      
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the Euler boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_NonReflecting(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
		
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Engine_Bleed(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Dielec(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                              CConfig *config, unsigned short val_marker);
	/*!
	 * \brief It performs the average value along a boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the average is evaluated.
	 */
	virtual void Mixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker);

	/*!
	 * \brief It performs the average value along a boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the average is evaluated.
	 */
	virtual void MPIMixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short marker_flag);

	/*!
	 * \brief it performs a mixed out average of the nodes of a boundary.
	 * \param[in] val_init_pressure -  initial pressure value
	 * \param[in] val_Averaged_Flux - flux averaged values.
     * \param[in] val_normal - normal vector.
     * \param[in] pressure_mix - value of the mixed-out avaraged pressure.
	 * \param[in] density_miz - value of the mixed-out avaraged density.
	 */
	virtual void MixedOut_Average (su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *pressure_mix, su2double *density_mix);

	/*!
	 * \brief it finds the root of an implicit equation that relates pressure and density.
	 * \param[in] pressure - pressure value
	 * \param[in] val_Averaged_Flux - flux averaged values.
     * \param[in] val_normal - normal vector.
     * \param[in] valfunc - Description of the numerical method.
	 * \param[in] density - value of the mixed-out avaraged density.
	 */
	virtual void MixedOut_Root_Function(su2double *pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *valfunc, su2double *density);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in]  c4k - Fourier transformation coefficients.
	 * \param[in]  nboundaryvertex - pithcwise ordered vertex.
	 */
	virtual void Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> > &c4k,signed long &nboundaryvertex);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in]  c2k - Fourier transformation coefficients.
	 * \param[in]  c3k - Fourier transformation coefficients.
	 * \param[in]  nboundaryvertex - pithcwise ordered vertex.
	 */
	virtual void Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> >& c2k,vector<std::complex<su2double> >& c3k,signed long& nboundaryvertex);

	/*!
	 * \brief A virtual member.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] intMarker - internal marker.
	 * \param[in] extMarker - external marker.
	 */
	 virtual void SetExtAveragedValue(CSolver *solver_container, unsigned short intMarker,  unsigned short extMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Density on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedDensity(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Pressure on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedPressure(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Enthalpy on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedEnthalpy(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Velocity on the surface <i>val_marker</i>.
	  */
	 virtual su2double* GetAveragedVelocity(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Grid Velocity on the surface <i>val_marker</i>.
	  */
	 virtual su2double* GetAveragedGridVelocity(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Entropy on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedEntropy(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Total Temperature on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedTotTemperature(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedTotPressure(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the flow angle on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetFlowAngle(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Mach Number on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedMach(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Normal Mach Number on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedNormalMach(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Mass flow on the surface <i>val_marker</i>..
	  */
	 virtual su2double GetMassFlow(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of Total Pressure Loss for turbomachinery performance.
	  */
	 virtual su2double GetTotalPressureLoss(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Kinetic energy loss for turbomachinery performance.
	  */
	 virtual su2double GetKineticEnergyLoss(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Total-total efficiency for turbomachinery performance.
	  */
	 virtual su2double GetTotalTotalEfficiency(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Total-static efficiency for turbomachinery performance.
	  */
	 virtual su2double GetTotalStaticEfficiency(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Eulerian Work for turbomachinery performance.
	  */
	 virtual su2double GetEulerianWork(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Inlet Total Enthalpy for turbomachinery performance.
	  */
	 virtual su2double GetTotalEnthalpyIn(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Inlet Flow Angle for turbomachinery performance.
	  */
	 virtual su2double GetFlowAngleIn(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Outlet Flow Angle for turbomachinery performance.
	  */
	 virtual su2double GetFlowAngleOut(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Inlet Mass Flow for turbomachinery performance.
	  */
	 virtual su2double GetMassFlowIn(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Outlet Mass FlowS for turbomachinery performance.
	  */
	 virtual su2double GetMassFlowOut(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Inlet Mach for turbomachinery performance.
	  */
	 virtual su2double GetMachIn(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Outlet Mach for turbomachinery performance.
	  */
	 virtual su2double GetMachOut(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the normal component of the Inlet Mach for turbomachinery performance.
	  */
	 virtual su2double GetNormalMachIn(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the normal component of the Outlet Mach for turbomachinery performance.
	  */
	 virtual su2double GetNormalMachOut(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Outlet Static Enthalpy for turbomachinery performance.
	  */
	 virtual su2double GetEnthalpyOut(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Outlet Isentropic Velocity for turbomachinery performance.
	  */
	 virtual su2double GetVelocityOutIs(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Outlet Pressure for turbomachinery performance.
	  */
	 virtual su2double GetPressureOut(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Pressure ratio for turbomachinery performance.
	  */
	 virtual su2double GetPressureRatio(unsigned short inMarkerTP);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Normal Velocity on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedNormalVelocity(unsigned short valMarker);

	 /*!
	  * \brief A virtual member.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Tangent Velocity on the surface <i>val_marker</i>.
	  */
	 virtual su2double GetAveragedTangVelocity(unsigned short valMarker);

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
	 */
	virtual void ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	 * \param[in] solver - solver containing the outlet information.
	 * \param[in] inMarker - marker related to the inlet.
	 * \param[in] outMarker - marker related to the outlet.
	 */
	virtual void TurboPerformance(CSolver *solver,  CConfig *config, unsigned short inMarker,  unsigned short outMarker, unsigned short Kind_TurboPerf, unsigned short inMarkerTP);

	/*!
	 * \brief A virtual member.
	 * \param[in] config - contains config file information.
	 */
	virtual void MPITurboPerformance(CConfig *config);
	/*!
	 * \brief A virtual member.
	 * \param[in] solver - solver containing the outlet information.
	 * \param[in] inMarker - marker related to the inlet.
	 * \param[in] outMarker - marker related to the outlet.
	 */
	virtual void StoreTurboPerformance(CSolver *solver, unsigned short inMarkerTP);


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
	virtual void SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config);
  
//	/*!
//	 * \brief A virtual member.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	virtual void SetSecondary_Gradient_GG(CGeometry *geometry, CConfig *config);
  
//	/*!
//	 * \brief A virtual member.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	virtual void SetSecondary_Gradient_LS(CGeometry *geometry, CConfig *config);
  
    /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config);
  
//  /*!
//	 * \brief A virtual member.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	virtual void Set_MPI_Secondary_Gradient(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPrimitive_Limiter_MPI(CGeometry *geometry, CConfig *config);
  
//  /*!
//	 * \brief A virtual member.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	virtual void SetSecondary_Limiter_MPI(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] iPoint - Index of the grid point.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetPreconditioner(CConfig *config, unsigned long iPoint);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	virtual void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                  CConfig *config, unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] StiffMatrix_Elem - Stiffness matrix of an element
	 */
	virtual void AddStiffMatrix(su2double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3 );
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                 CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	virtual void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \param[in] val_sensitivity - Value of the sensitivity coefficient.
	 */
	virtual void SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity);
    
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
	virtual void SetTotal_CDrag(su2double val_Total_CDrag);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CLift - Value of the total lift coefficient.
	 */
	virtual void SetTotal_CLift(su2double val_Total_CLift);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CT - Value of the total thrust coefficient.
	 */
	virtual void SetTotal_CT(su2double val_Total_CT);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_Total_CQ - Value of the total torque coefficient.
	 */
	virtual void SetTotal_CQ(su2double val_Total_CQ);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_Total_Heat - Value of the total heat load.
	 */
	virtual void SetTotal_HeatFlux(su2double val_Total_Heat);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_Total_MaxHeat - Value of the total heat load.
	 */
	virtual void SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat);
    
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCLift_Inv(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCLift_Visc(unsigned short val_marker);

    /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCMz_Inv(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCMz_Visc(unsigned short val_marker);
    
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CLift(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CDrag(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CSideForce(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CEff(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFx(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFy(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFz(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMx(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMy(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMz(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CLift_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CDrag_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CSideForce_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CEff_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFx_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFy_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFz_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMx_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMy_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMz_Inv(unsigned short val_marker);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCSideForce_Visc(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCDrag_Inv(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	virtual su2double GetInflow_MassFlow(unsigned short val_marker);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	virtual su2double GetExhaust_MassFlow(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the fan face pressure on the surface <i>val_marker</i>.
	 */
	virtual su2double GetInflow_Pressure(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the fan face mach on the surface <i>val_marker</i>.
	 */
	virtual su2double GetInflow_Mach(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCSideForce_Inv(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCEff_Inv(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual su2double GetCDrag_Visc(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CLift(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the sideforce coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CSideForce(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the efficiency coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CEff(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the thrust coefficient (force in the -x direction, inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CT(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the torque coefficient (moment in the -x direction, inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CQ(void);
    
    /*!
	 * \brief A virtual member.
	 * \return Value of the heat load (integrated heat flux).
	 */
	virtual su2double GetTotal_HeatFlux(void);
    
    /*!
	 * \brief A virtual member.
	 * \return Value of the heat load (integrated heat flux).
	 */
	virtual su2double GetTotal_MaxHeatFlux(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual su2double Get_PressureDrag(void);
    
    /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual su2double Get_ViscDrag(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the rotor Figure of Merit (FM) (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CMerit(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CEquivArea(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the difference of the presure and the target pressure.
	 */
	virtual su2double GetTotal_CpDiff(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the difference of the heat and the target heat.
	 */
	virtual su2double GetTotal_HeatFluxDiff(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the Free Surface coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CFreeSurface(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the FEA coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CFEA(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Near-Field Pressure coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CNearFieldOF(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	virtual void SetTotal_CEquivArea(su2double val_cequivarea);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_pressure - Value of the difference between pressure and the target pressure.
	 */
	virtual void SetTotal_CpDiff(su2double val_pressure);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_pressure - Value of the difference between heat and the target heat.
	 */
	virtual void SetTotal_HeatFluxDiff(su2double val_heat);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cfreesurface - Value of the Free Surface coefficient.
	 */
	virtual void SetTotal_CFreeSurface(su2double val_cfreesurface);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cfea - Value of the FEA coefficient.
	 */
	virtual void SetTotal_CFEA(su2double val_cfea);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
	 */
	virtual void SetTotal_CNearFieldOF(su2double val_cnearfieldpress);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CDrag(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the moment x coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CMx(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the moment y coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CMy(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the moment y coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CMz(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the force x coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CFx(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CFy(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_CFz(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the wave strength.
	 */
	virtual su2double GetTotal_CWave(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the wave strength.
	 */
	virtual su2double GetTotal_CHeat(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (inviscid contribution).
	 */
	virtual su2double GetAllBound_CLift_Inv(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual su2double GetAllBound_CDrag_Inv(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual su2double GetAllBound_CSideForce_Inv(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	virtual su2double GetAllBound_CEff_Inv(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMx_Inv(void);
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMy_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMz_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFx_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFy_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFz_Inv(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	virtual su2double GetAllBound_CLift_Visc(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	virtual su2double GetAllBound_CSideForce_Visc(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the drag coefficient (viscous contribution).
	 */
	virtual su2double GetAllBound_CDrag_Visc(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual su2double GetCPressure(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual su2double GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual void SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual su2double *GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the skin friction coefficient.
	 */
	virtual su2double GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual su2double GetHeatFlux(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual su2double GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual void SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the y plus.
	 */
	virtual su2double GetYPlus(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  virtual su2double GetStrainMag_Max(void);

  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  virtual su2double GetOmega_Max(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  virtual void SetStrainMag_Max(su2double val_strainmag_max);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  virtual void SetOmega_Max(su2double val_omega_max);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the adjoint density at the infinity.
	 */
	virtual su2double GetPsiRho_Inf(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the adjoint density at the infinity.
	 */
	virtual su2double* GetPsiRhos_Inf(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the adjoint energy at the infinity.
	 */
	virtual su2double GetPsiE_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the adjoint velocity vector at the infinity.
	 */
	virtual su2double GetPhi_Inf(unsigned short val_dim);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the geometrical sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_Sens_Geo(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the Mach sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_Sens_Mach(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the angle of attack sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_Sens_AoA(void);
    
	/*!
	 * \brief Set the total farfield pressure sensitivity coefficient.
	 * \return Value of the farfield pressure sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_Sens_Press(void);
    
	/*!
	 * \brief Set the total farfield temperature sensitivity coefficient.
	 * \return Value of the farfield temperature sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_Sens_Temp(void);

	/*!
	 * \author H. Kline
	 * \brief Set the total back pressure sensitivity coefficient.
	 * \return Value of the back pressure sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	virtual su2double GetTotal_Sens_BPress(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the density at the infinity.
	 */
	virtual su2double GetDensity_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable for the density.
	 * \return Value of the density at the infinity.
	 */
	virtual su2double GetDensity_Inf(unsigned short val_var);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the velocity at the infinity.
	 */
	virtual su2double GetModVelocity_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the density x energy at the infinity.
	 */
	virtual su2double GetDensity_Energy_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the pressure at the infinity.
	 */
	virtual su2double GetPressure_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the density x velocity at the infinity.
	 */
	virtual su2double GetDensity_Velocity_Inf(unsigned short val_dim);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \param[in] val_var - Index of the variable for the velocity.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	virtual su2double GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	virtual su2double GetVelocity_Inf(unsigned short val_dim);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the velocity at the infinity.
	 */
	virtual su2double *GetVelocity_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the viscosity at the infinity.
	 */
	virtual su2double GetViscosity_Inf(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the turbulent kinetic energy.
	 */
	virtual su2double GetTke_Inf(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the sensitivity coefficient.
	 */
	virtual su2double GetCSensitivity(unsigned short val_marker, unsigned long val_vertex);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetFreeSurface_Distance(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 * \return A pointer to an array containing a set of constants
	 */
	virtual su2double* GetConstants();

  /*!
   * \brief A virtual member.
   * \return average total pressure evaluated at an exit boundary marker
   */
  virtual su2double GetOneD_TotalPress(void);

  /*!
   * \brief A virtual member.
   * \param[in] val_exit_pt: value of the total average pressure at the exit.
   */
  virtual void SetOneD_TotalPress(su2double AveragePressure);

  /*!
   * \brief A virtual member.
   *\return average Mach number evaluated at an exit boundary marker
   */
  virtual su2double GetOneD_Mach(void);

  /*!
   * \brief A virtual member.
   * set average Mach number evaluated at an exit boundary marker
   */
  virtual void SetOneD_Mach(su2double AverageMach);
  
  /*!
   * \brief A virtual member.
   *\return average temperature evaluated at an exit boundary marker
   */
  virtual su2double GetOneD_Temp(void);
  
  /*!
   * \brief A virtual member.
   * set average temperature evaluated at an exit boundary marker
   */
  virtual void SetOneD_Temp(su2double AverageTemperature);
  
  /*!
   * \brief A virtual member.
   * \return average temperature evaluated at an exit boundary marker
   */
  virtual su2double GetOneD_MassFlowRate(void);
  
  /*!
   * \brief A virtual member.
   * set average temperature evaluated at an exit boundary marker
   */
  virtual void SetOneD_MassFlowRate(su2double MassFlowRate);
  
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged pressure at a marker.(same as area averaged pressure)
   */
  virtual su2double GetOneD_FluxAvgPress(void);
  
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged pressure at a marker. (same as area averaged pressure)
   */
  virtual void SetOneD_FluxAvgPress(su2double PressureRef);
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged density at a marker. (\f$ = (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) \f$)
   */
  virtual su2double GetOneD_FluxAvgDensity(void);
  
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged density at a marker.( \f$= (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) \f$)
   */
  virtual void SetOneD_FluxAvgDensity(su2double DensityRef);
  
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged velocity at a marker. = \f$ \sqrt ( \frac{\int((rho*u)*u^2dA)}{\int(rho*u*dA) }) \f$
   */
  virtual su2double GetOneD_FluxAvgVelocity(void);
  
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged velocity at a marker. = \f$ \sqrt (  \frac{\int((rho*u)*u^2dA)}{\int(rho*u*dA) }) \f$
   */
  virtual void SetOneD_FluxAvgVelocity(su2double VelocityRef);
  
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged enthalpy at a marker. =\f$ \frac{ \int(rho*u*h dA) }{ \int(rho *u *dA )} \f$
   */
  virtual su2double GetOneD_FluxAvgEntalpy(void);
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged enthalpy at a marker. =\f$ \frac{ \int(rho*u*h dA) }{ \int(rho *u *dA ) }\f$
   */
  virtual void SetOneD_FluxAvgEntalpy(su2double EnthalpyRef);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void GetSurface_Pressure(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] fea_config - Definition of the particular problem.
	 */
	virtual void SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry,
                           CGeometry **flow_geometry, CConfig *fea_config,
                           CConfig *flow_config, CNumerics *fea_numerics);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] solver1_geometry - Geometrical definition of the problem.
	 * \param[in] solver1_solution - Container vector with all the solutions.
	 * \param[in] solver1_config - Definition of the particular problem.
	 * \param[in] solver2_geometry - Geometrical definition of the problem.
	 * \param[in] solver2_solution - Container vector with all the solutions.
	 * \param[in] solver2_config - Definition of the particular problem.
	 */
	virtual void Copy_Zone_Solution(CSolver ***solver1_solution,
                                  CGeometry **solver1_geometry,
                                  CConfig *solver1_config,
                                  CSolver ***solver2_solution,
                                  CGeometry **solver2_geometry,
                                  CConfig *solver2_config);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] fea_config - Definition of the particular problem.
	 */
	virtual void SetInitialCondition(CGeometry **geometry,
                                   CSolver ***solver_container,
                                   CConfig *config, unsigned long ExtIter);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] flow_geometry - Geometrical definition of the problem.
	 * \param[in] flow_grid_movement - Geometrical definition of the problem.
	 * \param[in] flow_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void SetFlow_Displacement(CGeometry **flow_geometry,
                                    CVolumetricMovement *flow_grid_movement,
                                    CConfig *flow_config, CConfig *fea_config,
                                    CGeometry **fea_geometry,
                                    CSolver ***fea_solution);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void SetStruct_Displacement(CGeometry **fea_geometry,
            							CConfig *fea_config,
            							CSolver ***fea_solution);

	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void PredictStruct_Displacement(CGeometry **fea_geometry,
            								CConfig *fea_config,
            								CSolver ***fea_solution);

	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void ComputeAitken_Coefficient(CGeometry **fea_geometry,
            				  CConfig *fea_config,
            				  CSolver ***fea_solution,
            				  unsigned long iFSIIter);


	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void SetAitken_Relaxation(CGeometry **fea_geometry,
            						  CConfig *fea_config,
            						  CSolver ***fea_solution);

	/*!
	 * \brief A virtual member.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	virtual void Update_StructSolution(CGeometry **fea_geometry,
            						  CConfig *fea_config,
            						  CSolver ***fea_solution);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iter - Current external iteration number.
	 */
	virtual void LoadRestart(CGeometry **geometry, CSolver ***solver,
                           CConfig *config, int val_iter);
    
	/*!
	 * \brief Gauss method for solving a linear system.
	 * \param[in] A - Matrix Ax = b.
	 * \param[in] rhs - Right hand side.
	 * \param[in] nVar - Number of variables.
	 */
	void Gauss_Elimination(su2double** A, su2double* rhs, unsigned short nVar);
    
  /*!
  * \brief Get the number of Species present in the flow.
  */
	virtual unsigned short GetnSpecies(void);
  
  /*!
  * \brief A virtual member.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] solution - Container vector with all the solutions.
  */
	virtual void GetEngine_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  virtual void GetActuatorDisk_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);

  /*!
   * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
	 */
  virtual void SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                               CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
  * \brief Prepares and solves the aeroelastic equations.
  * \param[in] surface_movement - Surface movement classes of the problem.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] ExtIter - Physical iteration number.
  */
	void Aeroelastic(CSurfaceMovement *surface_movement, CGeometry *geometry, CConfig *config, unsigned long ExtIter);
    
  /*!
  * \brief Sets up the generalized eigenvectors and eigenvalues needed to solve the aeroelastic equations.
  * \param[in] PHI - Matrix of the generalized eigenvectors.
  * \param[in] lambda - The eigenvalues of the generalized eigensystem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetUpTypicalSectionWingModel(vector<vector<su2double> >& PHI, vector<su2double>& w, CConfig *config);
    
  /*!
  * \brief Solve the typical section wing model.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] Cl - Coefficient of lift at particular iteration.
  * \param[in] Cm - Moment coefficient about z-axis at particular iteration.
	* \param[in] config - Definition of the particular problem.
  * \param[in] val_Marker - Surface that is being monitored.
  * \param[in] displacements - solution of typical section wing model.
	*/
  
  void SolveTypicalSectionWingModel(CGeometry *geometry, su2double Cl, su2double Cm, CConfig *config, unsigned short val_Marker, vector<su2double>& displacements);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config_container - The particular config.
   */
  virtual void RegisterSolution(CGeometry *geometry, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config_container - The particular config.
   */
  virtual void RegisterOutput(CGeometry *geometry, CConfig *config);

   /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  virtual void SetAdjoint_Output(CGeometry *geometry, CConfig *config);

   /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  virtual void ExtractAdjoint_Solution(CGeometry *geometry,  CConfig *config);

  /*!
  * \brief A virtual member
  * \param[in] geometry - The geometrical definition of the problem.
  */
  virtual void RegisterObj_Func(CConfig *config);

  /*!
   * \brief  A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetSurface_Sensitivity(CGeometry *geometry, CConfig* config);

  /*!
   * \brief  A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetSensitivity(CGeometry *geometry, CConfig *config);

  virtual void SetAdj_ObjFunc(CGeometry *geometry, CConfig* config);

	/*!
	 * \brief A virtual member.
	 * \param[in] Set value of interest: 0 - Initial value, 1 - Current value.
	 */
	virtual void SetFSI_ConvValue(unsigned short val_index, su2double val_criteria);

	/*!
	 * \brief A virtual member.
	 * \param[in]  Value of interest: 0 - Initial value, 1 - Current value.
	 * \return Values to compare
	 */
	virtual su2double GetFSI_ConvValue(unsigned short val_index);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Compute_StiffMassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Compute_StiffMassDampMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void Compute_IntegrationConstants(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	virtual void SetSolution_time_n(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \return Value of the dynamic Aitken relaxation factor
	 */
	virtual su2double GetWAitken_Dyn(void);

	/*!
	 * \brief A virtual member.
	 * \return Value of the last Aitken relaxation factor in the previous time step.
	 */
	virtual su2double GetWAitken_Dyn_tn1(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] Value of the dynamic Aitken relaxation factor
	 */
	virtual void SetWAitken_Dyn(su2double waitk);

	/*!
	 * \brief A virtual member.
	 * \param[in] Value of the last Aitken relaxation factor in the previous time step.
	 */
	virtual void SetWAitken_Dyn_tn1(su2double waitk_tn1);

  /*!
   * \brief A virtual member.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] output - Reset the rhs vector.
   */
  virtual unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);

  /*!
   * \brief A virtual member.
   * \param[in] Value of freestream pressure.
   */
  virtual void SetPressure_Inf(su2double p_inf);

  /*!
   * \brief A virtual member.
   * \param[in] Value of freestream temperature.
   */
  virtual void SetTemperature_Inf(su2double t_inf);

  /*!
   * \brief A virtual member.
   * \param[in] kind_recording - Kind of AD recording.
   */
  virtual void SetRecording(CGeometry *geometry, CConfig *config, unsigned short kind_recording);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void RegisterVariables(CGeometry *geometry, CConfig *config, bool reset = false);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config);
};

/*!
 * \class CBaselineSolver
 * \brief Main class for defining a baseline solution from a restart file (for output).
 * \author F. Palacios, T. Economon.
 * \version 4.0.2 "Cardinal"
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
   * \param[in] solver - Container vector with all of the solvers.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iter - Current external iteration number.
	 */
	void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CBaselineSolver(void);
    
};

/*!
 * \class CEulerSolver
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
 */
class CEulerSolver : public CSolver {
protected:
	
  su2double
  Mach_Inf,	/*!< \brief Mach number at the infinity. */
	Density_Inf,	/*!< \brief Density at the infinity. */
	Energy_Inf,			/*!< \brief Energy at the infinity. */
  Temperature_Inf,			/*!< \brief Energy at the infinity. */
	Pressure_Inf,		/*!< \brief Pressure at the infinity. */
	*Velocity_Inf;		/*!< \brief Flow Velocity vector at the infinity. */
	
  su2double
  *CDrag_Inv,	/*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
	*CLift_Inv,			/*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
	*CSideForce_Inv,		/*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
	*CMx_Inv,			/*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
	*CMy_Inv,			/*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
	*CMz_Inv,			/*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
	*CFx_Inv,			/*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
	*CFy_Inv,			/*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
	*CFz_Inv,			/*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CLift_Inv, /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CDrag_Inv, /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSideForce_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
	*CEff_Inv,				/*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
	*CMerit_Inv,				/*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
	*CT_Inv,			/*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
	*CQ_Inv,			/*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
	*CEquivArea_Inv,				/*!< \brief Equivalent area (inviscid contribution) for each boundary. */
	*CNearFieldOF_Inv,				/*!< \brief Near field pressure (inviscid contribution) for each boundary. */
	**CPressure,		/*!< \brief Pressure coefficient for each boundary and vertex. */
	**CPressureTarget,		/*!< \brief Target Pressure coefficient for each boundary and vertex. */
	**HeatFlux,		/*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **HeatFluxTarget,		/*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **YPlus,		/*!< \brief Yplus for each boundary and vertex. */
  ***CharacPrimVar,		/*!< \brief Value of the characteristic variables at each boundary. */
	*ForceInviscid,		/*!< \brief Inviscid force for each boundary. */
	*MomentInviscid,	/*!< \brief Inviscid moment for each boundary. */
	*Inflow_MassFlow,	/*!< \brief Mass flow rate for each boundary. */
  *Bleed_MassFlow,	/*!< \brief Mass flow rate for each boundary. */
	*Exhaust_MassFlow,	/*!< \brief Mass flow rate for each boundary. */
	*Inflow_Pressure,	/*!< \brief Fan face pressure for each boundary. */
	*Inflow_Mach,	/*!< \brief Fan face mach number for each boundary. */
	*Inflow_Area,	/*!< \brief Boundary total area. */
  *Bleed_Pressure,	/*!< \brief Fan face pressure for each boundary. */
  *Bleed_Temperature,	/*!< \brief Fan face mach number for each boundary. */
  *Bleed_Area,	/*!< \brief Boundary total area. */
  *Exhaust_Area,	/*!< \brief Boundary total area. */
  *Exhaust_Pressure,	/*!< \brief Fan face pressure for each boundary. */
  *Exhaust_Temperature,	/*!< \brief Fan face mach number for each boundary. */
  Inflow_MassFlow_Total,	/*!< \brief Mass flow rate for each boundary. */
  Bleed_MassFlow_Total,	/*!< \brief Mass flow rate for each boundary. */
  Exhaust_MassFlow_Total,	/*!< \brief Mass flow rate for each boundary. */
	Inflow_Pressure_Total,	/*!< \brief Fan face pressure for each boundary. */
	Inflow_Mach_Total,	/*!< \brief Fan face mach number for each boundary. */
  Bleed_Pressure_Total,	/*!< \brief Fan face pressure for each boundary. */
  Bleed_Temperature_Total,	/*!< \brief Fan face mach number for each boundary. */
	InverseDesign;	/*!< \brief Inverse design functional for each boundary. */
	
  su2double
  AllBound_CDrag_Inv,	/*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
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
	
  su2double
  OneD_TotalPress, /*!< \brief average total pressure evaluated at an exit */
  OneD_Mach, /*!< \brief area average Mach evaluated at an exit */
  OneD_Temp, /*!< \brief area average Temperature evaluated at an exit */
  OneD_PressureRef, /*!< \brief area average Pressure evaluated at an exit */
  OneD_MassFlowRate, /*!< \brief Mass flow rate at an exit */
  OneD_DensityRef, /*!< \brief flux average density evaluated at an exit */
  OneD_EnthalpyRef, /*!< \brief flux average enthalpy evaluated at an exit */
  OneD_VelocityRef, /*!< \brief flux average velocity evaluated at an exit */
  Total_CDrag, /*!< \brief Total drag coefficient for all the boundaries. */
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
  Total_Heat,    /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat, /*!< \brief Maximum heat flux on all boundaries. */
	Total_CEquivArea,			/*!< \brief Total Equivalent Area coefficient for all the boundaries. */
	Total_CNearFieldOF,			/*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
  Total_CFreeSurface,			/*!< \brief Total Free Surface coefficient for all the boundaries. */
  Total_CpDiff,			/*!< \brief Total Equivalent Area coefficient for all the boundaries. */
	Total_HeatFluxDiff,			/*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_MassFlowRate;     /*!< \brief Total Mass Flow Rate on monitored boundaries. */
  su2double *Surface_CLift,   /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CDrag,          /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSideForce,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CEff,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,            /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,            /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,            /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,            /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,            /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz;            /*!< \brief z Moment coefficient for each monitoring surface. */
	su2double *iPoint_UndLapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*jPoint_UndLapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	su2double *SecondaryVar_i,	/*!< \brief Auxiliary vector for storing the solution at point i. */
	*SecondaryVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	su2double *PrimVar_i,	/*!< \brief Auxiliary vector for storing the solution at point i. */
	*PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	su2double **LowMach_Precontioner; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
	bool space_centered,  /*!< \brief True if space centered scheeme used. */
	euler_implicit,			/*!< \brief True if euler implicit scheme used. */
	least_squares;        /*!< \brief True if computing gradients by least squares. */
	su2double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	su2double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
  
  su2double *Primitive,		/*!< \brief Auxiliary nPrimVar vector. */
	*Primitive_i,				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
	*Primitive_j;				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */
  
  su2double *Secondary,		/*!< \brief Auxiliary nPrimVar vector. */
	*Secondary_i,				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
	*Secondary_j;				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  su2double Cauchy_Value,	/*!< \brief Summed value of the convergence indicator. */
	Cauchy_Func;			/*!< \brief Current value of the convergence indicator at one iteration. */
	unsigned short Cauchy_Counter;	/*!< \brief Number of elements of the Cauchy serial. */
	su2double *Cauchy_Serie;			/*!< \brief Complete Cauchy serial. */
	su2double Old_Func,	/*!< \brief Old value of the objective function (the function which is monitored). */
	New_Func;			/*!< \brief Current value of the objective function (the function which is monitored). */
  su2double AoA_old;  /*!< \brief Old value of the angle of attack (monitored). */

  CFluidModel  *FluidModel;  /*!< \brief fluid model used in the solver */
  su2double **AveragedVelocity,
  	  	 **AveragedNormal,
		 **AveragedGridVel,
  	  	  **AveragedFlux,
		  **TotalFlux,
		  *TotalArea,
		  *AveragedNormalVelocity,
		  *ExtAveragedNormalVelocity,
		  *AveragedTangVelocity,
		  *ExtAveragedTangVelocity,
		  *AveragedTangGridVelocity,
		  *AveragedMach,
		  *AveragedNormalMach,
		  *AveragedTangMach,
		  *AveragedEnthalpy,
		  *AveragedPressure,
		  *AveragedTotTemperature,
		  *AveragedTotPressure,
		  *ExtAveragedPressure,
		  *ExtAveragedTotTemperature,
		  *ExtAveragedTotPressure,
		  *AveragedDensity,
		  *ExtAveragedDensity,
		  *AveragedSoundSpeed,
		  *AveragedEntropy,
		  *MassFlow,
		  *FlowAngle;
  su2double *TotalStaticEfficiency,
  	  	  	*TotalTotalEfficiency,
			*KineticEnergyLoss,
			*TotalPressureLoss,
  	  	  	*MassFlowIn,
			*MassFlowOut,
			*FlowAngleIn,
			*FlowAngleOut,
			*EulerianWork,
			*TotalEnthalpyIn,
			*PressureRatio,
			*PressureOut,
			*EnthalpyOut,
			*MachIn,
			*MachOut,
			*NormalMachIn,
			*NormalMachOut,
			*VelocityOutIs;



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
  void Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config);
  
//  /*!
//	 * \brief Impose the send-receive boundary condition.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//  void Set_MPI_Secondary_Limiter(CGeometry *geometry, CConfig *config);

  /*!
	 * \brief Set the fluid solver nondimensionalization.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Compute the pressure at the infinity.
	 * \return Value of the pressure at the infinity.
	 */
	CFluidModel* GetFluidModel(void);


    /*!
	 * \brief Compute the density at the infinity.
	 * \return Value of the density at the infinity.
	 */
	su2double GetDensity_Inf(void);
    
	/*!
	 * \brief Compute 2-norm of the velocity at the infinity.
	 * \return Value of the 2-norm of the velocity at the infinity.
	 */
	su2double GetModVelocity_Inf(void);
    
	/*!
	 * \brief Compute the density multiply by energy at the infinity.
	 * \return Value of the density multiply by  energy at the infinity.
	 */
	su2double GetDensity_Energy_Inf(void);
    
	/*!
	 * \brief Compute the pressure at the infinity.
	 * \return Value of the pressure at the infinity.
	 */
	su2double GetPressure_Inf(void);

	/*!
	 * \brief Compute the density multiply by velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	su2double GetDensity_Velocity_Inf(unsigned short val_dim);
    
	/*!
	 * \brief Get the velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	su2double GetVelocity_Inf(unsigned short val_dim);
  
  /*!
	 * \brief Get the velocity at the infinity.
	 * \return Value of the velocity at the infinity.
	 */
	su2double *GetVelocity_Inf(void);
  
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
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Compute the extrapolated quantities, for MUSCL upwind 2nd reconstruction,
	 * in a more thermodynamic consistent way
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeConsExtrapolation(CConfig *config);
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
   * \brief Compute primitive variables and their gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh);
  
  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);

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
	void SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config);
  
//  /*!
//	 * \brief Compute the gradient of the primitive variables using Green-Gauss method,
//	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	void SetSecondary_Gradient_GG(CGeometry *geometry, CConfig *config);
  
//	/*!
//	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
//	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	void SetSecondary_Gradient_LS(CGeometry *geometry, CConfig *config);
    
  /*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Compute the limiter of the primitive variables.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config);
  
//  /*!
//	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
//	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	void Set_MPI_Secondary_Gradient(CGeometry *geometry, CConfig *config);
  
//	/*!
//	 * \brief Compute the limiter of the primitive variables.
//	 * \param[in] geometry - Geometrical definition of the problem.
//	 * \param[in] config - Definition of the particular problem.
//	 */
//	void SetSecondary_Limiter(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
	 * \param[in] iPoint - Index of the grid point
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPreconditioner(CConfig *config, unsigned long iPoint);
    
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
	 * \author: G.Gori, S.Vitale, M.Pini, A.Guardone, P.Colonna
	 *
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                      CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the interface boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config);
    
	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config);
  
  /*!
	 * \brief Impose the actuator disk boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_ActDisk_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config);
  
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
	 * \author: G.Gori, S.Vitale, M.Pini, A.Guardone, P.Colonna
	 *
	 * \brief Impose the boundary condition using characteristic recostruction.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

	/*!
	 * \author: G.Gori, S.Vitale, M.Pini, A.Guardone, P.Colonna
	 *
	 * \brief Impose the boundary condition using characteristic recostruction.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NonReflecting(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);


	/*!
	 * \brief Impose a subsonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose a supersonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose a supersonic outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                           CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                   CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
     
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                   CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the nacelle inflow boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
  
   /*!
    * \brief Impose the nacelle bleed boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
  void BC_Engine_Bleed(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                        CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Impose the ancelle exhaust boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] conv_numerics - Description of the numerical method.
     * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);

	/*!
	 * \brief It avarage the fluxes value along a boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void Mixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker);

	/*!
	 * \brief It avarage the fluxes value along a boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void MPIMixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short marker_flag);

	/*!
	 * \brief it performs a mixed out average of the nodes of a boundary.
	 * \param[in] val_init_pressure -  initial pressure value
	 * \param[in] val_Averaged_Flux - flux averaged values.
     * \param[in] val_normal - normal vector.
     * \param[in] pressure_mix - value of the mixed-out avaraged pressure.
	 * \param[in] density_miz - value of the mixed-out avaraged density.
	 */
	void MixedOut_Average (su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *pressure_mix, su2double *density_mix);

	/*!
	 * \brief it finds the root of an implicit equation that relates pressure and density.
	 * \param[in] pressure - pressure value
	 * \param[in] val_Averaged_Flux - flux averaged values.
     * \param[in] val_normal - normal vector.
     * \param[in] valfunc - Description of the numerical method.
	 * \param[in] density - value of the mixed-out avaraged density.
	 */
	void MixedOut_Root_Function(su2double *pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *valfunc, su2double *density);

	/*!
	 * \brief it performs a fourier transformation of a characteristic value.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in]  c4k - Fourier transformation coefficients.
	 * \param[in]  nboundaryvertex - pithcwise ordered vertex.
	 */
	void Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> > &c4k,signed long &nboundaryvertex);

	/*!
	 * \brief it performs a fourier transformation of a characteristic value.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 * \param[in]  c2k - Fourier transformation coefficients.
	 * \param[in]  c3k - Fourier transformation coefficients.
	 * \param[in]  nboundaryvertex - pithcwise ordered vertex.
	 */
	void Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> >& c2k,vector<std::complex<su2double> >& c3k,signed long& nboundaryvertex);

	/*!
	 * \brief A virtual member.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] intMarker - internal marker.
	 * \param[in] extMarker - external marker.
	 */
	 void SetExtAveragedValue(CSolver *solver_container, unsigned short intMarker,  unsigned short extMarker);

	 /*!
	  * \brief Provide the average density at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Density on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedDensity(unsigned short valMarker);

	 /*!
	  * \brief Provide the average pressure at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Pressure on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedTotPressure(unsigned short valMarker);

	 /*!
	  * \brief Provide Total Pressure Losses (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of Total Pressure Losses.
	  */
	 su2double GetTotalPressureLoss(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide Kinetic Energy Losses (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Kinetic Energy Losses.
	  */
	 su2double GetKineticEnergyLoss(unsigned short inMarkerTP);

	 /*!
		* \brief Provide Total-Total Efficiency (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Total-Total Efficiency.
	  */
	 su2double GetTotalTotalEfficiency(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide Total-Static Efficiency (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Total-Static Efficiency.
	  */
	 su2double GetTotalStaticEfficiency(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Eulerian Work (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Eulerian Work.
	  */
	 su2double GetEulerianWork(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Inlet Total Enthalpy (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Inlet Total Enthalpy.
	  */
	 su2double GetTotalEnthalpyIn(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Inlet Flow Angle (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Inlet Flow Angle.
	  */
	 su2double GetFlowAngleIn(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Outlet Flow Angle (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Outlet FLow Angle.
	  */
	 su2double GetFlowAngleOut(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Inlet Mass Flow (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Inlet Mass Flow.
	  */
	 su2double GetMassFlowIn(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Outlet Mass Flow (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Outlet Mass Flow.
	  */
	 su2double GetMassFlowOut(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Inlet Mach number (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Inlet Mach number.
	  */
	 su2double GetMachIn(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Outlet Mach number (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Outlet Mach number.
	  */
	 su2double GetMachOut(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the normal component of the Inlet Mach number (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the normal component of the Inlet Mach number.
	  */
	 su2double GetNormalMachIn(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the normal component of the Outlet Mach number (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the normal component of the Outlet Mach number.
	  */
	 su2double GetNormalMachOut(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Outlet Enthalpy (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Outlet Enthalpy.
	  */
	 su2double GetEnthalpyOut(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Isentropic Outlet Velocity (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Isentropic Outlet Velocity.
	  */
	 su2double GetVelocityOutIs(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide the Outlet Pressure (turbomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Outlet Pressure.
	  */
	 su2double GetPressureOut(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide Pressure Ratio (tubomachinery performance).
	  * \param[in] inMarkerTP - turboperformance marker.
	  * \return Value of the Pressure Ratio.
	  */
	 su2double GetPressureRatio(unsigned short inMarkerTP);

	 /*!
	  * \brief Provide Averaged Total Temperature at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Density on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedTotTemperature(unsigned short valMarker);

	 /*!
	  * \brief Provide the Average pressure at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Pressure on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedPressure(unsigned short valMarker);

	 /*!
	  * \brief Provide the MassFlow at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the MassFLow on the surface <i>val_marker</i>.
	  */
	 su2double GetMassFlow(unsigned short valMarker);

	 /*!
	  * \brief Provide the Flow Angle at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Flow Angle on the surface <i>val_marker</i>.
	  */
	 su2double GetFlowAngle(unsigned short valMarker);

	 /*!
	  * \brief Provide the Mach number at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Mach number on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedMach(unsigned short valMarker);

	 /*!
	  * \brief Provide the Normal Mach number at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Normal Mach number on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedNormalMach(unsigned short valMarker);

	 /*!
	  * \brief Provide the average enthalpy at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Enthalpy on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedEnthalpy(unsigned short valMarker);

	 /*!
	  * \brief Provide the average pressure at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Entropy on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedEntropy(unsigned short valMarker);

	 /*!
	  * \brief Provide the average pressure at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Components of the Average Velocity on the surface <i>val_marker</i>.
	  */
	 su2double* GetAveragedVelocity(unsigned short valMarker);

	 /*!
	  * \brief Provide the average grid velocity at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Grid Velocity on the surface <i>val_marker</i>.
	  */
	 su2double* GetAveragedGridVelocity(unsigned short valMarker);

	 /*!
	  * \brief Provide the Average Normal Velocity at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Normal Velocity on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedNormalVelocity(unsigned short valMarker);

	 /*!
	  * \brief Provide the Tangent Velocity at the boundary of interest.
	  * \param[in] val_marker - bound marker.
	  * \return Value of the Average Tangent Velocity on the surface <i>val_marker</i>.
	  */
	 su2double GetAveragedTangVelocity(unsigned short valMarker);

	 /*!
	 * \brief compare to values.
	 * \param[in] a - value 1.
	 * \param[in] b - value 2.
	 */
     static bool Compareval(std::vector<su2double> a,std::vector<su2double> b);

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
	void GetEngine_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
   * \brief Compute the Fan face Mach number.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  void GetActuatorDisk_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
	 * \brief Update the AoA and freestream velocity at the farfield.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
	 */
  void SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh, bool Output);
  
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
	 * \brief Compute turbomachinery performance.
	 * \param[in] solver - solver containing the outlet information.
	 * \param[in] inMarker - marker related to the inlet.
	 * \param[in] outMarker - marker related to the outlet.
	 */
	void TurboPerformance(CSolver *solver,  CConfig *config, unsigned short inMarker,  unsigned short outMarker, unsigned short Kind_TurboPerf , unsigned short inMarkerTP );

	/*!
	 * \brief Compute turbomachinery performance.
	 * \param[in] config - contains config file information.
	 */
	void MPITurboPerformance(CConfig *config);

	/*!
	 * \brief Compute turbomachinery performance.
	 * \param[in] solver - solver containing the outlet information.
	 * \param[in] inMarker - marker related to the inlet.
	 * \param[in] outMarker - marker related to the outlet.
	 */
	void StoreTurboPerformance(CSolver *solver,  unsigned short inMarkerTP );

	/*!
	 * \brief Provide the non dimensional lift coefficient (inviscid contribution).
	 * \param val_marker Surface where the coefficient is going to be computed.
	 * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCLift_Inv(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional z moment coefficient (inviscid contribution).
	 * \param val_marker Surface where the coefficient is going to be computed.
	 * \return Value of the z moment coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCMz_Inv(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional lift coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient on the surface <i>val_marker</i>.
	 */
	su2double GetSurface_CLift(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional drag coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient on the surface <i>val_marker</i>.
	 */
	su2double GetSurface_CDrag(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSideForce(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CEff(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFx(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFy(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFz(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional x moment coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
	 */
	su2double GetSurface_CMx(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional y moment coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
	 */
	su2double GetSurface_CMy(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional z moment coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
	 */
	su2double GetSurface_CMz(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CLift_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CDrag_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSideForce_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CEff_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFx_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFy_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFz_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMx_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMy_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMz_Inv(unsigned short val_marker);
  
	/*!
	 * \brief Provide the non dimensional drag coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCDrag_Inv(unsigned short val_marker);
    
	/*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	su2double GetInflow_MassFlow(unsigned short val_marker);
    
    /*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the mass flow rate on the surface <i>val_marker</i>.
	 */
	su2double GetExhaust_MassFlow(unsigned short val_marker);
    
	/*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the fan face pressure on the surface <i>val_marker</i>.
	 */
	su2double GetInflow_Pressure(unsigned short val_marker);
    
	/*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the fan face mach on the surface <i>val_marker</i>.
	 */
	su2double GetInflow_Mach(unsigned short val_marker);
    
	/*!
	 * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCSideForce_Inv(unsigned short val_marker);
    
	/*!
	 * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCEff_Inv(unsigned short val_marker);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
	 * \return Value of the sideforce coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CSideForce(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
	 * \return Value of the efficiency coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CEff(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CEquivArea(void);
  
  /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CpDiff(void);
  
  /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_HeatFluxDiff(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Near-Field pressure coefficient.
	 * \return Value of the NearField pressure coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CNearFieldOF(void);
    
	/*!
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_CEquivArea(su2double val_cequivarea);
  
  /*!
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_CpDiff(su2double val_pressure);
  
  /*!
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_HeatFluxDiff(su2double val_heat);
    
	/*!
	 * \brief Set the value of the Near-Field pressure oefficient.
	 * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
	 */
	void SetTotal_CNearFieldOF(su2double val_cnearfieldpress);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
	 * \param[in] val_Total_CLift - Value of the total lift coefficient.
	 */
	void SetTotal_CLift(su2double val_Total_CLift);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
	 * \return Value of the lift coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CLift(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
	 * \return Value of the drag coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CDrag(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
	 * \return Value of the moment x coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CMx(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
	 * \return Value of the moment y coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CMy(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
	 * \return Value of the moment z coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CMz(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
	 * \return Value of the force x coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CFx(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
	 * \return Value of the force y coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CFy(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
	 * \return Value of the force z coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CFz(void);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional thrust coefficient.
	 * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CT(void);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional thrust coefficient.
	 * \param[in] val_Total_CT - Value of the total thrust coefficient.
	 */
	void SetTotal_CT(su2double val_Total_CT);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional torque coefficient.
	 * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CQ(void);
    
    /*!
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	su2double GetTotal_HeatFlux(void);
    
    /*!
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	su2double GetTotal_MaxHeatFlux(void);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional torque coefficient.
	 * \param[in] val_Total_CQ - Value of the total torque coefficient.
	 */
	void SetTotal_CQ(su2double val_Total_CQ);
    
    /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Heat - Value of the heat load.
	 */
	void SetTotal_HeatFlux(su2double val_Total_Heat);
    
    /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Heat - Value of the heat load.
	 */
	void SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat);
    
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional rotor Figure of Merit.
	 * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CMerit(void);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
	 * \param[in] val_Total_CDrag - Value of the total drag coefficient.
	 */
	void SetTotal_CDrag(su2double val_Total_CDrag);
    
	/*!
	 * \brief Get the inviscid contribution to the lift coefficient.
	 * \return Value of the lift coefficient (inviscid contribution).
	 */
	su2double GetAllBound_CLift_Inv(void);
    
	/*!
	 * \brief Get the inviscid contribution to the drag coefficient.
	 * \return Value of the drag coefficient (inviscid contribution).
	 */
	su2double GetAllBound_CDrag_Inv(void);
    
	/*!
	 * \brief Get the inviscid contribution to the sideforce coefficient.
	 * \return Value of the sideforce coefficient (inviscid contribution).
	 */
	su2double GetAllBound_CSideForce_Inv(void);
    
	/*!
	 * \brief Get the inviscid contribution to the efficiency coefficient.
	 * \return Value of the efficiency coefficient (inviscid contribution).
	 */
	su2double GetAllBound_CEff_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMz_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFz_Inv(void);
  
	/*!
	 * \brief Provide the Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	su2double GetCPressure(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief Provide the Target Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	su2double GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief Set the value of the target Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
  void SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure);

  /*!
	 * \brief Value of the characteristic variables at the boundaries.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	su2double *GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex);
  
	/*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Free Surface coefficient.
	 * \return Value of the Free Surface coefficient (inviscid + viscous contribution).
	 */
	su2double GetTotal_CFreeSurface(void);
  
	/*!
	 * \brief Set the value of the Free Surface coefficient.
	 * \param[in] val_cfreesurface - Value of the Free Surface coefficient.
	 */
	void SetTotal_CFreeSurface(su2double val_cfreesurface);
  
	/*!
   * \brief Provide the averaged total pressure at a marker.
   */
	su2double GetOneD_TotalPress(void);
  
	/*!
   * \brief Set the value of averaged total pressure
   * \param[in] val_exit_pt - value of the averaged pressure
   */
	void SetOneD_TotalPress(su2double AveragePressure);
  
  /*!
   * \brief Provide the averaged Mach number at a marker.
   */
  su2double GetOneD_Mach(void);
  
  /*!
   * \brief Set the averaged Mach number at a marker.
   */
  void SetOneD_Mach(su2double AverageMach);
  
  /*!
   * \brief Provide the averaged Mach number at a marker.
   */
  su2double GetOneD_Temp(void);
  
  /*!
   * \brief Set the averaged Temperature at a marker.
   */
  void SetOneD_Temp(su2double AverageTemperature);
  
  /*!
   * \brief Provide the averaged Mach number at a marker.
   */
  su2double GetOneD_MassFlowRate(void);
  
  /*!
   * \brief Set the averaged Temperature at a marker.
   */
  void SetOneD_MassFlowRate(su2double MassFlowRate);
  
  /*!
   * \brief Get the flux averaged pressure at a marker.(same as area averaged pressure)
   */
  su2double GetOneD_FluxAvgPress(void);
  
  /*!
   * \brief Set the flux averaged pressure at a marker. (same as area averaged pressure)
   */
  void SetOneD_FluxAvgPress(su2double PressureRef);
  
  /*!
   * \brief Get the flux averaged density at a marker. ( = (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) )
   */
  su2double GetOneD_FluxAvgDensity(void);
  
  /*!
   * \brief Set the flux averaged density at a marker.( = (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) )
   */
  void SetOneD_FluxAvgDensity(su2double DensityRef);
  
  /*!
   * \brief Get the flux averaged velocity at a marker. = \f$ \sqrt ( \int((rho*u)*u^2dA)/\int(rho*u*dA) )\f$
   */
  su2double GetOneD_FluxAvgVelocity(void);
  
  /*!
   * \brief Set the flux averaged velocity at a marker. =\f$ sqrt ( \int((rho*u)*u^2dA)/\int(rho*u*dA) ) \f$
   */
  void SetOneD_FluxAvgVelocity(su2double VelocityRef);
  
  /*!
   * \brief Get the flux averaged enthalpy at a marker. = \f$ \int(rho*u*h dA) / \int(rho *u *dA ) \f$
   */
  su2double GetOneD_FluxAvgEntalpy(void);
  
  /*!
   * \brief Set the flux averaged enthalpy at a marker. =\f$ \int(rho*u*h dA) / \int(rho *u *dA ) \f$
   */
  void SetOneD_FluxAvgEntalpy(su2double EnthalpyRef);
  
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
	 * \brief Load a solution from a restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iter - Current external iteration number.
	 */
	void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter);
    
	/*!
	 * \brief Set the initial condition for the Euler Equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] ExtIter - External iteration.
	 */
	void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
  
	/*!
	 * \brief Recompute distance to the level set 0.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetFreeSurface_Distance(CGeometry *geometry, CConfig *config);
  

  /*!
   * \brief Set the freestream pressure.
   * \param[in] Value of freestream pressure.
   */
  void SetPressure_Inf(su2double p_inf);

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  void SetTemperature_Inf(su2double t_inf);
};

/*!
 * \class CNSSolver
 * \brief Main class for defining the Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
 */
class CNSSolver : public CEulerSolver {
private:
	su2double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;	/*!< \brief Turbulent kinetic energy at the infinity. */
	su2double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
	Prandtl_Turb;         /*!< \brief Turbulent Prandtl number. */
	su2double *CDrag_Visc,	/*!< \brief Drag coefficient (viscous contribution) for each boundary. */
	*CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for each boundary. */
	*CSideForce_Visc,		/*!< \brief Side force coefficient (viscous contribution) for each boundary. */
	*CMx_Visc,			/*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
	*CMy_Visc,			/*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
	*CMz_Visc,			/*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
	*CFx_Visc,			/*!< \brief Force x coefficient (viscous contribution) for each boundary. */
	*CFy_Visc,			/*!< \brief Force y coefficient (viscous contribution) for each boundary. */
	*CFz_Visc,			/*!< \brief Force z coefficient (viscous contribution) for each boundary. */
  *Surface_CLift_Visc,/*!< \brief Lift coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CDrag_Visc,/*!< \brief Drag coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CSideForce_Visc,/*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CEff_Visc,/*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFx_Visc,  /*!< \brief Force x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFy_Visc,  /*!< \brief Force y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFz_Visc,  /*!< \brief Force z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMx_Visc,  /*!< \brief Moment x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMy_Visc,  /*!< \brief Moment y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMz_Visc,  /*!< \brief Moment z coefficient (viscous contribution) for each monitoring surface. */
	*CEff_Visc,			/*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
	*CMerit_Visc,			/*!< \brief Rotor Figure of Merit (Viscous contribution) for each boundary. */
	*CT_Visc,		/*!< \brief Thrust coefficient (viscous contribution) for each boundary. */
	*CQ_Visc,		/*!< \brief Torque coefficient (viscous contribution) for each boundary. */
  *Heat_Visc,		/*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHeatFlux_Visc, /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
	**CSkinFriction;	/*!< \brief Skin friction coefficient for each boundary and vertex. */
	su2double *ForceViscous,	/*!< \brief Viscous force for each boundary. */
	*MomentViscous;			/*!< \brief Inviscid moment for each boundary. */
	su2double AllBound_CDrag_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
	AllBound_CLift_Visc,		/*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
	AllBound_CSideForce_Visc,		/*!< \brief Sideforce coefficient (viscous contribution) for all the boundaries. */
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
  AllBound_HeatFlux_Visc,		/*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHeatFlux_Visc; /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
  su2double StrainMag_Max, Omega_Max; /*!< \brief Maximum Strain Rate magnitude and Omega. */
  
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
	su2double GetViscosity_Inf(void);
  
  /*!
	 * \brief Get the turbulent kinetic energy at the infinity.
	 * \return Value of the turbulent kinetic energy at the infinity.
	 */
	su2double GetTke_Inf(void);
    
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
      
  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);

    /*!
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
    /*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
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
	su2double GetCLift_Visc(unsigned short val_marker);
    
    /*!
	 * \brief Get the non dimensional z moment coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCMz_Visc(unsigned short val_marker);
  
  /*!
	 * \brief Get the non dimensional sideforce coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCSideForce_Visc(unsigned short val_marker);
    
	/*!
	 * \brief Get the non dimensional drag coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	su2double GetCDrag_Visc(unsigned short val_marker);
    
	/*!
	 * \brief Get the total non dimensional lift coefficient (viscous contribution).
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	su2double GetAllBound_CLift_Visc(void);
  
  /*!
	 * \brief Get the total non dimensional sideforce coefficient (viscous contribution).
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	su2double GetAllBound_CSideForce_Visc(void);
    
	/*!
	 * \brief Get the total non dimensional drag coefficient (viscous contribution).
	 * \return Value of the drag coefficient (viscous contribution).
	 */
	su2double GetAllBound_CDrag_Visc(void);
    
	/*!
	 * \brief Compute the viscous residuals.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	su2double GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex);
    
	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	su2double GetHeatFlux(unsigned short val_marker, unsigned long val_vertex);
	
  /*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	su2double GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
	 * \brief Set the value of the target Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
  void SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat);
  
	/*!
	 * \brief Get the y plus.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the y plus.
	 */
	su2double GetYPlus(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Get the max Omega.
   * \return Value of the max Omega.
   */
  su2double GetOmega_Max(void);
  
  /*!
   * \brief Get the max Strain rate magnitude.
   * \return Value of the max Strain rate magnitude.
   */
  su2double GetStrainMag_Max(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  void SetStrainMag_Max(su2double val_strainmag_max);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  void SetOmega_Max(su2double val_omega_max);
  
};

/*!
 * \class CTurbSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 4.0.2 "Cardinal"
 */
class CTurbSolver : public CSolver {
protected:
	su2double *FlowPrimVar_i,  /*!< \brief Store the flow solution at point i. */
	*FlowPrimVar_j,         /*!< \brief Store the flow solution at point j. */
	*lowerlimit,            /*!< \brief contains lower limits for turbulence variables. */
	*upperlimit;            /*!< \brief contains upper limits for turbulence variables. */
	su2double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	su2double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */
    
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
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
  
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short iMesh);
  
	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                        CConfig *config, unsigned short iMesh, unsigned short iRKStep);
  
	/*!
	 * \brief Impose the Symmetry Plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
    /*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
  
};

/*!
 * \class CTurbSASolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 4.0.2 "Cardinal"
 */

class CTurbSASolver: public CTurbSolver {
private:
	su2double nu_tilde_Inf, nu_tilde_Engine;
	
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
	CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel);
    
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
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
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
    /*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
    /*!
	 * \brief Impose the engine inflow boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose the engine bleed boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Bleed(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                        CConfig *config, unsigned short val_marker);

	/*!
	 * \brief Impose the engine exhaust boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
    /*!
	 * \brief Impose the interface boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config);
    
	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config);
  
  /*!
	 * \brief Load a solution from a restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iter - Current external iteration number.
	 */
	void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter);

    
};

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Aranake.
 * \version 4.0.2 "Cardinal"
 */

class CTransLMSolver: public CTurbSolver {
private:
	su2double Intermittency_Inf, REth_Inf;
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
	 * \brief Restart residual and compute gradients.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                       CNumerics *numerics, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Compute the viscous residuals for the turbulent equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
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
	su2double *LinSysSolItmc;		/*!< \brief vector to store iterative solution of implicit linear system. */
	su2double *LinSysResItmc;		/*!< \brief vector to store iterative residual of implicit linear system. */
	su2double *rhsItmc;		/*!< \brief right hand side of implicit linear system. */
	CSysMatrix JacobianReth; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
	su2double *LinSysSolReth;		/*!< \brief vector to store iterative solution of implicit linear system. */
	su2double *LinSysResReth;		/*!< \brief vector to store iterative residual of implicit linear system. */
	su2double *rhsReth;		/*!< \brief right hand side of implicit linear system. */
};

/*!
 * \class CTurbSSTSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Campos, F. Palacios, T. Economon
 * \version 4.0.2 "Cardinal"
 */

class CTurbSSTSolver: public CTurbSolver {
private:
	su2double *constants,  /*!< \brief Constants for the model. */
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
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
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
    /*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                            unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
    
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Get the constants for the SST model.
	 * \return A pointer to an array containing a set of constants
	 */
	su2double* GetConstants();
    
};

/*!
 * \class CAdjEulerSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
 */
class CAdjEulerSolver : public CSolver {
protected:
	su2double PsiRho_Inf,	/*!< \brief PsiRho variable at the infinity. */
	PsiE_Inf,			/*!< \brief PsiE variable at the infinity. */
	*Phi_Inf;			/*!< \brief Phi vector at the infinity. */
	su2double *Sens_Mach, /*!< \brief Mach sensitivity coefficient for each boundary. */
	*Sens_AoA,			/*!< \brief Angle of attack sensitivity coefficient for each boundary. */
	*Sens_Geo,			/*!< \brief Shape sensitivity coefficient for each boundary. */
	*Sens_Press,			/*!< \brief Pressure sensitivity coefficient for each boundary. */
	*Sens_Temp,			/*!< \brief Temperature sensitivity coefficient for each boundary. */
	*Sens_BPress,     /*!< \brief Back pressure sensitivity coefficient for each boundary. */
	**CSensitivity;		/*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
	su2double Total_Sens_Mach;	/*!< \brief Total mach sensitivity coefficient for all the boundaries. */
	su2double Total_Sens_AoA;		/*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
	su2double Total_Sens_Geo;		/*!< \brief Total shape sensitivity coefficient for all the boundaries. */
	su2double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
	su2double Total_Sens_Temp;    /*!< \brief Total farfield sensitivity to temperature. */
	su2double Total_Sens_BPress;    /*!< \brief Total sensitivity to back pressure. */
	su2double *iPoint_UndLapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*jPoint_UndLapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	bool space_centered;  /*!< \brief True if space centered scheeme used. */
    su2double **Jacobian_Axisymmetric; /*!< \brief Storage for axisymmetric Jacobian. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
	su2double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	su2double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
  su2double *FlowPrimVar_i,	/*!< \brief Store the flow solution at point i. */
	*FlowPrimVar_j;        /*!< \brief Store the flow solution at point j. */

  su2double pnorm,
  Area_Monitored; /*!< \brief Store the total area of the monitored outflow surface (used for normalization in continuous adjoint outflow conditions) */
    
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
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iMesh, unsigned long Iteration);

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
	void GetEngine_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
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
	su2double GetPsiRho_Inf(void);
    
	/*!
	 * \brief Compute the adjoint energy at the infinity.
	 * \return Value of the adjoint energy at the infinity.
	 */
	su2double GetPsiE_Inf(void);
    
	/*!
	 * \brief Compute Phi (adjoint velocity) at the infinity.
	 * \param[in] val_dim - Index of the adjoint velocity vector.
	 * \return Value of the adjoint velocity vector at the infinity.
	 */
	su2double GetPhi_Inf(unsigned short val_dim);
    
	/*!
	 * \brief Compute the spatial integration using a centered scheme for the adjoint equations.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the interface adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the near-field adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the adjoint symmetry boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    

  /*!
   * \brief Impose the supersonic inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
	void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
	              unsigned short val_marker);
  
  /*!
   * \brief Impose the supersonic outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                           unsigned short val_marker);

	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Impose the engine inflow adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose the engine bleed adjoint boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Bleed(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                        CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Impose the engine exhaust boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);
    
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Compute the inviscid sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Smooth the inviscid sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Get the shape sensitivity coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the sensitivity coefficient.
	 */
	su2double GetCSensitivity(unsigned short val_marker, unsigned long val_vertex);
    
	/*!
	 * \brief Set the shape sensitivity coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \param[in] val_sensitivity - Value of the sensitivity coefficient.
	 */
	void SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity);
    
	/*!
	 * \brief Provide the total shape sensitivity coefficient.
	 * \return Value of the geometrical sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	su2double GetTotal_Sens_Geo(void);
    
	/*!
	 * \brief Set the total Mach number sensitivity coefficient.
	 * \return Value of the Mach sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	su2double GetTotal_Sens_Mach(void);
    
	/*!
	 * \brief Set the total angle of attack sensitivity coefficient.
	 * \return Value of the angle of attack sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	su2double GetTotal_Sens_AoA(void);
    
	/*!
	 * \brief Set the total farfield pressure sensitivity coefficient.
	 * \return Value of the farfield pressure sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	su2double GetTotal_Sens_Press(void);
    
	/*!
	 * \brief Set the total farfield temperature sensitivity coefficient.
	 * \return Value of the farfield temperature sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	su2double GetTotal_Sens_Temp(void);

	/*!
	 * \author H. Kline
	 * \brief Set the total Back pressure number sensitivity coefficient.
	 * \return Value of the Back sensitivity coefficient
	 *         (inviscid + viscous contribution).
	 */
	su2double GetTotal_Sens_BPress(void);
  
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
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
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
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iMesh, unsigned long Iteration);

    
	/*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition (heat flux).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
    /*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition (heat flux).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Compute the viscous sensitivity of the functional.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);
    
	/*!
	 * \brief Compute the viscous residuals for the adjoint equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
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
 * \author F. Palacios, A. Bueno.
 * \version 4.0.2 "Cardinal"
 */
class CAdjTurbSolver : public CSolver {
private:
	su2double PsiNu_Inf,	/*!< \brief PsiNu variable at the infinity. */
	*FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j;	/*!< \brief Store the flow solution at point j. */

	su2double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	su2double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
        
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
	CAdjTurbSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
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
	 * \brief Default destructor of the class.
	 */
	virtual ~CAdjTurbSolver(void);
    
	/*!
	 * \brief Impose the Navier-Stokes turbulent adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
  
  /*!
	 * \brief Impose an isothermal wall boundary condition (no-slip).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
  
	/*!
	 * \brief Impose the boundary condition to the far field using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Compute the viscous residuals for the turbulent adjoint equation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
};

/*! \class CPoissonSolver
 *  \brief Main class for defining the poisson potential solver.
 *  \author F. Palacios
 *  \version 4.0.2 "Cardinal"
 *  \date May 3, 2010.
 */
class CPoissonSolver : public CSolver {
private:
	su2double *Source_Vector;		  /*!< \brief Auxiliary vector for storing element source vector. */
  su2double **StiffMatrix_Elem; /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	su2double **StiffMatrix_Node;	/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
    
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CPoissonSolver(void);
    
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CPoissonSolver(CGeometry *geometry, CConfig *config);
    
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
	~CPoissonSolver(void);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Integrate the Poisson equation using a Galerkin method.
	 * \param[in] StiffMatrix_Elem - Element stiffness matrix
	 */
	void AddStiffMatrix(su2double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3);

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
	 * \brief Impose via the residual the Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
  void BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose via the residual the Neumann boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
  void BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
  
};

/*! \class CWaveSolver
 *  \brief Main class for defining the wave solver.
 *  \author F. Palacios
 *  \version 4.0.2 "Cardinal"
 *  \date May 3, 2010.
 */
class CWaveSolver : public CSolver {
private:
	su2double *CWave;	/*!< \brief Wave strength for each boundary. */
	su2double AllBound_CWave;	/*!< \brief Total wave strength for all the boundaries. */
	su2double Total_CWave; /*!< \brief Total wave strength for all the boundaries. */
    
    CSysMatrix StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	CSysMatrix StiffMatrixTime;	/*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
    
    su2double **StiffMatrix_Elem,			/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iRKStep);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose a Dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
	 * \brief Load a solution from a restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
	 */
	void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter);
    
	/*!
	 * \brief Compute the total wave strength coefficient.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Wave_Strength(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Build stiffness matrix in space.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSpace_Matrix(CGeometry *geometry,
                         CConfig   *config);

	/*!
	 * \brief Provide the total wave strength.
	 * \return Value of the wave strength.
	 */
	su2double GetTotal_CWave(void);
    
};

/*! \class CHeatSolver
 *  \brief Main class for defining the heat solver.
 *  \author F. Palacios
 *  \version 4.0.2 "Cardinal"
 *  \date May 3, 2010.
 */
class CHeatSolver : public CSolver {
private:
	su2double *CHeat;	     /*!< \brief Heat strength for each boundary. */
	su2double Total_CHeat; /*!< \brief Total Heat strength for all the boundaries. */
    
  CSysMatrix StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	CSysMatrix StiffMatrixTime;	 /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
    
  su2double **StiffMatrix_Elem; /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	su2double **StiffMatrix_Node;	 /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
    
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iRKStep);
  
  /*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition (heat flux).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                        unsigned short val_marker);
  
  /*!
	 * \brief Impose via the residual or brute force the Navier-Stokes adjoint boundary condition (heat flux).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
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
	 * \brief Provide the total heat strength.
	 * \return Value of the heat strength.
	 */
	su2double GetTotal_CHeat(void);
    
};

/*! \class CFEASolver
 *  \brief Main class for defining the FEA solver.
 *  \author F. Palacios, R. Sanchez.
 *  \version 4.0.2 "Cardinal"
 *  \date May 3, 2010.
 */
class CFEASolver : public CSolver {
private:
  
	su2double  Total_CFEA;			/*!< \brief Total FEA coefficient for all the boundaries. */
    CSysMatrix StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
	CSysMatrix StiffMatrixTime;	/*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */

    CSysMatrix MassMatrix; 		/*!< \brief Sparse structure for storing the mass matrix in Galerkin computations. */
	CSysMatrix DampMatrix;	/*!< \brief Sparse structure for storing the damping matrix in Galerkin computations. */

	CSysVector TimeRes_Aux;				/*!< \brief Auxiliary vector for adding mass and damping contributions to the residual. */
	CSysVector TimeRes;					/*!< \brief Vector for adding mass and damping contributions to the residual */
  
  su2double **StiffMatrix_Elem,			/*!< \brief Auxiliary matrices for storing elem to elem Stiffness Matrices. */
	**StiffMatrix_Node,					/*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
	**MassMatrix_Elem,					/*!< \brief Auxiliary matrices for storing elem to elem Mass Matrices. */
	**MassMatrix_Node,					/*!< \brief Auxiliary matrices for storing point to point Mass Matrices. */
	**MassMatrix_Node_Int,				/*!< \brief Auxiliary matrices for storing point to point Mass Matrices * a0. */
	**DampMatrix_Elem,					/*!< \brief Auxiliary matrices for storing elem to elem Damping Matrices. */
	**DampMatrix_Node,					/*!< \brief Auxiliary matrices for storing point to point Damping Matrices. */
	*DeadLoadVector_Elem,				/*!< \brief Auxiliary vector for storing point to point Dead Loads. */
	*DeadLoadVector_Node;				/*!< \brief Auxiliary vector for storing point to point Dead Loads. */

  su2double a_dt[8];			/*!< \brief Integration constants. */

  su2double WAitken_Dyn;			/*!< \brief Aitken's dynamic coefficient */
  su2double WAitken_Dyn_tn1;		/*!< \brief Aitken's dynamic coefficient in the previous iteration */

  su2double FSI_Conv[2];		/*!< \brief Values to check the convergence of the FSI problem. */



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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief Impose a displacement (constraint) boundary condition --> Clamped boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);
	/*!
	 * \brief Impose a displacement (constraint) boundary condition --> Clamped boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);
    
	/*!
	 * \brief Impose a displacement (constraint) boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker);
    
	/*!
	 * \brief Impose a load boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                     unsigned short val_marker);
    
	/*!
	 * \brief Impose a load boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                 unsigned short val_marker);
  
	/*!
	 * \brief Impose a load boundary condition in cartesian coordinates.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                 unsigned short val_marker);

	/*!
	 * \brief Impose a sine-wave load boundary condition in cartesian coordinates.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                 unsigned short val_marker);


  /*!
	 * \brief Impose a load boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Set residuals to zero.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
			unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Update the solution using an implicit solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);

	/*!
	 * \brief Update the solution using an implicit Newmark solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
    
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
	 * \brief Get the surface pressure from a file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void GetSurface_Pressure(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Set the the pressure load in the FEA solver.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] flow_solution - Container vector with all the solutions.
	 * \param[in] fea_config - Definition of the particular problem.
	 */
	void SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config, CNumerics *fea_numerics);
    
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
	su2double GetTotal_CFEA(void);
    
	/*!
	 * \brief Set the value of the FEA coefficient.
	 * \param[in] val_cfea - Value of the FEA coefficient.
	 */
	void SetTotal_CFEA(su2double val_cfea);

	/*!
	 * \brief Set the displacement for the nodes in the structural mesh
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_grid_movement - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] flow_geometry - Definition of the particular problem.
	 */
	void SetStruct_Displacement(CGeometry **fea_geometry,
                                CConfig *fea_config,
                                CSolver ***fea_solution);

	/*!
	 * \brief Predictor for structural displacements based on previous iterations
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_grid_movement - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] flow_geometry - Definition of the particular problem.
	 */
	void PredictStruct_Displacement(CGeometry **fea_geometry,
                                	CConfig *fea_config,
                                	CSolver ***fea_solution);

	/*!
	 * \brief Computation of Aitken's coefficient.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	void ComputeAitken_Coefficient(CGeometry **fea_geometry,
            				  CConfig *fea_config,
            				  CSolver ***fea_solution,
            				  unsigned long iFSIIter);

	/*!
	 * \brief Aitken's relaxation of the solution.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	void SetAitken_Relaxation(CGeometry **fea_geometry,
            				  CConfig *fea_config,
            				  CSolver ***fea_solution);

	/*!
	 * \brief Aitken's relaxation of the solution.
	 * \param[in] fea_geometry - Geometrical definition of the problem.
	 * \param[in] fea_config - Geometrical definition of the problem.
	 * \param[in] fea_geometry - Definition of the particular problem.
	 */
	void Update_StructSolution(CGeometry **fea_geometry,
            				  CConfig *fea_config,
            				  CSolver ***fea_solution);

	/*!
	 * \brief Get the value of the FSI convergence.
	 * \param[in] Set value of interest: 0 - Initial value, 1 - Current value.
	 */
	void SetFSI_ConvValue(unsigned short val_index, su2double val_criteria);

	/*!
	 * \brief Get the value of the FSI convergence.
	 * \param[in]  Value of interest: 0 - Initial value, 1 - Current value.
	 * \return Values to compare
	 */
	su2double GetFSI_ConvValue(unsigned short val_index);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Compute_StiffMassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Compute_StiffMassDampMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] solver - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Compute_IntegrationConstants(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config);

	/*!
	 * \brief Set the solution variables at time n to the current solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	void SetSolution_time_n(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Retrieve the value of the dynamic Aitken relaxation factor.
	 * \return Value of the dynamic Aitken relaxation factor.
	 */
	su2double GetWAitken_Dyn(void);

	/*!
	 * \brief Retrieve the value of the last Aitken relaxation factor in the previous time step.
	 * \return Value of the last Aitken relaxation factor in the previous time step.
	 */
	su2double GetWAitken_Dyn_tn1(void);

	/*!
	 * \brief Set the value of the dynamic Aitken relaxation factor
	 * \param[in] Value of the dynamic Aitken relaxation factor
	 */
	void SetWAitken_Dyn(su2double waitk);

	/*!
	 * \brief Set the value of the last Aitken relaxation factor in the current time step.
	 * \param[in] Value of the last Aitken relaxation factor in the current time step.
	 */
	void SetWAitken_Dyn_tn1(su2double waitk_tn1);

    
};

/*!
 * \class CAdjLevelSetSolver
 * \brief Main class for defining the level set solver.
 * \ingroup LevelSet_Model
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
 */
class CAdjLevelSetSolver : public CSolver {
protected:
	su2double *FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
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
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
	/*!
	 * \brief Compute the spatial integration using a upwind scheme.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term computation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the Navier-Stokes wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the Far Field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
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
 * \author F. Palacios
 * \version 4.0.2 "Cardinal"
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
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
	 * \param[in] numerics - Description of the numerical method.
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Source term integration.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
    
	/*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                          unsigned short val_marker);
    
	/*!
	 * \brief Impose the far-field boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                  unsigned short val_marker);
    
	/*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                   unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry plane boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method.
	 * \param[in] visc_numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                      unsigned short val_marker);
    
	/*!
	 * \brief Impose the dirichlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
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
 * \class CDiscAdjSolver
 * \brief Main class for defining the discrete adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring
 * \version 4.0.2 "Cardinal"
 */
class CDiscAdjSolver : public CSolver {
private:
  unsigned short KindDirect_Solver;
  CSolver *direct_solver;
  su2double *Sens_Mach, /*!< \brief Mach sensitivity coefficient for each boundary. */
  *Sens_AoA,			/*!< \brief Angle of attack sensitivity coefficient for each boundary. */
  *Sens_Geo,			/*!< \brief Shape sensitivity coefficient for each boundary. */
  *Sens_Press,			/*!< \brief Pressure sensitivity coefficient for each boundary. */
  *Sens_Temp,			/*!< \brief Temperature sensitivity coefficient for each boundary. */
  **CSensitivity;	/*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
  su2double Total_Sens_Mach;	/*!< \brief Total mach sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_AoA;		/*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Geo;		/*!< \brief Total shape sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
  su2double Total_Sens_Temp;    /*!< \brief Total farfield sensitivity to temperature. */
  su2double Total_Sens_BPress;    /*!< \brief Total sensitivity to outlet pressure. */
  su2double ObjFunc_Value;        /*!< \brief Value of the objective function. */
  su2double Mach, Alpha, Beta, Pressure, Temperature;
  unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CDiscAdjSolver(CGeometry *geometry, CConfig *config);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   */
  CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);

  /*!
   * \brief Performs the preprocessing of the adjoint AD-based solver.
   *        Registers all necessary variables on the tape. Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterSolution(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Performs the preprocessing of the adjoint AD-based solver.
   *        Registers all necessary variables that are output variables on the tape.
   *        Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterOutput(CGeometry *geometry, CConfig *config);

  /*!
  * \brief Sets the adjoint values of the output of the flow (+turb.) iteration
  *         before evaluation of the tape.
  * \param[in] geometry - The geometrical definition of the problem.
  * \param[in] config - The particular config.
  */
  void SetAdjoint_Output(CGeometry *geometry, CConfig *config);

  /*!
  * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
  *        after tape has been evaluated.
  * \param[in] geometry - The geometrical definition of the problem.
  * \param[in] config - The particular config.
  */
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config);

  /*!
  * \brief Register the objective function as output.
  * \param[in] geometry - The geometrical definition of the problem.
  */
  void RegisterObj_Func(CConfig *config);

  /*!
   * \brief Set the surface sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSurface_Sensitivity(CGeometry *geometry, CConfig* config);

  /*!
   * \brief Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the objective function.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetAdj_ObjFunc(CGeometry *geometry, CConfig* config);


  /*!
   * \brief Provide the total shape sensitivity coefficient.
   * \return Value of the geometrical sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_Geo(void);

  /*!
   * \brief Set the total Mach number sensitivity coefficient.
   * \return Value of the Mach sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_Mach(void);

  /*!
   * \brief Set the total angle of attack sensitivity coefficient.
   * \return Value of the angle of attack sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_AoA(void);

  /*!
   * \brief Set the total farfield pressure sensitivity coefficient.
   * \return Value of the farfield pressure sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_Press(void);

  /*!
   * \brief Set the total farfield temperature sensitivity coefficient.
   * \return Value of the farfield temperature sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_Temp(void);

  /*!
   * \author H. Kline
   * \brief Set the total Back pressure number sensitivity coefficient.
   * \return Value of the Back sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_BPress(void);

  /*!
   * \brief Get the shape sensitivity coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the sensitivity coefficient.
   */
  su2double GetCSensitivity(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetRecording(CGeometry *geometry, CConfig *config, unsigned short kind_recording);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void RegisterVariables(CGeometry *geometry, CConfig *config, bool reset = false);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config);
};
#include "solver_structure.inl"
