/*!
 * \file solver_structure.hpp
 * \brief Headers of the main subroutines for solving partial differential equations.
 *        The subroutines and functions are in the <i>solver_structure.cpp</i>,
 *        <i>solution_direct.cpp</i>, <i>solution_adjoint.cpp</i>, and
 *        <i>solution_linearized.cpp</i> files.
 * \author F. Palacios, T. Economon
 * \version 3.2.8.2 "eagle"
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
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "fluid_model.hpp"
#include "numerics_structure.hpp"
#include "variable_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "numerics_machine_learning.hpp"

using namespace std;

/*!
 * \class CSolver
 * \brief Main class for defining the PDE solution, it requires
 * a child class for each particular solver (Euler, Navier-Stokes, etc.)
 * \author F. Palacios
 * \version 3.2.8.2 "eagle"
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
	double Max_Delta_Time,	/*!< \brief Maximum value of the delta time for all the control volumes. */
	Min_Delta_Time;					/*!< \brief Minimum value of the delta time for all the control volumes. */
	double *Residual_RMS,	/*!< \brief Vector with the mean residual for each variable. */
  *Residual_Max,        /*!< \brief Vector with the maximal residual for each variable. */
	*Residual,						/*!< \brief Auxiliary nVar vector. */
	*Residual_i,					/*!< \brief Auxiliary nVar vector for storing the residual at point i. */
	*Residual_j;					/*!< \brief Auxiliary nVar vector for storing the residual at point j. */
  unsigned long *Point_Max; /*!< \brief Vector with the maximal residual for each variable. */
  double **Point_Max_Coord; /*!< \brief Vector with pointers to the coords of the maximal residual for each variable. */
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
   * \param[in] val_point - Value of the point index for the max residual.
   * \param[in] val_coord - Location (x,y,z) of the max residual point.
	 */
	void AddRes_Max(unsigned short val_var, double val_residual, unsigned long val_point, double* val_coord);
    
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
   * \brief Get the location of the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the location (x,y,z) of the biggest residual for the variable <i>val_var</i>.
   */
  double* GetPoint_Max_Coord(unsigned short val_var);
  
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
	virtual void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker);
  
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
	virtual void BC_IsothermalCatalytic_Wall(CGeometry *geometry,
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
	virtual void BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
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
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics *conv_numerics,
                                         CNumerics *visc_numerics,
                                         CConfig *config,
                                         unsigned short val_marker);
  
  /*!
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	virtual void BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics *conv_numerics,
                                            CNumerics *visc_numerics,
                                            CConfig *config,
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
	virtual void SetPreconditioner(CConfig *config, unsigned short iPoint);
    
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
	virtual void AddStiffMatrix(double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3 );
    
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
	 * \param[in] val_Total_Heat - Value of the total heat load.
	 */
	virtual void SetTotal_HeatFlux(double val_Total_Heat);
    
    /*!
	 * \brief A virtual member.
	 * \param[in] val_Total_MaxHeat - Value of the total heat load.
	 */
	virtual void SetTotal_MaxHeatFlux(double val_Total_MaxHeat);
    
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
	 * \return Value of the z moment coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCMz_Inv(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCMz_Visc(unsigned short val_marker);
    
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CLift(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CDrag(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CSideForce(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CEff(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CFx(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CFy(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CFz(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CMx(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CMy(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CMz(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CLift_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CDrag_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CSideForce_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CEff_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CFx_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CFy_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CFz_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CMx_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CMy_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  virtual double GetSurface_CMz_Inv(unsigned short val_marker);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	virtual double GetCSideForce_Visc(unsigned short val_marker);
    
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
	virtual double GetInflow_MassFlow(unsigned short val_marker);
    
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
	virtual double GetInflow_Pressure(unsigned short val_marker);
    
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the fan face mach on the surface <i>val_marker</i>.
	 */
	virtual double GetInflow_Mach(unsigned short val_marker);
    
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
	virtual double GetTotal_HeatFlux(void);
    
    /*!
	 * \brief A virtual member.
	 * \return Value of the heat load (integrated heat flux).
	 */
	virtual double GetTotal_MaxHeatFlux(void);
    
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
	 * \return Value of the difference of the presure and the target pressure.
	 */
	virtual double GetTotal_CpDiff(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the difference of the heat and the target heat.
	 */
	virtual double GetTotal_HeatFluxDiff(void);
  
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
	 * \param[in] val_pressure - Value of the difference between pressure and the target pressure.
	 */
	virtual void SetTotal_CpDiff(double val_pressure);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_pressure - Value of the difference between heat and the target heat.
	 */
	virtual void SetTotal_HeatFluxDiff(double val_heat);
  
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
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual double GetAllBound_CMx_Inv(void);
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual double GetAllBound_CMy_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual double GetAllBound_CMz_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual double GetAllBound_CFx_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual double GetAllBound_CFy_Inv(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual double GetAllBound_CFz_Inv(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	virtual double GetAllBound_CLift_Visc(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	virtual double GetAllBound_CSideForce_Visc(void);
  
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
	 * \return Value of the pressure coefficient.
	 */
	virtual double GetCPressureTarget(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual void SetCPressureTarget(unsigned short val_marker, unsigned short val_vertex, double val_pressure);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual double *GetCharacPrimVar(unsigned short val_marker, unsigned short val_vertex);
    
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
	virtual double GetHeatFlux(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	virtual double GetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	virtual void SetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex, double val_heat);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the y plus.
	 */
	virtual double GetYPlus(unsigned short val_marker, unsigned short val_vertex);

  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  virtual double GetStrainMag_Max(void);

  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  virtual double GetOmega_Max(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  virtual void SetStrainMag_Max(double val_strainmag_max);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  virtual void SetOmega_Max(double val_omega_max);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the adjoint density at the infinity.
	 */
	virtual double GetPsiRho_Inf(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the adjoint density at the infinity.
	 */
	virtual double* GetPsiRhos_Inf(void);
  
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
	 * \return Value of the velocity at the infinity.
	 */
	virtual double *GetVelocity_Inf(void);
    
	/*!
	 * \brief A virtual member.
	 * \return Value of the viscosity at the infinity.
	 */
	virtual double GetViscosity_Inf(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the turbulent kinetic energy.
	 */
	virtual double GetTke_Inf(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the sensitivity coefficient.
	 */
	virtual double GetCSensitivity(unsigned short val_marker, unsigned short val_vertex);
    
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
	virtual double* GetConstants();

  /*!
   * \brief A virtual member.
   * \return average total pressure evaluated at an exit boundary marker
   */
  virtual double GetOneD_TotalPress(void);

  /*!
   * \brief A virtual member.
   * \param[in] val_exit_pt: value of the total average pressure at the exit.
   */
  virtual void SetOneD_TotalPress(double AveragePressure);

  /*!
   * \brief A virtual member.
   *\return average Mach number evaluated at an exit boundary marker
   */
  virtual double GetOneD_Mach(void);

  /*!
   * \brief A virtual member.
   * set average Mach number evaluated at an exit boundary marker
   */
  virtual void SetOneD_Mach(double AverageMach);
  
  /*!
   * \brief A virtual member.
   *\return average temperature evaluated at an exit boundary marker
   */
  virtual double GetOneD_Temp(void);
  
  /*!
   * \brief A virtual member.
   * set average temperature evaluated at an exit boundary marker
   */
  virtual void SetOneD_Temp(double AverageTemperature);
  
  /*!
   * \brief A virtual member.
   * \return average temperature evaluated at an exit boundary marker
   */
  virtual double GetOneD_MassFlowRate(void);
  
  /*!
   * \brief A virtual member.
   * set average temperature evaluated at an exit boundary marker
   */
  virtual void SetOneD_MassFlowRate(double MassFlowRate);
  
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged pressure at a marker.(same as area averaged pressure)
   */
  virtual double GetOneD_FluxAvgPress(void);
  
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged pressure at a marker. (same as area averaged pressure)
   */
  virtual void SetOneD_FluxAvgPress(double PressureRef);
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged density at a marker. (\f$ = (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) \f$)
   */
  virtual double GetOneD_FluxAvgDensity(void);
  
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged density at a marker.( \f$= (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) \f$)
   */
  virtual void SetOneD_FluxAvgDensity(double DensityRef);
  
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged velocity at a marker. = \f$ \sqrt ( \frac{\int((rho*u)*u^2dA)}{\int(rho*u*dA) }) \f$
   */
  virtual double GetOneD_FluxAvgVelocity(void);
  
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged velocity at a marker. = \f$ \sqrt (  \frac{\int((rho*u)*u^2dA)}{\int(rho*u*dA) }) \f$
   */
  virtual void SetOneD_FluxAvgVelocity(double VelocityRef);
  
  /*!
   * \brief A virtual member.
   * \ Get the flux averaged enthalpy at a marker. =\f$ \frac{ \int(rho*u*h dA) }{ \int(rho *u *dA )} \f$
   */
  virtual double GetOneD_FluxAvgEntalpy(void);
  /*!
   * \brief A virtual member.
   * \ Set the flux averaged enthalpy at a marker. =\f$ \frac{ \int(rho*u*h dA) }{ \int(rho *u *dA ) }\f$
   */
  virtual void SetOneD_FluxAvgEntalpy(double EnthalpyRef);
  
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
                           CConfig *flow_config);
    
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
	void Gauss_Elimination(double** A, double* rhs, unsigned short nVar);
    
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
  void SetUpTypicalSectionWingModel(double (&PHI)[2][2],double (&lambda)[2], CConfig *config);
    
  /*!
  * \brief Solve the typical section wing model.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] Cl - Coefficient of lift at particular iteration.
  * \param[in] Cm - Moment coefficient about z-axis at particular iteration.
	* \param[in] config - Definition of the particular problem.
  * \param[in] val_Marker - Surface that is being monitored.
  * \param[in] displacements - solution of typical section wing model.
	*/
  void SolveTypicalSectionWingModel(CGeometry *geometry, double Cl, double Cm, CConfig *config, unsigned short val_Marker, double (&displacements)[4]);

};

/*!
 * \class CBaselineSolver
 * \brief Main class for defining a baseline solution from a restart file (for output).
 * \author F. Palacios, T. Economon.
 * \version 3.2.8.2 "eagle"
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
 * \version 3.2.8.2 "eagle"
 */
class CEulerSolver : public CSolver {
protected:
	
  double
  Mach_Inf,	/*!< \brief Mach number at the infinity. */
	Density_Inf,	/*!< \brief Density at the infinity. */
	Energy_Inf,			/*!< \brief Energy at the infinity. */
  Temperature_Inf,			/*!< \brief Energy at the infinity. */
	Pressure_Inf,		/*!< \brief Pressure at the infinity. */
	*Velocity_Inf;		/*!< \brief Flow Velocity vector at the infinity. */
	
  double
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
	
  double
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
	
  double
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
  double *Surface_CLift,   /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CDrag,          /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSideForce,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CEff,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,            /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,            /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,            /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,            /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,            /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz;            /*!< \brief z Moment coefficient for each monitoring surface. */
	double *iPoint_UndLapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*jPoint_UndLapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	double *SecondaryVar_i,	/*!< \brief Auxiliary vector for storing the solution at point i. */
	*SecondaryVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	double *PrimVar_i,	/*!< \brief Auxiliary vector for storing the solution at point i. */
	*PrimVar_j;			/*!< \brief Auxiliary vector for storing the solution at point j. */
	double **LowMach_Precontioner; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
	bool space_centered,  /*!< \brief True if space centered scheeme used. */
	euler_implicit,			/*!< \brief True if euler implicit scheme used. */
	least_squares;        /*!< \brief True if computing gradients by least squares. */
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
  
  double *Primitive,		/*!< \brief Auxiliary nPrimVar vector. */
	*Primitive_i,				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
	*Primitive_j;				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */
  
  double *Secondary,		/*!< \brief Auxiliary nPrimVar vector. */
	*Secondary_i,				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
	*Secondary_j;				/*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  double Cauchy_Value,	/*!< \brief Summed value of the convergence indicator. */
	Cauchy_Func;			/*!< \brief Current value of the convergence indicator at one iteration. */
	unsigned short Cauchy_Counter;	/*!< \brief Number of elements of the Cauchy serial. */
	double *Cauchy_Serie;			/*!< \brief Complete Cauchy serial. */
	double Old_Func,	/*!< \brief Old value of the objective function (the function which is monitored). */
	New_Func;			/*!< \brief Current value of the objective function (the function which is monitored). */
  double AoA_old;  /*!< \brief Old value of the angle of attack (monitored). */

  CFluidModel  *FluidModel;  /*!< \brief fluid model used in the solver */

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
	double GetDensity_Inf(void);
    
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
	 * \brief Compute the pressure at the infinity.
	 * \return Value of the pressure at the infinity.
	 */
	double GetPressure_Inf(void);

	/*!
	 * \brief Compute the density multiply by velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	double GetDensity_Velocity_Inf(unsigned short val_dim);
    
	/*!
	 * \brief Get the velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	double GetVelocity_Inf(unsigned short val_dim);
  
  /*!
	 * \brief Get the velocity at the infinity.
	 * \return Value of the velocity at the infinity.
	 */
	double *GetVelocity_Inf(void);
  
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
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
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
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
  
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
	 * \brief Provide the non dimensional lift coefficient (inviscid contribution).
	 * \param val_marker Surface where the coefficient is going to be computed.
	 * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCLift_Inv(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional z moment coefficient (inviscid contribution).
	 * \param val_marker Surface where the coefficient is going to be computed.
	 * \return Value of the z moment coefficient (inviscid contribution) on the surface <i>val_marker</i>.
	 */
	double GetCMz_Inv(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional lift coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the lift coefficient on the surface <i>val_marker</i>.
	 */
	double GetSurface_CLift(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional drag coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the drag coefficient on the surface <i>val_marker</i>.
	 */
	double GetSurface_CDrag(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CSideForce(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CEff(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CFx(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CFy(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CFz(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional x moment coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
	 */
	double GetSurface_CMx(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional y moment coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
	 */
	double GetSurface_CMy(unsigned short val_marker);
    
    /*!
	 * \brief Provide the non dimensional z moment coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
	 */
	double GetSurface_CMz(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CLift_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CDrag_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CSideForce_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CEff_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CFx_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CFy_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CFz_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CMx_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CMy_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  double GetSurface_CMz_Inv(unsigned short val_marker);
  
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
	double GetInflow_MassFlow(unsigned short val_marker);
    
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
	double GetInflow_Pressure(unsigned short val_marker);
    
	/*!
	 * \brief Provide the mass flow rate.
	 * \param val_marker Surface where the coeficient is going to be computed.
	 * \return Value of the fan face mach on the surface <i>val_marker</i>.
	 */
	double GetInflow_Mach(unsigned short val_marker);
    
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
	 * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_CpDiff(void);
  
  /*!
	 * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
	 * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
	 */
	double GetTotal_HeatFluxDiff(void);
    
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
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_CpDiff(double val_pressure);
  
  /*!
	 * \brief Set the value of the Equivalent Area coefficient.
	 * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
	 */
	void SetTotal_HeatFluxDiff(double val_heat);
    
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
	double GetTotal_HeatFlux(void);
    
    /*!
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	double GetTotal_MaxHeatFlux(void);
    
	/*!
	 * \brief Store the total (inviscid + viscous) non dimensional torque coefficient.
	 * \param[in] val_Total_CQ - Value of the total torque coefficient.
	 */
	void SetTotal_CQ(double val_Total_CQ);
    
    /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Heat - Value of the heat load.
	 */
	void SetTotal_HeatFlux(double val_Total_Heat);
    
    /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Heat - Value of the heat load.
	 */
	void SetTotal_MaxHeatFlux(double val_Total_MaxHeat);
    
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
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CMx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CMy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CMz_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CFx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CFy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CFz_Inv(void);
  
	/*!
	 * \brief Provide the Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	double GetCPressure(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
	 * \brief Provide the Target Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	double GetCPressureTarget(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
	 * \brief Set the value of the target Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
  void SetCPressureTarget(unsigned short val_marker, unsigned short val_vertex, double val_pressure);

  /*!
	 * \brief Value of the characteristic variables at the boundaries.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
	double *GetCharacPrimVar(unsigned short val_marker, unsigned short val_vertex);
  
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
   * \brief Provide the averaged total pressure at a marker.
   */
	double GetOneD_TotalPress(void);
  
	/*!
   * \brief Set the value of averaged total pressure
   * \param[in] val_exit_pt - value of the averaged pressure
   */
	void SetOneD_TotalPress(double AveragePressure);
  
  /*!
   * \brief Provide the averaged Mach number at a marker.
   */
  double GetOneD_Mach(void);
  
  /*!
   * \brief Set the averaged Mach number at a marker.
   */
  void SetOneD_Mach(double AverageMach);
  
  /*!
   * \brief Provide the averaged Mach number at a marker.
   */
  double GetOneD_Temp(void);
  
  /*!
   * \brief Set the averaged Temperature at a marker.
   */
  void SetOneD_Temp(double AverageTemperature);
  
  /*!
   * \brief Provide the averaged Mach number at a marker.
   */
  double GetOneD_MassFlowRate(void);
  
  /*!
   * \brief Set the averaged Temperature at a marker.
   */
  void SetOneD_MassFlowRate(double MassFlowRate);
  
  /*!
   * \brief Get the flux averaged pressure at a marker.(same as area averaged pressure)
   */
  double GetOneD_FluxAvgPress(void);
  
  /*!
   * \brief Set the flux averaged pressure at a marker. (same as area averaged pressure)
   */
  void SetOneD_FluxAvgPress(double PressureRef);
  
  /*!
   * \brief Get the flux averaged density at a marker. ( = (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) )
   */
  double GetOneD_FluxAvgDensity(void);
  
  /*!
   * \brief Set the flux averaged density at a marker.( = (gamma/(gamma-1)) / ( Pref*(href-1/2 uref^2) )
   */
  void SetOneD_FluxAvgDensity(double DensityRef);
  
  /*!
   * \brief Get the flux averaged velocity at a marker. = \f$ \sqrt ( \int((rho*u)*u^2dA)/\int(rho*u*dA) )\f$
   */
  double GetOneD_FluxAvgVelocity(void);
  
  /*!
   * \brief Set the flux averaged velocity at a marker. =\f$ sqrt ( \int((rho*u)*u^2dA)/\int(rho*u*dA) ) \f$
   */
  void SetOneD_FluxAvgVelocity(double VelocityRef);
  
  /*!
   * \brief Get the flux averaged enthalpy at a marker. = \f$ \int(rho*u*h dA) / \int(rho *u *dA ) \f$
   */
  double GetOneD_FluxAvgEntalpy(void);
  
  /*!
   * \brief Set the flux averaged enthalpy at a marker. =\f$ \int(rho*u*h dA) / \int(rho *u *dA ) \f$
   */
  void SetOneD_FluxAvgEntalpy(double EnthalpyRef);
  
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
  
};

/*!
 * \class CNSSolver
 * \brief Main class for defining the Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 * \version 3.2.8.2 "eagle"
 */
class CNSSolver : public CEulerSolver {
private:
	double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */
  double Tke_Inf;	/*!< \brief Turbulent kinetic energy at the infinity. */
	double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
	Prandtl_Turb;         /*!< \brief Turbulent Prandtl number. */
	double *CDrag_Visc,	/*!< \brief Drag coefficient (viscous contribution) for each boundary. */
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
	double *ForceViscous,	/*!< \brief Viscous force for each boundary. */
	*MomentViscous;			/*!< \brief Inviscid moment for each boundary. */
	double AllBound_CDrag_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
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
  double StrainMag_Max, Omega_Max; /*!< \brief Maximum Strain Rate magnitude and Omega. */
  
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
	 * \brief Get the turbulent kinetic energy at the infinity.
	 * \return Value of the turbulent kinetic energy at the infinity.
	 */
	double GetTke_Inf(void);
    
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
	double GetCLift_Visc(unsigned short val_marker);
    
    /*!
	 * \brief Get the non dimensional z moment coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the z moment coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	double GetCMz_Visc(unsigned short val_marker);
  
  /*!
	 * \brief Get the non dimensional sideforce coefficient (viscous contribution).
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
	 */
	double GetCSideForce_Visc(unsigned short val_marker);
    
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
	 * \brief Get the total non dimensional sideforce coefficient (viscous contribution).
	 * \return Value of the lift coefficient (viscous contribution).
	 */
	double GetAllBound_CSideForce_Visc(void);
    
	/*!
	 * \brief Get the total non dimensional drag coefficient (viscous contribution).
	 * \return Value of the drag coefficient (viscous contribution).
	 */
	double GetAllBound_CDrag_Visc(void);
    
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
	double GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex);
    
	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetHeatFlux(unsigned short val_marker, unsigned short val_vertex);
	
  /*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
	 * \brief Set the value of the target Pressure coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the pressure coefficient.
	 */
  void SetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex, double val_heat);
  
	/*!
	 * \brief Get the y plus.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the y plus.
	 */
	double GetYPlus(unsigned short val_marker, unsigned short val_vertex);
  
  /*!
   * \brief Get the max Omega.
   * \return Value of the max Omega.
   */
  double GetOmega_Max(void);
  
  /*!
   * \brief Get the max Strain rate magnitude.
   * \return Value of the max Strain rate magnitude.
   */
  double GetStrainMag_Max(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  void SetStrainMag_Max(double val_strainmag_max);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  void SetOmega_Max(double val_omega_max);
  
};

/*!
 * \class CTurbSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 3.2.8.2 "eagle"
 */
class CTurbSolver : public CSolver {
protected:
	double *FlowPrimVar_i,  /*!< \brief Store the flow solution at point i. */
	*FlowPrimVar_j,         /*!< \brief Store the flow solution at point j. */
	*lowerlimit,            /*!< \brief contains lower limits for turbulence variables. */
	*upperlimit;            /*!< \brief contains upper limits for turbulence variables. */
	double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */
    
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
 * \version 3.2.8.2 "eagle"
 */

class CTurbSASolver: public CTurbSolver {
private:
	double nu_tilde_Inf, nu_tilde_Engine;
	
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
                               CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the near-field boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                               CConfig *config, unsigned short val_marker);
  
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
 * \class CTurbMLSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author B. Tracey.
 * \version 3.2.8.2 "eagle"
 */

class CTurbMLSolver: public CTurbSolver {
private:
	double nu_tilde_Inf;
	
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbMLSolver(void);
  
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbMLSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbMLSolver(void);
  
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
	 * \brief Load a solution from a restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iter - Current external iteration number.
	 */
	void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter);
  
  double SAProduction, SADestruction, SACrossProduction, SASource, MLProduction, MLDestruction, MLCrossProduction, MLSource, SourceDiff;
  
};

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Aranake.
 * \version 3.2.8.2 "eagle"
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
 * \version 3.2.8.2 "eagle"
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
	double* GetConstants();
    
};

/*!
 * \class CAdjEulerSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 3.2.8.2 "eagle"
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
	double *iPoint_UndLapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*jPoint_UndLapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	bool space_centered;  /*!< \brief True if space centered scheeme used. */
    double **Jacobian_Axisymmetric; /*!< \brief Storage for axisymmetric Jacobian. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
  double *FlowPrimVar_i,	/*!< \brief Store the flow solution at point i. */
	*FlowPrimVar_j;        /*!< \brief Store the flow solution at point j. */

  double pnorm;
    
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
	void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                               unsigned short val_marker);
    
	/*!
	 * \brief Impose via the residual the near-field adjoint boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                               unsigned short val_marker);
    
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
 * \version 3.2.8.2 "eagle"
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
 * \version 3.2.8.2 "eagle"
 */
class CAdjTurbSolver : public CSolver {
private:
	double PsiNu_Inf,	/*!< \brief PsiNu variable at the infinity. */
	*FlowSolution_i,	/*!< \brief Store the flow solution at point i. */
	*FlowSolution_j;	/*!< \brief Store the flow solution at point j. */

	double Gamma;									/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;				/*!< \brief Fluids's Gamma - 1.0  . */
        
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

/*!
 * \class CLinEulerSolver
 * \brief Main class for defining the linearized Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 3.2.8.2 "eagle"
 */
class CLinEulerSolver : public CSolver {
private:
	double DeltaRho_Inf,	/*!< \brief Linearized density variable at the infinity. */
	DeltaE_Inf,				/*!< \brief Linearized energy at the infinity. */
	*DeltaVel_Inf;			/*!< \brief Linearized velocity vector at the infinity. */
	double *iPoint_UndLapl,	/*!< \brief Undivided Laplacians for centered scheme. */
	*jPoint_UndLapl;			/*!< \brief Undivided Laplacians for centered scheme. */
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
	CLinEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CLinEulerSolver(void);
    
	/*!
	 * \brief Compute the spatial integration using a centered scheme for the linearized equations.
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
	 * \brief Compute the undivided laplacian for the linearized solution.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
    
	/*!
	 * \brief Impose via the residual the linearized Euler wall boundary condition.
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
    
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

/*! \class CPoissonSolver
 *  \brief Main class for defining the poisson potential solver.
 *  \author F. Palacios
 *  \version 3.2.8.2 "eagle"
 *  \date May 3, 2010.
 */
class CPoissonSolver : public CSolver {
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
	void AddStiffMatrix(double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3);

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
 *  \version 3.2.8.2 "eagle"
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
	double GetTotal_CWave(void);
    
};

/*! \class CHeatSolver
 *  \brief Main class for defining the heat solver.
 *  \author F. Palacios
 *  \version 3.2.8.2 "eagle"
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
	double GetTotal_CHeat(void);
    
};

/*! \class CFEASolver
 *  \brief Main class for defining the FEA solver.
 *  \author F. Palacios
 *  \version 3.2.8.2 "eagle"
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
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iRKStep);
    
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
 * \class CAdjLevelSetSolver
 * \brief Main class for defining the level set solver.
 * \ingroup LevelSet_Model
 * \author F. Palacios
 * \version 3.2.8.2 "eagle"
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
 * \version 3.2.8.2 "eagle"
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
 * \class CTNE2EulerSolver
 * \brief Main class for defining the TNE2 Euler's flow solver.
 * \ingroup Euler_Equations
 * \author S. R. Copeland, F. Palacios
 * \version 2.0.6
 */
class CTNE2EulerSolver : public CSolver {
protected:

  unsigned short
  nSpecies;	                /*!< \brief Number of species in the gas mixture. */
  
  double
  Gamma,                    /*!< \brief Mixture Cp/Cv. */
	Gamma_Minus_One;	        /*!< \brief Mixture Cp/Cv - 1. */
	
  double
  Mach_Inf,       	        /*!< \brief Free stream Mach number. */
  *Density,                 /*!< \brief Free stream species density. */
  Energy_ve_Inf,            /*!< \brief Vib.-el. free stream energy. */
	Pressure_Inf,		          /*!< \brief Free stream pressure. */
	*Velocity_Inf,		        /*!< \brief Free stream flow velocity. */
  *MassFrac_Inf,            /*!< \brief Free stream species mass fraction. */
  Temperature_Inf,          /*!< \brief Trans.-rot. free stream temperature. */
  Temperature_ve_Inf;       /*!< \brief Vib.-el. free stream temperature. */

  double
  *lowerlimit,            /*!< \brief contains lower limits for conserved variables. */
	*upperlimit;            /*!< \brief contains upper limits for conserved variables. */

	double                    
  *CDrag_Inv,	              /*!< \brief Boundary invisc. Cd contribution. */
	*CLift_Inv,			          /*!< \brief Boundary invisc. Cl contribution. */
	*CSideForce_Inv,		      /*!< \brief Boundary invisc. sideforce contribution. */
	*CMx_Inv,			            /*!< \brief X-moment contribution on boundaries. */
	*CMy_Inv,			            /*!< \brief Y-moment contribution on boundaries. */
	*CMz_Inv,			            /*!< \brief Z-moment contribution on boundaries. */
	*CFx_Inv,			            /*!< \brief X-force contribution on boundaries. */
	*CFy_Inv,                 /*!< \brief Y-force contribution on boundaries. */
	*CFz_Inv,                 /*!< \brief Z-force contribution on boundaries. */
	*CEff_Inv,			        	/*!< \brief Invisc. Cl/Cd on boundaries. */
  *ForceInviscid,		        /*!< \brief Inviscid forces at domain boundaries. */
	*MomentInviscid,	        /*!< \brief Inviscid moments at domain boundaries. */
	**CPressure,		          /*!< \brief Cp at each node on domain boundaries . */
	**HeatFlux,		      /*!< \brief Cq at each node on domain boundaries . */
	**CharacPrimVar,		      /*!< \brief Value of the characteristic variables at the boundary . */
  AllBound_CDrag_Inv,	      /*!< \brief Sum of CDrag_Inv from all boundaries. */
	AllBound_CLift_Inv,			  /*!< \brief Sum of CLift_Inv from all boundaries. */
	AllBound_CSideForce_Inv,	/*!< \brief Sum of CSideForce_Inv from all boundaries. */
	AllBound_CMx_Inv,			    /*!< \brief Sum of CMx_Inv from all boundaries. */
	AllBound_CMy_Inv,			    /*!< \brief Sum of CMy_Inv from all boundaries. */
	AllBound_CMz_Inv,			    /*!< \brief Sum of CMz_Inv from all boundaries. */
	AllBound_CFx_Inv,			    /*!< \brief Sum of CFx_Inv from all boundaries. */
	AllBound_CFy_Inv,			    /*!< \brief Sum of CFy_Inv from all boundaries. */
	AllBound_CFz_Inv,			    /*!< \brief Sum of CFz_Inv from all boundaries. */
	AllBound_CEff_Inv;        /*!< \brief Sum of CEff_Inv from all boundaries. */
  
	double
  Total_CDrag,              /*!< \brief Total Cd. */
	Total_CLift,		          /*!< \brief Total Cl. */
	Total_CSideForce,		      /*!< \brief Total CSideForce. */
	Total_CMx,			          /*!< \brief Total CMx. */
	Total_CMy,                /*!< \brief Total CMy. */
	Total_CMz,                /*!< \brief Total CMz. */
	Total_CFx,                /*!< \brief Total CFx. */
	Total_CFy,                /*!< \brief Total CFy. */
	Total_CFz,                /*!< \brief Total CFz. */
	Total_CEff,               /*!< \brief Total CEff. */
  Total_Heat,                  /*!< \brief Total heat load. */
  Total_MaxHeat;               /*!< \brief Maximum heat flux on all boundaries. */

	double
  *PrimVar_i,	      /*!< \brief Vector for storing primitives at node i. */
	*PrimVar_j;			  /*!< \brief Vector for storing primitives at node j. */

	double
  **LowMach_Precontioner; /*!< \brief Matrix for storing the inverse of preconditioner. */
	
  unsigned long
  nMarker;				  /*!< \brief Total number of domain boundaries. */
  
	bool least_squares;     /*!< \brief Indicator for least-squares computed grads. */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CTNE2EulerSolver(void);
  
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTNE2EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CTNE2EulerSolver(void);
  
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
  void Set_MPI_Primitive(CGeometry *geometry, CConfig *config);
  
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
  
  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config);

  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config);

  /*!
	 * \brief Impose the send-receive boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config);
  
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
  void Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Set the fluid solver nondimensionalization.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
	 * \brief Set the maximum value of the eigenvalue.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Compute the density at the infinity.
	 * \return Value of the density at the infinity.
	 */
	double GetDensity_Inf(void);
  
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
	 * \brief Compute the pressure at the infinity.
	 * \return Value of the pressure at the infinity.
	 */
	double GetPressure_Inf(void);
  
	/*!
	 * \brief Compute the density multiply by velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the density multiply by the velocity at the infinity.
	 */
	double GetDensity_Velocity_Inf(unsigned short val_dim);
  
	/*!
	 * \brief Get the velocity at the infinity.
	 * \param[in] val_dim - Index of the velocity vector.
	 * \return Value of the velocity at the infinity.
	 */
	double GetVelocity_Inf(unsigned short val_dim);
  
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
   * \param[in] second_numerics - Description of the second numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                       CNumerics *second_numerics, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
	 */
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
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
  
  /*!
	 * \brief Compute the gradient of the primitive variables using a Least-Squares method,
	 *        and stores the result in the <i>Gradient_Primitive</i> variable.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, unsigned long val_Point);
  
	/*!
	 * \brief Compute the limiter of the primitive variables.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Compute slope limiter.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSolution_Limiter(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
	 * \param[in] iPoint - Index of the grid point
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetPreconditioner(CConfig *config, unsigned short iPoint);
   
	/*!
	 * \brief Impose via the residual the Euler wall boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] numerics - Description of the numerical method.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                     CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Impose the far-field boundary condition using characteristics.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Impose the symmetry boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Impose a subsonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
	/*!
	 * \brief Impose a supersonic inlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                           CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker);

  /*!
	 * \brief Impose the outlet boundary condition.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
   
	 */
	void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
	 * \brief Update the solution using an explicit Euler scheme.
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
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	double GetTotal_HeatFlux(void);
  
  /*!
	 * \brief Provide the total heat load.
	 * \return Value of the heat load (viscous contribution).
	 */
	double GetTotal_MaxHeatFlux(void);
  
  /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_Heat - Value of the heat load.
	 */
	void SetTotal_HeatFlux(double val_Total_Heat);
  
  /*!
	 * \brief Store the total heat load.
	 * \param[in] val_Total_HeatMax - Value of the heat load.
	 */
	void SetTotal_MaxHeatFlux(double val_Total_MaxHeat);
  
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
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CMx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CMy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CMz_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CFx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CFy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  double GetAllBound_CFz_Inv(void);
  
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
	 * \brief Load a direct flow solution for use with the adjoint solver.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_iZone - Current zone in the mesh.
	 */
	void GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone);
  
  /*!
	 * \brief Load the output data container with the variables to be written to the volume solution file.
   * \param[in] config - Definition of the particular problem.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] data_container - Container holding the output variable data.
   * \param[in] nOutput_Vars - Number of output variables being stored.
	 */
	void SetVolume_Output(CConfig *config, CGeometry *geometry, double **data_container, unsigned short nOutput_Vars);
  
};

/*!
 * \class CTNE2NSSolver
 * \brief Main class for defining the TNE2 Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author S. R. Copeland, F. Palacios
 * \version 2.0.6
 */
class CTNE2NSSolver : public CTNE2EulerSolver {
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
  *Heat_Visc,		/*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHeatFlux_Visc, /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
  
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
  AllBound_HeatFlux_Visc,		/*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHeatFlux_Visc; /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CTNE2NSSolver(void);
  
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTNE2NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CTNE2NSSolver(void);
  
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
	void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
	 * \brief Impose the symmetry boundary condition using the residual.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose a constant heat-flux condition at the wall.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_IsothermalCatalytic_Wall(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CNumerics *conv_numerics,
                                   CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker);
  
  /*!
	 * \brief Impose the Navier-Stokes boundary condition (strong).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_marker - Surface marker where the boundary condition is applied.
	 */
	void BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *conv_numerics,
                                      CNumerics *visc_numerics,
                                      CConfig *config,
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
	double GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex);
  
	/*!
	 * \brief Get the skin friction coefficient.
	 * \param[in] val_marker - Surface marker where the coefficient is computed.
	 * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
	 * \return Value of the heat transfer coefficient.
	 */
	double GetHeatFlux(unsigned short val_marker, unsigned short val_vertex);

};


/*!
 * \class CAdjEulerSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 * \version 2.0.6
 */
class CAdjTNE2EulerSolver : public CSolver {
protected:
  unsigned short
  nSpecies;
  
	double
  *PsiRho_Inf,	/*!< \brief Free-stream adjoint density. */
	*Phi_Inf,			/*!< \brief Phi vector at the infinity. */
  PsiE_Inf,			/*!< \brief PsiE variable at the infinity. */
  PsiEve_Inf;
  
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
	double *iPoint_UndLapl,	/*!< \brief Auxiliary variable for the undivided Laplacians. */
	*jPoint_UndLapl;			/*!< \brief Auxiliary variable for the undivided Laplacians. */
	bool space_centered;  /*!< \brief True if space centered scheeme used. */
	unsigned long nMarker;				/*!< \brief Total number of markers using the grid information. */
  double Gamma, Gamma_Minus_One;
  double **Jacobian_Axisymmetric;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjTNE2EulerSolver(void);
  
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CAdjTNE2EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CAdjTNE2EulerSolver(void);
  
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
	double* GetPsiRhos_Inf(void);
  
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
	 * \brief Update the solution using an explicit solver.
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
 * \version 2.0.6
 */
class CAdjTNE2NSSolver : public CAdjTNE2EulerSolver {
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CAdjTNE2NSSolver(void);
  
	/*!
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CAdjTNE2NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAdjTNE2NSSolver(void);
  
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


#include "solver_structure.inl"
