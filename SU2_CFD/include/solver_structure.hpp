/*!
 * \file solver_structure.hpp
 * \brief Headers of the main subroutines for solving partial differential equations.
 *        The subroutines and functions are in the <i>solver_structure.cpp</i>,
 *        <i>solution_direct.cpp</i>, <i>solution_adjoint.cpp</i>, and
 *        <i>solution_linearized.cpp</i> files.
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

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <set>
#include <stdlib.h>
#include <stdio.h>

#include "fluid_model.hpp"
#include "task_definition.hpp"
#include "numerics_structure.hpp"
#include "sgs_model.hpp"
#include "variables/CVariable.hpp"
#include "../../Common/include/gauss_structure.hpp"
#include "../../Common/include/element_structure.hpp"
#include "../../Common/include/fem_geometry_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/blas_structure.hpp"
#include "../../Common/include/graph_coloring_structure.hpp"
#include "../../Common/include/toolboxes/MMS/CVerificationSolution.hpp"

using namespace std;

/*!
 * \class CSolver
 * \brief Main class for defining the PDE solution, it requires
 * a child class for each particular solver (Euler, Navier-Stokes, etc.)
 * \author F. Palacios
 */
class CSolver {
protected:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
  unsigned short MGLevel;        /*!< \brief Multigrid level of this solver object. */
  unsigned short IterLinSolver;  /*!< \brief Linear solver iterations. */
  unsigned short nVar,          /*!< \brief Number of variables of the problem. */
  nPrimVar,                     /*!< \brief Number of primitive variables of the problem. */
  nPrimVarGrad,                 /*!< \brief Number of primitive variables of the problem in the gradient computation. */
  nSecondaryVar,                     /*!< \brief Number of primitive variables of the problem. */
  nSecondaryVarGrad,                 /*!< \brief Number of primitive variables of the problem in the gradient computation. */
  nVarGrad,                 /*!< \brief Number of variables for deallocating the LS Cvector. */
  nDim;                          /*!< \brief Number of dimensions of the problem. */
  unsigned long nPoint;          /*!< \brief Number of points of the computational grid. */
  unsigned long nPointDomain;   /*!< \brief Number of points of the computational grid. */
  su2double Max_Delta_Time,  /*!< \brief Maximum value of the delta time for all the control volumes. */
  Min_Delta_Time;          /*!< \brief Minimum value of the delta time for all the control volumes. */
  su2double *Residual_RMS,  /*!< \brief Vector with the mean residual for each variable. */
  *Residual_Max,        /*!< \brief Vector with the maximal residual for each variable. */
  *Residual,            /*!< \brief Auxiliary nVar vector. */
  *Residual_i,          /*!< \brief Auxiliary nVar vector for storing the residual at point i. */
  *Residual_j;          /*!< \brief Auxiliary nVar vector for storing the residual at point j. */
  su2double *Residual_BGS,  /*!< \brief Vector with the mean residual for each variable for BGS subiterations. */
  *Residual_Max_BGS;        /*!< \brief Vector with the maximal residual for each variable for BGS subiterations. */
  unsigned long *Point_Max; /*!< \brief Vector with the maximal residual for each variable. */
  unsigned long *Point_Max_BGS; /*!< \brief Vector with the maximal residual for each variable. */
  su2double **Point_Max_Coord; /*!< \brief Vector with pointers to the coords of the maximal residual for each variable. */
  su2double **Point_Max_Coord_BGS; /*!< \brief Vector with pointers to the coords of the maximal residual for each variable. */
  su2double *Solution,    /*!< \brief Auxiliary nVar vector. */
  *Solution_i,        /*!< \brief Auxiliary nVar vector for storing the solution at point i. */
  *Solution_j;        /*!< \brief Auxiliary nVar vector for storing the solution at point j. */
  su2double *Vector,  /*!< \brief Auxiliary nDim vector. */
  *Vector_i,      /*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point i. */
  *Vector_j;      /*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point j. */
  su2double *Res_Conv,  /*!< \brief Auxiliary nVar vector for storing the convective residual. */
  *Res_Visc,        /*!< \brief Auxiliary nVar vector for storing the viscous residual. */
  *Res_Sour,        /*!< \brief Auxiliary nVar vector for storing the viscous residual. */
  *Res_Conv_i,      /*!< \brief Auxiliary vector for storing the convective residual at point i. */
  *Res_Visc_i,      /*!< \brief Auxiliary vector for storing the viscous residual at point i. */
  *Res_Conv_j,      /*!< \brief Auxiliary vector for storing the convective residual at point j. */
  *Res_Visc_j;      /*!< \brief Auxiliary vector for storing the viscous residual at point j. */
  su2double **Jacobian_i,  /*!< \brief Auxiliary matrices for storing point to point Jacobians at point i. */
  **Jacobian_j;          /*!< \brief Auxiliary matrices for storing point to point Jacobians at point j. */
  su2double **Jacobian_ii,  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  **Jacobian_ij,        /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  **Jacobian_ji,        /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  **Jacobian_jj;        /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  su2double *iPoint_UndLapl,  /*!< \brief Auxiliary variable for the undivided Laplacians. */
  *jPoint_UndLapl;      /*!< \brief Auxiliary variable for the undivided Laplacians. */
  su2double **Smatrix,  /*!< \brief Auxiliary structure for computing gradients by least-squares */
  **Cvector;       /*!< \brief Auxiliary structure for computing gradients by least-squares */

  int *Restart_Vars;       /*!< \brief Auxiliary structure for holding the number of variables and points in a restart. */
  int Restart_ExtIter;     /*!< \brief Auxiliary structure for holding the external iteration offset from a restart. */
  passivedouble *Restart_Data; /*!< \brief Auxiliary structure for holding the data values from a restart. */
  unsigned short nOutputVariables;  /*!< \brief Number of variables to write. */

  unsigned long nMarker_InletFile;       /*!< \brief Auxiliary structure for holding the number of markers in an inlet profile file. */
  vector<string> Marker_Tags_InletFile;       /*!< \brief Auxiliary structure for holding the string names of the markers in an inlet profile file. */
  unsigned long *nRow_InletFile;       /*!< \brief Auxiliary structure for holding the number of rows for a particular marker in an inlet profile file. */
  unsigned long *nRowCum_InletFile;       /*!< \brief Auxiliary structure for holding the number of rows in cumulative storage format for a particular marker in an inlet profile file. */
  unsigned long maxCol_InletFile;       /*!< \brief Auxiliary structure for holding the maximum number of columns in all inlet marker profiles (for data structure size) */
  unsigned long *nCol_InletFile;       /*!< \brief Auxiliary structure for holding the number of columns for a particular marker in an inlet profile file. */
  passivedouble *Inlet_Data; /*!< \brief Auxiliary structure for holding the data values from an inlet profile file. */

  bool rotate_periodic;    /*!< \brief Flag that controls whether the periodic solution needs to be rotated for the solver. */
  bool implicit_periodic;  /*!< \brief Flag that controls whether the implicit system should be treated by the periodic BC comms. */
  
public:
  
  CSysVector<su2double> LinSysSol;    /*!< \brief vector to store iterative solution of implicit linear system. */
  CSysVector<su2double> LinSysRes;    /*!< \brief vector to store iterative residual of implicit linear system. */
  CSysVector<su2double> LinSysAux;    /*!< \brief vector to store iterative residual of implicit linear system. */
#ifndef CODI_FORWARD_TYPE
  CSysMatrix<passivedouble> Jacobian; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
  CSysSolve<passivedouble>  System;   /*!< \brief Linear solver/smoother. */
#else
  CSysMatrix<su2double> Jacobian;
  CSysSolve<su2double>  System;
#endif
  
  CSysMatrix<su2double> StiffMatrix; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations, and grid movement. */
  
  CSysVector<su2double> OutputVariables;    /*!< \brief vector to store the extra variables to be written. */
  string* OutputHeadingNames; /*< \brief vector of strings to store the headings for the exra variables */
  
  CVariable** node;  /*!< \brief Vector which the define the variables for each problem. */
  CVariable* node_infty; /*!< \brief CVariable storing the free stream conditions. */
  
  CVerificationSolution *VerificationSolution; /*!< \brief Verification solution class used within the solver. */

  /*!
   * \brief Constructor of the class.
   */
  CSolver(void);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSolver(void);
  
  /*!
   * \brief Routine to load a solver quantity into the data structures for MPI point-to-point communication and to launch non-blocking sends and recvs.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  void InitiateComms(CGeometry *geometry,
                     CConfig *config,
                     unsigned short commType);
  
  /*!
   * \brief Routine to complete the set of non-blocking communications launched by InitiateComms() and unpacking of the data in the solver class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  void CompleteComms(CGeometry *geometry,
                     CConfig *config,
                     unsigned short commType);
  
  /*!
   * \brief Routine to load a solver quantity into the data structures for MPI periodic communication and to launch non-blocking sends and recvs.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] val_periodic_index - Index for the periodic marker to be treated (first in a pair).
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  void InitiatePeriodicComms(CGeometry *geometry,
                             CConfig *config,
                             unsigned short val_periodic_index,
                             unsigned short commType);
  
  /*!
   * \brief Routine to complete the set of non-blocking periodic communications launched by InitiatePeriodicComms() and unpacking of the data in the solver class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] val_periodic_index - Index for the periodic marker to be treated (first in a pair).
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  void CompletePeriodicComms(CGeometry *geometry,
                             CConfig *config,
                             unsigned short val_periodic_index,
                             unsigned short commType);

  /*!
   * \brief Set number of linear solver iterations.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void SetIterLinSolver(unsigned short val_iterlinsolver);
  
  /*!
   * \brief Set the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void SetResidual_RMS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Communicate the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void SetResidual_BGS(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Set the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  virtual void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  virtual void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  virtual void SetNondimensionalization(CConfig *config, unsigned short iMesh);

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
   * \brief Adds the maximal residual, this is useful for the convergence history (overload).
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   * \param[in] val_point - Value of the point index for the max residual.
   * \param[in] val_coord - Location (x, y, z) of the max residual point.
   */
  void AddRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point, const su2double* val_coord);

  /*!
   * \brief Get the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  su2double GetRes_Max(unsigned short val_var);
  
  /*!
   * \brief Set the residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  void SetRes_BGS(unsigned short val_var, su2double val_residual);
  
  /*!
   * \brief Adds the residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  void AddRes_BGS(unsigned short val_var, su2double val_residual);
  
  /*!
   * \brief Get the residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  su2double GetRes_BGS(unsigned short val_var);
  
  /*!
   * \brief Set the maximal residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  void SetRes_Max_BGS(unsigned short val_var, su2double val_residual, unsigned long val_point);
  
  /*!
   * \brief Adds the maximal residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   * \param[in] val_point - Value of the point index for the max residual.
   * \param[in] val_coord - Location (x, y, z) of the max residual point.
   */
  void AddRes_Max_BGS(unsigned short val_var, su2double val_residual, unsigned long val_point, su2double* val_coord);
  
  /*!
   * \brief Get the maximal residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  su2double GetRes_Max_BGS(unsigned short val_var);
  
  /*!
   * \brief Get the residual for FEM structural analysis.
   * \param[in] val_var - Index of the variable.
   * \return Value of the residual for the variable in the position <i>val_var</i>.
   */
  virtual su2double GetRes_FEM(unsigned short val_var);
  
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
   * \brief Get the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  unsigned long GetPoint_Max_BGS(unsigned short val_var);
  
  /*!
   * \brief Get the location of the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the location (x, y, z) of the biggest residual for the variable <i>val_var</i>.
   */
  su2double* GetPoint_Max_Coord_BGS(unsigned short val_var);
  
  /*!
   * \brief Set Value of the residual due to the Geometric Conservation Law (GCL) for steady rotating frame problems.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRotatingFrame_GCL(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute the Green-Gauss gradient of the auxiliary variable.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SetAuxVar_Gradient_GG(CGeometry *geometry, CConfig *config);
  
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
  
  /*!
   * \brief Compute the pressure laplacian using in a incompressible solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] PressureLaplacian - Pressure laplacian.
   */
  void SetPressureLaplacian(CGeometry *geometry, CConfig *config, su2double *PressureLaplacian);
  
  /*!
   * \brief Set the old solution variables to the current solution value for Runge-Kutta iteration.
            It is a virtual function, because for the DG-FEM solver a different version is needed.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  virtual void Set_OldSolution(CGeometry *geometry);

  /*!
   * \brief Set the new solution variables to the current solution value for classical RK.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  virtual void Set_NewSolution(CGeometry *geometry);

  /*!
   * \brief Load the geometries at the previous time states n and nM1.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void Restart_OldGeometry(CGeometry *geometry, CConfig *config);
  
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
   * \param[in]     config          - Definition of the particular problem.
   * \param[in]     TimeSync        - The synchronization time.
   * \param[in,out] timeEvolved     - On input the time evolved before the time step,
                                      on output the time evolved after the time step.
   * \param[out]    syncTimeReached - Whether or not the synchronization time is reached.
   */
  virtual void CheckTimeSynchronization(CConfig         *config,
                                        const su2double TimeSync,
                                        su2double       &timeEvolved,
                                        bool            &syncTimeReached);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  virtual void ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                  CNumerics **numerics, CConfig *config,
                                  unsigned short iMesh);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  virtual void ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                         CNumerics **numerics, CConfig *config,
                                         unsigned short iMesh, unsigned short RunTime_EqSystem);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  virtual void ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                      CNumerics **numerics, CConfig *config,
                                      unsigned short iMesh, unsigned short RunTime_EqSystem);

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
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  virtual void Convective_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                           CConfig *config, unsigned short iMesh, unsigned short iRKStep);

	/*!
	 * \brief A virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] Output - boolean to determine whether to print output.
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
  virtual void Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] config - Definition of the particular problem.
    */
  virtual void Set_MPI_Nearfield(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] config - Definition of the particular problem.
    */
  virtual void Set_MPI_Interface(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Set_Heatflux_Areas(CGeometry *geometry, CConfig *config);
  
  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Evaluate_ObjFunc(CConfig *config);
  
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
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  
  virtual void BC_DispDir(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
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
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Damper(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                 unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                           CNumerics *numerics, CConfig *config);
  
  /*!
  * \brief Impose the interface state across sliding meshes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] solver_container - Container vector with all the solutions.
  * \param[in] conv_numerics - Description of the numerical method.
  * \param[in] visc_numerics - Description of the numerical method.
  * \param[in] config - Definition of the particular problem.
  */
  virtual void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                          CConfig *config, unsigned short val_marker, bool val_inlet_surface);
  
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
  virtual void BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container,
	                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief It computes Fourier transformation for the needed quantities along the pitch for each span in turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  virtual void PreprocessBC_Giles(CGeometry *geometry, CConfig *config, CNumerics *conv_numerics, unsigned short marker_flag);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Giles(CGeometry *geometry, CSolver **solver_container,
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
  virtual void BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
  virtual void BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
   * \param[in] geometry         - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics    - Description of the convective numerical method.
   * \param[in] visc_numerics    - Description of the viscous numerical method.
   * \param[in] config           - Definition of the particular problem.
   * \param[in] val_marker       - Surface marker where the boundary condition is applied.
   */
  virtual void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
  virtual void BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
    
	/*!
	 * \brief Impose the symmetry boundary condition using the residual.
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
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  virtual void BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker);
   
 /*!
   * \brief Get the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   */
  virtual su2double GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index);

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  virtual void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   * \param[in] component  - set value
   */
  virtual void SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component);

  /*!
   * \brief Get the number of outer states for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  virtual int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Set the number of outer states for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value      - number of outer states
   */
  virtual void SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   * \param[in] relaxation factor - relaxation factor for the change of the variables
   * \param[in] val_var           - value of the variable
   */
  virtual void SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   */
  virtual su2double GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var);

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
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  virtual void ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
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
   */
  virtual void ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
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
  virtual void Pressure_Forces(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Momentum_Forces(CGeometry *geometry, CConfig *config);
  
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
  virtual void Friction_Forces(CGeometry *geometry, CConfig *config);
    
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Buffet_Monitoring(CGeometry *geometry, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Heat_Fluxes(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
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
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetPrimitive_Limiter_MPI(CGeometry *geometry, CConfig *config);
  
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
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_CD(su2double val_Total_CD);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  virtual void SetTotal_CL(su2double val_Total_CL);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_NetThrust(su2double val_Total_NetThrust);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_Power(su2double val_Total_Power);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_SolidCD(su2double val_Total_SolidCD);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_ReverseFlow(su2double val_ReverseFlow);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_MFR(su2double val_Total_MFR);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_Prop_Eff(su2double val_Total_Prop_Eff);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_Poly_Eff(su2double val_Total_Poly_Eff);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_IDC(su2double val_Total_IDC);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_IDC_Mach(su2double val_Total_IDC_Mach);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_IDR(su2double val_Total_IDR);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  virtual void SetTotal_DC60(su2double val_Total_DC60);

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  virtual void SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  virtual void AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight);
  
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
  virtual su2double GetCL_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetCL_Visc(unsigned short val_marker);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CL(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CD(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CSF(unsigned short val_marker);
  
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
  virtual su2double GetSurface_CL_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CD_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CSF_Inv(unsigned short val_marker);
  
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
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CL_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CD_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CSF_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CEff_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFx_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFy_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFz_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMx_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMy_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMz_Visc(unsigned short val_marker);
    
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the buffet metric on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_Buffet_Metric(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CL_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CD_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CSF_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CEff_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFx_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFy_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CFz_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMx_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMy_Mnt(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_CMz_Mnt(unsigned short val_marker);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetCSF_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetCD_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  virtual su2double GetInflow_MassFlow(unsigned short val_marker);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  virtual void GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
   * \brief A virtual member.
   */
  virtual void GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  virtual void GetEllipticSpanLoad_Diff(CGeometry *geometry, CConfig *config);

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
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  virtual void SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output);

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
  virtual su2double GetCSF_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetCEff_Inv(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the integrated heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_HF_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the maximum heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetSurface_MaxHF_Visc(unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  virtual su2double GetCD_Visc(unsigned short val_marker);
  
  /*!
   * \author H. Kline
   * \brief Set the total "combo" objective (weighted sum of other values).
   * \param[in] ComboObj - Value of the combined objective.
   */
  virtual void SetTotal_ComboObj(su2double ComboObj);
  
  /*!
   * \author H. Kline
   * \brief Provide the total "combo" objective (weighted sum of other values).
   * \return Value of the "combo" objective values.
   */
  virtual su2double GetTotal_ComboObj(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the sideforce coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CSF(void);
  
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
   * \brief A virtual member.
   * \return Value of the average temperature.
   */
  virtual su2double GetTotal_AvgTemperature(void);
  
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
   * \return Value of the Aero drag (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_AeroCD(void);

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
   * \return Value of the FEA coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CFEA(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Near-Field Pressure coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CNearFieldOF(void);
  
  /*!
   * \author H. Kline
   * \brief Add to the value of the total 'combo' objective.
   * \param[in] val_obj - Value of the contribution to the 'combo' objective.
   */
  virtual void AddTotal_ComboObj(su2double val_obj);
  
  /*!
   * \brief A virtual member.
   * \return Value of the objective function for a reference geometry.
   */
  virtual su2double GetTotal_OFRefGeom(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the objective function for a reference node.
   */
  virtual su2double GetTotal_OFRefNode(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the objective function for the volume fraction.
   */
  virtual su2double GetTotal_OFVolFrac(void);

  /*!
   * \brief A virtual member.
   * \return Bool that defines whether the solution has an element-based file or not
   */
  virtual bool IsElementBased(void);

  /*!
   * \brief A virtual member.
   * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
   */
  virtual void SetTotal_CEquivArea(su2double val_cequivarea);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_aerocd - Value of the aero drag.
   */
  virtual void SetTotal_AeroCD(su2double val_aerocd);

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
   * \param[in] val_cfea - Value of the FEA coefficient.
   */
  virtual void SetTotal_CFEA(su2double val_cfea);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference geometry.
   */
  virtual void SetTotal_OFRefGeom(su2double val_ofrefgeom);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference node.
   */
  virtual void SetTotal_OFRefNode(su2double val_ofrefnode);

  /*!
   * \brief A virtual member.
   * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
   */
  virtual void SetTotal_CNearFieldOF(su2double val_cnearfieldpress);
  
  /*!
   * \brief A virtual member.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CL(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CD(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_NetThrust(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_Power(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_SolidCD(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_ReverseFlow(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_MFR(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_Prop_Eff(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_ByPassProp_Eff(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_Adiab_Eff(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_Poly_Eff(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_IDC(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_IDC_Mach(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_IDR(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_DC60(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the custom objective function.
   */
  virtual su2double GetTotal_Custom_ObjFunc(void);

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
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CoPx(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CoPy(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_CoPz(void);

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
  virtual su2double GetAllBound_CL_Inv(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CD_Inv(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CSF_Inv(void);
  
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
  virtual su2double GetAllBound_CoPx_Inv(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPy_Inv(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPz_Inv(void);

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
   * \return Value of the lift coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CL_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CD_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CSF_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CEff_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMx_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMy_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMz_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPx_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPy_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPz_Visc(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFx_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFy_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFz_Visc(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CL_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CD_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CSF_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CEff_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMx_Mnt(void);
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMy_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CMz_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPx_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPy_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CoPz_Mnt(void);

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFx_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFy_Mnt(void);
  
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  virtual su2double GetAllBound_CFz_Mnt(void);
    
  /*!
   * \brief A virtual member.
   * \return Value of the buffet metric.
   */
  virtual su2double GetTotal_Buffet_Metric(void);
  
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
  
  virtual void SetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual su2double *GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual void SetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  
  virtual void SetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual su2double GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  
  virtual su2double *GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  
  virtual su2double GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual unsigned long GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual void SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index);
  
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
   * \return Value of the pressure coefficient.
   */
  virtual su2double GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual void SetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex, su2double val_deltap);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual su2double GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  virtual void SetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex, su2double val_deltat);
  
  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total temperature is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is evaluated.
   * \return Value of the total temperature
   */
  virtual su2double GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total pressure is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is evaluated.
   * \return Value of the total pressure
   */
  virtual su2double GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the flow direction is evaluated
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is evaluated
   * \param[in] val_dim - The component of the flow direction unit vector to be evaluated
   * \return Component of a unit vector representing the flow direction.
   */
  virtual su2double GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim);
  
  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total temperature is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is set.
   * \param[in] val_ttotal - Value of the total temperature
   */
  virtual void SetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ttotal);
  
  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total pressure is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is set.
   * \param[in] val_ptotal - Value of the total pressure
   */
  virtual void SetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ptotal);
  
  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the flow direction is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is set.
   * \param[in] val_dim - The component of the flow direction unit vector to be set
   * \param[in] val_flowdir - Component of a unit vector representing the flow direction.
   */
  virtual void SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir);

  /*!
   * \brief A virtual member
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the turbulence variable (i.e. k is 0 in SST)
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  virtual void SetInlet_TurbVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_turb_var);

  /*!
   * \brief A virtual member
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  virtual void SetUniformInlet(CConfig* config, unsigned short iMarker);

  /*!
   * \brief A virtual member
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  virtual void SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief A virtual member
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  virtual su2double GetInletAtVertex(su2double *val_inlet,
                                     unsigned long val_inlet_point,
                                     unsigned short val_kind_marker,
                                     string val_marker,
                                     CGeometry *geometry,
                                     CConfig *config);

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  virtual void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the skin friction coefficient.
   */
  virtual su2double GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim);
  
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
   * \return Value of the buffet sensor.
   */
  virtual su2double GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex);
  
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
   * \brief Get the total back pressure sensitivity coefficient.
   * \return Value of the back pressure sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  virtual su2double GetTotal_Sens_BPress(void);

  /*!
   * \brief A virtual member.
   * \return Value of the density sensitivity.
   */
  virtual su2double GetTotal_Sens_Density(void);

  /*!
   * \brief A virtual member.
   * \return Value of the velocity magnitude sensitivity.
   */
  virtual su2double GetTotal_Sens_ModVel(void);

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
   * \return Value of nu tilde at the far-field.
   */
  virtual su2double GetNuTilde_Inf(void);

  /*!
   * \brief A virtual member.
   * \return Value of the turbulent kinetic energy.
   */
  virtual su2double GetTke_Inf(void);

  /*!
   * \brief A virtual member.
   * \return Value of the turbulent frequency.
   */
  virtual su2double GetOmega_Inf(void);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Young Modulus E
   */
  virtual su2double GetTotal_Sens_E(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity for the Poisson's ratio Nu
   */
  virtual su2double GetTotal_Sens_Nu(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the structural density sensitivity
   */
  virtual su2double GetTotal_Sens_Rho(unsigned short iVal);

  /*!
   * \brief A virtual member.
   * \return Value of the structural weight sensitivity
   */
  virtual su2double GetTotal_Sens_Rho_DL(unsigned short iVal);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField
   */
  virtual su2double GetTotal_Sens_EField(unsigned short iEField);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the FEA DV in the region iDVFEA
   */
  virtual su2double GetTotal_Sens_DVFEA(unsigned short iDVFEA);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Young Modulus E
   */
  virtual su2double GetGlobal_Sens_E(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Poisson's ratio Nu
   */
  virtual su2double GetGlobal_Sens_Nu(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the structural density sensitivity
   */
  virtual su2double GetGlobal_Sens_Rho(unsigned short iVal);

  /*!
   * \brief A virtual member.
   * \return Value of the structural weight sensitivity
   */
  virtual su2double GetGlobal_Sens_Rho_DL(unsigned short iVal);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField
   */
  virtual su2double GetGlobal_Sens_EField(unsigned short iEField);
  
  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the FEA DV in the region iDVFEA
   */
  virtual su2double GetGlobal_Sens_DVFEA(unsigned short iDVFEA);

  /*!
   * \brief A virtual member.
   * \return Value of the Young modulus from the adjoint solver
   */
  virtual su2double GetVal_Young(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Poisson's ratio from the adjoint solver
   */
  virtual su2double GetVal_Poisson(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the density for inertial effects, from the adjoint solver
   */
  virtual su2double GetVal_Rho(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the density for dead loads, from the adjoint solver
   */
  virtual su2double GetVal_Rho_DL(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Number of electric field variables from the adjoint solver
   */
  virtual unsigned short GetnEField(void);

  /*!
   * \brief A virtual member.
   * \return Number of design variables from the adjoint solver
   */
  virtual unsigned short GetnDVFEA(void);
  
  /*!
   * \brief A virtual member.
   */
  virtual void ReadDV(CConfig *config);

  /*!
   * \brief A virtual member.
   * \return Pointer to the values of the Electric Field
   */
  virtual su2double GetVal_EField(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Pointer to the values of the design variables
   */
  virtual su2double GetVal_DVFEA(unsigned short iVal);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the sensitivity coefficient.
   */
  virtual su2double GetCSensitivity(unsigned short val_marker, unsigned long val_vertex);
 
  /*!
   * \brief A virtual member.
   * \return A pointer to an array containing a set of constants
   */
  virtual su2double* GetConstants();
  
  /*!
   * \brief A virtual member.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_forcecoeff_history - Value of the force coefficient.
   */
  virtual void SetForceCoeff(su2double val_forcecoeff_history);

  /*!
   * \brief A virtual member.
   * \param[in] val_relaxcoeff_history - Value of the force coefficient.
   */
  virtual void SetRelaxCoeff(su2double val_relaxcoeff_history);

  /*!
   * \brief A virtual member.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  virtual void SetFSI_Residual(su2double val_FSI_residual);

  /*!
   * \brief A virtual member.
   * \param[out] val_forcecoeff_history - Value of the force coefficient.
   */
  virtual su2double GetForceCoeff();

  /*!
   * \brief A virtual member.
   * \param[out] val_relaxcoeff_history - Value of the relax coefficient.
   */
  virtual su2double GetRelaxCoeff();

  /*!
   * \brief A virtual member.
   * \param[out] val_FSI_residual - Value of the residual.
   */
  virtual su2double GetFSI_Residual();
  
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
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  virtual void SetInitialCondition(CGeometry **geometry,
                                   CSolver ***solver_container,
                                   CConfig *config, unsigned long ExtIter);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  virtual void ResetInitialCondition(CGeometry **geometry,
                                     CSolver ***solver_container,
                                     CConfig *config, unsigned long ExtIter);
  
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
                                         unsigned long iOuterIter);
  
  
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
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  virtual void LoadRestart(CGeometry **geometry, CSolver ***solver,
                           CConfig *config, int val_iter, bool val_update_geo);

  /*!
   * \brief Read a native SU2 restart file in ASCII format.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_SU2_Restart_ASCII(CGeometry *geometry, CConfig *config, string val_filename);

  /*!
   * \brief Read a native SU2 restart file in binary format.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_SU2_Restart_Binary(CGeometry *geometry, CConfig *config, string val_filename);

  /*!
   * \brief Read the metadata from a native SU2 restart file (ASCII or binary).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] adjoint - Boolean to identify the restart file of an adjoint run.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_SU2_Restart_Metadata(CGeometry *geometry, CConfig *config, bool adjoint_run, string val_filename);

  /*!
   * \brief Read a native SU2 inlet file in ASCII format.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_InletFile_ASCII(CGeometry *geometry, CConfig *config, string val_filename);

  /*!
   * \brief Load a inlet profile data from file into a particular solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_kind_solver - Solver container position.
   * \param[in] val_kind_marker - Kind of marker to apply the profiles.
   */
  void LoadInletProfile(CGeometry **geometry,
                        CSolver ***solver,
                        CConfig *config,
                        int val_iter,
                        unsigned short val_kind_solver,
                        unsigned short val_kind_marker);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_OFVolFrac(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Stiffness_Penalty(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   */
  virtual void LoadRestart_FSI(CGeometry *geometry, CConfig *config, int val_iter);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void RefGeom_Sensitivity(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void DE_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Stiffness_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] iElem - element parameter.
   * \param[out] iElem_iDe - ID of the Dielectric Elastomer region.
   */
  virtual unsigned short Get_iElem_iDe(unsigned long iElem);
  
  /*!
   * \brief A virtual member.
   * \param[in] i_DV - number of design variable.
   * \param[in] val_EField - value of the design variable.
   */
  virtual void Set_DV_Val(su2double val_EField, unsigned short i_DV);
  
  /*!
   * \brief A virtual member.
   * \param[in] i_DV - number of design variable.
   * \param[out] DV_Val - value of the design variable.
   */
  virtual su2double Get_DV_Val(unsigned short i_DV);
  
  /*!
   * \brief A virtual member.
   * \param[out] val_I - value of the objective function.
   */
  virtual su2double Get_val_I(void);
  
  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point i of the Mass Matrix.
   * \param[in] jPoint - Point j of the Mass Matrix.
   * \param[in] iVar - Variable i of the Mass Matrix submatrix.
   * \param[in] iVar - Variable j of the Mass Matrix submatrix.
   */
  virtual su2double Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar);
  
  /*!
   * \brief Gauss method for solving a linear system.
   * \param[in] A - Matrix Ax = b.
   * \param[in] rhs - Right hand side.
   * \param[in] nVar - Number of variables.
   */
  void Gauss_Elimination(su2double** A, su2double* rhs, unsigned short nVar);
  
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
   * \param[in] config - The particular config.
   */
  virtual void SetAdjoint_OutputMesh(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  virtual void ExtractAdjoint_Solution(CGeometry *geometry,  CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  virtual void ExtractAdjoint_Geometry(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  virtual void ExtractAdjoint_CrossTerm(CGeometry *geometry,  CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  virtual void ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry,  CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  virtual void ExtractAdjoint_CrossTerm_Geometry_Flow(CGeometry *geometry,  CConfig *config);
  
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
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Set_Prestretch(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Set_ReferenceGeometry(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Set_ElementProperties(CGeometry *geometry, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] CurrentTime - Current time step.
   * \param[in] RampTime - Time for application of the ramp.*
   * \param[in] config - Definition of the particular problem.
   */
  virtual su2double Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_MassRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  
  virtual void Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  
  virtual void Compute_DeadLoad(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  
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
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Compute_IntegrationConstants(CConfig *config);
  
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
   * \param[in] Value of the load increment for nonlinear structural analysis
   */
  virtual void SetLoad_Increment(su2double val_loadIncrement);
  
  /*!
   * \brief A virtual member.
   * \param[in] Value of the load increment for nonlinear structural analysis
   */
  virtual su2double GetLoad_Increment(void);

  /*!
   * \brief A virtual member.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] Value of freestream density.
   */
  virtual void SetDensity_Inf(su2double rho_inf);

  /*!
   * \brief A virtual member.
   * \param[in] val_dim - Index of the velocity vector.
   * \param[in] val_velocity - Value of the velocity.
   */
  virtual void SetVelocity_Inf(unsigned short val_dim, su2double val_velocity);

  /*!
   * \brief A virtual member.
   * \param[in] kind_recording - Kind of AD recording.
   */
  virtual void SetRecording(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] kind_recording - Kind of AD recording.
   */
  virtual void SetMesh_Recording(CGeometry **geometry, CVolumetricMovement *grid_movement,
      CConfig *config);
  
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
  
  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFreeStream_Solution(CConfig *config);

  /*!
   * \brief A virtual member.
   */
  virtual su2double* GetVecSolDOFs(void);

  /*!
   * \brief A virtual member.
   */
  virtual unsigned long GetnDOFsGlobal(void);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetTauWall_WF(CGeometry *geometry, CSolver** solver_container, CConfig* config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetNuTilde_WF(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                             CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void InitTurboContainers(CGeometry *geometry, CConfig *config);

  /*!
   * \brief virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the average is evaluated.
   */
  virtual void PreprocessAverage(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag);


  /*!
   * \brief virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the average is evaluated.
   */
  virtual void TurboAverageProcess(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag);

  /*!
   * \brief virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
   virtual void GatherInOutAverageValues(CConfig *config, CGeometry *geometry);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Density on the surface <i>val_marker</i>.
   */
  virtual su2double GetAverageDensity(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Pressure on the surface <i>val_marker</i>.
   */
  virtual su2double GetAveragePressure(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  virtual su2double* GetAverageTurboVelocity(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Nu on the surface <i>val_marker</i>.
   */
  virtual su2double GetAverageNu(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Kine on the surface <i>val_marker</i>.
   */
  virtual su2double GetAverageKine(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Omega on the surface <i>val_marker</i>.
   */
  virtual su2double GetAverageOmega(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Nu on the surface <i>val_marker</i>.
   */
  virtual su2double GetExtAverageNu(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Kine on the surface <i>val_marker</i>.
   */
  virtual su2double GetExtAverageKine(unsigned short valMarker, unsigned short iSpan);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Omega on the surface <i>val_marker</i>.
   */
  virtual su2double GetExtAverageOmega(unsigned short valMarker, unsigned short iSpan);


  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Density on the surface <i>val_marker</i>.
   */
  virtual void SetExtAverageDensity(unsigned short valMarker, unsigned short valSpan, su2double valDensity);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Pressure on the surface <i>val_marker</i>.
   */
  virtual void SetExtAveragePressure(unsigned short valMarker, unsigned short valSpan, su2double valPressure);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  virtual void SetExtAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan, unsigned short valIndex, su2double valTurboVelocity);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Nu on the surface <i>val_marker</i>.
   */
  virtual void SetExtAverageNu(unsigned short valMarker, unsigned short valSpan, su2double valNu);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Kine on the surface <i>val_marker</i>.
   */
  virtual void SetExtAverageKine(unsigned short valMarker, unsigned short valSpan, su2double valKine);

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Omega on the surface <i>val_marker</i>.
   */
  virtual void SetExtAverageOmega(unsigned short valMarker, unsigned short valSpan, su2double valOmega);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of inlet pressure.
   */
  virtual su2double GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet normal velocity.
   */
  virtual su2double* GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet density.
   */
  virtual su2double GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet pressure.
   */
  virtual su2double GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet normal velocity.
   */
  virtual su2double* GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetKineIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetNuIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetKineOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  virtual su2double GetNuOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetDensityIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetPressureIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetTurboVelocityIn(su2double* value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetDensityOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetPressureOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetTurboVelocityOut(su2double* value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetKineIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetOmegaIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetNuIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetKineOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetOmegaOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  virtual void SetNuOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFreeStream_TurboSolution(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   */
  virtual void SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                               CConfig *config, unsigned short iMesh);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetRoe_Dissipation(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] solver - Solver container
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetDES_LengthScale(CSolver** solver, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Routine that sets the flag controlling implicit treatment for periodic BCs.
   * \param[in] val_implicit_periodic - Flag controlling implicit treatment for periodic BCs.
   */
  void SetImplicitPeriodic(bool val_implicit_periodic);
  
  /*!
   * \brief Routine that sets the flag controlling solution rotation for periodic BCs.
   * \param[in] val_implicit_periodic - Flag controlling solution rotation for periodic BCs.
   */
  void SetRotatePeriodic(bool val_rotate_periodic);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  virtual void ComputeVerificationError(CGeometry *geometry, CConfig *config);
  
protected:
  /*!
   * \brief Allocate the memory for the verification solution, if necessary.
   * \param[in] nDim   - Number of dimensions of the problem.
   * \param[in] nVar   - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetVerificationSolution(unsigned short nDim, 
                               unsigned short nVar, 
                               CConfig        *config);
};

/*!
 * \class CBaselineSolver
 * \brief Main class for defining a baseline solution from a restart file (for output).
 * \author F. Palacios, T. Economon.
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
  CBaselineSolver(CGeometry *geometry, CConfig *config);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] nVar - Number of variables.
   * \param[in] field_names - Vector of variable names.
   */
  CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short nVar, vector<string> field_names);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CBaselineSolver(void);
  
  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);
  
  /*!
   * \brief Load a FSI solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   */
  void LoadRestart_FSI(CGeometry *geometry, CConfig *config, int val_iter);

  /*!
   * \brief Set the number of variables and string names from the restart file.
   * \param[in] config - Definition of the particular problem.
   */
  void SetOutputVariables(CGeometry *geometry, CConfig *config);

};

/*!
 * \class CBaselineSolver_FEM
 * \brief Main class for defining a baseline solution from a restart file for the DG-FEM solver output.
 * \author T. Economon.
 * \version 6.2.0 "Falcon"
 */
class CBaselineSolver_FEM : public CSolver {
protected:

  unsigned long nDOFsLocTot;    /*!< \brief Total number of local DOFs, including halos. */
  unsigned long nDOFsLocOwned;  /*!< \brief Number of owned local DOFs. */
  unsigned long nDOFsGlobal;    /*!< \brief Number of global DOFs. */

  unsigned long nVolElemTot;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned;  /*!< \brief Number of owned local volume elements. */
  CVolumeElementFEM *volElem;   /*!< \brief Array of the local volume elements, including halos. */

  vector<su2double> VecSolDOFs;    /*!< \brief Vector, which stores the solution variables in all the DOFs. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CBaselineSolver_FEM(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CBaselineSolver_FEM(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CBaselineSolver_FEM(void);

  /*!
   * \brief Set the number of variables and string names from the restart file.
   * \param[in] config - Definition of the particular problem.
   */
  void SetOutputVariables(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

  /*!
   * \brief Get a pointer to the vector of the solution degrees of freedom.
   * \return Pointer to the vector of the solution degrees of freedom.
   */
  su2double* GetVecSolDOFs(void);

};

/*!
 * \class CEulerSolver
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 */
class CEulerSolver : public CSolver {
protected:
  
  su2double
  Mach_Inf,  /*!< \brief Mach number at the infinity. */
  Density_Inf,  /*!< \brief Density at the infinity. */
  Energy_Inf,      /*!< \brief Energy at the infinity. */
  Temperature_Inf,      /*!< \brief Energy at the infinity. */
  Pressure_Inf,    /*!< \brief Pressure at the infinity. */
  *Velocity_Inf;    /*!< \brief Flow Velocity vector at the infinity. */
  
  su2double
  *CD_Inv,  /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Inv,      /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,    /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Inv,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Inv,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Inv,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,      /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,      /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,      /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Inv, /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Inv, /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Inv,        /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Inv,        /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Inv,      /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Inv,      /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  *CEquivArea_Inv,        /*!< \brief Equivalent area (inviscid contribution) for each boundary. */
  *CNearFieldOF_Inv,        /*!< \brief Near field pressure (inviscid contribution) for each boundary. */
  *CD_Mnt,  /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Mnt,      /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Mnt,    /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Mnt,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Mnt,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Mnt,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Mnt,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Mnt,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Mnt,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Mnt,      /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Mnt,      /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Mnt,      /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Mnt, /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Mnt, /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Mnt, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Mnt, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Mnt,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Mnt,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Mnt,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Mnt,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Mnt,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Mnt,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Mnt,        /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Mnt,        /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Mnt,      /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Mnt,      /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  *CEquivArea_Mnt,        /*!< \brief Equivalent area (inviscid contribution) for each boundary. */
  **CPressure,    /*!< \brief Pressure coefficient for each boundary and vertex. */
  **CPressureTarget,    /*!< \brief Target Pressure coefficient for each boundary and vertex. */
  **HeatFlux,    /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **HeatFluxTarget,    /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **YPlus,    /*!< \brief Yplus for each boundary and vertex. */
  ***CharacPrimVar,    /*!< \brief Value of the characteristic variables at each boundary. */
  ***DonorPrimVar,    /*!< \brief Value of the donor variables at each boundary. */
  *ForceInviscid,    /*!< \brief Inviscid force for each boundary. */
  *MomentInviscid,  /*!< \brief Inviscid moment for each boundary. */
  *ForceMomentum,    /*!< \brief Inviscid force for each boundary. */
  *MomentMomentum;  /*!< \brief Inviscid moment for each boundary. */
  su2double *Inflow_MassFlow,  /*!< \brief Mass flow rate for each boundary. */
  *Exhaust_MassFlow,  /*!< \brief Mass flow rate for each boundary. */
  *Inflow_Pressure,  /*!< \brief Fan face pressure for each boundary. */
  *Inflow_Mach,  /*!< \brief Fan face mach number for each boundary. */
  *Inflow_Area,  /*!< \brief Boundary total area. */
  *Exhaust_Area,  /*!< \brief Boundary total area. */
  *Exhaust_Pressure,  /*!< \brief Fan face pressure for each boundary. */
  *Exhaust_Temperature,  /*!< \brief Fan face mach number for each boundary. */
  Inflow_MassFlow_Total,  /*!< \brief Mass flow rate for each boundary. */
  Exhaust_MassFlow_Total,  /*!< \brief Mass flow rate for each boundary. */
  Inflow_Pressure_Total,  /*!< \brief Fan face pressure for each boundary. */
  Inflow_Mach_Total,  /*!< \brief Fan face mach number for each boundary. */
  InverseDesign;  /*!< \brief Inverse design functional for each boundary. */
  unsigned long **DonorGlobalIndex;    /*!< \brief Value of the donor global index. */
  su2double **ActDisk_DeltaP,    /*!< \brief Value of the Delta P. */
  **ActDisk_DeltaT;    /*!< \brief Value of the Delta T. */
  su2double **Inlet_Ptotal,    /*!< \brief Value of the Total P. */
  **Inlet_Ttotal,    /*!< \brief Value of the Total T. */
  ***Inlet_FlowDir;    /*!< \brief Value of the Flow Direction. */
  
  su2double
  AllBound_CD_Inv,  /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Inv,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Inv,      /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Inv,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Inv,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Inv,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Inv,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Inv,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Inv,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Inv,      /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Inv,      /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Inv,      /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Inv,      /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Inv,      /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Inv,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Inv,      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEquivArea_Inv,      /*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CNearFieldOF_Inv;      /*!< \brief Near-Field press coefficient (inviscid contribution) for all the boundaries. */
  
  su2double
  AllBound_CD_Mnt,  /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CL_Mnt,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CSF_Mnt,      /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CMx_Mnt,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CMy_Mnt,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CMz_Mnt,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CoPx_Mnt,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CoPy_Mnt,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CoPz_Mnt,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CFx_Mnt,      /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CFy_Mnt,      /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CFz_Mnt,      /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CEff_Mnt,      /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CMerit_Mnt,      /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
   AllBound_CT_Mnt,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
   AllBound_CQ_Mnt;      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */
  
  su2double
  Total_ComboObj, /*!< \brief Total 'combo' objective for all monitored boundaries */
  AoA_Prev, /*!< \brief Old value of the AoA for fixed lift mode. */
  Total_CD, /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CL,    /*!< \brief Total lift coefficient for all the boundaries. */
  Total_CL_Prev,    /*!< \brief Total lift coefficient for all the boundaries (fixed lift mode). */
  Total_SolidCD, /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CD_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_NetThrust, /*!< \brief Total drag coefficient for all the boundaries. */
  Total_Power, /*!< \brief Total drag coefficient for all the boundaries. */
  Total_ReverseFlow, /*!< \brief Total drag coefficient for all the boundaries. */
  Total_IDC,        /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDC_Mach,        /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDR,        /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_DC60,        /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_MFR,     /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Prop_Eff,     /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_ByPassProp_Eff,     /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Adiab_Eff,     /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Poly_Eff,     /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Custom_ObjFunc,        /*!< \brief Total custom objective function for all the boundaries. */
  Total_CSF,    /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CMx,      /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMx_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMy,      /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMy_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMz,      /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CMz_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CoPx,      /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CoPy,      /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CoPz,      /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CFx,      /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy,      /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz,      /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CEff,      /*!< \brief Total efficiency coefficient for all the boundaries. */
  Total_CMerit,      /*!< \brief Total rotor Figure of Merit for all the boundaries. */
  Total_CT,    /*!< \brief Total thrust coefficient for all the boundaries. */
  Total_CQ,    /*!< \brief Total torque coefficient for all the boundaries. */
  Total_Heat,    /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat, /*!< \brief Maximum heat flux on all boundaries. */
  Total_AeroCD,      /*!< \brief Total aero drag coefficient for all the boundaries. */
  Total_CEquivArea,      /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_CNearFieldOF,      /*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
  Total_CpDiff,      /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_HeatFluxDiff,      /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_MassFlowRate;     /*!< \brief Total Mass Flow Rate on monitored boundaries. */
  su2double *Surface_CL,   /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CD,          /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSF,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CEff,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,            /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,            /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,            /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,            /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,            /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz,            /*!< \brief z Moment coefficient for each monitoring surface. */
  *Surface_HF_Visc,            /*!< \brief Total (integrated) heat flux for each monitored surface. */
  *Surface_MaxHF_Visc;         /*!< \brief Maximum heat flux for each monitored surface. */
  
  su2double *SecondaryVar_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
  *SecondaryVar_j;      /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double *PrimVar_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
  *PrimVar_j;      /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double **LowMach_Precontioner; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
  unsigned long nMarker,        /*!< \brief Total number of markers using the grid information. */
  *nVertex;       /*!< \brief Store nVertex at each marker for deallocation */
  bool space_centered,  /*!< \brief True if space centered scheeme used. */
  euler_implicit,      /*!< \brief True if euler implicit scheme used. */
  least_squares;        /*!< \brief True if computing gradients by least squares. */
  su2double Gamma;                  /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;        /*!< \brief Fluids's Gamma - 1.0  . */
  
  su2double *Primitive,    /*!< \brief Auxiliary nPrimVar vector. */
  *Primitive_i,        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Primitive_j;        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */
  
  su2double *Secondary,    /*!< \brief Auxiliary nPrimVar vector. */
  *Secondary_i,        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Secondary_j;        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */
  
  su2double Cauchy_Value,  /*!< \brief Summed value of the convergence indicator. */
  Cauchy_Func;      /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter;  /*!< \brief Number of elements of the Cauchy serial. */
  su2double *Cauchy_Serie;      /*!< \brief Complete Cauchy serial. */
  su2double Old_Func,  /*!< \brief Old value of the objective function (the function which is monitored). */
  New_Func;      /*!< \brief Current value of the objective function (the function which is monitored). */
  su2double AoA_old;  /*!< \brief Old value of the angle of attack (monitored). */
  unsigned long AoA_Counter;
  bool AoA_FD_Change;
  unsigned long BCThrust_Counter;
  unsigned short nSpanWiseSections;  /*!< \brief Number of span-wise sections. */
  unsigned short nSpanMax; /*!< \brief Max number of maximum span-wise sections for all zones */
  unsigned short nMarkerTurboPerf;  /*!< \brief Number of turbo performance. */

  CFluidModel  *FluidModel;  /*!< \brief fluid model used in the solver */

  /*--- Turbomachinery Solver Variables ---*/
  su2double *** AverageFlux,
            ***SpanTotalFlux,
            ***AverageVelocity,
            ***AverageTurboVelocity,
            ***OldAverageTurboVelocity,
            ***ExtAverageTurboVelocity,
             **AveragePressure,
             **OldAveragePressure,
             **RadialEquilibriumPressure,
             **ExtAveragePressure,
             **AverageDensity,
             **OldAverageDensity,
             **ExtAverageDensity,
             **AverageNu,
             **AverageKine,
             **AverageOmega,
             **ExtAverageNu,
             **ExtAverageKine,
             **ExtAverageOmega;

  su2double  **DensityIn,
             **PressureIn,
             ***TurboVelocityIn,
             **DensityOut,
             **PressureOut,
             ***TurboVelocityOut,
             **KineIn,
             **OmegaIn,
             **NuIn,
             **KineOut,
             **OmegaOut,
             **NuOut;
  
  complex<su2double> ***CkInflow,
                     ***CkOutflow1,
                     ***CkOutflow2;

 /*--- End of Turbomachinery Solver Variables ---*/

  /* Sliding meshes variables */

  su2double ****SlidingState;
  int **SlidingStateNodes;

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
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh);

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
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);
  
  /*!
   * \brief Compute a pressure sensor switch.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute Ducros Sensor for Roe Dissipation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config);
  
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
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config);
  
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
   * \brief Compute the max eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Nearfield(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Interface(CGeometry *geometry, CConfig *config);

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  void Evaluate_ObjFunc(CConfig *config);
  
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
  * \brief Impose the interface state across sliding meshes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] solver_container - Container vector with all the solutions.
  * \param[in] conv_numerics - Description of the numerical method.
  * \param[in] visc_numerics - Description of the numerical method.
  * \param[in] config - Definition of the particular problem.
  */
  void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config);
    
  /*!
   * \brief Impose the engine inflow boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
  void BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
  void BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                         CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose the engine inflow boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
  void BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker, bool val_inlet_surface);
  
  /*!
   * \brief Impose the interface boundary condition using the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                             CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose the near-field boundary condition using the residual.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                             CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                   CNumerics *numerics, CConfig *config);
  
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
   * \brief Impose the boundary condition using characteristic recostruction.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container,
      CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief It computes Fourier transformation for the needed quantities along the pitch for each span in turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  void PreprocessBC_Giles(CGeometry *geometry, CConfig *config, CNumerics *conv_numerics,  unsigned short marker_flag);

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
  void BC_Giles(CGeometry *geometry, CSolver **solver_container,
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
   * \brief Impose a custom or verification boundary condition.
   * \param[in] geometry         - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics    - Description of the convective numerical method.
   * \param[in] visc_numerics    - Description of the viscous numerical method.
   * \param[in] config           - Definition of the particular problem.
   * \param[in] val_marker       - Surface marker where the boundary condition is applied.
   */
  void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
   * \brief Set the new solution variables to the current solution value for classical RK.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void Set_NewSolution(CGeometry *geometry);

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
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep);

  /*!
   * \brief Compute the Fan face Mach number.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  void GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);
  
  /*!
   * \brief Update the AoA and freestream velocity at the farfield.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                           CConfig *config, unsigned short iMesh, bool Output);
  
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
  void Pressure_Forces(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Momentum_Forces(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute turbomachinery performance.
   * \param[in] solver - solver containing the outlet information.
   * \param[in] inMarker - marker related to the inlet.
   * \param[in] outMarker - marker related to the outlet.
   */
  void TurboPerformance(CSolver *solver,  CConfig *config, unsigned short inMarker,  unsigned short outMarker, unsigned short Kind_TurboPerf , unsigned short inMarkerTP );
  
  /*!
   * \brief Compute turbomachinery performance.
   * \param[in] solver - solver containing the outlet information.
   * \param[in] inMarker - marker related to the inlet.
   * \param[in] outMarker - marker related to the outlet.
   */
  void StoreTurboPerformance(CSolver *solver,  unsigned short inMarkerTP );
  
 /*!
  * \brief Get the outer state for fluid interface nodes.
  * \param[in] val_marker - marker index
  * \param[in] val_vertex - vertex index
  * \param[in] val_state  - requested state component
  */
  su2double GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index);

  /*!
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCL_Inv(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CL(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF(unsigned short val_marker);
  
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
  su2double GetSurface_CL_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Inv(unsigned short val_marker);
  
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
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CL_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CEff_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFx_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFy_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFz_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMx_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMy_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMz_Mnt(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCD_Inv(unsigned short val_marker);
  
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
  su2double GetCSF_Inv(unsigned short val_marker);
  
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
  su2double GetTotal_CSF(void);
  
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
   * \brief Provide the total (inviscid + viscous) non dimensional aero CD.
   * \return Value of the Aero CD coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_AeroCD(void);

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
   * \author H. Kline
   * \brief Add to the value of the total 'combo' objective.
   * \param[in] val_obj - Value of the contribution to the 'combo' objective.
   */
  void AddTotal_ComboObj(su2double val_obj);
  
  /*!
   * \brief Set the value of the Equivalent Area coefficient.
   * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
   */
  void SetTotal_CEquivArea(su2double val_cequivarea);
  
  /*!
   * \brief Set the value of the Aero drag.
   * \param[in] val_cequivarea - Value of the aero drag.
   */
  void SetTotal_AeroCD(su2double val_aerocd);
  
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
   * \author H. Kline
   * \brief Set the total "combo" objective (weighted sum of other values).
   * \param[in] ComboObj - Value of the combined objective.
   */
  void SetTotal_ComboObj(su2double ComboObj);
  
  /*!
   * \author H. Kline
   * \brief Provide the total "combo" objective (weighted sum of other values).
   * \return Value of the "combo" objective values.
   */
  su2double GetTotal_ComboObj(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CL(void);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CD(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_NetThrust(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_Power(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_SolidCD(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_ReverseFlow(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_MFR(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_Prop_Eff(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_ByPassProp_Eff(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_Adiab_Eff(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_Poly_Eff(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_IDC(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_IDC_Mach(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_IDR(void);
 
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_DC60(void);
  
  /*!
   * \brief Provide the total custom objective function.
   * \return Value of the custom objective function.
   */
  su2double GetTotal_Custom_ObjFunc(void);

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
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CoPx(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CoPy(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CoPz(void);

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
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_CD(su2double val_Total_CD);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  void SetTotal_CL(su2double val_Total_CL);

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_NetThrust(su2double val_Total_NetThrust);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_Power(su2double val_Total_Power);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_SolidCD(su2double val_Total_SolidCD);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_ReverseFlow(su2double val_ReverseFlow);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_MFR(su2double val_Total_MFR);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_Prop_Eff(su2double val_Total_Prop_Eff);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_Poly_Eff(su2double val_Total_Poly_Eff);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_IDC(su2double val_Total_IDC);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_IDC_Mach(su2double val_Total_IDC_Mach);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_IDR(su2double val_Total_IDR);
  
  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_DC60(su2double val_Total_DC60);

  /*!
   * \brief Set the value of the custom objective function.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  void SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight);
  
  /*!
   * \brief Add the value of the custom objective function.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  void AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight);

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Inv(void);
  
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
  su2double GetAllBound_CoPx_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPy_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPz_Inv(void);
  
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
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CEff_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMx_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMy_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMz_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPx_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPy_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPz_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFx_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFy_Mnt(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFz_Mnt(void);
  
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
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value);
  
  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  su2double *GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value);
  
  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  su2double GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  unsigned long GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  su2double GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex, su2double val_deltap);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  su2double GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex, su2double val_deltat);
  
  /*!
   * \brief Value of the total temperature at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total temperature is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is evaluated.
   * \return Value of the total temperature
   */
  su2double GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the total pressure at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total pressure is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is evaluated.
   * \return Value of the total pressure
   */
  su2double GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A component of the unit vector representing the flow direction at an inlet boundary.
   * \param[in] val_marker - Surface marker where the flow direction is evaluated
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is evaluated
   * \param[in] val_dim - The component of the flow direction unit vector to be evaluated
   * \return Component of a unit vector representing the flow direction.
   */
  su2double GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim);
  
  /*!
   * \brief Set the value of the total temperature at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total temperature is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is set.
   * \param[in] val_ttotal - Value of the total temperature
   */
  void SetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ttotal);
  
  /*!
   * \brief Set the value of the total pressure at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total pressure is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is set.
   * \param[in] val_ptotal - Value of the total pressure
   */
  void SetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ptotal);
  
  /*!
   * \brief Set a component of the unit vector representing the flow direction at an inlet boundary.
   * \param[in] val_marker - Surface marker where the flow direction is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is set.
   * \param[in] val_dim - The component of the flow direction unit vector to be set
   * \param[in] val_flowdir - Component of a unit vector representing the flow direction.
   */
  void SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir);

  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(CConfig* config, unsigned short iMarker);

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the set of value imposed at an inlet.
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(su2double *val_inlet,
                             unsigned long val_inlet_point,
                             unsigned short val_kind_marker,
                             string val_marker,
                             CGeometry *geometry,
                             CConfig *config);

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config);

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
   * \brief Set the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex);
      
  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  void SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component);

  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  void SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value);

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex);
    
  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
  
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
  
  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(CConfig *config);

  /*!
   * \brief Initilize turbo containers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void InitTurboContainers(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_TurboSolution(CConfig *config);

  /*!
   * \brief It computes average quantities along the span for turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  void PreprocessAverage(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag);

  /*!
   * \brief It computes average quantities along the span for turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  void TurboAverageProcess(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag);

  /*!
   * \brief it performs a mixed out average of the nodes of a boundary.
   * \param[in] val_init_pressure -  initial pressure value
   * \param[in] val_Averaged_Flux - flux averaged values.
   * \param[in] val_normal - normal vector.
   * \param[in] pressure_mix - value of the mixed-out avaraged pressure.
   * \param[in] density_miz - value of the mixed-out avaraged density.
   */
  void MixedOut_Average (CConfig *config, su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double& pressure_mix, su2double& density_mix);

  /*!
   * \brief It gathers into the master node average quantities at inflow and outflow needed for turbomachinery analysis.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void GatherInOutAverageValues(CConfig *config, CGeometry *geometry);

  /*!
   * \brief it take a velocity in the cartesian reference of framework and transform into the turbomachinery frame of reference.
   * \param[in] cartesianVelocity - cartesian components of velocity vector.
   * \param[in] turboNormal - normal vector in the turbomachinery frame of reference.
   * \param[in] turboVelocity - velocity vector in the turbomachinery frame of reference.
   */
  void ComputeTurboVelocity(su2double *cartesianVelocity, su2double *turboNormal, su2double *turboVelocity, unsigned short marker_flag, unsigned short marker_kindturb);

  /*!
   * \brief it take a velocity in the cartesian reference of framework and transform into the turbomachinery frame of reference.
   * \param[in] cartesianVelocity - cartesian components of velocity vector.
   * \param[in] turboNormal - normal vector in the turbomachinery frame of reference.
   * \param[in] turboVelocity - velocity vector in the turbomachinery frame of reference.
   */
  void ComputeBackVelocity(su2double *turboVelocity, su2double *turboNormal, su2double *cartesianVelocity, unsigned short marker_flag, unsigned short marker_kindturb);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Density on the surface <i>val_marker</i>.
   */
  su2double GetAverageDensity(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average pressure at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Pressure on the surface <i>val_marker</i>.
   */
  su2double GetAveragePressure(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average turbo velocity average at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  su2double* GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Nu on the surface <i>val_marker</i>.
   */
  su2double GetAverageNu(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Kine on the surface <i>val_marker</i>.
   */
  su2double GetAverageKine(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Omega on the surface <i>val_marker</i>.
   */
  su2double GetAverageOmega(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Nu on the surface <i>val_marker</i>.
   */
  su2double GetExtAverageNu(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Kine on the surface <i>val_marker</i>.
   */
  su2double GetExtAverageKine(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Omega on the surface <i>val_marker</i>.
   */
  su2double GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan);

  /*!
   * \brief Set the external average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valDensity - value to set.
   */
  void SetExtAverageDensity(unsigned short valMarker, unsigned short valSpan, su2double valDensity);

  /*!
   * \brief Set the external average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valPressure - value to set.
   */
  void SetExtAveragePressure(unsigned short valMarker, unsigned short valSpan, su2double valPressure);

  /*!
   * \brief Set the external the average turbo velocity average at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  void SetExtAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan, unsigned short valIndex, su2double valTurboVelocity);

  /*!
   * \brief Set the external average turbulent Nu at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valNu - value to set.
   */
  void SetExtAverageNu(unsigned short valMarker, unsigned short valSpan, su2double valNu);

  /*!
   * \brief Set the external average turbulent Kine at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valKine - value to set.
   */
  void SetExtAverageKine(unsigned short valMarker, unsigned short valSpan, su2double valKine);

  /*!
   * \brief Set the external average turbulent Omega at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valOmega - value to set.
   */
  void SetExtAverageOmega(unsigned short valMarker, unsigned short valSpan, su2double valOmega);

  /*!
   * \brief Provide the inlet density to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the inlet pressure to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of inlet pressure.
   */
  su2double GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the inlet normal velocity to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet normal velocity.
   */
  su2double* GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the outlet density to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet density.
   */
  su2double GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the outlet pressure to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet pressure.
   */
  su2double GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the outlet normal velocity to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet normal velocity.
   */
  su2double* GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the inlet turbulent kei to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetKineIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the inlet turbulent omega to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the inlet turbulent nu to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetNuIn(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the outlet turbulent kei to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetKineOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the outlet turbulent omega to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Provide the outlet turbulent nu to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  su2double GetNuOut(unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set inlet density.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetDensityIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set inlet pressure.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetPressureIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set inlet normal velocity.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetTurboVelocityIn(su2double* value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set outlet density.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetDensityOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set outlet pressure.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetPressureOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set outlet normal velocity.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetTurboVelocityOut(su2double* value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set inlet turbulent kei.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetKineIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);
  /*!
   * \brief Set inlet turbulent omega.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetOmegaIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);
  /*!
   * \brief Set inlet turbulent Nu.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetNuIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan);

  /*!
   * \brief Set outlet turbulent kei.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetKineOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);
  /*!
   * \brief Set Outlet turbulent omega.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetOmegaOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);
  /*!
   * \brief Set outlet turbulent Nu.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  void SetNuOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan);
  
  /*!
   * \brief Compute the global error measures (L2, Linf) for verification cases.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVerificationError(CGeometry *geometry, CConfig *config);
};

/*!
 * \class CIncEulerSolver
 * \brief Main class for defining the incompressible Euler flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncEulerSolver : public CSolver {
protected:
  
  su2double
  Density_Inf,  /*!< \brief Density at the infinity. */
  Pressure_Inf,    /*!< \brief Pressure at the infinity. */
  *Velocity_Inf,    /*!< \brief Flow Velocity vector at the infinity. */
  Temperature_Inf;      /*!< \brief Temperature at infinity. */

  su2double
  *CD_Inv,  /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Inv,      /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,    /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Inv,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Inv,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Inv,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,      /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,      /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,      /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Inv, /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Inv, /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Inv,        /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Inv,        /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Inv,      /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Inv,      /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  *CD_Mnt,  /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Mnt,      /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Mnt,    /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Mnt,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Mnt,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Mnt,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Mnt,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Mnt,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Mnt,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Mnt,      /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Mnt,      /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Mnt,      /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Mnt, /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Mnt, /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Mnt, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Mnt, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Mnt,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Mnt,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Mnt,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Mnt,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Mnt,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Mnt,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Mnt,        /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Mnt,        /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Mnt,      /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Mnt,      /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  **CPressure,    /*!< \brief Pressure coefficient for each boundary and vertex. */
  **CPressureTarget,    /*!< \brief Target Pressure coefficient for each boundary and vertex. */
  **HeatFlux,    /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **HeatFluxTarget,    /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **YPlus,    /*!< \brief Yplus for each boundary and vertex. */
  ***CharacPrimVar,    /*!< \brief Value of the characteristic variables at each boundary. */
  *ForceInviscid,    /*!< \brief Inviscid force for each boundary. */
  *MomentInviscid,  /*!< \brief Inviscid moment for each boundary. */
  *ForceMomentum,    /*!< \brief Inviscid force for each boundary. */
  *MomentMomentum,  /*!< \brief Inviscid moment for each boundary. */
  InverseDesign;  /*!< \brief Inverse design functional for each boundary. */
  su2double **Inlet_Ptotal,    /*!< \brief Value of the Total P. */
  **Inlet_Ttotal,    /*!< \brief Value of the Total T. */
  ***Inlet_FlowDir;    /*!< \brief Value of the Flow Direction. */

  su2double
  AllBound_CD_Inv,  /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Inv,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Inv,      /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Inv,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Inv,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Inv,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Inv,      /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Inv,      /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Inv,      /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Inv,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Inv,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Inv,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Inv,      /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Inv,      /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Inv,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Inv;      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */


  su2double
  AllBound_CD_Mnt,  /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Mnt,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Mnt,      /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Mnt,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Mnt,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Mnt,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Mnt,      /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Mnt,      /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Mnt,      /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Mnt,      /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Mnt,      /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Mnt,      /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Mnt,      /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Mnt,      /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Mnt,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Mnt;      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */

  su2double
  AoA_Prev, /*!< \brief Old value of the AoA for fixed lift mode. */
  Total_ComboObj, /*!< \brief Total 'combo' objective for all monitored boundaries */
  Total_CD, /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CD_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CL,    /*!< \brief Total lift coefficient for all the boundaries. */
  Total_CL_Prev,    /*!< \brief Total lift coefficient for all the boundaries (fixed lift mode). */
  Total_CSF,    /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CMx,      /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMx_Prev,      /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMy,      /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMy_Prev,      /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMz,      /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CMz_Prev,      /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CoPx,      /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CoPy,      /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CoPz,      /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CFx,      /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy,      /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz,      /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CEff,      /*!< \brief Total efficiency coefficient for all the boundaries. */
  Total_CMerit,      /*!< \brief Total rotor Figure of Merit for all the boundaries. */
  Total_CT,    /*!< \brief Total thrust coefficient for all the boundaries. */
  Total_CQ,    /*!< \brief Total torque coefficient for all the boundaries. */
  Total_Heat,    /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat, /*!< \brief Maximum heat flux on all boundaries. */
  Total_CpDiff,      /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_HeatFluxDiff,      /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_Custom_ObjFunc,        /*!< \brief Total custom objective function for all the boundaries. */
  Total_MassFlowRate;     /*!< \brief Total Mass Flow Rate on monitored boundaries. */
  su2double *Surface_CL,   /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CD,          /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSF,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CEff,     /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,            /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,            /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,            /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,            /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,            /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz,            /*!< \brief z Moment coefficient for each monitoring surface. */
  *Surface_HF_Visc,     /*!< \brief Total (integrated) heat flux for each monitored surface. */
  *Surface_MaxHF_Visc;  /*!< \brief Maximum heat flux for each monitored surface. */

  su2double *SecondaryVar_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
  *SecondaryVar_j;      /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double *PrimVar_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
  *PrimVar_j;      /*!< \brief Auxiliary vector for storing the solution at point j. */
  unsigned long nMarker,        /*!< \brief Total number of markers using the grid information. */
  *nVertex;       /*!< \brief Store nVertex at each marker for deallocation */
  bool space_centered,  /*!< \brief True if space centered scheeme used. */
  euler_implicit,      /*!< \brief True if euler implicit scheme used. */
  least_squares;        /*!< \brief True if computing gradients by least squares. */
  su2double Gamma;                  /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;        /*!< \brief Fluids's Gamma - 1.0  . */
  
  su2double *Primitive,    /*!< \brief Auxiliary nPrimVar vector. */
  *Primitive_i,        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Primitive_j;        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  su2double Cauchy_Value,  /*!< \brief Summed value of the convergence indicator. */
  Cauchy_Func;      /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter;  /*!< \brief Number of elements of the Cauchy serial. */
  su2double *Cauchy_Serie;      /*!< \brief Complete Cauchy serial. */
  su2double Old_Func,  /*!< \brief Old value of the objective function (the function which is monitored). */
  New_Func;      /*!< \brief Current value of the objective function (the function which is monitored). */
  su2double AoA_old;  /*!< \brief Old value of the angle of attack (monitored). */
  unsigned long AoA_Counter;
  
  CFluidModel  *FluidModel;  /*!< \brief fluid model used in the solver */
  su2double **Preconditioner; /*!< \brief Auxiliary matrix for storing the low speed preconditioner. */

  /* Sliding meshes variables */

  su2double ****SlidingState;
  int **SlidingStateNodes;

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CIncEulerSolver(void);
  
  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIncEulerSolver(void);

  /*!
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh);

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
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  su2double GetPressure_Inf(void);
  
    /*!
   * \brief Get the temperature value at infinity.
   * \return Value of the temperature at infinity.
   */
  su2double GetTemperature_Inf(void);

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
   * \brief Set the velocity at infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \param[in] val_velocity - Value of the velocity.
   */
  void SetVelocity_Inf(unsigned short val_dim, su2double val_velocity);

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
   * \brief Source term integration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);
  
  /*!
   * \brief Compute a pressure sensor switch.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config);
  
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
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute the undivided laplacian for the solution, except the energy equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute the max eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  void Evaluate_ObjFunc(CConfig *config);

  /*!
   * \author: G.Gori, S.Vitale, M.Pini, A.Guardone, P.Colonna
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
   * \brief Impose a custom or verification boundary condition.
   * \param[in] geometry         - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics    - Description of the convective numerical method.
   * \param[in] visc_numerics    - Description of the viscous numerical method.
   * \param[in] config           - Definition of the particular problem.
   * \param[in] val_marker       - Surface marker where the boundary condition is applied.
   */
  void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
   * \brief Impose the interface state across sliding meshes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
   void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config);

  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                   CNumerics *numerics, CConfig *config);
  
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
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Pressure_Forces(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Momentum_Forces(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the solution using an implicit Euler scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCLift_Inv(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CL(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF(unsigned short val_marker);
  
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
  su2double GetSurface_CL_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Inv(unsigned short val_marker);
  
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
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CL_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CEff_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFx_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFy_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFz_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMx_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMy_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMz_Mnt(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional drag coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCD_Inv(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCSF_Inv(unsigned short val_marker);
  
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
  su2double GetTotal_CSF(void);
  
  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CEff(void);
  
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
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CLift - Value of the total lift coefficient.
   */
  void SetTotal_CLift(su2double val_Total_CLift);

  /*!
   * \brief Set the value of the custom objective function.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  void SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight);

  /*!
   * \brief Add the value of the custom objective function.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  void AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CL(void);

  /*!
   * \author H. Kline
   * \brief Set the total "combo" objective (weighted sum of other values).
   * \param[in] ComboObj - Value of the combined objective.
   */
  void SetTotal_ComboObj(su2double ComboObj);

  /*!
   * \author H. Kline
   * \brief Provide the total "combo" objective (weighted sum of other values).
   * \return Value of the "combo" objective values.
   */
  su2double GetTotal_ComboObj(void);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CD(void);
  
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
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CoPx(void);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CoPy(void);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CoPz(void);

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
   * \brief Provide the total custom objective function.
   * \return Value of the custom objective function.
   */
  su2double GetTotal_Custom_ObjFunc(void);

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CDrag - Value of the total drag coefficient.
   */
  void SetTotal_CD(su2double val_Total_CDrag);
  
  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Inv(void);
  
  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Inv(void);
  
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
  su2double GetAllBound_CoPx_Inv(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPy_Inv(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPz_Inv(void);

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
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CEff_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMx_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMy_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMz_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPx_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPy_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPz_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFx_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFy_Mnt(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFz_Mnt(void);

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
   * \brief Set the value of the max residual and BGS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);
  
  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
  
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

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  void SetDensity_Inf(su2double rho_inf);

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(CConfig *config);

  /*!
   * \brief Update the Beta parameter for the incompressible preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   */
  void SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                       CConfig *config, unsigned short iMesh);

  /*!
   * \brief Compute the preconditioner for low-Mach flows.
   * \param[in] iPoint - Index of the grid point
   * \param[in] config - Definition of the particular problem.
   */
  void SetPreconditioner(CConfig *config, unsigned long iPoint);
  
  /*!
   * \brief Value of the total temperature at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total temperature is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is evaluated.
   * \return Value of the total temperature
   */
  su2double GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the total pressure at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total pressure is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is evaluated.
   * \return Value of the total pressure
   */
  su2double GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief A component of the unit vector representing the flow direction at an inlet boundary.
   * \param[in] val_marker - Surface marker where the flow direction is evaluated
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is evaluated
   * \param[in] val_dim - The component of the flow direction unit vector to be evaluated
   * \return Component of a unit vector representing the flow direction.
   */
  su2double GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim);
  
  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(CConfig* config, unsigned short iMarker);
  
  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex);
  
  /*!
   * \brief Get the set of value imposed at an inlet.
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(su2double *val_inlet,
                             unsigned long val_inlet_point,
                             unsigned short val_kind_marker,
                             string val_marker,
                             CGeometry *geometry,
                             CConfig *config);
  
  /*!
   * \brief A virtual member.
   */
  void GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output);

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  void SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component);

  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  void SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value);

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Get the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   */
   su2double GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index);

  /*!
   * \brief Compute the global error measures (L2, Linf) for verification cases.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVerificationError(CGeometry *geometry, CConfig *config);
};

/*!
 * \class CNSSolver
 * \brief Main class for defining the Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 */
class CNSSolver : public CEulerSolver {
private:
  su2double Viscosity_Inf;  /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;  /*!< \brief Turbulent kinetic energy at the infinity. */
  su2double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb;         /*!< \brief Turbulent Prandtl number. */
  su2double *CD_Visc,  /*!< \brief Drag coefficient (viscous contribution) for each boundary. */
  *CL_Visc,    /*!< \brief Lift coefficient (viscous contribution) for each boundary. */
  *CSF_Visc,    /*!< \brief Side force coefficient (viscous contribution) for each boundary. */
  *CMx_Visc,      /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CMy_Visc,      /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CMz_Visc,      /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CoPx_Visc,      /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CoPy_Visc,      /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CoPz_Visc,      /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CFx_Visc,      /*!< \brief Force x coefficient (viscous contribution) for each boundary. */
  *CFy_Visc,      /*!< \brief Force y coefficient (viscous contribution) for each boundary. */
  *CFz_Visc,      /*!< \brief Force z coefficient (viscous contribution) for each boundary. */
  *Surface_CL_Visc,/*!< \brief Lift coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CD_Visc,/*!< \brief Drag coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CSF_Visc,/*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CEff_Visc,/*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFx_Visc,  /*!< \brief Force x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFy_Visc,  /*!< \brief Force y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFz_Visc,  /*!< \brief Force z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMx_Visc,  /*!< \brief Moment x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMy_Visc,  /*!< \brief Moment y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMz_Visc,  /*!< \brief Moment z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_Buffet_Metric,  /*!< \brief Integrated separation sensor for each monitoring surface. */
  *CEff_Visc,      /*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
  *CMerit_Visc,      /*!< \brief Rotor Figure of Merit (Viscous contribution) for each boundary. */
  *Buffet_Metric,    /*!< \brief Integrated separation sensor for each boundary. */
  *CT_Visc,    /*!< \brief Thrust coefficient (viscous contribution) for each boundary. */
  *CQ_Visc,    /*!< \brief Torque coefficient (viscous contribution) for each boundary. */
  *HF_Visc,    /*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHF_Visc, /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
  ***HeatConjugateVar,   /*!< \brief Conjugate heat transfer variables for each boundary and vertex. */
  ***CSkinFriction,  /*!< \brief Skin friction coefficient for each boundary and vertex. */
  **Buffet_Sensor;   /*!< \brief Separation sensor for each boundary and vertex. */
  su2double Total_Buffet_Metric;  /*!< \brief Integrated separation sensor for all the boundaries. */
  su2double *ForceViscous,  /*!< \brief Viscous force for each boundary. */
  *MomentViscous;      /*!< \brief Inviscid moment for each boundary. */
  su2double AllBound_CD_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
  AllBound_CL_Visc,    /*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
  AllBound_CSF_Visc,    /*!< \brief Sideforce coefficient (viscous contribution) for all the boundaries. */
  AllBound_CMx_Visc,      /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Visc,      /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Visc,      /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Visc,      /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Visc,      /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Visc,      /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Visc,      /*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CFx_Visc,      /*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Visc,      /*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Visc,      /*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Visc,      /*!< \brief Rotor Figure of Merit coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CT_Visc,    /*!< \brief Thrust coefficient (viscous contribution) for all the boundaries. */
  AllBound_CQ_Visc,    /*!< \brief Torque coefficient (viscous contribution) for all the boundaries. */
  AllBound_HF_Visc,    /*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHF_Visc; /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
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
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CL_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CEff_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFx_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFy_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFz_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMx_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMy_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMz_Visc(unsigned short val_marker);
    
  /*!
   * \brief Provide the buffet metric.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the buffet metric on the surface <i>val_marker</i>.
   */
  su2double GetSurface_Buffet_Metric(unsigned short val_marker);
  
  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CEff_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMx_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMy_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMz_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPx_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPy_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPz_Visc(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFx_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFy_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFz_Visc(void);
    
  /*!
   * \brief Get the buffet metric.
   * \return Value of the buffet metric.
   */
  su2double GetTotal_Buffet_Metric(void);
  
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
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);
    
  /*!
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  void Evaluate_ObjFunc(CConfig *config);
  
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
   * \brief Impose the Navier-Stokes boundary condition (strong) with values from a CHT coupling.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   */
  su2double GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   * \param[in] relaxation factor - relaxation factor for the change of the variables
   * \param[in] val_var           - value of the variable
   */
  void SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var);
  
  /*!
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(CGeometry *geometry, CConfig *config);
    
  /*!
   * \brief Compute the buffet sensor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Buffet_Monitoring(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Get the total heat flux.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the integrated heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetSurface_HF_Visc(unsigned short val_marker);
  
  /*!
   * \brief Get the maximum (per surface) heat flux.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the maximum heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetSurface_MaxHF_Visc(unsigned short val_marker);
  
  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCL_Visc(unsigned short val_marker);

  /*!
   * \brief Get the non dimensional sideforce coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCSF_Visc(unsigned short val_marker);
  
  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCD_Visc(unsigned short val_marker);
  
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
  su2double GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim);
  
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
   * \brief Get the value of the buffet sensor
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the buffet sensor.
   */
  su2double GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex);
  
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
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRoe_Dissipation(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Computes the wall shear stress (Tau_Wall) on the surface using a wall function.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetTauWall_WF(CGeometry *geometry, CSolver** solver_container, CConfig* config);
 
};

/*!
 * \class CIncNSSolver
 * \brief Main class for defining the incompressible Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncNSSolver : public CIncEulerSolver {
private:
  su2double Viscosity_Inf;  /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;  /*!< \brief Turbulent kinetic energy at the infinity. */
  su2double *CD_Visc,  /*!< \brief Drag coefficient (viscous contribution) for each boundary. */
  *CL_Visc,    /*!< \brief Lift coefficient (viscous contribution) for each boundary. */
  *CSF_Visc,    /*!< \brief Side force coefficient (viscous contribution) for each boundary. */
  *CMx_Visc,      /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CMy_Visc,      /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CMz_Visc,      /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CoPx_Visc,      /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CoPy_Visc,      /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CoPz_Visc,      /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CFx_Visc,      /*!< \brief Force x coefficient (viscous contribution) for each boundary. */
  *CFy_Visc,      /*!< \brief Force y coefficient (viscous contribution) for each boundary. */
  *CFz_Visc,      /*!< \brief Force z coefficient (viscous contribution) for each boundary. */
  *Surface_CL_Visc,/*!< \brief Lift coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CD_Visc,/*!< \brief Drag coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CSF_Visc,/*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CEff_Visc,/*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFx_Visc,  /*!< \brief Force x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFy_Visc,  /*!< \brief Force y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFz_Visc,  /*!< \brief Force z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMx_Visc,  /*!< \brief Moment x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMy_Visc,  /*!< \brief Moment y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMz_Visc,  /*!< \brief Moment z coefficient (viscous contribution) for each monitoring surface. */
  *CEff_Visc,      /*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
  *CMerit_Visc,      /*!< \brief Rotor Figure of Merit (Viscous contribution) for each boundary. */
  *CT_Visc,    /*!< \brief Thrust coefficient (viscous contribution) for each boundary. */
  *CQ_Visc,    /*!< \brief Torque coefficient (viscous contribution) for each boundary. */
  *HF_Visc,    /*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHF_Visc, /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
  ***HeatConjugateVar,   /*!< \brief Conjugate heat transfer variables for each boundary and vertex. */
  ***CSkinFriction;  /*!< \brief Skin friction coefficient for each boundary and vertex. */
  su2double *ForceViscous,  /*!< \brief Viscous force for each boundary. */
  *MomentViscous;      /*!< \brief Inviscid moment for each boundary. */
  su2double AllBound_CD_Visc, /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
  AllBound_CL_Visc,    /*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
  AllBound_CSF_Visc,    /*!< \brief Sideforce coefficient (viscous contribution) for all the boundaries. */
  AllBound_CMx_Visc,      /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Visc,      /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Visc,      /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Visc,      /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Visc,      /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Visc,      /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Visc,      /*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CFx_Visc,      /*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Visc,      /*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Visc,      /*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Visc,      /*!< \brief Rotor Figure of Merit coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CT_Visc,    /*!< \brief Thrust coefficient (viscous contribution) for all the boundaries. */
  AllBound_CQ_Visc,    /*!< \brief Torque coefficient (viscous contribution) for all the boundaries. */
  AllBound_HF_Visc,    /*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHF_Visc; /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
  su2double StrainMag_Max, Omega_Max; /*!< \brief Maximum Strain Rate magnitude and Omega. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CIncNSSolver(void);
    
  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);
    
  /*!
   * \brief Destructor of the class.
   */
  ~CIncNSSolver(void);
    
  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CL_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Visc(unsigned short val_marker);
    
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Visc(unsigned short val_marker);
    
  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CEff_Visc(unsigned short val_marker);
    
    /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFx_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFy_Visc(unsigned short val_marker);
    
  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CFz_Visc(unsigned short val_marker);
    
  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMx_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMy_Visc(unsigned short val_marker);
  
  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CMz_Visc(unsigned short val_marker);
  
  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CEff_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMx_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMy_Visc(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CMz_Visc(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPx_Visc(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPy_Visc(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CoPz_Visc(void);

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFx_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFy_Visc(void);
  
  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  su2double GetAllBound_CFz_Visc(void);
  
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
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);
  
  /*!
   * \brief Impose a no-slip condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose an isothermal temperature condition at the wall.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the (received) conjugate heat variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   * \param[in] relaxation factor - relaxation factor for the change of the variables
   * \param[in] val_var           - value of the variable
   */
  void SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   */
  su2double GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var);

  /*!
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Get the total heat flux.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the integrated heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetSurface_HF_Visc(unsigned short val_marker);
  
  /*!
   * \brief Get the maximum (per surface) heat flux.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the maximum heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetSurface_MaxHF_Visc(unsigned short val_marker);

  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCL_Visc(unsigned short val_marker);

  /*!
   * \brief Get the non dimensional sideforce coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the sideforce coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCSF_Visc(unsigned short val_marker);
  
  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCD_Visc(unsigned short val_marker);
  
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
  su2double GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim);
  
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
 */
class CTurbSolver : public CSolver {
protected:
  su2double *FlowPrimVar_i,  /*!< \brief Store the flow solution at point i. */
  *FlowPrimVar_j,         /*!< \brief Store the flow solution at point j. */
  *lowerlimit,            /*!< \brief contains lower limits for turbulence variables. */
  *upperlimit;            /*!< \brief contains upper limits for turbulence variables. */
  su2double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */
  su2double*** Inlet_TurbVars; /*!< \brief Turbulence variables at inlet profiles */
  unsigned long nMarker, /*!< \brief Total number of markers using the grid information. */
  *nVertex;              /*!< \brief Store nVertex at each marker for deallocation */
  
  /* Sliding meshes variables */

  su2double ****SlidingState;
  int **SlidingStateNodes;

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
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSolver(CGeometry* geometry, CConfig *config);
  
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
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                     unsigned short val_marker);

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                     unsigned short val_marker);

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Giles(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                     unsigned short val_marker);

  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                   CNumerics *numerics, CConfig *config);
 
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
  
  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);
  
 /*!
  * \brief Get the outer state for fluid interface nodes.
  * \param[in] val_marker - marker index
  * \param[in] val_vertex - vertex index
  * \param[in] val_state  - requested state component
  */
  su2double GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index);

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex);
      
  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  void SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component);

  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  void SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value);

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex);

  /*!
   * \brief Set custom turbulence variables at the vertex of an inlet.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the turbulence variable (i.e. k is 0 in SST)
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  void SetInlet_TurbVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_turb_var);
};

/*!
 * \class CTurbSASolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSASolver: public CTurbSolver {
private:
  su2double nu_tilde_Inf, nu_tilde_Engine, nu_tilde_ActDisk;
  
public:
  /*!
   * \brief Constructor of the class.
   */
  CTurbSASolver(void);
  
  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] FluidModel
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh);

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
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
  void BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
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
   * \brief Impose the fluid interface boundary condition using tranfer data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config);

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
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                        CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose an actuator disk outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                         CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker, bool val_inlet_surface);

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] solver - Solver container
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  void SetDES_LengthScale(CSolver** solver, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the set of value imposed at an inlet.
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(su2double *val_inlet,
                             unsigned long val_inlet_point,
                             unsigned short val_kind_marker,
                             string val_marker,
                             CGeometry *geometry,
                             CConfig *config);

  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(CConfig* config, unsigned short iMarker);

  /*!
   * \brief Get the value of nu tilde at the far-field.
   * \return Value of nu tilde at the far-field.
   */
  su2double GetNuTilde_Inf(void);

  /*!
   * \brief Compute nu tilde from the wall functions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void SetNuTilde_WF(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);
};

/*!
 * \class CTurbSSTSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Campos, F. Palacios, T. Economon
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
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
  void BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
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
  * \brief Impose the interface state across sliding meshes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] solver_container - Container vector with all the solutions.
  * \param[in] conv_numerics - Description of the numerical method.
  * \param[in] visc_numerics - Description of the numerical method.
  * \param[in] config - Definition of the particular problem.
  */
  void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config);

  /*!
   * \brief Get the constants for the SST model.
   * \return A pointer to an array containing a set of constants
   */
  su2double* GetConstants();
  
  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(CConfig *config);

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the set of value imposed at an inlet.
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(su2double *val_inlet,
                             unsigned long val_inlet_point,
                             unsigned short val_kind_marker,
                             string val_marker,
                             CGeometry *geometry,
                             CConfig *config);
  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(CConfig* config, unsigned short iMarker);

  /*!
   * \brief Get the value of the turbulent kinetic energy.
   * \return Value of the turbulent kinetic energy.
   */
  su2double GetTke_Inf(void);

  /*!
   * \brief Get the value of the turbulent frequency.
   * \return Value of the turbulent frequency.
   */
  su2double GetOmega_Inf(void);

};

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Aranake.
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] second_numerics - Description of the second numerical method.
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
  CSysMatrix<su2double> JacobianItmc; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
  su2double *LinSysSolItmc;    /*!< \brief vector to store iterative solution of implicit linear system. */
  su2double *LinSysResItmc;    /*!< \brief vector to store iterative residual of implicit linear system. */
  su2double *rhsItmc;    /*!< \brief right hand side of implicit linear system. */
  CSysMatrix<su2double> JacobianReth; /*!< \brief Complete sparse Jacobian structure for implicit computations. */
  su2double *LinSysSolReth;    /*!< \brief vector to store iterative solution of implicit linear system. */
  su2double *LinSysResReth;    /*!< \brief vector to store iterative residual of implicit linear system. */
  su2double *rhsReth;    /*!< \brief right hand side of implicit linear system. */
};

/*!
 * \class CAdjEulerSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 */
class CAdjEulerSolver : public CSolver {
protected:
  su2double PsiRho_Inf,  /*!< \brief PsiRho variable at the infinity. */
  PsiE_Inf,      /*!< \brief PsiE variable at the infinity. */
  *Phi_Inf;      /*!< \brief Phi vector at the infinity. */
  su2double *Sens_Mach, /*!< \brief Mach sensitivity coefficient for each boundary. */
  *Sens_AoA,      /*!< \brief Angle of attack sensitivity coefficient for each boundary. */
  *Sens_Geo,      /*!< \brief Shape sensitivity coefficient for each boundary. */
  *Sens_Press,      /*!< \brief Pressure sensitivity coefficient for each boundary. */
  *Sens_Temp,      /*!< \brief Temperature sensitivity coefficient for each boundary. */
  *Sens_BPress,     /*!< \brief Back pressure sensitivity coefficient for each boundary. */
  **CSensitivity,    /*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
  ***DonorAdjVar;    /*!< \brief Value of the donor variables at each boundary. */
  su2double Total_Sens_Mach;  /*!< \brief Total mach sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_AoA;    /*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Geo;    /*!< \brief Total shape sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
  su2double Total_Sens_Temp;    /*!< \brief Total farfield sensitivity to temperature. */
  su2double Total_Sens_BPress;    /*!< \brief Total sensitivity to back pressure. */
  bool space_centered;  /*!< \brief True if space centered scheeme used. */
  su2double **Jacobian_Axisymmetric; /*!< \brief Storage for axisymmetric Jacobian. */
  unsigned long nMarker;        /*!< \brief Total number of markers using the grid information. */
  su2double Gamma;                  /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;        /*!< \brief Fluids's Gamma - 1.0  . */
  su2double *FlowPrimVar_i,  /*!< \brief Store the flow solution at point i. */
  *FlowPrimVar_j;        /*!< \brief Store the flow solution at point j. */
  unsigned long **DonorGlobalIndex;    /*!< \brief Value of the donor global index. */

  su2double pnorm,
  Area_Monitored; /*!< \brief Store the total area of the monitored outflow surface (used for normalization in continuous adjoint outflow conditions) */
  
  unsigned long AoA_Counter;
  su2double ACoeff, ACoeff_inc, ACoeff_old;
  bool Update_ACoeff;
  
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
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Nearfield(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Interface(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config);

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
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  su2double *GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value);
  
  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  su2double GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  unsigned long GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex);
  
  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  void SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index);
  
  /*!
   * \brief Compute the sensor for higher order dissipation control in rotating problems.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config);
  
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
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                        CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose an actuator disk outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                         CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                  CConfig *config, unsigned short val_marker, bool val_inlet_surface);

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
   * \param[in] Output - boolean to determine whether to print output.
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
   * \brief Get the total Back pressure number sensitivity coefficient.
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

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

};

/*!
 * \class CAdjNSSolver
 * \brief Main class for defining the Navier-Stokes' adjoint flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
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
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] second_numerics - Description of the second numerical method.
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
 */
class CAdjTurbSolver : public CSolver {
private:
  su2double PsiNu_Inf,  /*!< \brief PsiNu variable at the infinity. */
  *FlowSolution_i,  /*!< \brief Store the flow solution at point i. */
  *FlowSolution_j;  /*!< \brief Store the flow solution at point j. */
  
  su2double Gamma;                  /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;        /*!< \brief Fluids's Gamma - 1.0  . */
  
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
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] second_numerics - Description of the second numerical method.
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
 *  \date May 3, 2010.
 */
class CPoissonSolver : public CSolver {
private:
  su2double *Source_Vector;      /*!< \brief Auxiliary vector for storing element source vector. */
  su2double **StiffMatrix_Elem; /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
  su2double **StiffMatrix_Node;  /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
  
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
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
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
 *  \date May 3, 2010.
 */
class CWaveSolver : public CSolver {
private:
  su2double *CWave;  /*!< \brief Wave strength for each boundary. */
  su2double AllBound_CWave;  /*!< \brief Total wave strength for all the boundaries. */
  su2double Total_CWave; /*!< \brief Total wave strength for all the boundaries. */
  
  CSysMatrix<su2double> StiffMatrixSpace; /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
  CSysMatrix<su2double> StiffMatrixTime;  /*!< \brief Sparse structure for storing the stiffness matrix in Galerkin computations. */
  
  su2double **StiffMatrix_Elem,      /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
  **StiffMatrix_Node;              /*!< \brief Auxiliary matrices for storing point to point Stiffness Matrices. */
  
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
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
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

/*! \class CHeatSolverFVM
 *  \brief Main class for defining the finite-volume heat solver.
 *  \author O. Burghardt
 *  \date January 19, 2018.
 */
class CHeatSolverFVM : public CSolver {
protected:
  unsigned short nVarFlow, nMarker, CurrentMesh;
  su2double *Heat_Flux, *Surface_HF, Total_HeatFlux, AllBound_HeatFlux,
            *AvgTemperature, Total_AvgTemperature, AllBound_AvgTemperature,
            *Primitive, *Primitive_Flow_i, *Primitive_Flow_j,
            *Surface_Areas, Total_HeatFlux_Areas, Total_HeatFlux_Areas_Monitor;
  su2double ***ConjugateVar, ***InterfaceVar;

public:

  /*!
   * \brief Constructor of the class.
   */
  CHeatSolverFVM(void);

  /*!
   * \brief Constructor of the class.
   */
  CHeatSolverFVM(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CHeatSolverFVM(void);

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                     unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

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
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                         CConfig *config, unsigned short iMesh);

  /*!
   * \brief Compute the undivided laplacian for the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config);

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


  void Set_Heatflux_Areas(CGeometry *geometry, CConfig *config);

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
   * \brief Impose a constant heat-flux condition at the wall.
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
  void BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the (received) conjugate heat variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                          unsigned short val_marker);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   */
  su2double GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var);

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   * \param[in] relaxation factor - relaxation factor for the change of the variables
   * \param[in] val_var           - value of the variable
   */
  void SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var);

  /*!
   * \brief Evaluate heat-flux related objectives.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Heat_Fluxes(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief Get value of the heat load (integrated heat flux).
   * \return Value of the heat load (integrated heat flux).
   */
  su2double GetTotal_HeatFlux(void);

  /*!
   * \brief Get value of the integral-averaged temperature.
   * \return Value of the integral-averaged temperature.
   */
  su2double GetTotal_AvgTemperature(void);

  /*!
   * \brief Update the solution using an implicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief Update the solution using an explicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh, unsigned long Iteration);

  /*!
   * \brief Set the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);

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
   * \brief Set the value of the max residual and BGS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);
};

/*! \class CFEASolver
 *  \brief Main class for defining a FEM solver for elastic structural problems.
 *  \author R. Sanchez.
 *  \date July 10, 2015.
 */
class CFEASolver : public CSolver {
private:
  
  su2double  Total_CFEA;        /*!< \brief Total FEA coefficient for all the boundaries. */
  /*!< We maintain the name to avoid defining a new function... */
  
  unsigned long nElement;
  unsigned short nMarker;
  int nFEA_Terms;
  
  bool element_based;             /*!< \brief Bool to determine if an element-based file is used. */

  su2double *GradN_X,
  *GradN_x;
  
  su2double **Jacobian_c_ij;      /*!< \brief Submatrix to store the constitutive term for node ij. */
  su2double **Jacobian_s_ij;      /*!< \brief Submatrix to store the stress contribution of node ij (diagonal). */
  su2double **MassMatrix_ij;      /*!< \brief Submatrix to store the term ij of the mass matrix. */
  su2double *Res_Stress_i;      /*!< \brief Submatrix to store the nodal stress contribution of node i. */
  
  su2double *Res_Ext_Surf;      /*!< \brief Auxiliary vector to store the surface load contribution to the residual */
  su2double *Res_Time_Cont;      /*!< \brief Auxiliary vector to store the surface load contribution to the residual */
  su2double *Res_FSI_Cont;      /*!< \brief Auxiliary vector to store the surface load contribution to the residual */
  
  su2double *Res_Dead_Load;      /*!< \brief Auxiliary vector to store the body load contribution to the residual */
  
  su2double *solutionPredictor;    /*!< \brief Auxiliary vector to store the solution predictor */
  
  su2double *Solution_Interm;    /*!< \brief Auxiliary vector to store the intermediate solution */
  
  su2double *SolRest;      /*!< \brief Auxiliary vector to restart the solution */
  
  su2double *nodeReactions;      /*!< \brief Auxiliary vector to store the reactions */
  
  su2double *normalVertex;      /*!< \brief Auxiliary vector to store the normals to a certain vertex */
  su2double **stressTensor;      /*!< \brief Auxiliary matrix to rebuild the stress tensor and compute reactions */
  
  su2double **mZeros_Aux;      /*!< \brief Submatrix to make zeros and impose clamped boundary conditions. */
  su2double **mId_Aux;        /*!< \brief Diagonal submatrix to impose clamped boundary conditions. */
  
  unsigned long *elProperties;  /*!< \brief Auxiliary vector to read the element properties from file */

  unsigned short *iElem_iDe;			/*!< \brief For DE cases, ID of the region considered for each iElem. */
  
  su2double a_dt[9];          /*!< \brief Integration constants. */
  
  su2double Conv_Ref[3];        /*!< \brief Reference values for convergence check: DTOL, RTOL, ETOL */
  su2double Conv_Check[3];      /*!< \brief Current values for convergence check: DTOL, RTOL, ETOL */
  su2double FSI_Conv[2];        /*!< \brief Values to check the convergence of the FSI problem. */
  
  su2double loadIncrement;      /*!< \brief Coefficient that determines the amount of load which is applied */
  
  su2double WAitken_Dyn;        /*!< \brief Aitken's dynamic coefficient */
  su2double WAitken_Dyn_tn1;      /*!< \brief Aitken's dynamic coefficient in the previous iteration */
  
  su2double PenaltyValue;      /*!< \brief Penalty value to maintain total stiffness constant */

  su2double Total_OFRefGeom;        /*!< \brief Total Objective Function: Reference Geometry. */
  su2double Total_OFRefNode;        /*!< \brief Total Objective Function: Reference Node. */
  su2double Total_OFVolFrac;        /*!< \brief Total Objective Function: Volume fraction (topology optimization). */

  su2double Global_OFRefGeom;        /*!< \brief Global Objective Function (added over time steps): Reference Geometry. */
  su2double Global_OFRefNode;        /*!< \brief Global Objective Function (added over time steps): Reference Node. */

  su2double Total_ForwardGradient;  /*!< \brief Vector of the total forward gradient. */
  
  su2double ForceCoeff;             /*!< \brief Load transfer coefficient . */
  su2double RelaxCoeff;             /*!< \brief Relaxation coefficient . */
  su2double FSI_Residual;           /*!< \brief FSI residual. */

public:
  
  CSysVector<su2double> TimeRes_Aux;      /*!< \brief Auxiliary vector for adding mass and damping contributions to the residual. */
  CSysVector<su2double> TimeRes;        /*!< \brief Vector for adding mass and damping contributions to the residual */
  CSysVector<su2double> LinSysReact;      /*!< \brief Vector to store the residual before applying the BCs */

  CSysVector<su2double> LinSysSol_Adj;   /*!< \brief Vector to store the solution of the adjoint problem */
  CSysVector<su2double> LinSysRes_Adj;   /*!< \brief Vector to store the residual of the adjoint problem */

  CSysMatrix<su2double> MassMatrix;       /*!< \brief Sparse structure for storing the mass matrix. */

  CElement*** element_container;   /*!< \brief Vector which the define the finite element structure for each problem. */
  CElementProperty** element_properties; /*!< \brief Vector which stores the properties of each element */

  
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
  virtual ~CFEASolver(void);

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
   * \brief Set the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
  
  /*!
   * \brief Reset the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);
  
  /*!
   * \brief Compute the time step for solving the FEM equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                    unsigned short iMesh, unsigned long Iteration);
  
  
  /*!
   * \brief Set a reference geometry for .
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_ReferenceGeometry(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute the stiffness matrix of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the stiffness matrix of the problem and the nodal stress terms at the same time (more efficient if full Newton Raphson).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the mass matrix of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the mass residual of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_MassRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);

  /*!
   * \brief Compute the nodal stress terms and add them to the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the stress at the nodes for output purposes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  
  void Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the dead loads.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_DeadLoad(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Initializes the matrices/residuals in the solution process (avoids adding over previous values).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_IntegrationConstants(CConfig *config);
  
  /*!
   * \brief Clamped boundary conditions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Enforce the solution to be 0 in the clamped nodes - Avoids accumulation of numerical error.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                       unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  
  void BC_DispDir(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
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
   * \brief Impose a load boundary condition normal to the boundary.
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
   * \brief Impose a damping load.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Damper(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                 unsigned short val_marker);

  /*!
   * \brief Required step for non conservative interpolation schemes where stresses are transferred instead of forces.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integrate_FSI_Loads(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the solution using an implicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Iterate using an implicit Newmark solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Update the solution using an implicit Newmark solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Iterate using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Update the solution using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Update the solution using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Update the solution using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Postprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                      unsigned short iMesh);
  
  /*!
   * \brief Routine to solve the Jacobian-Residual linearized system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Get the residual for FEM structural analysis.
   * \param[in] val_var - Index of the variable.
   * \return Value of the residual for the variable in the position <i>val_var</i>.
   */
  su2double GetRes_FEM(unsigned short val_var);
  
  /*!
   * \brief Provide the maximum Von Mises Stress for structural analysis.
   * \return Value of the maximum Von Mises Stress.
   */
  su2double GetTotal_CFEA(void);
  
  /*!
   * \brief Retrieve the value of the objective function for a reference geometry
   * \param[out] OFRefGeom - value of the objective function.
   */
  su2double GetTotal_OFRefGeom(void);
  
  /*!
   * \brief Retrieve the value of the objective function for a reference node
   * \param[out] OFRefNode - value of the objective function.
   */
  su2double GetTotal_OFRefNode(void);
  
  /*!
   * \brief Retrieve the value of the volume fraction objective function
   * \param[out] OFVolFrac - value of the objective function.
   */
  su2double GetTotal_OFVolFrac(void);

  /*!
   * \brief Determines whether there is an element-based file or not.
   * \return Bool that defines whether the solution has an element-based file or not
   */
  bool IsElementBased(void);

  /*!
   * \brief Set the value of the FEA coefficient.
   * \param[in] val_cfea - Value of the FEA coefficient.
   */
  void SetTotal_CFEA(su2double val_cfea);
  
  /*!
   * \brief Set the value of the objective function for a reference geometry.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference geometry.
   */
  void SetTotal_OFRefGeom(su2double val_ofrefgeom);
  
  /*!
   * \brief Set the value of the objective function for a reference node.
   * \param[in] val_ofrefnode - Value of the objective function for a reference node.
   */
  void SetTotal_OFRefNode(su2double val_ofrefnode);

  /*!
   * \brief Set the value of the force coefficient history for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_forcecoeff_history - Value of the force coefficient.
   */
  void SetForceCoeff(su2double val_forcecoeff_history);

  /*!
   * \brief Set the value of the FSI residual for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  void SetFSI_Residual(su2double val_FSI_residual);

  /*!
   * \brief Set the value of the FSI residual for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  void SetRelaxCoeff(su2double val_relaxcoeff_history);

  /*!
   * \brief Get the value of the force coefficient history for the history file.
   * \param[out] val_forcecoeff_history - Value of the force coefficient.
   */
  su2double GetForceCoeff(void);

  /*!
   * \brief Get the value of the relaxation coefficient history for the history file.
   * \param[out] val_relaxcoeff_history - Value of the relaxation coefficient.
   */
  su2double GetRelaxCoeff(void);

  /*!
   * \brief Get the value of the FSI residual for the history file.
   * \param[out] val_FSI_residual - Value of the residual.
   */
  su2double GetFSI_Residual(void);
  
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
                                 unsigned long iOuterIter);
  
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
   * \brief Compute the objective function for a reference geometry
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Compute the objective function for a reference node
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config);
  
  /*!
   * \brief Compute the objective function for a volume fraction
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFVolFrac(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief Compute the penalty due to the stiffness increase
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Stiffness_Penalty(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config);

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
  
  /*!
   * \brief Set the value of the load increment for nonlinear structural analysis
   * \param[in] Value of the coefficient
   */
  void SetLoad_Increment(su2double val_loadIncrement);
  
  /*!
   * \brief Get the value of the load increment for nonlinear structural analysis
   * \param[in] Value of the coefficient
   */
  su2double GetLoad_Increment(void);

  /*!
   * \brief Set a reference geometry for prestretched conditions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_Prestretch(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Retrieve the iDe index for DE computations
   * \param[in] iElem - element parameter.
   * \param[out] iElem_iDe - ID of the Dielectric Elastomer region.
   */
  unsigned short Get_iElem_iDe(unsigned long iElem);
  
  /*!
   * \brief Retrieve the Mass Matrix term (to add to the Jacobian of the adjoint problem)
   * \param[in] iPoint - Point i of the Mass Matrix.
   * \param[in] jPoint - Point j of the Mass Matrix.
   * \param[in] iVar - Variable i of the Mass Matrix submatrix.
   * \param[in] iVar - Variable j of the Mass Matrix submatrix.
   */
  su2double Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar);
  
  /*!
   * \brief Set the value of the max residual and BGS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);
  
  /*!
   * \brief Set container of element properties.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_ElementProperties(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Get multiplier for loads.
   * \param[in] CurrentTime - Current time step.
   * \param[in] config - Definition of the particular problem.
   */
  su2double Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, CConfig *config);
  
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
  
  /*!
   * \brief Filter the density field for topology optimization applications
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void FilterElementDensities(CGeometry *geometry, CConfig *config);

};

/*!
 * \class CTemplateSolver
 * \brief Main class for defining the template model solver.
 * \ingroup Template_Flow_Equation
 * \author F. Palacios
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
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
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
   * \param[in] second_numerics - Description of the second numerical method.
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
   * \brief Impose a custom or verification boundary condition.
   * \param[in] geometry         - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics    - Description of the convective numerical method.
   * \param[in] visc_numerics    - Description of the viscous numerical method.
   * \param[in] config           - Definition of the particular problem.
   * \param[in] val_marker       - Surface marker where the boundary condition is applied.
   */
  void BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
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
 */
class CDiscAdjSolver : public CSolver {
private:
  unsigned short KindDirect_Solver;
  CSolver *direct_solver;
  su2double **CSensitivity;  /*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
  su2double Total_Sens_Mach;  /*!< \brief Total mach sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_AoA;    /*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Geo;    /*!< \brief Total shape sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
  su2double Total_Sens_Temp;    /*!< \brief Total farfield sensitivity to temperature. */
  su2double Total_Sens_BPress;    /*!< \brief Total sensitivity to outlet pressure. */
  su2double Total_Sens_Density;    /*!< \brief Total sensitivity to initial density (incompressible). */
  su2double Total_Sens_ModVel;    /*!< \brief Total sensitivity to inlet velocity (incompressible). */
  su2double ObjFunc_Value;        /*!< \brief Value of the objective function. */
  su2double Mach, Alpha, Beta, Pressure, Temperature, BPressure, ModVel;
  unsigned long nMarker;        /*!< \brief Total number of markers using the grid information. */
  
  su2double *Solution_Geometry; /*!< \brief Auxiliary vector for the geometry solution (dimension nDim instead of nVar). */
  
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
   * \brief Destructor of the class.
   */
  ~CDiscAdjSolver(void);
  
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
   * \brief Sets the adjoint values of the output of the mesh deformation iteration
   *        before evaluation of the tape.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void SetAdjoint_OutputMesh(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
   *        after tape has been evaluated.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_Geometry(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Sets the adjoint values of the flow variables due to cross term contributions
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm(CGeometry *geometry,  CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry,  CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm_Geometry_Flow(CGeometry *geometry,  CConfig *config);
  
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
   * \brief Get the total Back pressure number sensitivity coefficient.
   * \return Value of the Back sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_BPress(void);

  /*!
   * \brief Get the total density sensitivity coefficient.
   * \return Value of the density sensitivity.
   */
  su2double GetTotal_Sens_Density(void);

  /*!
   * \brief Get the total velocity magnitude sensitivity coefficient.
   * \return Value of the velocity magnitude sensitivity.
   */
  su2double GetTotal_Sens_ModVel(void);

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
  void SetRecording(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetMesh_Recording(CGeometry **geometry, CVolumetricMovement *grid_movement,
      CConfig *config);
  
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
  
  /*!
   * \brief Update the dual-time derivatives.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);
  
  /*!
   * \brief Set the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);
};

/*!
 * \class CDiscAdjFEASolver
 * \brief Main class for defining the discrete adjoint solver for FE structural problems.
 * \ingroup Discrete_Adjoint
 * \author R. Sanchez
 */
class CDiscAdjFEASolver : public CSolver {
private:
  unsigned short KindDirect_Solver;
  CSolver *direct_solver;
  su2double *Sens_E,            /*!< \brief Young modulus sensitivity coefficient for each boundary. */
  *Sens_Nu,                     /*!< \brief Poisson's ratio sensitivity coefficient for each boundary. */
  *Sens_nL,                     /*!< \brief Normal pressure sensitivity coefficient for each boundary. */
  **CSensitivity;               /*!< \brief Shape sensitivity coefficient for each boundary and vertex. */

  su2double *Solution_Vel,      /*!< \brief Velocity componenent of the solution. */
  *Solution_Accel;              /*!< \brief Acceleration componenent of the solution. */

  su2double *SolRest;     /*!< \brief Auxiliary vector to restart the solution */

  su2double ObjFunc_Value;      /*!< \brief Value of the objective function. */
  su2double *normalLoads;       /*!< \brief Values of the normal loads for each marker iMarker_nL. */
  unsigned long nMarker;        /*!< \brief Total number of markers using the grid information. */
  unsigned long nMarker_nL;     /*!< \brief Total number of markers that have a normal load applied. */

  /*!< \brief Definition of element based sensitivities. */
  unsigned short nMPROP;        /*!< \brief Number of material properties */

  su2double *E_i,               /*!< \brief Values of the Young's Modulus. */
            *Nu_i,              /*!< \brief Values of the Poisson's ratio. */
            *Rho_i,             /*!< \brief Values of the density (for inertial effects). */
            *Rho_DL_i;          /*!< \brief Values of the density (for volume loading) */
  su2double *Local_Sens_E, *Global_Sens_E,            /*!< \brief Local and global sensitivity of the Young's modulus. */
            *Local_Sens_Nu, *Global_Sens_Nu,          /*!< \brief Local and global sensitivity of the Poisson ratio. */
            *Local_Sens_Rho, *Global_Sens_Rho,        /*!< \brief Local and global sensitivity of the density. */
            *Local_Sens_Rho_DL, *Global_Sens_Rho_DL;  /*!< \brief Local and global sensitivity of the volume load. */
  su2double *Total_Sens_E,       /*!< \brief Total sensitivity of the Young's modulus (time domain). */
            *Total_Sens_Nu,      /*!< \brief Local and global sensitivity of the Poisson ratio. */
            *Total_Sens_Rho,     /*!< \brief Local and global sensitivity of the density. */
            *Total_Sens_Rho_DL;  /*!< \brief Local and global sensitivity of the volume load. */

  bool de_effects;              /*!< \brief Determines if DE effects are considered. */
  unsigned short nEField;       /*!< \brief Number of electric field areas in the problem. */
  su2double *EField;            /*!< \brief Array that stores the electric field as design variables. */
  su2double *Local_Sens_EField, *Global_Sens_EField; /*!< \brief Local and global sensitivity of the Electric Field. */
  su2double *Total_Sens_EField;

  bool fea_dv;                  /*!< \brief Determines if the design variable we study is a FEA parameter. */
  unsigned short nDV;           /*!< \brief Number of design variables in the problem. */
  su2double *DV_Val;            /*!< \brief Value of the design variables. */
  su2double *Local_Sens_DV, *Global_Sens_DV;          /*!< \brief Local and global sensitivity of the Design Variable. */
  su2double *Total_Sens_DV;

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjFEASolver(void);
  
  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CDiscAdjFEASolver(CGeometry *geometry, CConfig *config);
  
  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   */
  CDiscAdjFEASolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEASolver(void);
  
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
   * \brief Sets the adjoint values of the structural variables due to cross term contributions
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm(CGeometry *geometry,  CConfig *config);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry,  CConfig *config);
  
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
   * \brief Provide the total Young's modulus sensitivity
   * \return Value of the total Young's modulus sensitivity
   *         (inviscid + viscous contribution).
   */
  su2double GetTotal_Sens_E(unsigned short iVal);
  
  /*!
   * \brief Set the total Poisson's ratio sensitivity.
   * \return Value of the Poisson's ratio sensitivity
   */
  su2double GetTotal_Sens_Nu(unsigned short iVal);
  
  /*!
   * \brief Get the total sensitivity for the structural density
   * \return Value of the structural density sensitivity
   */
  su2double GetTotal_Sens_Rho(unsigned short iVal);

  /*!
   * \brief Get the total sensitivity for the structural weight
   * \return Value of the structural weight sensitivity
   */
  su2double GetTotal_Sens_Rho_DL(unsigned short iVal);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField (time averaged)
   */
  su2double GetTotal_Sens_EField(unsigned short iEField);

  /*!
   * \brief A virtual member.
   * \return Value of the total sensitivity coefficient for the FEA DV in the region iDVFEA (time averaged)
   */
  su2double GetTotal_Sens_DVFEA(unsigned short iDVFEA);

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Young Modulus E
   */
  su2double GetGlobal_Sens_E(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the Mach sensitivity for the Poisson's ratio Nu
   */
  su2double GetGlobal_Sens_Nu(unsigned short iVal);
  
  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField
   */
  su2double GetGlobal_Sens_EField(unsigned short iEField);
  
  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the FEA DV in the region iDVFEA
   */
  su2double GetGlobal_Sens_DVFEA(unsigned short iDVFEA);

  /*!
   * \brief Get the total sensitivity for the structural density
   * \return Value of the structural density sensitivity
   */
  su2double GetGlobal_Sens_Rho(unsigned short iVal);

  /*!
   * \brief Get the total sensitivity for the structural weight
   * \return Value of the structural weight sensitivity
   */
  su2double GetGlobal_Sens_Rho_DL(unsigned short iVal);


  /*!
   * \brief Get the value of the Young modulus from the adjoint solver
   * \return Value of the Young modulus from the adjoint solver
   */
  su2double GetVal_Young(unsigned short iVal);
  
  /*!
   * \brief Get the value of the Poisson's ratio from the adjoint solver
   * \return Value of the Poisson's ratio from the adjoint solver
   */
  su2double GetVal_Poisson(unsigned short iVal);
  
  /*!
   * \brief Get the value of the density from the adjoint solver, for inertial effects
   * \return Value of the density from the adjoint solver
   */
  su2double GetVal_Rho(unsigned short iVal);
  
  /*!
   * \brief Get the value of the density from the adjoint solver, for dead loads
   * \return Value of the density for dead loads, from the adjoint solver
   */
  su2double GetVal_Rho_DL(unsigned short iVal);
  
  /*!
   * \brief Get the number of variables for the Electric Field from the adjoint solver
   * \return Number of electric field variables from the adjoint solver
   */
  unsigned short GetnEField(void);
  
  /*!
   * \brief Read the design variables for the adjoint solver
   */
  void ReadDV(CConfig *config);

  /*!
   * \brief Get the number of design variables from the adjoint solver,
   * \return Number of design variables from the adjoint solver
   */
  unsigned short GetnDVFEA(void);

  /*!
   * \brief Get the value of the Electric Field from the adjoint solver
   * \return Pointer to the values of the Electric Field
   */
  su2double GetVal_EField(unsigned short iVal);
  
  /*!
   * \brief Get the value of the design variables from the adjoint solver
   * \return Pointer to the values of the design variables
   */
  su2double GetVal_DVFEA(unsigned short iVal);

  /*!
   * \brief Set the value of the max residual and RMS residual.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   * \param[in] val_iterlinsolver - Number of linear iterations.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetRecording(CGeometry *geometry, CConfig *config);
  
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
  
  /*!
   * \brief Update the dual-time derivatives.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
  
  /*!
   * \brief Enforce the solution to be 0 in the clamped nodes - Avoids accumulation of numerical error.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
      unsigned short val_marker);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

};

/*!
 * \class CFEM_DG_EulerSolver
 * \brief Main class for defining the Euler Discontinuous Galerkin finite element flow solver.
 * \ingroup Euler_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 6.2.0 "Falcon"
 */
class CFEM_DG_EulerSolver : public CSolver {
protected:

  unsigned long nMarker; /*!< \brief Total number of markers using the grid information. */

  su2double Gamma;           /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One; /*!< \brief Fluids's Gamma - 1.0  . */

  CFluidModel  *FluidModel; /*!< \brief fluid model used in the solver */

  su2double
  Mach_Inf,	       /*!< \brief Mach number at infinity. */
  Density_Inf,	       /*!< \brief Density at infinity. */
  Energy_Inf,	       /*!< \brief Energy at infinity. */
  Temperature_Inf,     /*!< \brief Energy at infinity. */
  Pressure_Inf,	       /*!< \brief Pressure at infinity. */
  *Velocity_Inf;       /*!< \brief Flow velocity vector at infinity. */

  vector<su2double> ConsVarFreeStream; /*!< \brief Vector, which contains the free stream
                                        conservative variables. */
  su2double
  *CL_Inv,        /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CD_Inv,        /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,   /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,          /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,          /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,          /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,          /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,          /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,          /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CEff_Inv;         /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */

  su2double
  *Surface_CL_Inv,      /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Inv,      /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,        /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,        /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,        /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,        /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,        /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,        /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv;       /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each monitoring surface. */

  su2double
  AllBound_CL_Inv, 	  /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CD_Inv,      /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Inv, /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Inv, 	    /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Inv, 	    /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Inv, 	    /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Inv, 	    /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Inv, 	    /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Inv, 	    /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Inv; 	  /*!< \brief Total efficiency (Cl/Cd) (inviscid contribution) for all the boundaries. */

  su2double
  Total_CL, 	  /*!< \brief Total lift coefficient for all the boundaries. */
  Total_CD,      /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CSF, /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CFx, 	   /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy, 	   /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz, 	   /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CMx, 	   /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMy, 	   /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMz, 	   /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CEff; 	   /*!< \brief Total efficiency coefficient for all the boundaries. */

  su2double
  *Surface_CL,      /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CD,      /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSF, /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,        /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,        /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,        /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,        /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,        /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz,        /*!< \brief z Moment coefficient for each monitoring surface. */
  *Surface_CEff;       /*!< \brief Efficiency (Cl/Cd) for each monitoring surface. */

  su2double Cauchy_Value,         /*!< \brief Summed value of the convergence indicator. */
  Cauchy_Func; 	                  /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter;  /*!< \brief Number of elements of the Cauchy serial. */
  su2double *Cauchy_Serie; 	  /*!< \brief Complete Cauchy serial. */
  su2double Old_Func,             /*!< \brief Old value of the objective function (the function which is monitored). */
  New_Func; 	                  /*!< \brief Current value of the objective function (the function which is monitored). */


  unsigned long nDOFsLocTot;    /*!< \brief Total number of local DOFs, including halos. */
  unsigned long nDOFsLocOwned;  /*!< \brief Number of owned local DOFs. */
  unsigned long nDOFsGlobal;    /*!< \brief Number of global DOFs. */

  unsigned long nVolElemTot;    /*!< \brief Total number of local volume elements, including halos. */
  unsigned long nVolElemOwned;  /*!< \brief Number of owned local volume elements. */
  CVolumeElementFEM *volElem;   /*!< \brief Array of the local volume elements, including halos. */

  const unsigned long *nVolElemOwnedPerTimeLevel;    /*!< \brief Number of owned local volume elements
                                                                 per time level. Cumulative storage. */
  const unsigned long *nVolElemInternalPerTimeLevel; /*!< \brief Number of internal local volume elements per
                                                                 time level. Internal means that the solution
                                                                 data does not need to be communicated. */
  const unsigned long *nVolElemHaloPerTimeLevel;     /*!< \brief Number of halo volume elements
                                                                 per time level. Cumulative storage. */

  vector<vector<unsigned long> > ownedElemAdjLowTimeLevel; /*!< \brief List of owned elements per time level that are
                                                                       adjacent to elements of the lower time level. */
  vector<vector<unsigned long> > haloElemAdjLowTimeLevel;  /*!< \brief List of halo elements per time level that are
                                                                       adjacent to elements of the lower time level. */

  unsigned long nMeshPoints;    /*!< \brief Number of mesh points in the local part of the grid. */
  CPointFEM *meshPoints;        /*!< \brief Array of the points of the FEM mesh. */

  const unsigned long *nMatchingInternalFacesWithHaloElem;  /*!< \brief Number of local matching internal faces per time level
                                                                        between an owned and a halo element. Cumulative storage. */
  const unsigned long *nMatchingInternalFacesLocalElem;     /*!< \brief Number of local matching internal faces per time level
                                                                        between local elements. Cumulative storage. */

  CInternalFaceElementFEM *matchingInternalFaces;    /*!< \brief Array of the local matching internal faces. */
  CBoundaryFEM *boundaries;                          /*!< \brief Array of the boundaries of the FEM mesh. */

  unsigned short nStandardBoundaryFacesSol; /*!< \brief Number of standard boundary faces used for solution of the DG solver. */
  unsigned short nStandardElementsSol;      /*!< \brief Number of standard volume elements used for solution of the DG solver. */
  unsigned short nStandardMatchingFacesSol; /*!< \brief Number of standard matching internal faces used for solution of the DG solver. */

  const CFEMStandardBoundaryFace *standardBoundaryFacesSol; /*!< \brief Array that contains the standard boundary
                                                                        faces used for the solution of the DG solver. */
  const CFEMStandardElement *standardElementsSol; /*!< \brief Array that contains the standard volume elements
                                                              used for the solution of the DG solver. */
  const CFEMStandardInternalFace *standardMatchingFacesSol;  /*!< \brief Array that contains the standard matching
                                                                         internal faces used for the solution of
                                                                         the DG solver. */

  const su2double *timeCoefADER_DG;                        /*!< \brief The time coefficients in the iteration matrix of
                                                                       the ADER-DG predictor step. */
  const su2double *timeInterpolDOFToIntegrationADER_DG;    /*!< \brief The interpolation matrix between the time DOFs and
                                                                       the time integration points for ADER-DG. */
  const su2double *timeInterpolAdjDOFToIntegrationADER_DG; /*!< \brief The interpolation matrix between the time DOFs of adjacent
                                                                       elements of a higher time level and the time integration
                                                                       points for ADER-DG. */

  unsigned int sizeWorkArray;     /*!< \brief The size of the work array needed. */

  vector<su2double> TolSolADER;   /*!< \brief Vector, which stores the tolerances for the conserved
                                              variables in the ADER predictor step. */

  vector<su2double> VecSolDOFs;    /*!< \brief Vector, which stores the solution variables in the owned DOFs. */
  vector<su2double> VecSolDOFsNew; /*!< \brief Vector, which stores the new solution variables in the owned DOFs (needed for classical RK4 scheme). */
  vector<su2double> VecDeltaTime;  /*!< \brief Vector, which stores the time steps of the owned volume elements. */

  vector<su2double> VecSolDOFsPredictorADER; /*!< \brief Vector, which stores the ADER predictor solution in the owned
                                                         DOFs. These are both space and time DOFs. */

  vector<vector<su2double> > VecWorkSolDOFs; /*!< \brief Working double vector to store the conserved variables for
                                                         the DOFs for the different time levels. */

  vector<su2double> VecResDOFs;         /*!< \brief Vector, which stores the residuals in the owned DOFs. */
  vector<su2double> VecResFaces;        /*!< \brief Vector, which stores the residuals of the DOFs that
                                                    come from the faces, both boundary and internal. */
  vector<su2double> VecTotResDOFsADER;  /*!< \brief Vector, which stores the accumulated residuals of the
                                                    owned DOFs for the ADER corrector step. */


  vector<unsigned long> nEntriesResFaces; /*!< \brief Number of entries for the DOFs in the
                                                      residual of the faces. Cumulative storage. */
  vector<unsigned long> entriesResFaces;  /*!< \brief The corresponding entries in the residual of the faces. */

  vector<unsigned long> nEntriesResAdjFaces; /*!< \brief Number of entries for the DOFs in the residual of the faces,
                                                         where the face is adjacent to an element of lower time
                                                         level. Cumulative storage. */
  vector<unsigned long> entriesResAdjFaces;  /*!< \brief The corresponding entries in the residual of the faces. */

  vector<vector<unsigned long> > startLocResFacesMarkers; /*!< \brief The starting location in the residual of the
                                                                      faces for the time levels of the boundary
                                                                      markers. */

  vector<unsigned long> startLocResInternalFacesLocalElem; /*!< \brief The starting location in the residual of the
                                                                       faces for the time levels of internal faces
                                                                       between locally owned elements. */
  vector<unsigned long> startLocResInternalFacesWithHaloElem; /*!< \brief The starting location in the residual of the
                                                                          faces for the time levels of internal faces
                                                                          between an owned and a halo element. */

  bool symmetrizingTermsPresent;    /*!< \brief Whether or not symmetrizing terms are present in the
                                                discretization. */

  vector<unsigned long> nDOFsPerRank;                    /*!< \brief Number of DOFs per rank in
                                                                     cumulative storage format. */
  vector<vector<unsigned long> > nonZeroEntriesJacobian; /*!< \brief The ID's of the DOFs for the
                                                                     non-zero entries of the Jacobian
                                                                     for the locally owned DOFs. */

  int nGlobalColors;              /*!< \brief Number of global colors for the Jacobian computation. */

  vector<vector<unsigned long> > localDOFsPerColor;   /*!< \brief Double vector, which contains for every
                                                                  color the local DOFs. */
  vector<vector<int> > colorToIndEntriesJacobian;     /*!< \brief Double vector, which contains for every
                                                                  local DOF the mapping from the color to the
                                                                  entry in the Jacobian. A -1 indicates that
                                                                  the color does not contribute to the Jacobian
                                                                  of the DOF. */

  CBlasStructure *blasFunctions; /*!< \brief  Pointer to the object to carry out the BLAS functionalities. */

private:

#ifdef HAVE_MPI
  vector<vector<SU2_MPI::Request> > commRequests;  /*!< \brief Communication requests in the communication of the solution for all
                                                               time levels. These are both sending and receiving requests. */

  vector<vector<vector<unsigned long> > > elementsRecvMPIComm;  /*!< \brief Triple vector, which contains the halo elements
                                                                            for MPI communication for all time levels. */
  vector<vector<vector<unsigned long> > > elementsSendMPIComm;  /*!< \brief Triple vector, which contains the donor elements
                                                                            for MPI communication for all time levels. */

  vector<vector<int> > ranksRecvMPI; /*!< \brief Double vector, which contains the ranks from which the halo elements
                                                 are received for all time levels. */
  vector<vector<int> > ranksSendMPI; /*!< \brief Double vector, which contains the ranks to which the donor elements
                                                 are sent for all time levels. */

  vector<vector<vector<su2double> > > commRecvBuf;  /*!< \brief Receive buffers used to receive the solution data
                                                                in the communication pattern for all time levels. */
  vector<vector<vector<su2double> > > commSendBuf;  /*!< \brief Send buffers used to send the solution data
                                                                in the communication pattern for all time levels. */
#endif

  vector<vector<unsigned long> > elementsRecvSelfComm;  /*!< \brief Double vector, which contains the halo elements
                                                                    for self communication for all time levels. */
  vector<vector<unsigned long> > elementsSendSelfComm;  /*!< \brief Double vector, which contains the donor elements
                                                                    for self communication for all time levels. */

  vector<su2double> rotationMatricesPeriodicity;    /*!< \brief Vector, which contains the rotation matrices
                                                                for the rotational periodic transformations. */
  vector<vector<vector<unsigned long> > > halosRotationalPeriodicity; /*!< \brief Triple vector, which contains the indices
                                                                                  of halo elements for which a periodic
                                                                                  transformation must be applied for all
                                                                                  time levels. */

  vector<CTaskDefinition> tasksList; /*!< \brief List of tasks to be carried out in the computationally
                                                 intensive part of the solver. */
public:

  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_EulerSolver(void);

  /*!
   * \overload
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nDim - Dimension of the problem (2D or 3D).
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CFEM_DG_EulerSolver(CConfig *config, unsigned short val_nDim, unsigned short iMesh);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CFEM_DG_EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFEM_DG_EulerSolver(void);

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] config      - Definition of the particular problem.
   * \param[in] iMesh       - Index of the mesh in multigrid computations.
   * \param[in] writeOutput - Whether or not output must be written.
   */
  void SetNondimensionalization(CConfig        *config,
                                unsigned short iMesh,
                                const bool     writeOutput);
  using CSolver::SetNondimensionalization;

  /*!
   * \brief Get a pointer to the vector of the solution degrees of freedom.
   * \return Pointer to the vector of the solution degrees of freedom.
   */
  su2double* GetVecSolDOFs(void);

  /*!
   * \brief Get the global number of solution degrees of freedom for the calculation.
   * \return Global number of solution degrees of freedom
   */
  unsigned long GetnDOFsGlobal(void);

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
   * \brief Set the freestream pressure.
   * \param[in] Value of freestream pressure.
   */
  void SetPressure_Inf(su2double p_inf);

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  void SetTemperature_Inf(su2double t_inf);

  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container,
                           CConfig *config, unsigned long ExtIter);

  /*!
   * \brief Set the working solution of the first time level to the current
            solution. Used for Runge-Kutta type schemes.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void Set_OldSolution(CGeometry *geometry);

  /*!
   * \brief Set the new solution to the current solution for classical RK.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void Set_NewSolution(CGeometry *geometry);

  /*!
   * \brief Function to compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                    unsigned short iMesh, unsigned long Iteration);

  /*!
   * \brief Function, which checks whether or not the time synchronization point is reached
            when explicit time stepping is used.
   * \param[in]     config          - Definition of the particular problem.
   * \param[in]     TimeSync        - The synchronization time.
   * \param[in,out] timeEvolved     - On input the time evolved before the time step,
                                      on output the time evolved after the time step.
   * \param[out]    syncTimeReached - Whether or not the synchronization time is reached.
   */
  void CheckTimeSynchronization(CConfig         *config,
                                const su2double TimeSync,
                                su2double       &timeEvolved,
                                bool            &syncTimeReached);

  /*!
   * \brief Function, which processes the list of tasks to be executed by
            the DG solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                          CNumerics **numerics, CConfig *config,
                          unsigned short iMesh);

  /*!
   * \brief Function, to carry out the space time integration for ADER
            with time accurate local time stepping.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                 CNumerics **numerics, CConfig *config,
                                 unsigned short iMesh, unsigned short RunTime_EqSystem);

  /*!
   * \brief Function, which controls the computation of the spatial Jacobian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                              CNumerics **numerics, CConfig *config,
                              unsigned short iMesh, unsigned short RunTime_EqSystem);

  /*!
   * \brief Function, which determines the values of the tolerances in
            the predictor step of ADER-DG.
   */
  void TolerancesADERPredictorStep(void);

  /*!
   * \brief Function, carries out the predictor step of the ADER-DG
            time integration.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void ADER_DG_PredictorStep(CConfig             *config,
                             const unsigned long elemBeg,
                             const unsigned long elemEnd,
                             su2double           *workArray);

  /*!
   * \brief Function, which interpolates the predictor solution of ADER-DG
            to the time value that corresponds to iTime.
   * \param[in]  config            - Definition of the particular problem.
   * \param[in]  iTime             - Time index of the time integration point for the
                                     integration over the time slab in the corrector
                                     step of ADER-DG.
   * \param[in]  elemBeg           - Begin index of the element range to be computed. This
                                     range is for elements of the same time level.
   * \param[in]  elemEnd           - End index (not included) of the element range to be computed.
   * \param[in]  nAdjElem          - Number of elements of the next time level, which are
                                     adjacent to elements of the current time level.
   * \param[in]  adjElem           - The ID's of the adjacent elements.
   * \param[in]  secondPartTimeInt - Whether or not this is the second part of the
                                     time interval for the adjacent elements.
   * \param[out] solTimeLevel      - Array in which the interpolated solution for the
                                     time level considered must be stored.
   */
  void ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                           const unsigned short iTime,
                                           const unsigned long  elemBeg,
                                           const unsigned long  elemEnd,
                                           const unsigned long  nAdjElem,
                                           const unsigned long  *adjElem,
                                           const bool           secondPartTimeInt,
                                           su2double            *solTimeLevel);

  /*!
   * \brief Compute the artificial viscosity for shock capturing in DG. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  virtual void Shock_Capturing_DG(CConfig             *config,
                                  const unsigned long elemBeg,
                                  const unsigned long elemEnd,
                                  su2double           *workArray);

  /*!
   * \brief Compute the volume contributions to the spatial residual. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  virtual void Volume_Residual(CConfig             *config,
                               const unsigned long elemBeg,
                               const unsigned long elemEnd,
                               su2double           *workArray);

  /*!
   * \brief Function, which computes the spatial residual for the DG discretization.
   * \param[in]  timeLevel           - Time level of the time accurate local time stepping,
                                       if relevant.
   * \param[in]  config              - Definition of the particular problem.
   * \param[in]  numerics            - Description of the numerical method.
   * \param[in]  haloInfoNeededForBC - If true,  treat boundaries for which halo data is needed.
                                       If false, treat boundaries for which only owned data is needed.
   * \param[out] workArray           - Work array.
   */
  void Boundary_Conditions(const unsigned short timeLevel,
                           CConfig              *config,
                           CNumerics            **numerics,
                           const bool           haloInfoNeededForBC,
                           su2double            *workArray);

  /*!
   * \brief Compute the spatial residual for the given range of faces. It is a virtual
            function, because this function is overruled for Navier-Stokes.
   * \param[in]     config      - Definition of the particular problem.
   * \param[in]     indFaceBeg  - Starting index in the matching faces.
   * \param[in]     indFaceEnd  - End index in the matching faces.
   * \param[in,out] indResFaces - Index where to store the residuals in
                                  the vector of face residuals.
   * \param[in]     numerics    - Description of the numerical method.
   * \param[out]    workArray   - Work array.
   */
  virtual void ResidualFaces(CConfig             *config,
                             const unsigned long indFaceBeg,
                             const unsigned long indFaceEnd,
                             unsigned long       &indResFaces,
                             CNumerics           *numerics,
                             su2double           *workArray);

  /*!
   * \brief Function, which accumulates the space time residual of the ADER-DG
            time integration scheme for the owned elements.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - time level for which the residuals must be
                            accumulated.
   * \param[in] intPoint  - Index of the time integration point.
   */
  void AccumulateSpaceTimeResidualADEROwnedElem(CConfig             *config,
                                                const unsigned short timeLevel,
                                                const unsigned short intPoint);

  /*!
   * \brief Function, which accumulates the space time residual of the ADER-DG
            time integration scheme for the halo elements.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - time level for which the residuals must be
                            accumulated.
   * \param[in] intPoint  - Index of the time integration point.
   */
  void AccumulateSpaceTimeResidualADERHaloElem(CConfig             *config,
                                               const unsigned short timeLevel,
                                               const unsigned short intPoint);

  /*!
   * \brief Compute primitive variables and their gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iStep - Current step in the time accurate local time
                        stepping algorithm, if appropriate.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                     unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Impose via the residual the Euler wall boundary condition. It is a
            virtual function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Euler_Wall(CConfig                  *config,
                             const unsigned long      surfElemBeg,
                             const unsigned long      surfElemEnd,
                             const CSurfaceElementFEM *surfElem,
                             su2double                *resFaces,
                             CNumerics                *conv_numerics,
                             su2double                *workArray);
  using CSolver::BC_Euler_Wall;

  /*!
   * \brief Impose the far-field boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Far_Field(CConfig                  *config,
                            const unsigned long      surfElemBeg,
                            const unsigned long      surfElemEnd,
                            const CSurfaceElementFEM *surfElem,
                            su2double                *resFaces,
                            CNumerics                *conv_numerics,
                            su2double                *workArray);
  using CSolver::BC_Far_Field;

  /*!
   * \brief Impose the symmetry boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Sym_Plane(CConfig                  *config,
                            const unsigned long      surfElemBeg,
                            const unsigned long      surfElemEnd,
                            const CSurfaceElementFEM *surfElem,
                            su2double                *resFaces,
                            CNumerics                *conv_numerics,
                            su2double                *workArray);
  using CSolver::BC_Sym_Plane;

  /*!
   * \brief Impose the supersonic outlet boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Supersonic_Outlet(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray);
  using CSolver::BC_Supersonic_Outlet;

  /*!
   * \brief Impose the subsonic inlet boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Inlet(CConfig                  *config,
                        const unsigned long      surfElemBeg,
                        const unsigned long      surfElemEnd,
                        const CSurfaceElementFEM *surfElem,
                        su2double                *resFaces,
                        CNumerics                *conv_numerics,
                        unsigned short           val_marker,
                        su2double                *workArray);
  using CSolver::BC_Inlet;

  /*!
   * \brief Impose the outlet boundary condition.It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Outlet(CConfig                  *config,
                         const unsigned long      surfElemBeg,
                         const unsigned long      surfElemEnd,
                         const CSurfaceElementFEM *surfElem,
                         su2double                *resFaces,
                         CNumerics                *conv_numerics,
                         unsigned short           val_marker,
                         su2double                *workArray);
  using CSolver::BC_Outlet;

  /*!
   * \brief Impose a constant heat-flux condition at the wall. It is a virtual
            function, such that it can be overwritten for Navier-Stokes.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_HeatFlux_Wall(CConfig                  *config,
                                const unsigned long      surfElemBeg,
                                const unsigned long      surfElemEnd,
                                const CSurfaceElementFEM *surfElem,
                                su2double                *resFaces,
                                CNumerics                *conv_numerics,
                                unsigned short           val_marker,
                                su2double                *workArray);
  using CSolver::BC_HeatFlux_Wall;

  /*!
   * \brief Impose an isothermal condition at the wall. It is a virtual
            function, such that it can be overwritten for Navier-Stokes.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Isothermal_Wall(CConfig                  *config,
                                  const unsigned long      surfElemBeg,
                                  const unsigned long      surfElemEnd,
                                  const CSurfaceElementFEM *surfElem,
                                  su2double                *resFaces,
                                  CNumerics                *conv_numerics,
                                  unsigned short           val_marker,
                                  su2double                *workArray);
  using CSolver::BC_Isothermal_Wall;

  /*!
   * \brief Impose the boundary condition using characteristic reconstruction. It is
   *        a virtual function, such that it can be overwritten for Navier-Stokes.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Riemann(CConfig                  *config,
                          const unsigned long      surfElemBeg,
                          const unsigned long      surfElemEnd,
                          const CSurfaceElementFEM *surfElem,
                          su2double                *resFaces,
                          CNumerics                *conv_numerics,
                          unsigned short           val_marker,
                          su2double                *workArray);
  using CSolver::BC_Riemann;

  /*!
   * \brief Impose the user customized boundary condition. It is a virtual
            function, because for Navier-Stokes it is overwritten.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  virtual void BC_Custom(CConfig                  *config,
                         const unsigned long      surfElemBeg,
                         const unsigned long      surfElemEnd,
                         const CSurfaceElementFEM *surfElem,
                         su2double                *resFaces,
                         CNumerics                *conv_numerics,
                         su2double                *workArray);
  using CSolver::BC_Custom;

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
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                              unsigned short iRKStep);

  /*!
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetResidual_RMS_FEM(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the global error measures (L2, Linf) for verification cases.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVerificationError(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the solution for the ADER-DG scheme for the given range
            of elements.
   * \param[in] elemBeg - Begin index of the element range to be computed.
   * \param[in] elemEnd - End index (not included) of the element range to be computed.
   */
  void ADER_DG_Iteration(const unsigned long elemBeg,
                         const unsigned long elemEnd);

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Pressure_Forces(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

  /*!
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCL_Inv(unsigned short val_marker);

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
  su2double GetSurface_CL(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF(unsigned short val_marker);

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
  su2double GetSurface_CL_Inv(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CD_Inv(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  su2double GetSurface_CSF_Inv(unsigned short val_marker);

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
  su2double GetCD_Inv(unsigned short val_marker);

  /*!
   * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCSF_Inv(unsigned short val_marker);

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
  su2double GetTotal_CSF(void);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CEff(void);

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  void SetTotal_CL(su2double val_Total_CL);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CL(void);

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  su2double GetTotal_CD(void);

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
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  void SetTotal_CD(su2double val_Total_CD);

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  su2double GetAllBound_CL_Inv(void);

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  su2double GetAllBound_CD_Inv(void);

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  su2double GetAllBound_CSF_Inv(void);

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

protected:

  /*!
   * \brief Routine that initiates the non-blocking communication between ranks
            for the givem time level.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - The time level for which the communication must be
                            initiated.
   */
  void Initiate_MPI_Communication(CConfig *config,
                                  const unsigned short timeLevel);

  /*!
   * \brief Routine that initiates the reverse non-blocking communication
            between ranks.
   * \param[in] config    - Definition of the particular problem.
   * \param[in] timeLevel - The time level for which the reverse communication
                            must be initiated.
   */
  void Initiate_MPI_ReverseCommunication(CConfig *config,
                                         const unsigned short timeLevel);

  /*!
   * \brief Routine that completes the non-blocking communication between ranks.
   * \param[in] config              - Definition of the particular problem.
   * \param[in] timeLevel           - The time level for which the communication
                                      may be completed.
   * \param[in] commMustBeCompleted - Whether or not the communication must be completed.
   * \return  Whether or not the communication has been completed.
   */
  bool Complete_MPI_Communication(CConfig *config,
                                  const unsigned short timeLevel,
                                  const bool commMustBeCompleted);

  /*!
   * \brief Routine that completes the reverse non-blocking communication
            between ranks.
   * \param[in] config              - Definition of the particular problem.
   * \param[in] timeLevel           - The time level for which the communication
                                      may be completed.
   * \param[in] commMustBeCompleted - Whether or not the communication must be completed.
   * \return  Whether or not the communication has been completed.
   */
  bool Complete_MPI_ReverseCommunication(CConfig *config,
                                         const unsigned short timeLevel,
                                         const bool commMustBeCompleted);

  /*!
   * \brief Function, which computes the inviscid fluxes in face points.
   * \param[in]  config       - Definition of the particular problem.
   * \param[in]  nFaceSimul   - Number of faces that are treated simultaneously
                                to improve performance.
   * \param[in]  NPad         - Value of the padding parameter to obtain optimal
                                performance in the gemm computations.
   * \param[in]  nPoints      - Number of points per face for which the fluxes
                                must be computed.
   * \param[in]  normalsFace  - The normals in the points for the faces.
   * \param[in]  gridVelsFace - The grid velocities in the points for the faces.
   * \param[in]  solL         - Solution in the left state of the points.
   * \param[in]  solR         - Solution in the right state of the points.
   * \param[out] fluxes       - Inviscid fluxes in the points.
   * \param[in]  numerics     - Object, which contains the Riemann solver.
   */
  void ComputeInviscidFluxesFace(CConfig              *config,
                                 const unsigned short nFaceSimul,
                                 const unsigned short NPad,
                                 const unsigned long  nPoints,
                                 const su2double      *normalsFace[],
                                 const su2double      *gridVelsFace[],
                                 const su2double      *solL,
                                 const su2double      *solR,
                                 su2double            *fluxes,
                                 CNumerics            *numerics);

  /*!
   * \brief Function, which computes the inviscid fluxes in the face integration
            points of a chunk of matching internal faces.
   * \param[in]  config   - Definition of the particular problem.
   * \param[in]  lBeg     - Start index in matchingInternalFaces for which
                            the inviscid fluxes should be computed.
   * \param[in]  lEnd     - End index (not included) in matchingInternalFaces
                            for which the inviscid fluxes should be computed.
   * \param[in]  NPad     - Value of the padding parameter to obtain optimal
                            performance in the gemm computations.
   * \param[out] solIntL  - Solution in the left state of the integration points.
   * \param[out] solIntR  - Solution in the right state of the integration points.
   * \param[out] fluxes   - Inviscid fluxes in the integration points.
   * \param[in]  numerics - Object, which contains the Riemann solver.
   */
  void InviscidFluxesInternalMatchingFace(CConfig              *config,
                                          const unsigned long  lBeg,
                                          const unsigned long  lEnd,
                                          const unsigned short NPad,
                                          su2double            *solIntL,
                                          su2double            *solIntR,
                                          su2double            *fluxes,
                                          CNumerics            *numerics);
  /*!
   * \brief Function, which computes the left state of a boundary face.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of faces that are treated simultaneously
                              to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary elements for which the left state must be computed.
   * \param[out] solFace    - Temporary storage for the solution in the DOFs.
   * \param[out] solIntL    - Left states in the integration points of the face.
   */
  void LeftStatesIntegrationPointsBoundaryFace(CConfig                  *config,
                                               const unsigned short     nFaceSimul,
                                               const unsigned short     NPad,
                                               const CSurfaceElementFEM *surfElem,
                                               su2double                *solFace,
                                               su2double                *solIntL);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the inviscid wall boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of faces that are treated simultaneously
                              to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary elements for which the left state must
                              be computed.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Euler_Wall(CConfig                  *config,
                                 const unsigned short     nFaceSimul,
                                 const unsigned short     NPad,
                                 const CSurfaceElementFEM *surfElem,
                                 const su2double          *solIntL,
                                 su2double                *solIntR);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the inlet boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of fused faces, i.e. the number of faces
                              that are treated simultaneously to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary element for which the left state must be computed.
   * \param[in]  val_marker - Surface marker where the boundary condition is applied.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Inlet(CConfig                  *config,
                            const unsigned short     nFaceSimul,
                            const unsigned short     NPad,
                            const CSurfaceElementFEM *surfElem,
                            unsigned short           val_marker,
                            const su2double          *solIntL,
                            su2double                *solIntR);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the outlet boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of fused faces, i.e. the number of faces
                              that are treated simultaneously to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary element for which the left state must be computed.
   * \param[in]  val_marker - Surface marker where the boundary condition is applied.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Outlet(CConfig                  *config,
                             const unsigned short     nFaceSimul,
                             const unsigned short     NPad,
                             const CSurfaceElementFEM *surfElem,
                             unsigned short           val_marker,
                             const su2double          *solIntL,
                             su2double                *solIntR);

  /*!
   * \brief Function, which computes the boundary states in the integration points
            of the boundary face by applying the Riemann boundary conditions.
   * \param[in]  config     - Definition of the particular problem.
   * \param[in]  nFaceSimul - Number of fused faces, i.e. the number of faces
                              that are treated simultaneously to improve performance.
   * \param[in]  NPad       - Value of the padding parameter to obtain optimal
                              performance in the gemm computations.
   * \param[in]  surfElem   - Surface boundary element for which the left state must be computed.
   * \param[in]  val_marker - Surface marker where the boundary condition is applied.
   * \param[in]  solIntL    - Left states in the integration points of the face.
   * \param[out] solIntR    - Right states in the integration points of the face.
   */
  void BoundaryStates_Riemann(CConfig                  *config,
                              const unsigned short     nFaceSimul,
                              const unsigned short     NPad,
                              const CSurfaceElementFEM *surfElem,
                              unsigned short           val_marker,
                              const su2double          *solIntL,
                              su2double                *solIntR);
private:

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                   CVolumeElementFEM    *elem,
                                                   const su2double      *sol,
                                                   const unsigned short nSimul,
                                                   const unsigned short NPad,
                                                   su2double            *res,
                                                   su2double            *work);

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                   CVolumeElementFEM    *elem,
                                                   const su2double      *sol,
                                                   const unsigned short nSimul,
                                                   const unsigned short NPad,
                                                   su2double            *res,
                                                   su2double            *work);

  /*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                      CVolumeElementFEM    *elem,
                                                      const su2double      *sol,
                                                      const unsigned short nSimul,
                                                      const unsigned short NPad,
                                                      su2double            *res,
                                                      su2double            *work);

/*!
   * \brief Virtual function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  virtual void ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                      CVolumeElementFEM    *elem,
                                                      const su2double      *sol,
                                                      const unsigned short nSimul,
                                                      const unsigned short NPad,
                                                      su2double            *res,
                                                      su2double            *work);

  /*!
   * \brief Function, which computes the graph of the spatial discretization
            for the locally owned DOFs.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void DetermineGraphDOFs(const CMeshFEM *FEMGeometry,
                          CConfig        *config);

  /*!
   * \brief Function, which determines the meta data needed for the computation
            of the Jacobian of the spatial residual.
   * \param[in] DGGeometry     - Geometrical definition of the DG problem.
   * \param[in] colorLocalDOFs - Color of the locally stored DOFs.
   */
  void MetaDataJacobianComputation(const CMeshFEM    *FEMGeometry,
                                   const vector<int> &colorLocalDOFs);

  /*!
   * \brief Function, which sets up the list of tasks to be carried out in the
            computationally expensive part of the solver.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpTaskList(CConfig *config);

  /*!
   * \brief Function, which sets up the persistent communication of the flow
            variables in the DOFs.
   * \param[in] DGGeometry - Geometrical definition of the DG problem.
   * \param[in] config     - Definition of the particular problem.
   */
  void Prepare_MPI_Communication(const CMeshFEM *FEMGeometry,
                                 CConfig *config);

  /*!
   * \brief Function, which creates the final residual by summing up
            the contributions for the DOFs of the elements considered.
   * \param[in] timeLevel     - Time level of the elements for which the
                                final residual must be created.
   * \param[in] ownedElements - Whether owned or halo elements must be treated.
   */
  void CreateFinalResidual(const unsigned short timeLevel,
                           const bool ownedElements);

  /*!
   * \brief Function, which multiplies the residual by the inverse
            of the (lumped) mass matrix.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  useADER   - Whether or not the ADER residual must be multiplied.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void MultiplyResidualByInverseMassMatrix(CConfig             *config,
                                           const bool          useADER,
                                           const unsigned long elemBeg,
                                           const unsigned long elemEnd,
                                           su2double           *workArray);

  /*!
   * \brief Function, which computes the residual contribution from a boundary
            face in an inviscid computation when the boundary conditions have
            already been applied.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     nFaceSimul    - Number of faces that are treated simultaneously
                                    to improve performance.
   * \param[in]     NPad          - Value of the padding parameter to obtain optimal
                                    performance in the gemm computations.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     surfElem      - Surface boundary element for which the
                                    contribution to the residual must be computed.
   * \param[in]     solInt0       - Solution in the integration points of side 0.
                                    It is not const, because the array is used for
                                    temporary storage for the residual.
   * \param[in]     solInt1       - Solution in the integration points of side 1.
   * \param[out]    fluxes        - Temporary storage for the fluxes in the
                                    integration points.
   * \param[out]    resFaces      - Array to store the residuals of the face.
   * \param[in,out] indResFaces   - Index in resFaces, where the current residual
                                    should be stored. It is updated in the function
                                    for the next boundary element.
   */
  void ResidualInviscidBoundaryFace(CConfig                  *config,
                                    const unsigned short     nFaceSimul,
                                    const unsigned short     NPad,
                                    CNumerics                *conv_numerics,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *solInt0,
                                    su2double                *solInt1,
                                    su2double                *fluxes,
                                    su2double                *resFaces,
                                    unsigned long            &indResFaces);

protected:
  /*!
   * \brief Template function, which determines some meta data for the chunk of
            elements/faces that must be treated simulaneously.
   * \param[in]  elem       - Const pointer the volume or face elements for which
                              the meta data must be computed.
   * \param[in]  l          - Start index for the current chunk of elements/faces.
   * \param[in]  elemEnd    - End index (index not included) of the elements to be
                              treated in the residual computation from which this
                              function is called.
   * \param[in]  nElemSimul - Desired number of elements/faces that must be treated
                              simultaneously for optimal performance.
   * \param[in]  nPadMin    - Minimum number of the padding value in the gemm calls.
   * \param[out] lEnd       - Actual end index (not included) for this chunk of
                              elements.
   * \param[out] ind        - Index in the standard elements to which this chunk of
                              elements can be mapped.
   * \param[out] llEnd      - Actual number of elements/faces that are treated
                              simultaneously, llEnd = lEnd - l.
   * \param[out] NPad       - Actual padded N value in the gemm computations for
                              this chunk of elements.
   */
  template <class TElemType>
  void MetaDataChunkOfElem(const TElemType      *elem,
                           const unsigned long  l,
                           const unsigned long  elemEnd,
                           const unsigned short nElemSimul,
                           const unsigned short nPadMin,
                           unsigned long        &lEnd,
                           unsigned short       &ind,
                           unsigned short       &llEnd,
                           unsigned short       &NPad) {

    /* Determine the end index for this chunk of elements that must be
       treated simulaneously. The elements of this chunk must have the
       same standard element in order to make this work. */
    const unsigned long lEndMax = min(l+nElemSimul, elemEnd);

    ind = elem[l].indStandardElement;
    for(lEnd=l+1; lEnd<lEndMax; ++lEnd) {
      if(elem[lEnd].indStandardElement != ind) break;
    }

    /* Store the number of elements that are treated simultaneously in this chunk
       in llEnd and determine the padded N value in the gemm computations. */
    llEnd = lEnd - l;
    NPad  = llEnd*nVar;
    if( NPad%nPadMin ) NPad += nPadMin - (NPad%nPadMin);
  }
};

/*!
 * \class CFEM_DG_NSSolver
 * \brief Main class for defining the Navier-Stokes Discontinuous Galerkin finite element flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 6.2.0 "Falcon"
 */
class CFEM_DG_NSSolver : public CFEM_DG_EulerSolver {
private:
  su2double Viscosity_Inf; /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;       /*!< \brief Turbulent kinetic energy at the infinity. */
  su2double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb;            /*!< \brief Turbulent Prandtl number. */

  CSGSModel *SGSModel;     /*!< \brief LES Subgrid Scale model. */
  bool SGSModelUsed;       /*!< \brief Whether or not an LES Subgrid Scale model is used. */

  su2double
  *CL_Visc, 	            /*!< \brief Lift coefficient (viscous contribution) for each boundary. */
  *CD_Visc,              /*!< \brief Drag coefficient (viscous contribution) for each boundary. */
  *CSF_Visc,         /*!< \brief Side force coefficient (viscous contribution) for each boundary. */
  *CMx_Visc, 	            /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CMy_Visc, 	            /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CMz_Visc,          	    /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CFx_Visc,          	    /*!< \brief Force x coefficient (viscous contribution) for each boundary. */
  *CFy_Visc,          	    /*!< \brief Force y coefficient (viscous contribution) for each boundary. */
  *CFz_Visc, 	            /*!< \brief Force z coefficient (viscous contribution) for each boundary. */
  *CEff_Visc,               /*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
  *Surface_CL_Visc,      /*!< \brief Lift coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CD_Visc,      /*!< \brief Drag coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CSF_Visc, /*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CEff_Visc,       /*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFx_Visc,        /*!< \brief Force x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFy_Visc,        /*!< \brief Force y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFz_Visc,        /*!< \brief Force z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMx_Visc,        /*!< \brief Moment x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMy_Visc,        /*!< \brief Moment y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMz_Visc,        /*!< \brief Moment z coefficient (viscous contribution) for each monitoring surface. */
  *Heat_Visc,               /*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHeatFlux_Visc;        /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */

  su2double
  AllBound_CD_Visc,      /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
  AllBound_CL_Visc, 	    /*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
  AllBound_CSF_Visc, /*!< \brief Sideforce coefficient (viscous contribution) for all the boundaries. */
  AllBound_CMx_Visc, 	    /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Visc, 	    /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Visc, 	    /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Visc, 	    /*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CFx_Visc, 	    /*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Visc, 	    /*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Visc, 	    /*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_HeatFlux_Visc,   /*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHeatFlux_Visc; /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
  su2double StrainMag_Max, Omega_Max; /*!< \brief Maximum Strain Rate magnitude and Omega. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_NSSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_DG_NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DG_NSSolver(void);

  /*!
   * \brief Function to compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                    unsigned short iMesh, unsigned long Iteration);

  /*!
   * \brief Compute the artificial viscosity for shock capturing in DG.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void Shock_Capturing_DG(CConfig             *config,
                          const unsigned long elemBeg,
                          const unsigned long elemEnd,
                          su2double           *workArray);

  /*!
   * \brief Per-Olof Persson's method for capturing shock in DG
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void Shock_Capturing_DG_Persson(const unsigned long elemBeg,
                                  const unsigned long elemEnd,
                                  su2double           *workArray);

  /*!
   * \brief Compute the volume contributions to the spatial residual.
   * \param[in]  config    - Definition of the particular problem.
   * \param[in]  elemBeg   - Begin index of the element range to be computed.
   * \param[in]  elemEnd   - End index (not included) of the element range to be computed.
   * \param[out] workArray - Work array.
   */
  void Volume_Residual(CConfig             *config,
                       const unsigned long elemBeg,
                       const unsigned long elemEnd,
                       su2double           *workArray);

  /*!
   * \brief Compute the spatial residual for the given range of faces.
   * \param[in]     config      - Definition of the particular problem.
   * \param[in]     indFaceBeg  - Starting index in the matching faces.
   * \param[in]     indFaceEnd  - End index in the matching faces.
   * \param[in,out] indResFaces - Index where to store the residuals in
                                  the vector of face residuals.
   * \param[in]     numerics    - Description of the numerical method.
   * \param[out]    workArray   - Work array.
   */
  void ResidualFaces(CConfig             *config,
                     const unsigned long indFaceBeg,
                     const unsigned long indFaceEnd,
                     unsigned long       &indResFaces,
                     CNumerics           *numerics,
                     su2double           *workArray);

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray - Work array.
   */
  void BC_Euler_Wall(CConfig                  *config,
                     const unsigned long      surfElemBeg,
                     const unsigned long      surfElemEnd,
                     const CSurfaceElementFEM *surfElem,
                     su2double                *resFaces,
                     CNumerics                *conv_numerics,
                     su2double                *workArray);

  /*!
   * \brief Impose the far-field boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Far_Field(CConfig                  *config,
                    const unsigned long      surfElemBeg,
                    const unsigned long      surfElemEnd,
                    const CSurfaceElementFEM *surfElem,
                    su2double                *resFaces,
                    CNumerics                *conv_numerics,
                    su2double                *workArray);

  /*!
   * \brief Impose the symmetry boundary condition using the residual.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Sym_Plane(CConfig                  *config,
                    const unsigned long      surfElemBeg,
                    const unsigned long      surfElemEnd,
                    const CSurfaceElementFEM *surfElem,
                    su2double                *resFaces,
                    CNumerics                *conv_numerics,
                    su2double                *workArray);

 /*!
   * \brief Impose the supersonic outlet boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Supersonic_Outlet(CConfig                  *config,
                            const unsigned long      surfElemBeg,
                            const unsigned long      surfElemEnd,
                            const CSurfaceElementFEM *surfElem,
                            su2double                *resFaces,
                            CNumerics                *conv_numerics,
                            su2double                *workArray);

  /*!
   * \brief Impose the subsonic inlet boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Inlet(CConfig                  *config,
                const unsigned long      surfElemBeg,
                const unsigned long      surfElemEnd,
                const CSurfaceElementFEM *surfElem,
                su2double                *resFaces,
                CNumerics                *conv_numerics,
                unsigned short           val_marker,
                su2double                *workArray);

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Outlet(CConfig                  *config,
                 const unsigned long      surfElemBeg,
                 const unsigned long      surfElemEnd,
                 const CSurfaceElementFEM *surfElem,
                 su2double                *resFaces,
                 CNumerics                *conv_numerics,
                 unsigned short           val_marker,
                 su2double                *workArray);

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_HeatFlux_Wall(CConfig                  *config,
                        const unsigned long      surfElemBeg,
                        const unsigned long      surfElemEnd,
                        const CSurfaceElementFEM *surfElem,
                        su2double                *resFaces,
                        CNumerics                *conv_numerics,
                        unsigned short           val_marker,
                        su2double                *workArray);

  /*!
   * \brief Impose an isothermal condition at the wall.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Isothermal_Wall(CConfig                  *config,
                          const unsigned long      surfElemBeg,
                          const unsigned long      surfElemEnd,
                          const CSurfaceElementFEM *surfElem,
                          su2double                *resFaces,
                          CNumerics                *conv_numerics,
                          unsigned short           val_marker,
                          su2double                *workArray);

  /*!
   * \brief Impose the boundary condition using characteristic reconstruction.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[in]  val_marker    - Surface marker where the boundary condition is applied.
   * \param[out] workArray     - Work array.
   */
  void BC_Riemann(CConfig                  *config,
                  const unsigned long      surfElemBeg,
                  const unsigned long      surfElemEnd,
                  const CSurfaceElementFEM *surfElem,
                  su2double                *resFaces,
                  CNumerics                *conv_numerics,
                  unsigned short           val_marker,
                  su2double                *workArray);

  /*!
   * \brief Impose the user customized boundary condition.
   * \param[in]  config        - Definition of the particular problem.
   * \param[in]  surfElemBeg   - Start index in the list of surface elements.
   * \param[in]  surfElemEnd   - End index (not included) in the list of surface elements.
   * \param[in]  surfElem      - Array of surface elements for which the boundary
                                 conditions must be imposed.
   * \param[out] resFaces      - Array where the residual contribution from the
                                 surface elements must be stored.
   * \param[in]  conv_numerics - Description of the numerical method.
   * \param[out] workArray     - Work array.
   */
  void BC_Custom(CConfig                  *config,
                 const unsigned long      surfElemBeg,
                 const unsigned long      surfElemEnd,
                 const CSurfaceElementFEM *surfElem,
                 su2double                *resFaces,
                 CNumerics                *conv_numerics,
                 su2double                *workArray);

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
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCL_Visc(unsigned short val_marker);

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
  su2double GetCSF_Visc(unsigned short val_marker);

  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  su2double GetCD_Visc(unsigned short val_marker);

  /*!
   * \brief Get the total non dimensional lift coefficient (viscous contribution).
   * \return Value of the lift coefficient (viscous contribution).
   */
  su2double GetAllBound_CL_Visc(void);

  /*!
   * \brief Get the total non dimensional sideforce coefficient (viscous contribution).
   * \return Value of the lift coefficient (viscous contribution).
   */
  su2double GetAllBound_CSF_Visc(void);

  /*!
   * \brief Get the total non dimensional drag coefficient (viscous contribution).
   * \return Value of the drag coefficient (viscous contribution).
   */
  su2double GetAllBound_CD_Visc(void);

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

private:

  /*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                           CVolumeElementFEM    *elem,
                                           const su2double      *sol,
                                           const unsigned short nSimul,
                                           const unsigned short NPad,
                                           su2double            *res,
                                           su2double            *work);

/*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using an
            aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                           CVolumeElementFEM    *elem,
                                           const su2double      *sol,
                                           const unsigned short nSimul,
                                           const unsigned short NPad,
                                           su2double            *res,
                                           su2double            *work);
  /*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 2D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                              CVolumeElementFEM    *elem,
                                              const su2double      *sol,
                                              const unsigned short nSimul,
                                              const unsigned short NPad,
                                              su2double            *res,
                                              su2double            *work);

  /*!
   * \brief Function, which computes the spatial residual of the ADER-DG
            predictor step for the given volume element and solution using a
            non-aliased discretization in 3D.
   * \param[in]  config  - Definition of the particular problem.
   * \param[in]  elem    - Volume element for which the spatial residual of the
                           predictor step must be computed.
   * \param[in]  sol     - Solution for which the residual must be computed.
   * \param[in]  nSimul  - Number of entities (typically time integration points)
                           that are treated simultaneously.
   * \param[in]  NPad    - Padded N value in the matrix multiplications to
                           obtain better performance. The solution sol is stored
                           with this padded value to avoid a memcpy.
   * \param[out] res     - Residual of the spatial DOFs to be computed by this
                           function.
   * \param[out] work    - Work array.
   */
  void ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                              CVolumeElementFEM    *elem,
                                              const su2double      *sol,
                                              const unsigned short nSimul,
                                              const unsigned short NPad,
                                              su2double            *res,
                                              su2double            *work);
  /*!
   * \brief Function to compute the penalty terms in the integration
            points of a face.
   * \param[in]  indFaceChunk        - Index of the face in the chunk of fused faces.
   * \param[in]  nInt                - Number of integration points of the face.
   * \param[in]  NPad                - Value of the padding parameter to obtain optimal
                                       performance in the gemm computations.
   * \param[in]  solInt0             - Solution in the integration points of side 0.
   * \param[in]  solInt1             - Solution in the integration points of side 1.
   * \param[in]  viscosityInt0       - Viscosity in the integration points of side 0.
   * \param[in]  viscosityInt1       - Viscosity in the integration points of side 1.
   * \param[in]  kOverCvInt0         - Heat conductivity divided by Cv in the
                                       integration points of side 0.
   * \param[in]  kOverCvInt1         - Heat conductivity divided by Cv in the
                                       integration points of side 1.
   * \param[in]  ConstPenFace        - Penalty constant for this face.
   * \param[in]  lenScale0           - Length scale of the element of side 0.
   * \param[in]  lenScale1           - Length scale of the element of side 1.
   * \param[in]  metricNormalsFace   - Metric terms in the integration points, which
                                       contain the normals.
   * \param[out] penaltyFluxes       - Penalty fluxes in the integration points.
   */
  void PenaltyTermsFluxFace(const unsigned short indFaceChunk,
                            const unsigned short nInt,
                            const unsigned short NPad,
                            const su2double      *solInt0,
                            const su2double      *solInt1,
                            const su2double      *viscosityInt0,
                            const su2double      *viscosityInt1,
                            const su2double      *kOverCvInt0,
                            const su2double      *kOverCvInt1,
                            const su2double      ConstPenFace,
                            const su2double      lenScale0,
                            const su2double      lenScale1,
                            const su2double      *metricNormalsFace,
                                  su2double      *penaltyFluxes);

  /*!
   * \brief Function, which performs the treatment of the boundary faces for
            the Navier-Stokes equations for the most of the boundary conditions.
   * \param[in]     config                 - Definition of the particular problem.
   * \param[in]     conv_numerics          - Description of the numerical method.
   * \param[in]     nFaceSimul             - Number of faces that are treated simultaneously
                                             to improve performance.
   * \param[in]     NPad                   - Value of the padding parameter to obtain optimal
                                             performance in the gemm computations.
   * \param[in]     Wall_HeatFlux          - The value of the prescribed heat flux.
   * \param[in]     HeatFlux_Prescribed    - Whether or not the heat flux is prescribed by
                                             e.g. the boundary conditions.
   * \param[in]     Wall_Temperature       - The value of the prescribed wall temperature.
   * \param[in]     Temperature_Prescribed - Whether or not the temperature is precribed
                                             by e.g. the boundary conditions.
   * \param[in]     surfElem               - Surface boundary elements for which the
                                             residuals mut be computed.
   * \param[in]     solIntL                - Left states in the integration points of the face.
   * \param[in]     solIntR                - Right states in the integration points of the face.
   * \param[out]    workArray              - Storage for the local arrays.
   * \param[out]    resFaces               - Array to store the residuals of the face.
   * \param[in,out] indResFaces            - Index in resFaces, where the current residual
                                             should be stored. It is updated in the function
                                             for the next boundary element.
   * \param[in,out] wallModel              - Possible pointer to the wall model treatment.
                                             NULL pointer indicates no wall model treatment.
   */
  void ViscousBoundaryFacesBCTreatment(CConfig                  *config,
                                       CNumerics                *conv_numerics,
                                       const unsigned short     nFaceSimul,
                                       const unsigned short     NPad,
                                       const su2double          Wall_HeatFlux,
                                       const bool               HeatFlux_Prescribed,
                                       const su2double          Wall_Temperature,
                                       const bool               Temperature_Prescribed,
                                       const CSurfaceElementFEM *surfElem,
                                       const su2double          *solIntL,
                                       const su2double          *solIntR,
                                             su2double          *workArray,
                                             su2double          *resFaces,
                                             unsigned long      &indResFaces,
                                             CWallModel         *wallModel);

  /*!
   * \brief Function, which computes the viscous fluxes in the integration
            points for the boundary faces that must be treated simulaneously.
            This function uses the standard approach for computing the fluxes,
            i.e. no wall modeling.
   * \param[in]  config               - Definition of the particular problem.
   * \param[in]  nFaceSimul           - Number of faces that are treated simultaneously
                                        to improve performance.
   * \param[in]  NPad                 - Value of the padding parameter to obtain optimal
                                        performance in the gemm computations.
   * \param[in]  nInt                 - Number of integration points on the face.
   * \param[in]  nDOFsElem            - Number of DOFs of the adjacent element.
   * \param[in]  Wall_HeatFlux        - The value of the prescribed heat flux.
   * \param[in]  HeatFlux_Prescribed  - Whether or not the heat flux is prescribed by
                                        e.g. the boundary conditions.
   * \param[in]  derBasisElem         - Array, which contains the derivatives of the
                                        basis functions of the adjacent element
                                        in the integration points.
   * \param[in]  surfElem             - Surface boundary elements for which the
                                        viscous fluxes must be computed.
   * \param[in]  solIntL              - Left states in the integration points of the face.
   * \param[out] solElem              - Storage for the solution in the adjacent elements.
   * \param[out] gradSolInt           - Storage for the gradients of the solution in the
                                        integration points of the face.
   * \param[out] viscFluxes           - To be computed viscous fluxes in the
                                        integration points.
   * \param[out] viscosityInt         - To be computed viscosity in the integration points.
   * \param[out] kOverCvInt           - To be computed thermal conductivity in the
                                        integration points.
   */
  void ComputeViscousFluxesBoundaryFaces(CConfig                  *config,
                                         const unsigned short     nFaceSimul,
                                         const unsigned short     NPad,
                                         const unsigned short     nInt,
                                         const unsigned short     nDOFsElem,
                                         const su2double          Wall_HeatFlux,
                                         const bool               HeatFlux_Prescribed,
                                         const su2double          *derBasisElem,
                                         const CSurfaceElementFEM *surfElem,
                                         const su2double          *solIntL,
                                               su2double          *solElem,
                                               su2double          *gradSolInt,
                                               su2double          *viscFluxes,
                                               su2double          *viscosityInt,
                                               su2double          *kOverCvInt);

  /*!
   * \brief Function, which computes the viscous fluxes in the integration
            points for the boundary faces that must be treated simulaneously.
            The viscous fluxes are computed via a wall modeling approach.
   * \param[in]  config                 - Definition of the particular problem.
   * \param[in]  nFaceSimul             - Number of faces that are treated simultaneously
                                          to improve performance.
   * \param[in]  NPad                   - Value of the padding parameter to obtain optimal
                                          performance in the gemm computations.
   * \param[in]  nInt                   - Number of integration points on the face.
   * \param[in]  Wall_HeatFlux          - The value of the prescribed heat flux.
   * \param[in]  HeatFlux_Prescribed    - Whether or not the heat flux is prescribed by
                                          the boundary conditions.
   * \param[in]  Wall_Temperature       - The value of the prescribed wall temperature.
   * \param[in]  Temperature_Prescribed - Whether or not the temperature is precribed
                                          by  the boundary conditions
   * \param[in]  surfElem               - Surface boundary elements for which the
                                          viscous fluxes must be computed.
   * \param[in]  solIntL                - Left states in the integration points of the face.
   * \param[out] workArray              - Storage array
   * \param[out] viscFluxes             - To be computed viscous fluxes in the
                                          integration points.
   * \param[out] viscosityInt           - To be computed viscosity in the integration points.
   * \param[out] kOverCvInt             - To be computed thermal conductivity in the
                                          integration points.
   * \param[in,out] wallModel           - Pointer to the wall model treatment.
   */
  void WallTreatmentViscousFluxes(CConfig                  *config,
                                  const unsigned short     nFaceSimul,
                                  const unsigned short     NPad,
                                  const unsigned short     nInt,
                                  const su2double          Wall_HeatFlux,
                                  const bool               HeatFlux_Prescribed,
                                  const su2double          Wall_Temperature,
                                  const bool               Temperature_Prescribed,
                                  const CSurfaceElementFEM *surfElem,
                                  const su2double          *solIntL,
                                        su2double          *workArray,
                                        su2double          *viscFluxes,
                                        su2double          *viscosityInt,
                                        su2double          *kOverCvInt,
                                        CWallModel         *wallModel);

  /*!
   * \brief Function, which computes the residual contribution from a boundary
   face in a viscous computation when the boundary conditions have
   already been applied.
   * \param[in]     config        - Definition of the particular problem.
   * \param[in]     conv_numerics - Description of the numerical method.
   * \param[in]     nFaceSimul    - Number of fused faces, i.e. the number of faces
                                    that are treated simultaneously to improve performance.
   * \param[in]     NPad          - Value of the padding parameter to obtain optimal
                                    performance in the gemm computations.
   * \param[in]     surfElem      - Surface boundary elements for which the
                                    contribution to the residual must be computed.
   * \param[in]     solInt0       - Solution in the integration points of side 0.
   * \param[in]     solInt1       - Solution in the integration points of side 1.
   * \param[out]    paramFluxes   - Array used for temporary storage.
   * \param[out]    fluxes        - Temporary storage for the fluxes in the
                                    integration points.
   * \param[in,out] viscFluxes    - On input this array contains the viscous fluxes
                                    in the integration points. It is also used for
                                    temporary storage.
   * \param[in]     viscosityInt  - Temporary storage for the viscosity in the
                                    integration points.
   * \param[in]     kOverCvInt    - Temporary storage for the thermal conductivity
                                    over Cv in the integration points.
   * \param[out]    resFaces      - Array to store the residuals of the face.
   * \param[in,out] indResFaces   - Index in resFaces, where the current residual
                                    should be stored. It is updated in the function
                                    for the next boundary element.
   */
  void ResidualViscousBoundaryFace(CConfig                  *config,
                                   CNumerics                *conv_numerics,
                                   const unsigned short     nFaceSimul,
                                   const unsigned short     NPad,
                                   const CSurfaceElementFEM *surfElem,
                                   const su2double          *solInt0,
                                   const su2double          *solInt1,
                                   su2double                *paramFluxes,
                                   su2double                *fluxes,
                                   su2double                *viscFluxes,
                                   const su2double          *viscosityInt,
                                   const su2double          *kOverCvInt,
                                   su2double                *resFaces,
                                   unsigned long            &indResFaces);

  /*!
   * \brief Function to compute the symmetrizing terms in the integration
            points of a face.
   * \param[in]  indFaceChunk      - Index of the face in the chunk of fused faces.
   * \param[in]  nInt              - Number of integration points of the face.
   * \param[in]  NPad              - Value of the padding parameter to obtain optimal
                                     performance in the gemm computations.
   * \param[in]  solInt0           - Solution in the integration points of side 0.
   * \param[in]  solInt1           - Solution in the integration points of side 1.
   * \param[in]  viscosityInt0     - Viscosity in the integration points of side 0.
   * \param[in]  viscosityInt1     - Viscosity in the integration points of side 1.
   * \param[in]  kOverCvInt0       - Heat conductivity divided by Cv in the
                                     integration points of side 0.
   * \param[in]  kOverCvInt1       - Heat conductivity divided by Cv in the
                                     integration points of side 1.
   * \param[in]  metricNormalsFace - Metric terms in the integration points, which
                                     contain the normals.
   * \param[out] symmFluxes        - Symmetrizing fluxes in the integration points.
   */
  void SymmetrizingFluxesFace(const unsigned short indFaceChunk,
                              const unsigned short nInt,
                              const unsigned short NPad,
                              const su2double      *solInt0,
                              const su2double      *solInt1,
                              const su2double      *viscosityInt0,
                              const su2double      *viscosityInt1,
                              const su2double      *kOverCvInt0,
                              const su2double      *kOverCvInt1,
                              const su2double      *metricNormalsFace,
                                    su2double      *symmFluxes);

  /*!
   * \brief Function, which transforms the symmetrizing fluxes in the integration points
            such that they are suited to be multiplied by the parametric gradients of
            the basis functions.
   * \param[in]  indFaceChunk   - Index of the face in the chunk of fused faces.
   * \param[in]  nInt           - Number of integration points of the face.
   * \param[in]  NPad           - Value of the padding parameter to obtain optimal
                                  performance in the gemm computations.
   * \param[in]  halfTheta      - Half times the theta parameter in the symmetrizing terms.
   * \param[in]  symmFluxes     - Symmetrizing fluxes to be multiplied by the Cartesian
                                  gradients of the basis functions.
   * \param[in]  weights        - Integration weights of the integration points.
   * \param[in]  metricCoorFace - Derivatives of the parametric coordinates w.r.t. the
                                  Cartesian coordinates in the integration points of
                                  the face.
   * \param[out] paramFluxes    - Parametric fluxes in the integration points.
   */
  void TransformSymmetrizingFluxes(const unsigned short indFaceChunk,
                                   const unsigned short nInt,
                                   const unsigned short NPad,
                                   const su2double      halfTheta,
                                   const su2double      *symmFluxes,
                                   const su2double      *weights,
                                   const su2double      *metricCoorFace,
                                         su2double      *paramFluxes);

  /*!
   * \brief Function to compute the viscous normal fluxes in the integration points of a face.
   * \param[in]   adjVolElem          - Pointer to the adjacent volume.
   * \param[in]   indFaceChunk        - Index of the face in the chunk of fused faces.
   * \param[in]   nInt                - Number of integration points of the face.
   * \param[in]   NPad                - Value of the padding parameter to obtain optimal
                                        performance in the gemm computations.
   * \param[in]   Wall_HeatFlux       - The value of the prescribed heat flux.
   * \param[in]   HeatFlux_Prescribed - Whether or not the heat flux is prescribed by
                                        e.g. the boundary conditions.
   * \param[in]   solInt              - Solution in the integration points.
   * \param[in]   gradSolInt          - Gradient of the solution in the integration points.
   * \param[in]   metricCoorDerivFace - Metric terms in the integration points, which
                                        contain the derivatives of the parametric
                                        coordinates w.r.t. the Cartesian coordinates.
                                        Needed to compute the Cartesian gradients.
   * \param[in]   metricNormalsFace   - Metric terms in the integration points, which
                                        contain the normals.
   * \param[in]   wallDistanceInt     - Wall distances in the integration points of the face.
   * \param[out]  viscNormFluxes      - Viscous normal fluxes in the integration points.
   * \param[out]  viscosityInt        - Viscosity in the integration points, which is
                                        needed for other terms in the discretization.
   * \param[out]  kOverCvInt          - Thermal conductivity over Cv in the integration points,
                                        which is needed for other terms in the discretization.
   */
  void ViscousNormalFluxFace(const CVolumeElementFEM *adjVolElem,
                             const unsigned short    indFaceChunk,
                             const unsigned short    nInt,
                             const unsigned short    NPad,
                             const su2double         Wall_HeatFlux,
                             const bool              HeatFlux_Prescribed,
                             const su2double         *solInt,
                             const su2double         *gradSolInt,
                             const su2double         *metricCoorDerivFace,
                             const su2double         *metricNormalsFace,
                             const su2double         *wallDistanceInt,
                                   su2double         *viscNormFluxes,
                                   su2double         *viscosityInt,
                                   su2double         *kOverCvInt);

  /*!
   * \brief Function to compute the viscous normal flux in one integration point for a
            2D simulation.
   * \param[in]  sol            - Conservative variables.
   * \param[in]  solGradCart   - Cartesian gradients of the conservative variables.
   * \param[in]  normal        - Normal vector
   * \param[in]  HeatFlux      - Value of the prescribed heat flux. If not
                                 prescribed, this value should be zero.
   * \param[in]  factHeatFlux  - Multiplication factor for the heat flux. It is zero
                                 when the heat flux is prescribed and one when it has
                                 to be computed.
   * \param[in]  wallDist      - Distance to the nearest viscous wall, if appropriate.
   * \param[in   lenScale_LES  - LES length scale, if appropriate.
   * \param[out] Viscosity     - Total viscosity, to be computed.
   * \param[out] kOverCv       - Total thermal conductivity over Cv, to be computed.
   * \param[out] normalFlux    - Viscous normal flux, to be computed.
   */
  void ViscousNormalFluxIntegrationPoint_2D(const su2double *sol,
                                            const su2double solGradCart[4][2],
                                            const su2double *normal,
                                            const su2double HeatFlux,
                                            const su2double factHeatFlux,
                                            const su2double wallDist,
                                            const su2double lenScale_LES,
                                                  su2double &Viscosity,
                                                  su2double &kOverCv,
                                                  su2double *normalFlux);

  /*!
   * \brief Function to compute the viscous normal flux in one integration point for a
            3D simulation.
   * \param[in]  sol           - Conservative variables.
   * \param[in]  solGradCart   - Cartesian gradients of the conservative variables.
   * \param[in]  normal        - Normal vector
   * \param[in]  HeatFlux      - Value of the prescribed heat flux. If not
                                 prescribed, this value should be zero.
   * \param[in]  factHeatFlux  - Multiplication factor for the heat flux. It is zero
                                 when the heat flux is prescribed and one when it has
                                 to be computed.
   * \param[in]  wallDist      - Distance to the nearest viscous wall, if appropriate.
   * \param[in   lenScale_LES  - LES length scale, if appropriate.
   * \param[out] Viscosity     - Total viscosity, to be computed.
   * \param[out] kOverCv       - Total thermal conductivity over Cv, to be computed.
   * \param[out] normalFlux    - Viscous normal flux, to be computed.
   */
  void ViscousNormalFluxIntegrationPoint_3D(const su2double *sol,
                                            const su2double solGradCart[5][3],
                                            const su2double *normal,
                                            const su2double HeatFlux,
                                            const su2double factHeatFlux,
                                            const su2double wallDist,
                                            const su2double lenScale_LES,
                                                  su2double &Viscosity,
                                                  su2double &kOverCv,
                                                  su2double *normalFlux);
};

#include "solver_structure.inl"
