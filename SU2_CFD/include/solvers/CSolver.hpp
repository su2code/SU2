/*!
 * \file CSolver.hpp
 * \brief Headers of the CSolver class which is inherited by all of the other
 *        solvers
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <set>
#include <stdlib.h>
#include <stdio.h>

#include "../fluid_model.hpp"
#include "../task_definition.hpp"
#include "../numerics_structure.hpp"
#include "../sgs_model.hpp"
#include "../../../Common/include/fem_geometry_structure.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/linear_algebra/CSysMatrix.hpp"
#include "../../../Common/include/linear_algebra/CSysVector.hpp"
#include "../../../Common/include/linear_algebra/CSysSolve.hpp"
#include "../../../Common/include/grid_movement_structure.hpp"
#include "../../../Common/include/blas_structure.hpp"
#include "../../../Common/include/graph_coloring_structure.hpp"
#include "../../../Common/include/toolboxes/MMS/CVerificationSolution.hpp"
#include "../variables/CVariable.hpp"

using namespace std;

class CSolver {
protected:
  int rank,       /*!< \brief MPI Rank. */
  size;           /*!< \brief MPI Size. */
  bool adjoint;   /*!< \brief Boolean to determine whether solver is initialized as a direct or an adjoint solver. */
  unsigned short MGLevel;        /*!< \brief Multigrid level of this solver object. */
  unsigned short IterLinSolver;  /*!< \brief Linear solver iterations. */
  su2double ResLinSolver;        /*!< \brief Final linear solver residual. */
  su2double NonLinRes_Value,        /*!< \brief Summed value of the nonlinear residual indicator. */
  NonLinRes_Func;      /*!< \brief Current value of the nonlinear residual indicator at one iteration. */
  unsigned short NonLinRes_Counter;  /*!< \brief Number of elements of the nonlinear residual indicator series. */
  vector<su2double> NonLinRes_Series;      /*!< \brief Vector holding the nonlinear residual indicator series. */
  su2double Old_Func,  /*!< \brief Old value of the nonlinear residual indicator. */
  New_Func;      /*!< \brief Current value of the nonlinear residual indicator. */
  unsigned short nVar,           /*!< \brief Number of variables of the problem. */
  nPrimVar,                      /*!< \brief Number of primitive variables of the problem. */
  nPrimVarGrad,                  /*!< \brief Number of primitive variables of the problem in the gradient computation. */
  nSecondaryVar,                 /*!< \brief Number of primitive variables of the problem. */
  nSecondaryVarGrad,             /*!< \brief Number of primitive variables of the problem in the gradient computation. */
  nVarGrad,                      /*!< \brief Number of variables for deallocating the LS Cvector. */
  nDim;                          /*!< \brief Number of dimensions of the problem. */
  unsigned long nPoint;          /*!< \brief Number of points of the computational grid. */
  unsigned long nPointDomain;    /*!< \brief Number of points of the computational grid. */
  su2double Max_Delta_Time,  /*!< \brief Maximum value of the delta time for all the control volumes. */
  Min_Delta_Time;          /*!< \brief Minimum value of the delta time for all the control volumes. */
  su2double Max_CFL_Local;  /*!< \brief Maximum value of the CFL across all the control volumes. */
  su2double Min_CFL_Local;  /*!< \brief Minimum value of the CFL across all the control volumes. */
  su2double Avg_CFL_Local;  /*!< \brief Average value of the CFL across all the control volumes. */
  su2double *Residual_RMS,  /*!< \brief Vector with the mean residual for each variable. */
  *Residual_Max,            /*!< \brief Vector with the maximal residual for each variable. */
  *Residual,                /*!< \brief Auxiliary nVar vector. */
  *Residual_i,              /*!< \brief Auxiliary nVar vector for storing the residual at point i. */
  *Residual_j;              /*!< \brief Auxiliary nVar vector for storing the residual at point j. */
  su2double *Residual_BGS,  /*!< \brief Vector with the mean residual for each variable for BGS subiterations. */
  *Residual_Max_BGS;        /*!< \brief Vector with the maximal residual for each variable for BGS subiterations. */
  unsigned long *Point_Max;        /*!< \brief Vector with the maximal residual for each variable. */
  unsigned long *Point_Max_BGS;    /*!< \brief Vector with the maximal residual for each variable. */
  su2double **Point_Max_Coord;     /*!< \brief Vector with pointers to the coords of the maximal residual for each variable. */
  su2double **Point_Max_Coord_BGS; /*!< \brief Vector with pointers to the coords of the maximal residual for each variable. */
  su2double *Solution,    /*!< \brief Auxiliary nVar vector. */
  *Solution_i,            /*!< \brief Auxiliary nVar vector for storing the solution at point i. */
  *Solution_j;            /*!< \brief Auxiliary nVar vector for storing the solution at point j. */
  su2double *Vector,  /*!< \brief Auxiliary nDim vector. */
  *Vector_i,          /*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point i. */
  *Vector_j;          /*!< \brief Auxiliary nDim vector to do the reconstruction of the variables at point j. */
  su2double *Res_Conv,  /*!< \brief Auxiliary nVar vector for storing the convective residual. */
  *Res_Visc,            /*!< \brief Auxiliary nVar vector for storing the viscous residual. */
  *Res_Sour,            /*!< \brief Auxiliary nVar vector for storing the viscous residual. */
  *Res_Conv_i,          /*!< \brief Auxiliary vector for storing the convective residual at point i. */
  *Res_Visc_i,          /*!< \brief Auxiliary vector for storing the viscous residual at point i. */
  *Res_Conv_j,          /*!< \brief Auxiliary vector for storing the convective residual at point j. */
  *Res_Visc_j;          /*!< \brief Auxiliary vector for storing the viscous residual at point j. */
  su2double **Jacobian_i,   /*!< \brief Auxiliary matrices for storing point to point Jacobians at point i. */
  **Jacobian_j;             /*!< \brief Auxiliary matrices for storing point to point Jacobians at point j. */
  su2double **Jacobian_ii,  /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  **Jacobian_ij,            /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  **Jacobian_ji,            /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  **Jacobian_jj;            /*!< \brief Auxiliary matrices for storing point to point Jacobians. */
  su2double *iPoint_UndLapl,  /*!< \brief Auxiliary variable for the undivided Laplacians. */
  *jPoint_UndLapl;            /*!< \brief Auxiliary variable for the undivided Laplacians. */
  su2double **Smatrix,        /*!< \brief Auxiliary structure for computing gradients by least-squares */
  **Cvector;                  /*!< \brief Auxiliary structure for computing gradients by least-squares */

  int *Restart_Vars;                /*!< \brief Auxiliary structure for holding the number of variables and points in a restart. */
  int Restart_ExtIter;              /*!< \brief Auxiliary structure for holding the external iteration offset from a restart. */
  passivedouble *Restart_Data;      /*!< \brief Auxiliary structure for holding the data values from a restart. */
  unsigned short nOutputVariables;  /*!< \brief Number of variables to write. */

  unsigned long nMarker,                 /*!< \brief Total number of markers using the grid information. */
  *nVertex;                              /*!< \brief Store nVertex at each marker for deallocation */

  bool rotate_periodic;    /*!< \brief Flag that controls whether the periodic solution needs to be rotated for the solver. */
  bool implicit_periodic;  /*!< \brief Flag that controls whether the implicit system should be treated by the periodic BC comms. */

  bool dynamic_grid;       /*!< \brief Flag that determines whether the grid is dynamic (moving or deforming + grid velocities). */

  su2double ***VertexTraction;          /*- Temporary, this will be moved to a new postprocessing structure once in place -*/
  su2double ***VertexTractionAdjoint;   /*- Also temporary -*/

  string SolverName;      /*!< \brief Store the name of the solver for output purposes. */

  /*!
   * \brief Pure virtual function, all derived solvers MUST implement a method returning their "nodes".
   * \note Don't forget to call SetBaseClassPointerToNodes() in the constructor of the derived CSolver.
   * \return Nodes of the solver, upcast to their base class (CVariable).
   */
  virtual CVariable* GetBaseClassPointerToNodes() = 0;

  /*!
   * \brief Call this method to set "base_nodes" after the "nodes" variable of the derived solver is instantiated.
   * \note One could set base_nodes directly if it were not private but that could lead to confusion
   */
  inline void SetBaseClassPointerToNodes() { base_nodes = GetBaseClassPointerToNodes(); }

private:

  /*--- Private to prevent use by derived solvers, each solver MUST have its own "nodes" member of the
   most derived type possible, e.g. CEulerVariable has nodes of CEulerVariable* and not CVariable*.
   This variable is to avoid two virtual functions calls per call i.e. CSolver::GetNodes() returns
   directly instead of calling GetBaseClassPointerToNodes() or doing something equivalent. ---*/
  CVariable* base_nodes;  /*!< \brief Pointer to CVariable to allow polymorphic access to solver nodes. */

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
  string* OutputHeadingNames;               /*!< \brief vector of strings to store the headings for the exra variables */

  CVerificationSolution *VerificationSolution; /*!< \brief Verification solution class used within the solver. */

  vector<string> fields;
  /*!
   * \brief Constructor of the class.
   */
  CSolver(bool mesh_deform_mode = false);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSolver(void);

  /*!
   * \brief Allow outside access to the nodes of the solver, containing conservatives, primitives, etc.
   * \return Nodes of the solver.
   */
  inline CVariable* GetNodes() {
    assert(base_nodes!=nullptr && "CSolver::base_nodes was not set properly, see brief for CSolver::SetBaseClassPointerToNodes()");
    return base_nodes;
  }

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
  inline void SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

  /*!
   * \brief Set the final linear solver residual.
   * \param[in] val_reslinsolver - Value of final linear solver residual.
   */
  inline void SetResLinSolver(su2double val_reslinsolver) { ResLinSolver = val_reslinsolver; }

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
   * \brief Move the mesh in time
   */
  inline virtual void SetDualTime_Mesh(void){ }

  /*!
   * \brief Store the BGS solution in the previous subiteration in the corresponding vector.
   */
  void UpdateSolution_BGS(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void SetNondimensionalization(CConfig *config, unsigned short iMesh) { }

  /*!
   * \brief Get information whether the initialization is an adjoint solver or not.
   * \return <code>TRUE</code> means that it is an adjoint solver.
   */
  inline bool GetAdjoint(void) const { return adjoint; }

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline virtual CFluidModel* GetFluidModel(void) const { return NULL;}

  /*!
   * \brief Get number of linear solver iterations.
   * \return Number of linear solver iterations.
   */
  inline unsigned short GetIterLinSolver(void) const { return IterLinSolver; }

  /*!
   * \brief Get the final linear solver residual.
   * \return Value of final linear solver residual.
   */
  inline su2double GetResLinSolver(void) const { return ResLinSolver; }

  /*!
   * \brief Get the value of the maximum delta time.
   * \return Value of the maximum delta time.
   */
  inline su2double GetMax_Delta_Time(void) const { return Max_Delta_Time; }

  /*!
   * \brief Get the value of the minimum delta time.
   * \return Value of the minimum delta time.
   */
  inline su2double GetMin_Delta_Time(void) const { return Min_Delta_Time; }

  /*!
   * \brief Get the value of the maximum local CFL number.
   * \return Value of the maximum local CFL number.
   */
  inline su2double GetMax_CFL_Local(void) const { return Max_CFL_Local; }

  /*!
   * \brief Get the value of the minimum local CFL number.
   * \return Value of the minimum local CFL number.
   */
  inline su2double GetMin_CFL_Local(void) const { return Min_CFL_Local; }

  /*!
   * \brief Get the value of the average local CFL number.
   * \return Value of the average local CFL number.
   */
  inline su2double GetAvg_CFL_Local(void) const { return Avg_CFL_Local; }

  /*!
   * \brief Get the number of variables of the problem.
   */
  inline unsigned short GetnVar(void) const { return nVar; }

  /*!
   * \brief Get the number of variables of the problem.
   */
  inline unsigned short GetnPrimVar(void) const { return nPrimVar; }

  /*!
   * \brief Get the number of variables of the problem.
   */
  inline unsigned short GetnPrimVarGrad(void) const { return nPrimVarGrad; }

  /*!
   * \brief Get the number of variables of the problem.
   */
  inline unsigned short GetnSecondaryVar(void) const { return nSecondaryVar; }

  /*!
   * \brief Get the number of variables of the problem.
   */
  inline unsigned short GetnSecondaryVarGrad(void) const { return nSecondaryVarGrad; }

  /*!
   * \brief Get the number of variables of the problem.
   */
  inline unsigned short GetnOutputVariables(void) const { return nOutputVariables; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  inline virtual void SetResidual_DualTime(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CConfig *config,
                                           unsigned short iRKStep,
                                           unsigned short iMesh,
                                           unsigned short RunTime_EqSystem) { }

  /*!
   * \brief Set the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  inline void SetRes_RMS(unsigned short val_var, su2double val_residual) { Residual_RMS[val_var] = val_residual; }

  /*!
   * \brief Adds the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  inline void AddRes_RMS(unsigned short val_var, su2double val_residual) { Residual_RMS[val_var] += val_residual; }

  /*!
   * \brief Get the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  inline su2double GetRes_RMS(unsigned short val_var) const { return Residual_RMS[val_var]; }

  /*!
   * \brief Set the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  inline void SetRes_Max(unsigned short val_var,
                         su2double val_residual,
                         unsigned long val_point) {
    Residual_Max[val_var] = val_residual; Point_Max[val_var] = val_point;
  }

  /*!
   * \brief Adds the maximal residual, this is useful for the convergence history (overload).
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   * \param[in] val_point - Value of the point index for the max residual.
   * \param[in] val_coord - Location (x, y, z) of the max residual point.
   */
  inline void AddRes_Max(unsigned short val_var,
                         su2double val_residual,
                         unsigned long val_point,
                         const su2double* val_coord) {
    if (val_residual > Residual_Max[val_var]) {
      Residual_Max[val_var] = val_residual;
      Point_Max[val_var] = val_point;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Point_Max_Coord[val_var][iDim] = val_coord[iDim];
    }
  }


  /*!
   * \brief Get the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  inline su2double GetRes_Max(unsigned short val_var) const { return Residual_Max[val_var]; }

  /*!
   * \brief Set the residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  inline void SetRes_BGS(unsigned short val_var, su2double val_residual) { Residual_BGS[val_var] = val_residual; }

  /*!
   * \brief Adds the residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  inline void AddRes_BGS(unsigned short val_var, su2double val_residual) { Residual_BGS[val_var] += val_residual; }

  /*!
   * \brief Get the residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  inline su2double GetRes_BGS(unsigned short val_var) const { return Residual_BGS[val_var]; }

  /*!
   * \brief Set the maximal residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   */
  inline void SetRes_Max_BGS(unsigned short val_var,
                             su2double val_residual,
                             unsigned long val_point) {
    Residual_Max_BGS[val_var] = val_residual; Point_Max_BGS[val_var] = val_point;
  }

  /*!
   * \brief Adds the maximal residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_residual - Value of the residual to store in the position <i>val_var</i>.
   * \param[in] val_point - Value of the point index for the max residual.
   * \param[in] val_coord - Location (x, y, z) of the max residual point.
   */
  inline void AddRes_Max_BGS(unsigned short val_var,
                             su2double val_residual,
                             unsigned long val_point,
                             const su2double* val_coord) {
    if (val_residual > Residual_Max_BGS[val_var]) {
    Residual_Max_BGS[val_var] = val_residual;
    Point_Max_BGS[val_var] = val_point;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Point_Max_Coord_BGS[val_var][iDim] = val_coord[iDim];
    }
  }


  /*!
   * \brief Get the maximal residual for BGS subiterations.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  inline su2double GetRes_Max_BGS(unsigned short val_var) const { return Residual_Max_BGS[val_var]; }

  /*!
   * \brief Get the residual for FEM structural analysis.
   * \param[in] val_var - Index of the variable.
   * \return Value of the residual for the variable in the position <i>val_var</i>.
   */
  inline virtual su2double GetRes_FEM(unsigned short val_var) const { return 0.0; }

  /*!
   * \brief Get the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  inline unsigned long GetPoint_Max(unsigned short val_var) const { return Point_Max[val_var]; }

  /*!
   * \brief Get the location of the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the location (x, y, z) of the biggest residual for the variable <i>val_var</i>.
   */
  inline su2double* GetPoint_Max_Coord(unsigned short val_var) const { return Point_Max_Coord[val_var]; }

  /*!
   * \brief Get the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Value of the biggest residual for the variable in the position <i>val_var</i>.
   */
  inline unsigned long GetPoint_Max_BGS(unsigned short val_var) const { return Point_Max_BGS[val_var]; }

  /*!
   * \brief Get the location of the maximal residual, this is useful for the convergence history.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the location (x, y, z) of the biggest residual for the variable <i>val_var</i>.
   */
  inline su2double* GetPoint_Max_Coord_BGS(unsigned short val_var) const { return Point_Max_Coord_BGS[val_var]; }

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
   * \brief Add External to Solution vector.
   */
  void Add_External_To_Solution();

  /*!
   * \brief Add the current Solution vector to External.
   */
  void Add_Solution_To_External();

  /*!
   * \brief Update a given cross-term with relaxation and the running total (External).
   * \param[in] config - Definition of the particular problem.
   * \param[in,out] cross_term - The cross-term being updated.
   */
  void Update_Cross_Term(CConfig *config, su2passivematrix &cross_term);

  /*!
   * \brief Compute the Green-Gauss gradient of the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction = false);

  /*!
   * \brief Compute the Least Squares gradient of the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction = false);

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
  inline virtual void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief Set the old solution variables to the current solution value for Runge-Kutta iteration.
            It is a virtual function, because for the DG-FEM solver a different version is needed.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  inline virtual void Set_OldSolution(CGeometry *geometry) { base_nodes->Set_OldSolution(); }

  /*!
   * \brief Set the new solution variables to the current solution value for classical RK.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  inline virtual void Set_NewSolution(CGeometry *geometry) { }

  /*!
   * \brief Load the geometries at the previous time states n and nM1.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void Restart_OldGeometry(CGeometry *geometry, CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  inline virtual void SetTime_Step(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CConfig *config,
                                   unsigned short iMesh,
                                   unsigned long Iteration) { }

  /*!
   * \brief A virtual member.
   * \param[in]     config          - Definition of the particular problem.
   * \param[in]     TimeSync        - The synchronization time.
   * \param[in,out] timeEvolved     - On input the time evolved before the time step,
                                      on output the time evolved after the time step.
   * \param[out]    syncTimeReached - Whether or not the synchronization time is reached.
   */
  inline virtual void CheckTimeSynchronization(CConfig         *config,
                                               const su2double TimeSync,
                                               su2double       &timeEvolved,
                                               bool            &syncTimeReached) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void ProcessTaskList_DG(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics **numerics,
                                         CConfig *config,
                                         unsigned short iMesh) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void ADER_SpaceTimeIntegration(CGeometry *geometry,
                                                CSolver **solver_container,
                                                CNumerics **numerics,
                                                CConfig *config,
                                                unsigned short iMesh,
                                                unsigned short RunTime_EqSystem) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                             CNumerics **numerics, CConfig *config,
                                             unsigned short iMesh, unsigned short RunTime_EqSystem) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void Postprocessing(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CConfig *config,
                                     unsigned short iMesh) { }

  /*!
   * \brief A virtual member, overloaded.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   *
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void Postprocessing(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CConfig *config,
                                     CNumerics **numerics,
                                     unsigned short iMesh) { }
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  inline virtual void Centered_Residual(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *numerics,
                                        CConfig *config,
                                        unsigned short iMesh,
                                        unsigned short iRKStep) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void Upwind_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics,
                                      CConfig *config,
                                      unsigned short iMesh) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  inline virtual void Convective_Residual(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *numerics,
                                          CConfig *config,
                                          unsigned short iMesh,
                                          unsigned short iRKStep) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  inline virtual void Preprocessing(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CConfig *config,
                                    unsigned short iMesh,
                                    unsigned short iRKStep,
                                    unsigned short RunTime_EqSystem,
                                    bool Output) { }

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
  inline virtual void Preprocessing(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CConfig *config,
                                    CNumerics **numerics,
                                    unsigned short iMesh,
                                    unsigned long Iteration,
                                    unsigned short RunTime_EqSystem,
                                    bool Output) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Set_MPI_ActDisk(CSolver **solver_container,
                                      CGeometry *geometry,
                                      CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Set_Heatflux_Areas(CGeometry *geometry, CConfig *config) { }

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Evaluate_ObjFunc(CConfig *config) {};

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Euler_Wall(CGeometry      *geometry,
                                    CSolver        **solver_container,
                                    CNumerics      *conv_numerics,
                                    CNumerics      *visc_numerics,
                                    CConfig        *config,
                                    unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_Clamped(CGeometry *geometry,
                                 CNumerics *numerics,
                                 CConfig *config,
                                 unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_Clamped_Post(CGeometry *geometry,
                                      CNumerics *numerics,
                                      CConfig *config,
                                      unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_DispDir(CGeometry *geometry,
                                 CNumerics *numerics,
                                 CConfig *config,
                                 unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_Normal_Displacement(CGeometry *geometry,
                                             CNumerics *numerics,
                                             CConfig *config,
                                             unsigned short val_marker) { }


  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Normal_Load(CGeometry *geometry,
                                     CNumerics *numerics,
                                     CConfig *config,
                                     unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_Dir_Load(CGeometry *geometry,
                                  CNumerics *numerics,
                                  CConfig *config,
                                  unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_Sine_Load(CGeometry *geometry,
                                   CNumerics *numerics,
                                   CConfig *config,
                                   unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Damper(CGeometry *geometry,
                                CNumerics *numerics,
                                CConfig *config,
                                unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  inline virtual void BC_Deforming(CGeometry *geometry,
                                   CNumerics *numerics,
                                   CConfig *config,
                                   unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Interface_Boundary(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics *numerics,
                                            CConfig *config,
                                            unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_NearField_Boundary(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics *numerics,
                                            CConfig *config,
                                            unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void BC_Periodic(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *numerics,
                                  CConfig *config) { }

  /*!
   * \brief Impose the interface state across sliding meshes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void BC_Fluid_Interface(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics *conv_numerics,
                                         CNumerics *visc_numerics,
                                         CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_ActDisk_Inlet(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_ActDisk_Outlet(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *conv_numerics,
                                        CNumerics *visc_numerics,
                                        CConfig *config,
                                        unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  inline virtual void BC_ActDisk(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker,
                                 bool val_inlet_surface) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Isothermal_Wall(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics *conv_numerics,
                                         CNumerics *visc_numerics,
                                         CConfig *config,
                                         unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_HeatFlux_Wall(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Far_Field(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CNumerics *conv_numerics,
                                   CNumerics *visc_numerics,
                                   CConfig *config,
                                   unsigned short val_marker) { }

  /*!
   * \brief Impose via the residual the Euler boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Sym_Plane(CGeometry      *geometry,
                                   CSolver        **solver_container,
                                   CNumerics      *conv_numerics,
                                   CNumerics      *visc_numerics,
                                   CConfig        *config,
                                   unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Riemann(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_TurboRiemann(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *conv_numerics,
                                      CNumerics *visc_numerics,
                                      CConfig *config,
                                      unsigned short val_marker) { }

  /*!
   * \brief It computes Fourier transformation for the needed quantities along the pitch for each span in turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  inline virtual void PreprocessBC_Giles(CGeometry *geometry,
                                         CConfig *config,
                                         CNumerics *conv_numerics,
                                         unsigned short marker_flag) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Giles(CGeometry *geometry,
                               CSolver **solver_container,
                               CNumerics *conv_numerics,
                               CNumerics *visc_numerics,
                               CConfig *config,
                               unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Inlet(CGeometry *geometry,
                               CSolver **solver_container,
                               CNumerics *conv_numerics,
                               CNumerics *visc_numerics,
                               CConfig *config,
                               unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Inlet_Turbo(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *visc_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) { }
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Inlet_MixingPlane(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *conv_numerics,
                                           CNumerics *visc_numerics,
                                           CConfig *config,
                                           unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Supersonic_Inlet(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *conv_numerics,
                                          CNumerics *visc_numerics,
                                          CConfig *config,
                                          unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Supersonic_Outlet(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *conv_numerics,
                                           CNumerics *visc_numerics,
                                           CConfig *config,
                                           unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry         - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics    - Description of the convective numerical method.
   * \param[in] visc_numerics    - Description of the viscous numerical method.
   * \param[in] config           - Definition of the particular problem.
   * \param[in] val_marker       - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Custom(CGeometry *geometry,
                                CSolver **solver_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics,
                                CConfig *config,
                                unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Outlet(CGeometry *geometry,
                                CSolver **solver_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics,
                                CConfig *config,
                                unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Engine_Inflow(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Engine_Exhaust(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *conv_numerics,
                                        CNumerics *visc_numerics,
                                        CConfig *config,
                                        unsigned short val_marker) { }

  /*!
   * \brief Impose the symmetry boundary condition using the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Dielec(CGeometry *geometry,
                                CSolver **solver_container,
                                CNumerics *numerics,
                                CConfig *config,
                                unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_Electrode(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CNumerics *numerics,
                                   CConfig *config,
                                   unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline virtual void BC_ConjugateHeat_Interface(CGeometry *geometry,
                                                 CSolver **solver_container,
                                                 CNumerics *numerics,
                                                 CConfig *config,
                                                 unsigned short val_marker) { }

 /*!
   * \brief Get the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   * \param[in] donor_index- index of the donor node to get
   */
  inline virtual su2double GetSlidingState(unsigned short val_marker,
                                           unsigned long val_vertex,
                                           unsigned short val_state,
                                           unsigned long donor_index) const { return 0; }

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  inline virtual void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){}

  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   * \param[in] donor_index- index of the donor node to set
   * \param[in] component  - set value
   */
  inline virtual void SetSlidingState(unsigned short val_marker,
                                      unsigned long val_vertex,
                                      unsigned short val_state,
                                      unsigned long donor_index,
                                      su2double component){ }

  /*!
   * \brief Get the number of outer states for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  inline virtual int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief Set the number of outer states for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value      - number of outer states
   */
  inline virtual void SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value) { }

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   * \param[in] relaxation factor - relaxation factor for the change of the variables
   * \param[in] val_var           - value of the variable
   */
  inline virtual void SetConjugateHeatVariable(unsigned short val_marker,
                                               unsigned long val_vertex,
                                               unsigned short pos_var,
                                               su2double relaxation_factor,
                                               su2double val_var) { }

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   */
  inline virtual su2double GetConjugateHeatVariable(unsigned short val_marker,
                                                    unsigned long val_vertex,
                                                    unsigned short pos_var) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  inline virtual void ExplicitRK_Iteration(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CConfig *config,
                                           unsigned short iRKStep) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  inline virtual void ClassicalRK4_Iteration(CGeometry *geometry,
                                             CSolver **solver_container,
                                             CConfig *config,
                                             unsigned short iRKStep) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ExplicitEuler_Iteration(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ImplicitEuler_Iteration(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ComputeUnderRelaxationFactor(CSolver **solver_container, CConfig *config) { }

  /*!
   * \brief Adapt the CFL number based on the local under-relaxation parameters
   *        computed for each nonlinear iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver_container - Container vector with all the solutions.
   */
  void AdaptCFLNumber(CGeometry **geometry, CSolver ***solver_container, CConfig *config);

  /*!
   * \brief Reset the local CFL adaption variables
   */
  void ResetCFLAdapt();

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ImplicitNewmark_Iteration(CGeometry *geometry,
                                                CSolver **solver_container,
                                                CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ImplicitNewmark_Update(CGeometry *geometry,
                                             CSolver **solver_container,
                                             CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ImplicitNewmark_Relaxation(CGeometry *geometry,
                                                 CSolver **solver_container,
                                                 CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void GeneralizedAlpha_Iteration(CGeometry *geometry,
                                                 CSolver **solver_container,
                                                 CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void GeneralizedAlpha_UpdateDisp(CGeometry *geometry,
                                                  CSolver **solver_container,
                                                  CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void GeneralizedAlpha_UpdateSolution(CGeometry *geometry,
                                                      CSolver **solver_container,
                                                      CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void GeneralizedAlpha_UpdateLoads(CGeometry *geometry,
                                                   CSolver **solver_container,
                                                   CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Pressure_Forces(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Momentum_Forces(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Friction_Forces(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Buffet_Monitoring(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Heat_Fluxes(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  inline virtual void SetPrimitive_Gradient_GG(CGeometry *geometry,
                                               CConfig *config,
                                               bool reconstruction = false) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  inline virtual void SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction = false) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetPrimitive_Limiter_MPI(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Index of the grid point.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetPreconditioner(CConfig *config, unsigned long iPoint) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  inline virtual void Viscous_Residual(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *numerics,
                                       CConfig *config,
                                       unsigned short iMesh,
                                       unsigned short iRKStep) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void Source_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics,
                                      CNumerics *second_numerics,
                                      CConfig *config,
                                      unsigned short iMesh) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  inline virtual void Source_Template(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics,
                                      CConfig *config,
                                      unsigned short iMesh) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \param[in] val_sensitivity - Value of the sensitivity coefficient.
   */
  inline virtual void SetCSensitivity(unsigned short val_marker,
                                      unsigned long val_vertex,
                                      su2double val_sensitivity) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetForceProj_Vector(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetIntBoundary_Jump(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_CD(su2double val_Total_CD) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  inline virtual void SetTotal_CL(su2double val_Total_CL) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_NetThrust(su2double val_Total_NetThrust) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_Power(su2double val_Total_Power) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_SolidCD(su2double val_Total_SolidCD) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_MFR(su2double val_Total_MFR) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_IDC(su2double val_Total_IDC) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_IDR(su2double val_Total_IDR) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline virtual void SetTotal_DC60(su2double val_Total_DC60) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  inline virtual void SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  inline virtual void AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CT - Value of the total thrust coefficient.
   */
  inline virtual void SetTotal_CT(su2double val_Total_CT) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_CQ - Value of the total torque coefficient.
   */
  inline virtual void SetTotal_CQ(su2double val_Total_CQ) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_Heat - Value of the total heat load.
   */
  inline virtual void SetTotal_HeatFlux(su2double val_Total_Heat) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_Total_MaxHeat - Value of the total heat load.
   */
  inline virtual void SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Inviscid_Sensitivity(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *numerics,
                                           CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Smooth_Sensitivity(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CNumerics *numerics,
                                         CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Viscous_Sensitivity(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *numerics,
                                          CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCL_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCL_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CL(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CD(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CSF(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CEff(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFx(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFy(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFz(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMx(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMy(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMz(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CL_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CD_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CSF_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CEff_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFx_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFy_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFz_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMx_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMy_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMz_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CL_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CD_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CSF_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CEff_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFx_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFy_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFz_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMx_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMy_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMz_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the buffet metric on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_Buffet_Metric(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CL_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CD_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CSF_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CEff_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFx_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFy_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CFz_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMx_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMy_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_CMz_Mnt(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCSF_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCD_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetInflow_MassFlow(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  inline virtual void GetPower_Properties(CGeometry *geometry,
                                          CConfig *config,
                                          unsigned short iMesh,
                                          bool Output) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  inline virtual void GetOutlet_Properties(CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short iMesh,
                                           bool Output) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  inline virtual void SetFarfield_AoA(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CConfig *config,
                                      unsigned short iMesh,
                                      bool Output) { }

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] convergence - boolean for whether the solution is converged
   * \return boolean for whether the Fixed C_L mode is converged to target C_L
   */
  inline virtual bool FixedCL_Convergence(CConfig *config, bool convergence) { return false; }

  /*!
   * \brief A virtual member.
   * \return boolean for whether the Fixed C_L mode is currently in finite-differencing mode
   */
  inline virtual bool GetStart_AoA_FD(void) const { return false; }

  /*!
   * \brief A virtual member.
   * \return boolean for whether the Fixed C_L mode is currently in finite-differencing mode
   */
  inline virtual bool GetEnd_AoA_FD(void) const { return false; }

  /*!
   * \brief A virtual member.
   * \return value for the last iteration that the AoA was updated
   */
  inline virtual unsigned long GetIter_Update_AoA(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return value of the AoA before most recent update
   */
  inline virtual su2double GetPrevious_AoA(void) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return value of CL Driver control command (AoA_inc)
   */
  inline virtual su2double GetAoA_inc(void) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  inline virtual void SetActDisk_BCThrust(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config,
                                          unsigned short iMesh,
                                          bool Output) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetExhaust_MassFlow(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the fan face pressure on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetInflow_Pressure(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the fan face mach on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetInflow_Mach(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCSF_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCEff_Inv(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the integrated heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_HF_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the heat flux is computed.
   * \return Value of the maximum heat flux (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetSurface_MaxHF_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetCD_Visc(unsigned short val_marker) const { return 0; }

  /*!
   * \author H. Kline
   * \brief Set the total "combo" objective (weighted sum of other values).
   * \param[in] ComboObj - Value of the combined objective.
   */
  inline virtual void SetTotal_ComboObj(su2double ComboObj) {}

  /*!
   * \author H. Kline
   * \brief Provide the total "combo" objective (weighted sum of other values).
   * \return Value of the "combo" objective values.
   */
  inline virtual su2double GetTotal_ComboObj(void) const { return 0;}

  /*!
   * \brief A virtual member.
   * \return Value of the sideforce coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CSF() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CEff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the thrust coefficient (force in the -x direction, inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CT() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the torque coefficient (moment in the -x direction, inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CQ() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the heat load (integrated heat flux).
   */
  inline virtual su2double GetTotal_HeatFlux() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the heat load (integrated heat flux).
   */
  inline virtual su2double GetTotal_MaxHeatFlux() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the average temperature.
   */
  inline virtual su2double GetTotal_AvgTemperature() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the rotor Figure of Merit (FM) (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CMerit() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CEquivArea() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Aero drag (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_AeroCD() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the difference of the presure and the target pressure.
   */
  inline virtual su2double GetTotal_CpDiff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the difference of the heat and the target heat.
   */
  inline virtual su2double GetTotal_HeatFluxDiff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the FEA coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CFEA() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Near-Field Pressure coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CNearFieldOF() const { return 0; }

  /*!
   * \author H. Kline
   * \brief Add to the value of the total 'combo' objective.
   * \param[in] val_obj - Value of the contribution to the 'combo' objective.
   */
  inline virtual void AddTotal_ComboObj(su2double val_obj) {}

  /*!
   * \brief A virtual member.
   * \return Value of the objective function for a reference geometry.
   */
  inline virtual su2double GetTotal_OFRefGeom() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the objective function for a reference node.
   */
  inline virtual su2double GetTotal_OFRefNode() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the objective function for the volume fraction.
   */
  inline virtual su2double GetTotal_OFVolFrac() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the objective function for the structural compliance.
   */
  inline virtual su2double GetTotal_OFCompliance() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Bool that defines whether the solution has an element-based file or not
   */
  inline virtual bool IsElementBased(void) const { return false; }

  /*!
   * \brief A virtual member.
   * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
   */
  inline virtual void SetTotal_CEquivArea(su2double val_cequivarea) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_aerocd - Value of the aero drag.
   */
  inline virtual void SetTotal_AeroCD(su2double val_aerocd) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_pressure - Value of the difference between pressure and the target pressure.
   */
  inline virtual void SetTotal_CpDiff(su2double val_pressure) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_pressure - Value of the difference between heat and the target heat.
   */
  inline virtual void SetTotal_HeatFluxDiff(su2double val_heat) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_cfea - Value of the FEA coefficient.
   */
  inline virtual void SetTotal_CFEA(su2double val_cfea) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference geometry.
   */
  inline virtual void SetTotal_OFRefGeom(su2double val_ofrefgeom) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference node.
   */
  inline virtual void SetTotal_OFRefNode(su2double val_ofrefnode) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
   */
  inline virtual void SetTotal_CNearFieldOF(su2double val_cnearfieldpress) { }

  /*!
   * \brief A virtual member.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CL() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CD() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_NetThrust() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Power() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_SolidCD() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_ReverseFlow() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_MFR() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Prop_Eff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_ByPassProp_Eff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Adiab_Eff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Poly_Eff() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_IDC() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_IDC_Mach() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_IDR() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_DC60() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the custom objective function.
   */
  inline virtual su2double GetTotal_Custom_ObjFunc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CMx() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CMy() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CMz() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CoPx() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CoPy() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CoPz() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the force x coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CFx() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the force y coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CFy() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the force y coefficient (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_CFz() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CL_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CD_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CSF_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CEff_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMx_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMy_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMz_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPx_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPy_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPz_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFx_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFy_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFz_Inv() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CL_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CD_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CSF_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CEff_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMx_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMy_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMz_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPx_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPy_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPz_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFx_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFy_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFz_Visc() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CL_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CD_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CSF_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CEff_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMx_Mnt() const { return 0; }
  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMy_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CMz_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPx_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPy_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CoPz_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFx_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFy_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline virtual su2double GetAllBound_CFz_Mnt() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the buffet metric.
   */
  inline virtual su2double GetTotal_Buffet_Metric() const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double GetCPressure(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual void SetCPressureTarget(unsigned short val_marker,
                                         unsigned long val_vertex,
                                         su2double val_pressure) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */

  inline virtual void SetCharacPrimVar(unsigned short val_marker,
                                       unsigned long val_vertex,
                                       unsigned short val_var,
                                       su2double val_value) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.

   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double *GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual void SetDonorPrimVar(unsigned short val_marker,
                                      unsigned long val_vertex,
                                      unsigned short val_var,
                                      su2double val_value) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */

  inline virtual void SetDonorAdjVar(unsigned short val_marker,
                                     unsigned long val_vertex,
                                     unsigned short val_var,
                                     su2double val_value) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double GetDonorPrimVar(unsigned short val_marker,
                                           unsigned long val_vertex,
                                           unsigned short val_var) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */

  inline virtual su2double *GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex) const {
    return nullptr;
  }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */

  inline virtual su2double GetDonorAdjVar(unsigned short val_marker,
                                          unsigned long val_vertex,
                                          unsigned short val_var) const {
    return 0;
  }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual unsigned long GetDonorGlobalIndex(unsigned short val_marker,
                                                   unsigned long val_vertex) const {
    return 0;
  }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual void SetDonorGlobalIndex(unsigned short val_marker,
                                          unsigned long val_vertex,
                                          unsigned long val_index) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double *GetCharacPrimVar(unsigned short val_marker,
                                             unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual void SetActDisk_DeltaP(unsigned short val_marker,
                                        unsigned long val_vertex,
                                        su2double val_deltap) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual su2double GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex) { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual void SetActDisk_DeltaT(unsigned short val_marker,
                                        unsigned long val_vertex,
                                        su2double val_deltat) { }

  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total temperature is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is evaluated.
   * \return Value of the total temperature
   */
  inline virtual su2double GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total pressure is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is evaluated.
   * \return Value of the total pressure
   */
  inline virtual su2double GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the flow direction is evaluated
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is evaluated
   * \param[in] val_dim - The component of the flow direction unit vector to be evaluated
   * \return Component of a unit vector representing the flow direction.
   */
  inline virtual su2double GetInlet_FlowDir(unsigned short val_marker,
                                            unsigned long val_vertex,
                                            unsigned short val_dim) const { return 0; }

  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total temperature is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is set.
   * \param[in] val_ttotal - Value of the total temperature
   */
  inline virtual void SetInlet_Ttotal(unsigned short val_marker,
                                      unsigned long val_vertex,
                                      su2double val_ttotal) { }

  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the total pressure is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is set.
   * \param[in] val_ptotal - Value of the total pressure
   */
  inline virtual void SetInlet_Ptotal(unsigned short val_marker,
                                      unsigned long val_vertex,
                                      su2double val_ptotal) { }

  /*!
   * \brief A virtual member
   * \param[in] val_marker - Surface marker where the flow direction is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is set.
   * \param[in] val_dim - The component of the flow direction unit vector to be set
   * \param[in] val_flowdir - Component of a unit vector representing the flow direction.
   */
  inline virtual void SetInlet_FlowDir(unsigned short val_marker,
                                       unsigned long val_vertex,
                                       unsigned short val_dim,
                                       su2double val_flowdir) { }

  /*!
   * \brief A virtual member
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the turbulence variable (i.e. k is 0 in SST)
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  inline virtual void SetInlet_TurbVar(unsigned short val_marker,
                                       unsigned long val_vertex,
                                       unsigned short val_dim,
                                       su2double val_turb_var) { }

  /*!
   * \brief A virtual member
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  inline virtual void SetUniformInlet(CConfig* config, unsigned short iMarker) {};

  /*!
   * \brief A virtual member
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  inline virtual void SetInletAtVertex(su2double *val_inlet,
                                       unsigned short iMarker,
                                       unsigned long iVertex) { };

  /*!
   * \brief A virtual member
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  inline virtual su2double GetInletAtVertex(su2double *val_inlet,
                                            unsigned long val_inlet_point,
                                            unsigned short val_kind_marker,
                                            string val_marker,
                                            CGeometry *geometry,
                                            CConfig *config) const { return 0; }

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  inline virtual void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the skin friction coefficient.
   */
  inline virtual su2double GetCSkinFriction(unsigned short val_marker,
                                            unsigned long val_vertex,
                                            unsigned short val_dim) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the heat transfer coefficient.
   */
  inline virtual su2double GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the heat transfer coefficient.
   */
  inline virtual su2double GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline virtual void SetHeatFluxTarget(unsigned short val_marker,
                                        unsigned long val_vertex,
                                        su2double val_heat) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the buffet sensor.
   */
  inline virtual su2double GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the y plus.
   */
  inline virtual su2double GetYPlus(unsigned short val_marker, unsigned long val_vertex) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  inline virtual su2double GetStrainMag_Max(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  inline virtual su2double GetOmega_Max(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the StrainMag_Max
   */
  inline virtual void SetStrainMag_Max(su2double val_strainmag_max) { }

  /*!
   * \brief A virtual member.
   * \return Value of the Omega_Max
   */
  inline virtual void SetOmega_Max(su2double val_omega_max) { }

  /*!
   * \brief A virtual member.
   * \return Value of the adjoint density at the infinity.
   */
  inline virtual su2double GetPsiRho_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the adjoint energy at the infinity.
   */
  inline virtual su2double GetPsiE_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_dim - Index of the adjoint velocity vector.
   * \return Value of the adjoint velocity vector at the infinity.
   */
  inline virtual su2double GetPhi_Inf(unsigned short val_dim) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the geometrical sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Sens_Geo() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Mach sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Sens_Mach() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the angle of attack sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Sens_AoA() const { return 0; }

  /*!
   * \brief Set the total farfield pressure sensitivity coefficient.
   * \return Value of the farfield pressure sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Sens_Press() const { return 0; }

  /*!
   * \brief Set the total farfield temperature sensitivity coefficient.
   * \return Value of the farfield temperature sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Sens_Temp() const { return 0; }

  /*!
   * \author H. Kline
   * \brief Get the total back pressure sensitivity coefficient.
   * \return Value of the back pressure sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline virtual su2double GetTotal_Sens_BPress() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the density sensitivity.
   */
  inline virtual su2double GetTotal_Sens_Density() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the velocity magnitude sensitivity.
   */
  inline virtual su2double GetTotal_Sens_ModVel() const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the density at the infinity.
   */
  inline virtual su2double GetDensity_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the velocity at the infinity.
   */
  inline virtual su2double GetModVelocity_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the density x energy at the infinity.
   */
  inline virtual su2double GetDensity_Energy_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the pressure at the infinity.
   */
  inline virtual su2double GetPressure_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_dim - Index of the adjoint velocity vector.
   * \return Value of the density x velocity at the infinity.
   */
  inline virtual su2double GetDensity_Velocity_Inf(unsigned short val_dim) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_dim - Index of the velocity vector.
   * \return Value of the velocity at the infinity.
   */
  inline virtual su2double GetVelocity_Inf(unsigned short val_dim) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the velocity at the infinity.
   */
  inline virtual su2double *GetVelocity_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the viscosity at the infinity.
   */
  inline virtual su2double GetViscosity_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of nu tilde at the far-field.
   */
  inline virtual su2double GetNuTilde_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the turbulent kinetic energy.
   */
  inline virtual su2double GetTke_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the turbulent frequency.
   */
  inline virtual su2double GetOmega_Inf(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Young Modulus E
   */
  inline virtual su2double GetTotal_Sens_E(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity for the Poisson's ratio Nu
   */
  inline virtual su2double GetTotal_Sens_Nu(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the structural density sensitivity
   */
  inline virtual su2double GetTotal_Sens_Rho(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the structural weight sensitivity
   */
  inline virtual su2double GetTotal_Sens_Rho_DL(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField
   */
  inline virtual su2double GetTotal_Sens_EField(unsigned short iEField) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the FEA DV in the region iDVFEA
   */
  inline virtual su2double GetTotal_Sens_DVFEA(unsigned short iDVFEA) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Young Modulus E
   */
  inline virtual su2double GetGlobal_Sens_E(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Poisson's ratio Nu
   */
  inline virtual su2double GetGlobal_Sens_Nu(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the structural density sensitivity
   */
  inline virtual su2double GetGlobal_Sens_Rho(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the structural weight sensitivity
   */
  inline virtual su2double GetGlobal_Sens_Rho_DL(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the Electric Field in the region iEField
   */
  inline virtual su2double GetGlobal_Sens_EField(unsigned short iEField) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the sensitivity coefficient for the FEA DV in the region iDVFEA
   */
  inline virtual su2double GetGlobal_Sens_DVFEA(unsigned short iDVFEA) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Young modulus from the adjoint solver
   */
  inline virtual su2double GetVal_Young(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the Poisson's ratio from the adjoint solver
   */
  inline virtual su2double GetVal_Poisson(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the density for inertial effects, from the adjoint solver
   */
  inline virtual su2double GetVal_Rho(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Value of the density for dead loads, from the adjoint solver
   */
  inline virtual su2double GetVal_Rho_DL(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Number of electric field variables from the adjoint solver
   */
  inline virtual unsigned short GetnEField(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Number of design variables from the adjoint solver
   */
  inline virtual unsigned short GetnDVFEA(void) const { return 0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void ReadDV(CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \return Pointer to the values of the Electric Field
   */
  inline virtual su2double GetVal_EField(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \return Pointer to the values of the design variables
   */
  inline virtual su2double GetVal_DVFEA(unsigned short iVal) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the sensitivity coefficient.
   */
  inline virtual su2double GetCSensitivity(unsigned short val_marker,
                                           unsigned long val_vertex) const {
    return 0;
  }

  /*!
   * \brief A virtual member.
   * \return A pointer to an array containing a set of constants
   */
  inline virtual su2double* GetConstants() const { return NULL; }

  /*!
   * \brief A virtual member.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_forcecoeff_history - Value of the force coefficient.
   */
  inline virtual void SetForceCoeff(su2double val_forcecoeff_history) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_relaxcoeff_history - Value of the force coefficient.
   */
  inline virtual void SetRelaxCoeff(su2double val_relaxcoeff_history) { }

  /*!
   * \brief A virtual member.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  inline virtual void SetFSI_Residual(su2double val_FSI_residual) { }

  /*!
   * \brief A virtual member.
   * \param[out] val_forcecoeff_history - Value of the force coefficient.
   */
  inline virtual su2double GetForceCoeff(void) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[out] val_relaxcoeff_history - Value of the relax coefficient.
   */
  inline virtual su2double GetRelaxCoeff(void) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[out] val_FSI_residual - Value of the residual.
   */
  inline virtual su2double GetFSI_Residual(void) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  inline virtual void SetInitialCondition(CGeometry **geometry,
                                          CSolver ***solver_container,
                                          CConfig *config,
                                          unsigned long TimeIter) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  inline virtual void ResetInitialCondition(CGeometry **geometry,
                                            CSolver ***solver_container,
                                            CConfig *config,
                                            unsigned long TimeIter) { }

  /*!
   * \brief A virtual member.
   * \param[in] fea_geometry - Geometrical definition of the problem.
   * \param[in] fea_config - Geometrical definition of the problem.
   * \param[in] fea_geometry - Definition of the particular problem.
   */
  inline virtual void PredictStruct_Displacement(CGeometry **fea_geometry,
                                                 CConfig *fea_config,
                                                 CSolver ***fea_solution) { }

  /*!
   * \brief A virtual member.
   * \param[in] fea_geometry - Geometrical definition of the problem.
   * \param[in] fea_config - Geometrical definition of the problem.
   * \param[in] fea_geometry - Definition of the particular problem.
   */
  inline virtual void ComputeAitken_Coefficient(CGeometry **fea_geometry,
                                                CConfig *fea_config,
                                                CSolver ***fea_solution,
                                                unsigned long iOuterIter) { }


  /*!
   * \brief A virtual member.
   * \param[in] fea_geometry - Geometrical definition of the problem.
   * \param[in] fea_config - Geometrical definition of the problem.
   * \param[in] fea_geometry - Definition of the particular problem.
   */
  inline virtual void SetAitken_Relaxation(CGeometry **fea_geometry,
                                           CConfig *fea_config,
                                           CSolver ***fea_solution) { }

  /*!
   * \brief A virtual member.
   * \param[in] fea_geometry - Geometrical definition of the problem.
   * \param[in] fea_config - Geometrical definition of the problem.
   * \param[in] fea_geometry - Definition of the particular problem.
   */
  inline virtual void Update_StructSolution(CGeometry **fea_geometry,
                                            CConfig *fea_config,
                                            CSolver ***fea_solution) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  inline virtual void LoadRestart(CGeometry **geometry,
                                  CSolver ***solver,
                                  CConfig *config,
                                  int val_iter,
                                  bool val_update_geo) { }

  /*!
   * \brief Read a native SU2 restart file in ASCII format.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_SU2_Restart_ASCII(CGeometry *geometry,
                              CConfig *config,
                              string val_filename);

  /*!
   * \brief Read a native SU2 restart file in binary format.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_SU2_Restart_Binary(CGeometry *geometry,
                               CConfig *config,
                               string val_filename);

  /*!
   * \brief Read the metadata from a native SU2 restart file (ASCII or binary).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] adjoint - Boolean to identify the restart file of an adjoint run.
   * \param[in] val_filename - String name of the restart file.
   */
  void Read_SU2_Restart_Metadata(CGeometry *geometry,
                                 CConfig *config,
                                 bool adjoint_run,
                                 string val_filename);

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
  inline virtual void Compute_OFRefGeom(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_OFRefNode(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_OFVolFrac(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_OFCompliance(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Stiffness_Penalty(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics **numerics_container,
                                        CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   */
  inline virtual void LoadRestart_FSI(CGeometry *geometry,
                                      CConfig *config,
                                      int val_iter) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void RefGeom_Sensitivity(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config){ }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void DE_Sensitivity(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics **numerics_container,
                                     CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Stiffness_Sensitivity(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics **numerics_container,
                                            CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] iElem - element parameter.
   * \param[out] iElem_iDe - ID of the Dielectric Elastomer region.
   */
  inline virtual unsigned short Get_iElem_iDe(unsigned long iElem) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] i_DV - number of design variable.
   * \param[in] val_EField - value of the design variable.
   */
  inline virtual void Set_DV_Val(su2double val_EField, unsigned short i_DV) { }

  /*!
   * \brief A virtual member.
   * \param[in] i_DV - number of design variable.
   * \param[out] DV_Val - value of the design variable.
   */
  inline virtual su2double Get_DV_Val(unsigned short i_DV) const { return 0.0; }

  /*!
   * \brief Gauss method for solving a linear system.
   * \param[in] A - Matrix Ax = b.
   * \param[in] rhs - Right hand side.
   * \param[in] nVar - Number of variables.
   */
  void Gauss_Elimination(su2double** A,
                         su2double* rhs,
                         unsigned short nVar);

  /*!
   * \brief Prepares and solves the aeroelastic equations.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] TimeIter - Physical iteration number.
   */
  void Aeroelastic(CSurfaceMovement *surface_movement,
                   CGeometry *geometry,
                   CConfig *config,
                   unsigned long TimeIter);

  /*!
   * \brief Sets up the generalized eigenvectors and eigenvalues needed to solve the aeroelastic equations.
   * \param[in] PHI - Matrix of the generalized eigenvectors.
   * \param[in] w - The eigenvalues of the generalized eigensystem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpTypicalSectionWingModel(vector<vector<su2double> >& PHI,
                                    vector<su2double>& w,
                                    CConfig *config);

  /*!
   * \brief Solve the typical section wing model.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] Cl - Coefficient of lift at particular iteration.
   * \param[in] Cm - Moment coefficient about z-axis at particular iteration.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_Marker - Surface that is being monitored.
   * \param[in] displacements - solution of typical section wing model.
   */
  void SolveTypicalSectionWingModel(CGeometry *geometry,
                                    su2double Cl, su2double Cm,
                                    CConfig *config,
                                    unsigned short val_Marker,
                                    vector<su2double>& displacements);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config_container - The particular config.
   */
  inline virtual void RegisterSolution(CGeometry *geometry_container, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config_container - The particular config.
   */
  inline virtual void RegisterOutput(CGeometry *geometry_container, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  inline virtual void SetAdjoint_Output(CGeometry *geometry, CConfig *config){}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  inline virtual void SetAdjoint_OutputMesh(CGeometry *geometry, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  inline virtual void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  inline virtual void ExtractAdjoint_Geometry(CGeometry *geometry, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  inline virtual void ExtractAdjoint_CrossTerm(CGeometry *geometry, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  inline virtual void ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] solver_container - The solver container holding all solutions.
   * \param[in] config - The particular config.
   */
  inline virtual void ExtractAdjoint_CrossTerm_Geometry_Flow(CGeometry *geometry, CConfig *config) {}

  /*!
   * \brief A virtual member
   * \param[in] geometry - The geometrical definition of the problem.
   */
  inline virtual void RegisterObj_Func(CConfig *config){}

  /*!
   * \brief  A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetSurface_Sensitivity(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member. Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The solver container holding all terms of the solution.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config){ }

  /*!
   * \brief A virtual member. Extract and set the derivative of objective function.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetAdj_ObjFunc(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] Set value of interest: 0 - Initial value, 1 - Current value.
   */
  inline virtual void SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { };

  /*!
   * \brief A virtual member.
   * \param[in]  Value of interest: 0 - Initial value, 1 - Current value.
   * \return Values to compare
   */
  inline virtual su2double GetFSI_ConvValue(unsigned short val_index)const  { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] CurrentTime - Current time step.
   * \param[in] RampTime - Time for application of the ramp.*
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual su2double Compute_LoadCoefficient(su2double CurrentTime,
                                                   su2double RampTime,
                                                   CConfig *config) { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_StiffMatrix(CGeometry *geometry,
                                          CNumerics **numerics,
                                          CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_StiffMatrix_NodalStressRes(CGeometry *geometry,
                                                         CNumerics **numerics,
                                                         CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_MassMatrix(CGeometry *geometry,
                                         CNumerics **numerics,
                                         CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_MassRes(CGeometry *geometry,
                                      CNumerics **numerics,
                                      CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_NodalStressRes(CGeometry *geometry,
                                             CNumerics **numerics,
                                             CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */

  inline virtual void Compute_NodalStress(CGeometry *geometry,
                                          CNumerics **numerics,
                                          CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Compute_DeadLoad(CGeometry *geometry,
                                       CNumerics **numerics,
                                       CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void Solve_System(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \return Value of the dynamic Aitken relaxation factor
   */
  inline virtual su2double GetWAitken_Dyn(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \return Value of the last Aitken relaxation factor in the previous time step.
   */
  inline virtual su2double GetWAitken_Dyn_tn1(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] Value of the dynamic Aitken relaxation factor
   */
  inline virtual void SetWAitken_Dyn(su2double waitk) {  }

  /*!
   * \brief A virtual member.
   * \param[in] Value of the last Aitken relaxation factor in the previous time step.
   */
  inline virtual void SetWAitken_Dyn_tn1(su2double waitk_tn1) {  }

  /*!
   * \brief A virtual member.
   * \param[in] Value of the load increment for nonlinear structural analysis
   */
  inline virtual void SetLoad_Increment(su2double val_loadIncrement) {  }

  /*!
   * \brief A virtual member.
   * \param[in] Value of the load increment for nonlinear structural analysis
   */
  inline virtual su2double GetLoad_Increment() const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   */
  inline virtual unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                                      CConfig *config,
                                                      bool Output) { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] Value of freestream pressure.
   */
  inline virtual void SetPressure_Inf(su2double p_inf){}

  /*!
   * \brief A virtual member.
   * \param[in] Value of freestream temperature.
   */
  inline virtual void SetTemperature_Inf(su2double t_inf){}

  /*!
   * \brief A virtual member.
   * \param[in] Value of freestream density.
   */
  inline virtual void SetDensity_Inf(su2double rho_inf){}

  /*!
   * \brief A virtual member.
   * \param[in] val_dim - Index of the velocity vector.
   * \param[in] val_velocity - Value of the velocity.
   */
  inline virtual void SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) { }

  /*!
   * \brief A virtual member.
   * \param[in] kind_recording - Kind of AD recording.
   */
  inline virtual void SetRecording(CGeometry *geometry, CConfig *config){}

  /*!
   * \brief A virtual member.
   * \param[in] kind_recording - Kind of AD recording.
   */
  inline virtual void SetMesh_Recording(CGeometry **geometry,
                                        CVolumetricMovement *grid_movement,
                                        CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reset - If true reset variables to their initial values.
   */
  inline virtual void RegisterVariables(CGeometry *geometry,
                                        CConfig *config,
                                        bool reset = false) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetFreeStream_Solution(CConfig *config) { }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double* GetVecSolDOFs(void) { return NULL; }

  /*!
   * \brief A virtual member.
   */
  inline virtual unsigned long GetnDOFsGlobal(void) const { return 0; }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetTauWall_WF(CGeometry *geometry,
                                    CSolver** solver_container,
                                    CConfig* config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetNuTilde_WF(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config,
                                    unsigned short val_marker) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void InitTurboContainers(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the average is evaluated.
   */
  inline virtual void PreprocessAverage(CSolver **solver,
                                        CGeometry *geometry,
                                        CConfig *config,
                                        unsigned short marker_flag) { }

  /*!
   * \brief virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the average is evaluated.
   */
  inline virtual void TurboAverageProcess(CSolver **solver,
                                          CGeometry *geometry,
                                          CConfig *config,
                                          unsigned short marker_flag) { }

  /*!
   * \brief virtual member.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  inline virtual void GatherInOutAverageValues(CConfig *config, CGeometry *geometry) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Density on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetAverageDensity(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Pressure on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetAveragePressure(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  inline virtual su2double* GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan) const { return nullptr; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Nu on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetAverageNu(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Kine on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetAverageKine(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Omega on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetAverageOmega(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Nu on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetExtAverageNu(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Kine on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetExtAverageKine(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Omega on the surface <i>val_marker</i>.
   */
  inline virtual su2double GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Density on the surface <i>val_marker</i>.
   */
  inline virtual void SetExtAverageDensity(unsigned short valMarker,
                                           unsigned short valSpan,
                                           su2double valDensity) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Pressure on the surface <i>val_marker</i>.
   */
  inline virtual void SetExtAveragePressure(unsigned short valMarker,
                                            unsigned short valSpan,
                                            su2double valPressure) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  inline virtual void SetExtAverageTurboVelocity(unsigned short valMarker,
                                                 unsigned short valSpan,
                                                 unsigned short valIndex,
                                                 su2double valTurboVelocity) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Nu on the surface <i>val_marker</i>.
   */
  inline virtual void SetExtAverageNu(unsigned short valMarker,
                                      unsigned short valSpan,
                                      su2double valNu) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Kine on the surface <i>val_marker</i>.
   */
  inline virtual void SetExtAverageKine(unsigned short valMarker,
                                        unsigned short valSpan,
                                        su2double valKine) { }

  /*!
   * \brief A virtual member.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Omega on the surface <i>val_marker</i>.
   */
  inline virtual void SetExtAverageOmega(unsigned short valMarker,
                                         unsigned short valSpan,
                                         su2double valOmega) { }

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of inlet pressure.
   */
  inline virtual su2double GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet normal velocity.
   */
  inline virtual su2double* GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan) const {return NULL;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet density.
   */
  inline virtual su2double GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet pressure.
   */
  inline virtual su2double GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet normal velocity.
   */
  inline virtual su2double* GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan) const {return NULL;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetKineIn(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetNuIn(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetKineOut(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline virtual su2double GetNuOut(unsigned short inMarkerTP, unsigned short valSpan) const {return 0;}

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetDensityIn(su2double value,
                                   unsigned short inMarkerTP,
                                   unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetPressureIn(su2double value,
                                    unsigned short inMarkerTP,
                                    unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetTurboVelocityIn(su2double *value,
                                         unsigned short inMarkerTP,
                                         unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetDensityOut(su2double value,
                                    unsigned short inMarkerTP,
                                    unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetPressureOut(su2double value,
                                     unsigned short inMarkerTP,
                                     unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetTurboVelocityOut(su2double *value,
                                          unsigned short inMarkerTP,
                                          unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetKineIn(su2double value,
                                unsigned short inMarkerTP,
                                unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetOmegaIn(su2double value,
                                 unsigned short inMarkerTP,
                                 unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetNuIn(su2double value,
                              unsigned short inMarkerTP,
                              unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetKineOut(su2double value,
                                 unsigned short inMarkerTP,
                                 unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetOmegaOut(su2double value,
                                  unsigned short inMarkerTP,
                                  unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline virtual void SetNuOut(su2double value,
                               unsigned short inMarkerTP,
                               unsigned short valSpan) { }

  /*!
   * \brief A virtual member.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetFreeStream_TurboSolution(CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   */
  inline virtual void SetBeta_Parameter(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CConfig *config,
                                        unsigned short iMesh) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetRoe_Dissipation(CGeometry *geometry, CConfig *config) {}

  /*!
   * \brief A virtual member.
   * \param[in] solver - Solver container
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetDES_LengthScale(CSolver** solver,
                                         CGeometry *geometry,
                                         CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   * \param[in] referenceCoord - Determine if the mesh is deformed from the reference or from the current coordinates.
   */
  inline virtual void DeformMesh(CGeometry **geometry,
                                 CNumerics **numerics,
                                 CConfig *config) { }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   * \param[in] referenceCoord - Determine if the mesh is deformed from the reference or from the current coordinates.
   */
  inline virtual void SetMesh_Stiffness(CGeometry **geometry,
                                        CNumerics **numerics,
                                        CConfig *config) { }

  /*!
   * \brief Routine that sets the flag controlling implicit treatment for periodic BCs.
   * \param[in] val_implicit_periodic - Flag controlling implicit treatment for periodic BCs.
   */
  inline void SetImplicitPeriodic(bool val_implicit_periodic) { implicit_periodic = val_implicit_periodic; }

  /*!
   * \brief Routine that sets the flag controlling solution rotation for periodic BCs.
   * \param[in] val_implicit_periodic - Flag controlling solution rotation for periodic BCs.
   */
  inline void SetRotatePeriodic(bool val_rotate_periodic) { rotate_periodic = val_rotate_periodic; }

  /*!
   * \brief Retrieve the solver name for output purposes.
   * \param[out] val_solvername - Name of the solver.
   */
  inline string GetSolverName(void) {return SolverName;}

  /*!
   * \brief Get the solution fields.
   * \return A vector containing the solution fields.
   */
  inline vector<string> GetSolutionFields(){return fields;}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  inline virtual void ComputeVerificationError(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief Initialize the vertex traction containers at the vertices.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  inline void InitVertexTractionContainer(CGeometry *geometry, CConfig *config){

    unsigned long iVertex;
    unsigned short iMarker;

    VertexTraction = new su2double** [nMarker];
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      VertexTraction[iMarker] = new su2double* [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        VertexTraction[iMarker][iVertex] = new su2double [nDim]();
      }
    }
  }

  /*!
   * \brief Initialize the adjoint vertex traction containers at the vertices.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  inline void InitVertexTractionAdjointContainer(CGeometry *geometry, CConfig *config){

    unsigned long iVertex;
    unsigned short iMarker;

    VertexTractionAdjoint = new su2double** [nMarker];
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      VertexTractionAdjoint[iMarker] = new su2double* [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        VertexTractionAdjoint[iMarker][iVertex] = new su2double [nDim]();
      }
    }
  }

  /*!
   * \brief Compute the tractions at the vertices.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVertexTractions(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set the adjoints of the vertex tractions.
   * \param[in] iMarker  - Index of the marker
   * \param[in] iVertex  - Index of the relevant vertex
   * \param[in] iDim     - Dimension
   */
  inline su2double GetVertexTractions(unsigned short iMarker, unsigned long iVertex, unsigned short iDim) const {
    return VertexTraction[iMarker][iVertex][iDim];
  }

  /*!
   * \brief Register the vertex tractions as output.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void RegisterVertexTractions(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the adjoints of the vertex tractions.
   * \param[in] iMarker  - Index of the marker
   * \param[in] iVertex  - Index of the relevant vertex
   * \param[in] iDim     - Dimension
   * \param[in] val_adjoint - Value received for the adjoint (from another solver)
   */
  inline void StoreVertexTractionsAdjoint(unsigned short iMarker,
                                          unsigned long iVertex,
                                          unsigned short iDim,
                                          su2double val_adjoint){
    VertexTractionAdjoint[iMarker][iVertex][iDim] = val_adjoint;
  }

  /*!
   * \brief Set the adjoints of the vertex tractions to the AD structure.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void SetVertexTractionsAdjoint(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Get minimun volume in the mesh
   * \return
   */
  virtual su2double GetMinimum_Volume() const { return 0.0; }

  /*!
   * \brief Get maximum volume in the mesh
   * \return
   */
  virtual su2double GetMaximum_Volume() const { return 0.0; }

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
