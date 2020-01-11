/*!
 * \file CFEASolver.hpp
 * \brief Headers of the CFEASolver class
 * \author R. Sanchez.
 * \version 7.0.0 "Blackbird"
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

#include "CSolver.hpp"
#include "../variables/CFEABoundVariable.hpp"

/*! \class CFEASolver
 *  \brief Main class for defining a FEM solver for elastic structural problems.
 *  \author R. Sanchez.
 *  \date July 10, 2015.
 */
class CFEASolver : public CSolver {
private:
  
  su2double  Total_CFEA;        /*!< \brief Total FEA coefficient for all the boundaries. */
  /*!< We maintain the name to avoid defining a new function... */
  
  int nFEA_Terms; 
  bool topol_filter_applied;    /*!< \brief True if density filtering has been performed. */

  su2double *GradN_X,
  *GradN_x;
  
  su2double **Jacobian_c_ij;      /*!< \brief Submatrix to store the constitutive term for node ij. */
  su2double **Jacobian_s_ij;      /*!< \brief Submatrix to store the stress contribution of node ij (diagonal). */
  su2double **MassMatrix_ij;      /*!< \brief Submatrix to store the term ij of the mass matrix. */

  su2double *Res_Ext_Surf;      /*!< \brief Auxiliary vector to store the surface load contribution to the residual */
  su2double *Res_Time_Cont;     /*!< \brief Auxiliary vector to store the surface load contribution to the residual */
  su2double *Res_FSI_Cont;      /*!< \brief Auxiliary vector to store the surface load contribution to the residual */
  
  su2double *Res_Dead_Load;     /*!< \brief Auxiliary vector to store the body load contribution to the residual */
  
  su2double *solutionPredictor;  /*!< \brief Auxiliary vector to store the solution predictor */
  
  su2double *Solution_Interm;    /*!< \brief Auxiliary vector to store the intermediate solution */
  
  su2double *SolRest;      /*!< \brief Auxiliary vector to restart the solution */
  
  su2double *nodeReactions;      /*!< \brief Auxiliary vector to store the reactions */
  
  su2double *normalVertex;       /*!< \brief Auxiliary vector to store the normals to a certain vertex */
  su2double **stressTensor;      /*!< \brief Auxiliary matrix to rebuild the stress tensor and compute reactions */
  
  unsigned long *elProperties;   /*!< \brief Auxiliary vector to read the element properties from file */

  unsigned short *iElem_iDe;	 /*!< \brief For DE cases, ID of the region considered for each iElem. */
  
  su2double a_dt[9];             /*!< \brief Integration constants. */
  
  su2double Conv_Ref[3];        /*!< \brief Reference values for convergence check: DTOL, RTOL, ETOL */
  su2double Conv_Check[3];      /*!< \brief Current values for convergence check: DTOL, RTOL, ETOL */
  su2double FSI_Conv[2];        /*!< \brief Values to check the convergence of the FSI problem. */
  
  su2double loadIncrement;      /*!< \brief Coefficient that determines the amount of load which is applied */
  
  su2double WAitken_Dyn;        /*!< \brief Aitken's dynamic coefficient */
  su2double WAitken_Dyn_tn1;    /*!< \brief Aitken's dynamic coefficient in the previous iteration */
  
  su2double PenaltyValue;       /*!< \brief Penalty value to maintain total stiffness constant */

  su2double Total_OFRefGeom;        /*!< \brief Total Objective Function: Reference Geometry. */
  su2double Total_OFRefNode;        /*!< \brief Total Objective Function: Reference Node. */
  su2double Total_OFVolFrac;        /*!< \brief Total Objective Function: Volume fraction (topology optimization). */
  su2double Total_OFCompliance;     /*!< \brief Total Objective Function: Compliance (topology optimization). */

  su2double Global_OFRefGeom;       /*!< \brief Global Objective Function (added over time steps): Reference Geometry. */
  su2double Global_OFRefNode;       /*!< \brief Global Objective Function (added over time steps): Reference Node. */

  su2double Total_ForwardGradient;  /*!< \brief Vector of the total forward gradient. */
  
  su2double ForceCoeff;             /*!< \brief Load transfer coefficient . */
  su2double RelaxCoeff;             /*!< \brief Relaxation coefficient . */
  su2double FSI_Residual;           /*!< \brief FSI residual. */

protected:

  bool element_based;            /*!< \brief Bool to determine if an element-based file is used. */

  unsigned long nElement;       /*!< \brief Number of elements. */
  unsigned long IterLinSol;     /*!< \brief Number of iterations of the linear solver. */

  su2double **mZeros_Aux;       /*!< \brief Submatrix to make zeros and impose clamped boundary conditions. */
  su2double **mId_Aux;          /*!< \brief Diagonal submatrix to impose clamped boundary conditions. */

  su2double *Res_Stress_i;      /*!< \brief Submatrix to store the nodal stress contribution of node i. */

  CVariable* nodes = nullptr;   /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

public:
  
  CSysVector<su2double> TimeRes_Aux;    /*!< \brief Auxiliary vector for adding mass and damping contributions to the residual. */
  CSysVector<su2double> TimeRes;        /*!< \brief Vector for adding mass and damping contributions to the residual */
  CSysVector<su2double> LinSysReact;    /*!< \brief Vector to store the residual before applying the BCs */

  CSysVector<su2double> LinSysSol_Adj;  /*!< \brief Vector to store the solution of the adjoint problem */
  CSysVector<su2double> LinSysRes_Adj;  /*!< \brief Vector to store the residual of the adjoint problem */

  CSysMatrix<su2double> MassMatrix;     /*!< \brief Sparse structure for storing the mass matrix. */

  CElement*** element_container;   /*!< \brief Vector which the define the finite element structure for each problem. */
  CProperty** element_properties;  /*!< \brief Vector which stores the properties of each element */

  
  /*!
   * \brief Constructor of the class.
   */
  CFEASolver(bool mesh_deform_mode = false);
  
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
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter);
  
  /*!
   * \brief Reset the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter);
  
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
   * \brief Get the value of the reference coordinate to set on the element structure.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] indexNode - Index of the node.
   * \param[in] iDim - Dimension required.
   */
  inline su2double Get_ValCoord(CGeometry *geometry, unsigned long indexNode, unsigned short iDim) {return geometry->node[indexNode]->GetCoord(iDim);}
  
  /*!
   * \brief Compute the stiffness matrix of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the stiffness matrix of the problem and the nodal stress terms at the same time (more efficient if full Newton Raphson).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the mass matrix of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_MassMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the mass residual of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_MassRes(CGeometry *geometry, CNumerics **numerics, CConfig *config);

  /*!
   * \brief Compute the nodal stress terms and add them to the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the stress at the nodes for output purposes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  
  void Compute_NodalStress(CGeometry *geometry, CNumerics **numerics, CConfig *config);
  
  /*!
   * \brief Compute the dead loads.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_DeadLoad(CGeometry *geometry, CNumerics **numerics, CConfig *config);
  
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
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Clamped(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Enforce the solution to be 0 in the clamped nodes - Avoids accumulation of numerical error.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Clamped_Post(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  
  void BC_DispDir(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose a displacement (constraint) boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Normal_Displacement(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  
  /*!
   * \brief Impose a load boundary condition normal to the boundary.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Normal_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose a load boundary condition in cartesian coordinates.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Dir_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose a sine-wave load boundary condition in cartesian coordinates.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Sine_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);
  
  /*!
   * \brief Impose a damping load.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Damper(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Deformable boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */

  void BC_Deforming(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker);

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
  void Solve_System(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Get the residual for FEM structural analysis.
   * \param[in] val_var - Index of the variable.
   * \return Value of the residual for the variable in the position <i>val_var</i>.
   */
  inline su2double GetRes_FEM(unsigned short val_var) { return Conv_Check[val_var]; }
  
  /*!
   * \brief Provide the maximum Von Mises Stress for structural analysis.
   * \return Value of the maximum Von Mises Stress.
   */
  inline su2double GetTotal_CFEA() { return Total_CFEA; }
  
  /*!
   * \brief Retrieve the value of the objective function for a reference geometry
   * \param[out] OFRefGeom - value of the objective function.
   */
  inline su2double GetTotal_OFRefGeom(void){ return Total_OFRefGeom; }
  
  /*!
   * \brief Retrieve the value of the objective function for a reference node
   * \param[out] OFRefNode - value of the objective function.
   */
  inline su2double GetTotal_OFRefNode(void){ return Total_OFRefNode; }
  
  /*!
   * \brief Retrieve the value of the volume fraction objective function
   * \param[out] OFVolFrac - value of the objective function.
   */
  inline su2double GetTotal_OFVolFrac(void){ return Total_OFVolFrac; }
  
  /*!
   * \brief Retrieve the value of the structural compliance objective function
   * \return Value of the objective function.
   */
  inline su2double GetTotal_OFCompliance(void){ return Total_OFCompliance; }

  /*!
   * \brief Determines whether there is an element-based file or not.
   * \return Bool that defines whether the solution has an element-based file or not
   */
  inline bool IsElementBased(void){ return element_based; }

  /*!
   * \brief Set the value of the FEA coefficient.
   * \param[in] val_cfea - Value of the FEA coefficient.
   */
  inline void SetTotal_CFEA(su2double val_cfea) { Total_CFEA = val_cfea; }
  
  /*!
   * \brief Set the value of the objective function for a reference geometry.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference geometry.
   */
  inline void SetTotal_OFRefGeom(su2double val_ofrefgeom) { Total_OFRefGeom = val_ofrefgeom; }
  
  /*!
   * \brief Set the value of the objective function for a reference node.
   * \param[in] val_ofrefnode - Value of the objective function for a reference node.
   */
  inline void SetTotal_OFRefNode(su2double val_ofrefnode) { Total_OFRefNode = val_ofrefnode; }

  /*!
   * \brief Set the value of the force coefficient history for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_forcecoeff_history - Value of the force coefficient.
   */
  inline void SetForceCoeff(su2double val_forcecoeff_history) { ForceCoeff = val_forcecoeff_history; }

  /*!
   * \brief Set the value of the FSI residual for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  inline void SetFSI_Residual(su2double val_FSI_residual) { FSI_Residual = val_FSI_residual; }

  /*!
   * \brief Set the value of the FSI residual for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  inline void SetRelaxCoeff(su2double val_relaxcoeff_history){ RelaxCoeff = val_relaxcoeff_history;}

  /*!
   * \brief Get the value of the force coefficient history for the history file.
   * \param[out] val_forcecoeff_history - Value of the force coefficient.
   */
  inline su2double GetForceCoeff(void) { return ForceCoeff; }

  /*!
   * \brief Get the value of the relaxation coefficient history for the history file.
   * \param[out] val_relaxcoeff_history - Value of the relaxation coefficient.
   */
  inline su2double GetRelaxCoeff(void) { return RelaxCoeff; }

  /*!
   * \brief Get the value of the FSI residual for the history file.
   * \param[out] val_FSI_residual - Value of the residual.
   */
  inline su2double GetFSI_Residual(void) { return FSI_Residual; }
  
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
   * \brief Compute the compliance objective function
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFCompliance(CGeometry *geometry, CSolver **solver_container, CConfig *config);

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
  inline void SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { FSI_Conv[val_index] = val_criteria; }
  
  /*!
   * \brief Get the value of the FSI convergence.
   * \param[in]  Value of interest: 0 - Initial value, 1 - Current value.
   * \return Values to compare
   */
  inline su2double GetFSI_ConvValue(unsigned short val_index){ return FSI_Conv[val_index]; }
  
  /*!
   * \brief Retrieve the value of the dynamic Aitken relaxation factor.
   * \return Value of the dynamic Aitken relaxation factor.
   */
  inline su2double GetWAitken_Dyn(void) { return WAitken_Dyn; }
  
  /*!
   * \brief Retrieve the value of the last Aitken relaxation factor in the previous time step.
   * \return Value of the last Aitken relaxation factor in the previous time step.
   */
  inline su2double GetWAitken_Dyn_tn1(void) { return WAitken_Dyn_tn1; }
  
  /*!
   * \brief Set the value of the dynamic Aitken relaxation factor
   * \param[in] Value of the dynamic Aitken relaxation factor
   */
  inline void SetWAitken_Dyn(su2double waitk) { WAitken_Dyn = waitk; }
  
  /*!
   * \brief Set the value of the last Aitken relaxation factor in the current time step.
   * \param[in] Value of the last Aitken relaxation factor in the current time step.
   */
  inline void SetWAitken_Dyn_tn1(su2double waitk_tn1) { WAitken_Dyn_tn1 = waitk_tn1; }
  
  /*!
   * \brief Set the value of the load increment for nonlinear structural analysis
   * \param[in] Value of the coefficient
   */
  inline void SetLoad_Increment(su2double val_loadIncrement) { loadIncrement = val_loadIncrement; }
  
  /*!
   * \brief Get the value of the load increment for nonlinear structural analysis
   * \param[in] Value of the coefficient
   */
  inline su2double GetLoad_Increment(void) { return loadIncrement; }

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
  inline unsigned short Get_iElem_iDe(unsigned long iElem){ return iElem_iDe[iElem]; }
  
  /*!
   * \brief Retrieve the Mass Matrix term (to add to the Jacobian of the adjoint problem)
   * \param[in] iPoint - Point i of the Mass Matrix.
   * \param[in] jPoint - Point j of the Mass Matrix.
   * \param[in] iVar - Variable i of the Mass Matrix submatrix.
   * \param[in] iVar - Variable j of the Mass Matrix submatrix.
   */
  inline su2double Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar){
    return MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar); }
  
  
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
   * \param[in] reset - Not used by this class ATM.
   */
  void RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Filter the density field for topology optimization applications
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void FilterElementDensities(CGeometry *geometry, CConfig *config);

};

#include "../solver_inlines/CFEASolver.inl"