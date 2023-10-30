/*!
 * \file CFEASolver.hpp
 * \brief Finite element solver for elasticity problems.
 * \author R. Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFEASolverBase.hpp"

/*!
 * \class CFEASolver
 * \ingroup Elasticity_Equations
 * \brief Main class for defining a FEM solver for elastic structural problems.
 * \author R. Sanchez.
 */
class CFEASolver : public CFEASolverBase {
protected:

  unsigned long omp_chunk_size;     /*!< \brief Chunk size used in light point loops. */

  su2double Total_CFEA;             /*!< \brief Total FEA coefficient for all the boundaries. */

  unsigned short *iElem_iDe = nullptr; /*!< \brief For DE cases, ID of the region considered for each iElem. */

  su2double a_dt[9];                /*!< \brief Time integration constants. */

  su2double Conv_Check[3];          /*!< \brief Current values for convergence check: UTOL, RTOL, ETOL. */
  su2double FSI_Conv[2];            /*!< \brief Values to check the convergence of the FSI problem. */

  unsigned long idxIncrement = 0;   /*!< \brief Index of the current load increment */
  su2double loadIncrement = 1.0;    /*!< \brief Coefficient that determines the amount of load which is applied. */

  su2double WAitken_Dyn;            /*!< \brief Aitken's dynamic coefficient. */
  su2double WAitken_Dyn_tn1;        /*!< \brief Aitken's dynamic coefficient in the previous iteration. */

  su2double PenaltyValue;           /*!< \brief Penalty value to maintain total stiffness constant. */

  su2double Total_OFRefGeom;        /*!< \brief Total Objective Function: Reference Geometry. */
  su2double Total_OFRefNode;        /*!< \brief Total Objective Function: Reference Node. */
  su2double Total_OFVolFrac;        /*!< \brief Total Objective Function: Volume fraction (topology optimization). */
  su2double Total_OFDiscreteness;   /*!< \brief Total Objective Function: Discreteness (topology optimization). */
  su2double Total_OFCompliance;     /*!< \brief Total Objective Function: Compliance (topology optimization). */
  su2double Total_OFStressPenalty;  /*!< \brief Total Objective Function: Stress penalty. */

  su2double Global_OFRefGeom;       /*!< \brief Global Objective Function (added over time steps): Reference Geometry. */
  su2double Global_OFRefNode;       /*!< \brief Global Objective Function (added over time steps): Reference Node. */

  su2double Total_ForwardGradient;  /*!< \brief Vector of the total forward gradient. */

  su2double ForceCoeff;             /*!< \brief Load transfer coefficient . */
  su2double RelaxCoeff;             /*!< \brief Relaxation coefficient . */
  su2double FSI_Residual;           /*!< \brief FSI residual. */

  CSysVector<su2double> TimeRes_Aux;  /*!< \brief Auxiliary vector for adding mass and damping contributions to the residual. */
  CSysVector<su2double> TimeRes;      /*!< \brief Vector for adding mass and damping contributions to the residual */
  CSysVector<su2double> LinSysReact;  /*!< \brief Vector to store the residual before applying the BCs */

#ifndef CODI_FORWARD_TYPE
  CSysMatrix<su2mixedfloat> MassMatrix;   /*!< \brief Sparse structure for storing the mass matrix. */
#else
  CSysMatrix<su2double> MassMatrix;
#endif

  CProperty** element_properties = nullptr; /*!< \brief Vector which stores the properties of each element */

#ifdef HAVE_OMP
  vector<GridColor<> > ElemColoring;   /*!< \brief Element colors. */
  bool LockStrategy = false;           /*!< \brief Whether to use an OpenMP lock to guard updates of the Jacobian. */
  vector<omp_lock_t> UpdateLocks;      /*!< \brief Locks that may be used to protect accesses to CSysMatrix/Vector in element loops. */
#else
  array<DummyGridColor<>,1> ElemColoring;      /*--- Behaves like a normal integer type. ---*/
  static constexpr bool LockStrategy = false;  /*--- Lock strategy is never needed for MPI-only. ---*/
  DummyVectorOfLocks UpdateLocks;
#endif

  bool element_based;          /*!< \brief Bool to determine if an element-based file is used. */
  bool topol_filter_applied;   /*!< \brief True if density filtering has been performed. */
  bool initial_calc = true;    /*!< \brief Becomes false after first call to Preprocessing. */

  /*!
   * \brief The highest level in the variable hierarchy this solver can safely use,
   * CVariable is the common denominator between the FEA and Mesh deformation variables.
   */
  CVariable* nodes = nullptr;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

  /*!
   * \brief Get the element container index and number of nodes of a given VTK type.
   * \param[in] VTK_Type - Type of element.
   * \param[out] EL_KIND - Element container index.
   * \param[out] nNodes - Number of nodes.
   */
  static void GetElemKindAndNumNodes(unsigned short VTK_Type, int& EL_KIND, unsigned short& nNodes) {
    switch (VTK_Type) {
      case TRIANGLE:      nNodes = 3; EL_KIND = EL_TRIA;  break;
      case QUADRILATERAL: nNodes = 4; EL_KIND = EL_QUAD;  break;
      case TETRAHEDRON:   nNodes = 4; EL_KIND = EL_TETRA; break;
      case PYRAMID:       nNodes = 5; EL_KIND = EL_PYRAM; break;
      case PRISM:         nNodes = 6; EL_KIND = EL_PRISM; break;
      case HEXAHEDRON:    nNodes = 8; EL_KIND = EL_HEXA;  break;
      default: assert(false); nNodes = 0; EL_KIND = -(1<<30); break;
    }
  }

  /*!
   * \brief Actions required to initialize the supporting variables for hybrid parallel execution.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void HybridParallelInitialization(CGeometry* geometry);

  /*!
   * \brief Set container of element properties.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_ElementProperties(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set a reference geometry for .
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_ReferenceGeometry(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Set a reference geometry for prestretched conditions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_Prestretch(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Mitigation for an issue with Dirichlet boundary conditions and MPI,
   * some ranks do not get enough of the markers to cover their halo points.
   * This breaks the symmetry of the global matrix as columns are not fully eliminated.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] markers - List of essential BC markers.
   */
  void Set_VertexEliminationSchedule(CGeometry *geometry, const vector<unsigned short>& markers);

  /*!
   * \brief Compute constants for time integration.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_IntegrationConstants(const CConfig *config);

  /*!
   * \brief Write the forward mode gradient to file.
   * \param[in] config - Definition of the particular problem.
   * \param[in] newFile - If true, start file from scratch, else append values.
   * \param[in] fun - Value of the objective function.
   * \param[in] fun_avg - Time-averaged value of the objective function.
   * \param[in] der - Value of the derivative.
   * \param[in] der_avg - Time-averaged value of the derivative.
   */
  void OutputForwardModeGradient(const CConfig *config,
                                 bool newFile,
                                 su2double fun,
                                 su2double fun_avg,
                                 su2double der,
                                 su2double der_avg) const;

  /*!
   * \brief Compute the objective function for a reference geometry
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFRefGeom(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the objective function for a reference node
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFRefNode(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the objective function for a volume fraction
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFVolFrac(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the compliance objective function
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_OFCompliance(CGeometry *geometry, const CConfig *config);

public:
  /*!
   * \brief Constructor of the class.
   */
  CFEASolver(LINEAR_SOLVER_MODE mesh_deform_mode = LINEAR_SOLVER_MODE::STANDARD);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEASolver(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEASolver(void) override;

  /*!
   * \brief Set residuals to zero.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Preprocessing(CGeometry *geometry,
                     CSolver **solver_container,
                     CConfig *config,
                     CNumerics **numerics,
                     unsigned short iMesh,
                     unsigned long Iteration,
                     unsigned short RunTime_EqSystem,
                     bool Output) override;

  /*!
   * \brief Set the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) override;

  /*!
   * \brief Get the value of the reference coordinate to set on the element structure.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] indexNode - Index of the node.
   * \param[in] iDim - Dimension required.
   */
  inline virtual su2double Get_ValCoord(const CGeometry *geometry,
                                        unsigned long indexNode,
                                        unsigned short iDim) const {
    return geometry->nodes->GetCoord(indexNode, iDim);
  }

  /*!
   * \brief Compute the stiffness matrix of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_StiffMatrix(CGeometry *geometry,
                           CNumerics **numerics,
                           const CConfig *config) final;

  /*!
   * \brief Compute the stiffness matrix of the problem and the nodal stress terms at the same time.
   * \note This is more efficient for full Newton Raphson.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_StiffMatrix_NodalStressRes(CGeometry *geometry,
                                          CNumerics **numerics,
                                          const CConfig *config) final;

  /*!
   * \brief Compute the mass matrix of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_MassMatrix(const CGeometry *geometry,
                          CNumerics **numerics,
                          const CConfig *config) final;

  /*!
   * \brief Compute the mass residual of the problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_MassRes(const CGeometry *geometry,
                       CNumerics **numerics,
                       const CConfig *config) final;

  /*!
   * \brief Compute the nodal stress terms and add them to the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_NodalStressRes(CGeometry *geometry,
                              CNumerics **numerics,
                              const CConfig *config) final;

  /*!
   * \brief Compute the stress at the nodes for output purposes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_NodalStress(CGeometry *geometry,
                           CNumerics **numerics,
                           const CConfig *config);

  /*!
   * \brief Compute the dead loads.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Compute_DeadLoad(CGeometry *geometry,
                        CNumerics **numerics,
                        const CConfig *config) final;

  /*!
   * \brief Clamped boundary conditions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Clamped(CGeometry *geometry,
                  const CConfig *config,
                  unsigned short val_marker) final;

  /*!
   * \brief Enforce the solution to be 0 in the clamped nodes - Avoids accumulation of numerical error.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Clamped_Post(CGeometry *geometry,
                       const CConfig *config,
                       unsigned short val_marker) final;

  /*!
   * \brief Symmetry boundary conditions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Sym_Plane(CGeometry *geometry,
                    const CConfig *config,
                    unsigned short val_marker) final;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_DispDir(CGeometry *geometry,
                  const CConfig *config,
                  unsigned short val_marker) final;

  /*!
   * \brief Impose a load boundary condition normal to the boundary.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Normal_Load(CGeometry *geometry,
                      const CConfig *config,
                      unsigned short val_marker) final;

  /*!
   * \brief Impose a load boundary condition in cartesian coordinates.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Dir_Load(CGeometry *geometry,
                   const CConfig *config,
                   unsigned short val_marker) final;

  /*!
   * \brief Impose a sine-wave load boundary condition in cartesian coordinates.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_Sine_Load(CGeometry *geometry,
                           const CConfig *config,
                           unsigned short val_marker) final { }

  /*!
   * \brief Impose a damping load.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Damper(CGeometry *geometry,
                 const CConfig *config,
                 unsigned short val_marker) final;

  /*!
   * \brief Iterate using an implicit Newmark solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Numerical methods.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitNewmark_Iteration(const CGeometry *geometry,
                                 CNumerics **numerics,
                                 const CConfig *config) final;

  /*!
   * \brief Update the solution using an implicit Newmark solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitNewmark_Update(const CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitNewmark_Relaxation(const CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief Iterate using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Numerical methods.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_Iteration(const CGeometry *geometry,
                                  CNumerics **numerics,
                                  const CConfig *config) final;

  /*!
   * \brief Update the solution using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_UpdateDisp(const CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief Update the solution using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_UpdateSolution(const CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief Update the solution using an implicit Generalized Alpha solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void GeneralizedAlpha_UpdateLoads(const CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief Postprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] numerics - Implementation of numerical method.
   * \param[in] of_comp_mode - Mode to compute just the objective function.
   */
  void Postprocessing(CGeometry *geometry, CConfig *config, CNumerics **numerics, bool of_comp_mode) final;

  /*!
   * \brief Routine to solve the Jacobian-Residual linearized system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Solve_System(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Get the residual for FEM structural analysis.
   * \param[in] val_var - Index of the variable.
   * \return Value of the residual for the variable in the position <i>val_var</i>.
   */
  inline su2double GetRes_FEM(unsigned short val_var) const final { return Conv_Check[val_var]; }

  /*!
   * \brief Provide the maximum Von Mises Stress for structural analysis.
   */
  inline su2double GetTotal_CFEA(void) const final { return Total_CFEA; }

  /*!
   * \brief Retrieve the value of the objective function for a reference geometry
   */
  inline su2double GetTotal_OFRefGeom(void) const final { return Total_OFRefGeom; }

  /*!
   * \brief Retrieve the value of the objective function for a reference node
   */
  inline su2double GetTotal_OFRefNode(void) const final { return Total_OFRefNode; }

  /*!
   * \brief Retrieve the value of the volume fraction objective function
   */
  inline su2double GetTotal_OFVolFrac(void) const final { return Total_OFVolFrac; }

  /*!
   * \brief Retrieve the value of the discreteness objective function
   */
  inline su2double GetTotal_OFDiscreteness(void) const final { return Total_OFDiscreteness; }

  /*!
   * \brief Retrieve the value of the structural compliance objective function
   */
  inline su2double GetTotal_OFCompliance(void) const final { return Total_OFCompliance; }

  /*!
   * \brief Retrieve the value of the stress penalty objective function
   */
  inline su2double GetTotal_OFStressPenalty(void) const final { return Total_OFStressPenalty; }

  /*!
   * \brief Compute the objective function.
   * \param[in] config - Definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void Evaluate_ObjFunc(const CConfig *config, CSolver**) final {
    Total_ComboObj = 0.0;
    switch (config->GetKind_ObjFunc()) {
      case REFERENCE_GEOMETRY:
        Total_ComboObj = GetTotal_OFRefGeom();
        break;
      case REFERENCE_NODE:
        Total_ComboObj = GetTotal_OFRefNode();
        break;
      case TOPOL_COMPLIANCE:
        Total_ComboObj = GetTotal_OFCompliance();
        break;
      case VOLUME_FRACTION:
        Total_ComboObj = GetTotal_OFVolFrac();
        break;
      case TOPOL_DISCRETENESS:
        Total_ComboObj = GetTotal_OFDiscreteness();
        break;
      case STRESS_PENALTY:
        Total_ComboObj = GetTotal_OFStressPenalty();
        break;
      case CUSTOM_OBJFUNC:
        Total_ComboObj = Total_Custom_ObjFunc;
        break;
    }
  }

  /*!
   * \brief Determines whether there is an element-based file or not.
   * \return Bool that defines whether the solution has an element-based file or not
   */
  inline bool IsElementBased(void) const final { return element_based; }

  /*!
   * \brief Set the value of the FEA coefficient.
   * \param[in] val_cfea - Value of the FEA coefficient.
   */
  inline void SetTotal_CFEA(su2double val_cfea) final { Total_CFEA = val_cfea; }

  /*!
   * \brief Set the value of the objective function for a reference geometry.
   * \param[in] val_ofrefgeom - Value of the objective function for a reference geometry.
   */
  inline void SetTotal_OFRefGeom(su2double val_ofrefgeom) final { Total_OFRefGeom = val_ofrefgeom; }

  /*!
   * \brief Set the value of the objective function for a reference node.
   * \param[in] val_ofrefnode - Value of the objective function for a reference node.
   */
  inline void SetTotal_OFRefNode(su2double val_ofrefnode) final { Total_OFRefNode = val_ofrefnode; }

  /*!
   * \brief Set the value of the force coefficient history for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_forcecoeff - Value of the force coefficient.
   */
  inline void SetForceCoeff(su2double val_forcecoeff) final { ForceCoeff = val_forcecoeff; }

  /*!
   * \brief Set the value of the FSI residual for the history file.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_FSI_residual - Value of the residual.
   */
  inline void SetFSI_Residual(su2double val_FSI_residual) final { FSI_Residual = val_FSI_residual; }

  /*!
   * \brief Set the value of the FSI relaxation factor.
   * \param[in] iBGS - Number of BGS iteration.
   * \param[in] val_relaxcoeff - Value of the relaxation factor.
   */
  inline void SetRelaxCoeff(su2double val_relaxcoeff) final { RelaxCoeff = val_relaxcoeff; }

  /*!
   * \brief Get the value of the force coefficient history for the history file.
   * \param[out] val_forcecoeff_history - Value of the force coefficient.
   */
  inline su2double GetForceCoeff(void) const final { return ForceCoeff; }

  /*!
   * \brief Get the value of the relaxation coefficient history for the history file.
   * \param[out] val_relaxcoeff_history - Value of the relaxation coefficient.
   */
  inline su2double GetRelaxCoeff(void) const final { return RelaxCoeff; }

  /*!
   * \brief Get the value of the FSI residual for the history file.
   * \param[out] val_FSI_residual - Value of the residual.
   */
  inline su2double GetFSI_Residual(void) const final { return FSI_Residual; }

  /*!
   * \brief Predictor for structural displacements based on previous iterations
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Configuration of the problem.
   */
  void PredictStruct_Displacement(CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief Computation of Aitken's coefficient.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iOuterIter - Current outer iteration.
   */
  void ComputeAitken_Coefficient(CGeometry *geometry,
                                 const CConfig *config,
                                 unsigned long iOuterIter) final;

  /*!
   * \brief Aitken's relaxation of the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetAitken_Relaxation(CGeometry *geometry, const CConfig *config) final;

  /*!
   * \brief Compute the penalty due to the stiffness increase
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Stiffness_Penalty(CGeometry *geometry,
                         CNumerics **numerics_container,
                         CConfig *config);

  /*!
   * \brief Get the value of the FSI convergence.
   * \param[in] Set value of interest: 0 - Initial value, 1 - Current value.
   */
  inline void SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) final { FSI_Conv[val_index] = val_criteria; }

  /*!
   * \brief Get the value of the FSI convergence.
   * \param[in]  Value of interest: 0 - Initial value, 1 - Current value.
   * \return Values to compare
   */
  inline su2double GetFSI_ConvValue(unsigned short val_index) const final { return FSI_Conv[val_index]; }

  /*!
   * \brief Store the value of the last Aitken relaxation factor in the current time step.
   */
  inline void SetWAitken_Dyn_tn1() final { WAitken_Dyn_tn1 = WAitken_Dyn; }

  /*!
   * \brief Set the value of the load increment for nonlinear structural analysis
   * \param[in] Value of the coefficient
   */
  inline void SetLoad_Increment(unsigned long iInc, su2double loadInc) final {
    idxIncrement = iInc;
    loadIncrement = loadInc;
  }

  /*!
   * \brief Get the value of the load increment for nonlinear structural analysis
   * \param[in] Value of the coefficient
   */
  inline su2double GetLoad_Increment(void) const final { return loadIncrement; }

  /*!
   * \brief Retrieve the iDe index for DE computations
   * \param[in] iElem - element parameter.
   * \param[out] iElem_iDe - ID of the Dielectric Elastomer region.
   */
  inline unsigned short Get_iElem_iDe(unsigned long iElem) const final { return iElem_iDe[iElem]; }

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

  /*!
   * \brief Get multiplier for loads.
   * \param[in] CurrentTime - Current time step.
   * \param[in] config - Definition of the particular problem.
   */
  su2double Compute_LoadCoefficient(su2double CurrentTime,
                                    su2double RampTime,
                                    const CConfig *config) final;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reset - Not used by this class ATM.
   */
  void RegisterVariables(CGeometry *geometry,
                         CConfig *config,
                         bool reset) override;

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
  void FilterElementDensities(CGeometry *geometry, const CConfig *config);

};
