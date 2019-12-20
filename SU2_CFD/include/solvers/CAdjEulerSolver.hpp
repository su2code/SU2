/*!
 * \file CAdjEulerSolver.hpp
 * \brief Headers of the CAdjEulerSolver class
 * \author F. Palacios
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
#include "../variables/CAdjEulerVariable.hpp"

/*!
 * \class CAdjEulerSolver
 * \brief Main class for defining the Euler's adjoint flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios
 */
class CAdjEulerSolver : public CSolver {
protected:
  su2double
  PsiRho_Inf,     /*!< \brief PsiRho variable at the infinity. */
  PsiE_Inf,       /*!< \brief PsiE variable at the infinity. */
  *Phi_Inf;       /*!< \brief Phi vector at the infinity. */
  su2double
  *Sens_Mach,     /*!< \brief Mach sensitivity coefficient for each boundary. */
  *Sens_AoA,      /*!< \brief Angle of attack sensitivity coefficient for each boundary. */
  *Sens_Geo,      /*!< \brief Shape sensitivity coefficient for each boundary. */
  *Sens_Press,         /*!< \brief Pressure sensitivity coefficient for each boundary. */
  *Sens_Temp,          /*!< \brief Temperature sensitivity coefficient for each boundary. */
  *Sens_BPress,        /*!< \brief Back pressure sensitivity coefficient for each boundary. */
  **CSensitivity,      /*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
  ***DonorAdjVar;      /*!< \brief Value of the donor variables at each boundary. */
  su2double Total_Sens_Mach;     /*!< \brief Total mach sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_AoA;      /*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Geo;      /*!< \brief Total shape sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
  su2double Total_Sens_Temp;     /*!< \brief Total farfield sensitivity to temperature. */
  su2double Total_Sens_BPress;   /*!< \brief Total sensitivity to back pressure. */
  bool space_centered;           /*!< \brief True if space centered scheeme used. */
  su2double **Jacobian_Axisymmetric; /*!< \brief Storage for axisymmetric Jacobian. */
  su2double Gamma;                   /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;         /*!< \brief Fluids's Gamma - 1.0  . */
  su2double *FlowPrimVar_i,          /*!< \brief Store the flow solution at point i. */
  *FlowPrimVar_j;                    /*!< \brief Store the flow solution at point j. */
  unsigned long **DonorGlobalIndex;  /*!< \brief Value of the donor global index. */

  su2double pnorm,
  Area_Monitored; /*!< \brief Store the total area of the monitored outflow surface (used for normalization in continuous adjoint outflow conditions) */
  
  su2double ACoeff, ACoeff_inc, ACoeff_old;
  bool Update_ACoeff;

  CAdjEulerVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

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
  void BC_Euler_Wall(CGeometry      *geometry,
                     CSolver        **solver_container,
                     CNumerics      *conv_numerics,
                     CNumerics      *visc_numerics,
                     CConfig        *config,
                     unsigned short val_marker) override;

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
  void BC_Sym_Plane(CGeometry      *geometry,
                    CSolver        **solver_container,
                    CNumerics      *conv_numerics,
                    CNumerics      *visc_numerics,
                    CConfig        *config,
                    unsigned short val_marker) override;
  
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
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter);

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