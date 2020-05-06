/*!
 * \file CIncEulerSolver.hpp
 * \brief Headers of the CIncEulerSolver class
 * \author F. Palacios, T. Economon, T. Albring
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include "../variables/CIncEulerVariable.hpp"

/*!
 * \class CIncEulerSolver
 * \brief Main class for defining the incompressible Euler flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncEulerSolver : public CSolver {
protected:

  su2double
  Density_Inf,      /*!< \brief Density at the infinity. */
  Pressure_Inf,     /*!< \brief Pressure at the infinity. */
  *Velocity_Inf,    /*!< \brief Flow Velocity vector at the infinity. */
  Temperature_Inf;  /*!< \brief Temperature at infinity. */

  su2double
  *CD_Inv,       /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Inv,       /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,      /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,      /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,      /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,      /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Inv,     /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Inv,     /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Inv,     /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,      /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,      /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,      /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Inv,    /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Inv,    /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Inv,   /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv,  /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Inv,          /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Inv,        /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Inv,      /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Inv,      /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  *CD_Mnt,      /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Mnt,      /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Mnt,     /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Mnt,     /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Mnt,     /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Mnt,     /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Mnt,    /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Mnt,    /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Mnt,    /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Mnt,     /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Mnt,     /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Mnt,     /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Mnt,    /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Mnt,    /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Mnt,   /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Mnt,  /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Mnt,   /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Mnt,   /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Mnt,   /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Mnt,   /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Mnt,   /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Mnt,   /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Mnt,          /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Mnt,        /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Mnt,            /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Mnt,            /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  **CPressure,        /*!< \brief Pressure coefficient for each boundary and vertex. */
  **CPressureTarget,  /*!< \brief Target Pressure coefficient for each boundary and vertex. */
  **HeatFlux,         /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **HeatFluxTarget,   /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **YPlus,            /*!< \brief Yplus for each boundary and vertex. */
  ***CharacPrimVar,   /*!< \brief Value of the characteristic variables at each boundary. */
  *ForceInviscid,     /*!< \brief Inviscid force for each boundary. */
  *MomentInviscid,    /*!< \brief Inviscid moment for each boundary. */
  *ForceMomentum,     /*!< \brief Inviscid force for each boundary. */
  *MomentMomentum,    /*!< \brief Inviscid moment for each boundary. */
  InverseDesign;      /*!< \brief Inverse design functional for each boundary. */
  su2double
  **Inlet_Ptotal,    /*!< \brief Value of the Total P. */
  **Inlet_Ttotal,    /*!< \brief Value of the Total T. */
  ***Inlet_FlowDir;  /*!< \brief Value of the Flow Direction. */

  su2double
  AllBound_CD_Inv,      /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Inv,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Inv,     /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Inv,     /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Inv,     /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Inv,     /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Inv,     /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Inv,     /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Inv,     /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Inv,    /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Inv,    /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Inv,    /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Inv,    /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Inv,  /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Inv,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Inv;      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */


  su2double
  AllBound_CD_Mnt,      /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Mnt,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Mnt,     /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Mnt,     /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Mnt,     /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Mnt,     /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Mnt,     /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Mnt,     /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Mnt,     /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Mnt,    /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Mnt,    /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Mnt,    /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Mnt,    /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Mnt,  /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Mnt,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Mnt;      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */

  su2double
  Total_ComboObj,        /*!< \brief Total 'combo' objective for all monitored boundaries */
  Total_CD,              /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CL,              /*!< \brief Total lift coefficient for all the boundaries. */
  Total_CSF,             /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CMx,             /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMy,             /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMz,             /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CoPx,            /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CoPy,            /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CoPz,            /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CFx,             /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy,             /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz,             /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CEff,            /*!< \brief Total efficiency coefficient for all the boundaries. */
  Total_CMerit,          /*!< \brief Total rotor Figure of Merit for all the boundaries. */
  Total_CT,              /*!< \brief Total thrust coefficient for all the boundaries. */
  Total_CQ,              /*!< \brief Total torque coefficient for all the boundaries. */
  Total_Heat,            /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat,         /*!< \brief Maximum heat flux on all boundaries. */
  Total_CpDiff,          /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_HeatFluxDiff,    /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_Custom_ObjFunc,  /*!< \brief Total custom objective function for all the boundaries. */
  Total_MassFlowRate;    /*!< \brief Total Mass Flow Rate on monitored boundaries. */
  su2double
  *Surface_CL,          /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CD,          /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSF,         /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CEff,        /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,         /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,         /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,         /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,         /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,         /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz,         /*!< \brief z Moment coefficient for each monitoring surface. */
  *Surface_HF_Visc,     /*!< \brief Total (integrated) heat flux for each monitored surface. */
  *Surface_MaxHF_Visc;  /*!< \brief Maximum heat flux for each monitored surface. */

  su2double *SecondaryVar_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
  *SecondaryVar_j;            /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double *PrimVar_i,       /*!< \brief Auxiliary vector for storing the solution at point i. */
  *PrimVar_j;                 /*!< \brief Auxiliary vector for storing the solution at point j. */
  bool space_centered,  /*!< \brief True if space centered scheeme used. */
  euler_implicit,       /*!< \brief True if euler implicit scheme used. */
  least_squares;        /*!< \brief True if computing gradients by least squares. */
  su2double Gamma;                  /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;        /*!< \brief Fluids's Gamma - 1.0  . */

  su2double *Primitive,  /*!< \brief Auxiliary nPrimVar vector. */
  *Primitive_i,          /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Primitive_j;          /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  CFluidModel  *FluidModel;   /*!< \brief fluid model used in the solver */
  su2double **Preconditioner; /*!< \brief Auxiliary matrix for storing the low speed preconditioner. */

  /* Sliding meshes variables */

  su2double ****SlidingState;
  int **SlidingStateNodes;

  CIncEulerVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

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
  ~CIncEulerSolver(void) override;

  /*!
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh) final;

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline CFluidModel* GetFluidModel(void) const final { return FluidModel;}

  /*!
   * \brief Compute the density at the infinity.
   * \return Value of the density at the infinity.
   */
  inline su2double GetDensity_Inf(void) const final { return Density_Inf; }

  /*!
   * \brief Compute 2-norm of the velocity at the infinity.
   * \return Value of the 2-norm of the velocity at the infinity.
   */
  inline su2double GetModVelocity_Inf(void) const final {
    su2double Vel2 = 0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    return sqrt(Vel2);
  }


  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline su2double GetPressure_Inf(void) const final { return Pressure_Inf; }

    /*!
   * \brief Get the temperature value at infinity.
   * \return Value of the temperature at infinity.
   */
  inline su2double GetTemperature_Inf(void) const { return Temperature_Inf; }

  /*!
   * \brief Compute the density multiply by velocity at the infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \return Value of the density multiply by the velocity at the infinity.
   */
  inline su2double GetDensity_Velocity_Inf(unsigned short val_dim) const final { return Density_Inf*Velocity_Inf[val_dim]; }

  /*!
   * \brief Get the velocity at the infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \return Value of the velocity at the infinity.
   */
  inline su2double GetVelocity_Inf(unsigned short val_dim) const final { return Velocity_Inf[val_dim]; }

  /*!
   * \brief Get the velocity at the infinity.
   * \return Value of the velocity at the infinity.
   */
  inline su2double *GetVelocity_Inf(void) final { return Velocity_Inf; }

  /*!
   * \brief Set the velocity at infinity.
   * \param[in] val_dim - Index of the velocity vector.
   * \param[in] val_velocity - Value of the velocity.
   */
  inline void SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) final {
    Velocity_Inf[val_dim] = val_velocity;
  }

  /*!
   * \brief Compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned long Iteration) override;

  /*!
   * \brief Compute the spatial integration using a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Centered_Residual(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics **numerics_container,
                        CConfig *config,
                        unsigned short iMesh,
                        unsigned short iRKStep) final;

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry *geometry,
                      CSolver **solver_container,
                      CNumerics **numerics_container,
                      CConfig *config,
                      unsigned short iMesh) final;

  /*!
   * \brief Source term integration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) final;

  /*!
   * \brief Source term integration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Template(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) final;

  /*!
   * \brief Compute primitive variables and their gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned short iRKStep,
                    unsigned short RunTime_EqSystem,
                    bool Output) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry *geometry,
                      CSolver **solver_container,
                      CConfig *config,
                      unsigned short iMesh) final;

  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                       CConfig *config,
                                       bool Output);

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
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetPrimitive_Gradient_GG(CGeometry *geometry,
                                CConfig *config,
                                bool reconstruction = false) final;

  /*!
   * \brief Compute the gradient of the primitive variables using a Least-Squares method,
   *        and stores the result in the <i>Gradient_Primitive</i> variable.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetPrimitive_Gradient_LS(CGeometry *geometry,
                                CConfig *config,
                                bool reconstruction = false) final;

  /*!
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) final;

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
  void Evaluate_ObjFunc(CConfig *config) final;

  /*!
   * \author: T. Kattmann
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
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
   * \brief Impose the far-field boundary condition using characteristics.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) final;

  /*!
   * \brief Impose the symmetry boundary condition using the residual.
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
                    unsigned short val_marker) final;

  /*!
   * \brief Impose a subsonic inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) final;

  /*!
   * \brief Impose a custom or verification boundary condition.
   * \param[in] geometry         - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics    - Description of the convective numerical method.
   * \param[in] visc_numerics    - Description of the viscous numerical method.
   * \param[in] config           - Definition of the particular problem.
   * \param[in] val_marker       - Surface marker where the boundary condition is applied.
   */
  void BC_Custom(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) final;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) final;

  /*!
   * \brief Impose the interface state across sliding meshes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
   void BC_Fluid_Interface(CGeometry *geometry,
                           CSolver **solver_container,
                           CNumerics *conv_numerics,
                           CNumerics *visc_numerics,
                           CConfig *config) final;

  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry *geometry,
                   CSolver **solver_container,
                   CNumerics *numerics,
                   CConfig *config) final;

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
  void ExplicitRK_Iteration(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep) final;

  /*!
   * \brief Update the solution using the explicit Euler scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ExplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) final;

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Pressure_Forces(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Compute the pressure forces and all the adimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Momentum_Forces(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Update the solution using an implicit Euler scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) final;

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over a nonlinear iteration for stability.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeUnderRelaxationFactor(CSolver **solver, CConfig *config) final;

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL(unsigned short val_marker) const final { return Surface_CL[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD(unsigned short val_marker) const final { return Surface_CD[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF(unsigned short val_marker) const final { return Surface_CSF[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff(unsigned short val_marker) const final { return Surface_CEff[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx(unsigned short val_marker) const final { return Surface_CFx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy(unsigned short val_marker) const final { return Surface_CFy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz(unsigned short val_marker) const final { return Surface_CFz[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx(unsigned short val_marker) const final { return Surface_CMx[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy(unsigned short val_marker) const final { return Surface_CMy[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz(unsigned short val_marker) const final { return Surface_CMz[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Inv(unsigned short val_marker) const final { return Surface_CL_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Inv(unsigned short val_marker) const final { return Surface_CD_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Inv(unsigned short val_marker) const final { return Surface_CSF_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Inv(unsigned short val_marker) const final { return Surface_CEff_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Inv(unsigned short val_marker) const final { return Surface_CFx_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Inv(unsigned short val_marker) const final { return Surface_CFy_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Inv(unsigned short val_marker) const final { return Surface_CFz_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Inv(unsigned short val_marker) const final { return Surface_CMx_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Inv(unsigned short val_marker) const final { return Surface_CMy_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Inv(unsigned short val_marker) const final { return Surface_CMz_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional lift coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CL_Mnt(unsigned short val_marker) const final { return Surface_CL_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CD_Mnt(unsigned short val_marker) const final { return Surface_CD_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CSF_Mnt(unsigned short val_marker) const final { return Surface_CSF_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional side-force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the side-force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CEff_Mnt(unsigned short val_marker) const final { return Surface_CEff_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional x force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFx_Mnt(unsigned short val_marker) const final { return Surface_CFx_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional y force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFy_Mnt(unsigned short val_marker) const final { return Surface_CFy_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional z force coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z force coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CFz_Mnt(unsigned short val_marker) const final { return Surface_CFz_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional x moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the x moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMx_Mnt(unsigned short val_marker) const final { return Surface_CMx_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional y moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the y moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMy_Mnt(unsigned short val_marker) const final { return Surface_CMy_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional z moment coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the z moment coefficient on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_CMz_Mnt(unsigned short val_marker) const final { return Surface_CMz_Mnt[val_marker]; }

  /*!
   * \brief Provide the non dimensional drag coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the drag coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Inv(unsigned short val_marker) const final { return CD_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional sideforce coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the sideforce coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCSF_Inv(unsigned short val_marker) const final { return CSF_Inv[val_marker]; }

  /*!
   * \brief Provide the non dimensional efficiency coefficient (inviscid contribution).
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the efficiency coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCEff_Inv(unsigned short val_marker) const final { return CEff_Inv[val_marker]; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CSF() const final { return Total_CSF; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CEff() const final { return Total_CEff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
   * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CpDiff() const final { return Total_CpDiff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional Equivalent Area coefficient.
   * \return Value of the Equivalent Area coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_HeatFluxDiff() const final { return Total_HeatFluxDiff; }

  /*!
   * \brief Set the value of the Equivalent Area coefficient.
   * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
   */
  inline void SetTotal_CpDiff(su2double val_pressure) final { Total_CpDiff = val_pressure; }

  /*!
   * \brief Set the value of the Equivalent Area coefficient.
   * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
   */
  inline void SetTotal_HeatFluxDiff(su2double val_heat) final { Total_HeatFluxDiff = val_heat; }

  /*!
   * \brief Set the value of the custom objective function.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  inline void SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) final {
    Total_Custom_ObjFunc = val_total_custom_objfunc*val_weight;
  }

  /*!
   * \brief Add the value of the custom objective function.
   * \param[in] val_Total_Custom_ObjFunc - Value of the total custom objective function.
   * \param[in] val_weight - Value of the weight for the custom objective function.
   */
  inline void AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) final {
    Total_Custom_ObjFunc += val_total_custom_objfunc*val_weight;
  }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CL() const final { return Total_CL; }

  /*!
   * \author H. Kline
   * \brief Set the total "combo" objective (weighted sum of other values).
   * \param[in] ComboObj - Value of the combined objective.
   */
  inline void SetTotal_ComboObj(su2double ComboObj) final {Total_ComboObj = ComboObj; }

  /*!
   * \author H. Kline
   * \brief Provide the total "combo" objective (weighted sum of other values).
   * \return Value of the "combo" objective values.
   */
  inline su2double GetTotal_ComboObj() const final { return Total_ComboObj; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CD() const final { return Total_CD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMx() const final { return Total_CMx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMy() const final { return Total_CMy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMz() const final { return Total_CMz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x moment coefficient.
   * \return Value of the moment x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPx() const final { return Total_CoPx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y moment coefficient.
   * \return Value of the moment y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPy() const final { return Total_CoPy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z moment coefficient.
   * \return Value of the moment z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CoPz() const final { return Total_CoPz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional x force coefficient.
   * \return Value of the force x coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFx() const final { return Total_CFx; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional y force coefficient.
   * \return Value of the force y coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFy() const final { return Total_CFy; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional z force coefficient.
   * \return Value of the force z coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CFz() const final { return Total_CFz; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional thrust coefficient.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CT() const final { return Total_CT; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional thrust coefficient.
   * \param[in] val_Total_CT - Value of the total thrust coefficient.
   */
  inline void SetTotal_CT(su2double val_Total_CT) final { Total_CT = val_Total_CT; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional torque coefficient.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CQ() const final { return Total_CQ; }

  /*!
   * \brief Provide the total heat load.
   * \return Value of the heat load (viscous contribution).
   */
  inline su2double GetTotal_HeatFlux(void) const final { return Total_Heat; }

  /*!
   * \brief Provide the total heat load.
   * \return Value of the heat load (viscous contribution).
   */
  inline su2double GetTotal_MaxHeatFlux() const final { return Total_MaxHeat; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional torque coefficient.
   * \param[in] val_Total_CQ - Value of the total torque coefficient.
   */
  inline void SetTotal_CQ(su2double val_Total_CQ) final { Total_CQ = val_Total_CQ; }

  /*!
   * \brief Store the total heat load.
   * \param[in] val_Total_Heat - Value of the heat load.
   */
  inline void SetTotal_HeatFlux(su2double val_Total_Heat) final { Total_Heat = val_Total_Heat; }

  /*!
   * \brief Store the total heat load.
   * \param[in] val_Total_Heat - Value of the heat load.
   */
  inline void SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat) final { Total_MaxHeat = val_Total_MaxHeat; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional rotor Figure of Merit.
   * \return Value of the rotor efficiency coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CMerit() const final { return Total_CMerit; }

  /*!
   * \brief Provide the total custom objective function.
   * \return Value of the custom objective function.
   */
  inline su2double GetTotal_Custom_ObjFunc() const final { return Total_Custom_ObjFunc; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CDrag - Value of the total drag coefficient.
   */
  inline void SetTotal_CD(su2double val_Total_CD) final { Total_CD = val_Total_CD; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Inv() const final { return AllBound_CL_Inv; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Inv() const final { return AllBound_CD_Inv; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Inv() const final { return AllBound_CSF_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Inv() const final { return AllBound_CEff_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Inv() const final { return AllBound_CMx_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Inv() const final { return AllBound_CMy_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Inv() const final { return AllBound_CMz_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Inv() const final { return AllBound_CoPx_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Inv() const final { return AllBound_CoPy_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Inv() const final { return AllBound_CoPz_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Inv() const final { return AllBound_CFx_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Inv() const final { return AllBound_CFy_Inv; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Inv() const final { return AllBound_CFz_Inv; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Mnt() const final { return AllBound_CL_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Mnt() const final { return AllBound_CD_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the sideforce coefficient.
   * \return Value of the sideforce coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CSF_Mnt() const final { return AllBound_CSF_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CEff_Mnt() const final { return AllBound_CEff_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMx_Mnt() const final { return AllBound_CMx_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMy_Mnt() const final { return AllBound_CMy_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CMz_Mnt() const final { return AllBound_CMz_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPx_Mnt() const final { return AllBound_CoPx_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPy_Mnt() const final { return AllBound_CoPy_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CoPz_Mnt() const final { return AllBound_CoPz_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFx_Mnt() const final { return AllBound_CFx_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFy_Mnt() const final { return AllBound_CFy_Mnt; }

  /*!
   * \brief Get the inviscid contribution to the efficiency coefficient.
   * \return Value of the efficiency coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CFz_Mnt() const final { return AllBound_CFz_Mnt; }

  /*!
   * \brief Provide the Pressure coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double GetCPressure(unsigned short val_marker, unsigned long val_vertex) const final {
    return CPressure[val_marker][val_vertex];
  }

  /*!
   * \brief Provide the Target Pressure coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) const final {
    return CPressureTarget[val_marker][val_vertex];
  }

  /*!
   * \brief Set the value of the target Pressure coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetCPressureTarget(unsigned short val_marker,
                                 unsigned long val_vertex,
                                 su2double val_pressure) final {
    CPressureTarget[val_marker][val_vertex] = val_pressure;
  }

  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double *GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) const final {
    return CharacPrimVar[val_marker][val_vertex];
  }

  /*!
   * \brief Set the total residual adding the term that comes from the Dual Time Strategy.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SetResidual_DualTime(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep,
                            unsigned short iMesh,
                            unsigned short RunTime_EqSystem) final;

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
                   bool val_update_geo) final;

  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) final;

  /*!
   * \brief Set the freestream pressure.
   * \param[in] Value of freestream pressure.
   */
  inline void SetPressure_Inf(su2double p_inf) final { Pressure_Inf = p_inf; }

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetTemperature_Inf(su2double t_inf) final { Temperature_Inf = t_inf; }

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetDensity_Inf(su2double rho_inf) final { Density_Inf = rho_inf; }

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(CConfig *config) final;

  /*!
   * \brief Update the Beta parameter for the incompressible preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   */
  void SetBeta_Parameter(CGeometry *geometry,
                         CSolver **solver_container,
                         CConfig *config,
                         unsigned short iMesh) final;

  /*!
   * \brief Compute the preconditioner for low-Mach flows.
   * \param[in] iPoint - Index of the grid point
   * \param[in] config - Definition of the particular problem.
   */
  void SetPreconditioner(CConfig *config, unsigned long iPoint) final;

  /*!
   * \brief Value of the total temperature at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total temperature is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is evaluated.
   * \return Value of the total temperature
   */
  inline su2double GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) const final { return Inlet_Ttotal[val_marker][val_vertex]; }

  /*!
   * \brief Value of the total pressure at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total pressure is evaluated.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is evaluated.
   * \return Value of the total pressure
   */
  inline su2double GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) const final { return Inlet_Ptotal[val_marker][val_vertex]; }

  /*!
   * \brief A component of the unit vector representing the flow direction at an inlet boundary.
   * \param[in] val_marker - Surface marker where the flow direction is evaluated
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is evaluated
   * \param[in] val_dim - The component of the flow direction unit vector to be evaluated
   * \return Component of a unit vector representing the flow direction.
   */
  inline su2double GetInlet_FlowDir(unsigned short val_marker,
                                    unsigned long val_vertex,
                                    unsigned short val_dim) const final {
    return Inlet_FlowDir[val_marker][val_vertex][val_dim];
  }

  /*!
   * \brief Set a component of the unit vector representing the flow direction at an inlet boundary.
   * \param[in] val_marker - Surface marker where the flow direction is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is set.
   * \param[in] val_dim - The component of the flow direction unit vector to be set
   * \param[in] val_flowdir - Component of a unit vector representing the flow direction.
   */
  inline void SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir) final {
    /*--- Since this call can be accessed indirectly using python, do some error
     * checking to prevent segmentation faults ---*/
    if (val_marker >= nMarker)
      SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
    else if (Inlet_FlowDir == NULL || Inlet_FlowDir[val_marker] == NULL)
      SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
    else if (val_vertex >= nVertex[val_marker])
      SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
    else
      Inlet_FlowDir[val_marker][val_vertex][val_dim] = val_flowdir;
  }

  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(CConfig* config, unsigned short iMarker) final;

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(su2double *val_inlet,
                        unsigned short iMarker,
                        unsigned long iVertex) final;

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
                             CConfig *config) const final;

  /*!
   * \brief A virtual member.
   */
  void GetOutlet_Properties(CGeometry *geometry,
                            CConfig *config,
                            unsigned short iMesh,
                            bool Output) final;

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  inline void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex) final {
    int iVar;

    for( iVar = 0; iVar < nPrimVar+1; iVar++){
      if( SlidingState[val_marker][val_vertex][iVar] != NULL )
        delete [] SlidingState[val_marker][val_vertex][iVar];
    }

    for( iVar = 0; iVar < nPrimVar+1; iVar++)
      SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
  }




  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  inline void SetSlidingState(unsigned short val_marker,
                              unsigned long val_vertex,
                              unsigned short val_state,
                              unsigned long donor_index,
                              su2double component) final {
    SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
  }


  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  inline void SetnSlidingStates(unsigned short val_marker,
                                unsigned long val_vertex,
                                int value) final { SlidingStateNodes[val_marker][val_vertex] = value; }

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  inline int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex) const final{
    return SlidingStateNodes[val_marker][val_vertex];
  }

  /*!
   * \brief Get the outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] val_state  - requested state component
   * \param[in] donor_index- index of the donor node to get
   */
  inline su2double GetSlidingState(unsigned short val_marker,
                                   unsigned long val_vertex,
                                   unsigned short val_state,
                                   unsigned long donor_index) const final {
    return SlidingState[val_marker][val_vertex][val_state][donor_index];
  }

  /*!
   * \brief Compute the global error measures (L2, Linf) for verification cases.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVerificationError(CGeometry *geometry, CConfig *config) final;

};
