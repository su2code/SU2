/*!
 * \file CNEMOEulerSolver.hpp
 * \brief Headers of the CNEMOEulerSolver class
 * \author S. R. Copeland, F. Palacios, W. Maier.
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

#include "../variables/CNEMOEulerVariable.hpp"
#include "../fluid/CNEMOGas.hpp"
#include "CFVMFlowSolverBase.hpp"

/*!
 * \class CNEMOEulerSolver
 * \brief Main class for defining the NEMO Euler's flow solver.
 * \ingroup Euler_Equations
 * \author S. R. Copeland, F. Palacios, W. Maier.
 * \version 8.0.0 "Harrier"
 */
class CNEMOEulerSolver : public CFVMFlowSolverBase<CNEMOEulerVariable, ENUM_REGIME::COMPRESSIBLE> {
protected:

  su2double
  Prandtl_Lam = 0.0,              /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb = 0.0;             /*!< \brief Turbulent Prandtl number. */

  unsigned short nSpecies;        /*!< \brief Number of species in the gas mixture. */

  su2double
  Energy_ve_Inf,                  /*!< \brief Vib.-el. free stream energy. */
  Temperature_ve_Inf;             /*!< \brief Vib.-el. free stream temperature. */
  const su2double *MassFrac_Inf;  /*!< \brief Free stream species mass fraction. */

  su2double *Source;              /*!< \brief Auxiliary vector to store source terms. */

  unsigned long ErrorCounter = 0; /*!< \brief Counter for number of un-physical states. */

  su2double Global_Delta_Time = 0.0, /*!< \brief Time-step for TIME_STEPPING time marching strategy. */
  Global_Delta_UnstTimeND = 0.0;     /*!< \brief Unsteady time step for the dual time strategy. */

  CNEMOGas  *FluidModel;          /*!< \brief fluid model used in the solver */

  CNEMOEulerVariable* node_infty = nullptr;

  /*!
   * \brief Set the maximum value of the eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute a pressure sensor switch.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  inline void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief Set the value of undivided laplacian.
   * \param[in] val_und_lapl_i Undivided laplacian at point i.
   * \param[in] val_und_lapl_j Undivided laplacian at point j.
   */
  inline void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief Set reference values for pressure, forces, etc.
   */
  void SetReferenceValues(const CConfig& config) final;

public:
  CNEMOEulerSolver() = delete;

  /*!
   * \brief Contructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, const bool navier_stokes = false);

  /*!
   * \brief Destructor of the class.
   */
  ~CNEMOEulerSolver(void) override;

  /*!
   * \brief Compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container,
                    CConfig *config, unsigned short iMesh, unsigned long Iteration) final;

  /*!
   * \brief Compute the spatial integration using a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics,
                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) final;

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                       CConfig *config, unsigned short iMesh) final;

  /*!
   * \brief Recompute the extrapolated quantities, after MUSCL reconstruction,
   *        in a more thermodynamically consistent way.
   * \param[in] V - primitve variables.
   * \param[out] d*dU - reconstructed secondaryvariables.
   * \param[out] val_eves - reconstructed eve per species.
   * \param[out] val_cvves - reconstructed cvve per species.
   * \param[out] Gamma - reconstructed gamma.
   */
  static su2double ComputeConsistentExtrapolation(CNEMOGas *fluidmodel, unsigned short nSpecies, su2double *V,
                                                  su2double* dPdU, su2double* dTdU, su2double* dTvedU,
                                                  su2double* val_eves, su2double* val_cvves);
  /*!
   * \brief Source term integration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                       CConfig *config, unsigned short iMesh) final;

  /*!
   * \brief Preprocessing actions common to the Euler and NS solvers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void CommonPreprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                           unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                     unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Computes primitive variables.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  virtual unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                               CConfig *config, bool Output);

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over
   * a nonlinear iteration for stability.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeUnderRelaxationFactor(const CConfig *config) final;

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh);

  /*!
   * \brief Set all the conserved variables from the primitive vector..
   */
  void RecomputeConservativeVector(su2double *U, const su2double *V) const;

  /*!
   * \brief Check for unphysical points.
   * \return Boolean value of physical point
   */
  bool CheckNonPhys(const su2double *V) const;

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline CNEMOGas* GetFluidModel(void) const final { return FluidModel;}

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
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

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
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) final;

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
                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

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
                           CConfig *config, unsigned short val_marker) override;
  /*!
   * \brief Impose the supersonic outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics,
                            CConfig *config, unsigned short val_marker) override;
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
                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

  /*!
   * \brief Update the solution using an explicit Euler scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) final;

  /*!
   * \brief Update the solution using a general explicit RK scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Runge-Kutta step.
   */
  void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                            CConfig *config, unsigned short iRKStep) final;

  /*!
   * \brief Update the solution using the classical RK4 scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Runge-Kutta step.
   */
  void ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                              CConfig *config, unsigned short iRKStep) final;

  /*!
   * \brief Prepare an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PrepareImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) final;

  /*!
   * \brief Complete an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) final;

  /*!
   * \brief Print verification error to screen.
   * \param[in] config - Definition of the particular problem.
   */
  void PrintVerificationError(const CConfig* config) const final { }

  /*!
   * \brief Compute the Pressure sensor for NEMO schemes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPressureDiffusionSensor(CGeometry *geometry, CConfig *config);

};
