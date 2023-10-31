/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
 * \author F. Palacios, T. Economon
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

#include "CFVMFlowSolverBase.hpp"
#include "../variables/CEulerVariable.hpp"

/*!
 * \class CEulerSolver
 * \ingroup Euler_Equations
 * \brief Class for compressible inviscid flow problems, serves as base for Navier-Stokes/RANS.
 * \author F. Palacios
 */
class CEulerSolver : public CFVMFlowSolverBase<CEulerVariable, ENUM_REGIME::COMPRESSIBLE> {
protected:
  using BaseClass = CFVMFlowSolverBase<CEulerVariable, ENUM_REGIME::COMPRESSIBLE>;

  su2double
  Prandtl_Lam = 0.0,    /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb = 0.0;   /*!< \brief Turbulent Prandtl number. */

  su2double AllBound_CEquivArea_Inv=0.0; /*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
  vector<su2double> CEquivArea_Mnt;      /*!< \brief Equivalent area (inviscid contribution) for each boundary. */
  vector<su2double> CEquivArea_Inv;      /*!< \brief Equivalent area (inviscid contribution) for each boundary. */

  vector<su2double> Inflow_MassFlow;     /*!< \brief Mass flow rate for each boundary. */
  vector<su2double> Exhaust_MassFlow;    /*!< \brief Mass flow rate for each boundary. */
  vector<su2double> Inflow_Pressure;     /*!< \brief Fan face pressure for each boundary. */
  vector<su2double> Inflow_Mach;         /*!< \brief Fan face mach number for each boundary. */
  vector<su2double> Inflow_Area;         /*!< \brief Boundary total area. */
  vector<su2double> Exhaust_Area;        /*!< \brief Boundary total area. */
  vector<su2double> Exhaust_Pressure;    /*!< \brief Fan face pressure for each boundary. */
  vector<su2double> Exhaust_Temperature; /*!< \brief Fan face mach number for each boundary. */
  su2double
  Inflow_MassFlow_Total = 0.0,   /*!< \brief Mass flow rate for each boundary. */
  Exhaust_MassFlow_Total = 0.0,  /*!< \brief Mass flow rate for each boundary. */
  Inflow_Pressure_Total = 0.0,   /*!< \brief Fan face pressure for each boundary. */
  Inflow_Mach_Total = 0.0,       /*!< \brief Fan face mach number for each boundary. */
  InverseDesign = 0.0;           /*!< \brief Inverse design functional for each boundary. */
  vector<vector<unsigned long> > DonorGlobalIndex;  /*!< \brief Value of the donor global index. */
  vector<su2activematrix> DonorPrimVar;       /*!< \brief Value of the donor variables at each boundary. */
  vector<vector<su2double> > ActDisk_DeltaP;  /*!< \brief Value of the Delta P. */
  vector<vector<su2double> > ActDisk_DeltaT;  /*!< \brief Value of the Delta T. */

  su2activevector
  ActDisk_R;         /*!< \brief Value of the actuator disk Radius. */
  su2activematrix
  ActDisk_C,         /*!< \brief Value of the actuator disk Center. */
  ActDisk_Axis;      /*!< \brief Value of the actuator disk Axis. */
  vector<vector<su2double> > ActDisk_Fa; /*!< \brief Value of the actuator disk Axial Force per Unit Area. */
  vector<vector<su2double> > ActDisk_Fx; /*!< \brief Value of the actuator disk X component of the radial and tangential forces per Unit Area resultant. */
  vector<vector<su2double> > ActDisk_Fy; /*!< \brief Value of the actuator disk Y component of the radial and tangential forces per Unit Area resultant. */
  vector<vector<su2double> > ActDisk_Fz; /*!< \brief Value of the actuator disk Z component of the radial and tangential forces per Unit Area resultant. */

  su2double
  Total_CL_Prev = 0.0,        /*!< \brief Total lift coefficient for all the boundaries (fixed lift mode). */
  Total_SolidCD = 0.0,        /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CD_Prev = 0.0,        /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_NetThrust = 0.0,      /*!< \brief Total drag coefficient for all the boundaries. */
  Total_Power = 0.0,          /*!< \brief Total drag coefficient for all the boundaries. */
  Total_ReverseFlow = 0.0,    /*!< \brief Total drag coefficient for all the boundaries. */
  Total_IDC = 0.0,            /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDC_Mach = 0.0,       /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDR = 0.0,            /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_DC60 = 0.0,           /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_MFR = 0.0,            /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Prop_Eff = 0.0,       /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_ByPassProp_Eff = 0.0, /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Adiab_Eff = 0.0,      /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Poly_Eff = 0.0,       /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_CMx_Prev = 0.0,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMy_Prev = 0.0,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMz_Prev = 0.0,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_AeroCD = 0.0,         /*!< \brief Total aero drag coefficient for all the boundaries. */
  Total_CEquivArea = 0.0;     /*!< \brief Total Equivalent Area coefficient for all the boundaries. */

  su2double AoA_Prev,  /*!< \brief Old value of the angle of attack (monitored). */
  AoA_inc;
  bool Start_AoA_FD = false,  /*!< \brief Boolean for start of finite differencing for FixedCL mode */
  End_AoA_FD = false,         /*!< \brief Boolean for end of finite differencing for FixedCL mode */
  Update_AoA = false;         /*!< \brief Boolean to signal Angle of Attack Update */
  unsigned long Iter_Update_AoA = 0; /*!< \brief Iteration at which AoA was updated last */
  su2double dCL_dAlpha;              /*!< \brief Value of dCL_dAlpha used to control CL in fixed CL mode */
  unsigned long BCThrust_Counter;

  vector<CFluidModel*> FluidModel;   /*!< \brief fluid model used in the solver. */

  /*--- Turbomachinery Solver Variables ---*/

  vector<su2activematrix> AverageFlux;
  vector<su2activematrix> SpanTotalFlux;
  vector<su2activematrix> AverageVelocity;
  vector<su2activematrix> AverageTurboVelocity;
  vector<su2activematrix> OldAverageTurboVelocity;
  vector<su2activematrix> ExtAverageTurboVelocity;
  su2activematrix AveragePressure;
  su2activematrix OldAveragePressure;
  su2activematrix RadialEquilibriumPressure;
  su2activematrix ExtAveragePressure;
  su2activematrix AverageDensity;
  su2activematrix OldAverageDensity;
  su2activematrix ExtAverageDensity;
  su2activematrix AverageNu;
  su2activematrix AverageKine;
  su2activematrix AverageOmega;
  su2activematrix ExtAverageNu;
  su2activematrix ExtAverageKine;
  su2activematrix ExtAverageOmega;

  su2activematrix DensityIn;
  su2activematrix PressureIn;
  vector<su2activematrix> TurboVelocityIn;
  su2activematrix DensityOut;
  su2activematrix PressureOut;
  vector<su2activematrix> TurboVelocityOut;
  su2activematrix KineIn;
  su2activematrix OmegaIn;
  su2activematrix NuIn;
  su2activematrix KineOut;
  su2activematrix OmegaOut;
  su2activematrix NuOut;

  vector<su2matrix<complex<su2double> > > CkInflow, CkOutflow1, CkOutflow2;

  /*--- End of Turbomachinery Solver Variables ---*/

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
   * \brief Compute Ducros Sensor for Roe Dissipation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the Fan face Mach number.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  void GetPower_Properties(CGeometry *geometry, CConfig *config,
                           unsigned short iMesh, bool Output);

  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config);

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
   * \brief Read the actuator disk input file for the VARIABLE_LOAD type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void ReadActDisk_InputFile(CGeometry *geometry, CSolver **solver_container,
                           CConfig *config, unsigned short iMesh, bool Output);

  /*!
   * \brief Compute the max eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the undivided laplacian for the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUndivided_Laplacian(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the dissipation sensor for centered schemes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCentered_Dissipation_Sensor(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetRoe_Dissipation(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \return - The number of non-physical points.
   */
  virtual unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                               const CConfig *config);

  /*!
   * \brief Set gradients of coefficients for fixed CL mode
   * \param[in] config - Definition of the particular problem.
   */
  void SetCoefficient_Gradients(CConfig *config) const;

  /*!
   * \brief Instantiate a SIMD numerics object.
   * \param[in] solvers - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void InstantiateEdgeNumerics(const CSolver* const* solvers, const CConfig* config) final;

  /*!
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh);

  /*!
   * \brief Set reference values for pressure, forces, etc.
   */
  void SetReferenceValues(const CConfig& config) final;

public:
  CEulerSolver() = delete;

  /*!
   * \brief Main constructor of this class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Grid level.
   * \param[in] navier_stokes - True when the constructor is called by the derived class CNSSolver.
   */
  CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, const bool navier_stokes = false);

  /*!
   * \brief Destructor of the class.
   */
  ~CEulerSolver(void) override;

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline CFluidModel* GetFluidModel(void) const final { return FluidModel[omp_get_thread_num()]; }

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
                    unsigned long Iteration) final;

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
   * \brief Recompute the extrapolated quantities, after MUSCL reconstruction,
   *        in a more thermodynamically consistent way.
   * \note This method is static to improve the chances of it being used in a
   *       thread-safe manner.
   * \param[in,out] fluidModel - The fluid model.
   * \param[in] nDim - Number of physical dimensions.
   * \param[in,out] primitive - Primitive variables.
   * \param[out] secondary - Secondary variables.
   */
  static void ComputeConsistentExtrapolation(CFluidModel *fluidModel,
                                             unsigned short nDim,
                                             su2double *primitive,
                                             su2double *secondary);

  /*!
   * \brief Apply low Mach number correction to the primitives at two points,
   *        usually connected by an edge.
   * \note This method is static to improve the chances of it being used in a
   *       thread-safe manner.
   * \param[in,out] fluidModel - The fluid model.
   * \param[in] nDim - Number of physical dimensions.
   * \param[in,out] primitive_i - Primitive variables at point i.
   * \param[in,out] primitive_j - Primitive variables at point j.
   */
  static void LowMachPrimitiveCorrection(CFluidModel *fluidModel,
                                         unsigned short nDim,
                                         su2double *primitive_i,
                                         su2double *primitive_j);

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
                       unsigned short iMesh) override;

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
   * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iPoint - Index of the grid point.
   * \param[in] delta - Volume over delta t.
   * \param[in,out] preconditioner - The preconditioner matrix, must be allocated outside.
   */
  void SetPreconditioner(const CConfig *config, unsigned long iPoint,
                         su2double delta, su2activematrix& preconditioner) const;

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void Evaluate_ObjFunc(const CConfig *config, CSolver **solver) override;

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
   * \brief Impose the engine inflow boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
  void BC_ActDisk_Inlet(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) final;

  /*!
   * \brief Impose the engine exhaust boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Outlet(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *conv_numerics,
                         CNumerics *visc_numerics,
                         CConfig *config,
                         unsigned short val_marker) final;

  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  void BC_ActDisk(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker,
                  bool val_inlet_surface) final;

  /*!
   * \brief Impose an actuator disk with variable load boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  void BC_ActDisk_VariableLoad(CGeometry *geometry,
                               CSolver **solver_container,
                               CNumerics *conv_numerics,
                               CNumerics *visc_numerics,
                               CConfig *config,
                               unsigned short val_marker,
                               bool val_inlet_surface);

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
  void BC_Riemann(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker) final;

  /*!
   * \brief Impose the boundary condition using characteristic recostruction.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_TurboRiemann(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *conv_numerics,
                       CNumerics *visc_numerics,
                       CConfig *config,
                       unsigned short val_marker) final;

  /*!
   * \brief It computes Fourier transformation for the needed quantities along the pitch for each span in turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  void PreprocessBC_Giles(CGeometry *geometry,
                          CConfig *config,
                          CNumerics *conv_numerics,
                          unsigned short marker_flag) final;

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
  void BC_Giles(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
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
   * \brief Impose a supersonic inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Inlet(CGeometry *geometry,
                           CSolver **solver_container,
                           CNumerics *conv_numerics,
                           CNumerics *visc_numerics,
                           CConfig *config,
                           unsigned short val_marker) final;

  /*!
   * \brief Impose a supersonic outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Outlet(CGeometry *geometry,
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
   * \brief Impose the nacelle inflow boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Inflow(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) final;

  /*!
   * \brief Impose the ancelle exhaust boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Exhaust(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *conv_numerics,
                         CNumerics *visc_numerics,
                         CConfig *config,
                         unsigned short val_marker) final;

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
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ClassicalRK4_Iteration(CGeometry *geometry,
                              CSolver **solver_container,
                              CConfig *config,
                              unsigned short iRKStep) final;

  /*!
   * \brief Check for convergence of the Fixed CL mode to the target CL
   * \param[in] config - Definition of the particular problem.
   * \param[in] convergence - boolean for whether the solution is converged
   * \return boolean for whether the Fixed CL mode is converged to target CL
   */
  bool FixedCL_Convergence(CConfig *config, bool convergence) final;

  /*!
   * \brief Checking whether fixed CL mode in finite-differencing mode
   * \return boolean for whether the Fixed CL mode is currently in finite-differencing mode
   */
  inline bool GetStart_AoA_FD(void) const final { return Start_AoA_FD; }

  /*!
   * \brief Checking whether fixed CL mode in finite-differencing mode
   * \return boolean for whether the Fixed CL mode is currently in finite-differencing mode
   */
  inline bool GetEnd_AoA_FD(void) const final { return End_AoA_FD; }

    /*!
   * \brief Get the iteration of the last AoA update (Fixed CL Mode)
   * \return value for the last iteration that the AoA was updated
   */
  inline unsigned long GetIter_Update_AoA(void) const final { return Iter_Update_AoA; }

  /*!
   * \brief Get the AoA before the most recent update
   * \return value of the AoA before most recent update
   */
  inline su2double GetPrevious_AoA(void) const final { return AoA_Prev; }

  /*!
   * \brief Get the CL Driver's control command
   * \return value of CL Driver control command (AoA_inc)
   */
  inline su2double GetAoA_inc(void) const final { return AoA_inc; }

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
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  inline su2double GetInflow_MassFlow(unsigned short val_marker) const final { return Inflow_MassFlow[val_marker]; }

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  inline su2double GetExhaust_MassFlow(unsigned short val_marker) const final { return Exhaust_MassFlow[val_marker]; }

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the fan face pressure on the surface <i>val_marker</i>.
   */
  inline su2double GetInflow_Pressure(unsigned short val_marker) const final { return Inflow_Pressure[val_marker]; }

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the fan face mach on the surface <i>val_marker</i>.
   */
  inline su2double GetInflow_Mach(unsigned short val_marker) const final { return Inflow_Mach[val_marker]; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional aero CD.
   * \return Value of the Aero CD coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_AeroCD() const final { return Total_AeroCD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional Near-Field pressure coefficient.
   * \return Value of the NearField pressure coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CNearFieldOF() const final { return Total_CNearFieldOF; }

  /*!
   * \brief Set the value of the Aero drag.
   * \param[in] val_cequivarea - Value of the aero drag.
   */
  inline void SetTotal_AeroCD(su2double val_aerocd) final { Total_AeroCD = val_aerocd; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_NetThrust() const final { return Total_NetThrust; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Power() const final { return Total_Power; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_SolidCD() const final { return Total_SolidCD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_ReverseFlow() const final { return Total_ReverseFlow; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_MFR() const final { return Total_MFR; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Prop_Eff() const final { return Total_Prop_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_ByPassProp_Eff() const final { return Total_ByPassProp_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Adiab_Eff() const final { return Total_Adiab_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Poly_Eff() const final { return Total_Poly_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_IDC() const final { return Total_IDC; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_IDC_Mach() const final { return Total_IDC_Mach; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_IDR() const final { return Total_IDR; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_DC60() const final { return Total_DC60; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_NetThrust(su2double val_Total_NetThrust) final { Total_NetThrust = val_Total_NetThrust; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Power(su2double val_Total_Power) final { Total_Power = val_Total_Power; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_SolidCD(su2double val_Total_SolidCD) final { Total_SolidCD = val_Total_SolidCD; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) final { Total_ReverseFlow = val_Total_ReverseFlow; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_MFR(su2double val_Total_MFR) final { Total_MFR = val_Total_MFR; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) final { Total_Prop_Eff = val_Total_Prop_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) final { Total_ByPassProp_Eff = val_Total_ByPassProp_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) final { Total_Adiab_Eff = val_Total_Adiab_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) final { Total_Poly_Eff = val_Total_Poly_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_IDC(su2double val_Total_IDC) final { Total_IDC = val_Total_IDC; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) final { Total_IDC_Mach = val_Total_IDC_Mach; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_IDR(su2double val_Total_IDR) final { Total_IDR = val_Total_IDR; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_DC60(su2double val_Total_DC60) final { Total_DC60 = val_Total_DC60; }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline unsigned long GetDonorGlobalIndex(unsigned short val_marker,
                                           unsigned long val_vertex) const final {
    return DonorGlobalIndex[val_marker][val_vertex];
  }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetDonorGlobalIndex(unsigned short val_marker,
                                  unsigned long val_vertex,
                                  unsigned long val_index) final {
    DonorGlobalIndex[val_marker][val_vertex] = val_index;
  }

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config) final;

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
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(const CConfig *config) final;

  /*!
   * \brief Initilize turbo containers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void InitTurboContainers(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_TurboSolution(CConfig *config) final;

  /*!
   * \brief It computes average quantities along the span for turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  void PreprocessAverage(CSolver **solver,
                         CGeometry *geometry,
                         CConfig *config,
                         unsigned short marker_flag) final;

  /*!
   * \brief It computes average quantities along the span for turbomachinery analysis.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] marker_flag - Surface marker flag where the function is applied.
   */
  void TurboAverageProcess(CSolver **solver,
                           CGeometry *geometry,
                           CConfig *config,
                           unsigned short marker_flag) final;

  /*!
   * \brief it performs a mixed out average of the nodes of a boundary.
   * \param[in] val_init_pressure -  initial pressure value
   * \param[in] val_Averaged_Flux - flux averaged values.
   * \param[in] val_normal - normal vector.
   * \param[in] pressure_mix - value of the mixed-out avaraged pressure.
   * \param[in] density_miz - value of the mixed-out avaraged density.
   */
  void MixedOut_Average (CConfig *config,
                         su2double val_init_pressure,
                         const su2double *val_Averaged_Flux,
                         const su2double *val_normal,
                         su2double& pressure_mix,
                         su2double& density_mix);

  /*!
   * \brief It gathers into the master node average quantities at inflow and outflow needed for turbomachinery analysis.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void GatherInOutAverageValues(CConfig *config, CGeometry *geometry) final;

  /*!
   * \brief it take a velocity in the cartesian reference of framework and transform into the turbomachinery frame of reference.
   * \param[in] cartesianVelocity - cartesian components of velocity vector.
   * \param[in] turboNormal - normal vector in the turbomachinery frame of reference.
   * \param[in] turboVelocity - velocity vector in the turbomachinery frame of reference.
   */
  inline void ComputeTurboVelocity(const su2double *cartesianVelocity,
                                   const su2double *turboNormal,
                                   su2double *turboVelocity,
                                   unsigned short marker_flag,
                                   unsigned short kind_turb){

    if ((kind_turb == AXIAL && nDim == 3) || (kind_turb == CENTRIPETAL_AXIAL && marker_flag == OUTFLOW) || (kind_turb == AXIAL_CENTRIFUGAL && marker_flag == INFLOW) ){
      turboVelocity[2] =  turboNormal[0]*cartesianVelocity[0] + cartesianVelocity[1]*turboNormal[1];
      turboVelocity[1] =  turboNormal[0]*cartesianVelocity[1] - turboNormal[1]*cartesianVelocity[0];
      turboVelocity[0] = cartesianVelocity[2];
    }
    else{
      turboVelocity[0] =  turboNormal[0]*cartesianVelocity[0] + cartesianVelocity[1]*turboNormal[1];
      turboVelocity[1] =  turboNormal[0]*cartesianVelocity[1] - turboNormal[1]*cartesianVelocity[0];
      if (marker_flag == INFLOW){
        turboVelocity[0] *= -1.0;
        turboVelocity[1] *= -1.0;
      }
      if(nDim == 3)
        turboVelocity[2] = cartesianVelocity[2];
    }
  }

  /*!
   * \brief it take a velocity in the cartesian reference of framework and transform into the turbomachinery frame of reference.
   * \param[in] cartesianVelocity - cartesian components of velocity vector.
   * \param[in] turboNormal - normal vector in the turbomachinery frame of reference.
   * \param[in] turboVelocity - velocity vector in the turbomachinery frame of reference.
   */
  inline void ComputeBackVelocity(const su2double *turboVelocity,
                                  const su2double *turboNormal,
                                  su2double *cartesianVelocity,
                                  unsigned short marker_flag,
                                  unsigned short kind_turb){

    if ((kind_turb == AXIAL && nDim == 3) || (kind_turb == CENTRIPETAL_AXIAL && marker_flag == OUTFLOW) || (kind_turb == AXIAL_CENTRIFUGAL && marker_flag == INFLOW)){
      cartesianVelocity[0] = turboVelocity[2]*turboNormal[0] - turboVelocity[1]*turboNormal[1];
      cartesianVelocity[1] = turboVelocity[2]*turboNormal[1] + turboVelocity[1]*turboNormal[0];
      cartesianVelocity[2] = turboVelocity[0];
    }
    else{
      cartesianVelocity[0] =  turboVelocity[0]*turboNormal[0] - turboVelocity[1]*turboNormal[1];
      cartesianVelocity[1] =  turboVelocity[0]*turboNormal[1] + turboVelocity[1]*turboNormal[0];

      if (marker_flag == INFLOW){
        cartesianVelocity[0] *= -1.0;
        cartesianVelocity[1] *= -1.0;
      }

      if(nDim == 3)
        cartesianVelocity[2] = turboVelocity[2];
    }
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Density on the surface <i>val_marker</i>.
   */
  inline su2double GetAverageDensity(unsigned short valMarker, unsigned short valSpan) const final {
    return AverageDensity[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average pressure at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Pressure on the surface <i>val_marker</i>.
   */
  inline su2double GetAveragePressure(unsigned short valMarker, unsigned short valSpan) const final {
    return AveragePressure[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average turbo velocity average at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  inline const su2double* GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan) const final {
    return AverageTurboVelocity[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Nu on the surface <i>val_marker</i>.
   */
  inline su2double GetAverageNu(unsigned short valMarker, unsigned short valSpan) const final {
    return AverageNu[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Kine on the surface <i>val_marker</i>.
   */
  inline su2double GetAverageKine(unsigned short valMarker, unsigned short valSpan) const final {
    return AverageKine[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Omega on the surface <i>val_marker</i>.
   */
  inline su2double GetAverageOmega(unsigned short valMarker, unsigned short valSpan) const final {
    return AverageOmega[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Nu on the surface <i>val_marker</i>.
   */
  inline su2double GetExtAverageNu(unsigned short valMarker, unsigned short valSpan) const final {
    return ExtAverageNu[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Kine on the surface <i>val_marker</i>.
   */
  inline su2double GetExtAverageKine(unsigned short valMarker, unsigned short valSpan) const final {
    return ExtAverageKine[valMarker][valSpan];
  }

  /*!
   * \brief Provide the average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average turbulent Omega on the surface <i>val_marker</i>.
   */
  inline su2double GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan) const final {
    return ExtAverageOmega[valMarker][valSpan];
  }

  /*!
   * \brief Set the external average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valDensity - value to set.
   */
  inline void SetExtAverageDensity(unsigned short valMarker,
                                   unsigned short valSpan,
                                   su2double valDensity) final {
    ExtAverageDensity[valMarker][valSpan] = valDensity;
  }

  /*!
   * \brief Set the external average density at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valPressure - value to set.
   */
  inline void SetExtAveragePressure(unsigned short valMarker,
                                    unsigned short valSpan,
                                    su2double valPressure) final {
    ExtAveragePressure[valMarker][valSpan] = valPressure;
  }

  /*!
   * \brief Set the external the average turbo velocity average at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \return Value of the Average Total Pressure on the surface <i>val_marker</i>.
   */
  inline void SetExtAverageTurboVelocity(unsigned short valMarker,
                                         unsigned short valSpan,
                                         unsigned short valIndex,
                                         su2double valTurboVelocity) final {
    ExtAverageTurboVelocity[valMarker][valSpan][valIndex] = valTurboVelocity;
  }

  /*!
   * \brief Set the external average turbulent Nu at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valNu - value to set.
   */
  inline void SetExtAverageNu(unsigned short valMarker,
                              unsigned short valSpan,
                              su2double valNu) final {
    ExtAverageNu[valMarker][valSpan] = valNu;
  }

  /*!
   * \brief Set the external average turbulent Kine at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valKine - value to set.
   */
  inline void SetExtAverageKine(unsigned short valMarker,
                                unsigned short valSpan,
                                su2double valKine) final {
    ExtAverageKine[valMarker][valSpan] = valKine;
  }

  /*!
   * \brief Set the external average turbulent Omega at the boundary of interest.
   * \param[in] val_marker - bound marker.
   * \param[in] val_Span   - value of the Span.
   * \param[in] valOmega - value to set.
   */
  inline void SetExtAverageOmega(unsigned short valMarker,
                                 unsigned short valSpan,
                                 su2double valOmega) final {
    ExtAverageOmega[valMarker][valSpan] = valOmega;
  }

  /*!
   * \brief Provide the inlet density to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return DensityIn[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the inlet pressure to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of inlet pressure.
   */
  inline su2double GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return PressureIn[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the inlet normal velocity to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet normal velocity.
   */
  inline const su2double* GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return TurboVelocityIn[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the outlet density to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet density.
   */
  inline su2double GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return DensityOut[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the outlet pressure to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet pressure.
   */
  inline su2double GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return PressureOut[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the outlet normal velocity to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the outlet normal velocity.
   */
  inline const su2double* GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return TurboVelocityOut[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the inlet turbulent kei to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetKineIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return KineIn[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the inlet turbulent omega to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return OmegaIn[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the inlet turbulent nu to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetNuIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return NuIn[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the outlet turbulent kei to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetKineOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return KineOut[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the outlet turbulent omega to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return OmegaOut[inMarkerTP][valSpan];
  }

  /*!
   * \brief Provide the outlet turbulent nu to check convergence of conservative mixing-plane.
   * \param[in] inMarkerTP - bound marker.
   * \return Value of the inlet density.
   */
  inline su2double GetNuOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
    return NuOut[inMarkerTP][valSpan];
  }

  /*!
   * \brief Set inlet density.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetDensityIn(su2double value,
                           unsigned short inMarkerTP,
                           unsigned short valSpan) final {
    DensityIn[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set inlet pressure.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetPressureIn(su2double value,
                            unsigned short inMarkerTP,
                            unsigned short valSpan) final {
    PressureIn[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set inlet normal velocity.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetTurboVelocityIn(const su2double *value,
                                 unsigned short inMarkerTP,
                                 unsigned short valSpan) final {
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      TurboVelocityIn[inMarkerTP][valSpan][iDim] = value[iDim];
  }

  /*!
   * \brief Set outlet density.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetDensityOut(su2double value,
                            unsigned short inMarkerTP,
                            unsigned short valSpan) final {
    DensityOut[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set outlet pressure.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetPressureOut(su2double value,
                             unsigned short inMarkerTP,
                             unsigned short valSpan) final {
    PressureOut[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set outlet normal velocity.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetTurboVelocityOut(const su2double *value,
                                  unsigned short inMarkerTP,
                                  unsigned short valSpan) final {
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      TurboVelocityOut[inMarkerTP][valSpan][iDim] = value[iDim];
  }

  /*!
   * \brief Set inlet turbulent kei.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetKineIn(su2double value,
                        unsigned short inMarkerTP,
                        unsigned short valSpan) final {
    KineIn[inMarkerTP][valSpan] = value;
  }
  /*!
   * \brief Set inlet turbulent omega.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetOmegaIn(su2double value,
                        unsigned short inMarkerTP,
                        unsigned short valSpan) final {
    OmegaIn[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set inlet turbulent Nu.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetNuIn(su2double value,
                      unsigned short inMarkerTP,
                      unsigned short valSpan) final {
    NuIn[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set outlet turbulent kei.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetKineOut(su2double value,
                         unsigned short inMarkerTP,
                         unsigned short valSpan) final {
    KineOut[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set Outlet turbulent omega.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetOmegaOut(su2double value,
                          unsigned short inMarkerTP,
                          unsigned short valSpan) final {
    OmegaOut[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Set outlet turbulent Nu.
   * \param[in] value      - turboperformance value to set.
   * \param[in] inMarkerTP - turboperformance marker.
   */
  inline void SetNuOut(su2double value,
                       unsigned short inMarkerTP,
                       unsigned short valSpan) final {
    NuOut[inMarkerTP][valSpan] = value;
  }

  /*!
   * \brief Print verification error to screen.
   * \param[in] config - Definition of the particular problem.
   */
  void PrintVerificationError(const CConfig* config) const final;

  /*!
   * \brief The Euler and NS solvers support MPI+OpenMP.
   */
  inline bool GetHasHybridParallel() const final { return true; }

};
