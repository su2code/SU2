/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
 * \author F. Palacios, T. Economon
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
#include "../variables/CEulerVariable.hpp"

/*!
 * \class CSolver
 * \brief Main class for defining the PDE solution, it requires
 * a child class for each particular solver (Euler, Navier-Stokes, etc.)
 * \author F. Palacios
 */
class CEulerSolver : public CSolver {
protected:
  
  su2double
  Mach_Inf,         /*!< \brief Mach number at the infinity. */
  Density_Inf,      /*!< \brief Density at the infinity. */
  Energy_Inf,       /*!< \brief Energy at the infinity. */
  Temperature_Inf,  /*!< \brief Energy at the infinity. */
  Pressure_Inf,     /*!< \brief Pressure at the infinity. */
  *Velocity_Inf;    /*!< \brief Flow Velocity vector at the infinity. */
  
  su2double
  *CD_Inv,      /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Inv,      /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,     /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,     /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,     /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,     /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Inv,    /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Inv,    /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Inv,    /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,     /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,     /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,     /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Inv,   /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Inv,   /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Inv,  /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Inv, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Inv,  /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Inv,  /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Inv,  /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Inv,  /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Inv,  /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Inv,  /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Inv,         /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Inv,       /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Inv,           /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Inv,           /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  *CEquivArea_Inv,   /*!< \brief Equivalent area (inviscid contribution) for each boundary. */
  *CNearFieldOF_Inv, /*!< \brief Near field pressure (inviscid contribution) for each boundary. */
  *CD_Mnt,           /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Mnt,           /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Mnt,          /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Mnt,          /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Mnt,          /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Mnt,          /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Mnt,         /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Mnt,         /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Mnt,         /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Mnt,          /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Mnt,          /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Mnt,          /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
  *Surface_CL_Mnt,   /*!< \brief Lift coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CD_Mnt,   /*!< \brief Drag coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CSF_Mnt,  /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CEff_Mnt, /*!< \brief Side-force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFx_Mnt,  /*!< \brief x Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFy_Mnt,  /*!< \brief y Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CFz_Mnt,  /*!< \brief z Force coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMx_Mnt,  /*!< \brief x Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMy_Mnt,  /*!< \brief y Moment coefficient (inviscid contribution) for each monitoring surface. */
  *Surface_CMz_Mnt,  /*!< \brief z Moment coefficient (inviscid contribution) for each monitoring surface. */
  *CEff_Mnt,         /*!< \brief Efficiency (Cl/Cd) (inviscid contribution) for each boundary. */
  *CMerit_Mnt,       /*!< \brief Rotor Figure of Merit (inviscid contribution) for each boundary. */
  *CT_Mnt,           /*!< \brief Thrust coefficient (force in -x direction, inviscid contribution) for each boundary. */
  *CQ_Mnt,           /*!< \brief Torque coefficient (moment in -x direction, inviscid contribution) for each boundary. */
  *CEquivArea_Mnt,   /*!< \brief Equivalent area (inviscid contribution) for each boundary. */
  **CPressure,       /*!< \brief Pressure coefficient for each boundary and vertex. */
  **CPressureTarget, /*!< \brief Target Pressure coefficient for each boundary and vertex. */
  **HeatFlux,        /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **HeatFluxTarget,  /*!< \brief Heat transfer coefficient for each boundary and vertex. */
  **YPlus,           /*!< \brief Yplus for each boundary and vertex. */
  ***CharacPrimVar,  /*!< \brief Value of the characteristic variables at each boundary. */
  ***DonorPrimVar,   /*!< \brief Value of the donor variables at each boundary. */
  *ForceInviscid,    /*!< \brief Inviscid force for each boundary. */
  *MomentInviscid,   /*!< \brief Inviscid moment for each boundary. */
  *ForceMomentum,    /*!< \brief Inviscid force for each boundary. */
  *MomentMomentum;   /*!< \brief Inviscid moment for each boundary. */
  su2double
  *Inflow_MassFlow,  /*!< \brief Mass flow rate for each boundary. */
  *Exhaust_MassFlow, /*!< \brief Mass flow rate for each boundary. */
  *Inflow_Pressure,  /*!< \brief Fan face pressure for each boundary. */
  *Inflow_Mach,      /*!< \brief Fan face mach number for each boundary. */
  *Inflow_Area,      /*!< \brief Boundary total area. */
  *Exhaust_Area,     /*!< \brief Boundary total area. */
  *Exhaust_Pressure,      /*!< \brief Fan face pressure for each boundary. */
  *Exhaust_Temperature,   /*!< \brief Fan face mach number for each boundary. */
  Inflow_MassFlow_Total,  /*!< \brief Mass flow rate for each boundary. */
  Exhaust_MassFlow_Total, /*!< \brief Mass flow rate for each boundary. */
  Inflow_Pressure_Total,  /*!< \brief Fan face pressure for each boundary. */
  Inflow_Mach_Total,      /*!< \brief Fan face mach number for each boundary. */
  InverseDesign;          /*!< \brief Inverse design functional for each boundary. */
  unsigned long
  **DonorGlobalIndex;    /*!< \brief Value of the donor global index. */
  su2double
  **ActDisk_DeltaP,      /*!< \brief Value of the Delta P. */
  **ActDisk_DeltaT;      /*!< \brief Value of the Delta T. */
  su2double
  **Inlet_Ptotal,    /*!< \brief Value of the Total P. */
  **Inlet_Ttotal,    /*!< \brief Value of the Total T. */
  ***Inlet_FlowDir;  /*!< \brief Value of the Flow Direction. */
  
  su2double
  AllBound_CD_Inv,           /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Inv,           /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Inv,          /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Inv,          /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Inv,          /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Inv,          /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Inv,         /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Inv,         /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Inv,         /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Inv,          /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Inv,          /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Inv,          /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Inv,         /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Inv,       /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Inv,           /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Inv,           /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEquivArea_Inv,   /*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CNearFieldOF_Inv; /*!< \brief Near-Field press coefficient (inviscid contribution) for all the boundaries. */
  
  su2double
  AllBound_CD_Mnt,      /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Mnt,      /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Mnt,     /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Mnt,     /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Mnt,     /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Mnt,     /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Mnt,    /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Mnt,    /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Mnt,    /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Mnt,     /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Mnt,     /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Mnt,     /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Mnt,    /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Mnt,  /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Mnt,      /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Mnt;      /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */
  
  su2double
  Total_ComboObj,       /*!< \brief Total 'combo' objective for all monitored boundaries */
  Total_CD,             /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CL,             /*!< \brief Total lift coefficient for all the boundaries. */
  Total_CL_Prev,        /*!< \brief Total lift coefficient for all the boundaries (fixed lift mode). */
  Total_SolidCD,        /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CD_Prev,        /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_NetThrust,      /*!< \brief Total drag coefficient for all the boundaries. */
  Total_Power,          /*!< \brief Total drag coefficient for all the boundaries. */
  Total_ReverseFlow,    /*!< \brief Total drag coefficient for all the boundaries. */
  Total_IDC,            /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDC_Mach,       /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDR,            /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_DC60,           /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_MFR,            /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Prop_Eff,       /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_ByPassProp_Eff, /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Adiab_Eff,      /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Poly_Eff,       /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Custom_ObjFunc, /*!< \brief Total custom objective function for all the boundaries. */
  Total_CSF,      /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CMx,      /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMx_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMy,      /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMy_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMz,      /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CMz_Prev, /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CoPx,     /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CoPy,     /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CoPz,     /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CFx,      /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy,      /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz,      /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CEff,     /*!< \brief Total efficiency coefficient for all the boundaries. */
  Total_CMerit,   /*!< \brief Total rotor Figure of Merit for all the boundaries. */
  Total_CT,       /*!< \brief Total thrust coefficient for all the boundaries. */
  Total_CQ,       /*!< \brief Total torque coefficient for all the boundaries. */
  Total_Heat,     /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat,  /*!< \brief Maximum heat flux on all boundaries. */
  Total_AeroCD,   /*!< \brief Total aero drag coefficient for all the boundaries. */
  Total_CEquivArea,   /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_CNearFieldOF, /*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
  Total_CpDiff,       /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_HeatFluxDiff, /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_MassFlowRate; /*!< \brief Total Mass Flow Rate on monitored boundaries. */
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
  
  su2double
  *SecondaryVar_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
  *SecondaryVar_j;  /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double
  *PrimVar_i,      /*!< \brief Auxiliary vector for storing the solution at point i. */
  *PrimVar_j;      /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double **LowMach_Precontioner; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
  bool space_centered,  /*!< \brief True if space centered scheeme used. */
  euler_implicit,       /*!< \brief True if euler implicit scheme used. */
  least_squares;        /*!< \brief True if computing gradients by least squares. */
  su2double Gamma;                  /*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;        /*!< \brief Fluids's Gamma - 1.0  . */
  
  su2double *Primitive,    /*!< \brief Auxiliary nPrimVar vector. */
  *Primitive_i,            /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Primitive_j;            /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */
  
  su2double *Secondary,    /*!< \brief Auxiliary nPrimVar vector. */
  *Secondary_i,            /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Secondary_j;             /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  su2double AoA_Prev, /*!< \brief Old value of the angle of attack (monitored). */
  AoA_inc;
  bool Start_AoA_FD,  /*!< \brief Boolean for start of finite differencing for FixedCL mode */
  End_AoA_FD,         /*!< \brief Boolean for end of finite differencing for FixedCL mode */
  Update_AoA;         /*!< \brief Boolean to signal Angle of Attack Update */
  unsigned long Iter_Update_AoA; /*!< \brief Iteration at which AoA was updated last */
  su2double dCL_dAlpha;          /*!< \brief Value of dCL_dAlpha used to control CL in fixed CL mode */
  unsigned long BCThrust_Counter;
  unsigned short nSpanWiseSections;  /*!< \brief Number of span-wise sections. */
  unsigned short nSpanMax;           /*!< \brief Max number of maximum span-wise sections for all zones */
  unsigned short nMarkerTurboPerf;   /*!< \brief Number of turbo performance. */

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

  CEulerVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

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
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction = false);
  
  /*!
   * \brief Compute the gradient of the primitive variables using a Least-Squares method,
   *        and stores the result in the <i>Gradient_Primitive</i> variable.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction = false);
  
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
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  void Evaluate_ObjFunc(CConfig *config);
  
  /*!
   * \author: T. Kattmann
   *
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
  void BC_Sym_Plane(CGeometry      *geometry, 
                    CSolver        **solver_container, 
                    CNumerics      *conv_numerics, 
                    CNumerics      *visc_numerics, 
                    CConfig        *config, 
                    unsigned short val_marker) override;
  
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
   * \brief Check for convergence of the Fixed CL mode to the target CL
   * \param[in] config - Definition of the particular problem.
   * \param[in] convergence - boolean for whether the solution is converged
   * \return boolean for whether the Fixed CL mode is converged to target CL
   */
  bool FixedCL_Convergence(CConfig *config, bool convergence);

  /*!
   * \brief Checking whether fixed CL mode in finite-differencing mode 
   * \return boolean for whether the Fixed CL mode is currently in finite-differencing mode
   */
  bool GetStart_AoA_FD(void);

  /*!
   * \brief Checking whether fixed CL mode in finite-differencing mode 
   * \return boolean for whether the Fixed CL mode is currently in finite-differencing mode
   */
  bool GetEnd_AoA_FD(void);

    /*!
   * \brief Get the iteration of the last AoA update (Fixed CL Mode)
   * \return value for the last iteration that the AoA was updated
   */
  unsigned long GetIter_Update_AoA();

  /*!
   * \brief Get the AoA before the most recent update
   * \return value of the AoA before most recent update
   */
  su2double GetPrevious_AoA();

  /*!
   * \brief Get the CL Driver's control command
   * \return value of CL Driver control command (AoA_inc)
   */
  su2double GetAoA_inc();

  /*!
   * \brief Set gradients of coefficients for fixed CL mode
   * \param[in] config - Definition of the particular problem.
   */
  void SetCoefficient_Gradients(CConfig *config);
  
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
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over a nonlinear iteration for stability.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeUnderRelaxationFactor(CSolver **solver, CConfig *config);
  
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
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter);
  
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
  void MixedOut_Average (CConfig *config, su2double val_init_pressure, const su2double *val_Averaged_Flux,
                         const su2double *val_normal, su2double& pressure_mix, su2double& density_mix);

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
  void ComputeTurboVelocity(const su2double *cartesianVelocity, const su2double *turboNormal, su2double *turboVelocity,
                            unsigned short marker_flag, unsigned short marker_kindturb);

  /*!
   * \brief it take a velocity in the cartesian reference of framework and transform into the turbomachinery frame of reference.
   * \param[in] cartesianVelocity - cartesian components of velocity vector.
   * \param[in] turboNormal - normal vector in the turbomachinery frame of reference.
   * \param[in] turboVelocity - velocity vector in the turbomachinery frame of reference.
   */
  void ComputeBackVelocity(const su2double *turboVelocity, const su2double *turboNormal, su2double *cartesianVelocity,
                           unsigned short marker_flag, unsigned short marker_kindturb);

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

#include "../solver_inlines/CEulerSolver.inl"