/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
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
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

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
   * \brief Compute the density multiply by energy at the infinity.
   * \return Value of the density multiply by  energy at the infinity.
   */
  inline su2double GetDensity_Energy_Inf(void) const final { return Density_Inf*Energy_Inf; }

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline su2double GetPressure_Inf(void) const final { return Pressure_Inf; }

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
  inline su2double *GetVelocity_Inf(void) const final { return Velocity_Inf; }

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
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Centered_Residual(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *numerics,
                         CConfig *config,
                         unsigned short iMesh,
                         unsigned short iRKStep) final;

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) final;

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
  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CNumerics *second_numerics,
                       CConfig *config, unsigned short iMesh) override;

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
                                       bool Output) override;

  /*!
   * \brief Compute a pressure sensor switch.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Compute Ducros Sensor for Roe Dissipation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config) final;

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
   * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
   * \param[in] iPoint - Index of the grid point
   * \param[in] config - Definition of the particular problem.
   */
  void SetPreconditioner(CConfig *config, unsigned long iPoint) final;

  /*!
   * \brief Compute the undivided laplacian for the solution, except the energy equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Compute the max eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) final;

  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) final;

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  void Evaluate_ObjFunc(CConfig *config) override;

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
   * \brief Impose the interface boundary condition using the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Interface_Boundary(CGeometry *geometry,
                             CSolver **solver_container,
                             CNumerics *numerics,
                             CConfig *config,
                             unsigned short val_marker) final;

  /*!
   * \brief Impose the near-field boundary condition using the residual.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_NearField_Boundary(CGeometry *geometry,
                             CSolver **solver_container,
                             CNumerics *numerics,
                             CConfig *config,
                             unsigned short val_marker) final;

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
   * \brief Set the new solution variables to the current solution value for classical RK.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  inline void Set_NewSolution(CGeometry *geometry) final { nodes->SetSolution_New(); }

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
   * \brief Compute the Fan face Mach number.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  void GetPower_Properties(CGeometry *geometry,
                           CConfig *config,
                           unsigned short iMesh,
                           bool Output) final;

  /*!
   * \brief Update the AoA and freestream velocity at the farfield.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void SetActDisk_BCThrust(CGeometry *geometry,
                           CSolver **solver_container,
                           CConfig *config,
                           unsigned short iMesh,
                           bool Output) final;

  /*!
   * \brief Update the AoA and freestream velocity at the farfield.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void SetFarfield_AoA(CGeometry *geometry,
                       CSolver **solver_container,
                       CConfig *config,
                       unsigned short iMesh,
                       bool Output) final;

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
  void ExplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) final;

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
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Inv(unsigned short val_marker) const final { return CL_Inv[val_marker]; }

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
  inline su2double GetTotal_CEquivArea() const final { return Total_CEquivArea; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional aero CD.
   * \return Value of the Aero CD coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_AeroCD() const final { return Total_AeroCD; }

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
   * \brief Provide the total (inviscid + viscous) non dimensional Near-Field pressure coefficient.
   * \return Value of the NearField pressure coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CNearFieldOF() const final { return Total_CNearFieldOF; }

  /*!
   * \author H. Kline
   * \brief Add to the value of the total 'combo' objective.
   * \param[in] val_obj - Value of the contribution to the 'combo' objective.
   */
  inline void AddTotal_ComboObj(su2double val_obj) final {Total_ComboObj +=val_obj;}

  /*!
   * \brief Set the value of the Equivalent Area coefficient.
   * \param[in] val_cequivarea - Value of the Equivalent Area coefficient.
   */
  inline void SetTotal_CEquivArea(su2double val_cequivarea) final { Total_CEquivArea = val_cequivarea; }

  /*!
   * \brief Set the value of the Aero drag.
   * \param[in] val_cequivarea - Value of the aero drag.
   */
  inline void SetTotal_AeroCD(su2double val_aerocd) final { Total_AeroCD = val_aerocd; }

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
   * \brief Set the value of the Near-Field pressure oefficient.
   * \param[in] val_cnearfieldpress - Value of the Near-Field pressure coefficient.
   */
  inline void SetTotal_CNearFieldOF(su2double val_cnearfieldpress) final { Total_CNearFieldOF = val_cnearfieldpress; }

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
   * \brief Provide the total (inviscid + viscous) non dimensional lift coefficient.
   * \return Value of the lift coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CL() const final { return Total_CL; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CD() const final { return Total_CD; }

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
   * \brief Provide the total custom objective function.
   * \return Value of the custom objective function.
   */
  inline su2double GetTotal_Custom_ObjFunc() const final { return Total_Custom_ObjFunc; }

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
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_CD(su2double val_Total_CD) final { Total_CD = val_Total_CD; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  inline void SetTotal_CL(su2double val_Total_CL) final { Total_CL = val_Total_CL; }

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
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetCharacPrimVar(unsigned short val_marker,
                               unsigned long val_vertex,
                               unsigned short val_var,
                               su2double val_value) final {
    CharacPrimVar[val_marker][val_vertex][val_var] = val_value;
  }

  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double *GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex) const final{
    return DonorPrimVar[val_marker][val_vertex];
  }

  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetDonorPrimVar(unsigned short val_marker,
                              unsigned long val_vertex,
                              unsigned short val_var,
                              su2double val_value) final {
    DonorPrimVar[val_marker][val_vertex][val_var] = val_value;
  }

  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double GetDonorPrimVar(unsigned short val_marker,
                                   unsigned long val_vertex,
                                   unsigned short val_var) const final {
    return DonorPrimVar[val_marker][val_vertex][val_var];
  }

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
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double GetActDisk_DeltaP(unsigned short val_marker,
                                     unsigned long val_vertex) const final {
    return ActDisk_DeltaP[val_marker][val_vertex];
  }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetActDisk_DeltaP(unsigned short val_marker,
                                unsigned long val_vertex,
                                su2double val_deltap) final { ActDisk_DeltaP[val_marker][val_vertex] = val_deltap; }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex) final {
    return ActDisk_DeltaT[val_marker][val_vertex];
  }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetActDisk_DeltaT(unsigned short val_marker,
                                unsigned long val_vertex,
                                su2double val_deltat) final {
    ActDisk_DeltaT[val_marker][val_vertex] = val_deltat;
  }

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
   * \brief Set the value of the total temperature at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total temperature is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total temperature is set.
   * \param[in] val_ttotal - Value of the total temperature
   */
  inline void SetInlet_Ttotal(unsigned short val_marker,
                              unsigned long val_vertex,
                              su2double val_ttotal) final {
    /*--- Since this call can be accessed indirectly using python, do some error
     * checking to prevent segmentation faults ---*/
    if (val_marker >= nMarker)
      SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
    else if (Inlet_Ttotal == NULL || Inlet_Ttotal[val_marker] == NULL)
      SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
    else if (val_vertex >= nVertex[val_marker])
      SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
    else
      Inlet_Ttotal[val_marker][val_vertex] = val_ttotal;
  }


  /*!
   * \brief Set the value of the total pressure at an inlet boundary.
   * \param[in] val_marker - Surface marker where the total pressure is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the total pressure is set.
   * \param[in] val_ptotal - Value of the total pressure
   */
  inline void SetInlet_Ptotal(unsigned short val_marker,
                              unsigned long val_vertex,
                              su2double val_ptotal) final {
    /*--- Since this call can be accessed indirectly using python, do some error
     * checking to prevent segmentation faults ---*/
    if (val_marker >= nMarker)
      SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
    else if (Inlet_Ptotal == NULL || Inlet_Ptotal[val_marker] == NULL)
      SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
    else if (val_vertex >= nVertex[val_marker])
      SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
    else
      Inlet_Ptotal[val_marker][val_vertex] = val_ptotal;
  }


  /*!
   * \brief Set a component of the unit vector representing the flow direction at an inlet boundary.
   * \param[in] val_marker - Surface marker where the flow direction is set.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the flow direction is set.
   * \param[in] val_dim - The component of the flow direction unit vector to be set
   * \param[in] val_flowdir - Component of a unit vector representing the flow direction.
   */
  inline void SetInlet_FlowDir(unsigned short val_marker,
                               unsigned long val_vertex,
                               unsigned short val_dim,
                               su2double val_flowdir) final {
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
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config) final;

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
  inline int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex) const final {
    return SlidingStateNodes[val_marker][val_vertex];
  }

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
  inline void SetPressure_Inf(su2double p_inf) final {Pressure_Inf = p_inf;}

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetTemperature_Inf(su2double t_inf) final {Temperature_Inf = t_inf;}

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(CConfig *config) final;

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
  inline su2double* GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan) const final {
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
  inline su2double* GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan) const final {
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
  inline su2double* GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan) const final {
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
  inline void SetTurboVelocityIn(su2double *value,
                                 unsigned short inMarkerTP,
                                 unsigned short valSpan) final {
    unsigned short iDim;

    for(iDim = 0; iDim < nDim; iDim++)
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
  inline void SetTurboVelocityOut(su2double *value,
                                  unsigned short inMarkerTP,
                                  unsigned short valSpan) final {
    unsigned short iDim;

    for(iDim = 0; iDim < nDim; iDim++)
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
   * \brief Compute the global error measures (L2, Linf) for verification cases.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config   - Definition of the particular problem.
   */
  void ComputeVerificationError(CGeometry *geometry, CConfig *config) final;

};
