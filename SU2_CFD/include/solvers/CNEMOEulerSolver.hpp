/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
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
#include "../variables/CNEMOEulerVariable.hpp"
#include "../include/fluid/CNEMOGas.hpp"

/*!
 * \class CNEMOEulerSolver
 * \brief Main class for defining the NEMO Euler's flow solver.
 * \ingroup Euler_Equations
 * \author S. R. Copeland, F. Palacios, W. Maier.
 * \version 6.1.0
 */
class CNEMOEulerSolver : public CSolver {
protected:

  unsigned short
  nSpecies;	                             /*!< \brief Number of species in the gas mixture. */
                  
  su2double                  
  Gamma,                                 /*!< \brief Mixture Cp/Cv. */
  Gamma_Minus_One;                       /*!< \brief Mixture Cp/Cv - 1. */
                  
  su2double                  
  Mach_Inf,                              /*!< \brief Free stream Mach number. */
  *Density,                              /*!< \brief Free stream species density. */
  Density_Inf,                           /*!< \brief Free stream sdensity. */
  Energy_ve_Inf,                         /*!< \brief Vib.-el. free stream energy. */
  Pressure_Inf,	                         /*!< \brief Free stream pressure. */
  *Velocity_Inf,	                       /*!< \brief Free stream flow velocity. */
  Temperature_Inf,                       /*!< \brief Trans.-rot. free stream temperature. */
  Temperature_ve_Inf;                    /*!< \brief Vib.-el. free stream temperature. */
  const su2double *MassFrac_Inf;        /*!< \brief Free stream species mass fraction. */
  

  su2double
  *lowerlimit,         /*!< \brief contains lower limits for conserved variables. */
  *upperlimit;         /*!< \brief contains upper limits for conserved variables. */

  su2double
  *CD_Inv,           /*!< \brief Drag coefficient (inviscid contribution) for each boundary. */
  *CL_Inv,           /*!< \brief Lift coefficient (inviscid contribution) for each boundary. */
  *CSF_Inv,          /*!< \brief Sideforce coefficient (inviscid contribution) for each boundary. */
  *CMx_Inv,          /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CMy_Inv,          /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CMz_Inv,          /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CoPx_Inv,         /*!< \brief x Moment coefficient (inviscid contribution) for each boundary. */
  *CoPy_Inv,         /*!< \brief y Moment coefficient (inviscid contribution) for each boundary. */
  *CoPz_Inv,         /*!< \brief z Moment coefficient (inviscid contribution) for each boundary. */
  *CFx_Inv,          /*!< \brief x Force coefficient (inviscid contribution) for each boundary. */
  *CFy_Inv,          /*!< \brief y Force coefficient (inviscid contribution) for each boundary. */
  *CFz_Inv,          /*!< \brief z Force coefficient (inviscid contribution) for each boundary. */
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
  *Inflow_MassFlow,       /*!< \brief Mass flow rate for each boundary. */
  *Exhaust_MassFlow,      /*!< \brief Mass flow rate for each boundary. */
  *Inflow_Pressure,       /*!< \brief Fan face pressure for each boundary. */
  *Inflow_Mach,           /*!< \brief Fan face mach number for each boundary. */
  *Inflow_Area,           /*!< \brief Boundary total area. */
  *Exhaust_Area,          /*!< \brief Boundary total area. */
  *Exhaust_Pressure,      /*!< \brief Fan face pressure for each boundary. */
  *Exhaust_Temperature,   /*!< \brief Fan face mach number for each boundary. */
  Inflow_MassFlow_Total,  /*!< \brief Mass flow rate for each boundary. */
  Exhaust_MassFlow_Total, /*!< \brief Mass flow rate for each boundary. */
  Inflow_Pressure_Total,  /*!< \brief Fan face pressure for each boundary. */
  Inflow_Mach_Total,      /*!< \brief Fan face mach number for each boundary. */
  InverseDesign;          /*!< \brief Inverse design functional for each boundary. */

  su2double *Source;   /*!< \brief Auxiliary vector to store source terms. */

  unsigned long **DonorGlobalIndex; /*!< \brief Value of the donor global index. */
  su2double **ActDisk_DeltaP,       /*!< \brief Value of the Delta P. */
  **ActDisk_DeltaT;                 /*!< \brief Value of the Delta T. */
  su2double **Inlet_Ptotal,         /*!< \brief Value of the Total P. */
  **Inlet_Ttotal,                   /*!< \brief Value of the Total T. */
  ***Inlet_FlowDir;                 /*!< \brief Value of the Flow Direction. */

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
  AllBound_CD_Mnt,     /*!< \brief Total drag coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CL_Mnt,     /*!< \brief Total lift coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CSF_Mnt,    /*!< \brief Total sideforce coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMx_Mnt,    /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Mnt,    /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Mnt,    /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Mnt,   /*!< \brief Total x moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Mnt,   /*!< \brief Total y moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Mnt,   /*!< \brief Total z moment coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFx_Mnt,    /*!< \brief Total x force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Mnt,    /*!< \brief Total y force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Mnt,    /*!< \brief Total z force coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Mnt,   /*!< \brief Efficient coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Mnt, /*!< \brief Rotor Figure of Merit (inviscid contribution) for all the boundaries. */
  AllBound_CT_Mnt,     /*!< \brief Total thrust coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CQ_Mnt;     /*!< \brief Total torque coefficient (inviscid contribution) for all the boundaries. */

  su2double
  Total_ComboObj,       /*!< \brief Total 'combo' objective for all monitored boundaries */
  AoA_Prev,             /*!< \brief Old value of the AoA for fixed lift mode. */
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
  Total_CSF,            /*!< \brief Total sideforce coefficient for all the boundaries. */
  Total_CMx,            /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CMx_Prev,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMy,            /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CMy_Prev,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMz,            /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CMz_Prev,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CoPx,           /*!< \brief Total x moment coefficient for all the boundaries. */
  Total_CoPy,           /*!< \brief Total y moment coefficient for all the boundaries. */
  Total_CoPz,           /*!< \brief Total z moment coefficient for all the boundaries. */
  Total_CFx,            /*!< \brief Total x force coefficient for all the boundaries. */
  Total_CFy,            /*!< \brief Total y force coefficient for all the boundaries. */
  Total_CFz,            /*!< \brief Total z force coefficient for all the boundaries. */
  Total_CEff,           /*!< \brief Total efficiency coefficient for all the boundaries. */
  Total_CMerit,         /*!< \brief Total rotor Figure of Merit for all the boundaries. */
  Total_CT,             /*!< \brief Total thrust coefficient for all the boundaries. */
  Total_CQ,             /*!< \brief Total torque coefficient for all the boundaries. */
  Total_Heat,           /*!< \brief Total heat load for all the boundaries. */
  Total_MaxHeat,        /*!< \brief Maximum heat flux on all boundaries. */
  Total_AeroCD,         /*!< \brief Total aero drag coefficient for all the boundaries. */
  Total_CEquivArea,     /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_CNearFieldOF,   /*!< \brief Total Near-Field Pressure coefficient for all the boundaries. */
  Total_CpDiff,         /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_HeatFluxDiff,   /*!< \brief Total Equivalent Area coefficient for all the boundaries. */
  Total_MassFlowRate;   /*!< \brief Total Mass Flow Rate on monitored boundaries. */
  su2double
  *Surface_CL,         /*!< \brief Lift coefficient for each monitoring surface. */
  *Surface_CD,         /*!< \brief Drag coefficient for each monitoring surface. */
  *Surface_CSF,        /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CEff,       /*!< \brief Side-force coefficient for each monitoring surface. */
  *Surface_CFx,        /*!< \brief x Force coefficient for each monitoring surface. */
  *Surface_CFy,        /*!< \brief y Force coefficient for each monitoring surface. */
  *Surface_CFz,        /*!< \brief z Force coefficient for each monitoring surface. */
  *Surface_CMx,        /*!< \brief x Moment coefficient for each monitoring surface. */
  *Surface_CMy,        /*!< \brief y Moment coefficient for each monitoring surface. */
  *Surface_CMz,        /*!< \brief z Moment coefficient for each monitoring surface. */
  *Surface_HF_Visc,    /*!< \brief Total (integrated) heat flux for each monitored surface. */
  *Surface_MaxHF_Visc; /*!< \brief Maximum heat flux for each monitored surface. */

  su2double *iPoint_UndLapl, /*!< \brief Auxiliary variable for the undivided Laplacians. */
  *jPoint_UndLapl;           /*!< \brief Auxiliary variable for the undivided Laplacians. */
  su2double *SecondaryVar_i, /*!< \brief Auxiliary vector for storing the solution at point i. */
  *SecondaryVar_j;           /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double *PrimVar_i,      /*!< \brief Auxiliary vector for storing the solution at point i. */
  *PrimVar_j;                /*!< \brief Auxiliary vector for storing the solution at point j. */
  su2double
  **LowMach_Precontioner; /*!< \brief Auxiliary vector for storing the inverse of Roe-turkel preconditioner. */
  unsigned long nMarker,  /*!< \brief Total number of markers using the grid information. */
  *nVertex;               /*!< \brief Store nVertex at each marker for deallocation */
  bool space_centered,    /*!< \brief True if space centered scheeme used. */
  euler_implicit,         /*!< \brief True if euler implicit scheme used. */
  least_squares;          /*!< \brief True if computing gradients by least squares. */

  su2double *Primitive,   /*!< \brief Auxiliary nPrimVar vector. */
  *Primitive_i,           /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Primitive_j;           /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  su2double *Secondary,   /*!< \brief Auxiliary nPrimVar vector. */
  *Secondary_i,           /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
  *Secondary_j;           /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

  su2double Cauchy_Value,        /*!< \brief Summed value of the convergence indicator. */
  Cauchy_Func;                   /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter; /*!< \brief Number of elements of the Cauchy serial. */
  su2double *Cauchy_Serie;       /*!< \brief Complete Cauchy serial. */
  su2double Old_Func,            /*!< \brief Old value of the objective function (the function which is monitored). */
  New_Func;                      /*!< \brief Current value of the objective function (the function which is monitored). */
  su2double AoA_old;             /*!< \brief Old value of the angle of attack (monitored). */
  unsigned long AoA_Counter;
  bool AoA_FD_Change;
  unsigned long BCThrust_Counter;
  unsigned short nSpanWiseSections; /*!< \brief Number of span-wise sections. */
  unsigned short nSpanMax;          /*!< \brief Max number of maximum span-wise sections for all zones */
  unsigned short nMarkerTurboPerf;  /*!< \brief Number of turbo performance. */

  CNEMOGas  *FluidModel;         /*!< \brief fluid model used in the solver */

  

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

  CNEMOEulerVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */
  CNEMOEulerVariable* node_infty;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }


public:

  /*!
     * \brief Constructor of the class.
     */
  CNEMOEulerSolver(void);

  /*!
     * \overload
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
  CNEMOEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
     * \brief Destructor of the class.
     */
  virtual ~CNEMOEulerSolver(void);

// /*!
//    * \brief Compute slope limiter.
//    * \param[in] geometry - Geometrical definition of the problem.
//    * \param[in] config - Definition of the particular problem.
//    */
// void SetSolution_Limiter(CGeometry *geometry, CConfig *config);

  /*!
     * \brief Set the maximum value of the eigenvalue.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
  void SetMax_Eigenvalue(CGeometry *geometry, CConfig *config);

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
                    unsigned long Iteration);

  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter);

  /*!
   * \brief Reset Node Infty for discrete adjoint
   */
  void ResetNodeInfty(su2double pressure_inf, const su2double *massfrac_inf, su2double *mvec_inf, su2double temperature_inf,
                      su2double temperature_ve_inf, CConfig *config);

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
   * \brief Compute the spatial integration using a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                         CConfig *config, unsigned short iMesh, unsigned short iRKStep);

  /*!
     * \brief Compute the spatial integration using a upwind scheme.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] solver - Description of the numerical method.
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
     * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
     * \param[in] config - Definition of the particular problem.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
     */
   void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) final;

  /*!
     * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
     * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
     */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief Compute the Green-Gauss gradient of the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  //void SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction = false) final;

  /*!
   * \brief Compute the Least Squares gradient of the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
 // void SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction = false) final; 

  /*!
   * \brief Compute slope limiter.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  //void SetSolution_Limiter(CGeometry *geometry, CConfig *config) final;

  /*!
     * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
     * \param[in] iPoint - Index of the grid point
     * \param[in] config - Definition of the particular problem.
     */
  void SetPreconditioner(CConfig *config, unsigned short iPoint);

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh) final;

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   */
  void Evaluate_ObjFunc(CConfig *config) override;

  /*!
     * \brief Impose via the residual the Euler wall boundary condition.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] numerics - Description of the numerical method.
     * \param[in] config - Definition of the particular problem.
     * \param[in] val_marker - Surface marker where the boundary condition is applied.
     */
  void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

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
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

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
  void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) override;

  /*!
     * \brief Update the solution using an explicit Euler scheme.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Runge-Kutta step.
     */
  void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                            CConfig *config, unsigned short iRKStep) override;

  /*!
     * \brief Update the solution using an implicit Euler scheme.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
     */
  void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) override;

  /*!
     * \brief Compute the pressure forces and all the adimensional coefficients.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
  void Pressure_Forces(CGeometry *geometry, CConfig *config) override;

  /*!
     * \brief Compute the Momentum forces and all the adimensional coefficients.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
  void Momentum_Forces(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Provide the non dimensional lift coefficient (inviscid contribution).
   * \param val_marker Surface where the coefficient is going to be computed.
   * \return Value of the lift coefficient (inviscid contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Inv(unsigned short val_marker) const final { return CL_Inv[val_marker]; }

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
  inline su2double GetTotal_HeatFlux() const final { return Total_Heat; }

  /*!
   * \brief Provide the total heat load.
   * \return Value of the heat load (viscous contribution).
   */
  inline su2double GetTotal_MaxHeatFlux() const final { return Total_MaxHeat; }

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
   * \brief Store the total (inviscid + viscous) non dimensional lift coefficient.
   * \param[in] val_Total_CL - Value of the total lift coefficient.
   */
  inline void SetTotal_CL(su2double val_Total_CL) final { Total_CL = val_Total_CL; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_CD(su2double val_Total_CD) final { Total_CD = val_Total_CD; }

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
   * \brief Provide the Pressure coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline su2double GetCPressure(unsigned short val_marker, unsigned long val_vertex) const final {
    return CPressure[val_marker][val_vertex];
  }
     
  //su2double GetCPressure(unsigned short val_marker, unsigned long val_vertex) override;

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
                            unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) override;

  /*!
     * \brief Load a direct flow solution for use with the adjoint solver.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] val_iZone - Current zone in the mesh.
     */
  void GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone);

  /*!
     * \brief Load the output data container with the variables to be written to the volume solution file.
   * \param[in] config - Definition of the particular problem.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] data_container - Container holding the output variable data.
   * \param[in] nOutput_Vars - Number of output variables being stored.
     */
  void SetVolume_Output(CConfig *config, CGeometry *geometry, su2double **data_container, unsigned short nOutput_Vars);

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
   * \brief Set the old solution variables to the current solution value for Runge-Kutta iteration.
            It is a virtual function, because for the DG-FEM solver a different version is needed.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  //inline void Set_OldSolution(CGeometry *geometry) final { nodes->Set_OldSolution(); }

};