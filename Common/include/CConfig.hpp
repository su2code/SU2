/*!
 * \file CConfig.hpp
 * \brief All the information about the definition of the physical problem.
 *        The subroutines and functions are in the <i>CConfig.cpp</i> file.
 * \author F. Palacios, T. Economon, B. Tracey
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

#include "parallelization/mpi_structure.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <array>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <assert.h>

#include "option_structure.hpp"
#include "containers/container_decorators.hpp"

#ifdef HAVE_CGNS
#include "cgnslib.h"
#endif

using namespace std;

/*!
 * \class CConfig
 * \brief Main class for defining the problem; basically this class reads the configuration file, and
 *        stores all the information.
 * \author F. Palacios
 */

class CConfig {
private:
  SU2_MPI::Comm SU2_Communicator; /*!< \brief MPI communicator of SU2.*/
  int rank, size;                 /*!< \brief MPI rank and size.*/
  bool base_config;
  SU2_COMPONENT Kind_SU2;        /*!< \brief Kind of SU2 software component.*/
  unsigned short Ref_NonDim;      /*!< \brief Kind of non dimensionalization.*/
  unsigned short Ref_Inc_NonDim;  /*!< \brief Kind of non dimensionalization.*/
  unsigned short Kind_AverageProcess;            /*!< \brief Kind of mixing process.*/
  unsigned short Kind_PerformanceAverageProcess; /*!< \brief Kind of mixing process.*/
  unsigned short Kind_MixingPlaneInterface;      /*!< \brief Kind of mixing process.*/
  unsigned short Kind_SpanWise;                  /*!< \brief Kind of span-wise section computation.*/
  unsigned short *Kind_TurboMachinery;           /*!< \brief Kind of turbomachynery architecture.*/
  unsigned short iZone, nZone;    /*!< \brief Number of zones in the mesh. */
  unsigned short nZoneSpecified;  /*!< \brief Number of zones that are specified in config file. */
  su2double Highlite_Area;        /*!< \brief Highlite area. */
  su2double Fan_Poly_Eff;         /*!< \brief Fan polytropic effeciency. */
  su2double MinLogResidual;       /*!< \brief Minimum value of the log residual. */
  su2double EA_ScaleFactor;       /*!< \brief Equivalent Area scaling factor */
  su2double AdjointLimit;         /*!< \brief Adjoint variable limit */
  string* ConvField;              /*!< \brief Field used for convergence check.*/
  string FluidName;              /*!< \brief name of the applied fluid. */

  string* WndConvField;              /*!< \brief Function where to apply the windowed convergence criteria for the time average of the unsteady (single zone) flow problem. */
  unsigned short nConvField;         /*!< \brief Number of fields used to monitor convergence.*/
  unsigned short nWndConvField;      /*!< \brief Number of fields used to monitor time convergence.*/
  unsigned short Wnd_Cauchy_Elems;   /*!< \brief Number of elements to evaluate in the time iteration  for convergence of the time average of the unsteady (single zone)´ flow problem.  */
  su2double Wnd_Cauchy_Eps;          /*!< \brief Epsilon used for the convergence of the time average of the unsteady (single zone)´ flow problem. */
  unsigned long Wnd_StartConv_Iter;  /*!< \brief Start convergence criteria at this iteration after Start_Iter_Wnd. */
  bool Wnd_Cauchy_Crit;              /*!< \brief True => Cauchy criterion is used for time average objective function in unsteady flows. */

  bool MG_AdjointFlow;              /*!< \brief MG with the adjoint flow problem */
  su2double *PressureLimits,
  *DensityLimits,
  *TemperatureLimits;             /*!< \brief Limits for the primitive variables */
  bool ActDisk_DoubleSurface;     /*!< \brief actuator disk double surface  */
  bool Engine_HalfModel;          /*!< \brief only half model is in the computational grid  */
  bool ActDisk_SU2_DEF;           /*!< \brief actuator disk double surface  */
  unsigned short nFFD_Iter;       /*!< \brief Iteration for the point inversion problem. */
  unsigned short FFD_Blending;    /*!< \brief Kind of FFD Blending function. */
  su2double FFD_Tol;              /*!< \brief Tolerance in the point inversion problem. */
  bool FFD_IntPrev;                       /*!< \brief Enables self-intersection prevention procedure within the FFD box. */
  unsigned short FFD_IntPrev_MaxIter;     /*!< \brief Amount of iterations for FFD box self-intersection prevention procedure. */
  unsigned short FFD_IntPrev_MaxDepth;    /*!< \brief Maximum recursion depth for FFD box self-intersection procedure. */
  bool ConvexityCheck;                    /*!< \brief Enables convexity check on all mesh elements. */
  unsigned short ConvexityCheck_MaxIter;  /*!< \brief Amount of iterations for convexity check in deformations. */
  unsigned short ConvexityCheck_MaxDepth; /*!< \brief Maximum recursion depth for convexity check in deformations.*/
  su2double Opt_RelaxFactor;              /*!< \brief Scale factor for the line search. */
  su2double Opt_LineSearch_Bound;         /*!< \brief Bounds for the line search. */
  su2double StartTime;
  unsigned short SmoothNumGrid;           /*!< \brief Smooth the numerical grid. */
  bool ContinuousAdjoint,   /*!< \brief Flag to know if the code is solving an adjoint problem. */
  Viscous,                  /*!< \brief Flag to know if the code is solving a viscous problem. */
  EquivArea,                /*!< \brief Flag to know if the code is going to compute and plot the equivalent area. */
  Engine,                   /*!< \brief Flag to know if the code is going to compute a problem with engine. */
  InvDesign_Cp,             /*!< \brief Flag to know if the code is going to compute and plot the inverse design. */
  InvDesign_HeatFlux,       /*!< \brief Flag to know if the code is going to compute and plot the inverse design. */
  Wind_Gust,                /*!< \brief Flag to know if there is a wind gust. */
  Turb_Fixed_Values,        /*!< \brief Flag to know if there are fixed values for turbulence quantities in one half-plane. */
  Aeroelastic_Simulation,   /*!< \brief Flag to know if there is an aeroelastic simulation. */
  Weakly_Coupled_Heat,      /*!< \brief Flag to know if a heat equation should be weakly coupled to the incompressible solver. */
  Rotating_Frame,           /*!< \brief Flag to know if there is a rotating frame. */
  PoissonSolver,            /*!< \brief Flag to know if we are solving  poisson forces  in plasma solver. */
  Low_Mach_Precon,          /*!< \brief Flag to know if we are using a low Mach number preconditioner. */
  Low_Mach_Corr,            /*!< \brief Flag to know if we are using a low Mach number correction. */
  GravityForce,             /*!< \brief Flag to know if the gravity force is incuded in the formulation. */
  VorticityConfinement,     /*!< \brief Flag to know if the Vorticity Confinement is included in the formulation. */
  SubsonicEngine,           /*!< \brief Engine intake subsonic region. */
  Frozen_Visc_Cont,         /*!< \brief Flag for cont. adjoint problem with/without frozen viscosity. */
  Frozen_Visc_Disc,         /*!< \brief Flag for disc. adjoint problem with/without frozen viscosity. */
  Frozen_Limiter_Disc,      /*!< \brief Flag for disc. adjoint problem with/without frozen limiter. */
  Inconsistent_Disc,        /*!< \brief Use an inconsistent (primal/dual) discrete adjoint formulation. */
  Sens_Remove_Sharp,        /*!< \brief Flag for removing or not the sharp edges from the sensitivity computation. */
  Hold_GridFixed,           /*!< \brief Flag hold fixed some part of the mesh during the deformation. */
  Axisymmetric,             /*!< \brief Flag for axisymmetric calculations */
  Integrated_HeatFlux;      /*!< \brief Flag for heat flux BC whether it deals with integrated values.*/
  su2double Buffet_k;       /*!< \brief Sharpness coefficient for buffet sensor.*/
  su2double Buffet_lambda;  /*!< \brief Offset parameter for buffet sensor.*/
  su2double Damp_Engine_Inflow;   /*!< \brief Damping factor for the engine inlet. */
  su2double Damp_Engine_Exhaust;  /*!< \brief Damping factor for the engine exhaust. */
  su2double Damp_Res_Restric,     /*!< \brief Damping factor for the residual restriction. */
  Damp_Correc_Prolong;            /*!< \brief Damping factor for the correction prolongation. */
  su2double Position_Plane;    /*!< \brief Position of the Near-Field (y coordinate 2D, and z coordinate 3D). */
  su2double WeightCd;          /*!< \brief Weight of the drag coefficient. */
  su2double dCD_dCL;           /*!< \brief Fixed Cl mode derivate . */
  su2double dCMx_dCL;          /*!< \brief Fixed Cl mode derivate. */
  su2double dCMy_dCL;          /*!< \brief Fixed Cl mode derivate. */
  su2double dCMz_dCL;          /*!< \brief Fixed Cl mode derivate. */
  su2double CL_Target;         /*!< \brief Fixed Cl mode Target Cl. */
  su2double Confinement_Param; /*!< \brief Confinement paramenter for Vorticity Confinement method. */
  TIME_MARCHING TimeMarching;        /*!< \brief Steady or unsteady (time stepping or dual time stepping) computation. */
  su2double FixAzimuthalLine;        /*!< \brief Fix an azimuthal line due to misalignments of the nearfield. */
  su2double **DV_Value;              /*!< \brief Previous value of the design variable. */
  su2double Venkat_LimiterCoeff;     /*!< \brief Limiter coefficient */
  unsigned long LimiterIter;         /*!< \brief Freeze the value of the limiter after a number of iterations */
  su2double AdjSharp_LimiterCoeff;   /*!< \brief Coefficient to identify the limit of a sharp edge. */
  unsigned short SystemMeasurements; /*!< \brief System of measurements. */
  ENUM_REGIME Kind_Regime;           /*!< \brief Kind of flow regime: in/compressible. */
  unsigned short *Kind_ObjFunc;      /*!< \brief Kind of objective function. */
  su2double *Weight_ObjFunc;         /*!< \brief Weight applied to objective function. */
  unsigned short Kind_SensSmooth;    /*!< \brief Kind of sensitivity smoothing technique. */
  unsigned short Continuous_Eqns;    /*!< \brief Which equations to treat continuously (Hybrid adjoint)*/
  unsigned short Discrete_Eqns;      /*!< \brief Which equations to treat discretely (Hybrid adjoint). */
  unsigned short *Design_Variable;   /*!< \brief Kind of design variable. */
  unsigned short nTimeInstances;     /*!< \brief Number of periodic time instances for  harmonic balance. */
  su2double HarmonicBalance_Period;  /*!< \brief Period of oscillation to be used with harmonic balance computations. */
  su2double Delta_UnstTime,          /*!< \brief Time step for unsteady computations. */
  Delta_UnstTimeND;                  /*!< \brief Time step for unsteady computations (non dimensional). */
  su2double Total_UnstTime,       /*!< \brief Total time for unsteady computations. */
  Total_UnstTimeND;               /*!< \brief Total time for unsteady computations (non dimensional). */
  su2double Current_UnstTime,     /*!< \brief Global time of the unsteady simulation. */
  Current_UnstTimeND;             /*!< \brief Global time of the unsteady simulation. */
  unsigned short nMarker_Euler,   /*!< \brief Number of Euler wall markers. */
  nMarker_FarField,               /*!< \brief Number of far-field markers. */
  nMarker_Custom,                 /*!< \brief Number of custom markers. */
  nMarker_SymWall,                /*!< \brief Number of symmetry wall markers. */
  nMarker_PerBound,               /*!< \brief Number of periodic boundary markers. */
  nMarker_MixingPlaneInterface,   /*!< \brief Number of mixing plane interface boundary markers. */
  nMarker_Turbomachinery,         /*!< \brief Number turbomachinery markers. */
  nMarker_TurboPerformance,       /*!< \brief Number of turboperformance markers. */
  nSpanWiseSections_User,         /*!< \brief Number of spanwise sections to compute 3D BC and Performance for turbomachinery   */
  nMarker_Shroud,                 /*!< \brief Number of shroud markers to set grid velocity to 0.*/
  nMarker_NearFieldBound,         /*!< \brief Number of near field boundary markers. */
  nMarker_ActDiskInlet,           /*!< \brief Number of actuator disk inlet markers. */
  nMarker_ActDiskOutlet,          /*!< \brief Number of actuator disk outlet markers. */
  nMarker_Deform_Mesh_Sym_Plane,  /*!< \brief Number of markers with symmetric deformation */
  nMarker_Deform_Mesh,            /*!< \brief Number of deformable markers at the boundary. */
  nMarker_Fluid_Load,             /*!< \brief Number of markers in which the flow load is computed/employed. */
  nMarker_Fluid_InterfaceBound,   /*!< \brief Number of fluid interface markers. */
  nMarker_CHTInterface,           /*!< \brief Number of conjugate heat transfer interface markers. */
  nMarker_Inlet,                  /*!< \brief Number of inlet flow markers. */
  nMarker_Inlet_Species,          /*!< \brief Number of inlet species markers. */
  nSpecies_per_Inlet,             /*!< \brief Number of species defined per inlet markers. */
  nMarker_Inlet_Turb,             /*!< \brief Number of inlet turbulent markers. */
  nTurb_Properties,               /*!< \brief Number of turbulent properties per inlet markers. */
  nMarker_Riemann,                /*!< \brief Number of Riemann flow markers. */
  nMarker_Giles,                  /*!< \brief Number of Giles flow markers. */
  nRelaxFactor_Giles,             /*!< \brief Number of relaxation factors for Giles markers. */
  nMarker_Supersonic_Inlet,       /*!< \brief Number of supersonic inlet flow markers. */
  nMarker_Supersonic_Outlet,      /*!< \brief Number of supersonic outlet flow markers. */
  nMarker_Outlet,                 /*!< \brief Number of outlet flow markers. */
  nMarker_Smoluchowski_Maxwell,   /*!< \brief Number of smoluchowski/maxwell wall boundaries. */
  nMarker_Isothermal,             /*!< \brief Number of isothermal wall boundaries. */
  nMarker_HeatFlux,               /*!< \brief Number of constant heat flux wall boundaries. */
  nMarker_HeatTransfer,           /*!< \brief Number of heat-transfer/convection wall boundaries. */
  nMarker_EngineExhaust,          /*!< \brief Number of nacelle exhaust flow markers. */
  nMarker_EngineInflow,           /*!< \brief Number of nacelle inflow flow markers. */
  nMarker_Clamped,                /*!< \brief Number of clamped markers in the FEM. */
  nMarker_Displacement,           /*!< \brief Number of displacement surface markers. */
  nMarker_Load,                   /*!< \brief Number of load surface markers. */
  nMarker_Damper,                 /*!< \brief Number of damper surface markers. */
  nMarker_Load_Dir,               /*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_Disp_Dir,               /*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_Load_Sine,              /*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_FlowLoad,               /*!< \brief Number of load surface markers. */
  nMarker_Internal,               /*!< \brief Number of internal flow markers. */
  nMarker_All,                    /*!< \brief Total number of markers using the grid information. */
  nMarker_Max,                    /*!< \brief Max number of number of markers using the grid information. */
  nMarker_CfgFile;                /*!< \brief Total number of markers using the config file (note that in
                                        parallel computations this number can be different from nMarker_All). */

  bool Inlet_From_File;         /*!< \brief True if the inlet profile is to be loaded from a file. */
  string Inlet_Filename;        /*!< \brief Filename specifying an inlet profile. */
  su2double Inlet_Matching_Tol; /*!< \brief Tolerance used when matching a point to a point from the inlet file. */
  string ActDisk_FileName;      /*!< \brief Filename specifying an actuator disk. */

  string *Marker_Euler,           /*!< \brief Euler wall markers. */
  *Marker_FarField,               /*!< \brief Far field markers. */
  *Marker_Custom,
  *Marker_SymWall,                /*!< \brief Symmetry wall markers. */
  *Marker_PerBound,               /*!< \brief Periodic boundary markers. */
  *Marker_PerDonor,               /*!< \brief Rotationally periodic boundary donor markers. */
  *Marker_MixingPlaneInterface,   /*!< \brief MixingPlane interface boundary markers. */
  *Marker_TurboBoundIn,           /*!< \brief Turbomachinery performance boundary markers. */
  *Marker_TurboBoundOut,          /*!< \brief Turbomachinery performance boundary donor markers. */
  *Marker_NearFieldBound,         /*!< \brief Near Field boundaries markers. */
  *Marker_Deform_Mesh,            /*!< \brief Deformable markers at the boundary. */
  *Marker_Deform_Mesh_Sym_Plane,  /*!< \brief Marker with symmetric deformation. */
  *Marker_Fluid_Load,             /*!< \brief Markers in which the flow load is computed/employed. */
  *Marker_Fluid_InterfaceBound,   /*!< \brief Fluid interface markers. */
  *Marker_CHTInterface,           /*!< \brief Conjugate heat transfer interface markers. */
  *Marker_ActDiskInlet,           /*!< \brief Actuator disk inlet markers. */
  *Marker_ActDiskOutlet,          /*!< \brief Actuator disk outlet markers. */
  *Marker_Inlet,                  /*!< \brief Inlet flow markers. */
  *Marker_Inlet_Species,          /*!< \brief Inlet species markers. */
  *Marker_Inlet_Turb,             /*!< \brief Inlet turbulent markers. */
  *Marker_Riemann,                /*!< \brief Riemann markers. */
  *Marker_Giles,                  /*!< \brief Giles markers. */
  *Marker_Shroud,                 /*!< \brief Shroud markers. */
  *Marker_Supersonic_Inlet,       /*!< \brief Supersonic inlet flow markers. */
  *Marker_Supersonic_Outlet,      /*!< \brief Supersonic outlet flow markers. */
  *Marker_Outlet,                 /*!< \brief Outlet flow markers. */
  *Marker_Smoluchowski_Maxwell,   /*!< \brief Smoluchowski/Maxwell wall markers. */
  *Marker_Isothermal,             /*!< \brief Isothermal wall markers. */
  *Marker_HeatFlux,               /*!< \brief Constant heat flux wall markers. */
  *Marker_HeatTransfer,           /*!< \brief Heat-transfer/convection markers. */
  *Marker_RoughWall,              /*!< \brief Constant heat flux wall markers. */
  *Marker_EngineInflow,           /*!< \brief Engine Inflow flow markers. */
  *Marker_EngineExhaust,          /*!< \brief Engine Exhaust flow markers. */
  *Marker_Clamped,                /*!< \brief Clamped markers. */
  *Marker_Displacement,           /*!< \brief Displacement markers. */
  *Marker_Load,                   /*!< \brief Load markers. */
  *Marker_Damper,                 /*!< \brief Damper markers. */
  *Marker_Load_Dir,               /*!< \brief Load markers defined in cartesian coordinates. */
  *Marker_Disp_Dir,               /*!< \brief Load markers defined in cartesian coordinates. */
  *Marker_Load_Sine,              /*!< \brief Sine-wave loaded markers defined in cartesian coordinates. */
  *Marker_FlowLoad,               /*!< \brief Flow Load markers. */
  *Marker_Internal,               /*!< \brief Internal flow markers. */
  *Marker_All_TagBound;           /*!< \brief Global index for markers using grid information. */

  su2double *Exhaust_Temperature_Target;     /*!< \brief Specified total temperatures for nacelle boundaries. */
  su2double *Exhaust_Pressure_Target;        /*!< \brief Specified total pressures for nacelle boundaries. */
  su2double *Inlet_Ttotal;                   /*!< \brief Specified total temperatures for inlet boundaries. */
  su2double *Riemann_Var1, *Riemann_Var2;    /*!< \brief Specified values for Riemann boundary. */
  su2double **Riemann_FlowDir;               /*!< \brief Specified flow direction vector (unit vector) for Riemann boundaries. */
  su2double *Giles_Var1, *Giles_Var2,
  *RelaxFactorAverage, *RelaxFactorFourier;  /*!< \brief Specified values for Giles BC. */
  su2double **Giles_FlowDir;                 /*!< \brief Specified flow direction vector (unit vector) for Giles BC. */
  su2double *Inlet_Ptotal;                   /*!< \brief Specified total pressures for inlet boundaries. */
  su2double **Inlet_FlowDir;                 /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
  su2double *Inlet_Temperature;              /*!< \brief Specified temperatures for a supersonic inlet boundaries. */
  su2double *Inlet_Pressure;                 /*!< \brief Specified static pressures for supersonic inlet boundaries. */
  su2double **Inlet_Velocity;                /*!< \brief Specified flow velocity vectors for supersonic inlet boundaries. */
  su2double **Inlet_SpeciesVal;              /*!< \brief Specified species vector for inlet boundaries. */
  su2double **Inlet_TurbVal;                 /*!< \brief Specified turbulent intensity and viscosity ratio for inlet boundaries. */
  su2double *EngineInflow_Target;            /*!< \brief Specified fan face targets for nacelle boundaries. */
  su2double *Inflow_Mach;                    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_Pressure;                /*!< \brief Specified fan face pressure for nacelle boundaries. */
  su2double *Inflow_MassFlow;                /*!< \brief Specified fan face massflow for nacelle boundaries. */
  su2double *Inflow_ReverseMassFlow;         /*!< \brief Specified fan face reverse massflow for nacelle boundaries. */
  su2double *Inflow_TotalPressure;           /*!< \brief Specified fan face total pressure for nacelle boundaries. */
  su2double *Inflow_Temperature;             /*!< \brief Specified fan face temperature for nacelle boundaries. */
  su2double *Inflow_TotalTemperature;        /*!< \brief Specified fan face total temperature for nacelle boundaries. */
  su2double *Inflow_RamDrag;                 /*!< \brief Specified fan face ram drag for nacelle boundaries. */
  su2double *Inflow_Force;                   /*!< \brief Specified force for nacelle boundaries. */
  su2double *Inflow_Power;                   /*!< \brief Specified power for nacelle boundaries. */
  su2double *Exhaust_Pressure;               /*!< \brief Specified exhaust pressure for nacelle boundaries. */
  su2double *Exhaust_Temperature;            /*!< \brief Specified exhaust temperature for nacelle boundaries. */
  su2double *Exhaust_MassFlow;               /*!< \brief Specified exhaust mass flow for nacelle boundaries. */
  su2double *Exhaust_TotalPressure;          /*!< \brief Specified exhaust total pressure for nacelle boundaries. */
  su2double *Exhaust_TotalTemperature;       /*!< \brief Specified exhaust total temperature for nacelle boundaries. */
  su2double *Exhaust_GrossThrust;            /*!< \brief Specified exhaust gross thrust for nacelle boundaries. */
  su2double *Exhaust_Force;                  /*!< \brief Specified exhaust force for nacelle boundaries. */
  su2double *Exhaust_Power;                  /*!< \brief Specified exhaust power for nacelle boundaries. */
  su2double *Engine_Power;                   /*!< \brief Specified engine power for nacelle boundaries. */
  su2double *Engine_Mach;                    /*!< \brief Specified engine mach for nacelle boundaries. */
  su2double *Engine_Force;                   /*!< \brief Specified engine force for nacelle boundaries. */
  su2double *Engine_NetThrust;               /*!< \brief Specified engine net thrust for nacelle boundaries. */
  su2double *Engine_GrossThrust;             /*!< \brief Specified engine gross thrust for nacelle boundaries. */
  su2double *Engine_Area;                    /*!< \brief Specified engine area for nacelle boundaries. */
  su2double *Outlet_Pressure;                /*!< \brief Specified back pressures (static) for outlet boundaries. */
  su2double *Isothermal_Temperature;         /*!< \brief Specified isothermal wall temperatures (static). */
  su2double *HeatTransfer_Coeff;             /*!< \brief Specified heat transfer coefficients. */
  su2double *HeatTransfer_WallTemp;          /*!< \brief Specified temperatures at infinity alongside heat transfer coefficients. */
  su2double *Heat_Flux;                      /*!< \brief Specified wall heat fluxes. */
  su2double *Roughness_Height;               /*!< \brief Equivalent sand grain roughness for the marker according to config file. */
  su2double *Displ_Value;                    /*!< \brief Specified displacement for displacement boundaries. */
  su2double *Load_Value;                     /*!< \brief Specified force for load boundaries. */
  su2double *Damper_Constant;                /*!< \brief Specified constant for damper boundaries. */
  su2double *Load_Dir_Value;                 /*!< \brief Specified force for load boundaries defined in cartesian coordinates. */
  su2double *Load_Dir_Multiplier;            /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double *Disp_Dir_Value;                 /*!< \brief Specified force for load boundaries defined in cartesian coordinates. */
  su2double *Disp_Dir_Multiplier;            /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double **Load_Dir;                      /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
  su2double **Disp_Dir;                      /*!< \brief Specified structural displacement direction (unit vector). */
  su2double *Load_Sine_Amplitude;            /*!< \brief Specified amplitude for a sine-wave load. */
  su2double *Load_Sine_Frequency;            /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double **Load_Sine_Dir;                 /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
  su2double *FlowLoad_Value;                 /*!< \brief Specified force for flow load boundaries. */
  su2double *ActDiskInlet_MassFlow;          /*!< \brief Specified inlet mass flow for actuator disk. */
  su2double *ActDiskInlet_Temperature;       /*!< \brief Specified inlet temperature for actuator disk. */
  su2double *ActDiskInlet_TotalTemperature;  /*!< \brief Specified inlet total temperature for actuator disk. */
  su2double *ActDiskInlet_Pressure;          /*!< \brief Specified inlet pressure for actuator disk. */
  su2double *ActDiskInlet_TotalPressure;     /*!< \brief Specified inlet total pressure for actuator disk. */
  su2double *ActDiskInlet_RamDrag;           /*!< \brief Specified inlet ram drag for actuator disk. */
  su2double *ActDiskInlet_Force;             /*!< \brief Specified inlet force for actuator disk. */
  su2double *ActDiskInlet_Power;             /*!< \brief Specified inlet power for actuator disk. */
  su2double *ActDiskOutlet_MassFlow;         /*!< \brief Specified outlet mass flow for actuator disk. */
  su2double *ActDiskOutlet_Temperature;      /*!< \brief Specified outlet temperature for actuator disk. */
  su2double *ActDiskOutlet_TotalTemperature; /*!< \brief Specified outlet total temperatur for actuator disk. */
  su2double *ActDiskOutlet_Pressure;         /*!< \brief Specified outlet pressure for actuator disk. */
  su2double *ActDiskOutlet_TotalPressure;    /*!< \brief Specified outlet total pressure for actuator disk. */
  su2double *ActDiskOutlet_GrossThrust;      /*!< \brief Specified outlet gross thrust for actuator disk. */
  su2double *ActDiskOutlet_Force;            /*!< \brief Specified outlet force for actuator disk. */
  su2double *ActDiskOutlet_Power;            /*!< \brief Specified outlet power for actuator disk. */
  su2double **ActDisk_PressJump,
  **ActDisk_TempJump,  **ActDisk_Omega;      /*!< \brief Specified deltas for actuator disk.*/
  su2double *ActDisk_DeltaPress;             /*!< \brief Specified pressure delta for actuator disk. */
  su2double *ActDisk_DeltaTemp;              /*!< \brief Specified temperature delta for actuator disk. */
  su2double *ActDisk_TotalPressRatio;        /*!< \brief Specified tot. pres. ratio for actuator disk. */
  su2double *ActDisk_TotalTempRatio;         /*!< \brief Specified tot. temp. ratio for actuator disk. */
  su2double *ActDisk_StaticPressRatio;       /*!< \brief Specified press. ratio for actuator disk. */
  su2double *ActDisk_StaticTempRatio;        /*!< \brief Specified temp. ratio for actuator disk. */
  su2double *ActDisk_Power;                  /*!< \brief Specified power for actuator disk. */
  su2double *ActDisk_MassFlow;               /*!< \brief Specified mass flow for actuator disk. */
  su2double *ActDisk_Mach;                   /*!< \brief Specified mach for actuator disk. */
  su2double *ActDisk_Force;                  /*!< \brief Specified force for actuator disk. */
  su2double *Outlet_MassFlow;                /*!< \brief Mass flow for outlet boundaries. */
  su2double *Outlet_Density;                 /*!< \brief Avg. density for outlet boundaries. */
  su2double *Outlet_Area;                    /*!< \brief Area for outlet boundaries. */
  su2double *Surface_MassFlow;               /*!< \brief Massflow at the boundaries. */
  su2double *Surface_Mach;                   /*!< \brief Mach number at the boundaries. */
  su2double *Surface_Temperature;            /*!< \brief Temperature at the boundaries. */
  su2double *Surface_Pressure;               /*!< \brief Pressure at the boundaries. */
  su2double *Surface_Density;                /*!< \brief Density at the boundaries. */
  su2double *Surface_Enthalpy;               /*!< \brief Enthalpy at the boundaries. */
  su2double *Surface_NormalVelocity;         /*!< \brief Normal velocity at the boundaries. */
  su2double *Surface_Uniformity;             /*!< \brief Integral measure of the streamwise uniformity (absolute) at the boundaries (non-dim). */
  su2double *Surface_SecondaryStrength;      /*!< \brief Integral measure of the strength of secondary flows (absolute) at the boundaries (non-dim). */
  su2double *Surface_SecondOverUniform;      /*!< \brief Integral measure of the strength of secondary flows (relative to streamwise) at the boundaries (non-dim). */
  su2double *Surface_MomentumDistortion;     /*!< \brief Integral measure of the streamwise uniformity (relative to plug flow) at the boundaries (non-dim). */
  su2double *Surface_TotalTemperature;       /*!< \brief Total temperature at the boundaries. */
  su2double *Surface_TotalPressure;          /*!< \brief Total pressure at the boundaries. */
  su2double *Surface_PressureDrop;           /*!< \brief Pressure drop between boundaries. */
  su2double* Surface_Species_0;              /*!< \brief Average Species_0 at the boundaries. */
  su2double* Surface_Species_Variance;       /*!< \brief Species Variance at the boundaries. */
  su2double *Surface_DC60;                   /*!< \brief Specified surface DC60 for nacelle boundaries. */
  su2double *Surface_IDC;                    /*!< \brief Specified IDC for nacelle boundaries. */
  su2double *Surface_IDC_Mach;               /*!< \brief Specified IDC mach for nacelle boundaries. */
  su2double *Surface_IDR;                    /*!< \brief Specified surface IDR for nacelle boundaries. */
  su2double *ActDisk_NetThrust;              /*!< \brief Specified net thrust for nacelle boundaries. */
  su2double *ActDisk_BCThrust;               /*!< \brief Specified bc thrust for nacelle boundaries. */
  su2double *ActDisk_BCThrust_Old;           /*!< \brief Specified old bc thrust for nacelle boundaries. */
  su2double *ActDisk_GrossThrust;            /*!< \brief Specified gross thrust for nacelle boundaries. */
  su2double *ActDisk_Area;                   /*!< \brief Specified area for nacelle boundaries. */
  su2double *ActDisk_ReverseMassFlow;        /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double **Periodic_RotCenter;            /*!< \brief Rotational center for each periodic boundary. */
  su2double **Periodic_RotAngles;            /*!< \brief Rotation angles for each periodic boundary. */
  su2double **Periodic_Translation;          /*!< \brief Translation vector for each periodic boundary. */
  string *Marker_CfgFile_TagBound;           /*!< \brief Global index for markers using config file. */
  unsigned short *Marker_All_KindBC,         /*!< \brief Global index for boundaries using grid information. */
  *Marker_CfgFile_KindBC;                    /*!< \brief Global index for boundaries using config file. */
  short *Marker_All_SendRecv;                /*!< \brief Information about if the boundary is sended (+), received (-). */
  short *Marker_All_PerBound;                /*!< \brief Global index for periodic bc using the grid information. */

  unsigned long ExtIter;            /*!< \brief Current external iteration number. */
  unsigned long ExtIter_OffSet;     /*!< \brief External iteration number offset. */
  unsigned long IntIter;            /*!< \brief Current internal iteration number. */
  unsigned long OuterIter;          /*!< \brief Current Outer iterations for multizone problems. */
  unsigned long InnerIter;          /*!< \brief Current inner iterations for multizone problems. */
  unsigned long TimeIter;           /*!< \brief Current time iterations for multizone problems. */
  long Unst_AdjointIter;            /*!< \brief Iteration number to begin the reverse time integration in the direct solver for the unsteady adjoint. */
  long Iter_Avg_Objective;          /*!< \brief Iteration the number of time steps to be averaged, counting from the back */
  su2double PhysicalTime;           /*!< \brief Physical time at the current iteration in the solver for unsteady problems. */

  unsigned short nLevels_TimeAccurateLTS;   /*!< \brief Number of time levels for time accurate local time stepping. */
  unsigned short nTimeDOFsADER_DG;          /*!< \brief Number of time DOFs used in the predictor step of ADER-DG. */
  su2double *TimeDOFsADER_DG;               /*!< \brief The location of the ADER-DG time DOFs on the interval [-1,1]. */
  unsigned short nTimeIntegrationADER_DG;   /*!< \brief Number of time integration points ADER-DG. */
  su2double *TimeIntegrationADER_DG;        /*!< \brief The location of the ADER-DG time integration points on the interval [-1,1]. */
  su2double *WeightsIntegrationADER_DG;     /*!< \brief The weights of the ADER-DG time integration points on the interval [-1,1]. */
  unsigned short nRKStep;                   /*!< \brief Number of steps of the explicit Runge-Kutta method. */
  su2double *RK_Alpha_Step;                 /*!< \brief Runge-Kutta beta coefficients. */

  unsigned short nQuasiNewtonSamples;  /*!< \brief Number of samples used in quasi-Newton solution methods. */
  bool UseVectorization;       /*!< \brief Whether to use vectorized numerics schemes. */
  bool NewtonKrylov;           /*!< \brief Use a coupled Newton method to solve the flow equations. */
  array<unsigned short,3> NK_IntParam{{20, 3, 2}}; /*!< \brief Integer parameters for NK method. */
  array<su2double,4> NK_DblParam{{-2.0, 0.1, -3.0, 1e-4}}; /*!< \brief Floating-point parameters for NK method. */

  unsigned short nMGLevels;    /*!< \brief Number of multigrid levels (coarse levels). */
  unsigned short nCFL;         /*!< \brief Number of CFL, one for each multigrid level. */
  su2double
  CFLRedCoeff_Turb,            /*!< \brief CFL reduction coefficient on the LevelSet problem. */
  CFLRedCoeff_AdjFlow,         /*!< \brief CFL reduction coefficient for the adjoint problem. */
  CFLRedCoeff_AdjTurb,         /*!< \brief CFL reduction coefficient for the adjoint turbulent problem. */
  CFLRedCoeff_Species,         /*!< \brief CFL reduction coefficient on the species problem. */
  CFLFineGrid,                 /*!< \brief CFL of the finest grid. */
  Max_DeltaTime,               /*!< \brief Max delta time. */
  Unst_CFL;                    /*!< \brief Unsteady CFL number. */

  /* Gradient smoothing options */
  su2double SmoothingEps1;          /*!< \brief Parameter for the identity part in gradient smoothing. */
  su2double SmoothingEps2;          /*!< \brief Parameter for the Laplace part in gradient smoothing. */
  bool SmoothGradient;              /*!< \brief Flag for enabling gradient smoothing. */
  bool SmoothSepDim;                /*!< \brief Flag for enabling separated calculation for every dimension. */
  bool SmoothOnSurface;             /*!< \brief Flag for assembling the system only on the surface. */
  bool SmoothDirichletSurfaceBound; /*!< \brief Flag for using zero Dirichlet boundary in the surface case. */
  ENUM_SOBOLEV_MODUS SmoothNumMode; /*!< \brief The mode in which the Sobolev smoothing solver is applied. */

  unsigned short  Kind_Grad_Linear_Solver,  /*!< Numerical method to smooth the gradient */
  Kind_Grad_Linear_Solver_Prec;             /*!< \brief Preconditioner of the linear solver. */
  su2double Grad_Linear_Solver_Error;       /*!< \brief Min error of the linear solver for the gradient smoothing. */
  unsigned long Grad_Linear_Solver_Iter; /*!< \brief Max iterations of the linear solver for the gradient smoothing. */

  bool ReorientElements;       /*!< \brief Flag for enabling element reorientation. */
  string CustomObjFunc;        /*!< \brief User-defined objective function. */
  string CustomOutputs;        /*!< \brief User-defined functions for outputs. */
  unsigned short nDV,                  /*!< \brief Number of design variables. */
  nObj, nObjW;                         /*! \brief Number of objective functions. */
  unsigned short* nDV_Value;           /*!< \brief Number of values for each design variable (might be different than 1 if we allow arbitrary movement). */
  unsigned short nFFDBox;              /*!< \brief Number of ffd boxes. */
  unsigned short nTurboMachineryKind;  /*!< \brief Number turbomachinery types specified. */
  unsigned short nParamDV;             /*!< \brief Number of parameters of the design variable. */
  string DV_Filename;                  /*!< \brief Filename for providing surface positions from an external parameterization. */
  string DV_Unordered_Sens_Filename;   /*!< \brief Filename of volume sensitivities in an unordered ASCII format. */
  string DV_Sens_Filename;             /*!< \brief Filename of surface sensitivities written to an unordered ASCII format. */
  unsigned short
  Sensitivity_FileFormat;             /*!< \brief Format of the input volume sensitivity files (SU2_DOT). */
  su2double **ParamDV;                /*!< \brief Parameters of the design variable. */
  su2double **CoordFFDBox;            /*!< \brief Coordinates of the FFD boxes. */
  unsigned short **DegreeFFDBox;      /*!< \brief Degree of the FFD boxes. */
  string *FFDTag;                     /*!< \brief Parameters of the design variable. */
  string *TagFFDBox;                  /*!< \brief Tag of the FFD box. */
  unsigned short GeometryMode;        /*!< \brief Gemoetry mode (analysis or gradient computation). */
  unsigned short MGCycle;             /*!< \brief Kind of multigrid cycle. */
  unsigned short FinestMesh;          /*!< \brief Finest mesh for the full multigrid approach. */
  unsigned short nFFD_Fix_IDir,
  nFFD_Fix_JDir, nFFD_Fix_KDir;       /*!< \brief Number of planes fixed in the FFD. */
  unsigned short nMG_PreSmooth,       /*!< \brief Number of MG pre-smooth parameters found in config file. */
  nMG_PostSmooth,                     /*!< \brief Number of MG post-smooth parameters found in config file. */
  nMG_CorrecSmooth;                   /*!< \brief Number of MG correct-smooth parameters found in config file. */
  short *FFD_Fix_IDir,
  *FFD_Fix_JDir, *FFD_Fix_KDir;       /*!< \brief Exact sections. */
  unsigned short *MG_PreSmooth,       /*!< \brief Multigrid Pre smoothing. */
  *MG_PostSmooth,                     /*!< \brief Multigrid Post smoothing. */
  *MG_CorrecSmooth;                   /*!< \brief Multigrid Jacobi implicit smoothing of the correction. */
  su2double *LocationStations;        /*!< \brief Airfoil sections in wing slicing subroutine. */

  ENUM_MULTIZONE Kind_MZSolver;    /*!< \brief Kind of multizone solver.  */
  INC_DENSITYMODEL Kind_DensityModel; /*!< \brief Kind of the density model for incompressible flows. */
  CHT_COUPLING Kind_CHT_Coupling;  /*!< \brief Kind of coupling method used at CHT interfaces. */
  VISCOSITYMODEL Kind_ViscosityModel; /*!< \brief Kind of the Viscosity Model*/
  MIXINGVISCOSITYMODEL Kind_MixingViscosityModel; /*!< \brief Kind of the mixing Viscosity Model*/
  CONDUCTIVITYMODEL Kind_ConductivityModel; /*!< \brief Kind of the Thermal Conductivity Model */
  CONDUCTIVITYMODEL_TURB Kind_ConductivityModel_Turb; /*!< \brief Kind of the Turbulent Thermal Conductivity Model */
  DIFFUSIVITYMODEL Kind_Diffusivity_Model; /*!< \brief Kind of the mass diffusivity Model */
  FREESTREAM_OPTION Kind_FreeStreamOption; /*!< \brief Kind of free stream option to choose if initializing with density or temperature  */
  MAIN_SOLVER Kind_Solver;         /*!< \brief Kind of solver: Euler, NS, Continuous adjoint, etc.  */
  LIMITER Kind_SlopeLimit,    /*!< \brief Global slope limiter. */
  Kind_SlopeLimit_Flow,         /*!< \brief Slope limiter for flow equations.*/
  Kind_SlopeLimit_Turb,         /*!< \brief Slope limiter for the turbulence equation.*/
  Kind_SlopeLimit_AdjTurb,      /*!< \brief Slope limiter for the adjoint turbulent equation.*/
  Kind_SlopeLimit_AdjFlow,      /*!< \brief Slope limiter for the adjoint equation.*/
  Kind_SlopeLimit_Heat,         /*!< \brief Slope limiter for the adjoint equation.*/
  Kind_SlopeLimit_Species;      /*!< \brief Slope limiter for the species equation.*/
  unsigned short Kind_FluidModel,  /*!< \brief Kind of the Fluid Model: Ideal, van der Waals, etc. */
  Kind_InitOption,                 /*!< \brief Kind of Init option to choose if initializing with Reynolds number or with thermodynamic conditions   */
  Kind_GridMovement,               /*!< \brief Kind of the static mesh movement. */
  *Kind_SurfaceMovement,           /*!< \brief Kind of the static mesh movement. */
  nKind_SurfaceMovement,           /*!< \brief Kind of the dynamic mesh movement. */
  Kind_Gradient_Method,            /*!< \brief Numerical method for computation of spatial gradients. */
  Kind_Gradient_Method_Recon,      /*!< \brief Numerical method for computation of spatial gradients used for upwind reconstruction. */
  Kind_Deform_Linear_Solver,             /*!< Numerical method to deform the grid */
  Kind_Deform_Linear_Solver_Prec,        /*!< \brief Preconditioner of the linear solver. */
  Kind_Linear_Solver,                    /*!< \brief Numerical solver for the implicit scheme. */
  Kind_Linear_Solver_Prec,               /*!< \brief Preconditioner of the linear solver. */
  Kind_AdjTurb_Linear_Solver,            /*!< \brief Numerical solver for the turbulent adjoint implicit scheme. */
  Kind_AdjTurb_Linear_Prec,              /*!< \brief Preconditioner of the turbulent adjoint linear solver. */
  Kind_DiscAdj_Linear_Solver,            /*!< \brief Linear solver for the discrete adjoint system. */
  Kind_DiscAdj_Linear_Prec,              /*!< \brief Preconditioner of the discrete adjoint linear solver. */
  Kind_TimeNumScheme,           /*!< \brief Global explicit or implicit time integration. */
  Kind_TimeIntScheme_Flow,      /*!< \brief Time integration for the flow equations. */
  Kind_TimeIntScheme_FEM_Flow,  /*!< \brief Time integration for the flow equations. */
  Kind_ADER_Predictor,          /*!< \brief Predictor step of the ADER-DG time integration scheme. */
  Kind_TimeIntScheme_AdjFlow,   /*!< \brief Time integration for the adjoint flow equations. */
  Kind_TimeIntScheme_Turb,      /*!< \brief Time integration for the turbulence model. */
  Kind_TimeIntScheme_AdjTurb,   /*!< \brief Time integration for the adjoint turbulence model. */
  Kind_TimeIntScheme_Species,   /*!< \brief Time integration for the species model. */
  Kind_TimeIntScheme_Heat,      /*!< \brief Time integration for the wave equations. */
  Kind_TimeStep_Heat,           /*!< \brief Time stepping method for the (fvm) heat equation. */
  n_Datadriven_files;
  ENUM_DATADRIVEN_METHOD Kind_DataDriven_Method;       /*!< \brief Method used for datset regression in data-driven fluid models. */

  su2double DataDriven_Relaxation_Factor; /*!< \brief Relaxation factor for Newton solvers in data-driven fluid models. */

  STRUCT_TIME_INT Kind_TimeIntScheme_FEA;    /*!< \brief Time integration for the FEA equations. */
  STRUCT_SPACE_ITE Kind_SpaceIteScheme_FEA;  /*!< \brief Iterative scheme for nonlinear structural analysis. */
  unsigned short
  Kind_TimeIntScheme_Radiation, /*!< \brief Time integration for the Radiation equations. */
  Kind_ConvNumScheme,           /*!< \brief Global definition of the convective term. */
  Kind_ConvNumScheme_Flow,      /*!< \brief Centered or upwind scheme for the flow equations. */
  Kind_ConvNumScheme_FEM_Flow,  /*!< \brief Finite element scheme for the flow equations. */
  Kind_ConvNumScheme_Heat,      /*!< \brief Centered or upwind scheme for the flow equations. */
  Kind_ConvNumScheme_AdjFlow,   /*!< \brief Centered or upwind scheme for the adjoint flow equations. */
  Kind_ConvNumScheme_Turb,      /*!< \brief Centered or upwind scheme for the turbulence model. */
  Kind_ConvNumScheme_AdjTurb,   /*!< \brief Centered or upwind scheme for the adjoint turbulence model. */
  Kind_ConvNumScheme_Species,   /*!< \brief Centered or upwind scheme for the species model. */
  Kind_ConvNumScheme_Template,  /*!< \brief Centered or upwind scheme for the level set equation. */
  Kind_FEM,                     /*!< \brief Finite element scheme for the flow equations. */
  Kind_FEM_Flow,                /*!< \brief Finite element scheme for the flow equations. */
  Kind_Matrix_Coloring;         /*!< \brief Type of matrix coloring for sparse Jacobian computation. */

  CENTERED
  Kind_Centered,                /*!< \brief Centered scheme. */
  Kind_Centered_Flow,           /*!< \brief Centered scheme for the flow equations. */
  Kind_Centered_AdjFlow,        /*!< \brief Centered scheme for the adjoint flow equations. */
  Kind_Centered_Turb,           /*!< \brief Centered scheme for the turbulence model. */
  Kind_Centered_AdjTurb,        /*!< \brief Centered scheme for the adjoint turbulence model. */
  Kind_Centered_Species,        /*!< \brief Centered scheme for the species model. */
  Kind_Centered_Heat,           /*!< \brief Centered scheme for the heat transfer model. */
  Kind_Centered_Template;       /*!< \brief Centered scheme for the template model. */


  FEM_SHOCK_CAPTURING_DG Kind_FEM_Shock_Capturing_DG; /*!< \brief Shock capturing method for the FEM DG solver. */
  BGS_RELAXATION Kind_BGS_RelaxMethod; /*!< \brief Kind of relaxation method for Block Gauss Seidel method in FSI problems. */
  bool ReconstructionGradientRequired; /*!< \brief Enable or disable a second gradient calculation for upwind reconstruction only. */
  bool LeastSquaresRequired;    /*!< \brief Enable or disable memory allocation for least-squares gradient methods. */
  bool Energy_Equation;         /*!< \brief Solve the energy equation for incompressible flows. */

  UPWIND
  Kind_Upwind,                  /*!< \brief Upwind scheme. */
  Kind_Upwind_Flow,             /*!< \brief Upwind scheme for the flow equations. */
  Kind_Upwind_AdjFlow,          /*!< \brief Upwind scheme for the adjoint flow equations. */
  Kind_Upwind_Turb,             /*!< \brief Upwind scheme for the turbulence model. */
  Kind_Upwind_AdjTurb,          /*!< \brief Upwind scheme for the adjoint turbulence model. */
  Kind_Upwind_Species,          /*!< \brief Upwind scheme for the species model. */
  Kind_Upwind_Heat,             /*!< \brief Upwind scheme for the heat transfer model. */
  Kind_Upwind_Template;         /*!< \brief Upwind scheme for the template model. */

  bool MUSCL,              /*!< \brief MUSCL scheme .*/
  MUSCL_Flow,              /*!< \brief MUSCL scheme for the flow equations.*/
  MUSCL_Turb,              /*!< \brief MUSCL scheme for the turbulence equations.*/
  MUSCL_Heat,              /*!< \brief MUSCL scheme for the (fvm) heat equation.*/
  MUSCL_AdjFlow,           /*!< \brief MUSCL scheme for the adj flow equations.*/
  MUSCL_AdjTurb;           /*!< \brief MUSCL scheme for the adj turbulence equations.*/
  bool MUSCL_Species;      /*!< \brief MUSCL scheme for the species equations.*/
  bool Use_Accurate_Jacobians;  /*!< \brief Use numerically computed Jacobians for AUSM+up(2) and SLAU(2). */
  bool EulerPersson;       /*!< \brief Boolean to determine whether this is an Euler simulation with Persson shock capturing. */
  bool FSI_Problem = false,/*!< \brief Boolean to determine whether the simulation is FSI or not. */
  Multizone_Problem;       /*!< \brief Boolean to determine whether we are solving a multizone problem. */
  unsigned short nID_DV;   /*!< \brief ID for the region of FEM when computed using direct differentiation. */

  bool AD_Mode;             /*!< \brief Algorithmic Differentiation support. */
  bool AD_Preaccumulation;  /*!< \brief Enable or disable preaccumulation in the AD mode. */
  STRUCT_COMPRESS Kind_Material_Compress;  /*!< \brief Determines if the material is compressible or incompressible (structural analysis). */
  STRUCT_MODEL Kind_Material;              /*!< \brief Determines the material model to be used (structural analysis). */
  STRUCT_DEFORMATION Kind_Struct_Solver;   /*!< \brief Determines the geometric condition (small or large deformations) for structural analysis. */
  unsigned short Kind_DV_FEA;              /*!< \brief Kind of Design Variable for FEA problems.*/

  unsigned short nTurbVar;          /*!< \brief Number of Turbulence variables, i.e. 1 for SA-types, 2 for SST. */
  TURB_MODEL Kind_Turb_Model;       /*!< \brief Turbulent model definition. */
  SPECIES_MODEL Kind_Species_Model; /*!< \brief Species model definition. */
  TURB_SGS_MODEL Kind_SGS_Model;    /*!< \brief LES SGS model definition. */
  TURB_TRANS_MODEL Kind_Trans_Model;  /*!< \brief Transition model definition. */
  TURB_TRANS_CORRELATION Kind_Trans_Correlation;  /*!< \brief Transition correlation model definition. */
  su2double hRoughness;             /*!< \brief RMS roughness for Transition model. */
  unsigned short Kind_ActDisk, Kind_Engine_Inflow,
  *Kind_Data_Riemann,
  *Kind_Data_Giles;                /*!< \brief Kind of inlet boundary treatment. */
  INLET_TYPE Kind_Inlet;
  INLET_TYPE *Kind_Inc_Inlet;
  INC_OUTLET_TYPE *Kind_Inc_Outlet;
  unsigned short nWall_Types;      /*!< \brief Number of wall treatment types listed. */
  unsigned short nInc_Inlet;       /*!< \brief Number of inlet boundary treatment types listed. */
  unsigned short nInc_Outlet;      /*!< \brief Number of inlet boundary treatment types listed. */
  su2double Inc_Inlet_Damping;     /*!< \brief Damping factor applied to the iterative updates to the velocity at a pressure inlet in incompressible flow. */
  su2double Inc_Outlet_Damping;    /*!< \brief Damping factor applied to the iterative updates to the pressure at a mass flow outlet in incompressible flow. */
  bool Inc_Inlet_UseNormal;        /*!< \brief Flag for whether to use the local normal as the flow direction for an incompressible pressure inlet. */
  su2double Linear_Solver_Error;   /*!< \brief Min error of the linear solver for the implicit formulation. */
  su2double Deform_Linear_Solver_Error;          /*!< \brief Min error of the linear solver for the implicit formulation. */
  su2double Linear_Solver_Smoother_Relaxation;   /*!< \brief Relaxation factor for iterative linear smoothers. */
  unsigned long Linear_Solver_Iter;              /*!< \brief Max iterations of the linear solver for the implicit formulation. */
  unsigned long Deform_Linear_Solver_Iter;       /*!< \brief Max iterations of the linear solver for the implicit formulation. */
  unsigned long Linear_Solver_Restart_Frequency; /*!< \brief Restart frequency of the linear solver for the implicit formulation. */
  unsigned long Linear_Solver_Prec_Threads;      /*!< \brief Number of threads per rank for ILU and LU_SGS preconditioners. */
  unsigned short Linear_Solver_ILU_n;            /*!< \brief ILU fill=in level. */
  su2double SemiSpan;                   /*!< \brief Wing Semi span. */
  su2double Roe_Kappa;                  /*!< \brief Relaxation of the Roe scheme. */
  su2double Relaxation_Factor_Adjoint;  /*!< \brief Relaxation coefficient for variable updates of adjoint solvers. */
  su2double Relaxation_Factor_CHT;      /*!< \brief Relaxation coefficient for the update of conjugate heat variables. */
  su2double AdjTurb_Linear_Error;       /*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  su2double EntropyFix_Coeff;           /*!< \brief Entropy fix coefficient. */
  unsigned short AdjTurb_Linear_Iter;   /*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  unsigned short nLocationStations,     /*!< \brief Number of section cuts to make when outputting mesh and cp . */
  nWingStations;                        /*!< \brief Number of section cuts to make when calculating internal volume. */
  su2double Kappa_1st_AdjFlow,  /*!< \brief Lax 1st order dissipation coefficient for adjoint flow equations (coarse multigrid levels). */
  Kappa_2nd_AdjFlow,            /*!< \brief JST 2nd order dissipation coefficient for adjoint flow equations. */
  Kappa_4th_AdjFlow,            /*!< \brief JST 4th order dissipation coefficient for adjoint flow equations. */
  Kappa_1st_Flow,           /*!< \brief Lax 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
  Kappa_2nd_Flow,           /*!< \brief JST 2nd order dissipation coefficient for flow equations. */
  Kappa_4th_Flow,           /*!< \brief JST 4th order dissipation coefficient for flow equations. */
  Cent_Jac_Fix_Factor,              /*!< \brief Multiply the dissipation contribution to the Jacobian of central schemes
                                                by this factor to make the global matrix more diagonal dominant. */
  Cent_Inc_Jac_Fix_Factor;          /*!< \brief Multiply the dissipation contribution to the Jacobian of incompressible central schemes */
  su2double Geo_Waterline_Location; /*!< \brief Location of the waterline. */

  su2double Min_Beta_RoeTurkel,     /*!< \brief Minimum value of Beta for the Roe-Turkel low Mach preconditioner. */
  Max_Beta_RoeTurkel;               /*!< \brief Maximum value of Beta for the Roe-Turkel low Mach preconditioner. */
  unsigned long GridDef_Nonlinear_Iter;  /*!< \brief Number of nonlinear increments for grid deformation. */
  unsigned short Deform_StiffnessType;   /*!< \brief Type of element stiffness imposed for FEA mesh deformation. */
  bool Deform_Mesh;                      /*!< \brief Determines whether the mesh will be deformed. */
  bool Deform_Output;                    /*!< \brief Print the residuals during mesh deformation to the console. */
  su2double Deform_Tol_Factor;       /*!< \brief Factor to multiply smallest volume for deform tolerance (0.001 default) */
  su2double Deform_Coeff;            /*!< \brief Deform coeffienct */
  su2double Deform_Limit;            /*!< \brief Deform limit */
  unsigned short FFD_Continuity;     /*!< \brief Surface continuity at the intersection with the FFD */
  unsigned short FFD_CoordSystem;    /*!< \brief Define the coordinates system */
  su2double Deform_ElasticityMod,    /*!< \brief Young's modulus for volume deformation stiffness model */
  Deform_PoissonRatio,               /*!< \brief Poisson's ratio for volume deformation stiffness model */
  Deform_StiffLayerSize;             /*!< \brief Size of the layer of highest stiffness for wall distance-based mesh stiffness */
  bool FFD_Symmetry_Plane;           /*!< \brief FFD symmetry plane. */

  su2double Mach;             /*!< \brief Mach number. */
  su2double Reynolds;         /*!< \brief Reynolds number. */
  su2double Froude;           /*!< \brief Froude number. */
  su2double Length_Reynolds;  /*!< \brief Reynolds length (dimensional). */
  su2double AoA,              /*!< \brief Angle of attack (just external flow). */
  iH, AoS, AoA_Offset,
  AoS_Offset, AoA_Sens;       /*!< \brief Angle of sideSlip (just external flow). */
  bool Fixed_CL_Mode;         /*!< \brief Activate fixed CL mode (external flow only). */
  bool Eval_dOF_dCX;          /*!< \brief Activate fixed CL mode (external flow only). */
  bool Discard_InFiles;       /*!< \brief Discard angle of attack in solution and geometry files. */
  su2double Target_CL;        /*!< \brief Specify a target CL instead of AoA (external flow only). */
  su2double Total_CM;         /*!< \brief Specify a Total CM instead of AoA (external flow only). */
  su2double Total_CD;         /*!< \brief Specify a target CD instead of AoA (external flow only). */
  su2double dCL_dAlpha;       /*!< \brief value of dCl/dAlpha. */
  su2double dCM_diH;          /*!< \brief value of dCM/dHi. */
  unsigned long Iter_Fixed_CM;          /*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Iter_Fixed_NetThrust;   /*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Iter_dCL_dAlpha;        /*!< \brief Number of iterations to evaluate dCL_dAlpha. */
  unsigned long Update_Alpha;           /*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Update_iH;              /*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Update_BCThrust;        /*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  su2double dNetThrust_dBCThrust;       /*!< \brief value of dNetThrust/dBCThrust. */
  bool Update_BCThrust_Bool;            /*!< \brief Boolean flag for whether to update the AoA for fixed lift mode on a given iteration. */
  bool Update_AoA;                      /*!< \brief Boolean flag for whether to update the AoA for fixed lift mode on a given iteration. */
  unsigned long Update_AoA_Iter_Limit;  /*!< \brief Limit on number of iterations between AoA updates for fixed lift mode. */
  bool Finite_Difference_Mode;        /*!< \brief Flag to run the finite difference mode in fixed Cl mode. */
  su2double ChargeCoeff;              /*!< \brief Charge coefficient (just for poisson problems). */
  unsigned short Cauchy_Func_Flow,    /*!< \brief Function where to apply the convergence criteria in the flow problem. */
  Cauchy_Func_AdjFlow,                /*!< \brief Function where to apply the convergence criteria in the adjoint problem. */
  Cauchy_Elems;                       /*!< \brief Number of elements to evaluate. */
  unsigned short Residual_Func_Flow;  /*!< \brief Equation to apply residual convergence to. */
  unsigned short Res_FEM_CRIT;        /*!< \brief Criteria to apply to the FEM convergence (absolute/relative). */
  unsigned long StartConv_Iter;       /*!< \brief Start convergence criteria at iteration. */
  su2double Cauchy_Eps;               /*!< \brief Epsilon used for the convergence. */
  bool Restart,                       /*!< \brief Restart solution (for direct, adjoint, and linearized problems).*/
  Read_Binary_Restart,                /*!< \brief Read binary SU2 native restart files.*/
  Wrt_Restart_Overwrite,              /*!< \brief Overwrite restart files or append iteration number.*/
  Wrt_Surface_Overwrite,              /*!< \brief Overwrite surface output files or append iteration number.*/
  Wrt_Volume_Overwrite,               /*!< \brief Overwrite volume output files or append iteration number.*/
  Restart_Flow;                       /*!< \brief Restart flow solution for adjoint and linearized problems. */
  unsigned short nMarker_Monitoring,  /*!< \brief Number of markers to monitor. */
  nMarker_Designing,                  /*!< \brief Number of markers for the objective function. */
  nMarker_GeoEval,                    /*!< \brief Number of markers for the objective function. */
  nMarker_ZoneInterface,              /*!< \brief Number of markers in the zone interface. */
  nMarker_Plotting,                   /*!< \brief Number of markers to plot. */
  nMarker_Analyze,                    /*!< \brief Number of markers to analyze. */
  nMarker_Moving,                     /*!< \brief Number of markers in motion (DEFORMING, MOVING_WALL). */
  nMarker_PyCustom,                   /*!< \brief Number of markers that are customizable in Python. */
  nMarker_DV,                         /*!< \brief Number of markers affected by the design variables. */
  nMarker_WallFunctions,              /*!< \brief Number of markers for which wall functions must be applied. */
  nMarker_StrongBC,                   /*!< \brief Number of markers for which a strong BC must be applied. */
  nMarker_SobolevBC;                  /*!< \brief Number of markers treaded in the gradient problem. */
  string *Marker_Monitoring,          /*!< \brief Markers to monitor. */
  *Marker_Designing,                  /*!< \brief Markers to design. */
  *Marker_GeoEval,                    /*!< \brief Markers to evaluate geometry. */
  *Marker_Plotting,                   /*!< \brief Markers to plot. */
  *Marker_Analyze,                    /*!< \brief Markers to analyze. */
  *Marker_ZoneInterface,              /*!< \brief Markers in the FSI interface. */
  *Marker_Moving,                     /*!< \brief Markers in motion (DEFORMING, MOVING_WALL). */
  *Marker_PyCustom,                   /*!< \brief Markers that are customizable in Python. */
  *Marker_DV,                         /*!< \brief Markers affected by the design variables. */
  *Marker_WallFunctions,              /*!< \brief Markers for which wall functions must be applied. */
  *Marker_StrongBC,                   /*!< \brief Markers for which a strong BC must be applied. */
  *Marker_SobolevBC;                  /*!< \brief Markers in the gradient solver */

  unsigned short nConfig_Files;       /*!< \brief Number of config files for multiphysics problems. */
  string *Config_Filenames;           /*!< \brief List of names for configuration files. */
  SST_OPTIONS *SST_Options;           /*!< \brief List of modifications/corrections/versions of SST turbulence model.*/
  SA_OPTIONS *SA_Options;             /*!< \brief List of modifications/corrections/versions of SA turbulence model.*/
  LM_OPTIONS *LM_Options;             /*!< \brief List of modifications/corrections/versions of SA turbulence model.*/
  unsigned short nSST_Options;        /*!< \brief Number of SST options specified. */
  unsigned short nSA_Options;         /*!< \brief Number of SA options specified. */
  unsigned short nLM_Options;         /*!< \brief Number of SA options specified. */
  WALL_FUNCTIONS  *Kind_WallFunctions;        /*!< \brief The kind of wall function to use for the corresponding markers. */
  unsigned short  **IntInfo_WallFunctions;    /*!< \brief Additional integer information for the wall function markers. */
  su2double       **DoubleInfo_WallFunctions; /*!< \brief Additional double information for the wall function markers. */
  unsigned short  *Marker_All_Monitoring,     /*!< \brief Global index for monitoring using the grid information. */
  *Marker_All_GeoEval,               /*!< \brief Global index for geometrical evaluation. */
  *Marker_All_Plotting,              /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_Analyze,               /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_ZoneInterface,         /*!< \brief Global index for FSI interface markers using the grid information. */
  *Marker_All_Turbomachinery,        /*!< \brief Global index for Turbomachinery markers using the grid information. */
  *Marker_All_TurbomachineryFlag,    /*!< \brief Global index for Turbomachinery markers flag using the grid information. */
  *Marker_All_MixingPlaneInterface,  /*!< \brief Global index for MixingPlane interface markers using the grid information. */
  *Marker_All_DV,                    /*!< \brief Global index for design variable markers using the grid information. */
  *Marker_All_Moving,                /*!< \brief Global index for moving surfaces using the grid information. */
  *Marker_All_Deform_Mesh,           /*!< \brief Global index for deformable markers at the boundary. */
  *Marker_All_Deform_Mesh_Sym_Plane, /*!< \brief Global index for markers with symmetric deformations. */
  *Marker_All_Fluid_Load,            /*!< \brief Global index for markers in which the flow load is computed/employed. */
  *Marker_All_PyCustom,              /*!< \brief Global index for Python customizable surfaces using the grid information. */
  *Marker_All_Designing,             /*!< \brief Global index for moving using the grid information. */
  *Marker_All_SobolevBC,             /*!< \brief Global index for boundary condition applied to gradient smoothing. */
  *Marker_CfgFile_Monitoring,            /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Designing,             /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_GeoEval,               /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Plotting,              /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_Analyze,               /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_ZoneInterface,         /*!< \brief Global index for FSI interface using the config information. */
  *Marker_CfgFile_Turbomachinery,        /*!< \brief Global index for Turbomachinery  using the config information. */
  *Marker_CfgFile_TurbomachineryFlag,    /*!< \brief Global index for Turbomachinery flag using the config information. */
  *Marker_CfgFile_MixingPlaneInterface,  /*!< \brief Global index for MixingPlane interface using the config information. */
  *Marker_CfgFile_Moving,             /*!< \brief Global index for moving surfaces using the config information. */
  *Marker_CfgFile_Deform_Mesh,        /*!< \brief Global index for deformable markers at the boundary. */
  *Marker_CfgFile_Deform_Mesh_Sym_Plane, /*!< \brief Global index for markers with symmetric deformations. */
  *Marker_CfgFile_Fluid_Load,         /*!< \brief Global index for markers in which the flow load is computed/employed. */
  *Marker_CfgFile_PyCustom,           /*!< \brief Global index for Python customizable surfaces using the config information. */
  *Marker_CfgFile_DV,                 /*!< \brief Global index for design variable markers using the config information. */
  *Marker_CfgFile_PerBound,           /*!< \brief Global index for periodic boundaries using the config information. */
  *Marker_CfgFile_SobolevBC;          /*!< \brief Global index for boundary condition applied to gradient smoothing using the config information. */
  string *PlaneTag;                   /*!< \brief Global index for the plane adaptation (upper, lower). */
  su2double *nBlades;                 /*!< \brief number of blades for turbomachinery computation. */
  unsigned short Geo_Description;     /*!< \brief Description of the geometry. */
  unsigned short Mesh_FileFormat;     /*!< \brief Mesh input format. */
  TAB_OUTPUT Tab_FileFormat;          /*!< \brief Format of the output files. */
  unsigned short output_precision;    /*!< \brief <ofstream>.precision(value) for SU2_DOT and HISTORY output */
  unsigned short ActDisk_Jump;        /*!< \brief Format of the output files. */
  unsigned long StartWindowIteration; /*!< \brief Starting Iteration for long time Windowing apporach . */
  unsigned short nCFL_AdaptParam;     /*!< \brief Number of CFL parameters provided in config. */
  bool CFL_Adapt;        /*!< \brief Use adaptive CFL number. */
  bool HB_Precondition;  /*!< \brief Flag to turn on harmonic balance source term preconditioning */
  su2double RefArea,     /*!< \brief Reference area for coefficient computation. */
  RefElemLength,         /*!< \brief Reference element length for computing the slope limiting epsilon. */
  RefSharpEdges,         /*!< \brief Reference coefficient for detecting sharp edges. */
  RefLength,             /*!< \brief Reference length for moment computation. */
  *RefOriginMoment_X,    /*!< \brief X Origin for moment computation. */
  *RefOriginMoment_Y,    /*!< \brief Y Origin for moment computation. */
  *RefOriginMoment_Z,    /*!< \brief Z Origin for moment computation. */
  *CFL_AdaptParam,       /*!< \brief Information about the CFL ramp. */
  *RelaxFactor_Giles,    /*!< \brief Information about the under relaxation factor for Giles BC. */
  *CFL,                  /*!< \brief CFL number. */
  DomainVolume;          /*!< \brief Volume of the computational grid. */
  unsigned short
  nRefOriginMoment_X,      /*!< \brief Number of X-coordinate moment computation origins. */
  nRefOriginMoment_Y,      /*!< \brief Number of Y-coordinate moment computation origins. */
  nRefOriginMoment_Z;      /*!< \brief Number of Z-coordinate moment computation origins. */
  unsigned short nMesh_Box_Size;
  short *Mesh_Box_Size;          /*!< \brief Array containing the number of grid points in the x-, y-, and z-directions for the analytic RECTANGLE and BOX grid formats. */
  string Mesh_FileName,          /*!< \brief Mesh input file. */
  Mesh_Out_FileName,             /*!< \brief Mesh output file. */
  Solution_FileName,             /*!< \brief Flow solution input file. */
  Solution_AdjFileName,          /*!< \brief Adjoint solution input file for drag functional. */
  Volume_FileName,               /*!< \brief Flow variables output file. */
  Conv_FileName,                 /*!< \brief Convergence history output file. */
  Breakdown_FileName,            /*!< \brief Breakdown output file. */
  Restart_FileName,              /*!< \brief Restart file for flow variables. */
  Restart_AdjFileName,           /*!< \brief Restart file for adjoint variables, drag functional. */
  Adj_FileName,                  /*!< \brief Output file with the adjoint variables. */
  ObjFunc_Grad_FileName,         /*!< \brief Gradient of the objective function. */
  ObjFunc_Value_FileName,        /*!< \brief Objective function. */
  SurfCoeff_FileName,            /*!< \brief Output file with the flow variables on the surface. */
  SurfAdjCoeff_FileName,         /*!< \brief Output file with the adjoint variables on the surface. */
  SurfSens_FileName,             /*!< \brief Output file for the sensitivity on the surface (discrete adjoint). */
  VolSens_FileName,              /*!< \brief Output file for the sensitivity in the volume (discrete adjoint). */
  ObjFunc_Hess_FileName,         /*!< \brief Hessian approximation obtained by the Sobolev smoothing solver. */
  *DataDriven_Method_FileNames;    /*!< \brief Dataset information for data-driven fluid models. */

  bool
  Wrt_Performance,           /*!< \brief Write the performance summary at the end of a calculation.  */
  Wrt_AD_Statistics,         /*!< \brief Write the tape statistics (discrete adjoint).  */
  Wrt_MeshQuality,           /*!< \brief Write the mesh quality statistics to the visualization files.  */
  Wrt_MultiGrid,             /*!< \brief Write the coarse grids to the visualization files.  */
  Wrt_Projected_Sensitivity, /*!< \brief Write projected sensitivities (dJ/dx) on surfaces to ASCII file. */
  Plot_Section_Forces;       /*!< \brief Write sectional forces for specified markers. */
  unsigned short
  Console_Output_Verb,  /*!< \brief Level of verbosity for console output */
  Kind_Average;         /*!< \brief Particular average for the marker analyze. */
  su2double Gamma,      /*!< \brief Ratio of specific heats of the gas. */
  Bulk_Modulus,         /*!< \brief Value of the bulk modulus for incompressible flows. */
  Beta_Factor,          /*!< \brief Value of the epsilon^2 multiplier for Beta for the incompressible preconditioner. */
  Gas_Constant,         /*!< \brief Specific gas constant. */
  Gas_ConstantND,       /*!< \brief Non-dimensional specific gas constant. */
  *Molecular_Weight;    /*!< \brief Molecular weight of an incompressible ideal gas (g/mol). */
  unsigned short nMolecular_Weight, /*!< \brief Number of species molecular weights. */
  nSpecific_Heat_Cp;              /*!< \brief Number of species specific heat constants at constant pressure. */
  su2double *Specific_Heat_Cp, /*!< \brief Specific heat at constant pressure. */
  Thermal_Expansion_Coeff,    /*!< \brief Thermal expansion coefficient. */
  Thermal_Expansion_CoeffND,  /*!< \brief Non-dimensional thermal expansion coefficient. */
  Inc_Density_Ref,       /*!< \brief Reference density for custom incompressible non-dim. */
  Inc_Velocity_Ref,      /*!< \brief Reference velocity for custom incompressible non-dim. */
  Inc_Temperature_Ref,   /*!< \brief Reference temperature for custom incompressible non-dim. */
  Inc_Density_Init,      /*!< \brief Initial density for incompressible flows. */
  Inc_Temperature_Init,  /*!< \brief Initial temperature for incompressible flows w/ heat transfer. */
  Heat_Flux_Ref,         /*!< \brief Reference heat flux for non-dim. */
  Gas_Constant_Ref,      /*!< \brief Reference specific gas constant. */
  Temperature_Critical,  /*!< \brief Critical Temperature for real fluid model.  */
  Pressure_Critical,     /*!< \brief Critical Pressure for real fluid model.  */
  Density_Critical,      /*!< \brief Critical Density for real fluid model.  */
  Acentric_Factor,       /*!< \brief Acentric Factor for real fluid model.  */
  *Mu_Constant,           /*!< \brief Constant viscosity for ConstantViscosity model.  */
  *Thermal_Conductivity_Constant,  /*!< \brief Constant thermal conductivity for ConstantConductivity model.  */
  *Mu_Ref,                /*!< \brief Reference viscosity for Sutherland model.  */
  *Mu_Temperature_Ref,    /*!< \brief Reference temperature for Sutherland model.  */
  *Mu_S;                  /*!< \brief Reference S for Sutherland model.  */
  unsigned short nMu_Constant,   /*!< \brief Number of species constant viscosities. */
  nMu_Ref,                       /*!< \brief Number of species reference constants for Sutherland model. */
  nMu_Temperature_Ref,           /*!< \brief Number of species reference temperature for Sutherland model. */
  nMu_S,                         /*!< \brief Number of species reference S for Sutherland model. */
  nThermal_Conductivity_Constant,/*!< \brief Number of species constant thermal conductivity. */
  nPrandtl_Lam,                  /*!< \brief Number of species laminar Prandtl number. */
  nPrandtl_Turb,                 /*!< \brief Number of species turbulent Prandtl number. */
  nConstant_Lewis_Number;       /*!< \brief Number of species Lewis Number. */
  su2double Diffusivity_Constant;   /*!< \brief Constant mass diffusivity for scalar transport.  */
  su2double Diffusivity_ConstantND; /*!< \brief Non-dim. constant mass diffusivity for scalar transport.  */
  su2double Schmidt_Number_Laminar;   /*!< \brief Laminar Schmidt number for mass diffusion.  */
  su2double Schmidt_Number_Turbulent; /*!< \brief Turbulent Schmidt number for mass diffusion.  */
  su2double *Constant_Lewis_Number;   /*!< \brief Different Lewis number for mass diffusion.  */
  array<su2double, N_POLY_COEFFS> CpPolyCoefficientsND{{0.0}};  /*!< \brief Definition of the non-dimensional temperature polynomial coefficients for specific heat Cp. */
  array<su2double, N_POLY_COEFFS> MuPolyCoefficientsND{{0.0}};  /*!< \brief Definition of the non-dimensional temperature polynomial coefficients for viscosity. */
  array<su2double, N_POLY_COEFFS> KtPolyCoefficientsND{{0.0}};  /*!< \brief Definition of the non-dimensional temperature polynomial coefficients for thermal conductivity. */
  su2double TurbIntensityAndViscRatioFreeStream[2]; /*!< \brief Freestream turbulent intensity and viscosity ratio for turbulence and transition models. */
  su2double Energy_FreeStream,     /*!< \brief Free-stream total energy of the fluid.  */
  ModVel_FreeStream,               /*!< \brief Magnitude of the free-stream velocity of the fluid.  */
  ModVel_FreeStreamND,             /*!< \brief Non-dimensional magnitude of the free-stream velocity of the fluid.  */
  Density_FreeStream,              /*!< \brief Free-stream density of the fluid. */
  Viscosity_FreeStream,            /*!< \brief Free-stream viscosity of the fluid.  */
  Tke_FreeStream,                  /*!< \brief Total turbulent kinetic energy of the fluid.  */
  Intermittency_FreeStream,        /*!< \brief Freestream intermittency (for sagt transition model) of the fluid.  */
  ReThetaT_FreeStream,             /*!< \brief Freestream Transition Momentum Thickness Reynolds Number (for LM transition model) of the fluid.  */
  NuFactor_FreeStream,             /*!< \brief Ratio of turbulent to laminar viscosity. */
  NuFactor_Engine,                 /*!< \brief Ratio of turbulent to laminar viscosity at the engine. */
  SecondaryFlow_ActDisk,           /*!< \brief Ratio of turbulent to laminar viscosity at the actuator disk. */
  Initial_BCThrust,                /*!< \brief Ratio of turbulent to laminar viscosity at the actuator disk. */
  Pressure_FreeStream,             /*!< \brief Total pressure of the fluid. */
  Pressure_Thermodynamic,          /*!< \brief Thermodynamic pressure of the fluid. */
  Temperature_FreeStream,          /*!< \brief Total temperature of the fluid.  */
  Temperature_ve_FreeStream;       /*!< \brief Total vibrational-electronic temperature of the fluid.  */
  unsigned short wallModel_MaxIter; /*!< \brief maximum number of iterations for the Newton method for the wall model */
  su2double wallModel_Kappa,        /*!< \brief von Karman constant kappa for turbulence wall modeling */
  wallModel_B,                      /*!< \brief constant B for turbulence wall modeling */
  wallModel_RelFac,                 /*!< \brief relaxation factor for the Newton method used in the wall model */
  wallModel_MinYplus;               /*!< \brief minimum Y+ value, below which the wall model is not used anymore */
  su2double *Prandtl_Lam,      /*!< \brief Laminar Prandtl number for the gas.  */
  *Prandtl_Turb,               /*!< \brief Turbulent Prandtl number for the gas.  */
  Length_Ref,                 /*!< \brief Reference length for non-dimensionalization. */
  Pressure_Ref,               /*!< \brief Reference pressure for non-dimensionalization.  */
  Temperature_Ref,            /*!< \brief Reference temperature for non-dimensionalization.*/
  Temperature_ve_Ref,         /*!< \brief Reference vibrational-electronic temperature for non-dimensionalization.*/
  Density_Ref,                /*!< \brief Reference density for non-dimensionalization.*/
  Velocity_Ref,               /*!< \brief Reference velocity for non-dimensionalization.*/
  Time_Ref,                   /*!< \brief Reference time for non-dimensionalization. */
  Viscosity_Ref,              /*!< \brief Reference viscosity for non-dimensionalization. */
  Thermal_Conductivity_Ref,   /*!< \brief Reference conductivity for non-dimensionalization. */
  Energy_Ref,                 /*!< \brief Reference viscosity for non-dimensionalization. */
  Wall_Temperature,           /*!< \brief Temperature at an isotropic wall in Kelvin. */
  Omega_Ref,                  /*!< \brief Reference angular velocity for non-dimensionalization. */
  Force_Ref,                  /*!< \brief Reference body force for non-dimensionalization. */
  Pressure_FreeStreamND,      /*!< \brief Farfield pressure value (external flow). */
  Pressure_ThermodynamicND,   /*!< \brief Farfield thermodynamic pressure value. */
  Temperature_FreeStreamND,   /*!< \brief Farfield temperature value (external flow). */
  Temperature_ve_FreeStreamND,/*!< \brief Farfield vibrational-electronic temperature value (external flow). */
  Density_FreeStreamND,       /*!< \brief Farfield density value (external flow). */
  Velocity_FreeStreamND[3],   /*!< \brief Farfield velocity values (external flow). */
  Energy_FreeStreamND,        /*!< \brief Farfield energy value (external flow). */
  Viscosity_FreeStreamND,     /*!< \brief Farfield viscosity value (external flow). */
  Tke_FreeStreamND,           /*!< \brief Farfield kinetic energy (external flow). */
  Omega_FreeStreamND,         /*!< \brief Specific dissipation (external flow). */
  Omega_FreeStream;           /*!< \brief Specific dissipation (external flow). */
  bool Variable_Density;      /*!< \brief Variable density for incompressible flow. */
  unsigned short nElectric_Constant;    /*!< \brief Number of different electric constants. */
  su2double *Electric_Constant;         /*!< \brief Dielectric constant modulus. */
  su2double Knowles_B,                  /*!< \brief Knowles material model constant B. */
  Knowles_N;                            /*!< \brief Knowles material model constant N. */
  bool DE_Effects;                      /*!< Application of DE effects to FE analysis */
  bool RefGeom, RefGeomSurf;            /*!< Read a reference geometry for optimization purposes. */
  unsigned long refNodeID;              /*!< \brief Global ID for the reference node (optimization). */
  string RefGeom_FEMFileName;           /*!< \brief File name for reference geometry. */
  unsigned short RefGeom_FileFormat;    /*!< \brief Mesh input format. */
  STRUCT_2DFORM Kind_2DElasForm;        /*!< \brief Kind of bidimensional elasticity solver. */
  unsigned short nIterFSI_Ramp;         /*!< \brief Number of FSI subiterations during which a ramp is applied. */
  unsigned short iInst;                 /*!< \brief Current instance value */
  su2double AitkenStatRelax;      /*!< \brief Aitken's relaxation factor (if set as static) */
  su2double AitkenDynMaxInit;     /*!< \brief Aitken's maximum dynamic relaxation factor for the first iteration */
  su2double AitkenDynMinInit;     /*!< \brief Aitken's minimum dynamic relaxation factor for the first iteration */
  bool RampAndRelease;            /*!< \brief option for ramp load and release */
  bool Sine_Load;                 /*!< \brief option for sine load */
  su2double Thermal_Diffusivity;  /*!< \brief Thermal diffusivity used in the heat solver. */
  su2double Cyclic_Pitch,         /*!< \brief Cyclic pitch for rotorcraft simulations. */
  Collective_Pitch;               /*!< \brief Collective pitch for rotorcraft simulations. */
  su2double Mach_Motion;          /*!< \brief Mach number based on mesh velocity and freestream quantities. */

  su2double Motion_Origin[3] = {0.0}, /*!< \brief Mesh motion origin. */
  Translation_Rate[3] = {0.0},        /*!< \brief Translational velocity of the mesh. */
  Rotation_Rate[3] = {0.0},           /*!< \brief Angular velocity of the mesh . */
  Pitching_Omega[3] = {0.0},          /*!< \brief Angular frequency of the mesh pitching. */
  Pitching_Ampl[3] = {0.0},           /*!< \brief Pitching amplitude. */
  Pitching_Phase[3] = {0.0},          /*!< \brief Pitching phase offset. */
  Plunging_Omega[3] = {0.0},          /*!< \brief Angular frequency of the mesh plunging. */
  Plunging_Ampl[3] = {0.0};           /*!< \brief Plunging amplitude. */
  su2double *MarkerMotion_Origin, /*!< \brief Mesh motion origin of marker. */
  *MarkerTranslation_Rate,        /*!< \brief Translational velocity of marker. */
  *MarkerRotation_Rate,           /*!< \brief Angular velocity of marker. */
  *MarkerPitching_Omega,          /*!< \brief Angular frequency of marker. */
  *MarkerPitching_Ampl,           /*!< \brief Pitching amplitude of marker. */
  *MarkerPitching_Phase,          /*!< \brief Pitching phase offset of marker. */
  *MarkerPlunging_Omega,          /*!< \brief Angular frequency of marker.. */
  *MarkerPlunging_Ampl;           /*!< \brief Plunging amplitude of marker. */

  unsigned short
  nMarkerMotion_Origin,           /*!< \brief Number of values provided for mesh motion origin of marker. */
  nMarkerTranslation,             /*!< \brief Number of values provided for translational velocity of marker. */
  nMarkerRotation_Rate,           /*!< \brief Number of values provided for angular velocity of marker. */
  nMarkerPitching_Omega,          /*!< \brief Number of values provided for angular frequency of marker. */
  nMarkerPitching_Ampl,           /*!< \brief Number of values provided for pitching amplitude of marker. */
  nMarkerPitching_Phase,          /*!< \brief Number of values provided for pitching phase offset of marker. */
  nMarkerPlunging_Omega,          /*!< \brief Number of values provided for angular frequency of marker. */
  nMarkerPlunging_Ampl,           /*!< \brief Number of values provided for plunging amplitude of marker. */
  nRough_Wall;                    /*!< \brief Number of rough walls. */
  su2double  *Omega_HB;           /*!< \brief Frequency for Harmonic Balance Operator (in rad/s). */
  unsigned short
  nOmega_HB,                      /*!< \brief Number of frequencies in Harmonic Balance Operator. */
  nMoveMotion_Origin,             /*!< \brief Number of motion origins. */
  *MoveMotion_Origin;             /*!< \brief Keeps track if we should move moment origin. */
  vector<vector<vector<su2double> > > Aeroelastic_np1, /*!< \brief Aeroelastic solution at time level n+1. */
  Aeroelastic_n,                  /*!< \brief Aeroelastic solution at time level n. */
  Aeroelastic_n1;                 /*!< \brief Aeroelastic solution at time level n-1. */
  su2double FlutterSpeedIndex,    /*!< \brief The flutter speed index. */
  PlungeNaturalFrequency,         /*!< \brief Plunging natural frequency for Aeroelastic. */
  PitchNaturalFrequency,          /*!< \brief Pitch natural frequency for Aeroelastic. */
  AirfoilMassRatio,               /*!< \brief The airfoil mass ratio for Aeroelastic. */
  CG_Location,                    /*!< \brief Center of gravity location for Aeroelastic. */
  RadiusGyrationSquared;          /*!< \brief The radius of gyration squared for Aeroelastic. */
  su2double *Aeroelastic_plunge,  /*!< \brief Value of plunging coordinate at the end of an external iteration. */
  *Aeroelastic_pitch;             /*!< \brief Value of pitching coordinate at the end of an external iteration. */
  unsigned short AeroelasticIter; /*!< \brief Solve the aeroelastic equations every given number of internal iterations. */
  unsigned short Gust_Type,   /*!< \brief Type of Gust. */
  Gust_Dir;                   /*!< \brief Direction of the gust */
  su2double Gust_WaveLength,  /*!< \brief The gust wavelength. */
  Gust_Periods,               /*!< \brief Number of gust periods. */
  Gust_Ampl,                  /*!< \brief Gust amplitude. */
  Gust_Begin_Time,            /*!< \brief Time at which to begin the gust. */
  Gust_Begin_Loc;             /*!< \brief Location at which the gust begins. */
  /*! \brief Maximal scalar product of the normed far-field velocity vector and a space coordinate where fixed turbulence quantities are set. */
  su2double Turb_Fixed_Values_MaxScalarProd;
  long Visualize_CV;          /*!< \brief Node number for the CV to be visualized */
  bool ExtraOutput;           /*!< \brief Check if extra output need. */
  bool Wall_Functions;           /*!< \brief Use wall functions with the turbulence model */
  long ExtraHeatOutputZone;      /*!< \brief Heat solver zone with extra screen output */
  bool DeadLoad;                 /*!< \brief Application of dead loads to the FE analysis */
  bool PseudoStatic;             /*!< \brief Application of dead loads to the FE analysis */
  bool SteadyRestart;            /*!< \brief Restart from a steady state for FSI problems. */
  su2double Newmark_beta,        /*!< \brief Parameter alpha for Newmark method. */
  Newmark_gamma;                 /*!< \brief Parameter delta for Newmark method. */
  unsigned short nIntCoeffs;     /*!< \brief Number of integration coeffs for structural calculations. */
  su2double *Int_Coeffs;         /*!< \brief Time integration coefficients for structural method. */
  unsigned short nElasticityMod, /*!< \brief Number of different values for the elasticity modulus. */
  nPoissonRatio,                    /*!< \brief Number of different values for the Poisson ratio modulus. */
  nMaterialDensity;                 /*!< \brief Number of different values for the Material density. */
  su2double *ElasticityMod,         /*!< \brief Value of the elasticity moduli. */
  *PoissonRatio,                    /*!< \brief Value of the Poisson ratios. */
  *MaterialDensity;                 /*!< \brief Value of the Material densities. */
  unsigned short nElectric_Field,   /*!< \brief Number of different values for the electric field in the membrane. */
  nDim_Electric_Field;              /*!< \brief Dimensionality of the problem. */
  unsigned short nDim_RefNode;      /*!< \brief Dimensionality of the vector . */
  su2double *Electric_Field_Mod,    /*!< \brief Values of the modulus of the electric field. */
  *Electric_Field_Dir;              /*!< \brief Direction of the electric field. */
  su2double *RefNode_Displacement;  /*!< \brief Displacement of the reference node. */
  bool Ramp_Load;                         /*!< \brief Apply the load with linear increases. */
  unsigned short Dynamic_LoadTransfer;    /*!< \brief Method for dynamic load transferring. */
  bool IncrementalLoad;                   /*!< \brief Apply the load in increments (for nonlinear structural analysis). */
  unsigned long IncLoad_Nincrements;      /*!< \brief Number of increments. */
  su2double Ramp_Time;                    /*!< \brief Time until the maximum load is applied. */
  bool Predictor,                         /*!< \brief Determines whether a predictor step is used. */
  Relaxation;                             /*!< \brief Determines whether a relaxation step is used. */
  unsigned short Pred_Order;              /*!< \brief Order of the predictor for FSI applications. */
  INTERFACE_INTERPOLATOR Kind_Interpolation; /*!< \brief type of interpolation to use for FSI applications. */
  bool ConservativeInterpolation;            /*!< \brief Conservative approach for non matching mesh interpolation. */
  unsigned short NumNearestNeighbors;        /*!< \brief Number of neighbors used for Nearest Neighbor interpolation. */
  RADIAL_BASIS Kind_RadialBasisFunction;     /*!< \brief type of radial basis function to use for radial basis FSI. */
  bool RadialBasisFunction_PolynomialOption; /*!< \brief Option of whether to include polynomial terms in Radial Basis Function Interpolation or not. */
  su2double RadialBasisFunction_Parameter;   /*!< \brief Radial basis function parameter (radius). */
  su2double RadialBasisFunction_PruneTol;    /*!< \brief Tolerance to prune the RBF interpolation matrix. */
  bool Prestretch;                           /*!< \brief Read a reference geometry for optimization purposes. */
  string Prestretch_FEMFileName;             /*!< \brief File name for reference geometry. */
  string FEA_FileName;              /*!< \brief File name for element-based properties. */
  bool FEAAdvancedMode;             /*!< \brief Determine if advanced features are used from the element-based FEA analysis (experimental). */
  su2double RefGeom_Penalty,        /*!< \brief Penalty weight value for the reference geometry objective function. */
  RefNode_Penalty,                  /*!< \brief Penalty weight value for the reference node objective function. */
  DV_Penalty;                       /*!< \brief Penalty weight to add a constraint to the total amount of stiffness. */
  array<su2double,2> StressPenaltyParam = {{1.0, 20.0}}; /*!< \brief Allowed stress and KS aggregation exponent. */
  unsigned long Nonphys_Points,     /*!< \brief Current number of non-physical points in the solution. */
  Nonphys_Reconstr;                 /*!< \brief Current number of non-physical reconstructions for 2nd-order upwinding. */
  su2double ParMETIS_tolerance;     /*!< \brief Load balancing tolerance for ParMETIS. */
  long ParMETIS_pointWgt;           /*!< \brief Load balancing weight given to points. */
  long ParMETIS_edgeWgt;            /*!< \brief Load balancing weight given to edges. */
  unsigned short DirectDiff;        /*!< \brief Direct Differentation mode. */
  bool DiscreteAdjoint;                /*!< \brief AD-based discrete adjoint mode. */
  su2double Const_DES;                 /*!< \brief Detached Eddy Simulation Constant. */
  WINDOW_FUNCTION Kind_WindowFct;      /*!< \brief Type of window (weight) function for objective functional. */
  unsigned short Kind_HybridRANSLES;   /*!< \brief Kind of Hybrid RANS/LES. */
  unsigned short Kind_RoeLowDiss;      /*!< \brief Kind of Roe scheme with low dissipation for unsteady flows. */

  unsigned short nSpanWiseSections; /*!< \brief number of span-wise sections */
  unsigned short nSpanMaxAllZones;  /*!< \brief number of maximum span-wise sections for all zones */
  unsigned short *nSpan_iZones;     /*!< \brief number of span-wise sections for each zones */
  bool turbMixingPlane;             /*!< \brief option for turbulent mixingplane */
  bool SpatialFourier;              /*!< \brief option for computing the fourier transforms for subsonic non-reflecting BC. */
  bool RampRotatingFrame;           /*!< \brief option for ramping up or down the Rotating Frame values */
  bool RampOutletPressure;          /*!< \brief option for ramping up or down the outlet pressure */
  su2double AverageMachLimit;           /*!< \brief option for turbulent mixingplane */
  su2double FinalRotation_Rate_Z;       /*!< \brief Final rotation rate Z if Ramp rotating frame is activated. */
  su2double FinalOutletPressure;        /*!< \brief Final outlet pressure if Ramp outlet pressure is activated. */
  su2double MonitorOutletPressure;      /*!< \brief Monitor outlet pressure if Ramp outlet pressure is activated. */
  array<su2double, N_POLY_COEFFS> cp_polycoeffs{{0.0}};  /*!< \brief Array for specific heat polynomial coefficients. */
  array<su2double, N_POLY_COEFFS> mu_polycoeffs{{0.0}};  /*!< \brief Array for viscosity polynomial coefficients. */
  array<su2double, N_POLY_COEFFS> kt_polycoeffs{{0.0}};  /*!< \brief Array for thermal conductivity polynomial coefficients. */
  bool Body_Force;                      /*!< \brief Flag to know if a body force is included in the formulation. */

  ENUM_STREAMWISE_PERIODIC Kind_Streamwise_Periodic; /*!< \brief Kind of Streamwise periodic flow (pressure drop or massflow) */
  bool Streamwise_Periodic_Temperature;              /*!< \brief Use real periodicity for Energy equation or otherwise outlet source term. */
  su2double Streamwise_Periodic_PressureDrop;        /*!< \brief Value of prescribed pressure drop [Pa] which results in an artificial body force vector. */
  su2double Streamwise_Periodic_TargetMassFlow;      /*!< \brief Value of prescribed massflow [kg/s] which results in an delta p and therefore an artificial body force vector. */
  su2double Streamwise_Periodic_OutletHeat;          /*!< /brief Heatflux boundary [W/m^2] imposed at streamwise periodic outlet. */

  su2double *FreeStreamTurboNormal;     /*!< \brief Direction to initialize the flow in turbomachinery computation */
  su2double Restart_Bandwidth_Agg;      /*!< \brief The aggregate of the bandwidth for writing binary restarts (to be averaged later). */
  su2double Max_Vel2;                   /*!< \brief The maximum velocity^2 in the domain for the incompressible preconditioner. */
  bool topology_optimization;           /*!< \brief If the structural solver should consider a variable density field to penalize element stiffness. */
  string top_optim_output_file;         /*!< \brief File to where the derivatives w.r.t. element densities will be written to. */
  su2double simp_exponent;              /*!< \brief Exponent for the density-based stiffness penalization of the SIMP method. */
  su2double simp_minimum_stiffness;     /*!< \brief Lower bound for the stiffness penalization of the SIMP method. */
  ENUM_FILTER_KERNEL* top_optim_kernels;   /*!< \brief The kernels to use. */
  unsigned short top_optim_nKernel;        /*!< \brief Number of kernels specified. */
  unsigned short top_optim_nKernelParams;  /*!< \brief Number of kernel parameters specified. */
  unsigned short top_optim_nRadius;        /*!< \brief Number of radius values specified. */
  unsigned short top_optim_search_lim;     /*!< \brief Limit the maximum "logical radius" considered during filtering. */
  su2double *top_optim_kernel_params;  /*!< \brief The kernel parameters. */
  su2double *top_optim_filter_radius;  /*!< \brief Radius of the filter(s) used on the design density for topology optimization. */
  ENUM_PROJECTION_FUNCTION top_optim_proj_type;  /*!< \brief The projection function used in topology optimization. */
  su2double top_optim_proj_param;      /*!< \brief The value of the parameter for the projection function. */
  bool HeatSource;                     /*!< \brief Flag to know if there is a volumetric heat source on the flow. */
  su2double ValHeatSource;             /*!< \brief Value of the volumetric heat source on the flow (W/m3). */
  su2double Heat_Source_Rot_Z;         /*!< \brief Rotation of the volumetric heat source on the Z axis. */
  RADIATION_MODEL Kind_Radiation;      /*!< \brief Kind of radiation model used. */
  P1_INIT Kind_P1_Init;                /*!< \brief Kind of initialization used in the P1 model. */
  su2double Absorption_Coeff,          /*!< \brief Absorption coefficient of the medium (radiation). */
  Scattering_Coeff;                    /*!< \brief Scattering coefficient of the medium (radiation). */
  unsigned short nMarker_Emissivity;   /*!< \brief Number of markers for which the emissivity is defined. */
  string *Marker_Emissivity;           /*!< \brief Wall markers with defined emissivity. */
  su2double *Wall_Emissivity;          /*!< \brief Emissivity of the wall. */
  bool Radiation;                      /*!< \brief Determines if a radiation model is incorporated. */
  su2double CFL_Rad;                   /*!< \brief CFL Number for the radiation solver. */

  array<su2double,5> default_cfl_adapt;  /*!< \brief Default CFL adapt param array for the COption class. */
  su2double vel_init[3], /*!< \brief initial velocity array for the COption class. */
  vel_inf[3],            /*!< \brief freestream velocity array for the COption class. */
  eng_cyl[7],            /*!< \brief engine box array for the COption class. */
  eng_val[5],            /*!< \brief engine box array values for the COption class. */
  jst_coeff[2],          /*!< \brief artificial dissipation (flow) array for the COption class. */
  ffd_coeff[3],          /*!< \brief artificial dissipation (flow) array for the COption class. */
  mixedout_coeff[3],     /*!< \brief default mixedout algorithm coefficients for the COption class. */
  rampRotFrame_coeff[3], /*!< \brief ramp rotating frame coefficients for the COption class. */
  rampOutPres_coeff[3],  /*!< \brief ramp outlet pressure coefficients for the COption class. */
  jst_adj_coeff[2],      /*!< \brief artificial dissipation (adjoint) array for the COption class. */
  mesh_box_length[3],    /*!< \brief mesh box length for the COption class. */
  mesh_box_offset[3],    /*!< \brief mesh box offset for the COption class. */
  geo_loc[2],            /*!< \brief SU2_GEO section locations array for the COption class. */
  distortion[2],         /*!< \brief SU2_GEO section locations array for the COption class. */
  ea_lim[3],             /*!< \brief equivalent area limit array for the COption class. */
  grid_fix[6],           /*!< \brief fixed grid (non-deforming region) array for the COption class. */
  htp_axis[2],           /*!< \brief HTP axis for the COption class. */
  ffd_axis[3],           /*!< \brief FFD axis for the COption class. */
  inc_crit[3],           /*!< \brief incremental criteria array for the COption class. */
  extrarelfac[2],        /*!< \brief extra relaxation factor for Giles BC in the COption class. */
  sineload_coeff[3],     /*!< \brief values for a sine load. */
  body_force[3],         /*!< \brief body force vector for the COption class. */
  nacelle_location[5],   /*!< \brief Location of the nacelle. */
  hs_axes[3],            /*!< \brief principal axes (x, y, z) of the ellipsoid containing the heat source. */
  hs_center[3];          /*!< \brief position of the center of the heat source. */

  UPWIND Riemann_Solver_FEM;         /*!< \brief Riemann solver chosen for the DG method. */
  su2double Quadrature_Factor_Straight;      /*!< \brief Factor applied during quadrature of elements with a constant Jacobian. */
  su2double Quadrature_Factor_Curved;        /*!< \brief Factor applied during quadrature of elements with a non-constant Jacobian. */
  su2double Quadrature_Factor_Time_ADER_DG;  /*!< \brief Factor applied during quadrature in time for ADER-DG. */
  su2double Theta_Interior_Penalty_DGFEM;    /*!< \brief Factor for the symmetrizing terms in the DG discretization of the viscous fluxes. */
  unsigned short byteAlignmentMatMul;        /*!< \brief Number of bytes in the vectorization direction for the matrix multiplication. Multipe of 64. */
  unsigned short sizeMatMulPadding;          /*!< \brief The matrix size in the vectorization direction padded to a multiple of 8. Computed from byteAlignmentMatMul. */
  bool Compute_Entropy;                      /*!< \brief Whether or not to compute the entropy in the fluid model. */
  bool Use_Lumped_MassMatrix_DGFEM;          /*!< \brief Whether or not to use the lumped mass matrix for DGFEM. */
  bool Jacobian_Spatial_Discretization_Only; /*!< \brief Flag to know if only the exact Jacobian of the spatial discretization must be computed. */
  bool Compute_Average;                      /*!< \brief Whether or not to compute averages for unsteady simulations in FV or DG solver. */
  unsigned short Comm_Level;                 /*!< \brief Level of MPI communications to be performed. */
  VERIFICATION_SOLUTION Kind_Verification_Solution; /*!< \brief Verification solution for accuracy assessment. */

  bool Time_Domain;              /*!< \brief Determines if the multizone problem is solved in time-domain */
  unsigned long nOuterIter,      /*!< \brief Determines the number of outer iterations in the multizone problem */
  nInnerIter,                    /*!< \brief Determines the number of inner iterations in each multizone block */
  nTimeIter,                     /*!< \brief Determines the number of time iterations in the multizone problem */
  nIter,                         /*!< \brief Determines the number of pseudo-time iterations in a single-zone problem */
  Restart_Iter;                  /*!< \brief Determines the restart iteration in the multizone problem */
  su2double Time_Step;           /*!< \brief Determines the time step for the multizone problem */
  su2double Max_Time;            /*!< \brief Determines the maximum time for the time-domain problems */

  unsigned long HistoryWrtFreq[3],    /*!< \brief Array containing history writing frequencies for timer iter, outer iter, inner iter */
                ScreenWrtFreq[3];     /*!< \brief Array containing screen writing frequencies for timer iter, outer iter, inner iter */
  OUTPUT_TYPE* VolumeOutputFiles;     /*!< \brief File formats to output */
  unsigned short nVolumeOutputFiles=0;/*!< \brief Number of File formats to output */
  unsigned short nVolumeOutputFrequencies; /*!< \brief Number of frequencies for the volume outputs */
  unsigned long *VolumeOutputFrequencies; /*!< \brief list containing the writing frequencies */

  bool Multizone_Mesh;            /*!< \brief Determines if the mesh contains multiple zones. */
  bool Wrt_ZoneConv;              /*!< \brief Write the convergence history of each individual zone to screen. */
  bool Wrt_ZoneHist;              /*!< \brief Write the convergence history of each individual zone to file. */
  bool SpecialOutput,             /*!< \brief Determines if the special output is written. */
  Wrt_ForcesBreakdown;            /*!< \brief Determines if the forces breakdown file is written. */
  string *ScreenOutput,           /*!< \brief Kind of the screen output. */
  *HistoryOutput, *VolumeOutput;  /*!< \brief Kind of the output printed to the history file. */
  unsigned short nScreenOutput,   /*!< \brief Number of screen output variables (max: 6). */
  nHistoryOutput, nVolumeOutput;  /*!< \brief Number of variables printed to the history file. */
  bool Multizone_Residual;        /*!< \brief Determines if memory should be allocated for the multizone residual. */
  SST_ParsedOptions sstParsedOptions; /*!< \brief Additional parameters for the SST turbulence model. */
  SA_ParsedOptions saParsedOptions;   /*!< \brief Additional parameters for the SA turbulence model. */
  LM_ParsedOptions lmParsedOptions;   /*!< \brief Additional parameters for the LM transition model. */
  su2double uq_delta_b;         /*!< \brief Parameter used to perturb eigenvalues of Reynolds Stress Matrix */
  unsigned short eig_val_comp;  /*!< \brief Parameter used to determine type of eigenvalue perturbation */
  su2double uq_urlx;            /*!< \brief Under-relaxation factor */
  bool uq_permute;              /*!< \brief Permutation of eigenvectors */

  unsigned long pastix_fact_freq;  /*!< \brief (Re-)Factorization frequency for PaStiX */
  unsigned short pastix_verb_lvl;  /*!< \brief Verbosity level for PaStiX */
  unsigned short pastix_fill_lvl;  /*!< \brief Fill level for PaStiX ILU */

  string caseName;                 /*!< \brief Name of the current case */

  unsigned long edgeColorGroupSize; /*!< \brief Size of the edge groups colored for OpenMP parallelization of edge loops. */

  INLET_SPANWISE_INTERP Kind_InletInterpolationFunction; /*!brief type of spanwise interpolation function to use for the inlet face. */
  INLET_INTERP_TYPE Kind_Inlet_InterpolationType;    /*!brief type of spanwise interpolation data to use for the inlet face. */
  bool PrintInlet_InterpolatedData;               /*!brief option for printing the interpolated data file. */

  /*--- libROM configure options ---*/
  bool libROM;                              /*!< \brief Toggle saving to libROM. */
  string libROMbase_FileName;               /*!< \brief Base filename for libROM file saving. */
  POD_KIND POD_Basis_Gen;                   /*!< \brief Type of POD basis generation (static or incremental). */
  unsigned short maxBasisDim,               /*!< \brief Maximum number of POD basis dimensions. */
  rom_save_freq;                            /*!< \brief Frequency of unsteady time steps to save. */

  unsigned short nSpecies = 0;              /*!< \brief Number of transported species equations (for NEMO and species transport)*/

  /* other NEMO configure options*/
  unsigned short nSpecies_Cat_Wall,         /*!< \brief No. of species for a catalytic wall. */
  nSpecies_inlet,                           /*!< \brief No. of species for NEMO inlet. */
  iWall_Catalytic,                          /*!< \brief Iterator over catalytic walls. */
  nWall_Catalytic;                          /*!< \brief No. of catalytic walls. */
  su2double *Gas_Composition,               /*!< \brief Initial mass fractions of flow [dimensionless]. */
  *Supercatalytic_Wall_Composition,         /*!< \brief Supercatalytic wall mass fractions [dimensionless]. */
  pnorm_heat;                               /*!< \brief pnorm for heat-flux. */
  bool frozen,                              /*!< \brief Flag for determining if mixture is frozen. */
  ionization,                               /*!< \brief Flag for determining if free electron gas is in the mixture. */
  vt_transfer_res_limit,                    /*!< \brief Flag for determining if residual limiting for source term VT-transfer is used. */
  monoatomic,                               /*!< \brief Flag for monoatomic mixture. */
  Supercatalytic_Wall;                      /*!< \brief Flag for supercatalytic wall. */
  string GasModel,                          /*!< \brief Gas Model. */
  *Wall_Catalytic;                          /*!< \brief Pointer to catalytic walls. */
  TRANSCOEFFMODEL   Kind_TransCoeffModel;   /*!< \brief Transport coefficient Model for NEMO solver. */
  su2double CatalyticEfficiency;            /*!< \brief Wall catalytic efficiency. */
  su2double *Inlet_MassFrac;                /*!< \brief Specified Mass fraction vectors for NEMO inlet boundaries. */
  su2double Inlet_Temperature_ve;           /*!< \brief Specified Tve for supersonic inlet boundaries (NEMO solver). */

  /*--- Additional species solver options ---*/
  bool Species_Clipping;           /*!< \brief Boolean that activates solution clipping for scalar transport. */
  su2double* Species_Clipping_Max; /*!< \brief Maximum value of clipping for scalar transport. */
  su2double* Species_Clipping_Min; /*!< \brief Minimum value of clipping for scalar transport. */
  unsigned short nSpecies_Clipping_Max, nSpecies_Clipping_Min; /*!< \brief Number of entries of SPECIES_CLIPPING_MIN/MAX */
  bool Species_StrongBC;           /*!< \brief Boolean whether strong BC's are used for in- outlet of the species solver. */
  su2double* Species_Init;         /*!< \brief Initial uniform value for scalar transport. */
  unsigned short nSpecies_Init;    /*!< \brief Number of entries of SPECIES_INIT */

  /*--- Additional flamelet solver options ---*/
  su2double flame_init[8];       /*!< \brief Initial solution parameters for flamelet solver.*/

  /*--- lookup table ---*/
  unsigned short n_scalars = 0;       /*!< \brief Number of transported scalars for flamelet LUT approach. */
  unsigned short n_lookups = 0;       /*!< \brief Number of lookup variables, for visualization only. */
  unsigned short n_table_sources = 0; /*!< \brief Number of transported scalar source terms for LUT. */
  unsigned short n_user_scalars = 0;  /*!< \brief Number of user defined (auxiliary) scalar transport equations. */
  unsigned short n_user_sources = 0;  /*!< \brief Number of source terms for user defined (auxiliary) scalar transport equations. */
  unsigned short n_control_vars = 0;  /*!< \brief Number of controlling variables (independent variables) for the LUT. */

  string* controlling_variable_names;
  string* cv_source_names;
  vector<string> table_scalar_names;  /*!< \brief Names of transported scalar variables. */
  string* lookup_names;         /*!< \brief Names of passive look-up variables. */
  string* user_scalar_names;          /*!< \brief Names of the user defined (auxiliary) transported scalars .*/
  string* user_source_names;          /*!< \brief Names of the source terms for the user defined transported scalars. */

  /*!
   * \brief Set the default values of config options not set in the config file using another config object.
   * \param config - Config object to use the default values from.
   */
  void SetDefaultFromConfig(CConfig *config);

  /*!
   * \brief Set default values for all options not yet set.
   */
  void SetDefault();

  /*--- all_options is a map containing all of the options. This is used during config file parsing
   to track the options which have not been set (so the default values can be used). Without this map
   there would be no list of all the config file options. ---*/

  map<string, bool> all_options;

  /*--- brief param is a map from the option name (config file string) to its decoder (the specific child
   class of COptionBase that turns the string into a value) ---*/

  map<string, COptionBase*> option_map;


  // All of the addXxxOptions take in the name of the option, and a reference to the field of that option
  // in the option structure. Depending on the specific type, it may take in a default value, and may
  // take in extra options. The addXxxOptions mostly follow the same pattern, so please see addDoubleOption
  // for detailed comments.
  //
  // List options are those that can be an unknown number of elements, and also take in a reference to
  // an integer. This integer will be populated with the number of elements of that type unmarshaled.
  //
  // Array options are those with a fixed number of elements.
  //
  // List and Array options should also be able to be specified with the string "NONE" indicating that there
  // are no elements. This allows the option to be present in a config file but left blank.

  /*!< \brief addDoubleOption creates a config file parser for an option with the given name whose
   value can be represented by a su2double.*/

  void addDoubleOption(const string& name, su2double & option_field, su2double default_value);

  void addStringOption(const string& name, string & option_field, string default_value);

  void addIntegerOption(const string& name, int & option_field, int default_value);

  void addUnsignedLongOption(const string& name, unsigned long & option_field, unsigned long default_value);

  void addUnsignedShortOption(const string& name, unsigned short & option_field, unsigned short default_value);

  void addLongOption(const string& name, long & option_field, long default_value);

  void addBoolOption(const string& name, bool & option_field, bool default_value);

  // enum types work differently than all of the others because there are a small number of valid
  // string entries for the type. One must also provide a list of all the valid strings of that type.
  template <class Tenum, class Tfield>
  void addEnumOption(const string name, Tfield& option_field, const map<string,Tenum>& enum_map, Tenum default_value);

  // input_size is the number of options read in from the config file
  template <class Tenum, class Tfield>
  void addEnumListOption(const string name, unsigned short& input_size, Tfield*& option_field, const map<string,Tenum>& enum_map);

  void addDoubleArrayOption(const string& name, const int size, su2double* option_field);

  void addUShortArrayOption(const string& name, const int size, unsigned short* option_field);

  void addDoubleListOption(const string& name, unsigned short & size, su2double * & option_field);

  void addShortListOption(const string& name, unsigned short & size, short * & option_field);

  void addUShortListOption(const string& name, unsigned short & size, unsigned short * & option_field);

  void addULongListOption(const string& name, unsigned short & size, unsigned long * & option_field);

  void addStringListOption(const string& name, unsigned short & num_marker, string* & option_field);

  void addConvectOption(const string& name, unsigned short & space_field, CENTERED & centered_field, UPWIND & upwind_field);

  void addConvectFEMOption(const string& name, unsigned short & space_field, unsigned short & fem_field);

  void addMathProblemOption(const string& name, bool & ContinuousAdjoint, const bool & ContinuousAdjoint_default,
                            bool & DiscreteAdjoint, const bool & DiscreteAdjoint_default,
                            bool & Restart_Flow, const bool & Restart_Flow_default);

  void addDVParamOption(const string& name, unsigned short & nDV_field, su2double** & paramDV, string* & FFDTag,
                        unsigned short* & design_variable);

  void addDVValueOption(const string& name, unsigned short* & nDVValue_field, su2double** & valueDV, unsigned short & nDV_field,  su2double** & paramDV,
                        unsigned short* & design_variable);

  void addFFDDefOption(const string& name, unsigned short & nFFD_field, su2double** & coordFFD, string* & FFDTag);

  void addFFDDegreeOption(const string& name, unsigned short & nFFD_field, unsigned short** & degreeFFD);

  void addStringDoubleListOption(const string& name, unsigned short & list_size, string * & string_field,
                                 su2double* & double_field);

  void addInletOption(const string& name, unsigned short & nMarker_Inlet, string * & Marker_Inlet,
                      su2double* & Ttotal, su2double* & Ptotal, su2double** & FlowDir);

  void addInletSpeciesOption(const string& name, unsigned short & nMarker_Inlet_Species, string * & Marker_Inlet_Species,
                             su2double** & inlet_species_val, unsigned short & nSpecies_per_Inlet);

  void addInletTurbOption(const string& name, unsigned short& nMarker_Inlet_Turb, string*& Marker_Inlet_Turb,
                          su2double** & Turb_Properties, unsigned short & nTurb_Properties);

  template <class Tenum>
  void addRiemannOption(const string name, unsigned short & nMarker_Riemann, string * & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                        su2double* & var1, su2double* & var2, su2double** & FlowDir);

  template <class Tenum>
  void addGilesOption(const string name, unsigned short & nMarker_Giles, string * & Marker_Giles, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                     su2double* & var1, su2double* & var2, su2double** & FlowDir, su2double* & relaxfactor1, su2double* & relaxfactor2);

  void addExhaustOption(const string& name, unsigned short & nMarker_Exhaust, string * & Marker_Exhaust,
                        su2double* & Ttotal, su2double* & Ptotal);

  void addPeriodicOption(const string & name, unsigned short & nMarker_PerBound,
                         string* & Marker_PerBound, string* & Marker_PerDonor,
                         su2double** & RotCenter, su2double** & RotAngles, su2double** & Translation);

  void addTurboPerfOption(const string & name, unsigned short & nMarker_TurboPerf,
                          string* & Marker_TurboBoundIn, string* & Marker_TurboBoundOut);

  void addActDiskOption(const string & name,
                        unsigned short & nMarker_ActDiskInlet, unsigned short & nMarker_ActDiskOutlet, string* & Marker_ActDiskInlet, string* & Marker_ActDiskOutlet,
                        su2double** & ActDisk_PressJump, su2double** & ActDisk_TempJump, su2double** & ActDisk_Omega);

  void addWallFunctionOption(const string &name,               unsigned short &list_size,
                             string* &string_field,            WALL_FUNCTIONS* &val_Kind_WF,
                             unsigned short** &val_IntInfo_WF, su2double** &val_DoubleInfo_WF);

  void addPythonOption(const string& name);

public:

  /*!
   * \brief Tags for the different fields in a restart file.
   */
  vector<string> fields;

  /*!
   * \brief Constructor of the class which reads the input file.
   */
  CConfig(char case_filename[MAX_STRING_SIZE], SU2_COMPONENT val_software, bool verb_high);

  /*!
   * \brief Constructor of the class which takes an istream buffer containing the config options.
   */
  CConfig(istream &case_buffer, SU2_COMPONENT val_software, bool verb_high);

  /*!
   * \brief Constructor of the class which reads the input file and uses default options from another config.
   */
  CConfig(CConfig * config, char case_filename[MAX_STRING_SIZE], SU2_COMPONENT val_software, unsigned short val_iZone, unsigned short val_nZone, bool verb_high);

  /*!
   * \brief Constructor of the class which reads the input file.
   */
  CConfig(char case_filename[MAX_STRING_SIZE], SU2_COMPONENT val_software);

  /*!
   * \brief Constructor of the class which reads the input file.
   */
  CConfig(char case_filename[MAX_STRING_SIZE], CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CConfig(void);

  /*!
  * \brief Initialize common fields of the config structure.
  */
  void Init();

  /*!
  * \brief Set the number of zones
  */
  void SetnZone();

  /*!
  * \brief Set the physical dimension of the problem
  */
  void SetnDim();

  /*!
  * \brief Print the header to screen
  * \param val_software - Kind of software component
  */
  void SetHeader(SU2_COMPONENT val_software) const;

  /*!
   * \brief Get the MPI communicator of SU2.
   * \return MPI communicator of SU2.
   */
  SU2_MPI::Comm GetMPICommunicator() const;

  /*!
   * \brief Set the MPI communicator for SU2.
   * \param[in] Communicator - MPI communicator for SU2.
   */
  void SetMPICommunicator(SU2_MPI::Comm Communicator);

  /*!
   * \brief Gets the number of zones in the mesh file.
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_format - Format of the file with the grid information.
   * \return Total number of zones in the grid file.
   */
  static unsigned short GetnZone(const string& val_mesh_filename, unsigned short val_format);

  /*!
   * \brief Gets the number of dimensions in the mesh file
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_format - Format of the file with the grid information.
   * \return Total number of domains in the grid file.
   */
  static unsigned short GetnDim(const string& val_mesh_filename, unsigned short val_format);

  /*!
   * \brief Initializes pointers to null
   */
  void SetPointersNull(void);

  /*!
   * \brief breaks an input line from the config file into a set of tokens
   * \param[in] str - the input line string
   * \param[out] option_name - the name of the option found at the beginning of the line
   * \param[out] option_value - the tokens found after the "=" sign on the line
   * \return false if the line is empty or a commment, true otherwise
   */
  bool TokenizeString(string & str, string & option_name, vector<string> & option_value);

  /*!
   * \brief Get reference origin for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin (in cartesians coordinates) for moment computation.
   */
  std::array<su2double,3> GetRefOriginMoment(unsigned short val_marker) const {
    std::array<su2double,3> RefOriginMoment{{0.0}};
    if(val_marker < nMarker_Monitoring) {
      RefOriginMoment[0] = RefOriginMoment_X[val_marker];
      RefOriginMoment[1] = RefOriginMoment_Y[val_marker];
      RefOriginMoment[2] = RefOriginMoment_Z[val_marker];
    }
    return RefOriginMoment;
  }

  /*!
   * \brief Get reference origin x-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin x-coordinate (in cartesians coordinates) for moment computation.
   */
  su2double GetRefOriginMoment_X(unsigned short val_marker) const { return RefOriginMoment_X[val_marker]; }

  /*!
   * \brief Get reference origin y-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin y-coordinate (in cartesians coordinates) for moment computation.
   */
  su2double GetRefOriginMoment_Y(unsigned short val_marker) const { return RefOriginMoment_Y[val_marker]; }

  /*!
   * \brief Get reference origin z-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin z-coordinate (in cartesians coordinates) for moment computation.
   */
  su2double GetRefOriginMoment_Z(unsigned short val_marker) const { return RefOriginMoment_Z[val_marker]; }

  /*!
   * \brief Set reference origin x-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val_origin - New x-coordinate of the mesh motion origin.
   */
  void SetRefOriginMoment_X(unsigned short val_marker, su2double val_origin) { RefOriginMoment_X[val_marker] = val_origin; }

  /*!
   * \brief Set reference origin y-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val_origin - New y-coordinate of the mesh motion origin.
   */
  void SetRefOriginMoment_Y(unsigned short val_marker, su2double val_origin) { RefOriginMoment_Y[val_marker] = val_origin; }

  /*!
   * \brief Set reference origin z-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val_origin - New z-coordinate of the mesh motion origin.
   */
  void SetRefOriginMoment_Z(unsigned short val_marker, su2double val_origin) { RefOriginMoment_Z[val_marker] = val_origin; }

  /*!
   * \brief Get index of the upper and lower horizontal plane.
   * \param[in] index - 0 means upper surface, and 1 means lower surface.
   * \return Index of the upper and lower surface.
   */
  string GetPlaneTag(unsigned short index) const { return PlaneTag[index]; }

  /*!
   * \brief Get the integration limits for the equivalent area computation.
   * \param[in] index - 0 means x_min, and 1 means x_max.
   * \return Integration limits for the equivalent area computation.
   */
  su2double GetEA_IntLimit(unsigned short index) const { return ea_lim[index]; }

  /*!
   * \brief Get the integration limits for the equivalent area computation.
   * \param[in] index - 0 means x_min, and 1 means x_max.
   * \return Integration limits for the equivalent area computation.
   */
  su2double GetEA_ScaleFactor(void) const { return EA_ScaleFactor; }

  /*!
   * \brief Get the limit value for the adjoint variables.
   * \return Limit value for the adjoint variables.
   */
  su2double GetAdjointLimit(void) const { return AdjointLimit; }

  /*!
   * \brief Get the coordinates where of the box where the grid is going to be deformed.
   * \return Coordinates where of the box where the grid is going to be deformed.
   */
  const su2double *GetHold_GridFixed_Coord(void) const { return grid_fix; }

  /*!
   * \brief Get the values of subsonic engine.
   * \return Values of subsonic engine.
   */
  const su2double *GetSubsonicEngine_Values(void) const { return eng_val; }

  /*!
   * \brief Get the cycle of a subsonic engine.
   * \return Cyl of a subsonic engine.
   */
  const su2double *GetSubsonicEngine_Cyl(void) const { return eng_cyl; }

  /*!
   * \brief Get the distortion rack.
   * \return Distortion rack.
   */
  const su2double *GetDistortionRack(void) const { return distortion; }

  /*!
   * \brief Get Description of the geometry to be analyzed
   */
  unsigned short GetGeo_Description(void) const { return Geo_Description; }

  /*!
   * \brief Creates a tecplot file to visualize the partition made by the DDC software.
   * \return <code>TRUE</code> if the partition is going to be plotted; otherwise <code>FALSE</code>.
   */
  bool GetExtraOutput(void) const { return ExtraOutput; }

  /*!
   * \brief Heat solver zone with extra screen output.
   * \return Heat solver zone with extra screen output.
   */
  long GetExtraHeatOutputZone(void) const { return ExtraHeatOutputZone; }

  /*!
   * \brief Get the value of the Mach number (velocity divided by speed of sound).
   * \return Value of the Mach number.
   */
  su2double GetMach(void) const { return Mach; }

  /*!
   * \brief Get the value of the Gamma of fluid (ratio of specific heats).
   * \return Value of the constant: Gamma
   */
  su2double GetGamma(void) const { return Gamma; }

  /*!
   * \brief Get the value of the Confinement Parameter.
   * \return Value of the constant: Confinement Parameter
   */
  su2double GetConfinement_Param(void) const { return Confinement_Param; }

  /*!
   * \brief Get the values of the CFL adaption parameters.
   * \return Value of CFL adaption parameter
   */
  su2double GetCFL_AdaptParam(unsigned short val_index) const { return CFL_AdaptParam[val_index]; }

  /*!
   * \brief Get the value of the CFL adaption flag.
   * \return <code>TRUE</code> if CFL adaption is active; otherwise <code>FALSE</code>.
   */
  bool GetCFL_Adapt(void) const { return CFL_Adapt; }

  /*!
   * \brief Get the value of the limits for the sections.
   * \return Value of the limits for the sections.
   */
  su2double GetStations_Bounds(unsigned short val_var) const { return geo_loc[val_var]; }

  /*!
   * \brief Get the value of the vector that connects the cartesian axis with a sherical or cylindrical one.
   * \return Coordinate of the Axis.
   */
  su2double GetFFD_Axis(unsigned short val_var) const { return ffd_axis[val_var]; }

  /*!
   * \brief Get the value of the bulk modulus.
   * \return Value of the bulk modulus.
   */
  su2double GetBulk_Modulus(void) const { return Bulk_Modulus; }

  /*!
   * \brief Get the epsilon^2 multiplier for Beta in the incompressible preconditioner.
   * \return Value of the epsilon^2 multiplier for Beta in the incompressible preconditioner.
   */
  su2double GetBeta_Factor(void) const { return Beta_Factor; }

  /*!
   * \brief Get the value of specific gas constant.
   * \return Value of the constant: Gamma
   */
  su2double GetGas_Constant(void) const { return Gas_Constant; }

  /*!
   * \brief Get the value of specific gas constant.
   * \return Value of the constant: Gamma
   */
  su2double GetGas_ConstantND(void) const { return Gas_ConstantND; }

  /*!
   * \brief Get the value of the molecular weight for an incompressible ideal gas (g/mol).
   * \return Value of the molecular weight for an incompressible ideal gas (g/mol).
   */
  su2double GetMolecular_Weight(unsigned short val_index = 0) const { return Molecular_Weight[val_index]; }

  /*!
   * \brief Get the value of specific heat at constant pressure.
   * \return Value of the constant: Cp
   */
  su2double GetSpecific_Heat_Cp(unsigned short val_index = 0) const { return Specific_Heat_Cp[val_index]; }

  /*!
   * \brief Get the non-dimensional value of specific heat at constant pressure.
   * \return Value of the non-dim. constant: Cp
   */
  su2double GetSpecific_Heat_CpND(unsigned short val_index = 0) const { return Specific_Heat_Cp[val_index] / Gas_Constant_Ref; }

  /*!
   * \brief Get the value of wall temperature.
   * \return Value of the constant: Temperature
   */
  su2double GetWallTemperature(void) const { return Wall_Temperature; }

    /*!
   * \brief Get the p-norm for heat-flux objective functions (adjoint problem).
   * \return Value of the heat flux p-norm
   */
  su2double GetPnormHeat(void) const { return pnorm_heat; }

  /*!
   * \brief Get the reference value for the specific gas constant.
   * \return Reference value for the specific gas constant.
   */
  su2double GetGas_Constant_Ref(void) const { return Gas_Constant_Ref; }

  /*!
   * \brief Get the reference value for the heat flux.
   * \return Reference value for the heat flux.
   */
  su2double GetHeat_Flux_Ref(void) const { return Heat_Flux_Ref; }

  /*!
   * \brief Get the value of the freestream temperature.
   * \return Freestream temperature.
   */
  su2double GetTemperature_FreeStream(void) const { return Temperature_FreeStream; }
  /*!
   * \brief Get the value of the freestream vibrational-electronic temperature.
   * \return Freestream vibe-el temperature.
   */
  su2double GetTemperature_ve_FreeStream(void) const { return Temperature_ve_FreeStream; }

  /*!
   * \brief Get the value of the freestream energy.
   * \return Freestream energy.
   */
  su2double GetEnergy_FreeStream(void) const { return Energy_FreeStream; }

  /*!
   * \brief Get the value of the freestream viscosity.
   * \return Freestream viscosity.
   */
  su2double GetViscosity_FreeStream(void) const { return Viscosity_FreeStream; }

  /*!
   * \brief Get the value of the freestream density.
   * \return Freestream density.
   */
  su2double GetDensity_FreeStream(void) const { return Density_FreeStream; }

  /*!
   * \brief Get the magnitude of the free-stream velocity of the fluid.
   * \return Magnitude of the free-stream velocity.
   */
  su2double GetModVel_FreeStream(void) const { return ModVel_FreeStream; }

  /*!
   * \brief Get the non-dimensional magnitude of the free-stream velocity of the fluid.
   * \return Non-dimensional magnitude of the free-stream velocity.
   */
  su2double GetModVel_FreeStreamND(void) const { return ModVel_FreeStreamND; }

  /*!
   * \brief Get the value of the laminar Prandtl number.
   * \return Laminar Prandtl number.
   */
  su2double GetPrandtl_Lam(unsigned short val_index = 0) const { return Prandtl_Lam[val_index]; }

  /*!
   * \brief Get the value of the turbulent Prandtl number.
   * \return Turbulent Prandtl number.
   */
  su2double GetPrandtl_Turb(unsigned short val_index = 0) const { return Prandtl_Turb[val_index]; }

  /*!
   * \brief Get the value of the von Karman constant kappa for turbulence wall modeling.
   * \return von Karman constant.
   */
  su2double GetwallModel_Kappa() const { return wallModel_Kappa; }

  /*!
   * \brief Get the value of the max. number of Newton iterations for turbulence wall modeling.
   * \return Max number of iterations.
   */
  unsigned short GetwallModel_MaxIter() const { return wallModel_MaxIter; }

  /*!
   * \brief Get the value of the relaxation factor for turbulence wall modeling.
   * \return Relaxation factor.
   */
  su2double GetwallModel_RelFac() const { return wallModel_RelFac; }

  /*!
   * \brief Get the value of the minimum Y+ value below which the wall function is deactivated.
   * \return Minimum Y+ value.
   */
  su2double GetwallModel_MinYPlus() const { return wallModel_MinYplus; }

  /*!
   * \brief Get the value of the wall model constant B for turbulence wall modeling.
   * \return Wall model constant B.
   */
  su2double GetwallModel_B() const { return wallModel_B; }

  /*!
   * \brief Get the value of the thermal diffusivity for solids.
   * \return Thermal conductivity (solid).
   */
  su2double GetThermalDiffusivity(void) const { return Thermal_Diffusivity; }

  /*!
   * \brief Get the value of the reference length for non-dimensionalization.
   *        This value should always be 1 internally, and is not user-specified.
   * \return Reference length for non-dimensionalization.
   */
  su2double GetLength_Ref(void) const { return Length_Ref; }

  /*!
   * \brief Get the value of the reference pressure for non-dimensionalization.
   * \return Reference pressure for non-dimensionalization.
   */
  su2double GetPressure_Ref(void) const { return Pressure_Ref; }

  /*!
   * \brief Get the value of the reference energy for non-dimensionalization.
   * \return Reference energy for non-dimensionalization.
   */
  su2double GetEnergy_Ref(void) const { return Energy_Ref; }

  /*!
   * \brief Get the value of the reference temperature for non-dimensionalization.
   * \return Reference temperature for non-dimensionalization.
   */
  su2double GetTemperature_Ref(void) const { return Temperature_Ref; }

  /*!
   * \brief Get the value of the reference temperature for non-dimensionalization.
   * \return Reference temperature for non-dimensionalization.
   */
  su2double GetTemperature_ve_Ref(void) const { return Temperature_ve_Ref; }

  /*!
   * \brief Get the value of the reference density for non-dimensionalization.
   * \return Reference density for non-dimensionalization.
   */
  su2double GetDensity_Ref(void) const { return Density_Ref; }

  /*!
   * \brief Get the value of the reference velocity for non-dimensionalization.
   * \return Reference velocity for non-dimensionalization.
   */
  su2double GetVelocity_Ref(void) const { return Velocity_Ref; }

  /*!
   * \brief Get the value of the reference time for non-dimensionalization.
   * \return Reference time for non-dimensionalization.
   */
  su2double GetTime_Ref(void) const { return Time_Ref; }

  /*!
   * \brief Get the value of the reference viscosity for non-dimensionalization.
   * \return Reference viscosity for non-dimensionalization.
   */
  su2double GetViscosity_Ref(void) const { return Viscosity_Ref; }

  /*!
   * \brief Get the value of the reference viscosity for non-dimensionalization.
   * \return Reference viscosity for non-dimensionalization.
   */
  su2double GetHighlite_Area(void) const { return Highlite_Area; }

  /*!
   * \brief Get the value of the reference viscosity for non-dimensionalization.
   * \return Reference viscosity for non-dimensionalization.
   */
  su2double GetFan_Poly_Eff(void) const { return Fan_Poly_Eff; }

  /*!
   * \brief Get the value of the reference thermal conductivity for non-dimensionalization.
   * \return Reference thermal conductivity for non-dimensionalization.
   */
  su2double GetThermal_Conductivity_Ref(void) const { return Thermal_Conductivity_Ref; }

  /*!
   * \brief Get the value of the reference angular velocity for non-dimensionalization.
   * \return Reference angular velocity for non-dimensionalization.
   */
  su2double GetOmega_Ref(void) const { return Omega_Ref; }

  /*!
   * \brief Get the value of the reference force for non-dimensionalization.
   * \return Reference force for non-dimensionalization.
   */
  su2double GetForce_Ref(void) const { return Force_Ref; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream pressure.
   * \return Non-dimensionalized freestream pressure.
   */
  su2double GetPressure_FreeStream(void) const { return Pressure_FreeStream; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream pressure.
   * \return Non-dimensionalized freestream pressure.
   */
  su2double GetPressure_FreeStreamND(void) const { return Pressure_FreeStreamND; }

  /*!
   * \brief Get the value of the thermodynamic pressure.
   * \return Thermodynamic pressure.
   */
  su2double GetPressure_Thermodynamic(void) const { return Pressure_Thermodynamic; }

  /*!
   * \brief Get the value of the non-dimensionalized thermodynamic pressure.
   * \return Non-dimensionalized thermodynamic pressure.
   */
  su2double GetPressure_ThermodynamicND(void) const { return Pressure_ThermodynamicND; }

  /*!
   * \brief Get the vector of the dimensionalized freestream velocity.
   * \return Dimensionalized freestream velocity vector.
   */
  su2double* GetVelocity_FreeStream(void) { return vel_inf; }
  const su2double* GetVelocity_FreeStream(void) const { return vel_inf; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream temperature.
   * \return Non-dimensionalized freestream temperature.
   */
  su2double GetTemperature_FreeStreamND(void) const { return Temperature_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized vibrational-electronic freestream temperature.
   * \return Non-dimensionalized vibrational-electronic freestream temperature.
   */
  su2double GetTemperature_ve_FreeStreamND(void) const { return Temperature_ve_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream density.
   * \return Non-dimensionalized freestream density.
   */
  su2double GetDensity_FreeStreamND(void) const { return Density_FreeStreamND; }

  /*!
   * \brief Get the vector of the non-dimensionalized freestream velocity.
   * \return Non-dimensionalized freestream velocity vector.
   */
  su2double* GetVelocity_FreeStreamND(void) { return Velocity_FreeStreamND; }
  const su2double* GetVelocity_FreeStreamND(void) const { return Velocity_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream energy.
   * \return Non-dimensionalized freestream energy.
   */
  su2double GetEnergy_FreeStreamND(void) const { return Energy_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetViscosity_FreeStreamND(void) const { return Viscosity_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetTke_FreeStreamND(void) const { return Tke_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetOmega_FreeStreamND(void) const { return Omega_FreeStreamND; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetTke_FreeStream(void) const { return Tke_FreeStream; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetOmega_FreeStream(void) const { return Omega_FreeStream; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream intermittency.
   * \return Non-dimensionalized freestream intermittency.
   */
  su2double GetIntermittency_FreeStream(void) const { return Intermittency_FreeStream; }

  /*!
   * \brief Get the value of the freestream momentum thickness Reynolds number.
   * \return Freestream momentum thickness Reynolds number.
   */
  su2double GetReThetaT_FreeStream() const { return ReThetaT_FreeStream; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream turbulence intensity.
   * \return Non-dimensionalized freestream intensity.
   */
  su2double GetTurbulenceIntensity_FreeStream(void) const { return TurbIntensityAndViscRatioFreeStream[0]; }

  /*!
   * \brief Get the value of the non-dimensionalized freestream turbulence intensity.
   * \return Non-dimensionalized freestream intensity.
   */
  su2double GetNuFactor_FreeStream(void) const { return NuFactor_FreeStream; }

  /*!
   * \brief Get the value of the non-dimensionalized engine turbulence intensity.
   * \return Non-dimensionalized engine intensity.
   */
  su2double GetNuFactor_Engine(void) const { return NuFactor_Engine; }

  /*!
   * \brief Get the value of the non-dimensionalized actuator disk turbulence intensity.
   * \return Non-dimensionalized actuator disk intensity.
   */
  su2double GetSecondaryFlow_ActDisk(void) const { return SecondaryFlow_ActDisk; }

  /*!
   * \brief Get the value of the non-dimensionalized actuator disk turbulence intensity.
   * \return Non-dimensionalized actuator disk intensity.
   */
  su2double GetInitial_BCThrust(void) const { return Initial_BCThrust; }

  /*!
   * \brief Get the value of the non-dimensionalized actuator disk turbulence intensity.
   * \return Non-dimensionalized actuator disk intensity.
   */
  void SetInitial_BCThrust(su2double val_bcthrust) { Initial_BCThrust = val_bcthrust; }

  /*!
   * \brief Get the value of the turbulent to laminar viscosity ratio.
   * \return Ratio of turbulent to laminar viscosity ratio.
   */
  su2double GetTurb2LamViscRatio_FreeStream(void) const { return TurbIntensityAndViscRatioFreeStream[1]; }

  /*!
   * \brief Get the value of the Reynolds length.
   * \return Reynolds length.
   */
  su2double GetLength_Reynolds(void) const { return Length_Reynolds; }

  /*!
   * \brief Get the reference area for non dimensional coefficient computation. If the value from the
   *        is 0 then, the code will compute the reference area using the projection of the shape into
   *        the z plane (3D) or the x plane (2D).
   * \return Value of the reference area for coefficient computation.
   */
  su2double GetRefArea(void) const { return RefArea; }

  /*!
   * \brief Get the thermal expansion coefficient.
   * \return Value of the thermal expansion coefficient.
   */
  su2double GetThermal_Expansion_Coeff(void) const { return Thermal_Expansion_Coeff; }

  /*!
   * \brief Get the non-dim. thermal expansion coefficient.
   * \return Value of the non-dim. thermal expansion coefficient.
   */
  su2double GetThermal_Expansion_CoeffND(void) const { return Thermal_Expansion_CoeffND; }

  /*!
   * \brief Set the thermal expansion coefficient.
   * \param[in] val_thermal_expansion - thermal expansion coefficient
   */
  void SetThermal_Expansion_Coeff(su2double val_thermal_expansion) { Thermal_Expansion_Coeff = val_thermal_expansion; }

  /*!
   * \brief Set the non-dim. thermal expansion coefficient.
   * \param[in] val_thermal_expansion - non-dim. thermal expansion coefficient
   */
  void SetThermal_Expansion_CoeffND(su2double val_thermal_expansionnd) { Thermal_Expansion_CoeffND = val_thermal_expansionnd; }

  /*!
   * \brief Get the value of the reference density for custom incompressible non-dimensionalization.
   * \return Reference density for custom incompressible non-dimensionalization.
   */
  su2double GetInc_Density_Ref(void) const { return Inc_Density_Ref; }

  /*!
   * \brief Get the value of the reference velocity for custom incompressible non-dimensionalization.
   * \return Reference velocity for custom incompressible non-dimensionalization.
   */
  su2double GetInc_Velocity_Ref(void) const { return Inc_Velocity_Ref; }

  /*!
   * \brief Get the value of the reference temperature for custom incompressible non-dimensionalization.
   * \return Reference temperature for custom incompressible non-dimensionalization.
   */
  su2double GetInc_Temperature_Ref(void) const { return Inc_Temperature_Ref; }

  /*!
   * \brief Get the value of the initial density for incompressible flows.
   * \return Initial density for incompressible flows.
   */
  su2double GetInc_Density_Init(void) const { return Inc_Density_Init; }

  /*!
   * \brief Get the value of the initial velocity for incompressible flows.
   * \return Initial velocity for incompressible flows.
   */
  const su2double* GetInc_Velocity_Init(void) const { return vel_init; }

  /*!
   * \brief Get the value of the initial temperature for incompressible flows.
   * \return Initial temperature for incompressible flows.
   */
  su2double GetInc_Temperature_Init(void) const { return Inc_Temperature_Init; }

  /*!
   * \brief Get the flag for activating species transport clipping.
   * \return Flag for species clipping.
   */
  bool GetSpecies_Clipping() const { return Species_Clipping; }

  /*!
   * \brief Get the maximum bound for scalar transport clipping
   * \return Maximum value for scalar clipping
   */
  su2double GetSpecies_Clipping_Max(unsigned short iVar) const { return Species_Clipping_Max[iVar]; }

  /*!
   * \brief Get the minimum bound for scalar transport clipping
   * \return Minimum value for scalar clipping
   */
  su2double GetSpecies_Clipping_Min(unsigned short iVar) const { return Species_Clipping_Min[iVar]; }

  /*!
   * \brief Get initial species value/concentration in the range [0,1].
   * \return Initial species value/concentration
   */
  const su2double* GetSpecies_Init() const { return Species_Init; }

  /*!
   * \brief Get the flag for using strong BC's for in- and outlets in the species solver.
   * \return Flag for strong BC's.
   */
  bool GetSpecies_StrongBC() const { return Species_StrongBC; }

  /*!
   * \brief Get the flame initialization.
   *        (x1,x2,x3) = flame offset.
   *        (x4,x5,x6) = flame normal, separating unburnt from burnt.
   *        (x7) = flame thickness, the length from unburnt to burnt conditions.
   *        (x8) = flame burnt thickness, the length to stay at burnt conditions.
   * \return Flame initialization for the flamelet model.
   */
  const su2double* GetFlameInit() const { return flame_init; }

  /*!
   * \brief Get the number of control variables for flamelet model.
   */
  unsigned short GetNControlVars() const { return n_control_vars; }

  /*!
   * \brief Get the number of total transported scalars for flamelet model.
   */
  unsigned short GetNScalars() const { return n_scalars; }

  /*!
   * \brief Get the number of user scalars for flamelet model.
   */
  unsigned short GetNUserScalars() const { return n_user_scalars; }

  /*!
   * \brief Get the name of a specific controlling variable.
   */
  const string& GetControllingVariableName(unsigned short i_cv) const {
    return controlling_variable_names[i_cv];
  }

  /*!
   * \brief Get the name of the source term variable for a specific controlling variable.
   */
  const string& GetControllingVariableSourceName(unsigned short i_cv) const {
    return cv_source_names[i_cv];
  }
  /*!
   * \brief Get the name of the user scalar.
   */
  const string& GetUserScalarName(unsigned short i_user_scalar) const {
    static const std::string none = "NONE";
    if (n_user_scalars > 0) return user_scalar_names[i_user_scalar]; else return none;
  }

  /*!
   * \brief Get the name of the user scalar source term.
   */
  const string& GetUserSourceName(unsigned short i_user_source) const {
    static const std::string none = "NONE";
    if (n_user_sources > 0) return user_source_names[i_user_source]; else return none;
  }

  /*!
   * \brief Get the number of transported scalars for combustion.
   */
  unsigned short GetNLookups() const { return n_lookups; }

  /*!
   * \brief Get the name of the variable that we want to retrieve from the lookup table.
   */
  const string& GetLookupName(unsigned short i_lookup) const { return lookup_names[i_lookup]; }

  /*!
   * \brief Get the Young's modulus of elasticity.
   * \return Value of the Young's modulus of elasticity.
   */
  su2double GetElasticyMod(unsigned short id_val) const { return ElasticityMod[id_val]; }

  /*!
   * \brief Decide whether to apply DE effects to the model.
   * \return <code>TRUE</code> if the DE effects are to be applied, <code>FALSE</code> otherwise.
   */
  bool GetDE_Effects(void) const { return DE_Effects; }

  /*!
   * \brief Get the number of different electric constants.
   * \return Value of the DE modulus.
   */
  unsigned short GetnElectric_Constant(void) const { return nElectric_Constant; }

  /*!
   * \brief Get the value of the DE modulus.
   * \return Value of the DE modulus.
   */
  su2double GetElectric_Constant(unsigned short iVar) const { return Electric_Constant[iVar]; }

  /*!
   * \brief Get the value of the B constant in the Knowles material model.
   * \return Value of the B constant in the Knowles material model.
   */
  su2double GetKnowles_B(void) const { return Knowles_B; }

  /*!
   * \brief Get the value of the N constant in the Knowles material model.
   * \return Value of the N constant in the Knowles material model.
   */
  su2double GetKnowles_N(void) const { return Knowles_N; }

  /*!
   * \brief Get the kind of design variable for FEA.
   * \return Value of the DE voltage.
   */
  unsigned short GetDV_FEA(void) const { return Kind_DV_FEA; }

  /*!
   * \brief Get the ID of the reference node.
   * \return Number of FSI subiters.
   */
  unsigned long GetRefNode_ID(void) const { return refNodeID; }

  /*!
   * \brief Get the values for the reference node displacement.
   * \param[in] val_coeff - Index of the displacement.
   */
  su2double GetRefNode_Displacement(unsigned short val_coeff) const { return RefNode_Displacement[val_coeff]; }

  /*!
   * \brief Get the penalty weight value for the objective function.
   * \return  Penalty weight value for the reference geometry objective function.
   */
  su2double GetRefNode_Penalty(void) const { return RefNode_Penalty; }

  /*!
   * \brief Decide whether it's necessary to read a reference geometry.
   */
  bool GetRefGeom(void) const { return RefGeom; }

  /*!
   * \brief Consider only the surface of the reference geometry.
   */
  bool GetRefGeomSurf(void) const { return RefGeomSurf; }

  /*!
   * \brief Get the name of the file with the reference geometry of the structural problem.
   * \return Name of the file with the reference geometry of the structural problem.
   */
  string GetRefGeom_FEMFileName(void) const { return RefGeom_FEMFileName; }

  /*!
   * \brief Get the format of the reference geometry file.
   * \return Format of the reference geometry file.
   */
  unsigned short GetRefGeom_FileFormat(void) const { return RefGeom_FileFormat; }

  /*!
   * \brief Formulation for 2D elasticity (plane stress - strain)
   * \return Flag to 2D elasticity model.
   */
  STRUCT_2DFORM GetElas2D_Formulation() const { return Kind_2DElasForm; }

  /*!
   * \brief Decide whether it's necessary to read a reference geometry.
   * \return <code>TRUE</code> if it's necessary to read a reference geometry, <code>FALSE</code> otherwise.
   */
  bool GetPrestretch(void) const { return Prestretch; }

  /*!
   * \brief Get the name of the file with the element properties for structural problems.
   * \return Name of the file with the element properties of the structural problem.
   */
  string GetFEA_FileName(void) const { return FEA_FileName; }

  /*!
   * \brief Determine if advanced features are used from the element-based FEA analysis (experimental feature).
   * \return <code>TRUE</code> is experimental, <code>FALSE</code> is the default behaviour.
   */
  inline bool GetAdvanced_FEAElementBased(void) const { return FEAAdvancedMode; }

  /*!
   * \brief Get the name of the file with the reference geometry of the structural problem.
   * \return Name of the file with the reference geometry of the structural problem.
   */
  string GetPrestretch_FEMFileName(void) const { return Prestretch_FEMFileName; }

  /*!
   * \brief Get the Poisson's ratio.
   * \return Value of the Poisson's ratio.
   */
  su2double GetPoissonRatio(unsigned short id_val) const { return PoissonRatio[id_val]; }

  /*!
   * \brief Get the Material Density.
   * \return Value of the Material Density.
   */
  su2double GetMaterialDensity(unsigned short id_val) const { return MaterialDensity[id_val]; }

  /*!
   * \brief Compressibility/incompressibility of the solids analysed using the structural solver.
   * \return Compressible or incompressible.
   */
  STRUCT_COMPRESS GetMaterialCompressibility(void) const { return Kind_Material_Compress; }

  /*!
   * \brief Compressibility/incompressibility of the solids analysed using the structural solver.
   * \return Compressible or incompressible.
   */
  STRUCT_MODEL GetMaterialModel(void) const { return Kind_Material; }

  /*!
   * \brief Geometric conditions for the structural solver.
   * \return Small or large deformation structural analysis.
   */
  STRUCT_DEFORMATION GetGeometricConditions(void) const { return Kind_Struct_Solver; }

  /*!
   * \brief Get the reference length for computing moment (the default value is 1).
   * \return Reference length for moment computation.
   */
  su2double GetRefLength(void) const { return RefLength; }

  /*!
   * \brief Get the reference element length for computing the slope limiting epsilon.
   * \return Reference element length for slope limiting epsilon.
   */
  su2double GetRefElemLength(void) const { return RefElemLength; }

  /*!
   * \brief Get the reference coefficient for detecting sharp edges.
   * \return Reference coefficient for detecting sharp edges.
   */
  su2double GetRefSharpEdges(void) const { return RefSharpEdges; }

  /*!
   * \brief Get the volume of the whole domain using the fine grid, this value is common for all the grids
   *        in the multigrid method.
   * \return Volume of the whole domain.
   */
  su2double GetDomainVolume(void) const { return DomainVolume; }

  /*!
   * \brief In case the <i>RefArea</i> is equal to 0 then, it is necessary to compute a reference area,
   *        with this function we set the value of the reference area.
   * \param[in] val_area - Value of the reference area for non dimensional coefficient computation.
   */
  void SetRefArea(su2double val_area) { RefArea = val_area; }

  /*!
   * \brief In case the <i>SemiSpan</i> is equal to 0 then, it is necessary to compute the max y distance,
   *        with this function we set the value of the semi span.
   * \param[in] val_semispan - Value of the semispan.
   */
  void SetSemiSpan(su2double val_semispan) { SemiSpan = val_semispan; }

  /*!
   * \brief Set the value of the domain volume computed on the finest grid.
   * \note This volume do not include the volume of the body that is being simulated.
   * \param[in] val_volume - Value of the domain volume computed on the finest grid.
   */
  void SetDomainVolume(su2double val_volume) { DomainVolume = val_volume; }

  /*!
   * \brief Set the finest mesh in a multigrid strategy.
   * \note If we are using a Full Multigrid Strategy or a start up with finest grid, it is necessary
   *       to change several times the finest grid.
   * \param[in] val_finestmesh - Index of the finest grid.
   */
  void SetFinestMesh(unsigned short val_finestmesh) { FinestMesh = val_finestmesh; }

  /*!
   * \brief Set the kind of time integration scheme.
   * \note If we are solving different equations it will be necessary to change several
   *       times the kind of time integration, to choose the right scheme.
   * \param[in] val_kind_timeintscheme - Kind of time integration scheme.
   */
  void SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme) { Kind_TimeNumScheme = val_kind_timeintscheme; }

  /*!
   * \brief Set the parameters of the convective numerical scheme.
   * \note The parameters will change because we are solving different kind of equations.
   * \param[in] val_kind_convnumscheme - Center or upwind scheme.
   * \param[in] val_kind_centered - If centered scheme, kind of centered scheme (JST, etc.).
   * \param[in] val_kind_upwind - If upwind scheme, kind of upwind scheme (Roe, etc.).
   * \param[in] val_kind_slopelimit - If upwind scheme, kind of slope limit.
   * \param[in] val_muscl - Define if we apply a MUSCL scheme or not.
   * \param[in] val_kind_fem - If FEM, what kind of FEM discretization.
   */
  void SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme, CENTERED val_kind_centered,
                             UPWIND val_kind_upwind, LIMITER val_kind_slopelimit,
                             bool val_muscl,  unsigned short val_kind_fem);

  /*!
   * \brief Get the value of limiter coefficient.
   * \return Value of the limiter coefficient.
   */
  su2double GetVenkat_LimiterCoeff(void) const { return Venkat_LimiterCoeff; }

  /*!
   * \brief Freeze the value of the limiter after a number of iterations.
   * \return Number of iterations.
   */
  unsigned long GetLimiterIter(void) const { return LimiterIter; }

  /*!
   * \brief Get the value of sharp edge limiter.
   * \return Value of the sharp edge limiter coefficient.
   */
  su2double GetAdjSharp_LimiterCoeff(void) const { return AdjSharp_LimiterCoeff; }

  /*!
   * \brief Get the Reynolds number. Dimensionless number that gives a measure of the ratio of inertial forces
   *        to viscous forces and consequently quantifies the relative importance of these two types of forces
   *        for given flow condition.
   * \return Value of the Reynolds number.
   */
  su2double GetReynolds(void) const { return Reynolds; }

  /*!
   * \brief Get the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  su2double GetFroude(void) const { return Froude; }

  /*!
   * \brief Set the Froude number for free surface problems.
   * \param[in] val_froude - Value of the Froude number.
   */
  void SetFroude(su2double val_froude) { Froude = val_froude; }

  /*!
   * \brief Set the Mach number.
   * \param[in] val_mach - Value of the Mach number.
   */
  void SetMach(su2double val_mach) { Mach = val_mach; }

  /*!
   * \brief Set the Reynolds number.
   * \param[in] val_reynolds - Value of the Reynolds number.
   */
  void SetReynolds(su2double val_reynolds) { Reynolds = val_reynolds; }

  /*!
   * \brief Set the reference length for nondimensionalization.
   * \param[in] val_length_ref - Value of the reference length.
   */
  void SetLength_Ref(su2double val_length_ref) { Length_Ref = val_length_ref; }

  /*!
   * \brief Set the reference velocity for nondimensionalization.
   * \param[in] val_velocity_ref - Value of the reference velocity.
   */
  void SetVelocity_Ref(su2double val_velocity_ref) { Velocity_Ref = val_velocity_ref; }

  /*!
   * \brief Set the reference pressure for nondimensionalization.
   * \param[in] val_pressure_ref - Value of the reference pressure.
   */
  void SetPressure_Ref(su2double val_pressure_ref) { Pressure_Ref = val_pressure_ref; }

  /*!
   * \brief Set the reference pressure for nondimensionalization.
   * \param[in] val_density_ref - Value of the reference pressure.
   */
  void SetDensity_Ref(su2double val_density_ref) { Density_Ref = val_density_ref; }

  /*!
   * \brief Set the reference temperature for nondimensionalization.
   * \param[in] val_temperature_ref - Value of the reference temperature.
   */
  void SetTemperature_Ref(su2double val_temperature_ref) { Temperature_Ref = val_temperature_ref; }

  /*!
   * \brief Set the reference temperature.
   * \param[in] val_temperature_ve_ref - Value of the reference temperature.
   */
  void SetTemperature_ve_Ref(su2double val_temperature_ve_ref) { Temperature_ve_Ref = val_temperature_ve_ref; }

  /*!
   * \brief Set the reference time for nondimensionalization.
   * \param[in] val_time_ref - Value of the reference time.
   */
  void SetTime_Ref(su2double val_time_ref) { Time_Ref = val_time_ref; }

  /*!
   * \brief Set the reference energy for nondimensionalization.
   * \param[in] val_energy_ref - Value of the reference energy.
   */
  void SetEnergy_Ref(su2double val_energy_ref) { Energy_Ref = val_energy_ref; }

  /*!
   * \brief Set the reference Omega for nondimensionalization.
   * \param[in] val_omega_ref - Value of the reference omega.
   */
  void SetOmega_Ref(su2double val_omega_ref) { Omega_Ref = val_omega_ref; }

  /*!
   * \brief Set the reference Force for nondimensionalization.
   * \param[in] val_force_ref - Value of the reference Force.
   */
  void SetForce_Ref(su2double val_force_ref) { Force_Ref = val_force_ref; }

  /*!
   * \brief Set the reference gas-constant for nondimensionalization.
   * \param[in] val_gas_constant_ref - Value of the reference gas-constant.
   */
  void SetGas_Constant_Ref(su2double val_gas_constant_ref) { Gas_Constant_Ref = val_gas_constant_ref; }

  /*!
   * \brief Set the gas-constant.
   * \param[in] val_gas_constant - Value of the gas-constant.
   */
  void SetGas_Constant(su2double val_gas_constant) { Gas_Constant = val_gas_constant; }

  /*!
   * \brief Set the heat flux reference value.
   * \return Value of the reference heat flux.
   */
  void SetHeat_Flux_Ref(su2double val_heat_flux_ref) { Heat_Flux_Ref = val_heat_flux_ref; }

  /*!
   * \brief Set the reference viscosity for nondimensionalization.
   * \param[in] val_viscosity_ref - Value of the reference viscosity.
   */
  void SetViscosity_Ref(su2double val_viscosity_ref) { Viscosity_Ref = val_viscosity_ref; }

  /*!
   * \brief Set the reference conductivity for nondimensionalization.
   * \param[in] val_conductivity_ref - Value of the reference conductivity.
   */
  void SetConductivity_Ref(su2double val_conductivity_ref) { Thermal_Conductivity_Ref = val_conductivity_ref; }

  /*!
   * \brief Set the nondimensionalized freestream pressure.
   * \param[in] val_pressure_freestreamnd - Value of the nondimensionalized freestream pressure.
   */
  void SetPressure_FreeStreamND(su2double val_pressure_freestreamnd) { Pressure_FreeStreamND = val_pressure_freestreamnd; }

  /*!
   * \brief Set the freestream pressure.
   * \param[in] val_pressure_freestream - Value of the freestream pressure.
   */
  void SetPressure_FreeStream(su2double val_pressure_freestream) { Pressure_FreeStream = val_pressure_freestream; }

  /*!
   * \brief Set the non-dimensionalized thermodynamic pressure for low Mach problems.
   * \return Value of the non-dimensionalized thermodynamic pressure.
   */
  void SetPressure_ThermodynamicND(su2double val_pressure_thermodynamicnd) { Pressure_ThermodynamicND = val_pressure_thermodynamicnd; }

  /*!
   * \brief Set the thermodynamic pressure for low Mach problems.
   * \return Value of the thermodynamic pressure.
   */
  void SetPressure_Thermodynamic(su2double val_pressure_thermodynamic) { Pressure_Thermodynamic = val_pressure_thermodynamic; }

  /*!
   * \brief Set the nondimensionalized freestream density.
   * \param[in] val_density_freestreamnd - Value of the nondimensionalized freestream density.
   */
  void SetDensity_FreeStreamND(su2double val_density_freestreamnd) { Density_FreeStreamND = val_density_freestreamnd; }

  /*!
   * \brief Set the freestream density.
   * \param[in] val_density_freestream - Value of the freestream density.
   */
  void SetDensity_FreeStream(su2double val_density_freestream) { Density_FreeStream = val_density_freestream; }

  /*!
   * \brief Set the freestream viscosity.
   * \param[in] val_viscosity_freestream - Value of the freestream viscosity.
   */
  void SetViscosity_FreeStream(su2double val_viscosity_freestream) { Viscosity_FreeStream = val_viscosity_freestream; }

  /*!
   * \brief Set the magnitude of the free-stream velocity.
   * \param[in] val_modvel_freestream - Magnitude of the free-stream velocity.
   */
  void SetModVel_FreeStream(su2double val_modvel_freestream) { ModVel_FreeStream = val_modvel_freestream; }

  /*!
   * \brief Set the non-dimensional magnitude of the free-stream velocity.
   * \param[in] val_modvel_freestreamnd - Non-dimensional magnitude of the free-stream velocity.
   */
  void SetModVel_FreeStreamND(su2double val_modvel_freestreamnd) { ModVel_FreeStreamND = val_modvel_freestreamnd; }

  /*!
   * \brief Set the freestream temperature.
   * \param[in] val_temperature_freestream - Value of the freestream temperature.
   */
  void SetTemperature_FreeStream(su2double val_temperature_freestream) { Temperature_FreeStream = val_temperature_freestream; }

  /*!
   * \brief Set the non-dimensional freestream temperature.
   * \param[in] val_temperature_freestreamnd - Value of the non-dimensional freestream temperature.
   */
  void SetTemperature_FreeStreamND(su2double val_temperature_freestreamnd) { Temperature_FreeStreamND = val_temperature_freestreamnd; }

  /*!
   * \brief Set the freestream vibrational-electronic temperature.
   * \param[in] val_temperature_ve_freestream - Value of the freestream vibrational-electronic temperature.
   */
  void SetTemperature_ve_FreeStream(su2double val_temperature_ve_freestream) { Temperature_ve_FreeStream = val_temperature_ve_freestream; }

  /*!
   * \brief Set the non-dimensional freestream vibrational-electronic temperature.
   * \param[in] val_temperature_ve_freestreamnd - Value of the non-dimensional freestream vibrational-electronic temperature.
   */
  void SetTemperature_ve_FreeStreamND(su2double val_temperature_ve_freestreamnd) { Temperature_ve_FreeStreamND = val_temperature_ve_freestreamnd; }

  /*!
   * \brief Set the non-dimensional gas-constant.
   * \param[in] val_gas_constantnd - Value of the non-dimensional gas-constant.
   */
  void SetGas_ConstantND(su2double val_gas_constantnd) { Gas_ConstantND = val_gas_constantnd; }

  /*!
   * \brief Set the free-stream velocity.
   * \param[in] val_velocity_freestream - Value of the free-stream velocity component.
   * \param[in] val_dim - Value of the current dimension.
   */
  void SetVelocity_FreeStream(su2double val_velocity_freestream, unsigned short val_dim) { vel_inf[val_dim] = val_velocity_freestream; }

  /*!
   * \brief Set the non-dimensional free-stream velocity.
   * \param[in] val_velocity_freestreamnd - Value of the non-dimensional free-stream velocity component.
   * \param[in] val_dim - Value of the current dimension.
   */
  void SetVelocity_FreeStreamND(su2double val_velocity_freestreamnd, unsigned short val_dim) { Velocity_FreeStreamND[val_dim] = val_velocity_freestreamnd; }

  /*!
   * \brief Set the non-dimensional free-stream viscosity.
   * \param[in] val_viscosity_freestreamnd - Value of the non-dimensional free-stream viscosity.
   */
  void SetViscosity_FreeStreamND(su2double val_viscosity_freestreamnd) { Viscosity_FreeStreamND = val_viscosity_freestreamnd; }

  /*!
   * \brief Set the non-dimensional freestream turbulent kinetic energy.
   * \param[in] val_tke_freestreamnd - Value of the non-dimensional freestream turbulent kinetic energy.
   */
  void SetTke_FreeStreamND(su2double val_tke_freestreamnd) { Tke_FreeStreamND = val_tke_freestreamnd; }

  /*!
   * \brief Set the non-dimensional freestream specific dissipation rate omega.
   * \param[in] val_omega_freestreamnd - Value of the non-dimensional freestream specific dissipation rate omega.
   */
  void SetOmega_FreeStreamND(su2double val_omega_freestreamnd) { Omega_FreeStreamND = val_omega_freestreamnd; }

  /*!
   * \brief Set the freestream turbulent kinetic energy.
   * \param[in] val_tke_freestream - Value of the freestream turbulent kinetic energy.
   */
  void SetTke_FreeStream(su2double val_tke_freestream) { Tke_FreeStream = val_tke_freestream; }

  /*!
   * \brief Set the freestream specific dissipation rate omega.
   * \param[in] val_omega_freestream - Value of the freestream specific dissipation rate omega.
   */
  void SetOmega_FreeStream(su2double val_omega_freestream) { Omega_FreeStream = val_omega_freestream; }

  /*!
   * \brief Set the freestream momentum thickness Reynolds number.
   * \param[in] val_ReThetaT_freestream - Value of the freestream momentum thickness Reynolds number.
   */
  void SetReThetaT_FreeStream(su2double val_ReThetaT_freestream) { ReThetaT_FreeStream = val_ReThetaT_freestream; }

  /*!
   * \brief Set the non-dimensional freestream energy.
   * \param[in] val_energy_freestreamnd - Value of the non-dimensional freestream energy.
   */
  void SetEnergy_FreeStreamND(su2double val_energy_freestreamnd) { Energy_FreeStreamND = val_energy_freestreamnd; }

  /*!
   * \brief Set the freestream energy.
   * \param[in] val_energy_freestream - Value of the freestream energy.
   */
  void SetEnergy_FreeStream(su2double val_energy_freestream) { Energy_FreeStream = val_energy_freestream; }

  /*!
   * \brief Set the thermal diffusivity for solids.
   * \param[in] val_thermal_diffusivity - Value of the thermal diffusivity.
   */
  void SetThermalDiffusivity(su2double val_thermal_diffusivity) { Thermal_Diffusivity = val_thermal_diffusivity; }

  /*!
   * \brief Set the non-dimensional total time for unsteady simulations.
   * \param[in] val_total_unsttimend - Value of the non-dimensional total time.
   */
  void SetTotal_UnstTimeND(su2double val_total_unsttimend) { Total_UnstTimeND = val_total_unsttimend; }

  /*!
   * \brief Get the angle of attack of the body. This is the angle between a reference line on a lifting body
   *        (often the chord line of an airfoil) and the vector representing the relative motion between the
   *        lifting body and the fluid through which it is moving.
   * \return Value of the angle of attack.
   */
  su2double GetAoA(void) const { return AoA; }

  /*!
   * \brief Get the off set angle of attack of the body. The solution and the geometry
   *        file are able to modifity the angle of attack in the config file
   * \return Value of the off set angle of attack.
   */
  su2double GetAoA_Offset(void) const { return AoA_Offset; }

  /*!
   * \brief Get the off set sideslip angle of the body. The solution and the geometry
   *        file are able to modifity the angle of attack in the config file
   * \return Value of the off set sideslip angle.
   */
  su2double GetAoS_Offset(void) const { return AoS_Offset; }

  /*!
   * \brief Get the functional sensitivity with respect to changes in the angle of attack.
   * \return Value of the angle of attack.
   */
  su2double GetAoA_Sens(void) const { return AoA_Sens; }

  /*!
   * \brief Set the angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoA(su2double val_AoA) { AoA = val_AoA; }

  /*!
   * \brief Set the off set angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoA_Offset(su2double val_AoA_offset) { AoA_Offset = val_AoA_offset; }

  /*!
   * \brief Set the off set sideslip angle.
   * \param[in] val_AoA - Value of the off set sideslip angle.
   */
  void SetAoS_Offset(su2double val_AoS_offset) { AoS_Offset = val_AoS_offset; }

  /*!
   * \brief Set the angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoA_Sens(su2double val_AoA_sens) { AoA_Sens = val_AoA_sens; }

  /*!
   * \brief Set the angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoS(su2double val_AoS) { AoS = val_AoS; }

  /*!
   * \brief Get the angle of sideslip of the body. It relates to the rotation of the aircraft centerline from
   *        the relative wind.
   * \return Value of the angle of sideslip.
   */
  su2double GetAoS(void) const { return AoS; }

  /*!
   * \brief Get the charge coefficient that is used in the poissonal potential simulation.
   * \return Value of the charge coefficient.
   */
  su2double GetChargeCoeff(void) const { return ChargeCoeff; }

  /*!
   * \brief Get the number of multigrid levels.
   * \return Number of multigrid levels (without including the original grid).
   */
  unsigned short GetnMGLevels(void) const { return nMGLevels; }

  /*!
   * \brief Set the number of multigrid levels.
   * \param[in] val_nMGLevels - Index of the mesh were the CFL is applied
   */
  void SetMGLevels(unsigned short val_nMGLevels) {
    nMGLevels = val_nMGLevels;
    if (MGCycle == FULLMG_CYCLE) {
      SetFinestMesh(val_nMGLevels);
    }
  }

  /*!
   * \brief Get the index of the finest grid.
   * \return Index of the finest grid in a multigrid strategy, this is 0 unless we are
   performing a Full multigrid.
   */
  unsigned short GetFinestMesh(void) const { return FinestMesh; }

  /*!
   * \brief Get the kind of multigrid (V or W).
   * \note This variable is used in a recursive way to perform the different kind of cycles
   * \return 0 or 1 depending of we are dealing with a V or W cycle.
   */
  unsigned short GetMGCycle(void) const { return MGCycle; }

  /*!
   * \brief Get the king of evaluation in the geometrical module.
   * \return 0 or 1 depending of we are dealing with a V or W cycle.
   */
  unsigned short GetGeometryMode(void) const { return GeometryMode; }

  /*!
   * \brief Get the Courant Friedrich Levi number for each grid.
   * \param[in] val_mesh - Index of the mesh were the CFL is applied.
   * \return CFL number for each grid.
   */
  su2double GetCFL(unsigned short val_mesh) const { return CFL[val_mesh]; }

  /*!
   * \brief Get the Courant Friedrich Levi number for each grid.
   * \param[in] val_mesh - Index of the mesh were the CFL is applied.
   * \return CFL number for each grid.
   */
  void SetCFL(unsigned short val_mesh, su2double val_cfl) { CFL[val_mesh] = val_cfl; }

  /*!
   * \brief Get the Courant Friedrich Levi number for unsteady simulations.
   * \return CFL number for unsteady simulations.
   */
  su2double GetUnst_CFL(void) const { return Unst_CFL; }

  /*!
   * \brief Get information about element reorientation
   * \return    <code>TRUE</code> means that elements can be reoriented if suspected unhealthy
   */
  bool GetReorientElements(void) const { return ReorientElements; }

  /*!
   * \brief Get the Courant Friedrich Levi number for unsteady simulations.
   * \return CFL number for unsteady simulations.
   */
  su2double GetMax_DeltaTime(void) const { return Max_DeltaTime; }

  /*!
   * \brief Get a parameter of the particular design variable.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \param[in] val_param - Index of the parameter that we want to read.
   * \return Design variable parameter.
   */
  su2double GetParamDV(unsigned short val_dv, unsigned short val_param) const { return ParamDV[val_dv][val_param]; }

  /*!
   * \brief Get the coordinates of the FFD corner points.
   * \param[in] val_ffd - Index of the FFD box.
   * \param[in] val_coord - Index of the coordinate that we want to read.
   * \return Value of the coordinate.
   */
  su2double GetCoordFFDBox(unsigned short val_ffd, unsigned short val_index) const { return CoordFFDBox[val_ffd][val_index]; }

  /*!
   * \brief Get the degree of the FFD corner points.
   * \param[in] val_ffd - Index of the FFD box.
   * \param[in] val_degree - Index (I,J,K) to obtain the degree.
   * \return Value of the degree in a particular direction.
   */
  unsigned short GetDegreeFFDBox(unsigned short val_ffd, unsigned short val_index) const { return DegreeFFDBox[val_ffd][val_index]; }

  /*!
   * \brief Get the FFD Tag of a particular design variable.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \return Name of the FFD box.
   */
  string GetFFDTag(unsigned short val_dv) const { return FFDTag[val_dv]; }

  /*!
   * \brief Get the FFD Tag of a particular FFD box.
   * \param[in] val_ffd - Number of the FFD box that we want to read.
   * \return Name of the FFD box.
   */
  string GetTagFFDBox(unsigned short val_ffd) const { return TagFFDBox[val_ffd]; }

  /*!
   * \brief Get the number of design variables.
   * \return Number of the design variables.
   */
  unsigned short GetnDV(void) const { return nDV; }

  /*!
   * \brief Get the number of design variables.
   * \return Number of the design variables.
   */
  unsigned short GetnDV_Value(unsigned short iDV) const { return nDV_Value[iDV]; }

  /*!
   * \brief Get the total number of design variables.
   */
  unsigned short GetnDV_Total(void) const {
    if (!nDV_Value) return 0;
    unsigned short sum = 0;
    for (unsigned short iDV = 0; iDV < nDV; iDV++) {
      sum += nDV_Value[iDV];
    }
    return sum;
  }

  /*!
   * \brief Get the number of FFD boxes.
   * \return Number of FFD boxes.
   */
  unsigned short GetnFFDBox(void) const { return nFFDBox; }

  /*!
   * \brief Get the required continuity level at the surface intersection with the FFD
   * \return Continuity level at the surface intersection.
   */
  unsigned short GetFFD_Continuity(void) const { return FFD_Continuity; }

  /*!
   * \brief Get the coordinate system that we are going to use to define the FFD
   * \return Coordinate system (cartesian, spherical, etc).
   */
  unsigned short GetFFD_CoordSystem(void) const { return FFD_CoordSystem; }

  /*!
   * \brief Get the kind of FFD Blending function.
   * \return Kind of FFD Blending function.
   */
  unsigned short GetFFD_Blending(void) const { return FFD_Blending;}

  /*!
   * \brief Get the kind BSpline Order in i,j,k direction.
   * \return The kind BSpline Order in i,j,k direction.
   */
  const su2double* GetFFD_BSplineOrder() const { return ffd_coeff;}

  /*!
   * \brief Get the number of Runge-Kutta steps.
   * \return Number of Runge-Kutta steps.
   */
  unsigned short GetnRKStep(void) const { return nRKStep; }

  /*!
   * \brief Get the number of time levels for time accurate local time stepping.
   * \return Number of time levels.
   */
  unsigned short GetnLevels_TimeAccurateLTS(void) const { return nLevels_TimeAccurateLTS; }

  /*!
   * \brief Set the number of time levels for time accurate local time stepping.
   * \param[in] val_nLevels - The number of time levels to be set.
   */
  void SetnLevels_TimeAccurateLTS(unsigned short val_nLevels) { nLevels_TimeAccurateLTS = val_nLevels;}

  /*!
   * \brief Get the number time DOFs for ADER-DG.
   * \return Number of time DOFs used in ADER-DG.
   */
  unsigned short GetnTimeDOFsADER_DG(void) const { return nTimeDOFsADER_DG; }

  /*!
   * \brief Get the location of the time DOFs for ADER-DG on the interval [-1..1].
   * \return The location of the time DOFs used in ADER-DG.
   */
  const su2double *GetTimeDOFsADER_DG(void) const { return TimeDOFsADER_DG; }

  /*!
   * \brief Get the number time integration points for ADER-DG.
   * \return Number of time integration points used in ADER-DG.
   */
  unsigned short GetnTimeIntegrationADER_DG(void) const { return nTimeIntegrationADER_DG; }

  /*!
   * \brief Get the location of the time integration points for ADER-DG on the interval [-1..1].
   * \return The location of the time integration points used in ADER-DG.
   */
  const su2double *GetTimeIntegrationADER_DG(void) const { return TimeIntegrationADER_DG; }

  /*!
   * \brief Get the weights of the time integration points for ADER-DG.
   * \return The weights of the time integration points used in ADER-DG.
   */
  const su2double *GetWeightsIntegrationADER_DG(void) const { return WeightsIntegrationADER_DG; }

  /*!
   * \brief Get the total number of boundary markers of the local process including send/receive domains.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_All(void) const { return nMarker_All; }

  /*!
   * \brief Get the total number of boundary markers in the config file.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_CfgFile(void) const { return nMarker_CfgFile; }

  /*!
   * \brief Get the number of Euler boundary markers.
   * \return Number of Euler boundary markers.
   */
  unsigned short GetnMarker_Euler(void) const { return nMarker_Euler; }

  /*!
   * \brief Get the number of symmetry boundary markers.
   * \return Number of symmetry boundary markers.
   */
  unsigned short GetnMarker_SymWall(void) const { return nMarker_SymWall; }

  /*!
   * \brief Get the total number of boundary markers in the cfg plus the possible send/receive domains.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_Max(void) const { return nMarker_Max; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_EngineInflow(void) const { return nMarker_EngineInflow; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_EngineExhaust(void) const { return nMarker_EngineExhaust; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_NearFieldBound(void) const { return nMarker_NearFieldBound; }

  /*!
   * \brief Get the total number of deformable markers at the boundary.
   * \return Total number of deformable markers at the boundary.
   */
  unsigned short GetnMarker_Deform_Mesh(void) const { return nMarker_Deform_Mesh; }

  /*!
   * \brief Get the total number of markers in which the flow load is computed/employed.
   * \return Total number of markers in which the flow load is computed/employed.
   */
  unsigned short GetnMarker_Fluid_Load(void) const { return nMarker_Fluid_Load; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_Fluid_InterfaceBound(void) const { return nMarker_Fluid_InterfaceBound; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_ActDiskInlet(void) const { return nMarker_ActDiskInlet; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_ActDiskOutlet(void) const { return nMarker_ActDiskOutlet; }

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_Outlet(void) const { return nMarker_Outlet; }

  /*!
   * \brief Get the total number of monitoring markers.
   * \return Total number of monitoring markers.
   */
  unsigned short GetnMarker_Monitoring(void) const { return nMarker_Monitoring; }

  /*!
   * \brief Get the total number of DV markers.
   * \return Total number of DV markers.
   */
  unsigned short GetnMarker_DV(void) const { return nMarker_DV; }

  /*!
   * \brief Get the total number of moving markers.
   * \return Total number of moving markers.
   */
  unsigned short GetnMarker_Moving(void) const { return nMarker_Moving; }

  /*!
   * \brief Get the total number of markers for gradient treatment.
   * \return Total number of markers for gradient treatment.
   */
  unsigned short GetnMarker_SobolevBC(void) const { return nMarker_SobolevBC; }

  /*!
   * \brief Get the total number of Python customizable markers.
   * \return Total number of Python customizable markers.
   */
  unsigned short GetnMarker_PyCustom(void) const { return nMarker_PyCustom; }

  /*!
   * \brief Get the total number of moving markers.
   * \return Total number of moving markers.
   */
  unsigned short GetnMarker_Analyze(void) const { return nMarker_Analyze; }

  /*!
   * \brief Get the total number of periodic markers.
   * \return Total number of periodic markers.
   */
  unsigned short GetnMarker_Periodic(void) const { return nMarker_PerBound; }

  /*!
   * \brief Get the total (local) number of heat flux markers.
   * \return Total number of heat flux markers.
   */
  unsigned short GetnMarker_HeatFlux(void) const { return nMarker_HeatFlux; }

  /*!
   * \brief Get the total number of rough markers.
   * \return Total number of heat flux markers.
   */
  unsigned short GetnRoughWall(void) const { return nRough_Wall; }

  /*!
   * \brief Get the total number of objectives in kind_objective list
   * \return Total number of objectives in kind_objective list
   */
  unsigned short GetnObj(void) const { return nObj;}

  /*!
   * \brief Stores the number of marker in the simulation.
   * \param[in] val_nmarker - Number of markers of the problem.
   */
  void SetnMarker_All(unsigned short val_nmarker) { nMarker_All = val_nmarker; }

  /*!
   * \brief Get the starting direct iteration number for the unsteady adjoint (reverse time integration).
   * \return Starting direct iteration number for the unsteady adjoint.
   */
  long GetUnst_AdjointIter(void) const { return Unst_AdjointIter; }

  /*!
   * \brief Number of iterations to average (reverse time integration).
   * \return Starting direct iteration number for the unsteady adjoint.
   */
  unsigned long GetIter_Avg_Objective(void) const { return Iter_Avg_Objective ; }

  /*!
   * \brief Retrieves the number of periodic time instances for Harmonic Balance.
   * \return Number of periodic time instances for Harmonic Balance.
   */
  unsigned short GetnTimeInstances(void) const { return nTimeInstances; }

  /*!
   * \brief Retrieves the period of oscillations to be used with Harmonic Balance.
   * \return Period for Harmonic Balance.
   */
  su2double GetHarmonicBalance_Period(void) const { return HarmonicBalance_Period; }

  /*!
   * \brief Set the current external iteration number.
   * \param[in] val_iter - Current external iteration number.
   */
  void SetExtIter_OffSet(unsigned long val_iter) { ExtIter_OffSet = val_iter; }

  /*!
   * \brief Set the current FSI iteration number.
   * \param[in] val_iter - Current FSI iteration number.
   */
  void SetOuterIter(unsigned long val_iter) { OuterIter = val_iter; }

  /*!
   * \brief Set the current FSI iteration number.
   * \param[in] val_iter - Current FSI iteration number.
   */
  void SetInnerIter(unsigned long val_iter) { InnerIter = val_iter; }

  /*!
   * \brief Set the current time iteration number.
   * \param[in] val_iter - Current FSI iteration number.
   */
  void SetTimeIter(unsigned long val_iter) { TimeIter = val_iter; }

  /*!
   * \brief Get the current time iteration number.
   * \param[in] val_iter - Current time iterationnumber.
   */
  unsigned long GetTimeIter() const { return TimeIter; }

  /*!
   * \brief Get the current internal iteration number.
   * \return Current external iteration.
   */
  unsigned long GetExtIter_OffSet(void) const { return ExtIter_OffSet; }

  /*!
   * \brief Get the current FSI iteration number.
   * \return Current FSI iteration.
   */
  unsigned long GetOuterIter(void) const { return OuterIter; }

  /*!
   * \brief Get the current FSI iteration number.
   * \return Current FSI iteration.
   */
  unsigned long GetInnerIter(void) const { return InnerIter; }

  /*!
   * \brief Set the current physical time.
   * \param[in] val_t - Current physical time.
   */
  void SetPhysicalTime(su2double val_t) { PhysicalTime = val_t; }

  /*!
   * \brief Get the current physical time.
   * \return Current physical time.
   */
  su2double GetPhysicalTime(void) const { return PhysicalTime; }

  /*!
   * \brief Get information about writing the performance summary at the end of a calculation.
   * \return <code>TRUE</code> means that the performance summary will be written at the end of a calculation.
   */
  bool GetWrt_Performance(void) const { return Wrt_Performance; }

  /*!
   * \brief Get information about the computational graph (e.g. memory usage) when using AD in reverse mode.
   * \return <code>TRUE</code> means that the tape statistics will be written after each recording.
   */
  bool GetWrt_AD_Statistics(void) const { return Wrt_AD_Statistics; }

  /*!
   * \brief Get information about writing the mesh quality metrics to the visualization files.
   * \return <code>TRUE</code> means that the mesh quality metrics will be written to the visualization files.
   */
  bool GetWrt_MeshQuality(void) const { return Wrt_MeshQuality; }

  /*!
   * \brief Write coarse grids to the visualization files.
   */
  bool GetWrt_MultiGrid(void) const { return Wrt_MultiGrid; }

  /*!
   * \brief Get information about writing projected sensitivities on surfaces to an ASCII file with rows as x, y, z, dJ/dx, dJ/dy, dJ/dz for each vertex.
   * \return <code>TRUE</code> means that projected sensitivities on surfaces in an ASCII file with rows as x, y, z, dJ/dx, dJ/dy, dJ/dz for each vertex will be written.
   */
  bool GetWrt_Projected_Sensitivity(void) const { return Wrt_Projected_Sensitivity; }

  /*!
   * \brief Get information about the format for the input volume sensitvities.
   * \return Format of the input volume sensitivities.
   */
  unsigned short GetSensitivity_Format(void) const { return Sensitivity_FileFormat; }

  /*!
   * \brief Get information about writing sectional force files.
   * \return <code>TRUE</code> means that sectional force files will be written for specified markers.
   */
  bool GetPlot_Section_Forces(void) const { return Plot_Section_Forces; }

  /*!
   * \brief Get the alpha (convective) coefficients for the Runge-Kutta integration scheme.
   * \param[in] val_step - Index of the step.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double Get_Alpha_RKStep(unsigned short val_step) const { return RK_Alpha_Step[val_step]; }

  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_All_TagBound(unsigned short val_marker) const { return Marker_All_TagBound[val_marker]; }

  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_ActDiskInlet_TagBound(unsigned short val_marker) const { return Marker_ActDiskInlet[val_marker]; }

  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_ActDiskOutlet_TagBound(unsigned short val_marker) const { return Marker_ActDiskOutlet[val_marker]; }

  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Outlet_TagBound(unsigned short val_marker) const { return Marker_Outlet[val_marker]; }

  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_EngineInflow_TagBound(unsigned short val_marker) const { return Marker_EngineInflow[val_marker]; }

  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_EngineExhaust_TagBound(unsigned short val_marker) const { return Marker_EngineExhaust[val_marker]; }

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Monitoring_TagBound(unsigned short val_marker) const { return Marker_Monitoring[val_marker]; }

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_HeatFlux_TagBound(unsigned short val_marker) const { return Marker_HeatFlux[val_marker]; }

  /*!
   * \brief Get the tag if the iMarker defined in the geometry file.
   * \param[in] val_tag - Value of the tag in which we are interested.
   * \return Value of the marker <i>val_marker</i> that is in the geometry file
   *         for the surface that has the tag.
   */
  short GetMarker_All_TagBound(string val_tag)  {
    for (unsigned short iMarker = 0; iMarker < nMarker_All; iMarker++) {
      if (val_tag == Marker_All_TagBound[iMarker]) return iMarker;
    }
    return -1;
  }

  /*!
   * \brief Get the kind of boundary for each marker.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \return Kind of boundary for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_All_KindBC(unsigned short val_marker) const { return Marker_All_KindBC[val_marker]; }

  /*!
   * \brief Set the value of the boundary <i>val_boundary</i> (read from the config file)
   *        for the marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_boundary - Kind of boundary read from config file.
   */
  void SetMarker_All_KindBC(unsigned short val_marker, unsigned short val_boundary) { Marker_All_KindBC[val_marker] = val_boundary; }

  /*!
   * \brief Set the value of the index <i>val_index</i> (read from the geometry file) for
   *        the marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_index - Index of the surface read from geometry file.
   */
  void SetMarker_All_TagBound(unsigned short val_marker, string val_index) { Marker_All_TagBound[val_marker] = val_index; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be monitored <i>val_monitoring</i>
   *        (read from the config file).
   * \note This is important for non dimensional coefficient computation.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be monitored.
   */
  void SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring) { Marker_All_Monitoring[val_marker] = val_monitoring; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be monitored <i>val_monitoring</i>
   *        (read from the config file).
   * \note This is important for non dimensional coefficient computation.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be monitored.
   */
  void SetMarker_All_GeoEval(unsigned short val_marker, unsigned short val_geoeval) { Marker_All_GeoEval[val_marker] = val_geoeval; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be designed <i>val_designing</i>
   *        (read from the config file).
   * \note This is important for non dimensional coefficient computation.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be designed.
   */
  void SetMarker_All_Designing(unsigned short val_marker, unsigned short val_designing) { Marker_All_Designing[val_marker] = val_designing; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be plot <i>val_plotting</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_plotting - 0 or 1 depending if the the marker is going to be plot.
   */
  void SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting) { Marker_All_Plotting[val_marker] = val_plotting; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be plot <i>val_plotting</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_plotting - 0 or 1 depending if the the marker is going to be plot.
   */
  void SetMarker_All_Analyze(unsigned short val_marker, unsigned short val_analyze) { Marker_All_Analyze[val_marker] = val_analyze; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is part of the FSI interface <i>val_plotting</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_plotting - 0 or 1 depending if the the marker is part of the FSI interface.
   */
  void SetMarker_All_ZoneInterface(unsigned short val_marker, unsigned short val_fsiinterface) { Marker_All_ZoneInterface[val_marker] = val_fsiinterface; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is part of the Turbomachinery (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_turboperf - 0 if not part of Turbomachinery or greater than 1 if it is part.
   */
  void SetMarker_All_Turbomachinery(unsigned short val_marker, unsigned short val_turbo) { Marker_All_Turbomachinery[val_marker] = val_turbo; }

  /*!
   * \brief Set a flag to the marker <i>val_marker</i> part of the Turbomachinery (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_turboperflag - 0 if is not part of the Turbomachinery, flag INFLOW or OUTFLOW if it is part.
   */
  void SetMarker_All_TurbomachineryFlag(unsigned short val_marker, unsigned short val_turboflag) { Marker_All_TurbomachineryFlag[val_marker] = val_turboflag; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is part of the MixingPlane interface (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_turboperf - 0 if not part of the MixingPlane interface or greater than 1 if it is part.
   */
  void SetMarker_All_MixingPlaneInterface(unsigned short val_marker, unsigned short val_mixpla_interface) { Marker_All_MixingPlaneInterface[val_marker] = val_mixpla_interface; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be affected by design variables <i>val_moving</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_DV - 0 or 1 depending if the the marker is affected by design variables.
   */
  void SetMarker_All_DV(unsigned short val_marker, unsigned short val_DV) { Marker_All_DV[val_marker] = val_DV; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be moved <i>val_moving</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_moving - 0 or 1 depending if the the marker is going to be moved.
   */
  void SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving) { Marker_All_Moving[val_marker] = val_moving; }

  /*!
   * \brief Set if a marker how <i>val_marker</i> is going to be applied in gradient treatment.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_sobolev - 0 or 1 depending if the marker is selected.
   */
  void SetMarker_All_SobolevBC(unsigned short val_marker, unsigned short val_sobolev) { Marker_All_SobolevBC[val_marker] = val_sobolev; }

  /*!
   * \brief Set if a marker <i>val_marker</i> allows deformation at the boundary.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_interface - 0 or 1 depending if the the marker is or not a DEFORM_MESH marker.
   */
  void SetMarker_All_Deform_Mesh(unsigned short val_marker, unsigned short val_deform) { Marker_All_Deform_Mesh[val_marker] = val_deform; }

  /*!
   * \brief Set if a marker <i>val_marker</i> allows deformation at the boundary.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_interface - 0 or 1 depending if the the marker is or not a DEFORM_MESH_SYM_PLANE marker.
   */
  void SetMarker_All_Deform_Mesh_Sym_Plane(unsigned short val_marker, unsigned short val_deform) { Marker_All_Deform_Mesh_Sym_Plane[val_marker] = val_deform; }

  /*!
   * \brief Set if a in marker <i>val_marker</i> the flow load will be computed/employed.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_interface - 0 or 1 depending if the the marker is or not a Fluid_Load marker.
   */
  void SetMarker_All_Fluid_Load(unsigned short val_marker, unsigned short val_interface) { Marker_All_Fluid_Load[val_marker] = val_interface; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be customized in Python <i>val_PyCustom</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_PyCustom - 0 or 1 depending if the the marker is going to be customized in Python.
   */
  void SetMarker_All_PyCustom(unsigned short val_marker, unsigned short val_PyCustom) { Marker_All_PyCustom[val_marker] = val_PyCustom; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be periodic <i>val_perbound</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_perbound - Index of the surface with the periodic boundary.
   */
  void SetMarker_All_PerBound(unsigned short val_marker, short val_perbound) { Marker_All_PerBound[val_marker] = val_perbound; }

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be sent or receive <i>val_index</i>
   *        from another domain.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \param[in] val_index - Index of the surface read from geometry file.
   */
  void SetMarker_All_SendRecv(unsigned short val_marker, short val_index) { Marker_All_SendRecv[val_marker] = val_index; }

  /*!
   * \brief Get the send-receive information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return If positive, the information is sended to that domain, in case negative
   *         the information is receive from that domain.
   */
  short GetMarker_All_SendRecv(unsigned short val_marker) const { return Marker_All_SendRecv[val_marker]; }

  /*!
   * \brief Get an internal index that identify the periodic boundary conditions.
   * \param[in] val_marker - Value of the marker that correspond with the periodic boundary.
   * \return The internal index of the periodic boundary condition.
   */
  short GetMarker_All_PerBound(unsigned short val_marker) const { return Marker_All_PerBound[val_marker]; }

  /*!
   * \brief Get the monitoring information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be monitored.
   * \return 0 or 1 depending if the marker is going to be monitored.
   */
  unsigned short GetMarker_All_Monitoring(unsigned short val_marker) const { return Marker_All_Monitoring[val_marker]; }

  /*!
   * \brief Get the monitoring information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be monitored.
   * \return 0 or 1 depending if the marker is going to be monitored.
   */
  unsigned short GetMarker_All_GeoEval(unsigned short val_marker) const { return Marker_All_GeoEval[val_marker]; }

  /*!
   * \brief Get the design information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be monitored.
   * \return 0 or 1 depending if the marker is going to be monitored.
   */
  unsigned short GetMarker_All_Designing(unsigned short val_marker) const { return Marker_All_Designing[val_marker]; }

  /*!
   * \brief Get the plotting information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is going to be plotted.
   */
  unsigned short GetMarker_All_Plotting(unsigned short val_marker) const { return Marker_All_Plotting[val_marker]; }

  /*!
   * \brief Get the plotting information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is going to be plotted.
   */
  unsigned short GetMarker_All_Analyze(unsigned short val_marker) const { return Marker_All_Analyze[val_marker]; }

  /*!
   * \brief Get the FSI interface information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is part of the FSI interface.
   */
  unsigned short GetMarker_All_ZoneInterface(unsigned short val_marker) const { return Marker_All_ZoneInterface[val_marker]; }

  /*!
   * \brief Get the MixingPlane interface information for a marker <i>val_marker</i>.
   * \param[in] val_marker value of the marker on the grid.
   * \return 0 if is not part of the MixingPlane Interface and greater than 1 if it is part.
   */
  unsigned short GetMarker_All_MixingPlaneInterface(unsigned short val_marker) const { return Marker_All_MixingPlaneInterface[val_marker]; }

  /*!
   * \brief Get the Turbomachinery information for a marker <i>val_marker</i>.
   * \param[in] val_marker value of the marker on the grid.
   * \return 0 if is not part of the Turbomachinery and greater than 1 if it is part.
   */
  unsigned short GetMarker_All_Turbomachinery(unsigned short val_marker) const { return Marker_All_Turbomachinery[val_marker]; }

  /*!
   * \brief Get the Turbomachinery flag information for a marker <i>val_marker</i>.
   * \param[in] val_marker value of the marker on the grid.
   * \return 0 if is not part of the Turbomachinery, flag INFLOW or OUTFLOW if it is part.
   */
  unsigned short GetMarker_All_TurbomachineryFlag(unsigned short val_marker) const { return Marker_All_TurbomachineryFlag[val_marker]; }

  /*!
   * \brief Get the number of FSI interface markers <i>val_marker</i>.
   * \param[in] void.
   * \return Number of markers belonging to the FSI interface.
   */
  unsigned short GetMarker_n_ZoneInterface(void) const { return nMarker_ZoneInterface; }

  /*!
   * \brief Get the DV information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be affected by design variables.
   * \return 0 or 1 depending if the marker is going to be affected by design variables.
   */
  unsigned short GetMarker_All_DV(unsigned short val_marker) const { return Marker_All_DV[val_marker]; }

  /*!
   * \brief Get the motion information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is going to be moved.
   */
  unsigned short GetMarker_All_Moving(unsigned short val_marker) const { return Marker_All_Moving[val_marker]; }

  /*!
   * \brief Get the information if gradient treatment uses a marker <i>val_marker</i>.
   * \param[in] val_marker
   * \return 0 or 1 depending if the marker is going to be selected.
   */
  unsigned short GetMarker_All_SobolevBC(unsigned short val_marker) const { return Marker_All_SobolevBC[val_marker]; }

  /*!
   * \brief Get whether marker <i>val_marker</i> is a DEFORM_MESH marker
   * \param[in] val_marker - 0 or 1 depending if the the marker belongs to the DEFORM_MESH subset.
   * \return 0 or 1 depending if the marker belongs to the DEFORM_MESH subset.
   */
  unsigned short GetMarker_All_Deform_Mesh(unsigned short val_marker) const { return Marker_All_Deform_Mesh[val_marker]; }

  /*!
   * \brief Get whether marker <i>val_marker</i> is a DEFORM_MESH_SYM_PLANE marker
   * \param[in] val_marker - 0 or 1 depending if the the marker belongs to the DEFORM_MESH_SYM_PLANE subset.
   * \return 0 or 1 depending if the marker belongs to the DEFORM_MESH_SYM_PLANE subset.
   */
  unsigned short GetMarker_All_Deform_Mesh_Sym_Plane(unsigned short val_marker) const { return Marker_All_Deform_Mesh_Sym_Plane[val_marker]; }

  /*!
   * \brief Get whether marker <i>val_marker</i> is a Fluid_Load marker
   * \param[in] val_marker - 0 or 1 depending if the the marker belongs to the Fluid_Load subset.
   * \return 0 or 1 depending if the marker belongs to the Fluid_Load subset.
   */
  unsigned short GetMarker_All_Fluid_Load(unsigned short val_marker) const { return Marker_All_Fluid_Load[val_marker]; }

  /*!
   * \brief Get the Python customization for a marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \return 0 or 1 depending if the marker is going to be customized in Python.
   */
  unsigned short GetMarker_All_PyCustom(unsigned short val_marker) const { return Marker_All_PyCustom[val_marker];}

  /*!
   * \brief Get the airfoil sections in the slicing process.
   * \param[in] val_section - Index of the section.
   * \return Coordinate of the airfoil to slice.
   */
  su2double GetLocationStations(unsigned short val_section) const { return LocationStations[val_section]; }

  /*!
   * \brief Get the defintion of the nacelle location.
   * \param[in] val_index - Index of the section.
   * \return Coordinate of the nacelle location.
   */
  su2double GetNacelleLocation(unsigned short val_index) const { return nacelle_location[val_index]; }

  /*!
   * \brief Get the number of pre-smoothings in a multigrid strategy.
   * \param[in] val_mesh - Index of the grid.
   * \return Number of smoothing iterations.
   */
  unsigned short GetMG_PreSmooth(unsigned short val_mesh) const {
    if (nMG_PreSmooth == 0) return 1;
    return MG_PreSmooth[val_mesh];
  }

  /*!
   * \brief Get the number of post-smoothings in a multigrid strategy.
   * \param[in] val_mesh - Index of the grid.
   * \return Number of smoothing iterations.
   */
  unsigned short GetMG_PostSmooth(unsigned short val_mesh) const {
    if (nMG_PostSmooth == 0) return 0;
    return MG_PostSmooth[val_mesh];
  }

  /*!
   * \brief Get the number of implicit Jacobi smoothings of the correction in a multigrid strategy.
   * \param[in] val_mesh - Index of the grid.
   * \return Number of implicit smoothing iterations.
   */
  unsigned short GetMG_CorrecSmooth(unsigned short val_mesh) const {
    if (nMG_CorrecSmooth == 0) return 0;
    return MG_CorrecSmooth[val_mesh];
  }

  /*!
   * \brief plane of the FFD (I axis) that should be fixed.
   * \param[in] val_index - Index of the arrray with all the planes in the I direction that should be fixed.
   * \return Index of the plane that is going to be freeze.
   */
  short GetFFD_Fix_IDir(unsigned short val_index) const { return FFD_Fix_IDir[val_index]; }

  /*!
   * \brief plane of the FFD (J axis) that should be fixed.
   * \param[in] val_index - Index of the arrray with all the planes in the J direction that should be fixed.
   * \return Index of the plane that is going to be freeze.
   */
  short GetFFD_Fix_JDir(unsigned short val_index) const { return FFD_Fix_JDir[val_index]; }

  /*!
   * \brief plane of the FFD (K axis) that should be fixed.
   * \param[in] val_index - Index of the arrray with all the planes in the K direction that should be fixed.
   * \return Index of the plane that is going to be freeze.
   */
  short GetFFD_Fix_KDir(unsigned short val_index) const { return FFD_Fix_KDir[val_index]; }

  /*!
   * \brief Get the number of planes to fix in the I direction.
   * \return Number of planes to fix in the I direction.
   */
  unsigned short GetnFFD_Fix_IDir(void) const { return nFFD_Fix_IDir; }

  /*!
   * \brief Get the number of planes to fix in the J direction.
   * \return Number of planes to fix in the J direction.
   */
  unsigned short GetnFFD_Fix_JDir(void) const { return nFFD_Fix_JDir; }

  /*!
   * \brief Get the number of planes to fix in the K direction.
   * \return Number of planes to fix in the K direction.
   */
  unsigned short GetnFFD_Fix_KDir(void) const { return nFFD_Fix_KDir; }

  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  MAIN_SOLVER GetKind_Solver(void) const { return Kind_Solver; }

  /*!
   * \brief Return true if a fluid solver is in use.
   */
  bool GetFluidProblem(void) const {
    switch (Kind_Solver) {
      case MAIN_SOLVER::EULER : case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
      case MAIN_SOLVER::INC_EULER : case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
      case MAIN_SOLVER::NEMO_EULER : case MAIN_SOLVER::NEMO_NAVIER_STOKES:
      case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
        return true;
      default:
        return false;
    }
  }

  /*!
   * \brief Return true if a structural solver is in use.
   */
  bool GetStructuralProblem(void) const {
    return (Kind_Solver == MAIN_SOLVER::FEM_ELASTICITY) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM);
  }

  /*!
   * \brief Return true if a heat solver is in use.
   */
  bool GetHeatProblem(void) const {
    return (Kind_Solver == MAIN_SOLVER::HEAT_EQUATION) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_HEAT);
  }

  /*!
   * \brief Return true if a high order FEM solver is in use.
   */
  bool GetFEMSolver(void) const {
    switch (Kind_Solver) {
      case MAIN_SOLVER::FEM_EULER: case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_RANS: case MAIN_SOLVER::FEM_LES:
      case MAIN_SOLVER::DISC_ADJ_FEM_EULER: case MAIN_SOLVER::DISC_ADJ_FEM_NS: case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
        return true;
      default:
        return false;
    }
  }

  /*!
   * \brief Return true if a NEMO solver is in use.
   */
  bool GetNEMOProblem(void) const {
    switch (Kind_Solver) {
      case MAIN_SOLVER::NEMO_EULER : case MAIN_SOLVER::NEMO_NAVIER_STOKES:
        return true;
      default:
        return false;
    }
  }

   /*!
   * \brief Return true if an AUSM method is in use.
   */
  bool GetAUSMMethod(void) const {
    switch (Kind_Upwind_Flow) {
      case UPWIND::AUSM : case UPWIND::AUSMPLUSUP: case UPWIND::AUSMPLUSUP2: case UPWIND::AUSMPLUSM:
        return true;
      default:
        return false;
    }
  }

  /*!
   * \brief Kind of Multizone Solver.
   * \return Governing equation that we are solving.
   */
  ENUM_MULTIZONE GetKind_MZSolver(void) const { return Kind_MZSolver; }

  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  ENUM_REGIME GetKind_Regime(void) const { return Kind_Regime; }

  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  unsigned short GetSystemMeasurements(void) const { return SystemMeasurements; }

  /*!
   * \brief Gas model that we are using.
   * \return Gas model that we are using.
   */
  string GetGasModel(void) const {return GasModel;}

  /*!
   * \brief Get the transport coefficient model.
   * \return Index of transport coefficient model.
   */
  TRANSCOEFFMODEL GetKind_TransCoeffModel(void) const { return Kind_TransCoeffModel; }

  /*!
   * \brief Get the total number of heat flux markers.
   * \return Total number of heat flux markers.
   */
  unsigned short GetnWall_Catalytic(void) const { return nWall_Catalytic; }

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetWall_Catalytic_TagBound(unsigned short val_marker) const { return Wall_Catalytic[val_marker]; }

  /*!
   * \brief Get wall catalytic efficiency.
   * \return wall catalytic efficiency value.
   */
  su2double GetCatalytic_Efficiency(void) const { return CatalyticEfficiency; }

  /*!
   * \brief Fluid model that we are using.
   * \return Fluid model that we are using.
   */
  unsigned short GetKind_FluidModel(void) const { return Kind_FluidModel; }

  /*!
   * \brief Datadriven method for EoS evaluation.
   */
  ENUM_DATADRIVEN_METHOD GetKind_DataDriven_Method(void) const { return Kind_DataDriven_Method; }

  /*!
   * \brief Get name of the input file for the data-driven fluid model interpolation method.
   * \return Name of the input file for the interpolation method.
   */
  const string* GetDataDriven_FileNames(void) const { return DataDriven_Method_FileNames; }

  /*!
   * \brief Get number of listed look-up table or multi-layer perceptron input files.
   * \return Number of listed data-driven method input files.
   */
  unsigned short GetNDataDriven_Files(void) const { return n_Datadriven_files; }

  /*!
   * \brief Get Newton solver relaxation factor for data-driven fluid models.
   * \return Newton solver relaxation factor.
   */
  su2double GetRelaxation_DataDriven(void) const { return DataDriven_Relaxation_Factor; }

  /*!
   * \brief Returns the name of the fluid we are using in CoolProp.
   */
  string GetFluid_Name(void) const { return FluidName; }

  /*!
   * \brief Option to define the density model for incompressible flows.
   * \return Density model option
   */
  INC_DENSITYMODEL GetKind_DensityModel() const { return Kind_DensityModel; }

  /*!
   * \brief Selection of variable density option for incompressible flows.
   * \return Flag for variable density for incompressible flows.
   */
  bool GetVariable_Density_Model() const { return Variable_Density; }

  /*!
   * \brief Flag for whether to solve the energy equation for incompressible flows.
   * \return Flag for energy equation
   */
  bool GetEnergy_Equation(void) const { return Energy_Equation; }

  /*!
   * \brief free stream option to initialize the solution
   * \return free stream option
   */
  FREESTREAM_OPTION GetKind_FreeStreamOption() const { return Kind_FreeStreamOption; }

  /*!
   * \brief free stream option to initialize the solution
   * \return free stream option
   */
  unsigned short GetKind_InitOption(void) const { return Kind_InitOption; }
  /*!
   * \brief Get the value of the critical pressure.
   * \return Critical pressure.
   */
  su2double GetPressure_Critical(void) const { return Pressure_Critical; }

  /*!
   * \brief Get the value of the critical temperature.
   * \return Critical temperature.
   */
  su2double GetTemperature_Critical(void) const { return Temperature_Critical; }

  /*!
   * \brief Get the value of the critical pressure.
   * \return Critical pressure.
   */
  su2double GetAcentric_Factor(void) const { return Acentric_Factor; }

  /*!
   * \brief Get the value of the viscosity model.
   * \return Viscosity model.
   */
  VISCOSITYMODEL GetKind_ViscosityModel() const { return Kind_ViscosityModel; }

  /*!
   * \brief Get the value of the mixing model for viscosity.
   * \return Mixing Viscosity model.
   */
  MIXINGVISCOSITYMODEL GetKind_MixingViscosityModel() const { return Kind_MixingViscosityModel; }

  /*!
   * \brief Get the value of the thermal conductivity model.
   * \return Conductivity model.
   */
  CONDUCTIVITYMODEL GetKind_ConductivityModel() const { return Kind_ConductivityModel; }

  /*!
   * \brief Get the value of the turbulent thermal conductivity model.
   * \return Turbulent conductivity model.
   */
  CONDUCTIVITYMODEL_TURB GetKind_ConductivityModel_Turb() const { return Kind_ConductivityModel_Turb; }

  /*!
   * \brief Get the value of the mass diffusivity model.
   * \return Mass diffusivity model.
   */
  DIFFUSIVITYMODEL GetKind_Diffusivity_Model(void) const { return Kind_Diffusivity_Model; }

  /*!
   * \brief Get the value of the constant viscosity.
   * \return Constant viscosity.
   */
  su2double GetMu_Constant(unsigned short val_index = 0) const { return Mu_Constant[val_index]; }

  /*!
   * \brief Get the value of the non-dimensional constant viscosity.
   * \return Non-dimensional constant viscosity.
   */
  su2double GetMu_ConstantND(unsigned short val_index = 0) const { return Mu_Constant[val_index] / Viscosity_Ref; }

  /*!
   * \brief Get the value of the thermal conductivity.
   * \return Thermal conductivity.
   */
  su2double GetThermal_Conductivity_Constant(unsigned short val_index = 0) const {
    return Thermal_Conductivity_Constant[val_index];
  }

  /*!
   * \brief Get the value of the non-dimensional thermal conductivity.
   * \return Non-dimensional thermal conductivity.
   */
  su2double GetThermal_Conductivity_ConstantND(unsigned short val_index = 0) const {
    return Thermal_Conductivity_Constant[val_index] / Thermal_Conductivity_Ref;
  }

  /*!
   * \brief Get the value of the constant mass diffusivity for scalar transport.
   * \return Constant mass diffusivity.
   */
  su2double GetDiffusivity_Constant(void) const { return Diffusivity_Constant; }

  /*!
   * \brief Get the value of the non-dimensional constant mass diffusivity.
   * \return Non-dimensional constant mass diffusivity.
   */
  su2double GetDiffusivity_ConstantND(void) const { return Diffusivity_ConstantND; }

  /*!
   * \brief Get the value of the laminar Schmidt number for scalar transport.
   * \return Laminar Schmidt number for scalar transport.
   */
  su2double GetSchmidt_Number_Laminar(void) const { return Schmidt_Number_Laminar; }

  /*!
   * \brief Get the value of the turbulent Schmidt number for scalar transport.
   * \return Turbulent Schmidt number for scalar transport.
   */
  su2double GetSchmidt_Number_Turbulent(void) const { return Schmidt_Number_Turbulent; }

  /*!
   * \brief Get the value of the Lewis number for each species.
   * \return Lewis Number.
   */
  su2double GetConstant_Lewis_Number(unsigned short val_index = 0) const { return Constant_Lewis_Number[val_index]; }

  /*!
   * \brief Get the value of the reference viscosity for Sutherland model.
   * \return The reference viscosity.
   */
  su2double GetMu_Ref(unsigned short val_index = 0) const { return Mu_Ref[val_index]; }

  /*!
   * \brief Get the value of the non-dimensional reference viscosity for Sutherland model.
   * \return The non-dimensional reference viscosity.
   */
  su2double GetMu_RefND(unsigned short val_index = 0) const { return Mu_Ref[val_index] / Viscosity_Ref; }

  /*!
   * \brief Get the value of the reference temperature for Sutherland model.
   * \return The reference temperature.
   */
  su2double GetMu_Temperature_Ref(unsigned short val_index = 0) const { return Mu_Temperature_Ref[val_index]; }

  /*!
   * \brief Get the value of the non-dimensional reference temperature for Sutherland model.
   * \return The non-dimensional reference temperature.
   */
  su2double GetMu_Temperature_RefND(unsigned short val_index = 0) const {
    return Mu_Temperature_Ref[val_index] / Temperature_Ref;
  }

  /*!
   * \brief Get the value of the reference S for Sutherland model.
   * \return The reference S.
   */
  su2double GetMu_S(unsigned short val_index = 0) const { return Mu_S[val_index]; }

  /*!
   * \brief Get the value of the non-dimensional reference S for Sutherland model.
   * \return The non-dimensional reference S.
   */
  su2double GetMu_SND(unsigned short val_index = 0) const { return Mu_S[val_index] / Temperature_Ref; }

  /*!
   * \brief Get the number of coefficients in the temperature polynomial models.
   * \return The the number of coefficients in the temperature polynomial models.
   */
  unsigned short GetnPolyCoeffs(void) const { return N_POLY_COEFFS; }

  /*!
   * \brief Get the temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for specific heat Cp.
   */
  su2double GetCp_PolyCoeff(unsigned short val_index) const { return cp_polycoeffs[val_index]; }

  /*!
   * \brief Get the temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for specific heat Cp.
   */
  su2double GetCp_PolyCoeffND(unsigned short val_index) const { return CpPolyCoefficientsND[val_index]; }

  /*!
   * \brief Get the temperature polynomial coefficient for viscosity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for viscosity.
   */
  su2double GetMu_PolyCoeff(unsigned short val_index) const { return mu_polycoeffs[val_index]; }

  /*!
   * \brief Get the temperature polynomial coefficient for viscosity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Non-dimensional temperature polynomial coefficient for viscosity.
   */
  su2double GetMu_PolyCoeffND(unsigned short val_index) const { return MuPolyCoefficientsND[val_index]; }

  /*!
   * \brief Get the temperature polynomial coefficients for viscosity.
   * \return Non-dimensional temperature polynomial coefficients for viscosity.
   */
  const su2double* GetMu_PolyCoeffND(void) const { return MuPolyCoefficientsND.data(); }

  /*!
   * \brief Get the temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for thermal conductivity.
   */
  su2double GetKt_PolyCoeff(unsigned short val_index) const { return kt_polycoeffs[val_index]; }

  /*!
   * \brief Get the temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Non-dimensional temperature polynomial coefficient for thermal conductivity.
   */
  su2double GetKt_PolyCoeffND(unsigned short val_index) const { return KtPolyCoefficientsND[val_index]; }

  /*!
   * \brief Get the temperature polynomial coefficients for thermal conductivity.
   * \return Non-dimensional temperature polynomial coefficients for thermal conductivity.
   */
  const su2double* GetKt_PolyCoeffND(void) const { return KtPolyCoefficientsND.data(); }

  /*!
   * \brief Set the temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_coeff - Temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   */
  void SetCp_PolyCoeffND(su2double val_coeff, unsigned short val_index) { CpPolyCoefficientsND[val_index] = val_coeff; }

  /*!
   * \brief Set the temperature polynomial coefficient for viscosity.
   * \param[in] val_coeff - Non-dimensional temperature polynomial coefficient for viscosity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   */
  void SetMu_PolyCoeffND(su2double val_coeff, unsigned short val_index) { MuPolyCoefficientsND[val_index] = val_coeff; }

  /*!
   * \brief Set the temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_coeff - Non-dimensional temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   */
  void SetKt_PolyCoeffND(su2double val_coeff, unsigned short val_index) { KtPolyCoefficientsND[val_index] = val_coeff; }

  /*!
   * \brief Set the value of the non-dimensional constant mass diffusivity.
   */
  void SetDiffusivity_ConstantND(su2double diffusivity_const) { Diffusivity_ConstantND = diffusivity_const; }

  /*!
   * \brief Get the kind of method for computation of spatial gradients used for viscous and source terms.
   * \return Numerical method for computation of spatial gradients used for viscous and source terms.
   */
  unsigned short GetKind_Gradient_Method(void) const { return Kind_Gradient_Method; }

  /*!
   * \brief Get the kind of method for computation of spatial gradients used for upwind reconstruction.
   * \return Numerical method for computation of spatial gradients used for upwind reconstruction.
   */
  unsigned short GetKind_Gradient_Method_Recon(void) const { return Kind_Gradient_Method_Recon; }

  /*!
   * \brief Get flag for whether a second gradient calculation is required for upwind reconstruction alone.
   * \return <code>TRUE</code> means that a second gradient will be calculated for upwind reconstruction.
   */
  bool GetReconstructionGradientRequired(void) const { return ReconstructionGradientRequired; }

  /*!
   * \brief Get flag for whether a least-squares gradient method is being applied.
   * \return <code>TRUE</code> means that a least-squares gradient method is being applied.
   */
  bool GetLeastSquaresRequired(void) const { return LeastSquaresRequired; }

  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Linear_Solver(void) const { return Kind_Linear_Solver; }


  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Linear_Solver_Prec(void) const { return Kind_Linear_Solver_Prec; }

  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Deform_Linear_Solver(void) const { return Kind_Deform_Linear_Solver; }

  /*!
   * \brief Get min error of the linear solver for the implicit formulation.
   * \return Min error of the linear solver for the implicit formulation.
   */
  su2double GetLinear_Solver_Error(void) const { return Linear_Solver_Error; }

  /*!
   * \brief Get min error of the linear solver for the implicit formulation.
   * \return Min error of the linear solver for the implicit formulation.
   */
  su2double GetDeform_Linear_Solver_Error(void) const { return Deform_Linear_Solver_Error; }

  /*!
   * \brief Get max number of iterations of the linear solver for the implicit formulation.
   * \return Max number of iterations of the linear solver for the implicit formulation.
   */
  unsigned long GetLinear_Solver_Iter(void) const { return Linear_Solver_Iter; }

  /*!
   * \brief Get max number of iterations of the linear solver for the implicit formulation.
   * \return Max number of iterations of the linear solver for the implicit formulation.
   */
  unsigned long GetDeform_Linear_Solver_Iter(void) const { return Deform_Linear_Solver_Iter; }

  /*!
   * \brief Get the ILU fill-in level for the linear solver.
   * \return Fill in level of the ILU preconditioner for the linear solver.
   */
  unsigned short GetLinear_Solver_ILU_n(void) const { return Linear_Solver_ILU_n; }

  /*!
   * \brief Get restart frequency of the linear solver for the implicit formulation.
   * \return Restart frequency of the linear solver for the implicit formulation.
   */
  unsigned long GetLinear_Solver_Restart_Frequency(void) const { return Linear_Solver_Restart_Frequency; }

  /*!
   * \brief Get the relaxation factor for iterative linear smoothers.
   * \return Relaxation factor.
   */
  su2double GetLinear_Solver_Smoother_Relaxation(void) const { return Linear_Solver_Smoother_Relaxation; }

  /*!
   * \brief Get the relaxation factor for solution updates of adjoint solvers.
   */
  su2double GetRelaxation_Factor_Adjoint(void) const { return Relaxation_Factor_Adjoint; }

  /*!
   * \brief Get the relaxation coefficient of the CHT coupling.
   * \return relaxation coefficient of the CHT coupling.
   */
  su2double GetRelaxation_Factor_CHT(void) const { return Relaxation_Factor_CHT; }

  /*!
   * \brief Get the number of samples used in quasi-Newton methods.
   */
  unsigned short GetnQuasiNewtonSamples(void) const { return nQuasiNewtonSamples; }

  /*!
   * \brief Get whether to use vectorized numerics (if available).
   */
  bool GetUseVectorization(void) const { return UseVectorization; }

  /*!
   * \brief Get whether to use a Newton-Krylov method.
   */
  bool GetNewtonKrylov(void) const { return NewtonKrylov; }

  /*!
   * \brief Get Newton-Krylov integer parameters.
   */
  array<unsigned short,3> GetNewtonKrylovIntParam(void) const { return NK_IntParam; }

  /*!
   * \brief Get Newton-Krylov floating-point parameters.
   */
  array<su2double,4> GetNewtonKrylovDblParam(void) const { return NK_DblParam; }

  /*!
   * \brief Get the relaxation coefficient of the linear solver for the implicit formulation.
   * \return relaxation coefficient of the linear solver for the implicit formulation.
   */
  su2double GetRoe_Kappa(void) const { return Roe_Kappa; }

  /*!
   * \brief Get the wing semi span.
   * \return value of the wing semi span.
   */
  su2double GetSemiSpan(void) const { return SemiSpan; }

  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_AdjTurb_Linear_Solver(void) const { return Kind_AdjTurb_Linear_Solver; }

  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_AdjTurb_Linear_Prec(void) const { return Kind_AdjTurb_Linear_Prec; }

  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_DiscAdj_Linear_Solver(void) const { return Kind_DiscAdj_Linear_Solver; }

  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_DiscAdj_Linear_Prec(void) const { return Kind_DiscAdj_Linear_Prec; }

  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Deform_Linear_Solver_Prec(void) const { return Kind_Deform_Linear_Solver_Prec; }

  /*!
   * \brief Set the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  void SetKind_AdjTurb_Linear_Prec(unsigned short val_kind_prec) { Kind_AdjTurb_Linear_Prec = val_kind_prec; }

  /*!
   * \brief Get min error of the linear solver for the implicit formulation.
   * \return Min error of the linear solver for the implicit formulation.
   */
  su2double GetAdjTurb_Linear_Error(void) const { return AdjTurb_Linear_Error; }

  /*!
   * \brief Get the entropy fix.
   * \return Vaule of the entropy fix.
   */
  su2double GetEntropyFix_Coeff(void) const { return EntropyFix_Coeff; }

  /*!
   * \brief Get max number of iterations of the linear solver for the implicit formulation.
   * \return Max number of iterations of the linear solver for the implicit formulation.
   */
  unsigned short GetAdjTurb_Linear_Iter(void) const { return AdjTurb_Linear_Iter; }

  /*!
   * \brief Get CFL reduction factor for adjoint turbulence model.
   * \return CFL reduction factor.
   */
  su2double GetCFLRedCoeff_AdjTurb(void) const { return CFLRedCoeff_AdjTurb; }

  /*!
   * \brief Get the number of nonlinear increments for mesh deformation.
   * \return Number of nonlinear increments for mesh deformation.
   */
  unsigned long GetGridDef_Nonlinear_Iter(void) const { return GridDef_Nonlinear_Iter; }

  /*!
   * \brief Get information about whether the mesh will be deformed using pseudo linear elasticity.
   * \return <code>TRUE</code> means that grid deformation is active.
   */
  bool GetDeform_Mesh(void) const { return Deform_Mesh; }

  /*!
   * \brief Get information about writing grid deformation residuals to the console.
   * \return <code>TRUE</code> means that grid deformation residuals will be written to the console.
   */
  bool GetDeform_Output(void) const { return Deform_Output; }

  /*!
   * \brief Get factor to multiply smallest volume for deform tolerance.
   * \return Factor to multiply smallest volume for deform tolerance.
   */
  su2double GetDeform_Coeff(void) const { return Deform_Coeff; }

  /*!
   * \brief Get limit for the volumetric deformation.
   * \return Distance to the surface to be deformed.
   */
  su2double GetDeform_Limit(void) const { return Deform_Limit; }

  /*!
   * \brief Get Young's modulus for deformation (constant stiffness deformation)
   */
  su2double GetDeform_ElasticityMod(void) const { return Deform_ElasticityMod; }

  /*!
   * \brief Get Poisson's ratio for deformation (constant stiffness deformation)
   * \
   */
  su2double GetDeform_PoissonRatio(void) const { return Deform_PoissonRatio; }

  /*!
   * \brief Get the type of stiffness to impose for FEA mesh deformation.
   * \return type of stiffness to impose for FEA mesh deformation.
   */
  unsigned short GetDeform_Stiffness_Type(void) const { return Deform_StiffnessType; }

  /*!
   * \brief Get the size of the layer of highest stiffness for wall distance-based mesh stiffness.
   */
  su2double GetDeform_StiffLayerSize(void) const { return Deform_StiffLayerSize; }

  /*!
   * \brief Define the FFD box with a symetry plane.
   * \return <code>TRUE</code> if there is a symmetry plane in the FFD; otherwise <code>FALSE</code>.
   */
  bool GetFFD_Symmetry_Plane(void) const { return FFD_Symmetry_Plane; }

  /*!
   * \brief Get the kind of SU2 software component.
   * \return Kind of the SU2 software component.
   */
  SU2_COMPONENT GetKind_SU2(void) const { return Kind_SU2; }

  /*!
   * \brief Get the kind of non-dimensionalization.
   * \return Kind of non-dimensionalization.
   */
  unsigned short GetRef_NonDim(void) const { return Ref_NonDim; }

  /*!
   * \brief Get the kind of incompressible non-dimensionalization.
   * \return Kind of incompressible non-dimensionalization.
   */
  unsigned short GetRef_Inc_NonDim(void) const { return Ref_Inc_NonDim; }

  /*!
   * \brief Set the kind of SU2 software component.
   * \return Kind of the SU2 software component.
   */
  void SetKind_SU2(SU2_COMPONENT val_kind_su2) { Kind_SU2 = val_kind_su2 ; }

  /*!
   * \brief Get the number of Turbulence Variables.
   * \return Number of Turbulence Variables.
   */
  unsigned short GetnTurbVar(void) const { return nTurbVar; }

  /*!
   * \brief Get the kind of the turbulence model.
   * \return Kind of the turbulence model.
   */
  TURB_MODEL GetKind_Turb_Model(void) const { return Kind_Turb_Model; }

  /*!
   * \brief Get the kind of the transition model.
   * \return Kind of the transion model.
   */
  TURB_TRANS_MODEL GetKind_Trans_Model(void) const { return Kind_Trans_Model; }

  /*!
   * \brief Get the kind of the transition correlations.
   * \return Kind of the transition correlation.
   */
  TURB_TRANS_CORRELATION GetKind_Trans_Correlation(void) const { return Kind_Trans_Correlation; }

  /*!
   * \brief Get RMS roughness for Transtion model from config
   * \return Value of roughness.
   */
  su2double GethRoughness(void) const { return hRoughness; }

  /*!
   * \brief Get the kind of the species model.
   * \return Kind of the species model.
   */
  SPECIES_MODEL GetKind_Species_Model(void) const { return Kind_Species_Model; }

  /*!
   * \brief Get the kind of the subgrid scale model.
   * \return Kind of the subgrid scale model.
   */
  TURB_SGS_MODEL GetKind_SGS_Model(void) const { return Kind_SGS_Model; }

  /*!
   * \brief Get the kind of time integration method.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of time integration method.
   */
  unsigned short GetKind_TimeIntScheme(void) const { return Kind_TimeNumScheme; }

  /*!
   * \brief Get the kind of convective numerical scheme.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of the convective scheme.
   */
  unsigned short GetKind_ConvNumScheme(void) const { return Kind_ConvNumScheme; }

  /*!
   * \brief Get kind of center scheme for the convective terms.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of center scheme for the convective terms.
   */
  CENTERED GetKind_Centered(void) const { return Kind_Centered; }

  /*!
   * \brief Get kind of upwind scheme for the convective terms.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of upwind scheme for the convective terms.
   */
  UPWIND GetKind_Upwind(void) const { return Kind_Upwind; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL(void) const { return MUSCL; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Flow(void) const { return MUSCL_Flow; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Heat(void) const { return MUSCL_Heat; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Turb(void) const { return MUSCL_Turb; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Species(void) const { return MUSCL_Species; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_AdjFlow(void) const { return MUSCL_AdjFlow; }

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_AdjTurb(void) const { return MUSCL_AdjTurb; }

  /*!
   * \brief Get whether to "Use Accurate Jacobians" for AUSM+up(2) and SLAU(2).
   * \return yes/no.
   */
  bool GetUse_Accurate_Jacobians(void) const { return Use_Accurate_Jacobians; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the flow equations.
   */
  unsigned short GetKind_TimeIntScheme_Flow(void) const { return Kind_TimeIntScheme_Flow; }

  /*!
   * \brief Get the kind of scheme (aliased or non-aliased) to be used in the
   *        predictor step of ADER-DG.
   * \return Kind of scheme used in the predictor step of ADER-DG.
   */
  unsigned short GetKind_ADER_Predictor(void) const { return Kind_ADER_Predictor; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the plasma equations.
   */
  unsigned short GetKind_TimeIntScheme_Heat(void) const { return Kind_TimeIntScheme_Heat; }

  /*!
   * \brief Get the kind of time stepping
   *        for the heat equation.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of time stepping for the heat equation.
   */
  unsigned short GetKind_TimeStep_Heat(void) const { return Kind_TimeStep_Heat; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the plasma equations.
   */
  STRUCT_TIME_INT GetKind_TimeIntScheme_FEA(void) const { return Kind_TimeIntScheme_FEA; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the radiation equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the radiation equations.
   */
  unsigned short GetKind_TimeIntScheme_Radiation(void) const { return Kind_TimeIntScheme_Radiation; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the template equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the plasma equations.
   */
  unsigned short GetKind_TimeIntScheme_Template(void);

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the plasma equations.
   */
  STRUCT_SPACE_ITE GetKind_SpaceIteScheme_FEA(void) const { return Kind_SpaceIteScheme_FEA; }

  /*!
   * \brief Get the kind of convective numerical scheme for the flow
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_ConvNumScheme_Flow(void) const { return Kind_ConvNumScheme_Flow; }

  /*!
   * \brief Get the kind of convective numerical scheme for the flow
   *        equations (finite element).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_ConvNumScheme_FEM_Flow(void) const { return Kind_ConvNumScheme_FEM_Flow; }

  /*!
   * \brief Get the kind of convective numerical scheme for the template
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_ConvNumScheme_Template(void) const { return Kind_ConvNumScheme_Template; }

  /*!
   * \brief Get the kind of center convective numerical scheme for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the flow equations.
   */
  CENTERED GetKind_Centered_Flow(void) const { return Kind_Centered_Flow; }

  /*!
   * \brief Get the kind of center convective numerical scheme for the plasma equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_Centered_Template(void);

  /*!
   * \brief Get the kind of upwind convective numerical scheme for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the flow equations.
   */
  UPWIND GetKind_Upwind_Flow(void) const { return Kind_Upwind_Flow; }

  /*!
   * \brief Get the kind of finite element convective numerical scheme for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of finite element convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_FEM_Flow(void) const { return Kind_FEM_Flow; }

  /*!
   * \brief Get the kind of shock capturing method in FEM DG solver.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of shock capturing method in FEM DG solver.
   */
  FEM_SHOCK_CAPTURING_DG GetKind_FEM_DG_Shock(void) const { return Kind_FEM_Shock_Capturing_DG; }

  /*!
   * \brief Get the kind of matrix coloring used for the sparse Jacobian computation.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of matrix coloring used.
   */
  unsigned short GetKind_Matrix_Coloring(void) const { return Kind_Matrix_Coloring; }

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients.
   */
  LIMITER GetKind_SlopeLimit(void) const { return Kind_SlopeLimit; }

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the flow equations.
   */
  LIMITER GetKind_SlopeLimit_Flow(void) const { return Kind_SlopeLimit_Flow; }

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the turbulent equation.
   */
  LIMITER GetKind_SlopeLimit_Turb(void) const { return Kind_SlopeLimit_Turb; }

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the species equation.
   */
  LIMITER GetKind_SlopeLimit_Species() const { return Kind_SlopeLimit_Species; }

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the adjoint turbulent equation.
   */
  LIMITER GetKind_SlopeLimit_AdjTurb(void) const { return Kind_SlopeLimit_AdjTurb; }

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the adjoint flow equation.
   */
  LIMITER GetKind_SlopeLimit_AdjFlow(void) const { return Kind_SlopeLimit_AdjFlow; }

  /*!
   * \brief Value of the calibrated constant for the Lax method (center scheme).
   * \note This constant is used in coarse levels and with first order methods.
   * \return Calibrated constant for the Lax method.
   */
  su2double GetKappa_1st_Flow(void) const { return Kappa_1st_Flow; }

  /*!
   * \brief Value of the calibrated constant for the JST method (center scheme).
   * \return Calibrated constant for the JST method for the flow equations.
   */
  su2double GetKappa_2nd_Flow(void) const { return Kappa_2nd_Flow; }

  /*!
   * \brief Value of the calibrated constant for the JST method (center scheme).
   * \return Calibrated constant for the JST method for the flow equations.
   */
  su2double GetKappa_4th_Flow(void) const { return Kappa_4th_Flow; }

  /*!
   * \brief Factor by which to multiply the dissipation contribution to Jacobians of central schemes.
   * \return The factor.
   */
  su2double GetCent_Jac_Fix_Factor(void) const { return Cent_Jac_Fix_Factor; }

  /*!
   * \brief Factor by which to multiply the dissipation contribution to Jacobians of incompressible central schemes.
   * \return The factor.
   */
  su2double GetCent_Inc_Jac_Fix_Factor(void) const { return Cent_Inc_Jac_Fix_Factor; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the adjoint flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the adjoint flow equations.
   */
  unsigned short GetKind_TimeIntScheme_AdjFlow(void) const { return Kind_TimeIntScheme_AdjFlow; }

  /*!
   * \brief Get the kind of convective numerical scheme for the adjoint flow
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the adjoint flow equations.
   */
  unsigned short GetKind_ConvNumScheme_AdjFlow(void) const { return Kind_ConvNumScheme_AdjFlow; }

  /*!
   * \brief Get the kind of center convective numerical scheme for the adjoint flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the adjoint flow equations.
   */
  CENTERED GetKind_Centered_AdjFlow(void) const { return Kind_Centered_AdjFlow; }

  /*!
   * \brief Get the kind of upwind convective numerical scheme for the adjoint flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the adjoint flow equations.
   */
  UPWIND GetKind_Upwind_AdjFlow(void) const { return Kind_Upwind_AdjFlow; }

  /*!
   * \brief Value of the calibrated constant for the high order method (center scheme).
   * \return Calibrated constant for the high order center method for the adjoint flow equations.
   */
  su2double GetKappa_2nd_AdjFlow(void) const { return Kappa_2nd_AdjFlow; }

  /*!
   * \brief Value of the calibrated constant for the high order method (center scheme).
   * \return Calibrated constant for the high order center method for the adjoint flow equations.
   */
  su2double GetKappa_4th_AdjFlow(void) const { return Kappa_4th_AdjFlow; }

  /*!
   * \brief Value of the calibrated constant for the low order method (center scheme).
   * \return Calibrated constant for the low order center method for the adjoint flow equations.
   */
  su2double GetKappa_1st_AdjFlow(void) const { return Kappa_1st_AdjFlow; }

  /*!
   * \brief Get the kind of integration scheme (implicit)
   *        for the turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the turbulence equations.
   */
  unsigned short GetKind_TimeIntScheme_Turb(void) const { return Kind_TimeIntScheme_Turb; }

  /*!
   * \brief Get the kind of convective numerical scheme for the turbulence
   *        equations (upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the turbulence equations.
   */
  unsigned short GetKind_ConvNumScheme_Turb(void) const { return Kind_ConvNumScheme_Turb; }

  /*!
   * \brief Get the kind of center convective numerical scheme for the turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the turbulence equations.
   */
  CENTERED GetKind_Centered_Turb(void) const { return Kind_Centered_Turb; }

  /*!
   * \brief Get the kind of upwind convective numerical scheme for the turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the turbulence equations.
   */
  UPWIND GetKind_Upwind_Turb(void) const { return Kind_Upwind_Turb; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the adjoint turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the adjoint turbulence equations.
   */
  unsigned short GetKind_TimeIntScheme_AdjTurb(void) const { return Kind_TimeIntScheme_AdjTurb; }

  /*!
   * \brief Get the kind of convective numerical scheme for the adjoint turbulence
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the adjoint turbulence equations.
   */
  unsigned short GetKind_ConvNumScheme_AdjTurb(void) const { return Kind_ConvNumScheme_AdjTurb; }

  /*!
   * \brief Get the kind of integration scheme (implicit)
   *        for the Species equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the Species equations.
   */
  unsigned short GetKind_TimeIntScheme_Species() const { return Kind_TimeIntScheme_Species; }

  /*!
   * \brief Get the kind of convective numerical scheme for the Species
   *        equations (upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the Species equations.
   */
  unsigned short GetKind_ConvNumScheme_Species() const { return Kind_ConvNumScheme_Species; }

  /*!
   * \brief Get the kind of center convective numerical scheme for the Species equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the Species equations.
   */
  CENTERED GetKind_Centered_Species() const { return Kind_Centered_Species; }

  /*!
   * \brief Get the kind of upwind convective numerical scheme for the Species equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the Species equations.
   */
  UPWIND GetKind_Upwind_Species() const { return Kind_Upwind_Species; }

  /*!
   * \brief Returns true if bounded scalar mode is on for species transport.
   */
  bool GetBounded_Species() const { return (Kind_Upwind_Species == UPWIND::BOUNDED_SCALAR); }

  /*!
   * \brief Returns true if bounded scalar mode is on for turbulence transport.
   */
  bool GetBounded_Turb() const { return (Kind_Upwind_Turb == UPWIND::BOUNDED_SCALAR); }

  /*!
   * \brief Returns true if bounded scalar mode is used for any equation.
   */
  bool GetBounded_Scalar() const { return GetBounded_Species() || GetBounded_Turb(); }

  /*!
   * \brief Get the kind of convective numerical scheme for the heat equation.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the heat equation.
   */
  unsigned short GetKind_ConvNumScheme_Heat(void) const { return Kind_ConvNumScheme_Heat; }

  /*!
   * \brief Get the kind of center convective numerical scheme for the adjoint turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the adjoint turbulence equations.
   */
  CENTERED GetKind_Centered_AdjTurb(void) const { return Kind_Centered_AdjTurb; }

  /*!
   * \brief Get the kind of upwind convective numerical scheme for the adjoint turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the adjoint turbulence equations.
   */
  UPWIND GetKind_Upwind_AdjTurb(void) const { return Kind_Upwind_AdjTurb; }

  /*!
   * \brief Provides information about the way in which the turbulence will be treated by the
   *        cont. adjoint method.
   * \return <code>FALSE</code> means that the adjoint turbulence equations will be used.
   */
  bool GetFrozen_Visc_Cont(void) const { return Frozen_Visc_Cont; }

  /*!
   * \brief Provides information about the way in which the turbulence will be treated by the
   *        disc. adjoint method.
   * \return <code>FALSE</code> means that the adjoint turbulence equations will be used.
   */
  bool GetFrozen_Visc_Disc(void) const { return Frozen_Visc_Disc; }

  /*!
   * \brief Provides information about using an inconsistent (primal/dual) discrete adjoint formulation
   * \return <code>FALSE</code> means that the adjoint use the same numerical methods than the primal problem.
   */
  bool GetInconsistent_Disc(void) const { return Inconsistent_Disc; }

  /*!
   * \brief Provides information about the way in which the limiter will be treated by the
   *        disc. adjoint method.
   * \return <code>FALSE</code> means that the limiter computation is included.
   */
  bool GetFrozen_Limiter_Disc(void) const { return Frozen_Limiter_Disc; }

  /*!
   * \brief Provides information about if the sharp edges are going to be removed from the sensitivity.
   * \return <code>FALSE</code> means that the sharp edges will be removed from the sensitivity.
   */
  bool GetSens_Remove_Sharp(void) const { return Sens_Remove_Sharp; }

  /*!
   * \brief Get the kind of inlet boundary condition treatment (total conditions or mass flow).
   * \return Kind of inlet boundary condition.
   */
  INLET_TYPE GetKind_Inlet(void) const { return Kind_Inlet; }

  /*!
   * \brief Check if the inlet profile(s) are specified in an input file
   * \return True if an input file is to be used for the inlet profile(s)
   */
  bool GetInlet_Profile_From_File(void) const { return Inlet_From_File; }

  /*!
   * \brief Get name of the input file for the specified inlet profile.
   * \return Name of the input file for the specified inlet profile.
   */
  string GetInlet_FileName(void) const { return Inlet_Filename; }

  /*!
   * \brief Get name of the input file for the specified actuator disk.
   * \return Name of the input file for the specified actuator disk.
   */
  string GetActDisk_FileName(void) const { return ActDisk_FileName; }

  /*!
   * \brief Get the tolerance used for matching two points on a specified inlet
   * \return Tolerance used for matching a point to a specified inlet
   */
  su2double GetInlet_Profile_Matching_Tolerance(void) const { return Inlet_Matching_Tol; }

  /*!
   * \brief Get the type of incompressible inlet from the list.
   * \return Kind of the incompressible inlet.
   */
  INLET_TYPE GetKind_Inc_Inlet(const string& val_marker) const;

  /*!
   * \brief Get the total number of types in Kind_Inc_Inlet list
   * \return Total number of types in Kind_Inc_Inlet list
   */
  unsigned short GetnInc_Inlet(void) const { return nInc_Inlet;}

  /*!
   * \brief Flag for whether the local boundary normal is used as the flow direction for an incompressible pressure inlet.
   * \return <code>FALSE</code> means the prescribed flow direction is used.
   */
  bool GetInc_Inlet_UseNormal(void) const { return Inc_Inlet_UseNormal;}

  /*!
   * \brief Get the type of incompressible outlet from the list.
   * \return Kind of the incompressible outlet.
   */
  INC_OUTLET_TYPE GetKind_Inc_Outlet(const string& val_marker) const;

  /*!
   * \brief Get the damping factor applied to velocity updates at incompressible pressure inlets.
   * \return Damping factor applied to velocity updates at incompressible pressure inlets.
   */
  su2double GetInc_Inlet_Damping(void) const { return Inc_Inlet_Damping; }

  /*!
   * \brief Get the damping factor applied to pressure updates at incompressible mass flow outlet.
   * \return Damping factor applied to pressure updates at incompressible mass flow outlet.
   */
  su2double GetInc_Outlet_Damping(void) const { return Inc_Outlet_Damping; }

  /*!
   * \brief Get the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  unsigned short GetKind_AverageProcess(void) const { return Kind_AverageProcess; }

  /*!
   * \brief Get the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  unsigned short GetKind_PerformanceAverageProcess(void) const { return Kind_PerformanceAverageProcess; }

  /*!
   * \brief Set the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  void SetKind_AverageProcess(unsigned short new_AverageProcess) { Kind_AverageProcess = new_AverageProcess; }

  /*!
   * \brief Set the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  void SetKind_PerformanceAverageProcess(unsigned short new_AverageProcess) { Kind_PerformanceAverageProcess = new_AverageProcess; }

  /*!
   * \brief Get coeff for Rotating Frame Ramp.
   * \return coeff Ramp Rotating Frame.
   */
  su2double GetRampRotatingFrame_Coeff(unsigned short iCoeff) const { return rampRotFrame_coeff[iCoeff];}

  /*!
   * \brief Get Rotating Frame Ramp option.
   * \return Ramp Rotating Frame option.
   */
  bool GetRampRotatingFrame(void) const { return RampRotatingFrame;}

  /*!
   * \brief Get coeff for Outlet Pressure Ramp.
   * \return coeff Ramp Outlet Pressure.
   */
  su2double GetRampOutletPressure_Coeff(unsigned short iCoeff) const { return rampOutPres_coeff[iCoeff];}

  /*!
   * \brief Get final Outlet Pressure value for the ramp.
   * \return final Outlet Pressure value.
   */
  su2double GetFinalOutletPressure(void) const { return  FinalOutletPressure; }

  /*!
   * \brief Get final Outlet Pressure value for the ramp.
   * \return Monitor Outlet Pressure value.
   */
  su2double GetMonitorOutletPressure(void) const { return MonitorOutletPressure; }

  /*!
   * \brief Set Monitor Outlet Pressure value for the ramp.
   */
  void SetMonitotOutletPressure(su2double newMonPres) { MonitorOutletPressure = newMonPres;}

  /*!
   * \brief Get Outlet Pressure Ramp option.
   * \return Ramp Outlet pressure option.
   */
  bool GetRampOutletPressure(void) const { return RampOutletPressure;}

  /*!
   * \brief Get mixedout coefficients.
   * \return mixedout coefficient.
   */
  su2double GetMixedout_Coeff(unsigned short iCoeff) const { return mixedout_coeff[iCoeff];}

  /*!
   * \brief Get extra relaxation factor coefficients for the Giels BC.
   * \return mixedout coefficient.
   */
  su2double GetExtraRelFacGiles(unsigned short iCoeff) const { return extrarelfac[iCoeff];}

  /*!
   * \brief Get mach limit for average massflow-based procedure .
   * \return mach limit.
   */
  su2double GetAverageMachLimit(void) const { return AverageMachLimit;}

  /*!
   * \brief Get the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  unsigned short GetKind_MixingPlaneInterface(void) const { return Kind_MixingPlaneInterface;}

  /*!
   * \brief Get the kind of turbomachinery architecture.
   * \return Kind of turbomachinery architecture.
   */
  unsigned short GetKind_TurboMachinery(unsigned short val_iZone) const { return Kind_TurboMachinery[val_iZone]; }

  /*!
   * \brief Get the kind of turbomachinery architecture.
   * \return Kind of turbomachinery architecture.
   */
  unsigned short GetKind_SpanWise(void) const { return Kind_SpanWise; }

  /*!
   * \brief Verify if there is mixing plane interface specified from config file.
   * \return boolean.
   */
  bool GetBoolMixingPlaneInterface(void) const { return (nMarker_MixingPlaneInterface !=0);}

  /*!
   * \brief Verify if there is mixing plane interface specified from config file.
   * \return boolean.
   */
  bool GetBoolTurbMixingPlane(void) const { return turbMixingPlane;}

  /*!
   * \brief Verify if there is mixing plane interface specified from config file.
   * \return boolean.
   */
  bool GetSpatialFourier(void) const { return SpatialFourier;}

  /*!
   * \brief number mixing plane interface specified from config file.
   * \return number of bound.
   */
  unsigned short GetnMarker_MixingPlaneInterface(void) const { return nMarker_MixingPlaneInterface;}

  /*!
   * \brief Verify if there is Turbomachinery performance option specified from config file.
   * \return boolean.
   */
  bool GetBoolTurbomachinery(void) const { return (nMarker_Turbomachinery !=0);}

  /*!
   * \brief number Turbomachinery blades computed using the pitch information.
   * \return nBlades.
   */
  su2double GetnBlades(unsigned short val_iZone) const { return nBlades[val_iZone];}

  /*!
   * \brief number Turbomachinery blades computed using the pitch information.
   * \return nBlades.
   */
  void SetnBlades(unsigned short val_iZone, su2double nblades) { nBlades[val_iZone] = nblades;}

  /*!
   * \brief Verify if there is any Giles Boundary Condition option specified from config file.
   * \return boolean.
   */
  bool GetBoolGiles(void) const { return (nMarker_Giles!=0);}

  /*!
   * \brief Verify if there is any Riemann Boundary Condition option specified from config file.
   * \return boolean.
   */
  bool GetBoolRiemann(void) const { return (nMarker_Riemann!=0);}

  /*!
   * \brief number Turbomachinery performance option specified from config file.
   * \return number of bound.
   */
  unsigned short GetnMarker_Turbomachinery(void) const { return nMarker_Turbomachinery;}

  /*!
   * \brief Get number of shroud markers.
   * \return number of marker shroud.
   */
  unsigned short GetnMarker_Shroud(void) const { return nMarker_Shroud;}

  /*!
   * \brief Get the marker shroud.
   * \return marker shroud.
   */
  string GetMarker_Shroud(unsigned short val_marker) const { return Marker_Shroud[val_marker];}

  /*!
   * \brief number Turbomachinery performance option specified from config file.
   * \return number of bound.
   */
  unsigned short GetnMarker_TurboPerformance(void) const { return nMarker_TurboPerformance;}

  /*!
   * \brief number span-wise sections to compute 3D BC and performance for turbomachinery specified by the user.
   * \return number of span-wise sections.
   */
  unsigned short Get_nSpanWiseSections_User(void) const { return nSpanWiseSections_User;}

  /*!
   * \brief number span-wise sections to compute 3D BC and performance for turbomachinery.
   * \return number of span-wise sections.
   */
  unsigned short GetnSpanWiseSections(void) const { return nSpanWiseSections;}

  /*!
   * \brief set number of maximum span-wise sections among all zones .
   */
  void SetnSpanMaxAllZones(unsigned short val_nSpna_max) { nSpanMaxAllZones = val_nSpna_max;}

  /*!
   * \brief number span-wise sections to compute performance for turbomachinery.
   * \return number of max span-wise sections.
   */
  unsigned short GetnSpanMaxAllZones(void) const { return nSpanMaxAllZones;}

  /*!
   * \brief set number span-wise sections to compute 3D BC and performance for turbomachinery.
   */
  void SetnSpanWiseSections(unsigned short nSpan) { nSpanWiseSections = nSpan;}

  /*!
   * \brief set number span-wise sections to compute 3D BC and performance for turbomachinery.
   */
  unsigned short GetnSpan_iZones(unsigned short iZone) const { return nSpan_iZones[iZone];}

  /*!
   * \brief set number span-wise sections to compute 3D BC and performance for turbomachinery.
   */
  void SetnSpan_iZones(unsigned short nSpan, unsigned short iZone) { nSpan_iZones[iZone] = nSpan;}

  /*!
   * \brief get inlet bounds name for Turbomachinery performance calculation.
   * \return name of the bound.
   */
  string GetMarker_TurboPerf_BoundIn(unsigned short index) const { return Marker_TurboBoundIn[index];}

  /*!
   * \brief get outlet bounds name for Turbomachinery performance calculation.
   * \return name of the bound.
   */
  string GetMarker_TurboPerf_BoundOut(unsigned short index) const { return Marker_TurboBoundOut[index];}

  /*!
   * \brief get marker kind for Turbomachinery performance calculation.
   * \return kind index.
   */
  unsigned short GetKind_TurboPerf(unsigned short index);

  /*!
   * \brief get outlet bounds name for Turbomachinery performance calculation.
   * \return name of the bound.
   */
  string GetMarker_PerBound(unsigned short val_marker) const { return Marker_PerBound[val_marker];}

  /*!
   * \brief Get the kind of inlet boundary condition treatment (total conditions or mass flow).
   * \return Kind of inlet boundary condition.
   */
  unsigned short GetKind_Engine_Inflow(void) const { return Kind_Engine_Inflow; }

  /*!
   * \brief Get the kind of inlet boundary condition treatment (total conditions or mass flow).
   * \return Kind of inlet boundary condition.
   */
  unsigned short GetKind_ActDisk(void) const { return Kind_ActDisk; }

  /*!
   * \brief Set the kind of wall - rough or smooth.
   */
  void SetKindWall(string val_marker, unsigned short val_kindwall);

  /*!
   * \brief Get the number of sections.
   * \return Number of sections
   */
  unsigned short GetnLocationStations(void) const { return nLocationStations; }

  /*!
   * \brief Get the number of sections for computing internal volume.
   * \return Number of sections for computing internal volume.
   */
  unsigned short GetnWingStations(void) const { return nWingStations; }

  /*!
   * \brief Get the location of the waterline.
   * \return Z location of the waterline.
   */
  su2double GetGeo_Waterline_Location(void) const { return Geo_Waterline_Location; }

  /*!
   * \brief Provides information about the the nodes that are going to be moved on a deformation
   *        volumetric grid deformation.
   * \return <code>TRUE</code> means that only the points on the FFD box will be moved.
   */
  bool GetHold_GridFixed(void) const { return Hold_GridFixed; }

  /*!
   * \author H. Kline
   * \brief Get the kind of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \param[in] val_obj
   * \return Kind of objective function.
   */
  unsigned short GetKind_ObjFunc(unsigned short val_obj = 0) const { return Kind_ObjFunc[val_obj]; }

  /*!
   * \author H. Kline
   * \brief Get the weight of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \return Weight of objective function.
   */
  su2double GetWeight_ObjFunc(unsigned short val_obj) const { return Weight_ObjFunc[val_obj]; }

  /*!
   * \author H. Kline
   * \brief Set the weight of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \return Weight of objective function.
   */
  void SetWeight_ObjFunc(unsigned short val_obj, su2double val) { Weight_ObjFunc[val_obj] = val; }

  /*!
   * \brief Get the user expression for the custom objective function.
   */
  const string& GetCustomObjFunc() const { return CustomObjFunc; }

  /*!
   * \brief Get the user expressions for custom outputs.
   */
  const string& GetCustomOutputs() const { return CustomOutputs; }

  /*!
   * \brief Get the kind of sensitivity smoothing technique.
   * \return Kind of sensitivity smoothing technique.
   */
  unsigned short GetKind_SensSmooth(void) const { return Kind_SensSmooth; }

  /*!
   * \brief Provides information about the time integration, and change the write in the output
   *        files information about the iteration.
   * \return The kind of time integration: Steady state, time stepping method (unsteady) or
   *         dual time stepping method (unsteady).
   */
  TIME_MARCHING GetTime_Marching() const { return TimeMarching; }

  /*!
   * \brief Provides the number of species present in the gas mixture.
   * \return The number of species present in the gas mixture.
   */
  unsigned short GetnSpecies() const { return nSpecies; }

  /*!
   * \brief Provides the gas mass fractions of the flow.
   * \return Gas Mass fractions.
   */
  const su2double *GetGas_Composition(void) const { return Gas_Composition; }

  /*!
   * \brief Provides the gas mass fractions at the wall for supercat wall.
   * \return Supercat wall gas mass fractions.
   */
  const su2double *GetSupercatalytic_Wall_Composition(void) const { return Supercatalytic_Wall_Composition; }

  /*!
   * \brief Provides the restart information.
   * \return Restart information, if <code>TRUE</code> then the code will use the solution as restart.
   */
  bool GetRestart(void) const { return Restart; }

  /*!
   * \brief Flag for whether binary SU2 native restart files are read.
   * \return Flag for whether binary SU2 native restart files are read, if <code>TRUE</code> then the code will load binary restart files.
   */
  bool GetRead_Binary_Restart(void) const { return Read_Binary_Restart; }

  /*!
   * \brief Flag for whether restart solution files are overwritten.
   * \return Flag for overwriting. If Flag=false, iteration nr is appended to filename
   */
  bool GetWrt_Restart_Overwrite(void) const { return Wrt_Restart_Overwrite; }

    /*!
   * \brief Flag for whether visualization files are overwritten.
   * \return Flag for overwriting. If Flag=false, iteration nr is appended to filename
   */
  bool GetWrt_Surface_Overwrite(void) const { return Wrt_Surface_Overwrite; }

   /*!
   * \brief Flag for whether visualization files are overwritten.
   * \return Flag for overwriting. If Flag=false, iteration nr is appended to filename
   */
  bool GetWrt_Volume_Overwrite(void) const { return Wrt_Volume_Overwrite; }

  /*!
   * \brief Provides the number of varaibles.
   * \return Number of variables.
   */
  unsigned short GetnVar(void);

  /*!
   * \brief Provides the number of varaibles.
   * \return Number of variables.
   */
  unsigned short GetnZone(void) const { return nZone; }

  /*!
   * \brief Provides the number of varaibles.
   * \return Number of variables.
   */
  unsigned short GetiZone(void) const { return iZone; }

  /*!
   * \brief For some problems like adjoint or the linearized equations it
   *          is necessary to restart the flow solution.
   * \return Flow restart information, if <code>TRUE</code> then the code will restart the flow solution.
   */

  bool GetRestart_Flow(void) const { return Restart_Flow; }

  /*!
   * \brief Indicates whether the flow is frozen (chemistry deactivated).
   */
  bool GetFrozen(void) const { return frozen; }

  /*!
   * \brief Indicates whether electron gas is present in the gas mixture.
   */
  bool GetIonization(void) const { return ionization; }

  /*!
   * \brief Indicates whether the VT source residual is limited.
   */
  bool GetVTTransferResidualLimiting(void) const { return vt_transfer_res_limit; }

  /*!
   * \brief Indicates if mixture is monoatomic.
   */
  bool GetMonoatomic(void) const { return monoatomic; }

  /*!
   * \brief Indicates whether supercatalytic wall is used.
   */
  bool GetSupercatalytic_Wall(void) const { return Supercatalytic_Wall; }

  /*!
   * \brief Information about computing and plotting the equivalent area distribution.
   * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are computing the equivalent area.
   */
  bool GetEquivArea(void) const { return EquivArea; }

  /*!
   * \brief Information about computing and plotting the equivalent area distribution.
   * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are computing the equivalent area.
   */
  bool GetInvDesign_Cp(void) const { return InvDesign_Cp; }

  /*!
   * \brief Information about computing and plotting the equivalent area distribution.
   * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are computing the equivalent area.
   */
  bool GetInvDesign_HeatFlux(void) const { return InvDesign_HeatFlux; }

  /*!
   * \brief Get name of the input grid.
   * \return File name of the input grid.
   */
  string GetMesh_FileName(void) const { return Mesh_FileName; }

  /*!
   * \brief Get name of the output grid, this parameter is important for grid
   *        adaptation and deformation.
   * \return File name of the output grid.
   */
  string GetMesh_Out_FileName(void) const { return Mesh_Out_FileName; }

  /*!
   * \brief Get the name of the file with the solution of the flow problem.
   * \return Name of the file with the solution of the flow problem.
   */
  string GetSolution_FileName(void) const { return Solution_FileName; }

  /*!
   * \brief Get the name of the file with the solution of the adjoint flow problem
   *          with drag objective function.
   * \return Name of the file with the solution of the adjoint flow problem with
   *         drag objective function.
   */
  string GetSolution_AdjFileName(void) const { return Solution_AdjFileName; }

  /*!
   * \brief Get the format of the input/output grid.
   * \return Format of the input/output grid.
   */
  unsigned short GetMesh_FileFormat(void) const { return Mesh_FileFormat; }

  /*!
   * \brief Get the format of the output solution.
   * \return Format of the output solution.
   */
  TAB_OUTPUT GetTabular_FileFormat(void) const { return Tab_FileFormat; }

  /*!
   * \brief Get the output precision to be used in <ofstream>.precision(value) for history and SU2_DOT output.
   * \return Output precision.
   */
  unsigned short GetOutput_Precision(void) const { return output_precision; }

  /*!
   * \brief Get the format of the output solution.
   * \return Format of the output solution.
   */
  unsigned short GetActDisk_Jump(void) const { return ActDisk_Jump; }

  /*!
   * \brief Get the name of the file with the convergence history of the problem.
   * \return Name of the file with convergence history of the problem.
   */
  string GetConv_FileName(void) const { return Conv_FileName; }

  /*!
   * \brief Get the Starting Iteration for the windowing approach
   *        in Sensitivity Analysis for period-averaged outputs, which oscillate.
   * \return
   */
  unsigned long GetStartWindowIteration(void) const { return StartWindowIteration; }

  /*!
   * \brief Get Index of the window function used as weight in the cost functional
   * \return
   */
  WINDOW_FUNCTION GetKindWindow(void) const { return Kind_WindowFct; }

  /*!
   * \brief Get the name of the file with the forces breakdown of the problem.
   * \return Name of the file with forces breakdown of the problem.
   */
  string GetBreakdown_FileName(void) const { return Breakdown_FileName; }

  /*!
   * \brief Get the name of the file with the flow variables.
   * \return Name of the file with the primitive variables.
   */
  string GetVolume_FileName(void) const { return Volume_FileName; }

  /*!
   * \brief Add any numbers necessary to the filename (iteration number, zone ID ...)
   * \param[in] filename - the base filename.
   * \param[in] ext - the extension to be added.
   * \param[in] Iter - the current iteration
   * \return The new filename
   */
  string GetFilename(string filename, const string& ext, int Iter) const;

  /*!
   * \brief Add steady iteration number to the filename (does not overwrite previous files)
   * \param[in] filename - the base filename.
   * \param[in] inner_iter - the inner iterations
   * \param[in] outer_iter - the outer iterations
   * \return The new filename
   */
  string GetFilename_Iter(const string& filename_iter, unsigned long curInnerIter, unsigned long curOuterIter) const;

  /*!
   * \brief Append the zone index to the restart or the solution files.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultizone_FileName(string val_filename, int val_iZone, const string& ext) const;

  /*!
   * \brief Append the zone index to the restart or the solution files.
   * \param[in] val_filename - the base filename.
   * \param[in] val_iZone - the zone ID.
   * \param[in] ext - the filename extension.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultizone_HistoryFileName(string val_filename, int val_iZone, const string& ext) const;

  /*!
   * \brief Append the instance index to the restart or the solution files.
   * \param[in] val_filename - the base filename.
   * \param[in] val_iInst - the current instance.
   * \param[in] ext - the filename extension.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultiInstance_FileName(string val_filename, int val_iInst, const string& ext) const;

  /*!
   * \brief Append the instance index to the restart or the solution files.
   * \param[in] val_filename - the base filename.
   * \param[in] val_iInst - the current instance.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultiInstance_HistoryFileName(string val_filename, int val_iInst) const;

  /*!
   * \brief Get the name of the restart file for the flow variables.
   * \return Name of the restart file for the flow variables.
   */
  string GetRestart_FileName(void) const { return Restart_FileName; }

  /*!
   * \brief Get the name of the restart file for the adjoint variables (drag objective function).
   * \return Name of the restart file for the adjoint variables (drag objective function).
   */
  string GetRestart_AdjFileName(void) const { return Restart_AdjFileName; }

  /*!
   * \brief Get the name of the file with the adjoint variables.
   * \return Name of the file with the adjoint variables.
   */
  string GetAdj_FileName(void) const { return Adj_FileName; }

  /*!
   * \brief Get the name of the file with the gradient of the objective function.
   * \return Name of the file with the gradient of the objective function.
   */
  string GetObjFunc_Grad_FileName(void) const { return ObjFunc_Grad_FileName; }

  /*!
   * \brief Get the name of the file with the gradient of the objective function.
   * \return Name of the file with the gradient of the objective function.
   */
  string GetObjFunc_Value_FileName(void) const { return ObjFunc_Value_FileName; }

  /*!
   * \brief Get the name of the file with the surface information for the flow problem.
   * \return Name of the file with the surface information for the flow problem.
   */
  string GetSurfCoeff_FileName(void) const { return SurfCoeff_FileName; }

  /*!
   * \brief Get the name of the file with the surface information for the adjoint problem.
   * \return Name of the file with the surface information for the adjoint problem.
   */
  string GetSurfAdjCoeff_FileName(void) const { return SurfAdjCoeff_FileName; }

  /*!
   * \brief Get the name of the file with the surface sensitivity (discrete adjoint).
   * \return Name of the file with the surface sensitivity (discrete adjoint).
   */
  string GetSurfSens_FileName(void) const { return SurfSens_FileName; }

  /*!
   * \brief Get the name of the file with the volume sensitivity (discrete adjoint).
   * \return Name of the file with the volume sensitivity (discrete adjoint).
   */
  string GetVolSens_FileName(void) const { return VolSens_FileName; }

  /*!
   * \brief Augment the input filename with the iteration number for an unsteady file.
   * \param[in] val_filename - String value of the base filename.
   * \param[in] val_iter - Unsteady iteration number or time instance.
   * \param[in] ext - the filename extension.
   * \return Name of the file with the iteration number for an unsteady solution file.
   */
  string GetUnsteady_FileName(string val_filename, int val_iter, const string& ext) const;

  /*!
   * \brief Append the input filename string with the appropriate objective function extension.
   * \param[in] val_filename - String value of the base filename.
   * \return Name of the file with the appropriate objective function extension.
   */
  string GetObjFunc_Extension(string val_filename) const;

  /*!
   * \brief Get functional that is going to be used to evaluate the residual flow convergence.
   * \return Functional that is going to be used to evaluate the residual flow convergence.
   */
  unsigned short GetResidual_Func_Flow(void) const { return Residual_Func_Flow; }

  /*!
   * \brief Get functional that is going to be used to evaluate the flow convergence.
   * \return Functional that is going to be used to evaluate the flow convergence.
   */
  unsigned short GetCauchy_Func_Flow(void) const { return Cauchy_Func_Flow; }

  /*!
   * \brief Get functional that is going to be used to evaluate the adjoint flow convergence.
   * \return Functional that is going to be used to evaluate the adjoint flow convergence.
   */
  unsigned short GetCauchy_Func_AdjFlow(void) const { return Cauchy_Func_AdjFlow; }

  /*!
   * \brief Get the number of iterations that are considered in the Cauchy convergence criteria.
   * \return Number of elements in the Cauchy criteria.
   */
  unsigned short GetCauchy_Elems(void) const { return Cauchy_Elems; }

  /*!
   * \brief Get the number of iterations that are not considered in the convergence criteria.
   * \return Number of iterations before starting with the convergence criteria.
   */
  unsigned long GetStartConv_Iter(void) const { return StartConv_Iter; }

  /*!
   * \brief Get the value of convergence criteria for the Cauchy method in the direct,
   *        adjoint or linearized problem.
   * \return Value of the convergence criteria.
   */
  su2double GetCauchy_Eps(void) const { return Cauchy_Eps; }

  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetDelta_UnstTimeND(void) const { return Delta_UnstTimeND; }

  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetTotal_UnstTimeND(void) const { return Total_UnstTimeND; }

  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation.
   */
  su2double GetDelta_UnstTime(void) const { return Delta_UnstTime; }

  /*!
   * \brief Set the value of the unsteadty time step using the CFL number.
   * \param[in] val_delta_unsttimend - Value of the unsteady time step using CFL number.
   */
  void SetDelta_UnstTimeND(su2double val_delta_unsttimend) { Delta_UnstTimeND = val_delta_unsttimend; }

  /*!
   * \brief If we are performing an unsteady simulation, this is the
   *    value of max physical time for which we run the simulation
   * \return Value of the physical time in an unsteady simulation.
   */
  su2double GetTotal_UnstTime(void) const { return Total_UnstTime; }

  /*!
   * \brief If we are performing an unsteady simulation, this is the
   *    value of current time.
   * \return Value of the physical time in an unsteady simulation.
   */
  su2double GetCurrent_UnstTime(void) const { return Current_UnstTime; }

  /*!
   * \brief Divide the rectbles and hexahedron.
   * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
   */
  bool GetSubsonicEngine(void) const { return SubsonicEngine; }

  /*!
   * \brief Actuator disk defined with a double surface.
   * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
   */
  bool GetActDisk_DoubleSurface(void) const { return ActDisk_DoubleSurface; }

  /*!
   * \brief Only halg of the engine is in the compputational grid.
   * \return <code>TRUE</code> if the engine is complete; otherwise <code>FALSE</code>.
   */
  bool GetEngine_HalfModel(void) const { return Engine_HalfModel; }

  /*!
   * \brief Actuator disk defined with a double surface.
   * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
   */
  bool GetActDisk_SU2_DEF(void) const { return ActDisk_SU2_DEF; }

  /*!
   * \brief Value of the design variable step, we use this value in design problems.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \param[in] val_val - Value of the design variable that we want to read.
   * \return Design variable step.
   */
  su2double& GetDV_Value(unsigned short val_dv, unsigned short val_val = 0) { return DV_Value[val_dv][val_val]; }
  const su2double& GetDV_Value(unsigned short val_dv, unsigned short val_val = 0) const { return DV_Value[val_dv][val_val]; }

  /*!
   * \brief Set the value of the design variable step, we use this value in design problems.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \param[in] val_ind - value of initial deformation.
   * \param[in] val    - Value of the design variable.
   */
  void SetDV_Value(unsigned short val_dv, unsigned short val_ind, su2double val) { DV_Value[val_dv][val_ind] = val; }

  /*!
   * \brief Get information about the grid movement.
   * \return <code>TRUE</code> if there is a grid movement; otherwise <code>FALSE</code>.
   */
  bool GetGrid_Movement(void) const {
    return (Kind_GridMovement != NO_MOVEMENT) || (nKind_SurfaceMovement > 0);
  }

  /*!
   * \brief Get information about dynamic grids.
   * \return <code>TRUE</code> if there is a grid movement; otherwise <code>FALSE</code>.
   */
  bool GetDynamic_Grid(void) const { return GetGrid_Movement() || (Deform_Mesh && Time_Domain); }

  /*!
   * \brief Get information about the volumetric movement.
   * \return <code>TRUE</code> if there is a volumetric movement is required; otherwise <code>FALSE</code>.
   */
  bool GetVolumetric_Movement(void) const;

  /*!
   * \brief Get information about deforming markers.
   * \param[in] kind_movement - Kind of surface movement.
   * \return <code>TRUE</code> at least one surface of kind_movement moving; otherwise <code>FALSE</code>.
   */
  bool GetSurface_Movement(unsigned short kind_movement) const;

  /*!
   * \brief Set a surface movement marker.
   * \param[in] iMarker - Moving marker.
   * \param[in] kind_movement - Kind of surface movement.
   * \return <code>TRUE</code> at least one surface of kind_movement moving; otherwise <code>FALSE</code>.
   */
  void SetSurface_Movement(unsigned short iMarker, unsigned short kind_movement);

  /*!
   * \brief Get the type of dynamic mesh motion. Each zone gets a config file.
   * \return Type of dynamic mesh motion.
   */
  unsigned short GetKind_GridMovement() const { return Kind_GridMovement; }

  /*!
   * \brief Set the type of dynamic mesh motion.
   * \param[in] motion_Type - Specify motion type.
   */
  void SetKind_GridMovement(unsigned short motion_Type) { Kind_GridMovement = motion_Type; }

  /*!
   * \brief Get the type of surface motion.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving).
   * \return Type of surface motion.
   */
  unsigned short GetKind_SurfaceMovement(unsigned short iMarkerMoving) const { return Kind_SurfaceMovement[iMarkerMoving];}

  /*!
   * \brief Get the mach number based on the mesh velocity and freestream quantities.
   * \return Mach number based on the mesh velocity and freestream quantities.
   */
  su2double GetMach_Motion(void) const { return Mach_Motion; }

  /*!
   * \brief Get the mesh motion origin.
   * \param[in] iDim - spatial component
   * \return The mesh motion origin.
   */
  su2double GetMotion_Origin(unsigned short iDim) const { return Motion_Origin[iDim];}

  /*!
   * \brief Set the mesh motion origin.
   * \param[in] val - new value of the origin
   * \return The mesh motion origin.
   */
  void SetMotion_Origin(const su2double* val) { for (int iDim = 0; iDim < 3; iDim++) Motion_Origin[iDim] = val[iDim]; }

  /*!
   * \brief Get the mesh motion origin.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return The motion origin of the marker.
   */
  su2double GetMarkerMotion_Origin(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerMotion_Origin[3*iMarkerMoving + iDim];}

  /*!
   * \brief Set the mesh motion origin.
   * \param[in] val - new value of the origin
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   */
  void SetMarkerMotion_Origin(const su2double* val, unsigned short iMarkerMoving) {
    for (int iDim = 0; iDim < 3; iDim++) MarkerMotion_Origin[3*iMarkerMoving + iDim] = val[iDim];
  }

  /*!
   * \brief Get the translational velocity of the mesh.
   * \param[in] iDim - spatial component
   * \return Translational velocity of the mesh.
   */
  su2double GetTranslation_Rate(unsigned short iDim) const { return Translation_Rate[iDim];}

  /*!
   * \brief Set the translational velocity of the mesh.
   * \param[in] iDim - spatial component
   * \return Translational velocity of the mesh.
   */
  void SetTranslation_Rate(unsigned short iDim, su2double val) { Translation_Rate[iDim] = val;}

  /*!
   * \brief Get the translational velocity of the marker.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return Translational velocity of the marker.
   */
  su2double GetMarkerTranslationRate(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerTranslation_Rate[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the rotation rate of the mesh.
   * \param[in] iDim - spatial component
   * \return Translational velocity of the mesh.
   */
  su2double GetRotation_Rate(unsigned short iDim) const { return Rotation_Rate[iDim];}

  /*!
   * \brief Get the rotation rate of the mesh.
   * \param[in] iDim - spatial component
   * \param[in] val - new value of the rotation rate.
   * \return Translational velocity of the mesh.
   */
  void SetRotation_Rate(unsigned short iDim, su2double val) { Rotation_Rate[iDim] = val;}

  /*!
   * \brief Get the rotation rate of the marker.
   *  \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return Rotation velocity of the marker.
   */
  su2double GetMarkerRotationRate(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerRotation_Rate[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the pitching rate of the mesh.
   * \param[in] iDim - spatial component
   * \return Angular frequency of the mesh pitching.
   */
  su2double GetPitching_Omega(unsigned short iDim) const { return Pitching_Omega[iDim];}

  /*!
   * \brief Get pitching rate of the marker.
   * \param[in] iMarkerMoving - Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return  Angular frequency of the marker pitching.
   */
  su2double GetMarkerPitching_Omega(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerPitching_Omega[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the pitching amplitude of the mesh.
   * \param[in] iDim - spatial component
   * \return pitching amplitude of the mesh.
   */
  su2double GetPitching_Ampl(unsigned short iDim) const { return Pitching_Ampl[iDim];}

  /*!
   * \brief Get pitching amplitude of the marker.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return  pitching amplitude of the marker.
   */
  su2double GetMarkerPitching_Ampl(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerPitching_Ampl[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the pitching phase of the mesh.
   * \param[in] iDim - spatial component.
   * \return pitching phase of the mesh.
   */
  su2double GetPitching_Phase(unsigned short iDim) const { return Pitching_Phase[iDim];}

  /*!
   * \brief Get pitching phase of the marker.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving) \
   * \param[in] iDim - spatial component
   * \return pitching phase of the marker.
   */
  su2double GetMarkerPitching_Phase(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerPitching_Phase[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the plunging rate of the mesh.
   * \param[in] iDim - spatial component
   * \return Angular frequency of the mesh plunging.
   */
  su2double GetPlunging_Omega(unsigned short iDim) const { return Plunging_Omega[iDim];}

  /*!
   * \brief Get plunging rate of the marker.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return Angular frequency of the marker plunging.
   */
  su2double GetMarkerPlunging_Omega(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerPlunging_Omega[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the plunging amplitude of the mesh.
   * \param[in] iDim - spatial component
   * \return Plunging amplitude of the mesh.
   */
  su2double GetPlunging_Ampl(unsigned short iDim) const { return Plunging_Ampl[iDim];}

  /*!
   * \brief Get plunging amplitude of the marker.
   * \param[in] iMarkerMoving -  Index of the moving marker (as specified in Marker_Moving)
   * \param[in] iDim - spatial component
   * \return Plunging amplitude of the marker.
   */
  su2double GetMarkerPlunging_Ampl(unsigned short iMarkerMoving, unsigned short iDim) const { return MarkerPlunging_Ampl[3*iMarkerMoving + iDim];}

  /*!
   * \brief Get the angular velocity of the mesh about the z-axis.
   * \return Angular velocity of the mesh about the z-axis.
   */
  su2double GetFinalRotation_Rate_Z() const { return FinalRotation_Rate_Z;}

  /*!
   * \brief Set the angular velocity of the mesh about the z-axis.
   * \param[in] newRotation_Rate_Z - new rotation rate after computing the ramp value.
   */
  void SetRotation_Rate_Z(su2double newRotation_Rate_Z);

  /*!
   * \brief Get the Harmonic Balance frequency pointer.
   * \return Harmonic Balance Frequency pointer.
   */
  const su2double* GetOmega_HB(void) const { return  Omega_HB; }

  /*!
   * \brief Get if harmonic balance source term is to be preconditioned
   * \return yes or no to harmonic balance preconditioning
   */
  bool GetHB_Precondition(void) const { return HB_Precondition; }

  /*!
   * \brief Get if we should update the motion origin.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return yes or no to update motion origin.
   */
  unsigned short GetMoveMotion_Origin(unsigned short val_marker) const { return MoveMotion_Origin[val_marker]; }

  /*!
   * \brief Get the minimum value of Beta for Roe-Turkel preconditioner
   * \return the minimum value of Beta for Roe-Turkel preconditioner
   */
  su2double GetminTurkelBeta() const { return  Min_Beta_RoeTurkel; }

  /*!
   * \brief Get the minimum value of Beta for Roe-Turkel preconditioner
   * \return the minimum value of Beta for Roe-Turkel preconditioner
   */
  su2double GetmaxTurkelBeta() const { return  Max_Beta_RoeTurkel; }

  /*!
   * \brief Get information about the adibatic wall condition
   * \return <code>TRUE</code> if it is a adiabatic wall condition; otherwise <code>FALSE</code>.
   */
  bool GetAdiabaticWall(void);

  /*!
   * \brief Get information about the isothermal wall condition
   * \return <code>TRUE</code> if it is a isothermal wall condition; otherwise <code>FALSE</code>.
   */
  bool GetIsothermalWall(void);

  /*!
   * \brief Get information about the Low Mach Preconditioning
   * \return <code>TRUE</code> if we are using low Mach preconditioner; otherwise <code>FALSE</code>.
   */
  bool Low_Mach_Preconditioning(void) const { return Low_Mach_Precon; }

  /*!
   * \brief Get information about the Low Mach Correction
   * \return <code>TRUE</code> if we are using low Mach correction; otherwise <code>FALSE</code>.
   */
  bool Low_Mach_Correction(void) const { return Low_Mach_Corr; }

  /*!
   * \brief Get information about the poisson solver condition
   * \return <code>TRUE</code> if it is a poisson solver condition; otherwise <code>FALSE</code>.
   */
  bool GetPoissonSolver(void) const { return PoissonSolver; }

  /*!
   * \brief Get information about the gravity force.
   * \return <code>TRUE</code> if it uses the gravity force; otherwise <code>FALSE</code>.
   */
  bool GetGravityForce(void) const { return GravityForce; }

  /*!
   * \brief Get information about the Vorticity Confinement.
   * \return <code>TRUE</code> if it uses Vorticity Confinement; otherwise <code>FALSE</code>.
   */
  bool GetVorticityConfinement(void) const { return VorticityConfinement; }

  /*!
   * \brief Get information about the body force.
   * \return <code>TRUE</code> if it uses a body force; otherwise <code>FALSE</code>.
   */
  bool GetBody_Force(void) const { return Body_Force; }

  /*!
   * \brief Get a pointer to the body force vector.
   * \return A pointer to the body force vector.
   */
  const su2double* GetBody_Force_Vector(void) const { return body_force; }

  /*!
   * \brief Get information about the streamwise periodicity (None, Pressure_Drop, Massflow).
   * \return Driving force identification.
   */
  ENUM_STREAMWISE_PERIODIC GetKind_Streamwise_Periodic(void) const { return Kind_Streamwise_Periodic; }

  /*!
   * \brief Get information about the streamwise periodicity Energy equation handling.
   * \return Real periodic treatment of energy equation.
   */
  bool GetStreamwise_Periodic_Temperature(void) const { return Streamwise_Periodic_Temperature; }

  /*!
   * \brief Get the value of the artificial periodic outlet heat.
   * \return Heat value.
   */
  su2double GetStreamwise_Periodic_OutletHeat(void) const { return Streamwise_Periodic_OutletHeat; }

  /*!
   * \brief Get the value of the pressure delta from which body force vector is computed.
   * \return Delta Pressure for body force computation.
   */
  su2double GetStreamwise_Periodic_PressureDrop(void) const { return Streamwise_Periodic_PressureDrop; }

  /*!
   * \brief Set the value of the pressure delta from which body force vector is computed. Necessary for Restart metadata.
   */
  void SetStreamwise_Periodic_PressureDrop(su2double Streamwise_Periodic_PressureDrop_) { Streamwise_Periodic_PressureDrop = Streamwise_Periodic_PressureDrop_; }

  /*!
   * \brief Get the value of the massflow from which body force vector is computed.
   * \return Massflow for body force computation.
   */
  su2double GetStreamwise_Periodic_TargetMassFlow(void) const { return Streamwise_Periodic_TargetMassFlow; }

  /*!
   * \brief Get information about the volumetric heat source.
   * \return <code>TRUE</code> if it uses a volumetric heat source; otherwise <code>FALSE</code>.
   */
  inline bool GetHeatSource(void) const { return HeatSource; }

  /*!
   * \brief Get information about the volumetric heat source.
   * \return Value of the volumetric heat source
   */
  inline su2double GetHeatSource_Val(void) const {return ValHeatSource;}

  /*!
   * \brief Get the rotation angle of the volumetric heat source in axis Z.
   * \return Rotation (Z) of the volumetric heat source
   */
  inline su2double GetHeatSource_Rot_Z(void) const {return Heat_Source_Rot_Z;}

  /*!
   * \brief Set the rotation angle of the volumetric heat source in axis Z.
   * \param[in] val_rot - Rotation (Z) of the volumetric heat source
   */
  inline void SetHeatSource_Rot_Z(su2double val_rot) {Heat_Source_Rot_Z = val_rot;}

  /*!
   * \brief Get the position of the center of the volumetric heat source.
   * \return Pointer to the center of the ellipsoid that introduces a volumetric heat source.
   */
  inline const su2double* GetHeatSource_Center(void) const {return hs_center;}

  /*!
   * \brief Set the position of the center of the volumetric heat source.
   * \param[in] x_cent = X position of the center of the volumetric heat source.
   * \param[in] y_cent = Y position of the center of the volumetric heat source.
   * \param[in] z_cent = Z position of the center of the volumetric heat source.
   */
  inline void SetHeatSource_Center(su2double x_cent, su2double y_cent, su2double z_cent) {
    hs_center[0] = x_cent; hs_center[1] = y_cent; hs_center[2] = z_cent;
  }

  /*!
   * \brief Get the radius of the ellipsoid that introduces a volumetric heat source.
   * \return Pointer to the radii (x, y, z) of the ellipsoid that introduces a volumetric heat source.
   */
  inline const su2double* GetHeatSource_Axes(void) const {return hs_axes;}

  /*!
   * \brief Get information about the rotational frame.
   * \return <code>TRUE</code> if there is a rotational frame; otherwise <code>FALSE</code>.
   */
  bool GetRotating_Frame(void) const { return Rotating_Frame; }

  /*!
   * \brief Get information about the axisymmetric frame.
   * \return <code>TRUE</code> if there is a rotational frame; otherwise <code>FALSE</code>.
   */
  bool GetAxisymmetric(void) const { return Axisymmetric; }

  /*!
   * \brief Get information about there is a smoothing of the grid coordinates.
   * \return <code>TRUE</code> if there is smoothing of the grid coordinates; otherwise <code>FALSE</code>.
   */
  unsigned short GetSmoothNumGrid(void) const { return SmoothNumGrid; }

  /*!
   * \brief Subtract one to the index of the finest grid (full multigrid strategy).
   * \return Change the index of the finest grid.
   */
  void SubtractFinestMesh(void) { FinestMesh = FinestMesh-1; }

  /*!
   * \brief Obtain the kind of design variable.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \return Design variable identification.
   */
  unsigned short GetDesign_Variable(unsigned short val_dv) const { return Design_Variable[val_dv]; }

  /*!
   * \brief Get the buffet sensor sharpness coefficient.
   * \return Sharpness coefficient for buffet sensor.
   */
  su2double GetBuffet_k(void) const { return Buffet_k; }

  /*!
   * \brief Get the buffet sensor offset parameter.
   * \return Offset parameter for buffet sensor.
   */
  su2double GetBuffet_lambda(void) const { return Buffet_lambda; }

  /*!
   * \brief Get the index in the config information of the marker <i>val_marker</i>.
   * \note When we read the config file, it stores the markers in a particular vector.
   * \return Index in the config information of the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_TagBound(const string& val_marker) const;

  /*!
   * \brief Get the name in the config information of the marker number <i>val_marker</i>.
   * \note When we read the config file, it stores the markers in a particular vector.
   * \return Name of the marker in the config information of the marker <i>val_marker</i>.
   */
  string GetMarker_CfgFile_TagBound(unsigned short val_marker) const;

  /*!
   * \brief Get the boundary information (kind of boundary) in the config information of the marker <i>val_marker</i>.
   * \return Kind of boundary in the config information of the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_KindBC(const string& val_marker) const;

  /*!
   * \brief Get the monitoring information from the config definition for the marker <i>val_marker</i>.
   * \return Monitoring information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Monitoring(const string& val_marker) const;

  /*!
   * \brief Get the monitoring information from the config definition for the marker <i>val_marker</i>.
   * \return Monitoring information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_GeoEval(const string& val_marker) const;

  /*!
   * \brief Get the monitoring information from the config definition for the marker <i>val_marker</i>.
   * \return Monitoring information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Designing(const string& val_marker) const;

  /*!
   * \brief Get the plotting information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Plotting(const string& val_marker) const;

  /*!
   * \brief Get the plotting information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Analyze(const string& val_marker) const;

  /*!
   * \brief Get the multi-physics interface information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_ZoneInterface(const string& val_marker) const;

  /*!
   * \brief Get the TurboPerformance information from the config definition for the marker <i>val_marker</i>.
   * \return TurboPerformance information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Turbomachinery(const string& val_marker) const;

  /*!
   * \brief Get the TurboPerformance flag information from the config definition for the marker <i>val_marker</i>.
   * \return TurboPerformance flag information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_TurbomachineryFlag(const string& val_marker) const;

  /*!
   * \brief Get the MixingPlane interface information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_MixingPlaneInterface(const string& val_marker) const;

  /*!
   * \brief Get the DV information from the config definition for the marker <i>val_marker</i>.
   * \return DV information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_DV(const string& val_marker) const;

  /*!
   * \brief Get the motion information from the config definition for the marker <i>val_marker</i>.
   * \return Motion information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Moving(const string& val_marker) const;

  /*!
   * \brief Get the gradient boundary information from the config definition for the marker <i>val_marker</i>.
   * \return Gradient boundary information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_SobolevBC(const string& val_marker) const;

  /*!
   * \brief Get the DEFORM_MESH information from the config definition for the marker <i>val_marker</i>.
   * \return DEFORM_MESH information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Deform_Mesh(const string& val_marker) const;

  /*!
   * \brief Get the DEFORM_MESH_SYM_PLANE information from the config definition for the marker <i>val_marker</i>.
   * \return DEFORM_MESH_SYM_PLANE information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Deform_Mesh_Sym_Plane(const string& val_marker) const;

  /*!
   * \brief Get the Fluid_Load information from the config definition for the marker <i>val_marker</i>.
   * \return Fluid_Load information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Fluid_Load(const string& val_marker) const;

  /*!
   * \brief Get the Python customization information from the config definition for the marker <i>val_marker</i>.
   * \return Python customization information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_PyCustom(const string& val_marker) const;

  /*!
   * \brief Get the periodic information from the config definition of the marker <i>val_marker</i>.
   * \return Periodic information of the boundary in the config information of the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_PerBound(const string& val_marker) const;

  /*!
   * \brief  Get the name of the marker <i>val_marker</i>.
   * \return The interface which owns that marker <i>val_marker</i>.
   */
  unsigned short GetMarker_ZoneInterface(const string& val_marker) const;

  /*!
   * \brief  Get the name of the marker <i>val_iMarker</i>.
   * \return The name of the marker in the interface
   */
  string GetMarkerTag_ZoneInterface(unsigned short val_iMarker) const { return Marker_ZoneInterface[val_iMarker]; }

  /*!
   * \brief  Get the number of markers in the multizone interface.
   * \return The number markers in the multizone interface
   */
  unsigned short GetnMarker_ZoneInterface(void) const { return nMarker_ZoneInterface; }

  /*!
   * \brief Determines whether a marker with index iMarker is a solid boundary.
   * \param iMarker
   * \return <TRUE> it marker with index iMarker is a solid boundary.
   */
  bool GetSolid_Wall(unsigned short iMarker) const;

  /*!
   * \brief Determines whether a marker with index iMarker is a viscous no-slip boundary.
   * \param iMarker
   * \return <TRUE> it marker with index iMarker is a viscous no-slip boundary.
   */
  bool GetViscous_Wall(unsigned short iMarker) const;

  /*!
   * \brief Determines whether a marker with index iMarker is a catalytic boundary.
   * \param iMarker
   * \return <TRUE> it marker with index iMarker is a catalytic boundary.
   */
  bool GetCatalytic_Wall(unsigned short iMarker) const;

  /*!
   * \brief Determines if problem is adjoint.
   * \return true if Adjoint.
   */
  bool GetContinuous_Adjoint(void) const { return ContinuousAdjoint; }

  /*!
   * \brief Determines if problem is viscous.
   * \return true if Viscous.
   */
  bool GetViscous(void) const { return Viscous; }

  /*!
   * \brief Determines if problem has catalytic walls.
   * \return true if catalytic walls are present.
   */
  bool GetCatalytic(void) const { return nWall_Catalytic > 0; }

  /*!
   * \brief Provides the index of the solution in the container.
   * \param[in] val_eqsystem - Equation that is being solved.
   * \return Index on the solution container.
   */
  unsigned short GetContainerPosition(unsigned short val_eqsystem);

  /*!
   * \brief Value of the minimum residual value (log10 scale).
   * \return Value of the minimum residual value (log10 scale).
   */
  su2double GetMinLogResidual(void) const { return MinLogResidual; }

  /*!
   * \brief Value of the damping factor for the engine inlet bc.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Engine_Inflow(void) const { return Damp_Engine_Inflow; }

  /*!
   * \brief Value of the damping factor for the engine exhaust inlet bc.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Engine_Exhaust(void) const { return Damp_Engine_Exhaust; }

  /*!
   * \brief Value of the damping factor for the residual restriction.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Res_Restric(void) const { return Damp_Res_Restric; }

  /*!
   * \brief Value of the damping factor for the correction prolongation.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Correc_Prolong(void) const { return Damp_Correc_Prolong; }

  /*!
   * \brief Value of the position of the Near Field (y coordinate for 2D, and z coordinate for 3D).
   * \return Value of the Near Field position.
   */
  su2double GetPosition_Plane(void) const { return Position_Plane; }

  /*!
   * \brief Value of the weight of the drag coefficient in the Sonic Boom optimization.
   * \return Value of the weight of the drag coefficient in the Sonic Boom optimization.
   */
  su2double GetWeightCd(void) const { return WeightCd; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdNetThrust_dBCThrust(su2double val_dnetthrust_dbcthrust);

  /*!
   * \brief Value of the azimuthal line to fix due to a misalignments of the nearfield.
   * \return Azimuthal line to fix due to a misalignments of the nearfield.
   */
  su2double GetFixAzimuthalLine(void) const { return FixAzimuthalLine; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCD_dCL(void) const { return dCD_dCL; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCD_dCL(su2double val_dcd_dcl) { dCD_dCL = val_dcd_dcl; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCMx_dCL(void) const { return dCMx_dCL; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCMx_dCL(su2double val_dcmx_dcl) { dCMx_dCL = val_dcmx_dcl; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCMy_dCL(void) const { return dCMy_dCL; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCMy_dCL(su2double val_dcmy_dcl) { dCMy_dCL = val_dcmy_dcl; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCMz_dCL(void) const { return dCMz_dCL; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCMz_dCL(su2double val_dcmz_dcl) { dCMz_dCL = val_dcmz_dcl; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCL_dAlpha(su2double val_dcl_dalpha) { dCL_dAlpha = val_dcl_dalpha; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCM_diH(su2double val_dcm_dhi) { dCM_diH = val_dcm_dhi; }

  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetCL_Target(void) const { return CL_Target; }

  /*!
   * \brief Set the global parameters of each simulation for each runtime system.
   * \param[in] val_solver - Solver of the simulation.
   * \param[in] val_system - Runtime system that we are solving.
   */
  void SetGlobalParam(MAIN_SOLVER val_solver, unsigned short val_system);

  /*!
   * \brief Center of rotation for a rotational periodic boundary.
   */
  const su2double *GetPeriodicRotCenter(const string& val_marker) const;

  /*!
   * \brief Angles of rotation for a rotational periodic boundary.
   */
  const su2double *GetPeriodicRotAngles(const string& val_marker) const;

  /*!
   * \brief Translation vector for a translational periodic boundary.
   */
  const su2double *GetPeriodicTranslation(const string& val_marker) const;

  /*!
   * \brief Get the translation vector for a periodic transformation.
   * \param[in] val_index - Index corresponding to the periodic transformation.
   * \return The translation vector.
   */
  const su2double* GetPeriodic_Translation(unsigned short val_index ) const { return Periodic_Translation[val_index]; }

  /*!
   * \brief Get the rotationally periodic donor marker for boundary <i>val_marker</i>.
   * \return Periodic donor marker from the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_Periodic_Donor(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_NetThrust(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_Power(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_MassFlow(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_Mach(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_Force(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_BCThrust(const string& val_marker) const;

  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_BCThrust_Old(const string& val_marker) const;

  /*!
   * \brief Get the tip radius of th actuator disk.
   */
  su2double GetActDisk_Area(const string& val_marker) const;

  /*!
   * \brief Get the tip radius of th actuator disk.
   */
  su2double GetActDisk_ReverseMassFlow(const string& val_marker) const;

  /*!
   * \brief Get the thrust corffient of the actuator disk.
   */
  su2double GetActDisk_PressJump(const string& val_marker, unsigned short val_index) const;

  /*!
   * \brief Get the thrust corffient of the actuator disk.
   */
  su2double GetActDisk_TempJump(const string& val_marker, unsigned short val_index) const;

  /*!
   * \brief Get the rev / min of the actuator disk.
   */
  su2double GetActDisk_Omega(const string& val_marker, unsigned short val_index) const;

  /*!
   * \brief Get Actuator Disk Outlet for boundary <i>val_marker</i> (actuator disk inlet).
   * \return Actuator Disk Outlet from the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_ActDiskOutlet(const string& val_marker) const;

  /*!
   * \brief Get Actuator Disk Outlet for boundary <i>val_marker</i> (actuator disk inlet).
   * \return Actuator Disk Outlet from the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_EngineExhaust(const string& val_marker) const;

  /*!
   * \brief Get the internal index for a moving boundary <i>val_marker</i>.
   * \return Internal index for a moving boundary <i>val_marker</i>.
   */
  unsigned short GetMarker_Moving(const string& val_marker) const;

  /*!
   * \brief Get a bool for whether a marker is moving. <i>val_marker</i>.
   * \param[in] val_marker - Name of the marker to test.
   * \return True if the marker is a moving boundary <i>val_marker</i>.
   */
  inline bool GetMarker_Moving_Bool(const string& val_marker) const {
    return GetMarker_Moving(val_marker) < nMarker_Moving;
  }

  /*!
   * \brief Get the internal index for a DEFORM_MESH boundary <i>val_marker</i>.
   * \return Internal index for a DEFORM_MESH boundary <i>val_marker</i>.
   */
  unsigned short GetMarker_Deform_Mesh(const string& val_marker) const;

  /*!
   * \brief Get the internal index for a DEFORM_MESH_SYM_PLANE boundary <i>val_marker</i>.
   * \return Internal index for a DEFORM_MESH_SYM_PLANE boundary <i>val_marker</i>.
   */
  unsigned short GetMarker_Deform_Mesh_Sym_Plane(const string& val_marker) const;

  /*!
   * \brief Get a bool for whether the marker is deformed. <i>val_marker</i>.
   * \param[in] val_marker - Name of the marker to test.
   * \return True if the marker is a deforming boundary <i>val_marker</i>.
   */
  inline bool GetMarker_Deform_Mesh_Bool(const string& val_marker) const {
    return GetMarker_Deform_Mesh(val_marker) < nMarker_Deform_Mesh ||
        GetMarker_Deform_Mesh_Sym_Plane(val_marker) < nMarker_Deform_Mesh_Sym_Plane;
  }

  /*!
   * \brief Get the internal index for a Fluid_Load boundary <i>val_marker</i>.
   * \return Internal index for a Fluid_Load boundary <i>val_marker</i>.
   */
  unsigned short GetMarker_Fluid_Load(const string& val_marker) const;

  /*!
   * \brief Get the internal index for a gradient boundary condition <i>val_marker</i>.
   * \return Internal index for a gradient boundary  condition <i>val_marker</i>.
   */
  unsigned short GetMarker_SobolevBC(const string& val_marker) const;

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Moving_TagBound(unsigned short val_marker) const { return Marker_Moving[val_marker]; }

  /*!
   * \brief Get the name of the DEFORM_MESH boundary defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Deform_Mesh_TagBound(unsigned short val_marker) const { return Marker_Deform_Mesh[val_marker]; }

  /*!
   * \brief Get the name of the DEFORM_MESH_SYM_PLANE boundary defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Deform_Mesh_Sym_Plane_TagBound(unsigned short val_marker) const { return Marker_Deform_Mesh_Sym_Plane[val_marker]; }

  /*!
   * \brief Get the name of the Fluid_Load boundary defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Fluid_Load_TagBound(unsigned short val_marker) const { return Marker_Fluid_Load[val_marker]; }

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_PyCustom_TagBound(unsigned short val_marker) const { return Marker_PyCustom[val_marker]; }

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Analyze_TagBound(unsigned short val_marker) const { return Marker_Analyze[val_marker]; }

  /*!
   * \brief Get the total temperature at a nacelle boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total temperature.
   */
  su2double GetExhaust_Temperature_Target(const string& val_index) const;

  /*!
   * \brief Get the total temperature at an inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total temperature.
   */
  su2double GetInlet_Ttotal(const string& val_index) const;

  /*!
   * \brief Get the temperature at a supersonic inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet density.
   */
  su2double GetInlet_Temperature(const string& val_index) const;

  /*!
   * \brief Get the pressure at a supersonic inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet pressure.
   */
  su2double GetInlet_Pressure(const string& val_index) const;

  /*!
   * \brief Get the velocity vector at a supersonic inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet velocity vector.
   */
  const su2double* GetInlet_Velocity(const string& val_index) const;

  /*!
   * \brief Get the mass fraction vector for a NEMO inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet velocity vector.
   */
  const su2double* GetInlet_MassFrac() const { return Inlet_MassFrac; }

  /*!
   * \brief Get the Tve value for a NEMO inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet velocity vector.
   */
  su2double GetInlet_Temperature_ve() const { return Inlet_Temperature_ve; }

  /*!
   * \brief Get the total pressure at an inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total pressure.
   */
  su2double GetInlet_Ptotal(const string& val_index) const;

  /*!
   * \brief Set the total pressure at an inlet boundary.
   * \param[in] val_pressure - Pressure value at the inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   */
  void SetInlet_Ptotal(su2double val_pressure, const string& val_marker);

  /*!
   * \brief Get the species values at an inlet boundary
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet species values.
   */
  const su2double* GetInlet_SpeciesVal(const string& val_index) const;

  /*!
   * \brief Get the turbulent properties values at an inlet boundary
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet turbulent values.
   */
  const su2double* GetInlet_TurbVal(const string& val_index) const;

  /*!
   * \brief Get the total pressure at an nacelle boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total pressure.
   */
  su2double GetExhaust_Pressure_Target(const string& val_index) const;

  /*!
   * \brief Value of the CFL reduction in turbulence problems.
   * \return Value of the CFL reduction in turbulence problems.
   */
  su2double GetCFLRedCoeff_Turb(void) const { return CFLRedCoeff_Turb; }

  /*!
   * \brief Value of the CFL reduction in species problems.
   * \return Value of the CFL reduction in species problems.
   */
  su2double GetCFLRedCoeff_Species() const { return CFLRedCoeff_Species; }

  /*!
   * \brief Get the flow direction unit vector at an inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The flow direction vector.
   */
  const su2double* GetInlet_FlowDir(const string& val_index) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_Pressure(const string& val_index) const;

  /*!
   * \brief Set the back pressure (static) at an outlet boundary.
   * \param[in] val_pressure - Pressure value at the outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   */
  void SetOutlet_Pressure(su2double val_pressure, const string& val_marker);

  /*!
   * \brief Get the var 1 at Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return The var1
   */
  su2double GetRiemann_Var1(const string& val_marker) const;

  /*!
   * \brief Get the var 2 at Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return The var2
   */
  su2double GetRiemann_Var2(const string& val_marker) const;

  /*!
   * \brief Get the Flowdir at Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return The Flowdir
   */
  const su2double* GetRiemann_FlowDir(const string& val_marker) const;

  /*!
   * \brief Get Kind Data of Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return Kind data
   */
  unsigned short GetKind_Data_Riemann(const string& val_marker) const;

  /*!
   * \brief Get the var 1 for the Giels BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The var1
   */
  su2double GetGiles_Var1(const string& val_marker) const;

  /*!
   * \brief Get the var 2 for the Giles boundary.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The var2
   */
  su2double GetGiles_Var2(const string& val_marker) const;

  /*!
   * \brief Get the Flowdir for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The Flowdir
   */
  const su2double* GetGiles_FlowDir(const string& val_marker) const;

  /*!
   * \brief Get Kind Data for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return Kind data
   */
  unsigned short GetKind_Data_Giles(const string& val_marker) const;

  /*!
   * \brief Set the var 1 for Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   */
  void SetGiles_Var1(su2double newVar1, const string& val_marker);

  /*!
   * \brief Get the relax factor for the average component for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The relax factor for the average component
   */
  su2double GetGiles_RelaxFactorAverage(const string& val_marker) const;

  /*!
   * \brief Get the relax factor for the fourier component for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The relax factor for the fourier component
   */
  su2double GetGiles_RelaxFactorFourier(const string& val_marker) const;

  /*!
   * \brief Get the outlet pressure imposed as BC for internal flow.
   * \return outlet pressure
   */
  su2double GetPressureOut_BC() const;

  /*!
   * \brief Set the outlet pressure imposed as BC for internal flow.
   * \param[in] val_temp - New value of the outlet pressure.
   */
  void SetPressureOut_BC(su2double val_press);

  /*!
   * \brief Get the inlet velocity or pressure imposed for incompressible flow.
   * \return inlet velocity or pressure
   */
  su2double GetIncInlet_BC() const;

  /*!
   * \brief Set the inlet velocity or pressure imposed as BC for incompressible flow.
   * \param[in] val_in - New value of the inlet velocity or pressure.
   */
  void SetIncInlet_BC(su2double val_in);

  /*!
   * \brief Get the inlet temperature imposed as BC for incompressible flow.
   * \return inlet temperature
   */
  su2double GetIncTemperature_BC() const;

  /*!
   * \brief Set the inlet temperature imposed as BC for incompressible flow.
   * \param[in] val_temperature - New value of the inlet temperature.
   */
  void SetIncTemperature_BC(su2double val_temperature);

  /*!
   * \brief Get the outlet pressure imposed as BC for incompressible flow.
   * \return outlet pressure
   */
  su2double GetIncPressureOut_BC() const;

  /*!
   * \brief Set the outlet pressure imposed as BC for incompressible flow.
   * \param[in] val_pressure - New value of the outlet pressure.
   */
  void SetIncPressureOut_BC(su2double val_pressure);

  /*!
   * \brief Get the inlet total pressure imposed as BC for internal flow.
   * \return inlet total pressure
   */
  su2double GetTotalPressureIn_BC() const;

  /*!
   * \brief Get the inlet total temperature imposed as BC for internal flow.
   * \return inlet total temperature
   */
  su2double GetTotalTemperatureIn_BC() const;

  /*!
   * \brief Set the inlet total temperature imposed as BC for internal flow.
   * \param[in] val_temp - New value of the total temperature.
   */
  void SetTotalTemperatureIn_BC(su2double val_temp);

  /*!
   * \brief Get the inlet flow angle imposed as BC for internal flow.
   * \return inlet flow angle
   */
  su2double GetFlowAngleIn_BC() const;

  /*!
   * \brief Get the wall temperature (static) at an isothermal boundary.
   * \param[in] val_index - Index corresponding to the isothermal boundary.
   * \return The wall temperature.
   */
  su2double GetIsothermal_Temperature(const string& val_index) const;

  /*!
   * \brief Get the wall heat flux on a constant heat flux boundary.
   * \param[in] val_index - Index corresponding to the constant heat flux boundary.
   * \return The heat flux.
   */
  su2double GetWall_HeatFlux(const string& val_index) const;

  /*!
   * \brief Get the heat transfer coefficient on a heat transfer boundary.
   * \param[in] val_index - Index corresponding to the heat transfer boundary.
   * \return The heat transfer coefficient.
   */
  su2double GetWall_HeatTransfer_Coefficient(const string& val_index) const;

  /*!
   * \brief Get the temperature at inifinty on a heat transfer boundary.
   * \param[in] val_index - Index corresponding to the heat transfer boundary.
   * \return The temperature at infinity.
   */
  su2double GetWall_HeatTransfer_Temperature(const string& val_index) const;

  /*!
   * \brief Get the wall function treatment for the given boundary marker.
   * \param[in] val_marker - String of the viscous wall marker.
   * \return The type of wall function treatment.
   */
  WALL_FUNCTIONS GetWallFunction_Treatment(const string& val_marker) const;

  /*!
   * \brief Get the additional integer info for the wall function treatment
            for the given boundary marker.
   * \param[in] val_marker - String of the viscous wall marker.
   * \return Pointer to the integer info for the given marker.
   */
  const unsigned short* GetWallFunction_IntInfo(const string& val_marker) const;

  /*!
   * \brief Get the additional double info for the wall function treatment
            for the given boundary marker.
   * \param[in] val_marker - String of the viscous wall marker.
   * \return Pointer to the double info for the given marker.
   */
  const su2double* GetWallFunction_DoubleInfo(const string& val_marker) const;

  /*!
   * \brief Get the type of wall and roughness height on a wall boundary (Heatflux or Isothermal).
   * \param[in] val_index - Index corresponding to the boundary.
   * \return The wall type and roughness height.
   */
  pair<WALL_TYPE,su2double> GetWallRoughnessProperties(const string& val_marker) const;

  /*!
   * \brief Get the target (pressure, massflow, etc) at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \return Target (pressure, massflow, etc) .
   */
  su2double GetEngineInflow_Target(const string& val_marker) const;

  /*!
   * \brief Get the fan face Mach number at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The fan face Mach number.
   */
  su2double GetInflow_Mach(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow pressure.
   */
  su2double GetInflow_Pressure(const string& val_marker) const;

  /*!
   * \brief Get the mass flow rate at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine mass flow rate.
   */
  su2double GetInflow_MassFlow(const string& val_marker) const;

  /*!
   * \brief Get the percentage of reverse flow at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The percentage of reverse flow.
   */
  su2double GetInflow_ReverseMassFlow(const string& val_marker) const;

  /*!
   * \brief Get the percentage of reverse flow at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \return The percentage of reverse flow.
   */
  su2double GetInflow_ReverseMassFlow(unsigned short val_marker) const { return Inflow_ReverseMassFlow[val_marker]; }

  /*!
   * \brief Get the total pressure at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The total pressure.
   */
  su2double GetInflow_TotalPressure(const string& val_marker) const;

  /*!
   * \brief Get the temperature (static) at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow temperature.
   */
  su2double GetInflow_Temperature(const string& val_marker) const;

  /*!
   * \brief Get the total temperature at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow total temperature.
   */
  su2double GetInflow_TotalTemperature(const string& val_marker) const;

  /*!
   * \brief Get the ram drag at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow ram drag.
   */
  su2double GetInflow_RamDrag(const string& val_marker) const;

  /*!
   * \brief Get the force balance at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow force balance.
   */
  su2double GetInflow_Force(const string& val_marker) const;

  /*!
   * \brief Get the power at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow power.
   */
  su2double GetInflow_Power(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust pressure.
   */
  su2double GetExhaust_Pressure(const string& val_marker) const;

  /*!
   * \brief Get the temperature (static) at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust temperature.
   */
  su2double GetExhaust_Temperature(const string& val_marker) const;

  /*!
   * \brief Get the massflow at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust massflow.
   */
  su2double GetExhaust_MassFlow(const string& val_marker) const;

  /*!
   * \brief Get the total pressure at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust total pressure.
   */
  su2double GetExhaust_TotalPressure(const string& val_marker) const;

  /*!
   * \brief Get the total temperature at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The total temperature.
   */
  su2double GetExhaust_TotalTemperature(const string& val_marker) const;

  /*!
   * \brief Get the gross thrust at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return Gross thrust.
   */
  su2double GetExhaust_GrossThrust(const string& val_marker) const;

  /*!
   * \brief Get the force balance at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return Force balance.
   */
  su2double GetExhaust_Force(const string& val_marker) const;

  /*!
   * \brief Get the power at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return Power.
   */
  su2double GetExhaust_Power(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetInflow_Mach(unsigned short val_marker, su2double val_fanface_mach) { Inflow_Mach[val_marker] = val_fanface_mach; }

  /*!
   * \brief Set the fan face static pressure at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_pressure - Fan face static pressure.
   */
  void SetInflow_Pressure(unsigned short val_marker, su2double val_fanface_pressure) { Inflow_Pressure[val_marker] = val_fanface_pressure; }

  /*!
   * \brief Set the massflow at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_massflow - Massflow.
   */
  void SetInflow_MassFlow(unsigned short val_marker, su2double val_fanface_massflow) { Inflow_MassFlow[val_marker] = val_fanface_massflow; }

  /*!
   * \brief Set the reverse flow at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_reversemassflow - reverse flow.
   */
  void SetInflow_ReverseMassFlow(unsigned short val_marker, su2double val_fanface_reversemassflow) { Inflow_ReverseMassFlow[val_marker] = val_fanface_reversemassflow; }

  /*!
   * \brief Set the fan face total pressure at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_totalpressure - Fan face total pressure.
   */
  void SetInflow_TotalPressure(unsigned short val_marker, su2double val_fanface_totalpressure) { Inflow_TotalPressure[val_marker] = val_fanface_totalpressure; }

  /*!
   * \brief Set the fan face static temperature at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_pressure - Fan face static temperature.
   */
  void SetInflow_Temperature(unsigned short val_marker, su2double val_fanface_temperature) { Inflow_Temperature[val_marker] = val_fanface_temperature; }

  /*!
   * \brief Set the fan face total temperature at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_totaltemperature - Fan face total temperature.
   */
  void SetInflow_TotalTemperature(unsigned short val_marker, su2double val_fanface_totaltemperature) { Inflow_TotalTemperature[val_marker] = val_fanface_totaltemperature; }

  /*!
   * \brief Set the ram drag temperature at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_ramdrag - Ram drag value.
   */
  void SetInflow_RamDrag(unsigned short val_marker, su2double val_fanface_ramdrag) { Inflow_RamDrag[val_marker] = val_fanface_ramdrag; }

  /*!
   * \brief Set the force balance at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_force - Fan face force.
   */
  void SetInflow_Force(unsigned short val_marker, su2double val_fanface_force) { Inflow_Force[val_marker] = val_fanface_force; }

  /*!
   * \brief Set the power at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_force - Power.
   */
  void SetInflow_Power(unsigned short val_marker, su2double val_fanface_power) { Inflow_Power[val_marker] = val_fanface_power; }

  /*!
   * \brief Set the back pressure (static) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_pressure - Exhaust static pressure.
   */
  void SetExhaust_Pressure(unsigned short val_marker, su2double val_exhaust_pressure) { Exhaust_Pressure[val_marker] = val_exhaust_pressure; }

  /*!
   * \brief Set the temperature (static) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_temp - Exhaust static temperature.
   */
  void SetExhaust_Temperature(unsigned short val_marker, su2double val_exhaust_temp) { Exhaust_Temperature[val_marker] = val_exhaust_temp; }

  /*!
   * \brief Set the back pressure (static) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_temp - Exhaust static temperature.
   */
  void SetExhaust_MassFlow(unsigned short val_marker, su2double val_exhaust_massflow) { Exhaust_MassFlow[val_marker] = val_exhaust_massflow; }

  /*!
   * \brief Set the back pressure (total) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_totalpressure - Exhaust total pressure.
   */
  void SetExhaust_TotalPressure(unsigned short val_marker, su2double val_exhaust_totalpressure) { Exhaust_TotalPressure[val_marker] = val_exhaust_totalpressure; }

  /*!
   * \brief Set the total temperature at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_totaltemp - Exhaust total temperature.
   */
  void SetExhaust_TotalTemperature(unsigned short val_marker, su2double val_exhaust_totaltemp) { Exhaust_TotalTemperature[val_marker] = val_exhaust_totaltemp; }

  /*!
   * \brief Set the gross thrust at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_grossthrust - Exhaust gross thrust temperature.
   */
  void SetExhaust_GrossThrust(unsigned short val_marker, su2double val_exhaust_grossthrust) { Exhaust_GrossThrust[val_marker] = val_exhaust_grossthrust; }

  /*!
   * \brief Set the force balance at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_force - Exhaust force balance.
   */
  void SetExhaust_Force(unsigned short val_marker, su2double val_exhaust_force) { Exhaust_Force[val_marker] = val_exhaust_force; }

  /*!
   * \brief Set the power at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_power - Exhaust power.
   */
  void SetExhaust_Power(unsigned short val_marker, su2double val_exhaust_power) { Exhaust_Power[val_marker] = val_exhaust_power; }

  /*!
   * \brief Set the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_mach - Exhaust power.
   */
  void SetEngine_Mach(unsigned short val_marker, su2double val_engine_mach) { Engine_Mach[val_marker] = val_engine_mach; }

  /*!
   * \brief Set the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_force - Exhaust power.
   */
  void SetEngine_Force(unsigned short val_marker, su2double val_engine_force) { Engine_Force[val_marker] = val_engine_force; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_power - Exhaust power.
   */
  void SetEngine_Power(unsigned short val_marker, su2double val_engine_power) { Engine_Power[val_marker] = val_engine_power; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_netthrust - Exhaust power.
   */
  void SetEngine_NetThrust(unsigned short val_marker, su2double val_engine_netthrust) { Engine_NetThrust[val_marker] = val_engine_netthrust; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_grossthrust - Exhaust power.
   */
  void SetEngine_GrossThrust(unsigned short val_marker, su2double val_engine_grossthrust) { Engine_GrossThrust[val_marker] = val_engine_grossthrust; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_area - Exhaust power.
   */
  void SetEngine_Area(unsigned short val_marker, su2double val_engine_area) { Engine_Area[val_marker] = val_engine_area; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Mach(unsigned short val_marker) const { return Engine_Mach[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Force(unsigned short val_marker) const { return Engine_Force[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Power(unsigned short val_marker) const { return Engine_Power[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */

  su2double GetEngine_NetThrust(unsigned short val_marker) const { return Engine_NetThrust[val_marker]; }
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */

  su2double GetEngine_GrossThrust(unsigned short val_marker) const { return Engine_GrossThrust[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_marker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Area(unsigned short val_marker) const { return Engine_Area[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Temperature(unsigned short val_marker, su2double val_actdisk_temp) { ActDiskInlet_Temperature[val_marker] = val_actdisk_temp; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_TotalTemperature(unsigned short val_marker, su2double val_actdisk_totaltemp) { ActDiskInlet_TotalTemperature[val_marker] = val_actdisk_totaltemp; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Temperature(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_TotalTemperature(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Temperature(unsigned short val_marker, su2double val_actdisk_temp) { ActDiskOutlet_Temperature[val_marker] = val_actdisk_temp; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_TotalTemperature(unsigned short val_marker, su2double val_actdisk_totaltemp) { ActDiskOutlet_TotalTemperature[val_marker] = val_actdisk_totaltemp; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Temperature(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_TotalTemperature(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_MassFlow(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_MassFlow(unsigned short val_marker, su2double val_actdisk_massflow) { ActDiskInlet_MassFlow[val_marker] = val_actdisk_massflow; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_MassFlow(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_MassFlow(unsigned short val_marker, su2double val_actdisk_massflow) { ActDiskOutlet_MassFlow[val_marker] = val_actdisk_massflow; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Pressure(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_TotalPressure(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_DeltaPress(unsigned short val_marker) const { return ActDisk_DeltaPress[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_DeltaTemp(unsigned short val_marker) const { return ActDisk_DeltaTemp[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_TotalPressRatio(unsigned short val_marker) const { return ActDisk_TotalPressRatio[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_TotalTempRatio(unsigned short val_marker) const { return ActDisk_TotalTempRatio[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_StaticPressRatio(unsigned short val_marker) const { return ActDisk_StaticPressRatio[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_StaticTempRatio(unsigned short val_marker) const { return ActDisk_StaticTempRatio[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_NetThrust(unsigned short val_marker) const { return ActDisk_NetThrust[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_BCThrust(unsigned short val_marker) const { return ActDisk_BCThrust[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_BCThrust_Old(unsigned short val_marker) const { return ActDisk_BCThrust_Old[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_GrossThrust(unsigned short val_marker) const { return ActDisk_GrossThrust[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Area(unsigned short val_marker) const { return ActDisk_Area[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_ReverseMassFlow(unsigned short val_marker) const { return ActDisk_ReverseMassFlow[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_RamDrag(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Force(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Power(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Pressure(unsigned short val_marker, su2double val_actdisk_press) { ActDiskInlet_Pressure[val_marker] = val_actdisk_press; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_TotalPressure(unsigned short val_marker, su2double val_actdisk_totalpress) { ActDiskInlet_TotalPressure[val_marker] = val_actdisk_totalpress; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_DeltaPress(unsigned short val_marker, su2double val_actdisk_deltapress) { ActDisk_DeltaPress[val_marker] = val_actdisk_deltapress; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Power(unsigned short val_marker, su2double val_actdisk_power) { ActDisk_Power[val_marker] = val_actdisk_power; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_MassFlow(unsigned short val_marker, su2double val_actdisk_massflow) { ActDisk_MassFlow[val_marker] = val_actdisk_massflow; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Mach(unsigned short val_marker, su2double val_actdisk_mach) { ActDisk_Mach[val_marker] = val_actdisk_mach; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Force(unsigned short val_marker, su2double val_actdisk_force) { ActDisk_Force[val_marker] = val_actdisk_force; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_MassFlow(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetOutlet_MassFlow(unsigned short val_marker, su2double val_massflow) { Outlet_MassFlow[val_marker] = val_massflow; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_Density(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetOutlet_Density(unsigned short val_marker, su2double val_density) { Outlet_Density[val_marker] = val_density; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_Area(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetOutlet_Area(unsigned short val_marker, su2double val_area) { Outlet_Area[val_marker] = val_area; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_DC60(unsigned short val_marker, su2double val_surface_distortion) { Surface_DC60[val_marker] = val_surface_distortion; }

  /*!
   * \brief Set the massflow at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the mass flow.
   */
  void SetSurface_MassFlow(unsigned short val_marker, su2double val_surface_massflow) { Surface_MassFlow[val_marker] = val_surface_massflow; }

  /*!
   * \brief Set the mach number at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the mach number.
   */
  void SetSurface_Mach(unsigned short val_marker, su2double val_surface_mach) { Surface_Mach[val_marker] = val_surface_mach; }

  /*!
   * \brief Set the temperature at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the temperature.
   */
  void SetSurface_Temperature(unsigned short val_marker, su2double val_surface_temperature) { Surface_Temperature[val_marker] = val_surface_temperature; }

  /*!
   * \brief Set the pressure at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the pressure.
   */
  void SetSurface_Pressure(unsigned short val_marker, su2double val_surface_pressure) { Surface_Pressure[val_marker] = val_surface_pressure; }

  /*!
   * \brief Set the density at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_density - Value of the density.
   */
  void SetSurface_Density(unsigned short val_marker, su2double val_surface_density) { Surface_Density[val_marker] = val_surface_density; }

  /*!
   * \brief Set the enthalpy at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_density - Value of the density.
   */
  void SetSurface_Enthalpy(unsigned short val_marker, su2double val_surface_enthalpy) { Surface_Enthalpy[val_marker] = val_surface_enthalpy; }

  /*!
   * \brief Set the normal velocity at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_normalvelocity - Value of the normal velocity.
   */
  void SetSurface_NormalVelocity(unsigned short val_marker, su2double val_surface_normalvelocity) { Surface_NormalVelocity[val_marker] = val_surface_normalvelocity; }

  /*!
   * \brief Set the streamwise flow uniformity at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_streamwiseuniformity - Value of the streamwise flow uniformity.
   */
  void SetSurface_Uniformity(unsigned short val_marker, su2double val_surface_streamwiseuniformity) { Surface_Uniformity[val_marker] = val_surface_streamwiseuniformity; }

  /*!
   * \brief Set the secondary flow strength at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_secondarystrength - Value of the secondary flow strength.
   */
  void SetSurface_SecondaryStrength(unsigned short val_marker, su2double val_surface_secondarystrength) { Surface_SecondaryStrength[val_marker] = val_surface_secondarystrength; }

  /*!
   * \brief Set the relative secondary flow strength at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_secondaryoverstream - Value of the relative seondary flow strength.
   */
  void SetSurface_SecondOverUniform(unsigned short val_marker, su2double val_surface_secondaryoverstream) { Surface_SecondOverUniform[val_marker] = val_surface_secondaryoverstream; }

  /*!
   * \brief Set the momentum distortion at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_momentumdistortion - Value of the momentum distortion.
   */
  void SetSurface_MomentumDistortion(unsigned short val_marker, su2double val_surface_momentumdistortion) { Surface_MomentumDistortion[val_marker] = val_surface_momentumdistortion; }

  /*!
   * \brief Set the total temperature at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_totaltemperature - Value of the total temperature.
   */
  void SetSurface_TotalTemperature(unsigned short val_marker, su2double val_surface_totaltemperature) { Surface_TotalTemperature[val_marker] = val_surface_totaltemperature; }

  /*!
   * \brief Set the total pressure at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_totalpressure - Value of the total pressure.
   */
  void SetSurface_TotalPressure(unsigned short val_marker, su2double val_surface_totalpressure) { Surface_TotalPressure[val_marker] = val_surface_totalpressure; }

  /*!
   * \brief Set the pressure drop between two surfaces.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_pressuredrop - Value of the pressure drop.
   */
  void SetSurface_PressureDrop(unsigned short val_marker, su2double val_surface_pressuredrop) { Surface_PressureDrop[val_marker] = val_surface_pressuredrop; }

  /*!
   * \brief Set the average of species_0 at the surface.
   * \param[in] val_marker - Index corresponding to boundary.
   * \param[in] val_surface_species_0 - Value of avg species_0.
   */
  void SetSurface_Species_0(unsigned short val_marker, su2double val_surface_species_0) { Surface_Species_0[val_marker] = val_surface_species_0; }

  /*!
   * \brief Set the species variance at the surface.
   * \param[in] val_marker - Index corresponding to boundary.
   * \param[in] val_surface_species_variance - Value of the species variance.
   */
  void SetSurface_Species_Variance(unsigned short val_marker, su2double val_surface_species_variance) { Surface_Species_Variance[val_marker] = val_surface_species_variance; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_IDC(unsigned short val_marker, su2double val_surface_distortion) { Surface_IDC[val_marker] = val_surface_distortion; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_IDC_Mach(unsigned short val_marker, su2double val_surface_distortion) { Surface_IDC_Mach[val_marker] = val_surface_distortion; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_IDR(unsigned short val_marker, su2double val_surface_distortion) { Surface_IDR[val_marker] = val_surface_distortion; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_DeltaTemp(unsigned short val_marker, su2double val_actdisk_deltatemp) { ActDisk_DeltaTemp[val_marker] = val_actdisk_deltatemp; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_TotalPressRatio(unsigned short val_marker, su2double val_actdisk_pressratio) { ActDisk_TotalPressRatio[val_marker] = val_actdisk_pressratio; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_TotalTempRatio(unsigned short val_marker, su2double val_actdisk_tempratio) { ActDisk_TotalTempRatio[val_marker] = val_actdisk_tempratio; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_StaticPressRatio(unsigned short val_marker, su2double val_actdisk_pressratio) { ActDisk_StaticPressRatio[val_marker] = val_actdisk_pressratio; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_StaticTempRatio(unsigned short val_marker, su2double val_actdisk_tempratio) { ActDisk_StaticTempRatio[val_marker] = val_actdisk_tempratio; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_NetThrust(unsigned short val_marker, su2double val_actdisk_netthrust) { ActDisk_NetThrust[val_marker] = val_actdisk_netthrust; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust(const string& val_marker, su2double val_actdisk_bcthrust);

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust(unsigned short val_marker, su2double val_actdisk_bcthrust) { ActDisk_BCThrust[val_marker] = val_actdisk_bcthrust; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust_Old(const string& val_marker, su2double val_actdisk_bcthrust_old);

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust_Old(unsigned short val_marker, su2double val_actdisk_bcthrust_old) { ActDisk_BCThrust_Old[val_marker] = val_actdisk_bcthrust_old; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_GrossThrust(unsigned short val_marker, su2double val_actdisk_grossthrust) { ActDisk_GrossThrust[val_marker] = val_actdisk_grossthrust; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Area(unsigned short val_marker, su2double val_actdisk_area) { ActDisk_Area[val_marker] = val_actdisk_area; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_ReverseMassFlow(unsigned short val_marker, su2double val_actdisk_area) { ActDisk_ReverseMassFlow[val_marker] = val_actdisk_area; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_RamDrag(unsigned short val_marker, su2double val_actdisk_ramdrag) { ActDiskInlet_RamDrag[val_marker] = val_actdisk_ramdrag; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Force(unsigned short val_marker, su2double val_actdisk_force) { ActDiskInlet_Force[val_marker] = val_actdisk_force; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Power(unsigned short val_marker, su2double val_actdisk_power) { ActDiskInlet_Power[val_marker] = val_actdisk_power; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Power(unsigned short val_marker) const { return ActDisk_Power[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_MassFlow(unsigned short val_marker) const { return ActDisk_MassFlow[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Mach(unsigned short val_marker) const { return ActDisk_Mach[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Force(unsigned short val_marker) const { return ActDisk_Force[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_DC60(unsigned short val_marker) const { return Surface_DC60[val_marker]; }

  /*!
   * \brief Get the massflow at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The massflow.
   */
  su2double GetSurface_MassFlow(unsigned short val_marker) const { return Surface_MassFlow[val_marker]; }

  /*!
   * \brief Get the mach number at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The mach number.
   */
  su2double GetSurface_Mach(unsigned short val_marker) const { return Surface_Mach[val_marker]; }

  /*!
   * \brief Get the temperature at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The temperature.
   */
  su2double GetSurface_Temperature(unsigned short val_marker) const { return Surface_Temperature[val_marker]; }

  /*!
   * \brief Get the pressure at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The pressure.
   */
  su2double GetSurface_Pressure(unsigned short val_marker) const { return Surface_Pressure[val_marker]; }

  /*!
   * \brief Get the density at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The density.
   */
  su2double GetSurface_Density(unsigned short val_marker) const { return Surface_Density[val_marker]; }

  /*!
   * \brief Get the enthalpy at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The density.
   */
  su2double GetSurface_Enthalpy(unsigned short val_marker) const { return Surface_Enthalpy[val_marker]; }

  /*!
   * \brief Get the normal velocity at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The normal velocity.
   */
  su2double GetSurface_NormalVelocity(unsigned short val_marker) const { return Surface_NormalVelocity[val_marker]; }

  /*!
   * \brief Get the streamwise flow uniformity at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \return The streamwise flow uniformity.
   */
  su2double GetSurface_Uniformity(unsigned short val_marker) const { return Surface_Uniformity[val_marker]; }

  /*!
   * \brief Get the secondary flow strength at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \return The secondary flow strength.
   */
  su2double GetSurface_SecondaryStrength(unsigned short val_marker) const { return Surface_SecondaryStrength[val_marker]; }

  /*!
   * \brief Get the relative secondary flow strength at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \return The relative seondary flow strength.
   */
  su2double GetSurface_SecondOverUniform(unsigned short val_marker) const { return Surface_SecondOverUniform[val_marker]; }

  /*!
   * \brief Get the momentum distortion at the surface.
   * \param[in] val_marker - Index corresponding to the outlet boundary.
   * \return The momentum distortion.
   */
  su2double GetSurface_MomentumDistortion(unsigned short val_marker) const { return Surface_MomentumDistortion[val_marker]; }

  /*!
   * \brief Get the total temperature at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The total temperature.
   */
  su2double GetSurface_TotalTemperature(unsigned short val_marker) const { return Surface_TotalTemperature[val_marker]; }

  /*!
   * \brief Get the total pressure at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The total pressure.
   */
  su2double GetSurface_TotalPressure(unsigned short val_marker) const { return Surface_TotalPressure[val_marker]; }

  /*!
   * \brief Get the pressure drop between two surfaces.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The pressure drop.
   */
  su2double GetSurface_PressureDrop(unsigned short val_marker) const { return Surface_PressureDrop[val_marker]; }

  /*!
   * \brief Get avg species_0 at a boundary.
   * \param[in] val_index - Index corresponding to the boundary.
   * \return The avg species_0.
   */
  su2double GetSurface_Species_0(unsigned short val_marker) const { return Surface_Species_0[val_marker]; }

  /*!
   * \brief Get the species variance at a boundary.
   * \param[in] val_index - Index corresponding to the boundary.
   * \return The species variance.
   */
  su2double GetSurface_Species_Variance(unsigned short val_marker) const { return Surface_Species_Variance[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_IDC(unsigned short val_marker) const { return Surface_IDC[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_IDC_Mach(unsigned short val_marker) const { return Surface_IDC_Mach[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_IDR(unsigned short val_marker) const { return Surface_IDR[val_marker]; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Pressure(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_TotalPressure(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_GrossThrust(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Force(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Power(const string& val_marker) const;

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Pressure(unsigned short val_marker, su2double val_actdisk_press) { ActDiskOutlet_Pressure[val_marker] = val_actdisk_press; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_TotalPressure(unsigned short val_marker, su2double val_actdisk_totalpress) { ActDiskOutlet_TotalPressure[val_marker] = val_actdisk_totalpress; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_GrossThrust(unsigned short val_marker, su2double val_actdisk_grossthrust) { ActDiskOutlet_GrossThrust[val_marker] = val_actdisk_grossthrust; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Force(unsigned short val_marker, su2double val_actdisk_force) { ActDiskOutlet_Force[val_marker] = val_actdisk_force; }

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Power(unsigned short val_marker, su2double val_actdisk_power) { ActDiskOutlet_Power[val_marker] = val_actdisk_power; }

  /*!
   * \brief Get the displacement value at an displacement boundary.
   * \param[in] val_index - Index corresponding to the displacement boundary.
   * \return The displacement value.
   */
  su2double GetDispl_Value(const string& val_index) const;

  /*!
   * \brief Get the force value at an load boundary.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetLoad_Value(const string& val_index) const;

  /*!
   * \brief Get the constant value at a damper boundary.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The damper constant.
   */
  su2double GetDamper_Constant(const string& val_index) const;

  /*!
   * \brief Get the force value at a load boundary defined in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetLoad_Dir_Value(const string& val_index) const;

  /*!
   * \brief Get the force multiplier at a load boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load multiplier.
   */
  su2double GetLoad_Dir_Multiplier(const string& val_index) const;

  /*!
   * \brief Get the force value at a load boundary defined in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetDisp_Dir_Value(const string& val_index) const;

  /*!
   * \brief Get the force multiplier at a load boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load multiplier.
   */
  su2double GetDisp_Dir_Multiplier(const string& val_index) const;

  /*!
   * \brief Get the force direction at a loaded boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load direction.
   */
  const su2double* GetLoad_Dir(const string& val_index) const;

  /*!
   * \brief Get the force direction at a loaded boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load direction.
   */
  const su2double* GetDisp_Dir(const string& val_index) const;

  /*!
   * \brief Get the amplitude of the sine-wave at a load boundary defined in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetLoad_Sine_Amplitude(const string& val_index) const;

  /*!
   * \brief Get the frequency of the sine-wave at a load boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load frequency.
   */
  su2double GetLoad_Sine_Frequency(const string& val_index) const;

  /*!
   * \brief Get the force direction at a sine-wave loaded boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load direction.
   */
  const su2double* GetLoad_Sine_Dir(const string& val_index) const;

  /*!
   * \brief Get the force value at an load boundary.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetFlowLoad_Value(const string& val_index) const;

  /*!
   * \brief Cyclic pitch amplitude for rotor blades.
   * \return The specified cyclic pitch amplitude.
   */
  su2double GetCyclic_Pitch(void) const { return Cyclic_Pitch; }

  /*!
   * \brief Collective pitch setting for rotor blades.
   * \return The specified collective pitch setting.
   */
  su2double GetCollective_Pitch(void) const { return Collective_Pitch; }

  /*!
   * \brief Get name of the arbitrary mesh motion input file.
   * \return File name of the arbitrary mesh motion input file.
   */
  string GetDV_Filename(void) const { return DV_Filename; }

  /*!
   * \brief Get name of the unordered ASCII volume sensitivity file.
   * \return File name of the unordered ASCII volume sensitivity file.
   */
  string GetDV_Unordered_Sens_Filename(void) const { return DV_Unordered_Sens_Filename; }

  /*!
   * \brief Get name of the unordered ASCII surface sensitivity file.
   * \return File name of the unordered ASCII surface sensitivity file.
   */
  string GetDV_Sens_Filename(void) const { return DV_Sens_Filename; }

  /*!
   * \brief Set the config options.
   */
  void SetConfig_Options();

  /*!
   * \brief Set the config file parsing.
   */
  void SetConfig_Parsing(char case_filename[MAX_STRING_SIZE]);

  /*!
   * \brief Set the config file parsing.
   */
  void SetConfig_Parsing(istream &config_buffer);

  /*!
   * \brief Set the config file parsing.
   */
  bool SetRunTime_Parsing(char case_filename[MAX_STRING_SIZE]);

  /*!
   * \brief Config file postprocessing.
   */
  void SetPostprocessing(SU2_COMPONENT val_software, unsigned short val_izone, unsigned short val_nDim);

  /*!
   * \brief Config file markers processing.
   */
  void SetMarkers(SU2_COMPONENT val_software);

  /*!
   * \brief Config file output.
   */
  void SetOutput(SU2_COMPONENT val_software, unsigned short val_izone);

  /*!
   * \brief Value of Aeroelastic solution coordinate at time n+1.
   */
  vector<vector<su2double> > GetAeroelastic_np1(unsigned short iMarker) const { return Aeroelastic_np1[iMarker]; }

  /*!
   * \brief Value of Aeroelastic solution coordinate at time n.
   */
  vector<vector<su2double> > GetAeroelastic_n(unsigned short iMarker) const { return Aeroelastic_n[iMarker]; }

  /*!
   * \brief Value of Aeroelastic solution coordinate at time n-1.
   */
  vector<vector<su2double> > GetAeroelastic_n1(unsigned short iMarker) const { return Aeroelastic_n1[iMarker]; }

  /*!
   * \brief Value of Aeroelastic solution coordinate at time n+1.
   */
  void SetAeroelastic_np1(unsigned short iMarker, vector<vector<su2double> > solution) { Aeroelastic_np1[iMarker] = solution;}

  /*!
   * \brief Value of Aeroelastic solution coordinate at time n from time n+1.
   */
  void SetAeroelastic_n(void) { Aeroelastic_n = Aeroelastic_np1; }

  /*!
   * \brief Value of Aeroelastic solution coordinate at time n-1 from time n.
   */
  void SetAeroelastic_n1(void) { Aeroelastic_n1 = Aeroelastic_n; }

  /*!
   * \brief Aeroelastic Flutter Speed Index.
   */
  su2double GetAeroelastic_Flutter_Speed_Index(void) const { return FlutterSpeedIndex; }

  /*!
   * \brief Uncoupled Aeroelastic Frequency Plunge.
   */
  su2double GetAeroelastic_Frequency_Plunge(void) const { return PlungeNaturalFrequency; }

  /*!
   * \brief Uncoupled Aeroelastic Frequency Pitch.
   */
  su2double GetAeroelastic_Frequency_Pitch(void) const { return PitchNaturalFrequency; }

  /*!
   * \brief Aeroelastic Airfoil Mass Ratio.
   */
  su2double GetAeroelastic_Airfoil_Mass_Ratio(void) const { return AirfoilMassRatio; }

  /*!
   * \brief Aeroelastic center of gravity location.
   */
  su2double GetAeroelastic_CG_Location(void) const { return CG_Location; }

  /*!
   * \brief Aeroelastic radius of gyration squared.
   */
  su2double GetAeroelastic_Radius_Gyration_Squared(void) const { return RadiusGyrationSquared; }

  /*!
   * \brief Aeroelastic solve every x inner iteration.
   */
  unsigned short GetAeroelasticIter(void) const { return AeroelasticIter; }

  /*!
   * \brief Value of plunging coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Value of plunging coordinate.
   */
  su2double GetAeroelastic_plunge(unsigned short val_marker) const { return Aeroelastic_plunge[val_marker]; }

  /*!
   * \brief Value of pitching coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Value of pitching coordinate.
   */
  su2double GetAeroelastic_pitch(unsigned short val_marker) const { return Aeroelastic_pitch[val_marker]; }

  /*!
   * \brief Value of plunging coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val - value of plunging coordinate.
   */
  void SetAeroelastic_plunge(unsigned short val_marker, su2double val) { Aeroelastic_plunge[val_marker] = val; }

  /*!
   * \brief Value of pitching coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val - value of pitching coordinate.
   */
  void SetAeroelastic_pitch(unsigned short val_marker, su2double val) { Aeroelastic_pitch[val_marker] = val; }

  /*!
   * \brief Get information about the aeroelastic simulation.
   * \return <code>TRUE</code> if it is an aeroelastic case; otherwise <code>FALSE</code>.
   */
  bool GetAeroelastic_Simulation(void) const { return Aeroelastic_Simulation; }

  /*!
   * \brief Get information about the wind gust.
   * \return <code>TRUE</code> if there is a wind gust; otherwise <code>FALSE</code>.
   */
  bool GetWind_Gust(void) const { return Wind_Gust; }

  /*!
   * \brief Get the type of gust to simulate.
   * \return type of gust to use for the simulation.
   */
  unsigned short GetGust_Type(void) const { return Gust_Type; }

  /*!
   * \brief Get the gust direction.
   * \return the gust direction.
   */
  unsigned short GetGust_Dir(void) const { return Gust_Dir; }

  /*!
   * \brief Value of the gust wavelength.
   */
  su2double GetGust_WaveLength(void) const { return Gust_WaveLength; }

  /*!
   * \brief Value of the number of gust periods.
   */
  su2double GetGust_Periods(void) const { return Gust_Periods; }

  /*!
   * \brief Value of the gust amplitude.
   */
  su2double GetGust_Ampl(void) const { return Gust_Ampl; }

  /*!
   * \brief Value of the time at which to begin the gust.
   */
  su2double GetGust_Begin_Time(void) const { return Gust_Begin_Time; }

  /*!
   * \brief Value of the location ath which the gust begins.
   */
  su2double GetGust_Begin_Loc(void) const { return Gust_Begin_Loc; }

  /*!
   * \brief Get whether fixed values for turbulence quantities are applied.
   * \return <code>TRUE</code> if fixed values are applied; otherwise <code>FALSE</code>.
   */
  bool GetTurb_Fixed_Values(void) const { return Turb_Fixed_Values; }

  /*!
   * \brief Get shift of the upstream half-plane where fixed values for turbulence quantities are applied.
   * \details This half-plane is given by the condition that the dot product between the
   * coordinate vector and the normalized far-field velocity vector is less than what this
   * function returns.
   */
  su2double GetTurb_Fixed_Values_MaxScalarProd(void) const { return Turb_Fixed_Values_MaxScalarProd; }

  /*!
   * \brief Get the number of iterations to evaluate the parametric coordinates.
   * \return Number of iterations to evaluate the parametric coordinates.
   */
  unsigned short GetnFFD_Iter(void) const { return nFFD_Iter; }

  /*!
   * \brief Get the tolerance of the point inversion algorithm.
   * \return Tolerance of the point inversion algorithm.
   */
  su2double GetFFD_Tol(void) const { return FFD_Tol; }

  /*!
   * \brief Get information about whether to do a check on self-intersections within
      the FFD box based on value on the Jacobian determinant.
   * \param[out] FFD_IntPrev: <code>TRUE</code> if FFD intersection prevention is active; otherwise <code>FALSE</code>.
   * \param[out] FFD_IntPrev_MaxIter: Maximum number of iterations in the intersection prevention procedure.
   * \param[out] FFD_IntPrev_MaxDepth: Maximum recursion depth in the intersection prevention procedure.
   */
  tuple<bool, unsigned short, unsigned short> GetFFD_IntPrev(void) const {
    return make_tuple(FFD_IntPrev, FFD_IntPrev_MaxIter, FFD_IntPrev_MaxDepth);
  }

  /*!
   * \brief Get information about whether to do a check on convexity of the mesh elements.
   * \param[out] ConvexityCheck: <code>TRUE</code> if convexity check is active; otherwise <code>FALSE</code>.
   * \param[out] ConvexityCheck_MaxIter: Maximum number of iterations in the convexity check.
   * \param[out] ConvexityCheck_MaxDepth: Maximum recursion depth in the convexity check.
   */
  tuple<bool, unsigned short, unsigned short> GetConvexityCheck(void) const {
    return make_tuple(ConvexityCheck, ConvexityCheck_MaxIter, ConvexityCheck_MaxDepth);
  }

  /*!
   * \brief Get the scale factor for the line search.
   * \return Scale factor for the line search.
   */
  su2double GetOpt_RelaxFactor(void) const { return Opt_RelaxFactor; }

  /*!
   * \brief Get the bound for the line search.
   * \return Bound for the line search.
   */
  su2double GetOpt_LineSearch_Bound(void) const { return Opt_LineSearch_Bound; }

  /*!
   * \brief Set the scale factor for the line search.
   * \param[in] val_scale - scale of the deformation.
   */
  void SetOpt_RelaxFactor(su2double val_scale) { Opt_RelaxFactor = val_scale; }

  /*!
   * \brief Get the node number of the CV to visualize.
   * \return Node number of the CV to visualize.
   */
  long GetVisualize_CV(void) const { return Visualize_CV; }

  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetFixed_CL_Mode(void) const { return Fixed_CL_Mode; }

  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetEval_dOF_dCX(void) const { return Eval_dOF_dCX; }

  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetDiscard_InFiles(void) const { return Discard_InFiles; }

  /*!
   * \brief Get the value specified for the target CL.
   * \return Value of the target CL.
   */
  su2double GetTarget_CL(void) const { return Target_CL; }

  /*!
   * \brief Get the value for the lift curve slope for fixed CL mode.
   * \return Lift curve slope for fixed CL mode.
   */
  su2double GetdCL_dAlpha(void) const { return dCL_dAlpha; }

  /*!
   * \brief Number of iterations to evaluate dCL_dAlpha.
   * \return Number of iterations.
   */
  unsigned long GetIter_dCL_dAlpha(void) const { return Iter_dCL_dAlpha; }

  /*!
   * \brief Get the value of the damping coefficient for fixed CL mode.
   * \return Damping coefficient for fixed CL mode.
   */
  su2double GetdCM_diH(void) const { return dCM_diH; }

  /*!
   * \brief Get the value of iterations to re-evaluate the angle of attack.
   * \return Number of iterations.
   */
  unsigned long GetIter_Fixed_NetThrust(void) const { return Iter_Fixed_NetThrust; }

  /*!
   * \brief Get the value of the damping coefficient for fixed CL mode.
   * \return Damping coefficient for fixed CL mode.
   */
  su2double GetdNetThrust_dBCThrust(void) const { return dNetThrust_dBCThrust; }

  /*!
   * \brief Get the value of iterations to re-evaluate the angle of attack.
   * \return Number of iterations.
   */
  unsigned long GetUpdate_BCThrust(void) const { return Update_BCThrust; }

  /*!
   * \brief Set the value of the boolean for updating AoA in fixed lift mode.
   * \param[in] val_update - the bool for whether to update the AoA.
   */
  void SetUpdate_BCThrust_Bool(bool val_update) { Update_BCThrust_Bool = val_update; }

  /*!
   * \brief Set the value of the boolean for updating AoA in fixed lift mode.
   * \param[in] val_update - the bool for whether to update the AoA.
   */
  void SetUpdate_AoA(bool val_update) { Update_AoA = val_update; }

  /*!
   * \brief Get information about whether to update the AoA for fixed lift mode.
   * \return <code>TRUE</code> if we should update the AoA for fixed lift mode; otherwise <code>FALSE</code>.
   */
  bool GetUpdate_BCThrust_Bool(void) const { return Update_BCThrust_Bool; }

  /*!
   * \brief Get information about whether to update the AoA for fixed lift mode.
   * \return <code>TRUE</code> if we should update the AoA for fixed lift mode; otherwise <code>FALSE</code>.
   */
  bool GetUpdate_AoA(void) const { return Update_AoA; }

  /*!
   * \brief Get the maximum number of iterations between AoA updates for fixed C_L mode
   * \return Number of maximum iterations between AoA updates
   */
  unsigned long GetUpdate_AoA_Iter_Limit(void) const { return Update_AoA_Iter_Limit; }

  /*!
   * \brief Get whether at the end of finite differencing (Fixed CL mode)
   * \return boolean indicating end of finite differencing mode (Fixed CL mode)
   */
  bool GetFinite_Difference_Mode(void) const { return Finite_Difference_Mode; }

  /*!
   * \brief Set whether at the end of finite differencing (Fixed CL mode)
   */
  void SetFinite_Difference_Mode(bool val_fd_mode) { Finite_Difference_Mode = val_fd_mode; }

  /*!
   * \brief Set the current number of non-physical nodes in the solution.
   * \param[in] val_nonphys_points - current number of non-physical points.
   */
  void SetNonphysical_Points(unsigned long val_nonphys_points) { Nonphys_Points = val_nonphys_points; }

  /*!
   * \brief Get the current number of non-physical nodes in the solution.
   * \return Current number of non-physical points.
   */
  unsigned long GetNonphysical_Points(void) const { return Nonphys_Points; }

  /*!
   * \brief Set the current number of non-physical reconstructions for 2nd-order upwinding.
   * \param[in] val_nonphys_reconstr - current number of non-physical reconstructions for 2nd-order upwinding.
   */
  void SetNonphysical_Reconstr(unsigned long val_nonphys_reconstr) { Nonphys_Reconstr = val_nonphys_reconstr; }

  /*!
   * \brief Get the current number of non-physical reconstructions for 2nd-order upwinding.
   * \return Current number of non-physical reconstructions for 2nd-order upwinding.
   */
  unsigned long GetNonphysical_Reconstr(void) const { return Nonphys_Reconstr; }

  /*!
   * \brief Start the timer for profiling subroutines.
   * \param[in] val_start_time - the value of the start time.
   */
  void Tick(double *val_start_time);

  /*!
   * \brief Stop the timer for profiling subroutines and store results.
   * \param[in] val_start_time - the value of the start time.
   * \param[in] val_function_name - string for the name of the profiled subroutine.
   * \param[in] val_group_id - string for the name of the profiled subroutine.
   */
  void Tock(double val_start_time, const string& val_function_name, int val_group_id);

  /*!
   * \brief Write a CSV file containing the results of the profiling.
   */
  void SetProfilingCSV(void);

  /*!
   * \brief Start the timer for profiling subroutines.
   * \param[in] val_start_time - the value of the start time.
   */
  void GEMM_Tick(double *val_start_time) const;

  /*!
   * \brief Stop the timer for the GEMM profiling and store results.
   * \param[in] val_start_time - The value of the start time.
   * \param[in] M, N, K        - Matrix size of the GEMM call.
   */
  void GEMM_Tock(double val_start_time, int M, int N, int K) const;

  /*!
   * \brief Write a CSV file containing the results of the profiling.
   */
  void GEMMProfilingCSV(void);

  /*!
   * \brief Set freestream turbonormal for initializing solution.
   */
  void SetFreeStreamTurboNormal(const su2double* turboNormal);

  /*!
   * \brief Set freestream turbonormal for initializing solution.
   */
  const su2double* GetFreeStreamTurboNormal(void) const { return FreeStreamTurboNormal; }

  /*!
   * \brief Set multizone properties.
   */
  void SetMultizone(const CConfig *driver_config, const CConfig* const* config_container);

  /*!
   * \brief Get the verbosity level of the console output.
   * \return Verbosity level for the console output.
   */
  unsigned short GetConsole_Output_Verb(void) const { return Console_Output_Verb; }

  /*!
   * \brief Get the kind of marker analyze marker (area-averaged, mass flux averaged, etc).
   * \return Kind of average.
   */
  unsigned short GetKind_Average(void) const { return Kind_Average; }

  /*!
   * \brief Get the direct differentation method.
   * \return direct differentiation method.
   */
  unsigned short GetDirectDiff() const { return DirectDiff;}

  /*!
   * \brief Get the indicator whether we are solving an discrete adjoint problem.
   * \return the discrete adjoint indicator.
   */
  bool GetDiscrete_Adjoint(void) const { return DiscreteAdjoint; }

  /*!
   * \brief Get the number of subiterations while a ramp is applied.
   * \return Number of FSI subiters.
   */
  unsigned short GetnIterFSI_Ramp(void) const { return nIterFSI_Ramp; }

  /*!
   * \brief Get Aitken's relaxation parameter for static relaxation cases.
   * \return Aitken's relaxation parameters.
   */
  su2double GetAitkenStatRelax(void) const { return AitkenStatRelax; }

  /*!
   * \brief Get Aitken's maximum relaxation parameter for dynamic relaxation cases and first iteration.
   * \return Aitken's relaxation parameters.
   */
  su2double GetAitkenDynMaxInit(void) const { return AitkenDynMaxInit; }

  /*!
   * \brief Get Aitken's maximum relaxation parameter for dynamic relaxation cases and first iteration.
   * \return Aitken's relaxation parameters.
   */
  su2double GetAitkenDynMinInit(void) const { return AitkenDynMinInit; }

  /*!
   * \brief Decide whether to apply dead loads to the model.
   * \return <code>TRUE</code> if the dead loads are to be applied, <code>FALSE</code> otherwise.
   */
  bool GetDeadLoad(void) const { return DeadLoad; }

  /*!
   * \brief Identifies if the mesh is matching or not (temporary, while implementing interpolation procedures).
   * \return <code>TRUE</code> if the mesh is matching, <code>FALSE</code> otherwise.
   */
  bool GetPseudoStatic(void) const { return PseudoStatic; }

  /*!
   * \brief Identifies if we want to restart from a steady or an unsteady solution.
   * \return <code>TRUE</code> if we restart from steady state solution, <code>FALSE</code> otherwise.
   */
  bool GetSteadyRestart(void) const { return SteadyRestart; }

  /*!
   * \brief Get the current instance.
   * \return Current instance identifier.
   */
  unsigned short GetiInst(void) const { return iInst; }

  /*!
   * \brief Set the current instance.
   * \param[in] iInst - current instance identifier.
   */
  void SetiInst(unsigned short val_iInst) { iInst = val_iInst; }

  /*!
   * \brief Get Newmark alpha parameter.
   * \return Value of the Newmark alpha parameter.
   */
  su2double GetNewmark_beta(void) const { return Newmark_beta; }

  /*!
   * \brief Get Newmark delta parameter.
   * \return Value of the Newmark delta parameter.
   */
  su2double GetNewmark_gamma(void) const { return Newmark_gamma; }

  /*!
   * \brief Get the number of integration coefficients provided by the user.
   * \return Number of integration coefficients.
   */
  unsigned short GetnIntCoeffs(void) const { return nIntCoeffs; }

  /*!
   * \brief Get the number of different materials for the elasticity solver.
   * \return Number of different materials.
   */
  unsigned short GetnElasticityMat(void) const { return nElasticityMod; }

  /*!
   * \brief Get the integration coefficients for the Generalized Alpha - Newmark integration integration scheme.
   * \param[in] val_coeff - Index of the coefficient.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double Get_Int_Coeffs(unsigned short val_coeff) const { return Int_Coeffs[val_coeff]; }

  /*!
   * \brief Get the number of different values for the modulus of the electric field.
   * \return Number of different values for the modulus of the electric field.
   */
  unsigned short GetnElectric_Field(void) const { return nElectric_Field; }

  /*!
   * \brief Get the dimensionality of the electric field.
   * \return Number of integration coefficients.
   */
  unsigned short GetnDim_Electric_Field(void) const { return nDim_Electric_Field; }

  /*!
   * \brief Get the values for the electric field modulus.
   * \param[in] val_coeff - Index of the coefficient.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double Get_Electric_Field_Mod(unsigned short val_coeff) const { return Electric_Field_Mod[val_coeff]; }

  /*!
   * \brief Set the values for the electric field modulus.
   * \param[in] val_coeff - Index of the electric field.
   * \param[in] val_el_field - Value of the electric field.
   */
  void Set_Electric_Field_Mod(unsigned short val_coeff, su2double val_el_field) { Electric_Field_Mod[val_coeff] = val_el_field; }

  /*!
   * \brief Get the direction of the electric field in reference configuration.
   * \param[in] val_coeff - Index of the coefficient.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  const su2double* Get_Electric_Field_Dir(void) const { return Electric_Field_Dir; }

  /*!
   * \brief Check if the user wants to apply the load as a ramp.
   * \return    <code>TRUE</code> means that the load is to be applied as a ramp.
   */
  bool GetRamp_Load(void) const { return Ramp_Load; }

  /*!
   * \brief Get the maximum time of the ramp.
   * \return    Value of the max time while the load is linearly increased
   */
  su2double GetRamp_Time(void) const { return Ramp_Time; }

  /*!
   * \brief Check if the user wants to apply the load as a ramp.
   * \return  <code>TRUE</code> means that the load is to be applied as a ramp.
   */
  bool GetRampAndRelease_Load(void) const { return RampAndRelease; }

  /*!
   * \brief Check if the user wants to apply the load as a ramp.
   * \return  <code>TRUE</code> means that the load is to be applied as a ramp.
   */
  bool GetSine_Load(void) const { return Sine_Load; }

  /*!
   * \brief Get the sine load properties.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The pointer to the sine load values.
   */
  const su2double* GetLoad_Sine(void) const { return sineload_coeff; }

  /*!
   * \brief Get the kind of load transfer method we want to use for dynamic problems
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of transfer method for multiphysics problems
   */
  unsigned short GetDynamic_LoadTransfer(void) const { return Dynamic_LoadTransfer; }

  /*!
   * \brief Get the penalty weight value for the objective function.
   * \return  Penalty weight value for the reference geometry objective function.
   */
  su2double GetRefGeom_Penalty(void) const { return RefGeom_Penalty; }

  /*!
   * \brief Get the penalty weight value for the objective function.
   * \return  Penalty weight value for the reference geometry objective function.
   */
  su2double GetTotalDV_Penalty(void) const { return DV_Penalty; }

  /*!
   * \brief Get the maximum allowed VM stress and KS exponent for the stress penalty objective function.
   */
  array<su2double,2> GetStressPenaltyParam(void) const { return StressPenaltyParam; }

  /*!
   * \brief Get whether a predictor is used for FSI applications.
   * \return Bool: determines if predictor is used or not
   */
  bool GetPredictor(void) const { return Predictor; }

  /*!
   * \brief Get the order of the predictor for FSI applications.
   * \return Order of predictor
   */
  unsigned short GetPredictorOrder(void) const { return Pred_Order; }

  /*!
   * \brief Get boolean for using Persson's shock capturing method in Euler flow DG-FEM
   * \return Boolean for using Persson's shock capturing method in Euler flow DG-FEM
   */
  bool GetEulerPersson(void) const { return EulerPersson; }

  /*!
   * \brief Set boolean for using Persson's shock capturing method in Euler flow DG-FEM
   * \param[in] val_EulerPersson - Boolean for using Persson's shock capturing method in Euler flow DG-FEM
   */
  void SetEulerPersson(bool val_EulerPersson) { EulerPersson = val_EulerPersson; }

  /*!
   * \brief Get whether a relaxation parameter is used for FSI applications.
   * \return Bool: determines if relaxation parameter  is used or not
   */
  bool GetRelaxation(void) const { return Relaxation; }

  /*!
   * \brief Check if the simulation we are running is a FSI simulation
   * \return Value of the physical time in an unsteady simulation.
   */
  bool GetFSI_Simulation(void) const { return FSI_Problem || (nMarker_Fluid_Load > 0); }

  /*!
   * \brief Set that the simulation we are running is a multizone simulation
   * \param[in] MZ_problem - boolean that determines is Multizone_Problem is true/false.
   */
  void SetMultizone_Problem(bool MZ_problem) { Multizone_Problem = MZ_problem; }

  /*!
   * \brief Get whether the simulation we are running is a multizone simulation
   * \return Multizone_Problem - boolean that determines is Multizone_Problem is true/false.
   */
  bool GetMultizone_Problem(void) const { return Multizone_Problem; }

  /*!
   * \brief Get the ID for the FEA region that we want to compute the gradient for using direct differentiation
   * \return ID
   */
  unsigned short GetnID_DV(void) const { return nID_DV; }

  /*!
   * \brief Check if we want to apply an incremental load to the nonlinear structural simulation
   * \return <code>TRUE</code> means that the load is to be applied in increments.
   */
  bool GetIncrementalLoad(void) const { return IncrementalLoad; }

  /*!
   * \brief Get the number of increments for an incremental load.
   * \return Number of increments.
   */
  unsigned long GetNumberIncrements(void) const { return IncLoad_Nincrements; }

  /*!
   * \brief Get the value of the criteria for applying incremental loading.
   * \return Value of the log10 of the residual.
   */
  su2double GetIncLoad_Criteria(unsigned short val_var) const { return inc_crit[val_var]; }

  /*!
   * \brief Get the relaxation method chosen for the simulation
   * \return Value of the relaxation method
   */
  BGS_RELAXATION GetRelaxation_Method_BGS(void) const { return Kind_BGS_RelaxMethod; }

  /*!
   * \brief Get the kind of Riemann solver for the DG method (FEM flow solver).
   * \note This value is obtained from the config file, and it is constant during the computation.
   * \return Kind of Riemann solver for the DG method (FEM flow solver).
   */
  UPWIND GetRiemann_Solver_FEM(void) const { return Riemann_Solver_FEM; }

  /*!
   * \brief Get the factor applied during quadrature of straight elements.
   * \return The specified straight element quadrature factor.
   */
  su2double GetQuadrature_Factor_Straight(void) const { return Quadrature_Factor_Straight; }

  /*!
   * \brief Get the factor applied during quadrature of curved elements.
   * \return The specified curved element quadrature factor.
   */
  su2double GetQuadrature_Factor_Curved(void) const { return Quadrature_Factor_Curved; }

  /*!
   * \brief Get the factor applied during time quadrature for ADER-DG.
   * \return The specified ADER-DG time quadrature factor.
   */
  su2double GetQuadrature_Factor_Time_ADER_DG(void) const { return Quadrature_Factor_Time_ADER_DG; }

  /*!
   * \brief Function to make available the multiplication factor theta of the
   *        symmetrizing terms in the DG discretization of the viscous terms.
   * \return The specified factor for the DG discretization.
   */
  su2double GetTheta_Interior_Penalty_DGFEM(void) const { return Theta_Interior_Penalty_DGFEM; }

  /*!
   * \brief Function to make available the matrix size in vectorization in
            order to optimize the gemm performance.
   * \return The matrix size in this direction.
   */
  unsigned short GetSizeMatMulPadding(void) const { return sizeMatMulPadding; }

  /*!
   * \brief Function to make available whether or not the entropy must be computed.
   * \return The boolean whether or not the entropy must be computed.
   */
  bool GetCompute_Entropy(void) const { return Compute_Entropy; }

  /*!
   * \brief Function to make available whether or not the lumped mass matrix
            must be used for steady computations.
   * \return The boolean whether or not to use the lumped mass matrix.
   */
  bool GetUse_Lumped_MassMatrix_DGFEM(void) const { return Use_Lumped_MassMatrix_DGFEM; }

  /*!
   * \brief Function to make available whether or not only the exact Jacobian
   *        of the spatial discretization must be computed.
   * \return The boolean whether or not the Jacobian must be computed.
   */
  bool GetJacobian_Spatial_Discretization_Only(void) const { return Jacobian_Spatial_Discretization_Only; }

  /*!
   * \brief Get the interpolation method used for matching between zones.
   */
  INTERFACE_INTERPOLATOR GetKindInterpolation(void) const { return Kind_Interpolation; }

  /*!
   * \brief Get option of whether to use conservative interpolation between zones.
   */
  bool GetConservativeInterpolation(void) const { return ConservativeInterpolation && GetStructuralProblem(); }

  /*!
   * \brief Get the basis function to use for radial basis function interpolation for FSI.
   */
  RADIAL_BASIS GetKindRadialBasisFunction(void) const { return Kind_RadialBasisFunction; }

  /*!
   * \brief Get option of whether to use polynomial terms in Radial Basis Function interpolation.
   */
  bool GetRadialBasisFunctionPolynomialOption(void) const { return RadialBasisFunction_PolynomialOption; }

  /*!
   * \brief Get the basis function radius to use for radial basis function interpolation for FSI.
   */
  su2double GetRadialBasisFunctionParameter(void) const { return RadialBasisFunction_Parameter; }

  /*!
   * \brief Get the tolerance used to prune the interpolation matrix (making it sparser).
   */
  su2double GetRadialBasisFunctionPruneTol(void) const { return RadialBasisFunction_PruneTol; }

  /*!
   * \brief Get the number of donor points to use in Nearest Neighbor interpolation.
   */
  unsigned short GetNumNearestNeighbors(void) const { return NumNearestNeighbors; }

  /*!
   * \brief Get the kind of inlet face interpolation function to use.
   */
  inline INLET_SPANWISE_INTERP GetKindInletInterpolationFunction(void) const { return Kind_InletInterpolationFunction; }

  /*!
   * \brief Get the kind of inlet face interpolation data type.
   */
  inline INLET_INTERP_TYPE GetKindInletInterpolationType (void) const  { return Kind_Inlet_InterpolationType; }

  /*!
   * \brief Get whether to print inlet interpolated data or not.
   */
  bool GetPrintInlet_InterpolatedData(void) const { return PrintInlet_InterpolatedData; }

  /*!
   * \brief Get the amount of eigenvalue perturbation to be done
   * \return Value of the uq_delta_b parameter
   */
  su2double GetUQ_Delta_B(void) const { return uq_delta_b; }

  /*!
   * \brief Get the kind of eigenspace perturbation to be done
   * \return Value of the eig_val_comp
   */
  unsigned short GetEig_Val_Comp(void) const { return eig_val_comp; }

  /*!
   * \brief Get the underelaxation factor
   * \return Value of the uq_urlx parameter
   */
  su2double GetUQ_URLX(void) const { return uq_urlx; }

  /*!
   * \brief Get information about eigenspace perturbation
   * \return <code>TRUE</code> means eigenspace perterturbation will be used
   */
  bool GetUQ_Permute(void) const { return uq_permute; }

  /*!
   * \brief Get information about whether to use wall functions.
   * \return <code>TRUE</code> if wall functions are on; otherwise <code>FALSE</code>.
   */
  bool GetWall_Functions(void) const { return Wall_Functions; }

  /*!
   * \brief Get the AD support.
   */
  bool GetAD_Mode(void) const { return AD_Mode;}

  /*!
   * \brief Set the maximum velocity^2 in the domain for the incompressible preconditioner.
   * \param[in] Value of the maximum velocity^2 in the domain for the incompressible preconditioner.
   */
  void SetMax_Vel2(su2double val_max_vel2) { Max_Vel2 = val_max_vel2; }

  /*!
   * \brief Get the maximum velocity^2 in the domain for the incompressible preconditioner.
   * \return Value of the maximum velocity^2 in the domain for the incompressible preconditioner.
   */
  su2double GetMax_Vel2(void) const { return Max_Vel2; }

  /*!
   * \brief Set the sum of the bandwidth for writing binary restarts (to be averaged later).
   * \param[in] Sum of the bandwidth for writing binary restarts.
   */
  void SetRestart_Bandwidth_Agg(su2double val_restart_bandwidth_sum) { Restart_Bandwidth_Agg = val_restart_bandwidth_sum; }

  /*!
   * \brief Set the sum of the bandwidth for writing binary restarts (to be averaged later).
   * \return Sum of the bandwidth for writing binary restarts.
   */
  su2double GetRestart_Bandwidth_Agg(void) const { return Restart_Bandwidth_Agg; }

  /*!
   * \brief Get the Kind of Hybrid RANS/LES.
   * \return Value of Hybrid RANS/LES method.
   */
  unsigned short GetKind_HybridRANSLES(void) const { return Kind_HybridRANSLES; }

  /*!
   * \brief Get the Kind of Roe Low Dissipation Scheme for Unsteady flows.
   * \return Value of Low dissipation approach.
   */
  unsigned short GetKind_RoeLowDiss(void) const { return Kind_RoeLowDiss; }

  /*!
   * \brief Get the DES Constant.
   * \return Value of DES constant.
   */
  su2double GetConst_DES(void) const { return Const_DES; }

  /*!
   * \brief Get if AD preaccumulation should be performed.
   */
  bool GetAD_Preaccumulation(void) const { return AD_Preaccumulation;}

  /*!
   * \brief Get the heat equation.
   * \return YES if weakly coupled heat equation for inc. flow is enabled.
   */
  bool GetWeakly_Coupled_Heat(void) const { return Weakly_Coupled_Heat; }

  /*!
   * \brief Get the CHT couling method.
   * \return Kind of the method.
   */
  CHT_COUPLING GetKind_CHT_Coupling() const { return Kind_CHT_Coupling; }

  /*!
   * \brief Check if values passed to the BC_HeatFlux-Routine are already integrated.
   * \return YES if the passed values is the integrated heat flux over the marker's surface.
   */
  bool GetIntegrated_HeatFlux() const { return Integrated_HeatFlux; }

  /*!
   * \brief Get Compute Average.
   * \return YES if start computing averages
   */
  bool GetCompute_Average(void) const { return Compute_Average;}

  /*!
   * \brief Get the verification solution.
   * \return The verification solution to be used.
   */
  VERIFICATION_SOLUTION GetVerification_Solution(void) const { return Kind_Verification_Solution;}

  /*!
   * \brief Get topology optimization.
   */
  bool GetTopology_Optimization(void) const { return topology_optimization; }

  /*!
   * \brief Get name of output file for topology optimization derivatives.
   */
  string GetTopology_Optim_FileName(void) const { return top_optim_output_file; }

  /*!
   * \brief Get exponent for density-based stiffness penalization.
   */
  su2double GetSIMP_Exponent(void) const { return simp_exponent; }

  /*!
   * \brief Get lower bound for density-based stiffness penalization.
   */
  su2double GetSIMP_MinStiffness(void) const { return simp_minimum_stiffness; }

  /*!
   * \brief Number of kernels to use in filtering the design density field.
   */
  unsigned short GetTopology_Optim_Num_Kernels(void) const { return top_optim_nKernel; }

  /*!
   * \brief Get the i'th kernel to use, its parameter, and the radius.
   */
  void GetTopology_Optim_Kernel(const unsigned short iKernel, ENUM_FILTER_KERNEL &type,
                                su2double &param, su2double &radius) const {
    type = top_optim_kernels[iKernel];
    param = top_optim_kernel_params[iKernel];
    radius = top_optim_filter_radius[iKernel];
  }

  /*!
   * \brief Get the maximum "logical radius" (degree of neighborhood) to consider in the neighbor search.
   */
  unsigned short GetTopology_Search_Limit(void) const { return top_optim_search_lim; }

  /*!
   * \brief Get the type and parameter for the projection function used in topology optimization
   */
  void GetTopology_Optim_Projection(ENUM_PROJECTION_FUNCTION &type, su2double &param) const {
    type = top_optim_proj_type;  param = top_optim_proj_param;
  }

  /*!
   * \brief Get the filenames of the individual config files
   * \return File name of the config file for zone "index"
   */
  string GetConfigFilename(unsigned short index) const { return Config_Filenames[index]; }

  /*!
   * \brief Get the number of config files
   * \return Number of config filenames in CONFIG_LIST
   */
  unsigned short GetnConfigFiles(void) const { return nConfig_Files; }

  /*!
   * \brief Check if the multizone problem is solved for time domain.
   * \return YES if time-domain is considered.
   */
  bool GetTime_Domain(void) const { return Time_Domain; }

  /*!
   * \brief Get the number of inner iterations
   * \return Number of inner iterations on each multizone block
   */
  unsigned long GetnInner_Iter(void) const { return nInnerIter; }

  /*!
   * \brief Get the number of outer iterations
   * \return Number of outer iterations for the multizone problem
   */
  unsigned long GetnOuter_Iter(void) const { return nOuterIter; }

  /*!
   * \brief Get the number of time iterations
   * \return Number of time steps run
   */
  unsigned long GetnTime_Iter(void) const { return nTimeIter; }

  /*!
   * \brief Set the number of time iterations
   * \param[in] val_iter - Number of time steps run
   */
  void SetnTime_Iter(unsigned long val_iter) { nTimeIter = val_iter; }

  /*!
   * \brief Get the restart iteration
   * \return Iteration for the restart of multizone problems
   */
  unsigned long GetRestart_Iter(void) const { return Restart_Iter; }

  /*!
   * \brief Get the time step for multizone problems
   * \return Time step for multizone problems, it is set on all the zones
   */
  su2double GetTime_Step(void) const { return Time_Step; }

  /*!
   * \brief Get the maximum simulation time for time-domain problems
   * \return Simulation time for multizone problems, it is set on all the zones
   */
  su2double GetMax_Time(void) const { return Max_Time; }

  /*!
   * \brief Get the level of MPI communications to be performed.
   * \return Level of MPI communications.
   */
  unsigned short GetComm_Level(void) const { return Comm_Level; }

  /*!
   * \brief Check if the mesh read supports multiple zones.
   * \return YES if multiple zones can be contained in the mesh file.
   */
  bool GetMultizone_Mesh(void) const { return Multizone_Mesh; }

  /*!
   * \brief Check if the mesh read supports multiple zones.
   * \return YES if multiple zones can be contained in the mesh file.
   */
  bool GetMultizone_Residual(void) const { return Multizone_Residual; }

  /*!
   * \brief Get the Kind of Radiation model applied.
   * \return Kind of radiation model used.
   */
  RADIATION_MODEL GetKind_RadiationModel(void) const { return Kind_Radiation; }

  /*!
   * \brief Get the Kind of P1 initialization method applied.
   * \return Kind of P1 initialization method used.
   */
  P1_INIT GetKind_P1_Init(void) const { return Kind_P1_Init; }

  /*!
   * \brief Get the value of the absorption coefficient of the medium.
   * \return Value of the absorption coefficient of the medium.
   */
  su2double GetAbsorption_Coeff(void) const { return Absorption_Coeff; }

  /*!
   * \brief Get the value of the scattering coefficient of the medium.
   * \return Value of the scattering coefficient of the medium.
   */
  su2double GetScattering_Coeff(void) const { return Scattering_Coeff; }

  /*!
   * \brief Get the wall emissivity at a boundary.
   * \param[in] val_index - Index corresponding to the boundary.
   * \return The wall emissivity.
   */
  su2double GetWall_Emissivity(const string& val_index) const;

  /*!
   * \brief Get if boundary is strong or weak.
   * \param[in] val_index - Index corresponding to the boundary.
   * \return true if strong BC.
   */
  bool GetMarker_StrongBC(const string& val_index) const;

  /*!
   * \brief Get the value of the CFL condition for radiation solvers.
   * \return Value of the CFL condition for radiation solvers.
   */
  su2double GetCFL_Rad(void) const { return CFL_Rad; }

  /*!
   * \brief Determines if radiation needs to be incorporated to the analysis.
   * \return Radiation boolean
   */
  bool AddRadiation(void) const { return Radiation; }

  /*!
   * \brief Check if the convergence history of each individual zone is written to screen
   * \return YES if the zone convergence history of each individual zone must be written to screen
   */
  bool GetWrt_ZoneConv(void) const { return Wrt_ZoneConv; }

  /*!
   * \brief Check if the convergence history of each individual zone is written to file
   * \return YES if the zone convergence history of each individual zone must be written to file
   */
  bool GetWrt_ZoneHist(void) const { return Wrt_ZoneHist; }

  /*!
   * \brief Check if the special output is written
   * \return YES if the special output is written.
   */
  bool GetSpecial_Output(void) const { return SpecialOutput; }

  /*!
   * \brief Check if the forces breakdown file is written
   * \return YES if the forces breakdown file is written.
   */
  bool GetWrt_ForcesBreakdown(void) const { return Wrt_ForcesBreakdown; }

  /*!
   * \brief Get the number of grid points in the analytic RECTANGLE or BOX grid in the specified coordinate direction.
   * \return Number of grid points in the analytic RECTANGLE or BOX grid in the specified coordinate direction.
   */
  short GetMeshBoxSize(unsigned short val_iDim) const { return Mesh_Box_Size[val_iDim]; }

  /*!
   * \brief Get the length of the analytic RECTANGLE or BOX grid in the specified coordinate direction.
   * \return Length the analytic RECTANGLE or BOX grid in the specified coordinate direction.
   */
  su2double GetMeshBoxLength(unsigned short val_iDim) const { return mesh_box_length[val_iDim]; }

  /*!
   * \brief Get the offset from 0.0 of the analytic RECTANGLE or BOX grid in the specified coordinate direction.
   * \return Offset from 0.0 the analytic RECTANGLE or BOX grid in the specified coordinate direction.
   */
  su2double GetMeshBoxOffset(unsigned short val_iDim) const { return mesh_box_offset[val_iDim]; }

  /*!
   * \brief Get the number of screen output variables requested (maximum 6)
   */
  unsigned short GetnScreenOutput(void) const { return nScreenOutput; }

  /*!
   * \brief Get the screen output field iField
   */
  string GetScreenOutput_Field(unsigned short iField) const { return ScreenOutput[iField]; }

  /*!
   * \brief Get the number of history output variables requested
   */
  unsigned short GetnHistoryOutput(void) const { return nHistoryOutput; }

  /*!
   * \brief Get the history output field iField
   */
  string GetHistoryOutput_Field(unsigned short iField) const { return HistoryOutput[iField]; }

  /*!
   * \brief Get the number of history output variables requested
   */
  unsigned short GetnVolumeOutput(void) const { return nVolumeOutput; }

  /*!
   * \brief Get the history output field iField
   */
  string GetVolumeOutput_Field(unsigned short iField) const { return VolumeOutput[iField]; }

  /*!
  * \brief Get the convergence fields for monitoring
  * \param[in] iField - Index of the field
  * return Field name for monitoring convergence
  */
  string GetConv_Field(unsigned short iField) const { return ConvField[iField]; }

  /*!
   * \brief Get functional that is going to be used to evaluate the convergence of the windowed time average of the unsteady problem.
   * \param[in] iField - Index of the field
   * \return Field name for monitoring convergence
   */
  string GetWndConv_Field(unsigned short iField) const { return WndConvField[iField]; }

  /*!
   * \brief Get the number of iterations that are considered in the Cauchy convergence criteria for the windowed time average of the unsteady problem.
   * \return Number of elements in the Cauchy criteria windowed time average of the unsteady problem.
   */
  unsigned short GetWnd_Cauchy_Elems(void) const { return Wnd_Cauchy_Elems; }

  /*!
   * \brief Get the value of convergence criteria for the Cauchy method for the time averaged
   *        windowed objective functions for unsteady flows
   * \return Value of the convergence criteria.
   */
  su2double GetWnd_Cauchy_Eps(void) const { return Wnd_Cauchy_Eps; }

  /*!
   * \brief Get the number of iterations that are not considered in the convergence criteria for the windowed average output function
   * \return Number of iterations before starting with the convergence criteria for the windowed average output function.
   */
  unsigned long  GetWnd_StartConv_Iter(void) const { return Wnd_StartConv_Iter; }

  /*!
   * \brief Get the boolean value, whether the the Cauchy method for the time averaged
   *        windowed objective functions for unsteady flows is used or not.
   * \return Boolean value, if the criterion is used.
   */
  bool GetWnd_Cauchy_Crit(void) const { return Wnd_Cauchy_Crit; }

  /*!
  * \brief Get the number of convergence monitoring fields for time convergence monitoring.
  * return Number of convergence monitoring fields.
  */
  unsigned short GetnWndConv_Field() const { return nWndConvField; }

  /*!
  * \brief Get the number of convergence monitoring fields for inner convergence monitoring.
  * return Number of convergence monitoring fields.
  */
  unsigned short GetnConv_Field() const { return nConvField; }

  /*!
   * \brief Set the start time to track a phase of the code (preprocessing, compute, output).
   * \param[in] Value of the start time to track a phase of the code.
   */
  void Set_StartTime(su2double starttime) { StartTime = starttime; }

  /*!
   * \brief Get the start time to track a phase of the code (preprocessing, compute, output).
   * \return Value of the start time to track a phase of the code.
   */
  su2double Get_StartTime() const { return StartTime; }

  /*!
   * \brief GetHistory_Wrt_Freq_Inner
   * \return
   */
  unsigned long GetHistory_Wrt_Freq(unsigned short iter) const { return HistoryWrtFreq[iter];}

  /*!
   * \brief SetHistory_Wrt_Freq_Inner
   * \param[in] iter: index for Time (0), Outer (1), or Inner (2) iterations
   * \param[in] nIter: Number of iterations
   */
  void SetHistory_Wrt_Freq(unsigned short iter, unsigned long nIter) { HistoryWrtFreq[iter] = nIter;}

  /*!
   * \brief GetScreen_Wrt_Freq_Inner
   * \param[in] iter: index for Time (0), Outer (1), or Inner (2) iterations
   * \return
   */
  unsigned long GetScreen_Wrt_Freq(unsigned short iter) const { return ScreenWrtFreq[iter]; }

  /*!
   * \brief SetScreen_Wrt_Freq_Inner
   * \param[in] iter: index for Time (0), Outer (1), or Inner (2) iterations
   * \param[in] nIter: Number of iterations
   */
  void SetScreen_Wrt_Freq(unsigned short iter, unsigned long nIter) { ScreenWrtFreq[iter] = nIter; }

  /*!
   * \brief GetVolumeOutputFiles
   */
  const OUTPUT_TYPE* GetVolumeOutputFiles() const { return VolumeOutputFiles; }

  /*!
   * \brief GetnVolumeOutputFiles
   */
  unsigned short GetnVolumeOutputFiles() const { return nVolumeOutputFiles; }

  /*!
   * \brief GetVolumeOutputFrequency
   * \param[in] iFile: index of file number for which the writing frequency needs to be returned.
   */
  unsigned long GetVolumeOutputFrequency(unsigned short iFile) const { return VolumeOutputFrequencies[iFile]; }

  /*!
   * \brief Get the desired factorization frequency for PaStiX
   * \return Number of calls to 'Build' that trigger re-factorization.
   */
  unsigned long GetPastixFactFreq(void) const { return pastix_fact_freq; }

  /*!
   * \brief Get the desired level of verbosity for PaStiX
   * \return 0 - Quiet, 1 - During factorization and cleanup, 2 - Even more detail.
   */
  unsigned short GetPastixVerbLvl(void) const { return pastix_verb_lvl; }

  /*!
   * \brief Get the desired level of fill for the PaStiX ILU
   * \return Level of fill.
   */
  unsigned short GetPastixFillLvl(void) const { return pastix_fill_lvl; }

  /*!
   * \brief Check if an option is present in the config file
   * \param[in] - Name of the option
   * \return <TRUE> if option was set in the config file
   */
  bool OptionIsSet(string option) const { return all_options.find(option) == all_options.end(); }

  /*!
   * \brief Get the name of the current case
   * \return the case name
   */
  const string& GetCaseName() const { return caseName; }

  /*!
   * \brief Get the number of threads per rank to use for ILU and LU_SGS preconditioners.
   * \return Number of threads per rank.
   */
  unsigned long GetLinear_Solver_Prec_Threads(void) const { return Linear_Solver_Prec_Threads; }

  /*!
   * \brief Get the size of the edge groups colored for OpenMP parallelization of edge loops.
   */
  unsigned long GetEdgeColoringGroupSize(void) const { return edgeColorGroupSize; }

  /*!
   * \brief Get the ParMETIS load balancing tolerance.
   */
  passivedouble GetParMETIS_Tolerance() const { return SU2_TYPE::GetValue(ParMETIS_tolerance); }

  /*!
   * \brief Get the ParMETIS load balancing weight for points.
   */
  long GetParMETIS_PointWeight() const { return ParMETIS_pointWgt; }

  /*!
   * \brief Get the ParMETIS load balancing weight for edges
   */
  long GetParMETIS_EdgeWeight() const { return ParMETIS_edgeWgt; }

  /*!
   * \brief Find the marker index (if any) that is part of a given interface pair.
   * \param[in] iInterface - Number of the interface pair being tested, starting at 0.
   * \return -1 if (on this mpi rank) the zone defined by config is not part of the interface.
   */
  short FindInterfaceMarker(unsigned short iInterface) const;

  /*!
   * \brief Get whether or not to save solution data to libROM.
   * \return True if specified in config file.
   */
  bool GetSave_libROM(void) const {return libROM; }

  /*!
   * \brief Get the name of the file for libROM to save.
   * \return Filename prefix for libROM to save to (default: "su2").
   */
  string GetlibROMbase_FileName(void) const { return libROMbase_FileName; }

  /*!
   * \brief Static or incremental toggle for POD basis generation type.
   * \return Type of POD generation type
   */
  POD_KIND GetKind_PODBasis(void) const { return POD_Basis_Gen; }

  /*!
   * \brief Get maximum number of POD basis dimensions (default: 100).
   * \return Maximum number of POD basis vectors.
   */
  unsigned short GetMax_BasisDim(void) const { return maxBasisDim; }

  /*!
   * \brief Get frequency of unsteady time steps to save (default: 1).
   * \return Save frequency for unsteady time steps.
   */
  unsigned short GetRom_SaveFreq(void) const { return rom_save_freq; }

  /*!
   * \brief Check if the gradient smoothing is active
   * \return true means that smoothing is applied to the sensitivities
   */
  bool GetSmoothGradient(void) const {return SmoothGradient; }

  /*!
   * \brief Gets the factor epsilon in front of the Laplace term
   * \return epsilon
   */
  su2double GetSmoothingEps1(void) const { return SmoothingEps1; }

  /*!
   * \brief Gets the factor zeta in front of the identity term
   * \return zeta
   */
  su2double GetSmoothingEps2(void) const { return SmoothingEps2; }

  /*!
   * \brief Check if we split in the dimensions
   * \return true means that smoothing is for each dimension separate
   */
  bool GetSmoothSepDim(void) const { return SmoothSepDim; }

  /*!
   * \brief Check if we assemble the operator on the surface
   * \return true means that smoothing is done on the surface level
   */
  bool GetSmoothOnSurface(void) const { return SmoothOnSurface; }

  /*!
   * \brief Check if we use zero Dirichlet boundarys on the bound of the surface
   * \return true means that we use zero Dirichlet boundary
   */
  bool GetDirichletSurfaceBound(void) const { return SmoothDirichletSurfaceBound; }

  /*!
   * \brief The modus of operation for the Sobolev solver
   * \return returns on what level we operate
   */
  ENUM_SOBOLEV_MODUS GetSobMode(void) const { return SmoothNumMode; }

  /*!
   * \brief Get the name of the file with the hessian of the objective function.
   * \return Name of the file with the hessian of the objective function.
   */
  string GetObjFunc_Hess_FileName(void) const { return ObjFunc_Hess_FileName; }

  /*!
   * \brief Get min error of the linear solver for the gradient smoothing.
   * \return Min error of the linear solver for the gradient smoothing.
   */
  su2double GetGrad_Linear_Solver_Error(void) const { return Grad_Linear_Solver_Error; }

  /*!
   * \brief Get the kind of solver for the gradient smoothing.
   * \return Numerical solver for the gradient smoothing.
   */
  unsigned short GetKind_Grad_Linear_Solver(void) const { return Kind_Grad_Linear_Solver; }

  /*!
   * \brief Get the kind of preconditioner for the gradient smoothing.
   * \return Numerical preconditioner for the gradient smoothing.
   */
  unsigned short GetKind_Grad_Linear_Solver_Prec(void) const { return Kind_Grad_Linear_Solver_Prec; }

  /*!
   * \brief Get max number of iterations of the for the gradient smoothing.
   * \return Max number of iterations of the linear solver for the gradient smoothing.
   */
  unsigned long GetGrad_Linear_Solver_Iter(void) const { return Grad_Linear_Solver_Iter; }

  /*!
   * \brief Get parsed SST option data structure.
   * \return SST option data structure.
   */
  SST_ParsedOptions GetSSTParsedOptions() const { return sstParsedOptions; }

  /*!
   * \brief Get parsed SA option data structure.
   * \return SA option data structure.
   */
  SA_ParsedOptions GetSAParsedOptions() const { return saParsedOptions; }

  /*!
   * \brief Get parsed LM option data structure.
   * \return LM option data structure.
   */
  LM_ParsedOptions GetLMParsedOptions() const { return lmParsedOptions; }

};
