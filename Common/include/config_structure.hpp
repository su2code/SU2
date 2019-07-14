/*!
 * \file config_structure.hpp
 * \brief All the information about the definition of the physical problem.
 *        The subroutines and functions are in the <i>config_structure.cpp</i> file.
 * \author F. Palacios, T. Economon, B. Tracey
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "./mpi_structure.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <assert.h>

#include "./option_structure.hpp"
#include "./datatype_structure.hpp"

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
  int rank, size;
  unsigned short Kind_SU2; /*!< \brief Kind of SU2 software component.*/
  unsigned short Ref_NonDim; /*!< \brief Kind of non dimensionalization.*/
  unsigned short Ref_Inc_NonDim; /*!< \brief Kind of non dimensionalization.*/
  unsigned short Kind_AverageProcess; /*!< \brief Kind of mixing process.*/
  unsigned short Kind_PerformanceAverageProcess; /*!< \brief Kind of mixing process.*/
  unsigned short Kind_MixingPlaneInterface; /*!< \brief Kind of mixing process.*/
  unsigned short Kind_SpanWise; /*!< \brief Kind of span-wise section computation.*/
  unsigned short *Kind_TurboMachinery;  /*!< \brief Kind of turbomachynery architecture.*/
  unsigned short iZone, nZone; /*!< \brief Number of zones in the mesh. */
  unsigned short nZoneSpecified; /*!< \brief Number of zones that are specified in config file. */
  su2double Highlite_Area; /*!< \brief Highlite area. */
  su2double Fan_Poly_Eff; /*!< \brief Highlite area. */
  su2double OrderMagResidual; /*!< \brief Order of magnitude reduction. */
  su2double MinLogResidual; /*!< \brief Minimum value of the log residual. */
  su2double OrderMagResidualFSI; /*!< \brief Order of magnitude reduction. */
  su2double MinLogResidualFSI; /*!< \brief Minimum value of the log residual. */
  su2double OrderMagResidual_BGS_F; /*!< \brief Order of magnitude reduction. */
  su2double MinLogResidual_BGS_F; /*!< \brief Minimum value of the log residual. */
  su2double OrderMagResidual_BGS_S; /*!< \brief Order of magnitude reduction. */
  su2double MinLogResidual_BGS_S; /*!< \brief Minimum value of the log residual. */
  su2double Res_FEM_UTOL; 		/*!< \brief UTOL criteria for structural FEM. */
  su2double Res_FEM_RTOL; 		/*!< \brief RTOL criteria for structural FEM. */
  su2double Res_FEM_ETOL; 		/*!< \brief ETOL criteria for structural FEM. */
  su2double Res_FEM_ADJ;     /*!< \brief Convergence criteria for adjoint FEM. */
  su2double EA_ScaleFactor; /*!< \brief Equivalent Area scaling factor */
  su2double* EA_IntLimit; /*!< \brief Integration limits of the Equivalent Area computation */
  su2double AdjointLimit; /*!< \brief Adjoint variable limit */
  su2double* Obj_ChainRuleCoeff; /*!< \brief Array defining objective function for adjoint problem based on chain rule in terms of gradient w.r.t. density, velocity, pressure */
  bool MG_AdjointFlow; /*!< \brief MG with the adjoint flow problem */
  su2double* SubsonicEngine_Cyl; /*!< \brief Coordinates of the box subsonic region */
  su2double* SubsonicEngine_Values; /*!< \brief Values of the box subsonic region */
  su2double* Hold_GridFixed_Coord; /*!< \brief Coordinates of the box to hold fixed the nbumerical grid */
  su2double *DistortionRack;
  su2double *PressureLimits,
  *DensityLimits,
  *TemperatureLimits; /*!< \brief Limits for the primitive variables */
  bool ActDisk_DoubleSurface;  /*!< \brief actuator disk double surface  */
  bool Engine_HalfModel;  /*!< \brief only half model is in the computational grid  */
  bool ActDisk_SU2_DEF;  /*!< \brief actuator disk double surface  */
  unsigned short ConvCriteria;	/*!< \brief Kind of convergence criteria. */
  unsigned short nFFD_Iter; 	/*!< \brief Iteration for the point inversion problem. */
  unsigned short FFD_Blending; /*!< \brief Kind of FFD Blending function. */
  su2double* FFD_BSpline_Order; /*!< \brief BSpline order in i,j,k direction. */
  su2double FFD_Tol;  	/*!< \brief Tolerance in the point inversion problem. */
  su2double Opt_RelaxFactor;  	/*!< \brief Scale factor for the line search. */
  su2double Opt_LineSearch_Bound;  	/*!< \brief Bounds for the line search. */
  bool Write_Conv_FSI;			/*!< \brief Write convergence file for FSI problems. */
  bool ContinuousAdjoint,			/*!< \brief Flag to know if the code is solving an adjoint problem. */
  Viscous,                /*!< \brief Flag to know if the code is solving a viscous problem. */
  EquivArea,				/*!< \brief Flag to know if the code is going to compute and plot the equivalent area. */
  Engine,				/*!< \brief Flag to know if the code is going to compute a problem with engine. */
  InvDesign_Cp,				/*!< \brief Flag to know if the code is going to compute and plot the inverse design. */
  InvDesign_HeatFlux,				/*!< \brief Flag to know if the code is going to compute and plot the inverse design. */
  Grid_Movement,			/*!< \brief Flag to know if there is grid movement. */
  Wind_Gust,              /*!< \brief Flag to know if there is a wind gust. */
  Aeroelastic_Simulation, /*!< \brief Flag to know if there is an aeroelastic simulation. */
  Weakly_Coupled_Heat, /*!< \brief Flag to know if a heat equation should be weakly coupled to the incompressible solver. */
  Rotating_Frame,			/*!< \brief Flag to know if there is a rotating frame. */
  PoissonSolver,			/*!< \brief Flag to know if we are solving  poisson forces  in plasma solver. */
  Low_Mach_Precon,		/*!< \brief Flag to know if we are using a low Mach number preconditioner. */
  Low_Mach_Corr,			/*!< \brief Flag to know if we are using a low Mach number correction. */
  GravityForce,			/*!< \brief Flag to know if the gravity force is incuded in the formulation. */
  SmoothNumGrid,			/*!< \brief Smooth the numerical grid. */
  AdaptBoundary,			/*!< \brief Adapt the elements on the boundary. */
  SubsonicEngine,			/*!< \brief Engine intake subsonic region. */
  Frozen_Visc_Cont,			/*!< \brief Flag for cont. adjoint problem with/without frozen viscosity. */
  Frozen_Visc_Disc,			/*!< \brief Flag for disc. adjoint problem with/without frozen viscosity. */
  Frozen_Limiter_Disc,			/*!< \brief Flag for disc. adjoint problem with/without frozen limiter. */
  Inconsistent_Disc,      /*!< \brief Use an inconsistent (primal/dual) discrete adjoint formulation. */
  Sens_Remove_Sharp,			/*!< \brief Flag for removing or not the sharp edges from the sensitivity computation. */
  Hold_GridFixed,	/*!< \brief Flag hold fixed some part of the mesh during the deformation. */
  Axisymmetric, /*!< \brief Flag for axisymmetric calculations */
  Integrated_HeatFlux, /*!< \brief Flag for heat flux BC whether it deals with integrated values.*/
  Buffet_Monitoring;       /*!< \brief Flag for computing the buffet sensor.*/
  su2double Buffet_k;     /*!< \brief Sharpness coefficient for buffet sensor.*/
  su2double Buffet_lambda; /*!< \brief Offset parameter for buffet sensor.*/
  su2double Damp_Engine_Inflow;	/*!< \brief Damping factor for the engine inlet. */
  su2double Damp_Engine_Exhaust;	/*!< \brief Damping factor for the engine exhaust. */
  su2double Damp_Res_Restric,	/*!< \brief Damping factor for the residual restriction. */
  Damp_Correc_Prolong; /*!< \brief Damping factor for the correction prolongation. */
  su2double Position_Plane; /*!< \brief Position of the Near-Field (y coordinate 2D, and z coordinate 3D). */
  su2double WeightCd; /*!< \brief Weight of the drag coefficient. */
  su2double dCD_dCL; /*!< \brief Weight of the drag coefficient. */
  su2double dCMx_dCL; /*!< \brief Weight of the drag coefficient. */
  su2double dCMy_dCL; /*!< \brief Weight of the drag coefficient. */
  su2double dCMz_dCL; /*!< \brief Weight of the drag coefficient. */
  su2double dCD_dCMy; /*!< \brief Weight of the drag coefficient. */
  su2double CL_Target; /*!< \brief Weight of the drag coefficient. */
  su2double CM_Target; /*!< \brief Weight of the drag coefficient. */
  su2double *HTP_Min_XCoord, *HTP_Min_YCoord; /*!< \brief Identification of the HTP. */
  unsigned short Unsteady_Simulation;	/*!< \brief Steady or unsteady (time stepping or dual time stepping) computation. */
  unsigned short Dynamic_Analysis;	/*!< \brief Static or dynamic structural analysis. */
  unsigned short nStartUpIter;	/*!< \brief Start up iterations using the fine grid. */
  su2double FixAzimuthalLine; /*!< \brief Fix an azimuthal line due to misalignments of the nearfield. */
  su2double **DV_Value;		/*!< \brief Previous value of the design variable. */
  su2double Venkat_LimiterCoeff;				/*!< \brief Limiter coefficient */
  unsigned long LimiterIter;	/*!< \brief Freeze the value of the limiter after a number of iterations */
  su2double AdjSharp_LimiterCoeff;				/*!< \brief Coefficient to identify the limit of a sharp edge. */
  unsigned short SystemMeasurements; /*!< \brief System of measurements. */
  unsigned short Kind_Regime;  /*!< \brief Kind of adjoint function. */
  unsigned short *Kind_ObjFunc;  /*!< \brief Kind of objective function. */
  su2double *Weight_ObjFunc;    /*!< \brief Weight applied to objective function. */
  unsigned short Kind_SensSmooth; /*!< \brief Kind of sensitivity smoothing technique. */
  unsigned short Continuous_Eqns; /*!< \brief Which equations to treat continuously (Hybrid adjoint)*/
  unsigned short Discrete_Eqns; /*!< \brief Which equations to treat discretely (Hybrid adjoint). */
  unsigned short *Design_Variable; /*!< \brief Kind of design variable. */
  unsigned short Kind_Adaptation;	/*!< \brief Kind of numerical grid adaptation. */
  unsigned short nTimeInstances;  /*!< \brief Number of periodic time instances for  harmonic balance. */
  su2double HarmonicBalance_Period;		/*!< \brief Period of oscillation to be used with harmonic balance computations. */
  su2double New_Elem_Adapt;			/*!< \brief Elements to adapt in the numerical grid adaptation process. */
  su2double Delta_UnstTime,			/*!< \brief Time step for unsteady computations. */
  Delta_UnstTimeND;						/*!< \brief Time step for unsteady computations (non dimensional). */
  su2double Delta_DynTime,		/*!< \brief Time step for dynamic structural computations. */
  Total_DynTime,				/*!< \brief Total time for dynamic structural computations. */
  Current_DynTime;			/*!< \brief Global time of the dynamic structural computations. */
  su2double Total_UnstTime,						/*!< \brief Total time for unsteady computations. */
  Total_UnstTimeND;								/*!< \brief Total time for unsteady computations (non dimensional). */
  su2double Current_UnstTime,									/*!< \brief Global time of the unsteady simulation. */
  Current_UnstTimeND;									/*!< \brief Global time of the unsteady simulation. */
  unsigned short nMarker_Euler,	/*!< \brief Number of Euler wall markers. */
  nMarker_FarField,				/*!< \brief Number of far-field markers. */
  nMarker_Custom,
  nMarker_SymWall,				/*!< \brief Number of symmetry wall markers. */
  nMarker_PerBound,				/*!< \brief Number of periodic boundary markers. */
  nMarker_MixingPlaneInterface,				/*!< \brief Number of mixing plane interface boundary markers. */
  nMarker_Turbomachinery,				/*!< \brief Number turbomachinery markers. */
  nMarker_TurboPerformance,				/*!< \brief Number of turboperformance markers. */
  nSpanWiseSections_User,			/*!< \brief Number of spanwise sections to compute 3D BC and Performance for turbomachinery   */
  nMarker_Shroud,/*!< \brief Number of shroud markers to set grid velocity to 0.*/
  nMarker_NearFieldBound,				/*!< \brief Number of near field boundary markers. */
  nMarker_ActDiskInlet, nMarker_ActDiskOutlet,
  nMarker_InterfaceBound,				/*!< \brief Number of interface boundary markers. */
  nMarker_Fluid_InterfaceBound,				/*!< \brief Number of fluid interface markers. */
  nMarker_CHTInterface,     /*!< \brief Number of conjugate heat transfer interface markers. */
  nMarker_Dirichlet,				/*!< \brief Number of interface boundary markers. */
  nMarker_Inlet,					/*!< \brief Number of inlet flow markers. */
  nMarker_Riemann,					/*!< \brief Number of Riemann flow markers. */
  nMarker_Giles,					/*!< \brief Number of Giles flow markers. */
  nRelaxFactor_Giles,                                   /*!< \brief Number of relaxation factors for Giles markers. */
  nMarker_Supersonic_Inlet,					/*!< \brief Number of supersonic inlet flow markers. */
  nMarker_Supersonic_Outlet,					/*!< \brief Number of supersonic outlet flow markers. */
  nMarker_Outlet,					/*!< \brief Number of outlet flow markers. */
  nMarker_Isothermal,     /*!< \brief Number of isothermal wall boundaries. */
  nMarker_HeatFlux,       /*!< \brief Number of constant heat flux wall boundaries. */
  nMarker_EngineExhaust,					/*!< \brief Number of nacelle exhaust flow markers. */
  nMarker_EngineInflow,					/*!< \brief Number of nacelle inflow flow markers. */
  nMarker_Clamped,						/*!< \brief Number of clamped markers in the FEM. */
  nMarker_Displacement,					/*!< \brief Number of displacement surface markers. */
  nMarker_Load,					/*!< \brief Number of load surface markers. */
  nMarker_Damper,         /*!< \brief Number of damper surface markers. */
  nMarker_Load_Dir,					/*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_Disp_Dir,         /*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_Load_Sine,					/*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_FlowLoad,					/*!< \brief Number of load surface markers. */
  nMarker_Neumann,				/*!< \brief Number of Neumann flow markers. */
  nMarker_Internal,				/*!< \brief Number of Neumann flow markers. */
  nMarker_All,					/*!< \brief Total number of markers using the grid information. */
  nMarker_Max,					/*!< \brief Max number of number of markers using the grid information. */
  nMarker_CfgFile;					/*!< \brief Total number of markers using the config file
                             (note that using parallel computation this number can be different
                             from nMarker_All). */
  bool Inlet_From_File; /*!< \brief True if the inlet profile is to be loaded from a file. */
  string Inlet_Filename; /*!< \brief Filename specifying an inlet profile. */
  su2double Inlet_Matching_Tol; /*!< \brief Tolerance used when matching a point to a point from the inlet file. */
  string *Marker_Euler,			/*!< \brief Euler wall markers. */
  *Marker_FarField,				/*!< \brief Far field markers. */
  *Marker_Custom,
  *Marker_SymWall,				/*!< \brief Symmetry wall markers. */
  *Marker_PerBound,				/*!< \brief Periodic boundary markers. */
  *Marker_PerDonor,				/*!< \brief Rotationally periodic boundary donor markers. */
  *Marker_MixingPlaneInterface,				/*!< \brief MixingPlane interface boundary markers. */
  *Marker_TurboBoundIn,				/*!< \brief Turbomachinery performance boundary markers. */
  *Marker_TurboBoundOut,				/*!< \brief Turbomachinery performance boundary donor markers. */
  *Marker_NearFieldBound,				/*!< \brief Near Field boundaries markers. */
  *Marker_InterfaceBound,				/*!< \brief Interface boundaries markers. */
  *Marker_Fluid_InterfaceBound,				/*!< \brief Fluid interface markers. */
  *Marker_CHTInterface,         /*!< \brief Conjugate heat transfer interface markers. */
  *Marker_ActDiskInlet,
  *Marker_ActDiskOutlet,
  *Marker_Dirichlet,				/*!< \brief Interface boundaries markers. */
  *Marker_Inlet,					/*!< \brief Inlet flow markers. */
  *Marker_Riemann,					/*!< \brief Riemann markers. */
  *Marker_Giles,					/*!< \brief Giles markers. */
  *Marker_Shroud,                                       /*!< \brief Shroud markers. */
  *Marker_Supersonic_Inlet,					/*!< \brief Supersonic inlet flow markers. */
  *Marker_Supersonic_Outlet,					/*!< \brief Supersonic outlet flow markers. */
  *Marker_Outlet,					/*!< \brief Outlet flow markers. */
  *Marker_Isothermal,     /*!< \brief Isothermal wall markers. */
  *Marker_HeatFlux,       /*!< \brief Constant heat flux wall markers. */
  *Marker_EngineInflow,					/*!< \brief Engine Inflow flow markers. */
  *Marker_EngineExhaust,					/*!< \brief Engine Exhaust flow markers. */
  *Marker_Clamped,						/*!< \brief Clamped markers. */
  *Marker_Displacement,					/*!< \brief Displacement markers. */
  *Marker_Load,					/*!< \brief Load markers. */
  *Marker_Damper,         /*!< \brief Damper markers. */
  *Marker_Load_Dir,					/*!< \brief Load markers defined in cartesian coordinates. */
  *Marker_Disp_Dir,         /*!< \brief Load markers defined in cartesian coordinates. */
  *Marker_Load_Sine,					/*!< \brief Sine-wave loaded markers defined in cartesian coordinates. */
  *Marker_FlowLoad,					/*!< \brief Flow Load markers. */
  *Marker_Neumann,					/*!< \brief Neumann flow markers. */
  *Marker_Internal,					/*!< \brief Neumann flow markers. */
  *Marker_All_TagBound;				/*!< \brief Global index for markers using grid information. */
  su2double *Dirichlet_Value;    /*!< \brief Specified Dirichlet value at the boundaries. */
  su2double *Exhaust_Temperature_Target;    /*!< \brief Specified total temperatures for nacelle boundaries. */
  su2double *Exhaust_Pressure_Target;    /*!< \brief Specified total pressures for nacelle boundaries. */
  su2double *Inlet_Ttotal;    /*!< \brief Specified total temperatures for inlet boundaries. */
  su2double *Riemann_Var1, *Riemann_Var2;    /*!< \brief Specified values for Riemann boundary. */
  su2double **Riemann_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for Riemann boundaries. */
  su2double *Giles_Var1, *Giles_Var2, *RelaxFactorAverage, *RelaxFactorFourier;    /*!< \brief Specified values for Giles BC. */
  su2double **Giles_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for Giles BC. */
  su2double *Inlet_Ptotal;    /*!< \brief Specified total pressures for inlet boundaries. */
  su2double **Inlet_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
  su2double *Inlet_Temperature;    /*!< \brief Specified temperatures for a supersonic inlet boundaries. */
  su2double *Inlet_Pressure;    /*!< \brief Specified static pressures for supersonic inlet boundaries. */
  su2double **Inlet_Velocity;  /*!< \brief Specified flow velocity vectors for supersonic inlet boundaries. */
  su2double *EngineInflow_Target;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_Mach;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_Pressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_ReverseMassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_TotalPressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_Temperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_TotalTemperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_RamDrag;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_Force;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Inflow_Power;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_Pressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_Temperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_TotalPressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_TotalTemperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_GrossThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_Force;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Exhaust_Power;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Engine_Power;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Engine_Mach;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Engine_Force;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Engine_NetThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Engine_GrossThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Engine_Area;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Outlet_Pressure;    /*!< \brief Specified back pressures (static) for outlet boundaries. */
  su2double *Isothermal_Temperature; /*!< \brief Specified isothermal wall temperatures (static). */
  su2double *Heat_Flux;  /*!< \brief Specified wall heat fluxes. */
  su2double *Displ_Value;    /*!< \brief Specified displacement for displacement boundaries. */
  su2double *Load_Value;    /*!< \brief Specified force for load boundaries. */
  su2double *Damper_Constant;    /*!< \brief Specified constant for damper boundaries. */
  su2double *Load_Dir_Value;    /*!< \brief Specified force for load boundaries defined in cartesian coordinates. */
  su2double *Load_Dir_Multiplier;    /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double *Disp_Dir_Value;    /*!< \brief Specified force for load boundaries defined in cartesian coordinates. */
   su2double *Disp_Dir_Multiplier;    /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double **Load_Dir;  /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
  su2double **Disp_Dir;  /*!< \brief Specified structural displacement direction (unit vector). */
  su2double *Load_Sine_Amplitude;    /*!< \brief Specified amplitude for a sine-wave load. */
  su2double *Load_Sine_Frequency;    /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double **Load_Sine_Dir;  /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
  su2double *FlowLoad_Value;    /*!< \brief Specified force for flow load boundaries. */
  su2double *ActDiskInlet_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_Temperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_TotalTemperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_Pressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_TotalPressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_RamDrag;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_Force;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskInlet_Power;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_Temperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_TotalTemperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_Pressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_TotalPressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_GrossThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_Force;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDiskOutlet_Power;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double **ActDisk_PressJump, **ActDisk_TempJump,  **ActDisk_Omega;
  su2double *ActDisk_DeltaPress;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_DeltaTemp;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_TotalPressRatio;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_TotalTempRatio;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_StaticPressRatio;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_StaticTempRatio;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_Power;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_Mach;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_Force;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Outlet_MassFlow;    /*!< \brief Mass flow for outlet boundaries. */
  su2double *Outlet_Density;    /*!< \brief Avg. density for outlet boundaries. */
  su2double *Outlet_Area;    /*!< \brief Area for outlet boundaries. */
  su2double *Surface_MassFlow;    /*!< \brief Massflow at the boundaries. */
  su2double *Surface_Mach;    /*!< \brief Mach number at the boundaries. */
  su2double *Surface_Temperature;    /*!< \brief Temperature at the boundaries. */
  su2double *Surface_Pressure;    /*!< \brief Pressure at the boundaries. */
  su2double *Surface_Density;    /*!< \brief Density at the boundaries. */
  su2double *Surface_Enthalpy;    /*!< \brief Enthalpy at the boundaries. */
  su2double *Surface_NormalVelocity;    /*!< \brief Normal velocity at the boundaries. */
  su2double *Surface_Uniformity;  /*!< \brief Integral measure of the streamwise uniformity (absolute) at the boundaries (non-dim). */
  su2double *Surface_SecondaryStrength;     /*!< \brief Integral measure of the strength of secondary flows (absolute) at the boundaries (non-dim). */
  su2double *Surface_SecondOverUniform;   /*!< \brief Integral measure of the strength of secondary flows (relative to streamwise) at the boundaries (non-dim). */
  su2double *Surface_MomentumDistortion;    /*!< \brief Integral measure of the streamwise uniformity (relative to plug flow) at the boundaries (non-dim). */
  su2double *Surface_TotalTemperature;   /*!< \brief Total temperature at the boundaries. */
  su2double *Surface_TotalPressure;    /*!< \brief Total pressure at the boundaries. */
  su2double *Surface_PressureDrop;    /*!< \brief Pressure drop between boundaries. */
  su2double *Surface_DC60;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Surface_IDC;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Surface_IDC_Mach;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *Surface_IDR;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_NetThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_BCThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_BCThrust_Old;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_GrossThrust;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_Area;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double *ActDisk_ReverseMassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  su2double **Periodic_RotCenter;  /*!< \brief Rotational center for each periodic boundary. */
  su2double **Periodic_RotAngles;      /*!< \brief Rotation angles for each periodic boundary. */
  su2double **Periodic_Translation;      /*!< \brief Translation vector for each periodic boundary. */
  unsigned short nPeriodic_Index;     /*!< \brief Number of SEND_RECEIVE periodic transformations. */
  su2double **Periodic_Center;         /*!< \brief Rotational center for each SEND_RECEIVE boundary. */
  su2double **Periodic_Rotation;      /*!< \brief Rotation angles for each SEND_RECEIVE boundary. */
  su2double **Periodic_Translate;      /*!< \brief Translation vector for each SEND_RECEIVE boundary. */
  string *Marker_CfgFile_TagBound;			/*!< \brief Global index for markers using config file. */
  unsigned short *Marker_All_KindBC,			/*!< \brief Global index for boundaries using grid information. */
  *Marker_CfgFile_KindBC;		/*!< \brief Global index for boundaries using config file. */
  short *Marker_All_SendRecv;		/*!< \brief Information about if the boundary is sended (+), received (-). */
  short *Marker_All_PerBound;	/*!< \brief Global index for periodic bc using the grid information. */
  unsigned long nExtIter;			/*!< \brief Number of external iterations. */
  unsigned long ExtIter;			/*!< \brief Current external iteration number. */
  unsigned long ExtIter_OffSet;			/*!< \brief External iteration number offset. */
  unsigned long IntIter;			/*!< \brief Current internal iteration number. */
  unsigned long OuterIter;			/*!< \brief Current Outer Iteration for multizone problems. */
  unsigned long Unst_nIntIter;			/*!< \brief Number of internal iterations (Dual time Method). */
  unsigned long Dyn_nIntIter;			/*!< \brief Number of internal iterations (Newton-Raphson Method for nonlinear structural analysis). */
  long Unst_RestartIter;			/*!< \brief Iteration number to restart an unsteady simulation (Dual time Method). */
  long Unst_AdjointIter;			/*!< \brief Iteration number to begin the reverse time integration in the direct solver for the unsteady adjoint. */
  long Iter_Avg_Objective;			/*!< \brief Iteration the number of time steps to be averaged, counting from the back */
  long Dyn_RestartIter;                         /*!< \brief Iteration number to restart a dynamic structural analysis. */
  su2double PhysicalTime;                       /*!< \brief Physical time at the current iteration in the solver for unsteady problems. */
  unsigned short nLevels_TimeAccurateLTS;       /*!< \brief Number of time levels for time accurate local time stepping. */
  unsigned short nTimeDOFsADER_DG;              /*!< \brief Number of time DOFs used in the predictor step of ADER-DG. */
  su2double *TimeDOFsADER_DG;                   /*!< \brief The location of the ADER-DG time DOFs on the interval [-1,1]. */
  unsigned short nTimeIntegrationADER_DG;       /*!< \brief Number of time integration points ADER-DG. */
  su2double *TimeIntegrationADER_DG;            /*!< \brief The location of the ADER-DG time integration points on the interval [-1,1]. */
  su2double *WeightsIntegrationADER_DG;         /*!< \brief The weights of the ADER-DG time integration points on the interval [-1,1]. */
  unsigned short nRKStep;			/*!< \brief Number of steps of the explicit Runge-Kutta method. */
  su2double *RK_Alpha_Step;			/*!< \brief Runge-Kutta beta coefficients. */
  unsigned short nMGLevels;		/*!< \brief Number of multigrid levels (coarse levels). */
  unsigned short nCFL;			/*!< \brief Number of CFL, one for each multigrid level. */
  su2double
  CFLRedCoeff_Turb,		/*!< \brief CFL reduction coefficient on the LevelSet problem. */
  CFLRedCoeff_AdjFlow,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
  CFLRedCoeff_AdjTurb,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
  CFLFineGrid,		/*!< \brief CFL of the finest grid. */
  CFLSolid,       /*!< \brief CFL in (heat) solid solvers. */
  Max_DeltaTime,  		/*!< \brief Max delta time. */
  Unst_CFL;		/*!< \brief Unsteady CFL number. */
  bool ReorientElements;		/*!< \brief Flag for enabling element reorientation. */
  bool AddIndNeighbor;			/*!< \brief Include indirect neighbor in the agglomeration process. */
  unsigned short nDV,		/*!< \brief Number of design variables. */
  nObj, nObjW;              /*! \brief Number of objective functions. */
  unsigned short* nDV_Value;		/*!< \brief Number of values for each design variable (might be different than 1 if we allow arbitrary movement). */
  unsigned short nFFDBox;		/*!< \brief Number of ffd boxes. */
  unsigned short nGridMovement;		/*!< \brief Number of grid movement types specified. */
  unsigned short nTurboMachineryKind; 	/*!< \brief Number turbomachinery types specified. */
  unsigned short nParamDV;		/*!< \brief Number of parameters of the design variable. */
  string DV_Filename;      /*!< \brief Filename for providing surface positions from an external parameterization. */
  string DV_Unordered_Sens_Filename;      /*!< \brief Filename of volume sensitivities in an unordered ASCII format. */
  string DV_Sens_Filename;      /*!< \brief Filename of surface sensitivities written to an unordered ASCII format. */
  unsigned short Sensitivity_FileFormat; /*!< \brief Format of the input volume sensitivity files (SU2_DOT). */
  su2double **ParamDV;				/*!< \brief Parameters of the design variable. */
  su2double **CoordFFDBox;				/*!< \brief Coordinates of the FFD boxes. */
  unsigned short **DegreeFFDBox;	/*!< \brief Degree of the FFD boxes. */
  string *FFDTag;				/*!< \brief Parameters of the design variable. */
  string *TagFFDBox;				/*!< \brief Tag of the FFD box. */
  unsigned short GeometryMode;			/*!< \brief Gemoetry mode (analysis or gradient computation). */
  unsigned short MGCycle;			/*!< \brief Kind of multigrid cycle. */
  unsigned short FinestMesh;		/*!< \brief Finest mesh for the full multigrid approach. */
  unsigned short nFFD_Fix_IDir, nFFD_Fix_JDir, nFFD_Fix_KDir;                 /*!< \brief Number of planes fixed in the FFD. */
  unsigned short nMG_PreSmooth,                 /*!< \brief Number of MG pre-smooth parameters found in config file. */
  nMG_PostSmooth,                             /*!< \brief Number of MG post-smooth parameters found in config file. */
  nMG_CorrecSmooth;                           /*!< \brief Number of MG correct-smooth parameters found in config file. */
  short *FFD_Fix_IDir, *FFD_Fix_JDir, *FFD_Fix_KDir;	/*!< \brief Exact sections. */
  unsigned short *MG_PreSmooth,	/*!< \brief Multigrid Pre smoothing. */
  *MG_PostSmooth,					/*!< \brief Multigrid Post smoothing. */
  *MG_CorrecSmooth;					/*!< \brief Multigrid Jacobi implicit smoothing of the correction. */
  su2double *LocationStations;   /*!< \brief Airfoil sections in wing slicing subroutine. */
  su2double *NacelleLocation;   /*!< \brief Definition of the nacelle location. */
  unsigned short Kind_Solver,	/*!< \brief Kind of solver Euler, NS, Continuous adjoint, etc.  */
  *Kind_Solver_PerZone,  /*!< \brief Kind of solvers for each zone Euler, NS, Continuous adjoint, etc.  */
  Kind_MZSolver,         /*!< \brief Kind of multizone solver.  */
  Kind_FluidModel,			/*!< \brief Kind of the Fluid Model: Ideal or Van der Walls, ... . */
  Kind_ViscosityModel,			/*!< \brief Kind of the Viscosity Model*/
  Kind_ConductivityModel,			/*!< \brief Kind of the Thermal Conductivity Model*/
  Kind_ConductivityModel_Turb,      /*!< \brief Kind of the Turbulent Thermal Conductivity Model*/
  Kind_FreeStreamOption,			/*!< \brief Kind of free stream option to choose if initializing with density or temperature  */
  Kind_InitOption,			/*!< \brief Kind of Init option to choose if initializing with Reynolds number or with thermodynamic conditions   */
  Kind_GasModel,				/*!< \brief Kind of the Gas Model. */
  Kind_DensityModel,				/*!< \brief Kind of the density model for incompressible flows. */
  *Kind_GridMovement,    /*!< \brief Kind of the unsteady mesh movement. */
  Kind_Gradient_Method,		/*!< \brief Numerical method for computation of spatial gradients. */
  Kind_Deform_Linear_Solver, /*!< Numerical method to deform the grid */
  Kind_Deform_Linear_Solver_Prec,		/*!< \brief Preconditioner of the linear solver. */
  Kind_Linear_Solver,		/*!< \brief Numerical solver for the implicit scheme. */
  Kind_Linear_Solver_FSI_Struc,	 /*!< \brief Numerical solver for the structural part in FSI problems. */
  Kind_Linear_Solver_Prec,		/*!< \brief Preconditioner of the linear solver. */
  Kind_Linear_Solver_Prec_FSI_Struc,		/*!< \brief Preconditioner of the linear solver for the structural part in FSI problems. */
  Kind_AdjTurb_Linear_Solver,		/*!< \brief Numerical solver for the turbulent adjoint implicit scheme. */
  Kind_AdjTurb_Linear_Prec,		/*!< \brief Preconditioner of the turbulent adjoint linear solver. */
  Kind_DiscAdj_Linear_Solver, /*!< \brief Linear solver for the discrete adjoint system. */
  Kind_DiscAdj_Linear_Prec,  /*!< \brief Preconditioner of the discrete adjoint linear solver. */
  Kind_DiscAdj_Linear_Solver_FSI_Struc, /*!< \brief Linear solver for the discrete adjoint system in the structural side of FSI problems. */
  Kind_DiscAdj_Linear_Prec_FSI_Struc,   /*!< \brief Preconditioner of the discrete adjoint linear solver in the structural side of FSI problems. */
  Kind_SlopeLimit,				/*!< \brief Global slope limiter. */
  Kind_SlopeLimit_Flow,		/*!< \brief Slope limiter for flow equations.*/
  Kind_SlopeLimit_Turb,		/*!< \brief Slope limiter for the turbulence equation.*/
  Kind_SlopeLimit_AdjTurb,	/*!< \brief Slope limiter for the adjoint turbulent equation.*/
  Kind_SlopeLimit_AdjFlow,	/*!< \brief Slope limiter for the adjoint equation.*/
  Kind_TimeNumScheme,			/*!< \brief Global explicit or implicit time integration. */
  Kind_TimeIntScheme_Flow,	/*!< \brief Time integration for the flow equations. */
  Kind_TimeIntScheme_FEM_Flow,  /*!< \brief Time integration for the flow equations. */
  Kind_ADER_Predictor,          /*!< \brief Predictor step of the ADER-DG time integration scheme. */
  Kind_TimeIntScheme_AdjFlow,		/*!< \brief Time integration for the adjoint flow equations. */
  Kind_TimeIntScheme_Turb,	/*!< \brief Time integration for the turbulence model. */
  Kind_TimeIntScheme_AdjTurb,	/*!< \brief Time integration for the adjoint turbulence model. */
  Kind_TimeIntScheme_Heat,	/*!< \brief Time integration for the wave equations. */
  Kind_TimeStep_Heat, /*!< \brief Time stepping method for the (fvm) heat equation. */
  Kind_TimeIntScheme_FEA,	/*!< \brief Time integration for the FEA equations. */
  Kind_SpaceIteScheme_FEA,	/*!< \brief Iterative scheme for nonlinear structural analysis. */
  Kind_ConvNumScheme,			/*!< \brief Global definition of the convective term. */
  Kind_ConvNumScheme_Flow,	/*!< \brief Centered or upwind scheme for the flow equations. */
  Kind_ConvNumScheme_FEM_Flow,  /*!< \brief Finite element scheme for the flow equations. */
  Kind_ConvNumScheme_Heat,	/*!< \brief Centered or upwind scheme for the flow equations. */
  Kind_ConvNumScheme_AdjFlow,		/*!< \brief Centered or upwind scheme for the adjoint flow equations. */
  Kind_ConvNumScheme_Turb,	/*!< \brief Centered or upwind scheme for the turbulence model. */
  Kind_ConvNumScheme_AdjTurb,	/*!< \brief Centered or upwind scheme for the adjoint turbulence model. */
  Kind_ConvNumScheme_Template,	/*!< \brief Centered or upwind scheme for the level set equation. */
  Kind_Centered,				/*!< \brief Centered scheme. */
  Kind_Centered_Flow,			/*!< \brief Centered scheme for the flow equations. */
  Kind_Centered_AdjFlow,			/*!< \brief Centered scheme for the adjoint flow equations. */
  Kind_Centered_Turb,			/*!< \brief Centered scheme for the turbulence model. */
  Kind_Centered_AdjTurb,		/*!< \brief Centered scheme for the adjoint turbulence model. */
  Kind_Centered_Template,		/*!< \brief Centered scheme for the template model. */
  Kind_Upwind,				/*!< \brief Upwind scheme. */
  Kind_Upwind_Flow,			/*!< \brief Upwind scheme for the flow equations. */
  Kind_Upwind_AdjFlow,			/*!< \brief Upwind scheme for the adjoint flow equations. */
  Kind_Upwind_Turb,			/*!< \brief Upwind scheme for the turbulence model. */
  Kind_Upwind_AdjTurb,		/*!< \brief Upwind scheme for the adjoint turbulence model. */
  Kind_Upwind_Template,			/*!< \brief Upwind scheme for the template model. */
  Kind_FEM,                     /*!< \brief Finite element scheme for the flow equations. */
  Kind_FEM_Flow,                        /*!< \brief Finite element scheme for the flow equations. */
  Kind_FEM_DG_Shock,      /*!< \brief Shock capturing method for the FEM DG solver. */
  Kind_Matrix_Coloring,   /*!< \brief Type of matrix coloring for sparse Jacobian computation. */
  Kind_Solver_Fluid_FSI,		/*!< \brief Kind of solver for the fluid in FSI applications. */
  Kind_Solver_Struc_FSI,		/*!< \brief Kind of solver for the structure in FSI applications. */
  Kind_BGS_RelaxMethod;				/*!< \brief Kind of relaxation method for Block Gauss Seidel method in FSI problems. */
  bool Energy_Equation;         /*!< \brief Solve the energy equation for incompressible flows. */
  bool MUSCL,		/*!< \brief MUSCL scheme .*/
  MUSCL_Flow,		/*!< \brief MUSCL scheme for the flow equations.*/
  MUSCL_Turb,	 /*!< \brief MUSCL scheme for the turbulence equations.*/
  MUSCL_Heat,	 /*!< \brief MUSCL scheme for the (fvm) heat equation.*/
  MUSCL_AdjFlow,		/*!< \brief MUSCL scheme for the adj flow equations.*/
  MUSCL_AdjTurb, 	/*!< \brief MUSCL scheme for the adj turbulence equations.*/
  Use_Accurate_Jacobians;   /*!< \brief Use numerically computed Jacobians for AUSM+up(2) and SLAU(2). */
  bool EulerPersson;        /*!< \brief Boolean to determine whether this is an Euler simulation with Persson shock capturing. */
  bool FSI_Problem,			/*!< \brief Boolean to determine whether the simulation is FSI or not. */
  ZoneSpecific_Problem,   /*!< \brief Boolean to determine whether we wish to use zone-specific solvers. */
  Multizone_Problem;      /*!< \brief Boolean to determine whether we are solving a multizone problem. */
  unsigned short nID_DV;  /*!< \brief ID for the region of FEM when computed using direct differentiation. */
  bool AD_Mode;         /*!< \brief Algorithmic Differentiation support. */
  bool AD_Preaccumulation;   /*!< \brief Enable or disable preaccumulation in the AD mode. */
  unsigned short Kind_Material_Compress,	/*!< \brief Determines if the material is compressible or incompressible (structural analysis). */
  Kind_Material,			/*!< \brief Determines the material model to be used (structural analysis). */
  Kind_Struct_Solver,		/*!< \brief Determines the geometric condition (small or large deformations) for structural analysis. */
  Kind_DV_FEA;				/*!< \brief Kind of Design Variable for FEA problems.*/
  unsigned short Kind_Turb_Model;			/*!< \brief Turbulent model definition. */
  unsigned short Kind_SGS_Model;                        /*!< \brief LES SGS model definition. */
  unsigned short Kind_Trans_Model,			/*!< \brief Transition model definition. */
  Kind_ActDisk, Kind_Engine_Inflow, Kind_Inlet, *Kind_Inc_Inlet, *Kind_Inc_Outlet, *Kind_Data_Riemann, *Kind_Data_Giles;           /*!< \brief Kind of inlet boundary treatment. */
  unsigned short nInc_Inlet;  /*!< \brief Number of inlet boundary treatment types listed. */
  unsigned short nInc_Outlet;  /*!< \brief Number of inlet boundary treatment types listed. */
  su2double Inc_Inlet_Damping;  /*!< \brief Damping factor applied to the iterative updates to the velocity at a pressure inlet in incompressible flow. */
  su2double Inc_Outlet_Damping; /*!< \brief Damping factor applied to the iterative updates to the pressure at a mass flow outlet in incompressible flow. */
  bool Inc_Inlet_UseNormal;    /*!< \brief Flag for whether to use the local normal as the flow direction for an incompressible pressure inlet. */
  su2double Linear_Solver_Error;		/*!< \brief Min error of the linear solver for the implicit formulation. */
  su2double Deform_Linear_Solver_Error;    /*!< \brief Min error of the linear solver for the implicit formulation. */
  su2double Linear_Solver_Error_FSI_Struc;		/*!< \brief Min error of the linear solver for the implicit formulation in the structural side for FSI problems . */
  su2double Linear_Solver_Error_Heat;        /*!< \brief Min error of the linear solver for the implicit formulation in the fvm heat solver . */
  su2double Linear_Solver_Smoother_Relaxation;  /*!< \brief Relaxation factor for iterative linear smoothers. */
  unsigned long Linear_Solver_Iter;		/*!< \brief Max iterations of the linear solver for the implicit formulation. */
  unsigned long Deform_Linear_Solver_Iter;   /*!< \brief Max iterations of the linear solver for the implicit formulation. */
  unsigned long Linear_Solver_Iter_FSI_Struc;		/*!< \brief Max iterations of the linear solver for FSI applications and structural solver. */
  unsigned long Linear_Solver_Iter_Heat;       /*!< \brief Max iterations of the linear solver for the implicit formulation in the fvm heat solver. */
  unsigned long Linear_Solver_Restart_Frequency;   /*!< \brief Restart frequency of the linear solver for the implicit formulation. */
  unsigned short Linear_Solver_ILU_n;		/*!< \brief ILU fill=in level. */
  su2double SemiSpan;		/*!< \brief Wing Semi span. */
  su2double Roe_Kappa;		/*!< \brief Relaxation of the Roe scheme. */
  su2double Relaxation_Factor_Flow;		/*!< \brief Relaxation coefficient of the linear solver mean flow. */
  su2double Relaxation_Factor_Turb;		/*!< \brief Relaxation coefficient of the linear solver turbulence. */
  su2double Relaxation_Factor_AdjFlow;		/*!< \brief Relaxation coefficient of the linear solver adjoint mean flow. */
  su2double Relaxation_Factor_CHT;  /*!< \brief Relaxation coefficient for the update of conjugate heat variables. */
  su2double AdjTurb_Linear_Error;		/*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  su2double EntropyFix_Coeff;              /*!< \brief Entropy fix coefficient. */
  unsigned short AdjTurb_Linear_Iter;		/*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  su2double *Stations_Bounds;                  /*!< \brief Airfoil section limit. */
  unsigned short nLocationStations,      /*!< \brief Number of section cuts to make when outputting mesh and cp . */
  nWingStations;               /*!< \brief Number of section cuts to make when calculating internal volume. */
  su2double* Kappa_Flow,           /*!< \brief Numerical dissipation coefficients for the flow equations. */
  *Kappa_AdjFlow,                  /*!< \brief Numerical dissipation coefficients for the adjoint flow equations. */
  *Kappa_Heat;                    /*!< \brief Numerical dissipation coefficients for the (fvm) heat equation. */  
  su2double* FFD_Axis;       /*!< \brief Numerical dissipation coefficients for the adjoint equations. */
  su2double Kappa_1st_AdjFlow,	/*!< \brief JST 1st order dissipation coefficient for adjoint flow equations (coarse multigrid levels). */
  Kappa_2nd_AdjFlow,			/*!< \brief JST 2nd order dissipation coefficient for adjoint flow equations. */
  Kappa_4th_AdjFlow,			/*!< \brief JST 4th order dissipation coefficient for adjoint flow equations. */
  Kappa_1st_Flow,			/*!< \brief JST 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
  Kappa_2nd_Flow,			/*!< \brief JST 2nd order dissipation coefficient for flow equations. */
  Kappa_4th_Flow,			/*!< \brief JST 4th order dissipation coefficient for flow equations. */
  Kappa_2nd_Heat,     /*!< \brief 2nd order dissipation coefficient for heat equation. */
  Kappa_4th_Heat,     /*!< \brief 4th order dissipation coefficient for heat equation. */
  Cent_Jac_Fix_Factor;/*!< \brief Multiply the dissipation contribution to the Jacobian of central schemes by this factor to make the global matrix more diagonal dominant. */
  su2double Geo_Waterline_Location; /*!< \brief Location of the waterline. */
  
  su2double Min_Beta_RoeTurkel,		/*!< \brief Minimum value of Beta for the Roe-Turkel low Mach preconditioner. */
  Max_Beta_RoeTurkel;		/*!< \brief Maximum value of Beta for the Roe-Turkel low Mach preconditioner. */
  unsigned long GridDef_Nonlinear_Iter, /*!< \brief Number of nonlinear increments for grid deformation. */
  GridDef_Linear_Iter; /*!< \brief Number of linear smoothing iterations for grid deformation. */
  unsigned short Deform_Stiffness_Type; /*!< \brief Type of element stiffness imposed for FEA mesh deformation. */
  bool Deform_Output;  /*!< \brief Print the residuals during mesh deformation to the console. */
  su2double Deform_Tol_Factor; /*!< Factor to multiply smallest volume for deform tolerance (0.001 default) */
  su2double Deform_Coeff; /*!< Deform coeffienct */
  su2double Deform_Limit; /*!< Deform limit */
  unsigned short FFD_Continuity; /*!< Surface continuity at the intersection with the FFD */
  unsigned short FFD_CoordSystem; /*!< Define the coordinates system */
  su2double Deform_ElasticityMod, Deform_PoissonRatio; /*!< young's modulus and poisson ratio for volume deformation stiffness model */
  bool Visualize_Surface_Def;  /*!< \brief Flag to visualize the surface deformacion in SU2_DEF. */
  bool Visualize_Volume_Def; /*!< \brief Flag to visualize the volume deformation in SU2_DEF. */
  bool FFD_Symmetry_Plane;	/*!< \brief FFD symmetry plane. */
  su2double Mach;		/*!< \brief Mach number. */
  su2double Reynolds;	/*!< \brief Reynolds number. */
  su2double Froude;	/*!< \brief Froude number. */
  su2double Length_Reynolds;	/*!< \brief Reynolds length (dimensional). */
  su2double AoA,			/*!< \brief Angle of attack (just external flow). */
  iH, AoS, AoA_Offset, AoS_Offset, AoA_Sens;		/*!< \brief Angle of sideSlip (just external flow). */
  bool Fixed_CL_Mode;			/*!< \brief Activate fixed CL mode (external flow only). */
  bool Fixed_CM_Mode;			/*!< \brief Activate fixed CL mode (external flow only). */
  bool Eval_dOF_dCX;			/*!< \brief Activate fixed CL mode (external flow only). */
  bool Discard_InFiles; /*!< \brief Discard angle of attack in solution and geometry files. */
  su2double Target_CL;			/*!< \brief Specify a target CL instead of AoA (external flow only). */
  su2double Target_CM;			/*!< \brief Specify a target CL instead of AoA (external flow only). */
  su2double Total_CM;			/*!< \brief Specify a target CL instead of AoA (external flow only). */
  su2double Total_CD;			/*!< \brief Specify a target CL instead of AoA (external flow only). */
  su2double dCL_dAlpha;        /*!< \brief value of dCl/dAlpha. */
  su2double dCM_diH;        /*!< \brief value of dCM/dHi. */
  unsigned long Iter_Fixed_CL;			/*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Iter_Fixed_CM;			/*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Iter_Fixed_NetThrust;			/*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Iter_dCL_dAlpha;   /*!< \brief Number of iterations to evaluate dCL_dAlpha. */
  unsigned long Update_Alpha;			/*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Update_iH;			/*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  unsigned long Update_BCThrust;			/*!< \brief Iterations to re-evaluate the angle of attack (external flow only). */
  su2double dNetThrust_dBCThrust;        /*!< \brief value of dCl/dAlpha. */
  bool Update_BCThrust_Bool;			/*!< \brief Boolean flag for whether to update the AoA for fixed lift mode on a given iteration. */
  bool Update_AoA;			/*!< \brief Boolean flag for whether to update the AoA for fixed lift mode on a given iteration. */
  bool Update_HTPIncidence;			/*!< \brief Boolean flag for whether to update the AoA for fixed lift mode on a given iteration. */
  su2double ChargeCoeff;		/*!< \brief Charge coefficient (just for poisson problems). */
  unsigned short Cauchy_Func_Flow,	/*!< \brief Function where to apply the convergence criteria in the flow problem. */
  Cauchy_Func_AdjFlow,				/*!< \brief Function where to apply the convergence criteria in the adjoint problem. */
  Cauchy_Elems;						/*!< \brief Number of elements to evaluate. */
  unsigned short Residual_Func_Flow;	/*!< \brief Equation to apply residual convergence to. */
  unsigned short Res_FEM_CRIT;  /*!< \brief Criteria to apply to the FEM convergence (absolute/relative). */
  unsigned long StartConv_Iter;	/*!< \brief Start convergence criteria at iteration. */
  su2double Cauchy_Eps;	/*!< \brief Epsilon used for the convergence. */
  unsigned long Wrt_Sol_Freq,	/*!< \brief Writing solution frequency. */
  Wrt_Sol_Freq_DualTime,	/*!< \brief Writing solution frequency for Dual Time. */
  Wrt_Con_Freq,				/*!< \brief Writing convergence history frequency. */
  Wrt_Con_Freq_DualTime;				/*!< \brief Writing convergence history frequency. */
  bool Wrt_Unsteady;  /*!< \brief Write unsteady data adding header and prefix. */
  bool Wrt_Dynamic;  		/*!< \brief Write dynamic data adding header and prefix. */
  bool Restart,	/*!< \brief Restart solution (for direct, adjoint, and linearized problems).*/
  Wrt_Binary_Restart,	/*!< \brief Write binary SU2 native restart files.*/
  Read_Binary_Restart,	/*!< \brief Read binary SU2 native restart files.*/
  Restart_Flow;	/*!< \brief Restart flow solution for adjoint and linearized problems. */
  unsigned short nMarker_Monitoring,	/*!< \brief Number of markers to monitor. */
  nMarker_Designing,					/*!< \brief Number of markers for the objective function. */
  nMarker_GeoEval,					/*!< \brief Number of markers for the objective function. */
  nMarker_ZoneInterface, /*!< \brief Number of markers in the zone interface. */
  nMarker_Plotting,					/*!< \brief Number of markers to plot. */
  nMarker_Analyze,					/*!< \brief Number of markers to plot. */
  nMarker_Moving,               /*!< \brief Number of markers in motion (DEFORMING, MOVING_WALL, or FLUID_STRUCTURE). */
  nMarker_PyCustom,               /*!< \brief Number of markers that are customizable in Python. */
  nMarker_DV,               /*!< \brief Number of markers affected by the design variables. */
  nMarker_WallFunctions;    /*!< \brief Number of markers for which wall functions must be applied. */
  string *Marker_Monitoring,     /*!< \brief Markers to monitor. */
  *Marker_Designing,         /*!< \brief Markers to plot. */
  *Marker_GeoEval,         /*!< \brief Markers to plot. */
  *Marker_Plotting,          /*!< \brief Markers to plot. */
  *Marker_Analyze,          /*!< \brief Markers to plot. */
  *Marker_ZoneInterface,          /*!< \brief Markers in the FSI interface. */
  *Marker_Moving,            /*!< \brief Markers in motion (DEFORMING, MOVING_WALL, or FLUID_STRUCTURE). */
  *Marker_PyCustom,            /*!< \brief Markers that are customizable in Python. */
  *Marker_DV,            /*!< \brief Markers affected by the design variables. */
  *Marker_WallFunctions; /*!< \brief Markers for which wall functions must be applied. */
  unsigned short  nConfig_Files;          /*!< \brief Number of config files for multiphysics problems. */
  string *Config_Filenames;               /*!< \brief List of names for configuration files. */
  unsigned short  *Kind_WallFunctions;        /*!< \brief The kind of wall function to use for the corresponding markers. */
  unsigned short  **IntInfo_WallFunctions;    /*!< \brief Additional integer information for the wall function markers. */
  su2double       **DoubleInfo_WallFunctions; /*!< \brief Additional double information for the wall function markers. */
  unsigned short  *Marker_All_Monitoring,        /*!< \brief Global index for monitoring using the grid information. */
  *Marker_All_GeoEval,       /*!< \brief Global index for geometrical evaluation. */
  *Marker_All_Plotting,        /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_Analyze,        /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_ZoneInterface,        /*!< \brief Global index for FSI interface markers using the grid information. */
  *Marker_All_Turbomachinery,        /*!< \brief Global index for Turbomachinery markers using the grid information. */
  *Marker_All_TurbomachineryFlag,        /*!< \brief Global index for Turbomachinery markers flag using the grid information. */
  *Marker_All_MixingPlaneInterface,        /*!< \brief Global index for MixingPlane interface markers using the grid information. */    
  *Marker_All_DV,          /*!< \brief Global index for design variable markers using the grid information. */
  *Marker_All_Moving,          /*!< \brief Global index for moving surfaces using the grid information. */
  *Marker_All_PyCustom,                 /*!< \brief Global index for Python customizable surfaces using the grid information. */
  *Marker_All_Designing,         /*!< \brief Global index for moving using the grid information. */
  *Marker_CfgFile_Monitoring,     /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Designing,      /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_GeoEval,      /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Plotting,     /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_Analyze,     /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_ZoneInterface,     /*!< \brief Global index for FSI interface using the config information. */
  *Marker_CfgFile_Turbomachinery,     /*!< \brief Global index for Turbomachinery  using the config information. */
  *Marker_CfgFile_TurbomachineryFlag,     /*!< \brief Global index for Turbomachinery flag using the config information. */
  *Marker_CfgFile_MixingPlaneInterface,     /*!< \brief Global index for MixingPlane interface using the config information. */
  *Marker_CfgFile_Moving,       /*!< \brief Global index for moving surfaces using the config information. */
  *Marker_CfgFile_PyCustom,        /*!< \brief Global index for Python customizable surfaces using the config information. */
  *Marker_CfgFile_DV,       /*!< \brief Global index for design variable markers using the config information. */
  *Marker_CfgFile_PerBound;     /*!< \brief Global index for periodic boundaries using the config information. */
  string *PlaneTag;      /*!< \brief Global index for the plane adaptation (upper, lower). */
  su2double DualVol_Power;			/*!< \brief Power for the dual volume in the grid adaptation sensor. */
  su2double *nBlades;						/*!< \brief number of blades for turbomachinery computation. */
  unsigned short Analytical_Surface;	/*!< \brief Information about the analytical definition of the surface for grid adaptation. */
  unsigned short Geo_Description;	/*!< \brief Description of the geometry. */
  unsigned short Mesh_FileFormat;	/*!< \brief Mesh input format. */
  unsigned short Output_FileFormat;	/*!< \brief Format of the output files. */
  unsigned short ActDisk_Jump;	/*!< \brief Format of the output files. */
  bool CFL_Adapt;      /*!< \brief Adaptive CFL number. */
  bool HB_Precondition;    /*< \brief Flag to turn on harmonic balance source term preconditioning */
  su2double RefArea,		/*!< \brief Reference area for coefficient computation. */
  RefElemLength,				/*!< \brief Reference element length for computing the slope limiting epsilon. */
  RefSharpEdges,				/*!< \brief Reference coefficient for detecting sharp edges. */
  RefLength,			/*!< \brief Reference length for moment computation. */
  *RefOriginMoment,           /*!< \brief Origin for moment computation. */
  *RefOriginMoment_X,      /*!< \brief X Origin for moment computation. */
  *RefOriginMoment_Y,      /*!< \brief Y Origin for moment computation. */
  *RefOriginMoment_Z,      /*!< \brief Z Origin for moment computation. */
  *CFL_AdaptParam,      /*!< \brief Information about the CFL ramp. */
  *RelaxFactor_Giles,      /*!< \brief Information about the under relaxation factor for Giles BC. */
  *CFL,
  *HTP_Axis,      /*!< \brief Location of the HTP axis. */
  DomainVolume;		/*!< \brief Volume of the computational grid. */
  unsigned short nRefOriginMoment_X,    /*!< \brief Number of X-coordinate moment computation origins. */
  nRefOriginMoment_Y,           /*!< \brief Number of Y-coordinate moment computation origins. */
  nRefOriginMoment_Z;           /*!< \brief Number of Z-coordinate moment computation origins. */
  string Mesh_FileName,			/*!< \brief Mesh input file. */
  Mesh_Out_FileName,				/*!< \brief Mesh output file. */
  Solution_FlowFileName,			/*!< \brief Flow solution input file. */
  Solution_LinFileName,			/*!< \brief Linearized flow solution input file. */
  Solution_AdjFileName,			/*!< \brief Adjoint solution input file for drag functional. */
  Solution_FEMFileName,			/*!< \brief Solution input file for structural problem. */
  Solution_AdjFEMFileName,     /*!< \brief Adjoint solution input file for structural problem. */
  Flow_FileName,					/*!< \brief Flow variables output file. */
  Structure_FileName,					/*!< \brief Structure variables output file. */
  SurfStructure_FileName,					/*!< \brief Surface structure variables output file. */
  AdjStructure_FileName,         /*!< \brief Structure variables output file. */
  AdjSurfStructure_FileName,         /*!< \brief Surface structure variables output file. */
  SurfHeat_FileName,					/*!< \brief Surface structure variables output file. */
  Heat_FileName,					/*!< \brief Heat variables output file. */
  Residual_FileName,				/*!< \brief Residual variables output file. */
  Conv_FileName,					/*!< \brief Convergence history output file. */
  Breakdown_FileName,			    /*!< \brief Breakdown output file. */
  Conv_FileName_FSI,					/*!< \brief Convergence history output file. */
  Restart_FlowFileName,			/*!< \brief Restart file for flow variables. */
  Restart_HeatFileName,			/*!< \brief Restart file for heat variables. */
  Restart_AdjFileName,			/*!< \brief Restart file for adjoint variables, drag functional. */
  Restart_FEMFileName,			/*!< \brief Restart file for FEM elasticity. */
  Restart_AdjFEMFileName,      /*!< \brief Restart file for FEM elasticity. */
  Adj_FileName,					/*!< \brief Output file with the adjoint variables. */
  ObjFunc_Grad_FileName,			/*!< \brief Gradient of the objective function. */
  ObjFunc_Value_FileName,			/*!< \brief Objective function. */
  SurfFlowCoeff_FileName,			/*!< \brief Output file with the flow variables on the surface. */
  SurfAdjCoeff_FileName,			/*!< \brief Output file with the adjoint variables on the surface. */
  New_SU2_FileName,       		/*!< \brief Output SU2 mesh file converted from CGNS format. */
  SurfSens_FileName,			/*!< \brief Output file for the sensitivity on the surface (discrete adjoint). */
  VolSens_FileName;			/*!< \brief Output file for the sensitivity in the volume (discrete adjoint). */
  bool Low_MemoryOutput,      /*!< \brief Output less information for lower memory use */
  Wrt_Output,                 /*!< \brief Write any output files */
  Wrt_Vol_Sol,                /*!< \brief Write a volume solution file */
  Wrt_Srf_Sol,                /*!< \brief Write a surface solution file */
  Wrt_Csv_Sol,                /*!< \brief Write a surface comma-separated values solution file */
  Wrt_Crd_Sol,                /*!< \brief Write a binary file with the grid coordinates only. */
  Wrt_Residuals,              /*!< \brief Write residuals to solution file */
  Wrt_Surface,                /*!< \brief Write solution at each surface */
  Wrt_Limiters,              /*!< \brief Write residuals to solution file */
  Wrt_SharpEdges,              /*!< \brief Write residuals to solution file */
  Wrt_Halo,                   /*!< \brief Write rind layers in solution files */
  Wrt_Performance,            /*!< \brief Write the performance summary at the end of a calculation.  */
  Wrt_InletFile,                   /*!< \brief Write a template inlet profile file */
  Wrt_Slice,                   /*!< \brief Write 1D slice of a 2D cartesian solution */
  Wrt_Projected_Sensitivity,   /*!< \brief Write projected sensitivities (dJ/dx) on surfaces to ASCII file. */
  Plot_Section_Forces;       /*!< \brief Write sectional forces for specified markers. */
  unsigned short Console_Output_Verb,  /*!< \brief Level of verbosity for console output */
  Kind_Average;        /*!< \brief Particular average for the marker analyze. */
  unsigned short nPolyCoeffs; /*!< \brief Number of coefficients in temperature polynomial fits for fluid models. */
  su2double Gamma,			/*!< \brief Ratio of specific heats of the gas. */
  Bulk_Modulus,			/*!< \brief Value of the bulk modulus for incompressible flows. */
  Beta_Factor,			/*!< \brief Value of the epsilon^2 multiplier for Beta for the incompressible preconditioner. */
  Gas_Constant,     /*!< \brief Specific gas constant. */
  Gas_ConstantND,     /*!< \brief Non-dimensional specific gas constant. */
  Molecular_Weight,     /*!< \brief Molecular weight of an incompressible ideal gas (g/mol). */
  Specific_Heat_Cp,     /*!< \brief Specific heat at constant pressure. */
  Specific_Heat_CpND,     /*!< \brief Non-dimensional specific heat at constant pressure. */
  Specific_Heat_Cp_Solid, /*!< \brief Specific heat in solids. */
  Specific_Heat_Cv,     /*!< \brief Specific heat at constant volume. */
  Specific_Heat_CvND,     /*!< \brief Non-dimensional specific heat at constant volume. */
  Thermal_Expansion_Coeff,     /*!< \brief Thermal expansion coefficient. */
  Thermal_Expansion_CoeffND,     /*!< \brief Non-dimensional thermal expansion coefficient. */
  Inc_Density_Ref,    /*!< \brief Reference density for custom incompressible non-dim. */
  Inc_Velocity_Ref,    /*!< \brief Reference velocity for custom incompressible non-dim. */
  Inc_Temperature_Ref,    /*!< \brief Reference temperature for custom incompressible non-dim. */
  Inc_Density_Init,    /*!< \brief Initial density for incompressible flows. */
  *Inc_Velocity_Init,    /*!< \brief Initial velocity vector for incompressible flows. */
  Inc_Temperature_Init,    /*!< \brief Initial temperature for incompressible flows w/ heat transfer. */
  Heat_Flux_Ref,  /*!< \brief Reference heat flux for non-dim. */
  Gas_Constant_Ref, /*!< \brief Reference specific gas constant. */
  Temperature_Critical,   /*!< \brief Critical Temperature for real fluid model.  */
  Pressure_Critical,   /*!< \brief Critical Pressure for real fluid model.  */
  Density_Critical,   /*!< \brief Critical Density for real fluid model.  */
  Acentric_Factor,   /*!< \brief Acentric Factor for real fluid model.  */
  Mu_Constant,     /*!< \brief Constant viscosity for ConstantViscosity model.  */
  Mu_ConstantND,   /*!< \brief Non-dimensional constant viscosity for ConstantViscosity model.  */
  Kt_Constant,     /*!< \brief Constant thermal conductivity for ConstantConductivity model.  */
  Kt_ConstantND,   /*!< \brief Non-dimensional constant thermal conductivity for ConstantConductivity model.  */
  Mu_Ref,     /*!< \brief Reference viscosity for Sutherland model.  */
  Mu_RefND,   /*!< \brief Non-dimensional reference viscosity for Sutherland model.  */
  Mu_Temperature_Ref,     /*!< \brief Reference temperature for Sutherland model.  */
  Mu_Temperature_RefND,   /*!< \brief Non-dimensional reference temperature for Sutherland model.  */
  Mu_S,     /*!< \brief Reference S for Sutherland model.  */
  Mu_SND,   /*!< \brief Non-dimensional reference S for Sutherland model.  */
  *CpPolyCoefficients,   /*!< \brief Definition of the temperature polynomial coefficients for specific heat Cp. */
  *MuPolyCoefficients,   /*!< \brief Definition of the temperature polynomial coefficients for viscosity. */
  *KtPolyCoefficients,   /*!< \brief Definition of the temperature polynomial coefficients for thermal conductivity. */
  *CpPolyCoefficientsND,   /*!< \brief Definition of the non-dimensional temperature polynomial coefficients for specific heat Cp. */
  *MuPolyCoefficientsND,   /*!< \brief Definition of the non-dimensional temperature polynomial coefficients for viscosity. */
  *KtPolyCoefficientsND,   /*!< \brief Definition of the non-dimensional temperature polynomial coefficients for thermal conductivity. */
  Thermal_Conductivity_Solid, /*!< \brief Thermal conductivity in solids. */
  Thermal_Diffusivity_Solid, /*!< \brief Thermal diffusivity in solids. */
  Temperature_Freestream_Solid, /*!< \brief Temperature in solids at freestream conditions. */
  Density_Solid,      /*!< \brief Total density in solids. */  
  *Velocity_FreeStream,     /*!< \brief Free-stream velocity vector of the fluid.  */
  Energy_FreeStream,     /*!< \brief Free-stream total energy of the fluid.  */
  ModVel_FreeStream,     /*!< \brief Magnitude of the free-stream velocity of the fluid.  */
  ModVel_FreeStreamND,     /*!< \brief Non-dimensional magnitude of the free-stream velocity of the fluid.  */
  Density_FreeStream,     /*!< \brief Free-stream density of the fluid. */
  Viscosity_FreeStream,     /*!< \brief Free-stream viscosity of the fluid.  */
  Tke_FreeStream,     /*!< \brief Total turbulent kinetic energy of the fluid.  */
  Intermittency_FreeStream,     /*!< \brief Freestream intermittency (for sagt transition model) of the fluid.  */
  TurbulenceIntensity_FreeStream,     /*!< \brief Freestream turbulent intensity (for sagt transition model) of the fluid.  */
  Turb2LamViscRatio_FreeStream,          /*!< \brief Ratio of turbulent to laminar viscosity. */
  NuFactor_FreeStream,  /*!< \brief Ratio of turbulent to laminar viscosity. */
  NuFactor_Engine,  /*!< \brief Ratio of turbulent to laminar viscosity at the engine. */
  SecondaryFlow_ActDisk,  /*!< \brief Ratio of turbulent to laminar viscosity at the actuator disk. */
  Initial_BCThrust,  /*!< \brief Ratio of turbulent to laminar viscosity at the actuator disk. */
  Pressure_FreeStream,     /*!< \brief Total pressure of the fluid. */
  Pressure_Thermodynamic,     /*!< \brief Thermodynamic pressure of the fluid. */
  Temperature_FreeStream,  /*!< \brief Total temperature of the fluid.  */
  Temperature_ve_FreeStream,  /*!< \brief Total vibrational-electronic temperature of the fluid.  */
  *MassFrac_FreeStream, /*!< \brief Mixture mass fractions of the fluid. */
  Prandtl_Lam,      /*!< \brief Laminar Prandtl number for the gas.  */
  Prandtl_Turb,     /*!< \brief Turbulent Prandtl number for the gas.  */
  Length_Ref,       /*!< \brief Reference length for non-dimensionalization. */
  Pressure_Ref,     /*!< \brief Reference pressure for non-dimensionalization.  */
  Temperature_Ref,  /*!< \brief Reference temperature for non-dimensionalization.*/
  Density_Ref,      /*!< \brief Reference density for non-dimensionalization.*/
  Velocity_Ref,     /*!< \brief Reference velocity for non-dimensionalization.*/
  Time_Ref,                  /*!< \brief Reference time for non-dimensionalization. */
  Viscosity_Ref,              /*!< \brief Reference viscosity for non-dimensionalization. */
  Conductivity_Ref,           /*!< \brief Reference conductivity for non-dimensionalization. */
  Energy_Ref,                 /*!< \brief Reference viscosity for non-dimensionalization. */
  Wall_Temperature,           /*!< \brief Temperature at an isotropic wall in Kelvin. */
  Omega_Ref,                  /*!< \brief Reference angular velocity for non-dimensionalization. */
  Force_Ref,                  /*!< \brief Reference body force for non-dimensionalization. */
  Pressure_FreeStreamND,      /*!< \brief Farfield pressure value (external flow). */
  Pressure_ThermodynamicND,   /*!< \brief Farfield thermodynamic pressure value. */
  Temperature_FreeStreamND,   /*!< \brief Farfield temperature value (external flow). */
  Density_FreeStreamND,       /*!< \brief Farfield density value (external flow). */
  Velocity_FreeStreamND[3],   /*!< \brief Farfield velocity values (external flow). */
  Energy_FreeStreamND,        /*!< \brief Farfield energy value (external flow). */
  Viscosity_FreeStreamND,     /*!< \brief Farfield viscosity value (external flow). */
  Tke_FreeStreamND,           /*!< \brief Farfield kinetic energy (external flow). */
  Omega_FreeStreamND,         /*!< \brief Specific dissipation (external flow). */
  Omega_FreeStream;           /*!< \brief Specific dissipation (external flow). */
  unsigned short nElectric_Constant; /*!< \brief Number of different electric constants. */
  su2double *Electric_Constant;   /*!< \brief Dielectric constant modulus. */
  su2double Knowles_B,      /*!< \brief Knowles material model constant B. */
  Knowles_N;                /*!< \brief Knowles material model constant N. */
  bool DE_Effects; 						/*!< Application of DE effects to FE analysis */
  bool RefGeom; 						/*!< Read a reference geometry for optimization purposes. */
  unsigned long refNodeID;     /*!< \brief Global ID for the reference node (optimization). */
  string RefGeom_FEMFileName;    			/*!< \brief File name for reference geometry. */
  unsigned short RefGeom_FileFormat;	/*!< \brief Mesh input format. */
  unsigned short Kind_2DElasForm;			/*!< \brief Kind of bidimensional elasticity solver. */
  unsigned short nIterFSI;	  /*!< \brief Number of maximum number of subiterations in a FSI problem. */
  unsigned short nIterFSI_Ramp;  /*!< \brief Number of FSI subiterations during which a ramp is applied. */
  unsigned short iInst;       /*!< \brief Current instance value */
  su2double AitkenStatRelax;	/*!< \brief Aitken's relaxation factor (if set as static) */
  su2double AitkenDynMaxInit;	/*!< \brief Aitken's maximum dynamic relaxation factor for the first iteration */
  su2double AitkenDynMinInit;	/*!< \brief Aitken's minimum dynamic relaxation factor for the first iteration */
  bool RampAndRelease;        /*!< \brief option for ramp load and release */
  bool Sine_Load;             /*!< \brief option for sine load */
  su2double *SineLoad_Coeff;  /*!< \brief Stores the load coefficient */
  su2double Thermal_Diffusivity;			/*!< \brief Thermal diffusivity used in the heat solver. */
  su2double Cyclic_Pitch,     /*!< \brief Cyclic pitch for rotorcraft simulations. */
  Collective_Pitch;           /*!< \brief Collective pitch for rotorcraft simulations. */
  su2double Mach_Motion;			/*!< \brief Mach number based on mesh velocity and freestream quantities. */
  su2double *Motion_Origin_X, /*!< \brief X-coordinate of the mesh motion origin. */
  *Motion_Origin_Y,           /*!< \brief Y-coordinate of the mesh motion origin. */
  *Motion_Origin_Z,           /*!< \brief Z-coordinate of the mesh motion origin. */
  *Translation_Rate_X,        /*!< \brief Translational velocity of the mesh in the x-direction. */
  *Translation_Rate_Y,        /*!< \brief Translational velocity of the mesh in the y-direction. */
  *Translation_Rate_Z,        /*!< \brief Translational velocity of the mesh in the z-direction. */
  *Rotation_Rate_X,           /*!< \brief Angular velocity of the mesh about the x-axis. */
  *Rotation_Rate_Y,           /*!< \brief Angular velocity of the mesh about the y-axis. */
  *Rotation_Rate_Z,           /*!< \brief Angular velocity of the mesh about the z-axis. */
  *Pitching_Omega_X,          /*!< \brief Angular frequency of the mesh pitching about the x-axis. */
  *Pitching_Omega_Y,          /*!< \brief Angular frequency of the mesh pitching about the y-axis. */
  *Pitching_Omega_Z,          /*!< \brief Angular frequency of the mesh pitching about the z-axis. */
  *Pitching_Ampl_X,           /*!< \brief Pitching amplitude about the x-axis. */
  *Pitching_Ampl_Y,           /*!< \brief Pitching amplitude about the y-axis. */
  *Pitching_Ampl_Z,           /*!< \brief Pitching amplitude about the z-axis. */
  *Pitching_Phase_X,          /*!< \brief Pitching phase offset about the x-axis. */
  *Pitching_Phase_Y,          /*!< \brief Pitching phase offset about the y-axis. */
  *Pitching_Phase_Z,          /*!< \brief Pitching phase offset about the z-axis. */
  *Plunging_Omega_X,          /*!< \brief Angular frequency of the mesh plunging in the x-direction. */
  *Plunging_Omega_Y,          /*!< \brief Angular frequency of the mesh plunging in the y-direction. */
  *Plunging_Omega_Z,          /*!< \brief Angular frequency of the mesh plunging in the z-direction. */
  *Plunging_Ampl_X,           /*!< \brief Plunging amplitude in the x-direction. */
  *Plunging_Ampl_Y,           /*!< \brief Plunging amplitude in the y-direction. */
  *Plunging_Ampl_Z,           /*!< \brief Plunging amplitude in the z-direction. */
  *Omega_HB;                  /*!< \brief Frequency for Harmonic Balance Operator (in rad/s). */
  unsigned short nMotion_Origin_X,    /*!< \brief Number of X-coordinate mesh motion origins. */
  nMotion_Origin_Y,           /*!< \brief Number of Y-coordinate mesh motion origins. */
  nMotion_Origin_Z,           /*!< \brief Number of Z-coordinate mesh motion origins. */
  nTranslation_Rate_X,        /*!< \brief Number of Translational x-velocities for mesh motion. */
  nTranslation_Rate_Y,        /*!< \brief Number of Translational y-velocities for mesh motion. */
  nTranslation_Rate_Z,        /*!< \brief Number of Translational z-velocities for mesh motion. */
  nRotation_Rate_X,           /*!< \brief Number of Angular velocities about the x-axis for mesh motion. */
  nRotation_Rate_Y,           /*!< \brief Number of Angular velocities about the y-axis for mesh motion. */
  nRotation_Rate_Z,           /*!< \brief Number of Angular velocities about the z-axis for mesh motion. */
  nPitching_Omega_X,          /*!< \brief Number of Angular frequencies about the x-axis for pitching. */
  nPitching_Omega_Y,          /*!< \brief Number of Angular frequencies about the y-axis for pitching. */
  nPitching_Omega_Z,          /*!< \brief Number of Angular frequencies about the z-axis for pitching. */
  nPitching_Ampl_X,           /*!< \brief Number of Pitching amplitudes about the x-axis. */
  nPitching_Ampl_Y,           /*!< \brief Number of Pitching amplitudes about the y-axis. */
  nPitching_Ampl_Z,           /*!< \brief Number of Pitching amplitudes about the z-axis. */
  nPitching_Phase_X,          /*!< \brief Number of Pitching phase offsets about the x-axis. */
  nPitching_Phase_Y,          /*!< \brief Number of Pitching phase offsets about the y-axis. */
  nPitching_Phase_Z,          /*!< \brief Number of Pitching phase offsets about the z-axis. */
  nPlunging_Omega_X,          /*!< \brief Number of Angular frequencies in the x-direction for plunging. */
  nPlunging_Omega_Y,          /*!< \brief Number of Angular frequencies in the y-direction for plunging. */
  nPlunging_Omega_Z,          /*!< \brief Number of Angular frequencies in the z-direction for plunging. */
  nPlunging_Ampl_X,           /*!< \brief Number of Plunging amplitudes in the x-direction. */
  nPlunging_Ampl_Y,           /*!< \brief Number of Plunging amplitudes in the y-direction. */
  nPlunging_Ampl_Z,           /*!< \brief Number of Plunging amplitudes in the z-direction. */
  nOmega_HB,                /*!< \brief Number of frequencies in Harmonic Balance Operator. */
  nMoveMotion_Origin,         /*!< \brief Number of motion origins. */
  *MoveMotion_Origin;         /*!< \brief Keeps track if we should move moment origin. */
  vector<vector<vector<su2double> > > Aeroelastic_np1, /*!< \brief Aeroelastic solution at time level n+1. */
  Aeroelastic_n,              /*!< \brief Aeroelastic solution at time level n. */
  Aeroelastic_n1;             /*!< \brief Aeroelastic solution at time level n-1. */
  su2double FlutterSpeedIndex,/*!< \brief The flutter speed index. */
  PlungeNaturalFrequency,     /*!< \brief Plunging natural frequency for Aeroelastic. */
  PitchNaturalFrequency,      /*!< \brief Pitch natural frequency for Aeroelastic. */
  AirfoilMassRatio,           /*!< \brief The airfoil mass ratio for Aeroelastic. */
  CG_Location,                /*!< \brief Center of gravity location for Aeroelastic. */
  RadiusGyrationSquared;      /*!< \brief The radius of gyration squared for Aeroelastic. */
  su2double *Aeroelastic_plunge, /*!< \brief Value of plunging coordinate at the end of an external iteration. */
  *Aeroelastic_pitch;         /*!< \brief Value of pitching coordinate at the end of an external iteration. */
  unsigned short AeroelasticIter; /*!< \brief Solve the aeroelastic equations every given number of internal iterations. */
  unsigned short Gust_Type,	  /*!< \brief Type of Gust. */
  Gust_Dir;                   /*!< \brief Direction of the gust */
  su2double Gust_WaveLength,  /*!< \brief The gust wavelength. */
  Gust_Periods,               /*!< \brief Number of gust periods. */
  Gust_Ampl,                  /*!< \brief Gust amplitude. */
  Gust_Begin_Time,            /*!< \brief Time at which to begin the gust. */
  Gust_Begin_Loc;             /*!< \brief Location at which the gust begins. */
  long Visualize_CV;          /*!< \brief Node number for the CV to be visualized */
  bool ExtraOutput;
  bool Wall_Functions;         /*!< \brief Use wall functions with the turbulence model */
  long ExtraHeatOutputZone;   /*!< \brief Heat solver zone with extra screen output */
  bool DeadLoad; 	          	/*!< Application of dead loads to the FE analysis */
  bool PseudoStatic;    /*!< Application of dead loads to the FE analysis */
  bool SteadyRestart; 	      /*!< Restart from a steady state for FSI problems. */
  su2double Newmark_beta,		/*!< \brief Parameter alpha for Newmark method. */
  Newmark_gamma;				      /*!< \brief Parameter delta for Newmark method. */
  unsigned short nIntCoeffs;	/*!< \brief Number of integration coeffs for structural calculations. */
  su2double *Int_Coeffs;		  /*!< \brief Time integration coefficients for structural method. */
  unsigned short nElasticityMod,  /*!< \brief Number of different values for the elasticity modulus. */
  nPoissonRatio,                    /*!< \brief Number of different values for the Poisson ratio modulus. */
  nMaterialDensity;                 /*!< \brief Number of different values for the Material density. */
  su2double *ElasticityMod,         /*!< \brief Value of the elasticity moduli. */
  *PoissonRatio,                    /*!< \brief Value of the Poisson ratios. */
  *MaterialDensity;                 /*!< \brief Value of the Material densities. */
  unsigned short nElectric_Field,	/*!< \brief Number of different values for the electric field in the membrane. */
  nDim_Electric_Field;				/*!< \brief Dimensionality of the problem. */
  unsigned short nDim_RefNode;   /*!< \brief Dimensionality of the vector . */
  su2double *Electric_Field_Mod, 	/*!< \brief Values of the modulus of the electric field. */
  *Electric_Field_Dir;				/*!< \brief Direction of the electric field. */
  su2double *RefNode_Displacement;  /*!< \brief Displacement of the reference node. */
  bool Ramp_Load;				          /*!< \brief Apply the load with linear increases. */
  unsigned short Dynamic_LoadTransfer;  /*!< \brief Method for dynamic load transferring. */
  bool IncrementalLoad;		    /*!< \brief Apply the load in increments (for nonlinear structural analysis). */
  unsigned long IncLoad_Nincrements; /*!< \brief Number of increments. */
  su2double *IncLoad_Criteria;/*!< \brief Criteria for the application of incremental loading. */
  su2double Ramp_Time;			  /*!< \brief Time until the maximum load is applied. */
  bool Predictor,             /*!< \brief Determines whether a predictor step is used. */
  Relaxation;                 /*!< \brief Determines whether a relaxation step is used. */
  unsigned short Pred_Order;  /*!< \brief Order of the predictor for FSI applications. */
  unsigned short Kind_Interpolation; /*!\brief type of interpolation to use for FSI applications. */
  bool ConservativeInterpolation; /*!\brief Conservative approach for non matching mesh interpolation. */
  unsigned short Kind_RadialBasisFunction; /*!\brief type of radial basis function to use for radial basis FSI. */
  bool RadialBasisFunction_PolynomialOption; /*!\brief Option of whether to include polynomial terms in Radial Basis Function Interpolation or not. */
  su2double RadialBasisFunction_Parameter; /*!\brief Radial basis function parameter. */
  bool Prestretch;            /*!< Read a reference geometry for optimization purposes. */
  string Prestretch_FEMFileName;         /*!< \brief File name for reference geometry. */
  string FEA_FileName;         /*!< \brief File name for element-based properties. */
  su2double RefGeom_Penalty,        /*!< \brief Penalty weight value for the reference geometry objective function. */
  RefNode_Penalty,            /*!< \brief Penalty weight value for the reference node objective function. */
  DV_Penalty;                 /*!< \brief Penalty weight to add a constraint to the total amount of stiffness. */
  bool addCrossTerm;          /*!< \brief Evaluates the need to add the cross term when setting the adjoint output. */
  unsigned long Nonphys_Points, /*!< \brief Current number of non-physical points in the solution. */
  Nonphys_Reconstr;      /*!< \brief Current number of non-physical reconstructions for 2nd-order upwinding. */
  bool ParMETIS;      /*!< \brief Boolean for activating ParMETIS mode (while testing). */
  unsigned short DirectDiff; /*!< \brief Direct Differentation mode. */
  bool DiscreteAdjoint; /*!< \brief AD-based discrete adjoint mode. */
  unsigned long Wrt_Surf_Freq_DualTime;	/*!< \brief Writing surface solution frequency for Dual Time. */
  su2double Const_DES;   /*!< \brief Detached Eddy Simulation Constant. */
  unsigned short Kind_HybridRANSLES; /*!< \brief Kind of Hybrid RANS/LES. */
  unsigned short Kind_RoeLowDiss;    /*!< \brief Kind of Roe scheme with low dissipation for unsteady flows. */
  bool QCR;                   /*!< \brief Spalart-Allmaras with Quadratic Constitutive Relation, 2000 version (SA-QCR2000) . */
  su2double *default_vel_inf, /*!< \brief Default freestream velocity array for the COption class. */
  *default_eng_cyl,           /*!< \brief Default engine box array for the COption class. */
  *default_eng_val,           /*!< \brief Default engine box array values for the COption class. */
  *default_cfl_adapt,         /*!< \brief Default CFL adapt param array for the COption class. */
  *default_jst_coeff,         /*!< \brief Default artificial dissipation (flow) array for the COption class. */
  *default_ffd_coeff,         /*!< \brief Default artificial dissipation (flow) array for the COption class. */
  *default_mixedout_coeff,    /*!< \brief Default default mixedout algorithm coefficients for the COption class. */
  *default_rampRotFrame_coeff,/*!< \brief Default ramp rotating frame coefficients for the COption class. */
  *default_rampOutPres_coeff, /*!< \brief Default ramp outlet pressure coefficients for the COption class. */
  *default_jst_adj_coeff,      /*!< \brief Default artificial dissipation (adjoint) array for the COption class. */
  *default_ad_coeff_heat,     /*!< \brief Default artificial dissipation (heat) array for the COption class. */  
  *default_obj_coeff,         /*!< \brief Default objective array for the COption class. */
  *default_geo_loc,           /*!< \brief Default SU2_GEO section locations array for the COption class. */
  *default_distortion,        /*!< \brief Default SU2_GEO section locations array for the COption class. */
  *default_ea_lim,            /*!< \brief Default equivalent area limit array for the COption class. */
  *default_grid_fix,          /*!< \brief Default fixed grid (non-deforming region) array for the COption class. */
  *default_htp_axis,          /*!< \brief Default HTP axis for the COption class. */
  *default_ffd_axis,          /*!< \brief Default FFD axis for the COption class. */
  *default_inc_crit,          /*!< \brief Default incremental criteria array for the COption class. */
  *default_extrarelfac,       /*!< \brief Default extra relaxation factor for Giles BC in the COption class. */
  *default_sineload_coeff;    /*!< \brief Default values for a sine load. */
  unsigned short nSpanWiseSections; /*!< \brief number of span-wise sections */
  unsigned short nSpanMaxAllZones; /*!< \brief number of maximum span-wise sections for all zones */
  unsigned short *nSpan_iZones;  /*!< \brief number of span-wise sections for each zones */
  bool turbMixingPlane;   /*!< \brief option for turbulent mixingplane */
  bool SpatialFourier; /*!< \brief option for computing the fourier transforms for subsonic non-reflecting BC. */
  bool RampRotatingFrame;   /*!< \brief option for ramping up or down the Rotating Frame values */
  bool RampOutletPressure;  /*!< \brief option for ramping up or down the outlet pressure */
  su2double *Mixedout_Coeff; /*!< \brief coefficient for the  */
  su2double *RampRotatingFrame_Coeff; /*!< \brief coefficient for Rotating frame ramp */
  su2double *RampOutletPressure_Coeff; /*!< \brief coefficient for outlet pressure ramp */
  su2double AverageMachLimit;       /*!< \brief option for turbulent mixingplane */
  su2double *FinalRotation_Rate_Z; /*!< \brief Final rotation rate Z if Ramp rotating frame is activated. */
  su2double FinalOutletPressure; /*!< \brief Final outlet pressure if Ramp outlet pressure is activated. */
  su2double MonitorOutletPressure; /*!< \brief Monitor outlet pressure if Ramp outlet pressure is activated. */
  su2double *default_body_force;        /*!< \brief Default body force vector for the COption class. */
  su2double *default_nacelle_location;        /*!< \brief Location of the nacelle. */
  su2double *default_cp_polycoeffs;        /*!< \brief Array for specific heat polynomial coefficients. */
  su2double *default_mu_polycoeffs;        /*!< \brief Array for viscosity polynomial coefficients. */
  su2double *default_kt_polycoeffs;        /*!< \brief Array for thermal conductivity polynomial coefficients. */
  su2double *ExtraRelFacGiles; /*!< \brief coefficient for extra relaxation factor for Giles BC*/
  bool Body_Force;            /*!< \brief Flag to know if a body force is included in the formulation. */
  su2double *Body_Force_Vector;  /*!< \brief Values of the prescribed body force vector. */
  su2double *FreeStreamTurboNormal; /*!< \brief Direction to initialize the flow in turbomachinery computation */
  su2double Restart_Bandwidth_Agg; /*!< \brief The aggregate of the bandwidth for writing binary restarts (to be averaged later). */
  su2double Max_Vel2; /*!< \brief The maximum velocity^2 in the domain for the incompressible preconditioner. */
  bool topology_optimization; /*!< \brief If the structural solver should consider a variable density field to penalize element stiffness. */
  string top_optim_output_file; /*!< \brief File to where the derivatives w.r.t. element densities will be written to. */
  su2double simp_exponent; /*!< \brief Exponent for the density-based stiffness penalization of the SIMP method. */
  su2double simp_minimum_stiffness; /*!< \brief Lower bound for the stiffness penalization of the SIMP method. */
  unsigned short top_optim_nKernel, /*!< \brief Number of kernels specified. */
                *top_optim_kernels, /*!< \brief The kernels to use. */
                 top_optim_nKernelParams, /*!< \brief Number of kernel parameters specified. */
                 top_optim_nRadius; /*!< \brief Number of radius values specified. */
  su2double *top_optim_kernel_params, /*!< \brief The kernel parameters. */
            *top_optim_filter_radius; /*!< \brief Radius of the filter(s) used on the design density for topology optimization. */
  unsigned short top_optim_proj_type; /*!< \brief The projection function used in topology optimization. */
  su2double top_optim_proj_param;  /*!< \brief The value of the parameter for the projection function. */

  unsigned short Riemann_Solver_FEM;         /*!< \brief Riemann solver chosen for the DG method. */
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
  unsigned short Kind_Verification_Solution;  /*!< \brief Verification solution for accuracy assessment. */

  ofstream *ConvHistFile;       /*!< \brief Store the pointer to each history file */
  bool Time_Domain;             /*!< \brief Determines if the multizone problem is solved in time-domain */
  unsigned long Outer_Iter,    /*!< \brief Determines the number of outer iterations in the multizone problem */
  Inner_Iter,                   /*!< \brief Determines the number of inner iterations in each multizone block */
  Time_Iter,                    /*!< \brief Determines the number of time iterations in the multizone problem */
  Iter,                         /*!< \brief Determines the number of pseudo-time iterations in a single-zone problem */
  Restart_Iter;                 /*!< \brief Determines the restart iteration in the multizone problem */
  su2double Time_Step;          /*!< \brief Determines the time step for the multizone problem */
  su2double Max_Time;           /*!< \brief Determines the maximum time for the time-domain problems */
  bool Multizone_Mesh;          /*!< \brief Determines if the mesh contains multiple zones. */
  bool SinglezoneDriver;        /*!< \brief Determines if the single-zone driver is used. (TEMPORARY) */
  bool SpecialOutput,           /*!< \brief Determines if the special output is written. */
  Wrt_ForcesBreakdown;          /*!< \brief Determines if the forces breakdown file is written. */
  bool Multizone_Residual;      /*!< \brief Determines if memory should be allocated for the multizone residual. */
  
  bool using_uq;                /*!< \brief Using uncertainty quantification with SST model */
  su2double uq_delta_b;            /*!< \brief Parameter used to perturb eigenvalues of Reynolds Stress Matrix */
  unsigned short eig_val_comp;  /*!< \brief Parameter used to determine type of eigenvalue perturbation */
  su2double uq_urlx;            /*!< \brief Under-relaxation factor */
  bool uq_permute;              /*!< \brief Permutation of eigenvectors */

  /*--- all_options is a map containing all of the options. This is used during config file parsing
   to track the options which have not been set (so the default values can be used). Without this map
   there would be no list of all the config file options. ---*/
  
  map<string, bool> all_options;
  
  /*--- brief param is a map from the option name (config file string) to its decoder (the specific child
   class of COptionBase that turns the string into a value) ---*/
  
  map<string, COptionBase*> option_map;
  
  
  // All of the addXxxOptions take in the name of the option, and a refernce to the field of that option
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
  
  /*!<\brief addDoubleOption creates a config file parser for an option with the given name whose
   value can be represented by a su2double.*/
  
  void addDoubleOption(const string name, su2double & option_field, su2double default_value) {
    // Check if the key is already in the map. If this fails, it is coder error
    // and not user error, so throw.
    assert(option_map.find(name) == option_map.end());
    
    // Add this option to the list of all the options
    all_options.insert(pair<string, bool>(name, true));
    
    // Create the parser for a su2double option with a reference to the option_field and the desired
    // default value. This will take the string in the config file, convert it to a su2double, and
    // place that su2double in the memory location specified by the reference.
    COptionBase* val = new COptionDouble(name, option_field, default_value);
    
    // Create an association between the option name ("CFL") and the parser generated above.
    // During configuration, the parsing script will get the option name, and use this map
    // to find how to parse that option.
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addStringOption(const string name, string & option_field, string default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionString(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addIntegerOption(const string name, int & option_field, int default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionInt(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addUnsignedLongOption(const string name, unsigned long & option_field, unsigned long default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionULong(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addUnsignedShortOption(const string name, unsigned short & option_field, unsigned short default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionUShort(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addLongOption(const string name, long & option_field, long default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionLong(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addBoolOption(const string name, bool & option_field, bool default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionBool(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  // enum types work differently than all of the others because there are a small number of valid
  // string entries for the type. One must also provide a list of all the valid strings of that type.
  template <class Tenum>
  void addEnumOption(const string name, unsigned short & option_field, const map<string, Tenum> & enum_map, Tenum default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionEnum<Tenum>(name, enum_map, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
    return;
  }
  
  
  // input_size is the number of options read in from the config file
  template <class Tenum>
  void addEnumListOption(const string name, unsigned short & input_size, unsigned short * & option_field, const map<string, Tenum> & enum_map) {
    input_size = 0;
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionEnumList<Tenum>(name, enum_map, option_field, input_size);
    option_map.insert( pair<string, COptionBase*>(name, val) );
  }
  
  void addDoubleArrayOption(const string name, const int size, su2double * & option_field, su2double * default_value) {
    
    //  su2double * def = new su2double [size];
    //  for (int i = 0; i < size; i++) {
    //    def[i] = default_value[i];
    //  }
    
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionDoubleArray(name, size, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addDoubleListOption(const string name, unsigned short & size, su2double * & option_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionDoubleList(name, size, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addShortListOption(const string name, unsigned short & size, short * & option_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionShortList(name, size, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addUShortListOption(const string name, unsigned short & size, unsigned short * & option_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionUShortList(name, size, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addStringListOption(const string name, unsigned short & num_marker, string* & option_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionStringList(name, num_marker, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addConvectOption(const string name, unsigned short & space_field, unsigned short & centered_field, unsigned short & upwind_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionConvect(name, space_field, centered_field, upwind_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addConvectFEMOption(const string name, unsigned short & space_field, unsigned short & fem_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionFEMConvect(name, space_field, fem_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addMathProblemOption(const string name, bool & ContinuousAdjoint, const bool & ContinuousAdjoint_default,
                            bool & DiscreteAdjoint, const bool & DiscreteAdjoint_default,
                            bool & Restart_Flow, const bool & Restart_Flow_default) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionMathProblem(name, ContinuousAdjoint, ContinuousAdjoint_default, DiscreteAdjoint, DiscreteAdjoint_default, Restart_Flow, Restart_Flow_default);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addDVParamOption(const string name, unsigned short & nDV_field, su2double** & paramDV, string* & FFDTag,
                        unsigned short* & design_variable) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionDVParam(name, nDV_field, paramDV, FFDTag, design_variable);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addDVValueOption(const string name, unsigned short* & nDVValue_field, su2double** & valueDV, unsigned short & nDV_field,  su2double** & paramDV,
                        unsigned short* & design_variable) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionDVValue(name, nDVValue_field, valueDV, nDV_field, paramDV, design_variable);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addFFDDefOption(const string name, unsigned short & nFFD_field, su2double** & coordFFD, string* & FFDTag) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionFFDDef(name, nFFD_field, coordFFD, FFDTag);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addFFDDegreeOption(const string name, unsigned short & nFFD_field, unsigned short** & degreeFFD) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionFFDDegree(name, nFFD_field, degreeFFD);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addStringDoubleListOption(const string name, unsigned short & list_size, string * & string_field,
                                 su2double* & double_field) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionStringDoubleList(name, list_size, string_field, double_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addInletOption(const string name, unsigned short & nMarker_Inlet, string * & Marker_Inlet,
                      su2double* & Ttotal, su2double* & Ptotal, su2double** & FlowDir) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionInlet(name, nMarker_Inlet, Marker_Inlet, Ttotal, Ptotal, FlowDir);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  template <class Tenum>
  void addRiemannOption(const string name, unsigned short & nMarker_Riemann, string * & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                        su2double* & var1, su2double* & var2, su2double** & FlowDir) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionRiemann<Tenum>(name, nMarker_Riemann, Marker_Riemann, option_field, enum_map, var1, var2, FlowDir);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  template <class Tenum>
  void addGilesOption(const string name, unsigned short & nMarker_Giles, string * & Marker_Giles, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                     su2double* & var1, su2double* & var2, su2double** & FlowDir, su2double* & relaxfactor1, su2double* & relaxfactor2) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionGiles<Tenum>(name, nMarker_Giles, Marker_Giles, option_field, enum_map, var1, var2, FlowDir, relaxfactor1, relaxfactor2);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addExhaustOption(const string name, unsigned short & nMarker_Exhaust, string * & Marker_Exhaust,
                        su2double* & Ttotal, su2double* & Ptotal) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionExhaust(name, nMarker_Exhaust, Marker_Exhaust, Ttotal, Ptotal);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addPeriodicOption(const string & name, unsigned short & nMarker_PerBound,
                         string* & Marker_PerBound, string* & Marker_PerDonor,
                         su2double** & RotCenter, su2double** & RotAngles, su2double** & Translation) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionPeriodic(name, nMarker_PerBound, Marker_PerBound, Marker_PerDonor, RotCenter, RotAngles, Translation);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
 
  void addTurboPerfOption(const string & name, unsigned short & nMarker_TurboPerf,
                    string* & Marker_TurboBoundIn, string* & Marker_TurboBoundOut) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionTurboPerformance(name, nMarker_TurboPerf, Marker_TurboBoundIn, Marker_TurboBoundOut);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addActDiskOption(const string & name,
                        unsigned short & nMarker_ActDiskInlet, unsigned short & nMarker_ActDiskOutlet, string* & Marker_ActDiskInlet, string* & Marker_ActDiskOutlet,
                        su2double** & ActDisk_PressJump, su2double** & ActDisk_TempJump, su2double** & ActDisk_Omega) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionActDisk(name,
                                          nMarker_ActDiskInlet, nMarker_ActDiskOutlet, Marker_ActDiskInlet, Marker_ActDiskOutlet,
                                          ActDisk_PressJump, ActDisk_TempJump, ActDisk_Omega);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addWallFunctionOption(const string &name,               unsigned short &list_size,
                             string* &string_field,            unsigned short* &val_Kind_WF,
                             unsigned short** &val_IntInfo_WF, su2double** &val_DoubleInfo_WF) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionWallFunction(name, list_size, string_field, val_Kind_WF,
                                               val_IntInfo_WF, val_DoubleInfo_WF);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addPythonOption(const string name) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionPython(name);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
public:
  
  vector<string> fields; /*!< \brief Tags for the different fields in a restart file. */
  
  /*!
   * \brief Constructor of the class which reads the input file.
   */
  CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_nDim, bool verb_high);
  
  /*!
   * \brief Constructor of the class which reads the input file.
   */
  CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software);
  
  /*!
   * \brief Constructor of the class which reads the input file.
   */
  CConfig(char case_filename[MAX_STRING_SIZE], CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CConfig(void);
  
  /*!
   * \brief Get the MPI communicator of SU2.
   * \return MPI communicator of SU2.
   */
  SU2_MPI::Comm GetMPICommunicator();

  /*!
   * \brief Set the MPI communicator for SU2.
   * \param[in] Communicator - MPI communicator for SU2.
   */
  void SetMPICommunicator(SU2_MPI::Comm Communicator);

  /*!
   * \brief Gets the number of zones in the mesh file.
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_format - Format of the file with the grid information.
   * \param[in] config - Definition of the particular problem.
   * \return Total number of zones in the grid file.
   */
  static unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config);
  
  /*!
   * \brief Gets the number of dimensions in the mesh file
   * \param[in] val_mesh_filename - Name of the file with the grid information.
   * \param[in] val_format - Format of the file with the grid information.
   * \return Total number of domains in the grid file.
   */
  static unsigned short GetnDim(string val_mesh_filename, unsigned short val_format);
  
  /*!
   * \brief Initializes pointers to null
   */
  void SetPointersNull(void);
  
  /*!
   * \brief breaks an input line from the config file into a set of tokens
   * \param[in] str - the input line string
   * \param[out] option_name - the name of the option found at the beginning of the line
   * \param[out] option_value - the tokens found after the "=" sign on the line
   * \returns false if the line is empty or a commment, true otherwise
   */
  bool TokenizeString(string & str, string & option_name,
                      vector<string> & option_value);
  
  /*!
   * \brief Get reference origin for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin (in cartesians coordinates) for moment computation.
   */
  su2double *GetRefOriginMoment(unsigned short val_marker);
  
  /*!
   * \brief Get reference origin x-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin x-coordinate (in cartesians coordinates) for moment computation.
   */
  su2double GetRefOriginMoment_X(unsigned short val_marker);
  
  /*!
   * \brief Get reference origin y-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin y-coordinate (in cartesians coordinates) for moment computation.
   */
  su2double GetRefOriginMoment_Y(unsigned short val_marker);
  
  /*!
   * \brief Get reference origin z-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Reference origin z-coordinate (in cartesians coordinates) for moment computation.
   */
  su2double GetRefOriginMoment_Z(unsigned short val_marker);
  
  /*!
   * \brief Set reference origin x-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val_origin - New x-coordinate of the mesh motion origin.
   */
  void SetRefOriginMoment_X(unsigned short val_marker, su2double val_origin);
  
  /*!
   * \brief Set reference origin y-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val_origin - New y-coordinate of the mesh motion origin.
   */
  void SetRefOriginMoment_Y(unsigned short val_marker, su2double val_origin);
  
  /*!
   * \brief Set reference origin z-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val_origin - New z-coordinate of the mesh motion origin.
   */
  void SetRefOriginMoment_Z(unsigned short val_marker, su2double val_origin);
  
  /*!
   * \brief Get index of the upper and lower horizontal plane.
   * \param[in] index - 0 means upper surface, and 1 means lower surface.
   * \return Index of the upper and lower surface.
   */
  string GetPlaneTag(unsigned short index);
  
  /*!
   * \brief Get the integration limits for the equivalent area computation.
   * \param[in] index - 0 means x_min, and 1 means x_max.
   * \return Integration limits for the equivalent area computation.
   */
  su2double GetEA_IntLimit(unsigned short index);
  
  /*!
   * \brief Get the integration limits for the equivalent area computation.
   * \param[in] index - 0 means x_min, and 1 means x_max.
   * \return Integration limits for the equivalent area computation.
   */
  su2double GetEA_ScaleFactor(void);
  
  /*!
   * \brief Get the limit value for the adjoint variables.
   * \return Limit value for the adjoint variables.
   */
  su2double GetAdjointLimit(void);
  
  /*!
   * \brief Get the the coordinates where of the box where the grid is going to be deformed.
   * \return Coordinates where of the box where the grid is going to be deformed.
   */
  su2double *GetHold_GridFixed_Coord(void);
  
  /*!
   * \brief Get the the coordinates where of the box where a subsonic region is imposed.
   * \return Coordinates where of the box where the grid is going to be a subsonic region.
   */
  su2double *GetSubsonicEngine_Values(void);
  
  /*!
   * \brief Get the the coordinates where of the box where a subsonic region is imposed.
   * \return Coordinates where of the box where the grid is going to be a subsonic region.
   */
  su2double *GetSubsonicEngine_Cyl(void);
  
  /*!
   * \brief Get the the coordinates where of the box where a subsonic region is imposed.
   * \return Coordinates where of the box where the grid is going to be a subsonic region.
   */
  su2double *GetDistortionRack(void);
  
  /*!
   * \brief Get the power of the dual volume in the grid adaptation sensor.
   * \return Power of the dual volume in the grid adaptation sensor.
   */
  su2double GetDualVol_Power(void);
  
  /*!
   * \brief Get Information about if there is an analytical definition of the surface for doing the
   *        grid adaptation.
   * \return Definition of the surfaces. NONE implies that there isn't any analytical definition
   *         and it will use and interpolation.
   */
  unsigned short GetAnalytical_Surface(void);
  
  /*!
   * \brief Get Description of the geometry to be analyzed
   */
  unsigned short GetGeo_Description(void);
  
  /*!
   * \brief Creates a tecplot file to visualize the partition made by the DDC software.
   * \return <code>TRUE</code> if the partition is going to be plotted; otherwise <code>FALSE</code>.
   */
  bool GetExtraOutput(void);

  /*!
   * \brief Heat solver zone with extra screen output.
   * \return Heat solver zone with extra screen output.
   */
  long GetExtraHeatOutputZone(void);
  
  /*!
   * \brief Get the value of the Mach number (velocity divided by speed of sound).
   * \return Value of the Mach number.
   */
  su2double GetMach(void);
  
  /*!
   * \brief Get the value of the Gamma of fluid (ratio of specific heats).
   * \return Value of the constant: Gamma
   */
  su2double GetGamma(void);
  
  /*!
   * \brief Get the values of the CFL adapation.
   * \return Value of CFL adapation
   */
  su2double GetCFL_AdaptParam(unsigned short val_index);
  
  /*!
   * \brief Get the values of the CFL adapation.
   * \return Value of CFL adapation
   */
  bool GetCFL_Adapt(void);
  
  /*!
   * \brief Get the values of the CFL adapation.
   * \return Value of CFL adapation
   */
  su2double GetHTP_Axis(unsigned short val_index);
  
  /*!
   * \brief Get the value of the limits for the sections.
   * \return Value of the limits for the sections.
   */
  su2double GetStations_Bounds(unsigned short val_var);
  
  /*!
   * \brief Get the value of the vector that connects the cartesian axis with a sherical or cylindrical one.
   * \return Coordinate of the Axis.
   */
  su2double GetFFD_Axis(unsigned short val_var);
  
  /*!
   * \brief Get the value of the bulk modulus.
   * \return Value of the bulk modulus.
   */
  su2double GetBulk_Modulus(void);
  
  /*!
   * \brief Get the epsilon^2 multiplier for Beta in the incompressible preconditioner.
   * \return Value of the epsilon^2 multiplier for Beta in the incompressible preconditioner.
   */
  su2double GetBeta_Factor(void);
  
  /*!
   * \brief Get the value of specific gas constant.
   * \return Value of the constant: Gamma
   */
  su2double GetGas_Constant(void);
  
  /*!
   * \brief Get the value of specific gas constant.
   * \return Value of the constant: Gamma
   */
  su2double GetGas_ConstantND(void);
  
  /*!
   * \brief Get the value of the molecular weight for an incompressible ideal gas (g/mol).
   * \return Value of the molecular weight for an incompressible ideal gas (g/mol).
   */
  su2double GetMolecular_Weight(void);
  
  /*!
   * \brief Get the value of specific heat at constant pressure.
   * \return Value of the constant: Cp
   */
  su2double GetSpecific_Heat_Cp(void);

  /*!
   * \brief Get the value of the specific heat for solids.
   * \return Specific heat number (solid).
   */
  su2double GetSpecific_Heat_Cp_Solid(void);
  
  /*!
   * \brief Get the non-dimensional value of specific heat at constant pressure.
   * \return Value of the non-dim. constant: Cp
   */
  su2double GetSpecific_Heat_CpND(void);

  /*!
   * \brief Get the value of specific heat at constant volume.
   * \return Value of the constant: Cv
   */
  su2double GetSpecific_Heat_Cv(void);
  
  /*!
   * \brief Get the non-dimensional value of specific heat at constant volume.
   * \return Value of the non-dim. constant: Cv
   */
  su2double GetSpecific_Heat_CvND(void);

  /*!
   * \brief Get the coefficients of the Blottner viscosity model
   * \param[in] val_Species - Index of the species
   * \param[in] val_Coeff - Index of the coefficient (As, Bs, Cs)
   * \return Value of the Blottner coefficient
   */
  su2double GetBlottnerCoeff(unsigned short val_Species, unsigned short val_Coeff);
  
  /*!
   * \brief Get the p-norm for heat-flux objective functions (adjoint problem).
   * \return Value of the heat flux p-norm
   */
  su2double GetPnormHeat(void);
  
  /*!
   * \brief Get the value of wall temperature.
   * \return Value of the constant: Temperature
   */
  su2double GetWallTemperature(void);
  
  /*!
   * \brief Get the reference value for the specific gas constant.
   * \return Reference value for the specific gas constant.
   */
  su2double GetGas_Constant_Ref(void);
  
  /*!
   * \brief Get the reference value for the heat flux.
   * \return Reference value for the heat flux.
   */
  su2double GetHeat_Flux_Ref(void);

  /*!
   * \brief Get the value of the frestream temperature.
   * \return Freestream temperature.
   */
  su2double GetTemperature_FreeStream(void);
  
  /*!
   * \brief Get the value of the frestream temperature.
   * \return Freestream temperature.
   */
  su2double GetEnergy_FreeStream(void);
  
  /*!
   * \brief Get the value of the frestream temperature.
   * \return Freestream temperature.
   */
  su2double GetViscosity_FreeStream(void);
  
  /*!
   * \brief Get the value of the frestream temperature.
   * \return Freestream temperature.
   */
  su2double GetDensity_FreeStream(void);

  /*!
   * \brief Get the value of the solid density.
   * \return Solid density.
   */
  su2double GetDensity_Solid(void);
  
  /*!
   * \brief Get the value of the frestream temperature.
   * \return Freestream temperature.
   */
  su2double GetModVel_FreeStream(void);
  
  /*!
   * \brief Get the value of the frestream temperature.
   * \return Freestream temperature.
   */
  su2double GetModVel_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the frestream vibrational-electronic temperature.
   * \return Freestream temperature.
   */
  su2double GetTemperature_ve_FreeStream(void);
  
  /*!
   * \brief Get the value of the laminar Prandtl number.
   * \return Laminar Prandtl number.
   */
  su2double GetPrandtl_Lam(void);
  
  /*!
   * \brief Get the value of the turbulent Prandtl number.
   * \return Turbulent Prandtl number.
   */
  su2double GetPrandtl_Turb(void);

  /*!
   * \brief Get the value of the thermal conductivity for solids.
   * \return Thermal conductivity (solid).
   */
  su2double GetThermalConductivity_Solid(void);

  /*!
   * \brief Get the value of the thermal diffusivity for solids.
   * \return Thermal conductivity (solid).
   */
  su2double GetThermalDiffusivity_Solid(void);

  /*!
   * \brief Get the temperature in solids at freestream conditions.
   * \return Freestream temperature (solid).
   */
  su2double GetTemperature_Freestream_Solid(void);
  
  /*!
   * \brief Get the value of the reference length for non-dimensionalization.
   *        This value should always be 1 internally, and is not user-specified.
   * \return Reference length for non-dimensionalization.
   */
  su2double GetLength_Ref(void);
  
  /*!
   * \brief Get the value of the reference pressure for non-dimensionalization.
   * \return Reference pressure for non-dimensionalization.
   */
  su2double GetPressure_Ref(void);
  
  /*!
   * \brief Get the value of the reference pressure for non-dimensionalization.
   * \return Reference pressure for non-dimensionalization.
   */
  su2double GetEnergy_Ref(void);
  
  /*!
   * \brief Get the value of the reference temperature for non-dimensionalization.
   * \return Reference temperature for non-dimensionalization.
   */
  su2double GetTemperature_Ref(void);
  
  /*!
   * \brief Get the value of the reference density for non-dimensionalization.
   * \return Reference density for non-dimensionalization.
   */
  su2double GetDensity_Ref(void);
  
  /*!
   * \brief Get the value of the reference velocity for non-dimensionalization.
   * \return Reference velocity for non-dimensionalization.
   */
  su2double GetVelocity_Ref(void);
  
  /*!
   * \brief Get the value of the reference time for non-dimensionalization.
   * \return Reference time for non-dimensionalization.
   */
  su2double GetTime_Ref(void);
  
  /*!
   * \brief Get the value of the reference viscosity for non-dimensionalization.
   * \return Reference viscosity for non-dimensionalization.
   */
  su2double GetViscosity_Ref(void);
  
  /*!
   * \brief Get the value of the reference viscosity for non-dimensionalization.
   * \return Reference viscosity for non-dimensionalization.
   */
  su2double GetHighlite_Area(void);
  
  /*!
   * \brief Get the value of the reference viscosity for non-dimensionalization.
   * \return Reference viscosity for non-dimensionalization.
   */
  su2double GetFan_Poly_Eff(void);
  
  /*!
   * \brief Get the value of the reference conductivity for non-dimensionalization.
   * \return Reference conductivity for non-dimensionalization.
   */
  su2double GetConductivity_Ref(void);
  
  /*!
   * \brief Get the value of the reference angular velocity for non-dimensionalization.
   * \return Reference angular velocity for non-dimensionalization.
   */
  su2double GetOmega_Ref(void);
  
  /*!
   * \brief Get the value of the reference force for non-dimensionalization.
   * \return Reference force for non-dimensionalization.
   */
  su2double GetForce_Ref(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream pressure.
   * \return Non-dimensionalized freestream pressure.
   */
  su2double GetPressure_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream pressure.
   * \return Non-dimensionalized freestream pressure.
   */
  su2double GetPressure_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the thermodynamic pressure.
   * \return Thermodynamic pressure.
   */
  su2double GetPressure_Thermodynamic(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized thermodynamic pressure.
   * \return Non-dimensionalized thermodynamic pressure.
   */
  su2double GetPressure_ThermodynamicND(void);

  /*!
   * \brief Get the vector of the dimensionalized freestream velocity.
   * \return Dimensionalized freestream velocity vector.
   */
  su2double* GetVelocity_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream temperature.
   * \return Non-dimensionalized freestream temperature.
   */
  su2double GetTemperature_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream density.
   * \return Non-dimensionalized freestream density.
   */
  su2double GetDensity_FreeStreamND(void);
  
  /*!
   * \brief Get the vector of the non-dimensionalized freestream velocity.
   * \return Non-dimensionalized freestream velocity vector.
   */
  su2double* GetVelocity_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream energy.
   * \return Non-dimensionalized freestream energy.
   */
  su2double GetEnergy_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetViscosity_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetTke_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetOmega_FreeStreamND(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetTke_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream viscosity.
   * \return Non-dimensionalized freestream viscosity.
   */
  su2double GetOmega_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream intermittency.
   * \return Non-dimensionalized freestream intermittency.
   */
  su2double GetIntermittency_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream turbulence intensity.
   * \return Non-dimensionalized freestream intensity.
   */
  su2double GetTurbulenceIntensity_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized freestream turbulence intensity.
   * \return Non-dimensionalized freestream intensity.
   */
  su2double GetNuFactor_FreeStream(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized engine turbulence intensity.
   * \return Non-dimensionalized engine intensity.
   */
  su2double GetNuFactor_Engine(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized actuator disk turbulence intensity.
   * \return Non-dimensionalized actuator disk intensity.
   */
  su2double GetSecondaryFlow_ActDisk(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized actuator disk turbulence intensity.
   * \return Non-dimensionalized actuator disk intensity.
   */
  su2double GetInitial_BCThrust(void);
  
  /*!
   * \brief Get the value of the non-dimensionalized actuator disk turbulence intensity.
   * \return Non-dimensionalized actuator disk intensity.
   */
  void SetInitial_BCThrust(su2double val_bcthrust);
  
  /*!
   * \brief Get the value of the turbulent to laminar viscosity ratio.
   * \return Ratio of turbulent to laminar viscosity ratio.
   */
  su2double GetTurb2LamViscRatio_FreeStream(void);
  
  /*!
   * \brief Get the vector of free stream mass fraction values.
   * \return Ratio of species mass to mixture mass.
   */
  su2double* GetMassFrac_FreeStream(void);
  
  /*!
   * \brief Get the value of the Reynolds length.
   * \return Reynolds length.
   */
  su2double GetLength_Reynolds(void);
  
  /*!
   * \brief Get the start up iterations using the fine grid, this works only for multigrid problems.
   * \return Start up iterations using the fine grid.
   */
  unsigned short GetnStartUpIter(void);
  
  /*!
   * \brief Get the reference area for non dimensional coefficient computation. If the value from the
   *        is 0 then, the code will compute the reference area using the projection of the shape into
   *        the z plane (3D) or the x plane (2D).
   * \return Value of the reference area for coefficient computation.
   */
  su2double GetRefArea(void);
  
  /*!
   * \brief Get the wave speed.
   * \return Value of the wave speed.
   */
  su2double GetThermalDiffusivity(void);
  
  /*!
   * \brief Get the thermal expansion coefficient.
   * \return Value of the thermal expansion coefficient.
   */
  su2double GetThermal_Expansion_Coeff(void);

  /*!
   * \brief Get the non-dim. thermal expansion coefficient.
   * \return Value of the non-dim. thermal expansion coefficient.
   */
  su2double GetThermal_Expansion_CoeffND(void);

  /*!
   * \brief Set the thermal expansion coefficient.
   * \param[in] val_thermal_expansion - thermal expansion coefficient
   */
  void SetThermal_Expansion_Coeff(su2double val_thermal_expansion);

  /*!
   * \brief Set the non-dim. thermal expansion coefficient.
   * \param[in] val_thermal_expansion - non-dim. thermal expansion coefficient
   */
  void SetThermal_Expansion_CoeffND(su2double val_thermal_expansionnd);

  /*!
   * \brief Get the value of the reference density for custom incompressible non-dimensionalization.
   * \return Reference density for custom incompressible non-dimensionalization.
   */
  su2double GetInc_Density_Ref(void);

  /*!
   * \brief Get the value of the reference velocity for custom incompressible non-dimensionalization.
   * \return Reference velocity for custom incompressible non-dimensionalization.
   */
  su2double GetInc_Velocity_Ref(void);

  /*!
   * \brief Get the value of the reference temperature for custom incompressible non-dimensionalization.
   * \return Reference temperature for custom incompressible non-dimensionalization.
   */
  su2double GetInc_Temperature_Ref(void);

  /*!
   * \brief Get the value of the initial density for incompressible flows.
   * \return Initial density for incompressible flows.
   */
  su2double GetInc_Density_Init(void);

  /*!
   * \brief Get the value of the initial velocity for incompressible flows.
   * \return Initial velocity for incompressible flows.
   */
  su2double* GetInc_Velocity_Init(void);

  /*!
   * \brief Get the value of the initial temperature for incompressible flows.
   * \return Initial temperature for incompressible flows.
   */
  su2double GetInc_Temperature_Init(void);

  /*!
   * \brief Get the Young's modulus of elasticity.
   * \return Value of the Young's modulus of elasticity.
   */
  su2double GetElasticyMod(unsigned short id_val);
  
  /*!
    * \brief Decide whether to apply DE effects to the model.
    * \return <code>TRUE</code> if the DE effects are to be applied, <code>FALSE</code> otherwise.
    */
  
  bool GetDE_Effects(void);
  
  /*!
    * \brief Decide whether to predict the DE effects for the next time step.
    * \return <code>TRUE</code> if the DE effects are to be applied, <code>FALSE</code> otherwise.
    */
  
  bool GetDE_Predicted(void);
  
  /*!
   * \brief Get the number of different electric constants.
   * \return Value of the DE modulus.
   */
  unsigned short GetnElectric_Constant(void);

  /*!
   * \brief Get the value of the DE modulus.
   * \return Value of the DE modulus.
   */
  su2double GetElectric_Constant(unsigned short iVar);

  /*!
   * \brief Get the value of the B constant in the Knowles material model.
   * \return Value of the B constant in the Knowles material model.
   */
  su2double GetKnowles_B(void);

  /*!
   * \brief Get the value of the N constant in the Knowles material model.
   * \return Value of the N constant in the Knowles material model.
   */
  su2double GetKnowles_N(void);

  /*!
   * \brief Get the kind of design variable for FEA.
   * \return Value of the DE voltage.
   */
  unsigned short GetDV_FEA(void);

  /*!
   * \brief Get the ID of the reference node.
   * \return Number of FSI subiters.
   */
  unsigned long GetRefNode_ID(void);

  /*!
   * \brief Get the values for the reference node displacement.
   * \param[in] val_coeff - Index of the displacement.
   */
  su2double GetRefNode_Displacement(unsigned short val_coeff);

  /*!
   * \brief Get the penalty weight value for the objective function.
   * \return  Penalty weight value for the reference geometry objective function.
   */
  su2double GetRefNode_Penalty(void);

  /*!
    * \brief Decide whether it's necessary to read a reference geometry.
    * \return <code>TRUE</code> if it's necessary to read a reference geometry, <code>FALSE</code> otherwise.
    */

  bool GetRefGeom(void);

  /*!
   * \brief Get the name of the file with the reference geometry of the structural problem.
   * \return Name of the file with the reference geometry of the structural problem.
   */
  string GetRefGeom_FEMFileName(void);

  /*!
   * \brief Get the format of the reference geometry file.
   * \return Format of the reference geometry file.
   */
  unsigned short GetRefGeom_FileFormat(void);

    /*!
   * \brief Formulation for 2D elasticity (plane stress - strain)
   * \return Flag to 2D elasticity model.
   */
  unsigned short GetElas2D_Formulation(void);
  
  /*!
   * \brief Decide whether it's necessary to read a reference geometry.
   * \return <code>TRUE</code> if it's necessary to read a reference geometry, <code>FALSE</code> otherwise.
   */
  
  bool GetPrestretch(void);
  
  /*!
    * \brief Decide whether it's necessary to add the cross term for adjoint FSI.
    * \return <code>TRUE</code> if it's necessary to add the cross term, <code>FALSE</code> otherwise.
    */
  
  bool Add_CrossTerm(void);
  
  /*!
    * \brief Set the boolean addCrossTerm to true or false.
    */
  
  void Set_CrossTerm(bool needCrossTerm);

  /*!
   * \brief Get the name of the file with the element properties for structural problems.
   * \return Name of the file with the element properties of the structural problem.
   */
  string GetFEA_FileName(void);

  /*!
   * \brief Get the name of the file with the reference geometry of the structural problem.
   * \return Name of the file with the reference geometry of the structural problem.
   */
  string GetPrestretch_FEMFileName(void);
  
  /*!
   * \brief Get the Poisson's ratio.
   * \return Value of the Poisson's ratio.
   */
  su2double GetPoissonRatio(unsigned short id_val);
  
  /*!
   * \brief Get the Material Density.
   * \return Value of the Material Density.
   */
  su2double GetMaterialDensity(unsigned short id_val);
  
  /*!
   * \brief Compressibility/incompressibility of the solids analysed using the structural solver.
   * \return Compressible or incompressible.
   */
  unsigned short GetMaterialCompressibility(void);
  
  /*!
   * \brief Compressibility/incompressibility of the solids analysed using the structural solver.
   * \return Compressible or incompressible.
   */
  unsigned short GetMaterialModel(void);
  
  /*!
   * \brief Geometric conditions for the structural solver.
   * \return Small or large deformation structural analysis.
   */
  unsigned short GetGeometricConditions(void);
  
  /*!
   * \brief Get the reference length for computing moment (the default value is 1).
   * \return Reference length for moment computation.
   */
  su2double GetRefLength(void);
  
  /*!
   * \brief Get the reference element length for computing the slope limiting epsilon.
   * \return Reference element length for slope limiting epsilon.
   */
  su2double GetRefElemLength(void);
  
  /*!
   * \brief Get the reference coefficient for detecting sharp edges.
   * \return Reference coefficient for detecting sharp edges.
   */
  su2double GetRefSharpEdges(void);
  
  /*!
   * \brief Get the volume of the whole domain using the fine grid, this value is common for all the grids
   *        in the multigrid method.
   * \return Volume of the whole domain.
   */
  su2double GetDomainVolume(void);
  
  /*!
   * \brief In case the <i>RefArea</i> is equal to 0 then, it is necessary to compute a reference area,
   *        with this function we set the value of the reference area.
   * \param[in] val_area - Value of the reference area for non dimensional coefficient computation.
   */
  void SetRefArea(su2double val_area);
  
  /*!
   * \brief In case the <i>SemiSpan</i> is equal to 0 then, it is necessary to compute the max y distance,
   *        with this function we set the value of the semi span.
   * \param[in] val_semispan - Value of the semispan.
   */
  void SetSemiSpan(su2double val_semispan);
  
  /*!
   * \brief Set the value of the domain volume computed on the finest grid.
   * \note This volume do not include the volume of the body that is being simulated.
   * \param[in] val_volume - Value of the domain volume computed on the finest grid.
   */
  void SetDomainVolume(su2double val_volume);
  
  /*!
   * \brief Set the finest mesh in a multigrid strategy.
   * \note If we are using a Full Multigrid Strategy or a start up with finest grid, it is necessary
   *       to change several times the finest grid.
   * \param[in] val_finestmesh - Index of the finest grid.
   */
  void SetFinestMesh(unsigned short val_finestmesh);
  
  /*!
   * \brief Set the kind of time integration scheme.
   * \note If we are solving different equations it will be necessary to change several
   *       times the kind of time integration, to choose the right scheme.
   * \param[in] val_kind_timeintscheme - Kind of time integration scheme.
   */
  void SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme);
  
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
  void SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme, unsigned short val_kind_centered,
                             unsigned short val_kind_upwind,        unsigned short val_kind_slopelimit,
                             bool val_muscl,                        unsigned short val_kind_fem);

  /*!
   * \brief Get the value of limiter coefficient.
   * \return Value of the limiter coefficient.
   */
  su2double GetVenkat_LimiterCoeff(void);
  
  /*!
   * \brief Freeze the value of the limiter after a number of iterations.
   * \return Number of iterations.
   */
  unsigned long GetLimiterIter(void);
  
  /*!
   * \brief Get the value of sharp edge limiter.
   * \return Value of the sharp edge limiter coefficient.
   */
  su2double GetAdjSharp_LimiterCoeff(void);
  
  /*!
   * \brief Get the Reynolds number. Dimensionless number that gives a measure of the ratio of inertial forces
   *        to viscous forces and consequently quantifies the relative importance of these two types of forces
   *        for given flow condition.
   * \return Value of the Reynolds number.
   */
  su2double GetReynolds(void);
  
  /*!
   * \brief Get the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  su2double GetFroude(void);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetFroude(su2double val_froude);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetMach(su2double val_mach);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetReynolds(su2double val_reynolds);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetLength_Ref(su2double val_length_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetVelocity_Ref(su2double val_velocity_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetPressure_Ref(su2double val_pressure_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetDensity_Ref(su2double val_density_ref);
  
  /*!
   * \brief Set the reference temperature.
   * \return Value of the Froude number.
   */
  void SetTemperature_Ref(su2double val_temperature_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetTime_Ref(su2double val_time_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetEnergy_Ref(su2double val_energy_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetOmega_Ref(su2double val_omega_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetForce_Ref(su2double val_force_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetGas_Constant_Ref(su2double val_gas_constant_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetGas_Constant(su2double val_gas_constant);
  
  /*!
   * \brief Set the value of the specific heat at constant pressure (incompressible fluids with energy equation).
   * \param[in] val_specific_heat_cp - specific heat at constant pressure.
   */
  void SetSpecific_Heat_Cp(su2double val_specific_heat_cp);

  /*!
   * \brief Set the non-dimensional value of the specific heat at constant pressure (incompressible fluids with energy equation).
   * \param[in] val_specific_heat_cpnd - non-dim. specific heat at constant pressure.
   */
  void SetSpecific_Heat_CpND(su2double val_specific_heat_cpnd);

  /*!
   * \brief Set the value of the specific heat at constant volume (incompressible fluids with energy equation).
   * \param[in] val_specific_heat_cv - specific heat at constant volume.
   */
  void SetSpecific_Heat_Cv(su2double val_specific_heat_cv);

  /*!
   * \brief Set the non-dimensional value of the specific heat at constant volume (incompressible fluids with energy equation).
   * \param[in] val_specific_heat_cvnd - non-dim. specific heat at constant pressure.
   */
  void SetSpecific_Heat_CvND(su2double val_specific_heat_cvnd);

  /*!
   * \brief Set the heat flux reference value.
   * \return Value of the reference heat flux.
   */
  void SetHeat_Flux_Ref(su2double val_heat_flux_ref);

  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetViscosity_Ref(su2double val_viscosity_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetConductivity_Ref(su2double val_conductivity_ref);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetPressure_FreeStreamND(su2double val_pressure_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetPressure_FreeStream(su2double val_pressure_freestream);
  
  /*!
   * \brief Set the non-dimensionalized thermodynamic pressure for low Mach problems.
   * \return Value of the non-dimensionalized thermodynamic pressure.
   */
  void SetPressure_ThermodynamicND(su2double val_pressure_thermodynamicnd);
  
  /*!
   * \brief Set the thermodynamic pressure for low Mach problems.
   * \return Value of the thermodynamic pressure.
   */
  void SetPressure_Thermodynamic(su2double val_pressure_thermodynamic);

  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetDensity_FreeStreamND(su2double val_density_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetDensity_FreeStream(su2double val_density_freestream);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetViscosity_FreeStream(su2double val_viscosity_freestream);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetModVel_FreeStream(su2double val_modvel_freestream);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetModVel_FreeStreamND(su2double val_modvel_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetTemperature_FreeStream(su2double val_temperature_freestream);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetTemperature_FreeStreamND(su2double val_temperature_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetGas_ConstantND(su2double val_gas_constantnd);
  
  /*!
   * \brief Set the free-stream velocity.
   * \param[in] val_velocity_freestream - Value of the free-stream velocity component.
   * \param[in] val_dim - Value of the current dimension.
   */
  void SetVelocity_FreeStream(su2double val_velocity_freestream, unsigned short val_dim);

  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetVelocity_FreeStreamND(su2double val_velocity_freestreamnd, unsigned short val_dim);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetViscosity_FreeStreamND(su2double val_viscosity_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetTke_FreeStreamND(su2double val_tke_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetOmega_FreeStreamND(su2double val_omega_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetTke_FreeStream(su2double val_tke_freestream);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetOmega_FreeStream(su2double val_omega_freestream);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetEnergy_FreeStreamND(su2double val_energy_freestreamnd);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetEnergy_FreeStream(su2double val_energy_freestream);

  /*!
   * \brief Set the thermal diffusivity for solids.
   * \return Value of the Froude number.
   */
  void SetThermalDiffusivity_Solid(su2double val_thermal_diffusivity);
  
  /*!
   * \brief Set the Froude number for free surface problems.
   * \return Value of the Froude number.
   */
  void SetTotal_UnstTimeND(su2double val_total_unsttimend);
  
  /*!
   * \brief Get the angle of attack of the body. This is the angle between a reference line on a lifting body
   *        (often the chord line of an airfoil) and the vector representing the relative motion between the
   *        lifting body and the fluid through which it is moving.
   * \return Value of the angle of attack.
   */
  su2double GetAoA(void);
  
  /*!
   * \brief Get the off set angle of attack of the body. The solution and the geometry
   *        file are able to modifity the angle of attack in the config file
   * \return Value of the off set angle of attack.
   */
  su2double GetAoA_Offset(void);
  
  /*!
   * \brief Get the off set sideslip angle of the body. The solution and the geometry
   *        file are able to modifity the angle of attack in the config file
   * \return Value of the off set sideslip angle.
   */
  su2double GetAoS_Offset(void);
  
  /*!
   * \brief Get the functional sensitivity with respect to changes in the angle of attack.
   * \return Value of the angle of attack.
   */
  su2double GetAoA_Sens(void);
  
  /*!
   * \brief Set the angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoA(su2double val_AoA);
  
  /*!
   * \brief Set the off set angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoA_Offset(su2double val_AoA_offset);
  
  /*!
   * \brief Set the off set sideslip angle.
   * \param[in] val_AoA - Value of the off set sideslip angle.
   */
  void SetAoS_Offset(su2double val_AoS_offset);
  
  /*!
   * \brief Set the angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoA_Sens(su2double val_AoA_sens);
  
  /*!
   * \brief Set the angle of attack.
   * \param[in] val_AoA - Value of the angle of attack.
   */
  void SetAoS(su2double val_AoS);
  
  /*!
   * \brief Get the angle of sideslip of the body. It relates to the rotation of the aircraft centerline from
   *        the relative wind.
   * \return Value of the angle of sideslip.
   */
  su2double GetAoS(void);
  
  /*!
   * \brief Get the charge coefficient that is used in the poissonal potential simulation.
   * \return Value of the charge coefficient.
   */
  su2double GetChargeCoeff(void);
  
  /*!
   * \brief Get the number of multigrid levels.
   * \return Number of multigrid levels (without including the original grid).
   */
  unsigned short GetnMGLevels(void);
  
  /*!
   * \brief Set the number of multigrid levels.
   * \param[in] val_nMGLevels - Index of the mesh were the CFL is applied
   */
  void SetMGLevels(unsigned short val_nMGLevels);
  
  /*!
   * \brief Get the index of the finest grid.
   * \return Index of the finest grid in a multigrid strategy, this is 0 unless we are
   performing a Full multigrid.
   */
  unsigned short GetFinestMesh(void);
  
  /*!
   * \brief Get the kind of multigrid (V or W).
   * \note This variable is used in a recursive way to perform the different kind of cycles
   * \return 0 or 1 depending of we are dealing with a V or W cycle.
   */
  unsigned short GetMGCycle(void);
  
  /*!
   * \brief Get the king of evaluation in the geometrical module.
   * \return 0 or 1 depending of we are dealing with a V or W cycle.
   */
  unsigned short GetGeometryMode(void);
  
  /*!
   * \brief Get the Courant Friedrich Levi number for each grid.
   * \param[in] val_mesh - Index of the mesh were the CFL is applied.
   * \return CFL number for each grid.
   */
  su2double GetCFL(unsigned short val_mesh);

  /*!
   * \brief Get the Courant Friedrich Levi number for solid solvers.
   * \param[in] val_mesh - Index of the mesh were the CFL is applied.
   * \return CFL number for each grid.
   */
  su2double GetCFL_Solid(void);
  
  /*!
   * \brief Get the Courant Friedrich Levi number for each grid.
   * \param[in] val_mesh - Index of the mesh were the CFL is applied.
   * \return CFL number for each grid.
   */
  void SetCFL(unsigned short val_mesh, su2double val_cfl);

  /*!
   * \brief Get the Courant Friedrich Levi number for unsteady simulations.
   * \return CFL number for unsteady simulations.
   */
  su2double GetUnst_CFL(void);

  /*!
   * \brief Get information about element reorientation
   * \return 	<code>TRUE</code> means that elements can be reoriented if suspected unhealthy
   */
  bool GetReorientElements(void);
  
  /*!
   * \brief Get the Courant Friedrich Levi number for unsteady simulations.
   * \return CFL number for unsteady simulations.
   */
  su2double GetMax_DeltaTime(void);
  
  /*!
   * \brief Get a parameter of the particular design variable.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \param[in] val_param - Index of the parameter that we want to read.
   * \return Design variable parameter.
   */
  su2double GetParamDV(unsigned short val_dv, unsigned short val_param);
  
  /*!
   * \brief Get the coordinates of the FFD corner points.
   * \param[in] val_ffd - Index of the FFD box.
   * \param[in] val_coord - Index of the coordinate that we want to read.
   * \return Value of the coordinate.
   */
  su2double GetCoordFFDBox(unsigned short val_ffd, unsigned short val_index);
  
  /*!
   * \brief Get the degree of the FFD corner points.
   * \param[in] val_ffd - Index of the FFD box.
   * \param[in] val_degree - Index (I,J,K) to obtain the degree.
   * \return Value of the degree in a particular direction.
   */
  unsigned short GetDegreeFFDBox(unsigned short val_ffd, unsigned short val_index);
  
  /*!
   * \brief Get the FFD Tag of a particular design variable.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \return Name of the FFD box.
   */
  string GetFFDTag(unsigned short val_dv);
  
  /*!
   * \brief Get the FFD Tag of a particular FFD box.
   * \param[in] val_ffd - Number of the FFD box that we want to read.
   * \return Name of the FFD box.
   */
  string GetTagFFDBox(unsigned short val_ffd);
  
  /*!
   * \brief Get the number of design variables.
   * \return Number of the design variables.
   */
  unsigned short GetnDV(void);
  
  /*!
   * \brief Get the number of design variables.
   * \return Number of the design variables.
   */
  unsigned short GetnDV_Value(unsigned short iDV);
  
  /*!
   * \brief Get the number of FFD boxes.
   * \return Number of FFD boxes.
   */
  unsigned short GetnFFDBox(void);
  
  /*!
   * \brief Get the required continuity level at the surface intersection with the FFD
   * \return Continuity level at the surface intersection.
   */
  unsigned short GetFFD_Continuity(void);
  
  /*!
   * \brief Get the coordinate system that we are going to use to define the FFD
   * \return Coordinate system (cartesian, spherical, etc).
   */
  unsigned short GetFFD_CoordSystem(void);
  
  /*!
   * \brief Get the kind of FFD Blending function.
   * \return Kind of FFD Blending function.
   */
  unsigned short GetFFD_Blending(void);
  
  /*!
   * \brief Get the kind BSpline Order in i,j,k direction.
   * \return The kind BSpline Order in i,j,k direction.
   */
  su2double* GetFFD_BSplineOrder();
  
  /*!
   * \brief Get the number of Runge-Kutta steps.
   * \return Number of Runge-Kutta steps.
   */
  unsigned short GetnRKStep(void);

  /*!
   * \brief Get the number of time levels for time accurate local time stepping.
   * \return Number of time levels.
   */
  unsigned short GetnLevels_TimeAccurateLTS(void);

  /*!
   * \brief Set the number of time levels for time accurate local time stepping.
   * \param[in] val_nLevels - The number of time levels to be set.
   */
  void SetnLevels_TimeAccurateLTS(unsigned short val_nLevels);

  /*!
   * \brief Get the number time DOFs for ADER-DG.
   * \return Number of time DOFs used in ADER-DG.
   */
  unsigned short GetnTimeDOFsADER_DG(void);

  /*!
   * \brief Get the location of the time DOFs for ADER-DG on the interval [-1..1].
   * \return The location of the time DOFs used in ADER-DG.
   */
  su2double *GetTimeDOFsADER_DG(void);

  /*!
   * \brief Get the number time integration points for ADER-DG.
   * \return Number of time integration points used in ADER-DG.
   */
  unsigned short GetnTimeIntegrationADER_DG(void);

  /*!
   * \brief Get the location of the time integration points for ADER-DG on the interval [-1..1].
   * \return The location of the time integration points used in ADER-DG.
   */
  su2double *GetTimeIntegrationADER_DG(void);

  /*!
   * \brief Get the weights of the time integration points for ADER-DG.
   * \return The weights of the time integration points used in ADER-DG.
   */
  su2double *GetWeightsIntegrationADER_DG(void);

  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_All(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_Max(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_EngineInflow(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_EngineExhaust(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_NearFieldBound(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_InterfaceBound(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_Fluid_InterfaceBound(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_ActDiskInlet(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_ActDiskOutlet(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_Outlet(void);
  
  /*!
   * \brief Get the total number of monitoring markers.
   * \return Total number of monitoring markers.
   */
  unsigned short GetnMarker_Monitoring(void);
  
  /*!
   * \brief Get the total number of DV markers.
   * \return Total number of DV markers.
   */
  unsigned short GetnMarker_DV(void);
  
  /*!
   * \brief Get the total number of moving markers.
   * \return Total number of moving markers.
   */
  unsigned short GetnMarker_Moving(void);

  /*!
   * \brief Get the total number of Python customizable markers.
   * \return Total number of Python customizable markers.
   */
  unsigned short GetnMarker_PyCustom(void);
  
  /*!
   * \brief Get the total number of moving markers.
   * \return Total number of moving markers.
   */
  unsigned short GetnMarker_Analyze(void);

  /*!
   * \brief Get the total number of periodic markers.
   * \return Total number of periodic markers.
   */
  unsigned short GetnMarker_Periodic(void);

  /*!
   * \brief Get the total number of heat flux markers.
   * \return Total number of heat flux markers.
   */
  unsigned short GetnMarker_HeatFlux(void);
  
  /*!
   * \brief Get the total number of objectives in kind_objective list
   * \return Total number of objectives in kind_objective list
   */
  unsigned short GetnObj(void);
  
  /*!
   * \brief Stores the number of marker in the simulation.
   * \param[in] val_nmarker - Number of markers of the problem.
   */
  void SetnMarker_All(unsigned short val_nmarker);
  
  /*!
   * \brief Get the number of external iterations.
   * \return Number of external iterations.
   */
  unsigned long GetnExtIter(void);
  
  /*!
   * \brief Get the number of internal iterations.
   * \return Number of internal iterations.
   */
  unsigned long GetUnst_nIntIter(void);
  
  /*!
   * \brief Get the number of internal iterations for the Newton-Raphson Method in nonlinear structural applications.
   * \return Number of internal iterations.
   */
  unsigned long GetDyn_nIntIter(void);
  
  /*!
   * \brief Get the restart iteration number for unsteady simulations.
   * \return Restart iteration number for unsteady simulations.
   */
  long GetUnst_RestartIter(void);
  
  /*!
   * \brief Get the starting direct iteration number for the unsteady adjoint (reverse time integration).
   * \return Starting direct iteration number for the unsteady adjoint.
   */
  long GetUnst_AdjointIter(void);
  
  /*!
   * \brief Number of iterations to average (reverse time integration).
   * \return Starting direct iteration number for the unsteady adjoint.
   */
  unsigned long GetIter_Avg_Objective(void);
  
  /*!
   * \brief Get the restart iteration number for dynamic structural simulations.
   * \return Restart iteration number for dynamic structural simulations.
   */
  long GetDyn_RestartIter(void);

  /*!
   * \brief Retrieves the number of periodic time instances for Harmonic Balance.
   * \return: Number of periodic time instances for Harmonic Balance.
   */
  unsigned short GetnTimeInstances(void);
  
  /*!
   * \brief Retrieves the period of oscillations to be used with Harmonic Balance.
   * \return: Period for Harmonic Balance.
   */
  su2double GetHarmonicBalance_Period(void);
  
  /*!
   * \brief Set the number of external iterations.
   * \note This is important in no time depending methods, where only
   *       one external iteration is needed.
   * \param[in] val_niter - Set the number of external iterations.
   */
  void SetnExtIter(unsigned long val_niter);
  
  /*!
   * \brief Set the current external iteration number.
   * \param[in] val_iter - Current external iteration number.
   */
  void SetExtIter(unsigned long val_iter);
  
  /*!
   * \brief Set the current external iteration number.
   * \param[in] val_iter - Current external iteration number.
   */
  void SetExtIter_OffSet(unsigned long val_iter);
  
  /*!
   * \brief Set the current FSI iteration number.
   * \param[in] val_iter - Current FSI iteration number.
   */
  void SetOuterIter(unsigned long val_iter);
  
  /*!
   * \brief Set the current internal iteration number.
   * \param[in] val_iter - Current external iteration number.
   */
  void SetIntIter(unsigned long val_iter);
  
  /*!
   * \brief Get the current external iteration number.
   * \return Current external iteration.
   */
  unsigned long GetExtIter(void);
  
  /*!
   * \brief Get the current internal iteration number.
   * \return Current external iteration.
   */
  unsigned long GetExtIter_OffSet(void);
  
  /*!
   * \brief Get the current FSI iteration number.
   * \return Current FSI iteration.
   */
  unsigned long GetOuterIter(void);
  
  /*!
   * \brief Get the current internal iteration number.
   * \return Current internal iteration.
   */
  unsigned long GetIntIter(void);

  /*!
   * \brief Set the current physical time.
   * \param[in] val_t - Current physical time.
   */
  void SetPhysicalTime(su2double val_t);
  
  /*!
   * \brief Get the current physical time.
   * \return Current physical time.
   */
  su2double GetPhysicalTime(void);
  
  /*!
   * \brief Get the frequency for writing the solution file.
   * \return It writes the solution file with this frequency.
   */
  unsigned long GetWrt_Sol_Freq(void);
  
  /*!
   * \brief Get the frequency for writing the solution file in Dual Time.
   * \return It writes the solution file with this frequency.
   */
  unsigned long GetWrt_Sol_Freq_DualTime(void);
  
  /*!
   * \brief Get the frequency for writing the convergence file.
   * \return It writes the convergence file with this frequency.
   */
  unsigned long GetWrt_Con_Freq(void);
  
  /*!
   * \brief Set the frequency for writing the convergence file.
   * \return It writes the convergence file with this frequency.
   */
  void SetWrt_Con_Freq(unsigned long val_freq);

  /*!
   * \brief Get the frequency for writing the convergence file in Dual Time.
   * \return It writes the convergence file with this frequency.
   */
  unsigned long GetWrt_Con_Freq_DualTime(void);
  
  /*!
   * \brief Get information about writing unsteady headers and file extensions.
   * \return 	<code>TRUE</code> means that unsteady solution files will be written.
   */
  bool GetWrt_Unsteady(void);

  /*!
   * \brief Get information about writing output files.
   * \return <code>TRUE</code> means that output files will be written.
   */
  bool GetWrt_Output(void);

  /*!
   * \brief Get information about writing a volume solution file.
   * \return <code>TRUE</code> means that a volume solution file will be written.
   */
  bool GetWrt_Vol_Sol(void);
  
  /*!
   * \brief Get information about writing a volume solution file.
   * \return <code>TRUE</code> means that a volume solution file will be written.
   */
  bool GetLow_MemoryOutput(void);
  
  /*!
   * \brief Get information about writing a surface solution file.
   * \return <code>TRUE</code> means that a surface solution file will be written.
   */
  bool GetWrt_Srf_Sol(void);
  
  /*!
   * \brief Get information about writing a surface comma-separated values (CSV) solution file.
   * \return <code>TRUE</code> means that a surface comma-separated values (CSV) solution file will be written.
   */
  bool GetWrt_Csv_Sol(void);
  
  /*!
   * \brief Get information about writing a binary coordinates file.
   * \return <code>TRUE</code> means that a binary coordinates file will be written.
   */
  bool GetWrt_Crd_Sol(void);
  
  /*!
   * \brief Get information about writing residuals to volume solution file.
   * \return <code>TRUE</code> means that residuals will be written to the solution file.
   */
  bool GetWrt_Residuals(void);
  
  /*!
   * \brief Get information about writing residuals to volume solution file.
   * \return <code>TRUE</code> means that residuals will be written to the solution file.
   */
  bool GetWrt_Limiters(void);
  
  /*!
   * \brief Write solution at each surface.
   * \return <code>TRUE</code> means that the solution at each surface will be written.
   */
  bool GetWrt_Surface(void);
  
  /*!
   * \brief Get information about writing residuals to volume solution file.
   * \return <code>TRUE</code> means that residuals will be written to the solution file.
   */
  bool GetWrt_SharpEdges(void);
  
  /*!
   * \brief Get information about writing rind layers to the solution files.
   * \return <code>TRUE</code> means that rind layers will be written to the solution file.
   */
  bool GetWrt_Halo(void);

  /*!
   * \brief Get information about writing the performance summary at the end of a calculation.
   * \return <code>TRUE</code> means that the performance summary will be written at the end of a calculation.
   */
  bool GetWrt_Performance(void);
  
  /*!
   * \brief Get information about writing a template inlet profile file.
   * \return <code>TRUE</code> means that a template inlet profile file will be written.
   */
  bool GetWrt_InletFile(void);

  /*!
   * \brief Set information about writing a template inlet profile file.
   * \param[in] val_wrt_inletfile - flag for whether to write a template inlet profile file.
   */
  void SetWrt_InletFile(bool val_wrt_inletfile);

  /*!
   * \brief Get information about writing a 1D slice of a 2D cartesian solution.
   * \return <code>TRUE</code> means that a 1D slice of a 2D cartesian solution will be written.
   */
  bool GetWrt_Slice(void);

  /*!
   * \brief Get information about writing projected sensitivities on surfaces to an ASCII file with rows as x, y, z, dJ/dx, dJ/dy, dJ/dz for each vertex.
   * \return <code>TRUE</code> means that projected sensitivities on surfaces in an ASCII file with rows as x, y, z, dJ/dx, dJ/dy, dJ/dz for each vertex will be written.
   */
  bool GetWrt_Projected_Sensitivity(void);
  
  /*!
   * \brief Get information about the format for the input volume sensitvities.
   * \return Format of the input volume sensitivities.
   */
  unsigned short GetSensitivity_Format(void);
  
  /*!
   * \brief Get information about writing sectional force files.
   * \return <code>TRUE</code> means that sectional force files will be written for specified markers.
   */
  bool GetPlot_Section_Forces(void);
  
  /*!
   * \brief Get the alpha (convective) coefficients for the Runge-Kutta integration scheme.
   * \param[in] val_step - Index of the step.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double Get_Alpha_RKStep(unsigned short val_step);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_All_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_ActDiskInlet_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_ActDiskOutlet_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Outlet_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_EngineInflow_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_EngineExhaust_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Monitoring_TagBound(unsigned short val_marker);

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_HeatFlux_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the tag if the iMarker defined in the geometry file.
   * \param[in] val_tag - Value of the tag in which we are interested.
   * \return Value of the marker <i>val_marker</i> that is in the geometry file
   *         for the surface that has the tag.
   */
  short GetMarker_All_TagBound(string val_tag);
  
  /*!
   * \brief Get the kind of boundary for each marker.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \return Kind of boundary for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_All_KindBC(unsigned short val_marker);
  
  /*!
   * \brief Set the value of the boundary <i>val_boundary</i> (read from the config file)
   *        for the marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_boundary - Kind of boundary read from config file.
   */
  void SetMarker_All_KindBC(unsigned short val_marker, unsigned short val_boundary);
  
  /*!
   * \brief Set the value of the index <i>val_index</i> (read from the geometry file) for
   *        the marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_index - Index of the surface read from geometry file.
   */
  void SetMarker_All_TagBound(unsigned short val_marker, string val_index);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be monitored <i>val_monitoring</i>
   *        (read from the config file).
   * \note This is important for non dimensional coefficient computation.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be monitored.
   */
  void SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be monitored <i>val_monitoring</i>
   *        (read from the config file).
   * \note This is important for non dimensional coefficient computation.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be monitored.
   */
  void SetMarker_All_GeoEval(unsigned short val_marker, unsigned short val_geoeval);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be designed <i>val_designing</i>
   *        (read from the config file).
   * \note This is important for non dimensional coefficient computation.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be designed.
   */
  void SetMarker_All_Designing(unsigned short val_marker, unsigned short val_designing);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be plot <i>val_plotting</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_plotting - 0 or 1 depending if the the marker is going to be plot.
   */
  void SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be plot <i>val_plotting</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_plotting - 0 or 1 depending if the the marker is going to be plot.
   */
  void SetMarker_All_Analyze(unsigned short val_marker, unsigned short val_analyze);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is part of the FSI interface <i>val_plotting</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_plotting - 0 or 1 depending if the the marker is part of the FSI interface.
   */
  void SetMarker_All_ZoneInterface(unsigned short val_marker, unsigned short val_fsiinterface);
 
  /*!
   * \brief Set if a marker <i>val_marker</i> is part of the Turbomachinery (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_turboperf - 0 if not part of Turbomachinery or greater than 1 if it is part.
   */
  void SetMarker_All_Turbomachinery(unsigned short val_marker, unsigned short val_turbo);

  /*!
   * \brief Set a flag to the marker <i>val_marker</i> part of the Turbomachinery (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_turboperflag - 0 if is not part of the Turbomachinery, flag INFLOW or OUTFLOW if it is part.
   */
  void SetMarker_All_TurbomachineryFlag(unsigned short val_marker, unsigned short val_turboflag);

  /*!
   * \brief Set if a marker <i>val_marker</i> is part of the MixingPlane interface (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_turboperf - 0 if not part of the MixingPlane interface or greater than 1 if it is part.
   */
  void SetMarker_All_MixingPlaneInterface(unsigned short val_marker, unsigned short val_mixplan_interface);
   
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be affected by design variables <i>val_moving</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_DV - 0 or 1 depending if the the marker is affected by design variables.
   */
  void SetMarker_All_DV(unsigned short val_marker, unsigned short val_DV);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be moved <i>val_moving</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_moving - 0 or 1 depending if the the marker is going to be moved.
   */
  void SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving);

  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be customized in Python <i>val_PyCustom</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_PyCustom - 0 or 1 depending if the the marker is going to be customized in Python.
   */
  void SetMarker_All_PyCustom(unsigned short val_marker, unsigned short val_PyCustom);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be periodic <i>val_perbound</i>
   *        (read from the config file).
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_perbound - Index of the surface with the periodic boundary.
   */
  void SetMarker_All_PerBound(unsigned short val_marker, short val_perbound);
  
  /*!
   * \brief Set if a marker <i>val_marker</i> is going to be sent or receive <i>val_index</i>
   *        from another domain.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \param[in] val_index - Index of the surface read from geometry file.
   */
  void SetMarker_All_SendRecv(unsigned short val_marker, short val_index);
  
  /*!
   * \brief Get the send-receive information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return If positive, the information is sended to that domain, in case negative
   *         the information is receive from that domain.
   */
  short GetMarker_All_SendRecv(unsigned short val_marker);
  
  /*!
   * \brief Get an internal index that identify the periodic boundary conditions.
   * \param[in] val_marker - Value of the marker that correspond with the periodic boundary.
   * \return The internal index of the periodic boundary condition.
   */
  short GetMarker_All_PerBound(unsigned short val_marker);
  
  /*!
   * \brief Get the monitoring information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be monitored.
   * \return 0 or 1 depending if the marker is going to be monitored.
   */
  unsigned short GetMarker_All_Monitoring(unsigned short val_marker);
  
  /*!
   * \brief Get the monitoring information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be monitored.
   * \return 0 or 1 depending if the marker is going to be monitored.
   */
  unsigned short GetMarker_All_GeoEval(unsigned short val_marker);
  
  /*!
   * \brief Get the design information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be monitored.
   * \return 0 or 1 depending if the marker is going to be monitored.
   */
  unsigned short GetMarker_All_Designing(unsigned short val_marker);
  
  /*!
   * \brief Get the plotting information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is going to be plotted.
   */
  unsigned short GetMarker_All_Plotting(unsigned short val_marker);
  
  /*!
   * \brief Get the plotting information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is going to be plotted.
   */
  unsigned short GetMarker_All_Analyze(unsigned short val_marker);
  
  /*!
   * \brief Get the FSI interface information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is part of the FSI interface.
   */
  unsigned short GetMarker_All_ZoneInterface(unsigned short val_marker);
  
  /*!
	 * \brief Get the MixingPlane interface information for a marker <i>val_marker</i>.
	 * \param[in] val_marker value of the marker on the grid.
	 * \return 0 if is not part of the MixingPlane Interface and greater than 1 if it is part.
	 */
	unsigned short GetMarker_All_MixingPlaneInterface(unsigned short val_marker);

	/*!
	 * \brief Get the Turbomachinery information for a marker <i>val_marker</i>.
	 * \param[in] val_marker value of the marker on the grid.
	 * \return 0 if is not part of the Turbomachinery and greater than 1 if it is part.
	 */
	unsigned short GetMarker_All_Turbomachinery(unsigned short val_marker);

	/*!
	 * \brief Get the Turbomachinery flag information for a marker <i>val_marker</i>.
	 * \param[in] val_marker value of the marker on the grid.
	 * \return 0 if is not part of the Turbomachinery, flag INFLOW or OUTFLOW if it is part.
	 */
	unsigned short GetMarker_All_TurbomachineryFlag(unsigned short val_marker);

	/*!
   * \brief Get the number of FSI interface markers <i>val_marker</i>.
   * \param[in] void.
   * \return Number of markers belonging to the FSI interface.
   */
  unsigned short GetMarker_n_ZoneInterface(void);
  
  /*!
   * \brief Get the DV information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be affected by design variables.
   * \return 0 or 1 depending if the marker is going to be affected by design variables.
   */
  unsigned short GetMarker_All_DV(unsigned short val_marker);
  
  /*!
   * \brief Get the motion information for a marker <i>val_marker</i>.
   * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
   * \return 0 or 1 depending if the marker is going to be moved.
   */
  unsigned short GetMarker_All_Moving(unsigned short val_marker);

  /*!
   * \brief Get the Python customization for a marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \return 0 or 1 depending if the marker is going to be customized in Python.
   */
  unsigned short GetMarker_All_PyCustom(unsigned short val_marker);
  
  /*!
   * \brief Get the airfoil sections in the slicing process.
   * \param[in] val_section - Index of the section.
   * \return Coordinate of the airfoil to slice.
   */
  su2double GetLocationStations(unsigned short val_section);
  
  /*!
   * \brief Get the defintion of the nacelle location.
   * \param[in] val_index - Index of the section.
   * \return Coordinate of the nacelle location.
   */
  su2double GetNacelleLocation(unsigned short val_index);
  
  /*!
   * \brief Get the number of pre-smoothings in a multigrid strategy.
   * \param[in] val_mesh - Index of the grid.
   * \return Number of smoothing iterations.
   */
  unsigned short GetMG_PreSmooth(unsigned short val_mesh);
  
  /*!
   * \brief Get the number of post-smoothings in a multigrid strategy.
   * \param[in] val_mesh - Index of the grid.
   * \return Number of smoothing iterations.
   */
  unsigned short GetMG_PostSmooth(unsigned short val_mesh);
  
  /*!
   * \brief Get the number of implicit Jacobi smoothings of the correction in a multigrid strategy.
   * \param[in] val_mesh - Index of the grid.
   * \return Number of implicit smoothing iterations.
   */
  unsigned short GetMG_CorrecSmooth(unsigned short val_mesh);
  
  /*!
   * \brief plane of the FFD (I axis) that should be fixed.
   * \param[in] val_index - Index of the arrray with all the planes in the I direction that should be fixed.
   * \return Index of the plane that is going to be freeze.
   */
  short GetFFD_Fix_IDir(unsigned short val_index);
  
  /*!
   * \brief plane of the FFD (J axis) that should be fixed.
   * \param[in] val_index - Index of the arrray with all the planes in the J direction that should be fixed.
   * \return Index of the plane that is going to be freeze.
   */
  short GetFFD_Fix_JDir(unsigned short val_index);
  
  /*!
   * \brief plane of the FFD (K axis) that should be fixed.
   * \param[in] val_index - Index of the arrray with all the planes in the K direction that should be fixed.
   * \return Index of the plane that is going to be freeze.
   */
  short GetFFD_Fix_KDir(unsigned short val_index);
  
  /*!
   * \brief Get the number of planes to fix in the I direction.
   * \return Number of planes to fix in the I direction.
   */
  unsigned short GetnFFD_Fix_IDir(void);
  
  /*!
   * \brief Get the number of planes to fix in the J direction.
   * \return Number of planes to fix in the J direction.
   */
  unsigned short GetnFFD_Fix_JDir(void);
  
  /*!
   * \brief Get the number of planes to fix in the K direction.
   * \return Number of planes to fix in the K direction.
   */
  unsigned short GetnFFD_Fix_KDir(void);
  
  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  unsigned short GetKind_Solver(void);
  
  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  void SetKind_Solver(unsigned short val_solver);
  
  /*!
   * \brief Kind of Multizone Solver.
   * \return Governing equation that we are solving.
   */
  unsigned short GetKind_MZSolver(void);

  
  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  unsigned short GetKind_Regime(void);
  
  /*!
   * \brief Governing equations of the flow (it can be different from the run time equation).
   * \param[in] val_zone - Zone where the soler is applied.
   * \return Governing equation that we are solving.
   */
  unsigned short GetSystemMeasurements(void);
  
  /*!
   * \brief Gas model that we are using.
   * \return Gas model that we are using.
   */
  unsigned short GetKind_GasModel(void);
  
  /*!
   * \brief Fluid model that we are using.
   * \return Fluid model that we are using.
   */
  unsigned short GetKind_FluidModel(void);
  
  /*!
   * \brief Option to define the density model for incompressible flows.
   * \return Density model option
   */
  unsigned short GetKind_DensityModel(void);
  
  /*!
   * \brief Flag for whether to solve the energy equation for incompressible flows.
   * \return Flag for energy equation
   */
  bool GetEnergy_Equation(void);

  /*!
   * \brief free stream option to initialize the solution
   * \return free stream option
   */
  unsigned short GetKind_FreeStreamOption(void);
  
  /*!
   * \brief free stream option to initialize the solution
   * \return free stream option
   */
  unsigned short GetKind_InitOption(void);
  /*!
   * \brief Get the value of the critical pressure.
   * \return Critical pressure.
   */
  su2double GetPressure_Critical(void);
  
  /*!
   * \brief Get the value of the critical temperature.
   * \return Critical temperature.
   */
  su2double GetTemperature_Critical(void);
  
  /*!
   * \brief Get the value of the critical pressure.
   * \return Critical pressure.
   */
  su2double GetAcentric_Factor(void);
  
  /*!
   * \brief Get the value of the viscosity model.
   * \return Viscosity model.
   */
  unsigned short GetKind_ViscosityModel(void);
  
  /*!
   * \brief Get the value of the thermal conductivity model.
   * \return Conductivity model.
   */
  unsigned short GetKind_ConductivityModel(void);
  
  /*!
   * \brief Get the value of the turbulent thermal conductivity model.
   * \return Turbulent conductivity model.
   */
  unsigned short GetKind_ConductivityModel_Turb(void);
  
  /*!
   * \brief Get the value of the constant viscosity.
   * \return Constant viscosity.
   */
  su2double GetMu_Constant(void);

  /*!
   * \brief Get the value of the non-dimensional constant viscosity.
   * \return Non-dimensional constant viscosity.
   */
  su2double GetMu_ConstantND(void);

  /*!
   * \brief Get the value of the thermal conductivity.
   * \return Thermal conductivity.
   */
  su2double GetKt_Constant(void);
  
  /*!
   * \brief Get the value of the non-dimensional thermal conductivity.
   * \return Non-dimensional thermal conductivity.
   */
  su2double GetKt_ConstantND(void);
  
  /*!
   * \brief Get the value of the reference viscosity for Sutherland model.
   * \return The reference viscosity.
   */
  su2double GetMu_Ref(void);

  /*!
   * \brief Get the value of the non-dimensional reference viscosity for Sutherland model.
   * \return The non-dimensional reference viscosity.
   */
  su2double GetMu_RefND(void);
  
  /*!
   * \brief Get the value of the reference temperature for Sutherland model.
   * \return The reference temperature.
   */
  su2double GetMu_Temperature_Ref(void);

  /*!
   * \brief Get the value of the non-dimensional reference temperature for Sutherland model.
   * \return The non-dimensional reference temperature.
   */
  su2double GetMu_Temperature_RefND(void);
  
  /*!
   * \brief Get the value of the reference S for Sutherland model.
   * \return The reference S.
   */
  su2double GetMu_S(void);

  /*!
   * \brief Get the value of the non-dimensional reference S for Sutherland model.
   * \return The non-dimensional reference S.
   */
  su2double GetMu_SND(void);
  
  /*!
   * \brief Get the number of coefficients in the temperature polynomial models.
   * \return The the number of coefficients in the temperature polynomial models.
   */
  unsigned short GetnPolyCoeffs(void);
  
  /*!
   * \brief Get the temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for specific heat Cp.
   */
  su2double GetCp_PolyCoeff(unsigned short val_index);

  /*!
   * \brief Get the temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for specific heat Cp.
   */
  su2double GetCp_PolyCoeffND(unsigned short val_index);
  
  /*!
   * \brief Get the temperature polynomial coefficient for viscosity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for viscosity.
   */
  su2double GetMu_PolyCoeff(unsigned short val_index);
  
  /*!
   * \brief Get the temperature polynomial coefficient for viscosity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Non-dimensional temperature polynomial coefficient for viscosity.
   */
  su2double GetMu_PolyCoeffND(unsigned short val_index);
  
  /*!
   * \brief Get the temperature polynomial coefficients for viscosity.
   * \return Non-dimensional temperature polynomial coefficients for viscosity.
   */
  su2double* GetMu_PolyCoeffND(void);
  
  /*!
   * \brief Get the temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Temperature polynomial coefficient for thermal conductivity.
   */
  su2double GetKt_PolyCoeff(unsigned short val_index);
  
  /*!
   * \brief Get the temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   * \return Non-dimensional temperature polynomial coefficient for thermal conductivity.
   */
  su2double GetKt_PolyCoeffND(unsigned short val_index);
  
  /*!
   * \brief Get the temperature polynomial coefficients for thermal conductivity.
   * \return Non-dimensional temperature polynomial coefficients for thermal conductivity.
   */
  su2double* GetKt_PolyCoeffND(void);
  
  /*!
   * \brief Set the value of the non-dimensional constant viscosity.
   */
  void SetMu_ConstantND(su2double mu_const);
  
  /*!
   * \brief Set the value of the non-dimensional thermal conductivity.
   */
  void SetKt_ConstantND(su2double kt_const);
  
  /*!
   * \brief Set the value of the non-dimensional reference viscosity for Sutherland model.
   */
  void SetMu_RefND(su2double mu_ref);
  
  /*!
   * \brief Set the value of the non-dimensional reference temperature for Sutherland model.
   */
  void SetMu_Temperature_RefND(su2double mu_Tref);
  
  /*!
   * \brief Set the value of the non-dimensional S for Sutherland model.
   */
  void SetMu_SND(su2double mu_s);
  
  /*!
   * \brief Set the temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_coeff - Temperature polynomial coefficient for specific heat Cp.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   */
  void SetCp_PolyCoeffND(su2double val_coeff, unsigned short val_index);
  
  /*!
   * \brief Set the temperature polynomial coefficient for viscosity.
   * \param[in] val_coeff - Non-dimensional temperature polynomial coefficient for viscosity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   */
  void SetMu_PolyCoeffND(su2double val_coeff, unsigned short val_index);
  
  /*!
   * \brief Set the temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_coeff - Non-dimensional temperature polynomial coefficient for thermal conductivity.
   * \param[in] val_index - Index of the array with all polynomial coefficients.
   */
  void SetKt_PolyCoeffND(su2double val_coeff, unsigned short val_index);
  
  /*!
   * \brief Get the kind of method for computation of spatial gradients.
   * \return Numerical method for computation of spatial gradients.
   */
  unsigned short GetKind_Gradient_Method(void);
  
  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Linear_Solver(void);
  
  
  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Linear_Solver_Prec(void);
  
  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Deform_Linear_Solver(void);
  
  /*!
   * \brief Set the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  void SetKind_Deform_Linear_Solver_Prec(unsigned short val_kind_prec);
  
  /*!
   * \brief Set the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  void SetKind_Linear_Solver_Prec(unsigned short val_kind_prec);
  
  /*!
   * \brief Get min error of the linear solver for the implicit formulation.
   * \return Min error of the linear solver for the implicit formulation.
   */
  su2double GetLinear_Solver_Error(void);
  
  /*!
   * \brief Get min error of the linear solver for the implicit formulation.
   * \return Min error of the linear solver for the implicit formulation.
   */
  su2double GetDeform_Linear_Solver_Error(void);
  
  /*!
   * \brief Get max number of iterations of the linear solver for the implicit formulation.
   * \return Max number of iterations of the linear solver for the implicit formulation.
   */
  unsigned long GetLinear_Solver_Iter(void);
  
  /*!
   * \brief Get max number of iterations of the linear solver for the implicit formulation.
   * \return Max number of iterations of the linear solver for the implicit formulation.
   */
  unsigned long GetDeform_Linear_Solver_Iter(void);
  
  /*!
   * \brief Get the ILU fill-in level for the linear solver.
   * \return Fill in level of the ILU preconditioner for the linear solver.
   */
  unsigned short GetLinear_Solver_ILU_n(void);

  /*!
   * \brief Get restart frequency of the linear solver for the implicit formulation.
   * \return Restart frequency of the linear solver for the implicit formulation.
   */
  unsigned long GetLinear_Solver_Restart_Frequency(void);
  
  /*!
   * \brief Get the relaxation factor for iterative linear smoothers.
   * \return Relaxation factor.
   */
  su2double GetLinear_Solver_Smoother_Relaxation(void) const;
  
  /*!
   * \brief Get the relaxation coefficient of the linear solver for the implicit formulation.
   * \return relaxation coefficient of the linear solver for the implicit formulation.
   */
  su2double GetRelaxation_Factor_Flow(void);
  
  /*!
   * \brief Get the relaxation coefficient of the linear solver for the implicit formulation.
   * \return relaxation coefficient of the linear solver for the implicit formulation.
   */
  su2double GetRelaxation_Factor_AdjFlow(void);
  
  /*!
   * \brief Get the relaxation coefficient of the linear solver for the implicit formulation.
   * \return relaxation coefficient of the linear solver for the implicit formulation.
   */
  su2double GetRelaxation_Factor_Turb(void);

  /*!
   * \brief Get the relaxation coefficient of the CHT coupling.
   * \return relaxation coefficient of the CHT coupling.
   */
  su2double GetRelaxation_Factor_CHT(void);
  
  /*!
   * \brief Get the relaxation coefficient of the linear solver for the implicit formulation.
   * \return relaxation coefficient of the linear solver for the implicit formulation.
   */
  su2double GetRoe_Kappa(void);
  
  /*!
   * \brief Get the wing semi span.
   * \return value of the wing semi span.
   */
  su2double GetSemiSpan(void);
  
  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_AdjTurb_Linear_Solver(void);
  
  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_AdjTurb_Linear_Prec(void);
  
  /*!
   * \brief Get the kind of solver for the implicit solver.
   * \return Numerical solver for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_DiscAdj_Linear_Solver(void);
  
  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_DiscAdj_Linear_Prec(void);
  
  /*!
   * \brief Get the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  unsigned short GetKind_Deform_Linear_Solver_Prec(void);
  
  /*!
   * \brief Set the kind of preconditioner for the implicit solver.
   * \return Numerical preconditioner for implicit formulation (solving the linear system).
   */
  void SetKind_AdjTurb_Linear_Prec(unsigned short val_kind_prec);
  
  /*!
   * \brief Get min error of the linear solver for the implicit formulation.
   * \return Min error of the linear solver for the implicit formulation.
   */
  su2double GetAdjTurb_Linear_Error(void);
  
  /*!
   * \brief Get the entropy fix.
   * \return Vaule of the entropy fix.
   */
  su2double GetEntropyFix_Coeff(void);
  
  /*!
   * \brief Get max number of iterations of the linear solver for the implicit formulation.
   * \return Max number of iterations of the linear solver for the implicit formulation.
   */
  unsigned short GetAdjTurb_Linear_Iter(void);
  
  /*!
   * \brief Get CFL reduction factor for adjoint turbulence model.
   * \return CFL reduction factor.
   */
  su2double GetCFLRedCoeff_AdjTurb(void);
  
  /*!
   * \brief Get the number of linear smoothing iterations for mesh deformation.
   * \return Number of linear smoothing iterations for mesh deformation.
   */
  unsigned long GetGridDef_Linear_Iter(void);
  
  /*!
   * \brief Get the number of nonlinear increments for mesh deformation.
   * \return Number of nonlinear increments for mesh deformation.
   */
  unsigned long GetGridDef_Nonlinear_Iter(void);
  
  /*!
   * \brief Get information about writing grid deformation residuals to the console.
   * \return <code>TRUE</code> means that grid deformation residuals will be written to the console.
   */
  bool GetDeform_Output(void);
  
  /*!
   * \brief Get factor to multiply smallest volume for deform tolerance.
   * \return Factor to multiply smallest volume for deform tolerance.
   */
  su2double GetDeform_Coeff(void);
  
  /*!
   * \brief Get limit for the volumetric deformation.
   * \return Distance to the surface to be deformed.
   */
  su2double GetDeform_Limit(void);
  
  /*!
   * \brief Get Young's modulus for deformation (constant stiffness deformation)
   */
  su2double GetDeform_ElasticityMod(void);
  
  /*!
   * \brief Get Poisson's ratio for deformation (constant stiffness deformation)
   * \
   */
  su2double GetDeform_PoissonRatio(void);
  
  /*!
   * \brief Get the type of stiffness to impose for FEA mesh deformation.
   * \return type of stiffness to impose for FEA mesh deformation.
   */
  unsigned short GetDeform_Stiffness_Type(void);
  
  /*!
   * \brief Creates a tecplot file to visualize the volume deformation deformation made by the DEF software.
   * \return <code>TRUE</code> if the deformation is going to be plotted; otherwise <code>FALSE</code>.
   */
  bool GetVisualize_Volume_Def(void);
  
  /*!
   * \brief Creates a teot file to visualize the surface deformation deformation made by the DEF software.
   * \return <code>TRUE</code> if the deformation is going to be plotted; otherwise <code>FALSE</code>.
   */
  bool GetVisualize_Surface_Def(void);
  
  /*!
   * \brief Define the FFD box with a symetry plane.
   * \return <code>TRUE</code> if there is a symmetry plane in the FFD; otherwise <code>FALSE</code>.
   */
  bool GetFFD_Symmetry_Plane(void);
  
  /*!
   * \brief Get the kind of SU2 software component.
   * \return Kind of the SU2 software component.
   */
  unsigned short GetKind_SU2(void);
  
  /*!
   * \brief Get the kind of non-dimensionalization.
   * \return Kind of non-dimensionalization.
   */
  unsigned short GetRef_NonDim(void);
  
  /*!
   * \brief Get the kind of incompressible non-dimensionalization.
   * \return Kind of incompressible non-dimensionalization.
   */
  unsigned short GetRef_Inc_NonDim(void);

  /*!
   * \brief Get the kind of SU2 software component.
   * \return Kind of the SU2 software component.
   */
  void SetKind_SU2(unsigned short val_kind_su2);
  
  /*!
   * \brief Get the kind of the turbulence model.
   * \return Kind of the turbulence model.
   */
  unsigned short GetKind_Turb_Model(void);
  
  /*!
   * \brief Get the kind of the transition model.
   * \return Kind of the transion model.
   */
  unsigned short GetKind_Trans_Model(void);

  /*!
   * \brief Get the kind of the subgrid scale model.
   * \return Kind of the subgrid scale model.
   */
  unsigned short GetKind_SGS_Model(void);

  /*!
   * \brief Get the kind of adaptation technique.
   * \return Kind of adaptation technique.
   */
  unsigned short GetKind_Adaptation(void);
  
  /*!
   * \brief Get the number of new elements added in the adaptation process.
   * \return percentage of new elements that are going to be added in the adaptation.
   */
  su2double GetNew_Elem_Adapt(void);
  
  /*!
   * \brief Get the kind of time integration method.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of time integration method.
   */
  unsigned short GetKind_TimeIntScheme(void);
  
  /*!
   * \brief Get the kind of convective numerical scheme.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of the convective scheme.
   */
  unsigned short GetKind_ConvNumScheme(void);
  
  /*!
   * \brief Get kind of center scheme for the convective terms.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of center scheme for the convective terms.
   */
  unsigned short GetKind_Centered(void);
  
  /*!
   * \brief Get kind of upwind scheme for the convective terms.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return Kind of upwind scheme for the convective terms.
   */
  unsigned short GetKind_Upwind(void);
  
  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL(void);
  
  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Flow(void);
  
  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Heat(void);

  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_Turb(void);
  
  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_AdjFlow(void);
  
  /*!
   * \brief Get if the upwind scheme used MUSCL or not.
   * \note This is the information that the code will use, the method will
   *       change in runtime depending of the specific equation (direct, adjoint,
   *       linearized) that is being solved.
   * \return MUSCL scheme.
   */
  bool GetMUSCL_AdjTurb(void);
  
  /*!
   * \brief Get whether to "Use Accurate Jacobians" for AUSM+up(2) and SLAU(2).
   * \return yes/no.
   */
  inline bool GetUse_Accurate_Jacobians(void) { return Use_Accurate_Jacobians; }

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the flow equations.
   */
  unsigned short GetKind_TimeIntScheme_Flow(void);

  /*!
   * \brief Get the kind of scheme (aliased or non-aliased) to be used in the
   *        predictor step of ADER-DG.
   * \return Kind of scheme used in the predictor step of ADER-DG.
   */
  unsigned short GetKind_ADER_Predictor(void);

  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the plasma equations.
   */
  unsigned short GetKind_TimeIntScheme_Heat(void);
  
  /*!
   * \brief Get the kind of time stepping
   *        for the heat equation.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of time stepping for the heat equation.
   */
  unsigned short GetKind_TimeStep_Heat(void);
  
  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the plasma equations.
   */
  unsigned short GetKind_TimeIntScheme_FEA(void);
  
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
  unsigned short GetKind_SpaceIteScheme_FEA(void);
  
  /*!
   * \brief Get the kind of convective numerical scheme for the flow
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_ConvNumScheme_Flow(void);

  /*!
   * \brief Get the kind of convective numerical scheme for the flow
   *        equations (finite element).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_ConvNumScheme_FEM_Flow(void);

  /*!
   * \brief Get the kind of convective numerical scheme for the template
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_ConvNumScheme_Template(void);
  
  /*!
   * \brief Get the kind of center convective numerical scheme for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_Centered_Flow(void);
  
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
  unsigned short GetKind_Upwind_Flow(void);

  /*!
   * \brief Get the kind of finite element convective numerical scheme for the flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of finite element convective numerical scheme for the flow equations.
   */
  unsigned short GetKind_FEM_Flow(void);

  /*!
   * \brief Get the kind of shock capturing method in FEM DG solver.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of shock capturing method in FEM DG solver.
   */
  unsigned short GetKind_FEM_DG_Shock(void);

  /*!
   * \brief Get the kind of matrix coloring used for the sparse Jacobian computation.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of matrix coloring used.
   */
  unsigned short GetKind_Matrix_Coloring(void);

  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients.
   */
  unsigned short GetKind_SlopeLimit(void);
  
  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the flow equations.
   */
  unsigned short GetKind_SlopeLimit_Flow(void);
  
  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the turbulent equation.
   */
  unsigned short GetKind_SlopeLimit_Turb(void);
  
  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the adjoint turbulent equation.
   */
  unsigned short GetKind_SlopeLimit_AdjTurb(void);
  
  /*!
   * \brief Get the method for limiting the spatial gradients.
   * \return Method for limiting the spatial gradients solving the adjoint flow equation.
   */
  unsigned short GetKind_SlopeLimit_AdjFlow(void);
  
  /*!
   * \brief Value of the calibrated constant for the Lax method (center scheme).
   * \note This constant is used in coarse levels and with first order methods.
   * \return Calibrated constant for the Lax method.
   */
  su2double GetKappa_1st_Flow(void);
  
  /*!
   * \brief Value of the calibrated constant for the JST method (center scheme).
   * \return Calibrated constant for the JST method for the flow equations.
   */
  su2double GetKappa_2nd_Flow(void);
  
  /*!
   * \brief Value of the calibrated constant for the JST method (center scheme).
   * \return Calibrated constant for the JST method for the flow equations.
   */
  su2double GetKappa_4th_Flow(void);

  /*!
   * \brief Value of the calibrated constant for the JST method (center scheme).
   * \return Calibrated constant for the JST-like method for the heat equations.
   */
  su2double GetKappa_2nd_Heat(void);

  /*!
   * \brief Value of the calibrated constant for the JST-like method (center scheme).
   * \return Calibrated constant for the JST-like method for the heat equation.
   */
  su2double GetKappa_4th_Heat(void);
  
  /*!
   * \brief Factor by which to multiply the dissipation contribution to Jacobians of central schemes.
   * \return The factor.
   */
  inline su2double GetCent_Jac_Fix_Factor(void) { return Cent_Jac_Fix_Factor; }
  
  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the adjoint flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the adjoint flow equations.
   */
  unsigned short GetKind_TimeIntScheme_AdjFlow(void);
  
  /*!
   * \brief Get the kind of convective numerical scheme for the adjoint flow
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the adjoint flow equations.
   */
  unsigned short GetKind_ConvNumScheme_AdjFlow(void);
  
  /*!
   * \brief Get the kind of center convective numerical scheme for the adjoint flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the adjoint flow equations.
   */
  unsigned short GetKind_Centered_AdjFlow(void);
  
  /*!
   * \brief Get the kind of upwind convective numerical scheme for the adjoint flow equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the adjoint flow equations.
   */
  unsigned short GetKind_Upwind_AdjFlow(void);
  
  /*!
   * \brief Value of the calibrated constant for the high order method (center scheme).
   * \return Calibrated constant for the high order center method for the adjoint flow equations.
   */
  su2double GetKappa_2nd_AdjFlow(void);
  
  /*!
   * \brief Value of the calibrated constant for the high order method (center scheme).
   * \return Calibrated constant for the high order center method for the adjoint flow equations.
   */
  su2double GetKappa_4th_AdjFlow(void);
  
  /*!
   * \brief Value of the calibrated constant for the low order method (center scheme).
   * \return Calibrated constant for the low order center method for the adjoint flow equations.
   */
  su2double GetKappa_1st_AdjFlow(void);
  
  /*!
   * \brief Get the kind of integration scheme (implicit)
   *        for the turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the turbulence equations.
   */
  unsigned short GetKind_TimeIntScheme_Turb(void);
  
  /*!
   * \brief Get the kind of convective numerical scheme for the turbulence
   *        equations (upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the turbulence equations.
   */
  unsigned short GetKind_ConvNumScheme_Turb(void);
  
  /*!
   * \brief Get the kind of center convective numerical scheme for the turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the turbulence equations.
   */
  unsigned short GetKind_Centered_Turb(void);
  
  /*!
   * \brief Get the kind of upwind convective numerical scheme for the turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the turbulence equations.
   */
  unsigned short GetKind_Upwind_Turb(void);
  
  /*!
   * \brief Get the kind of integration scheme (explicit or implicit)
   *        for the adjoint turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of integration scheme for the adjoint turbulence equations.
   */
  unsigned short GetKind_TimeIntScheme_AdjTurb(void);
  
  /*!
   * \brief Get the kind of convective numerical scheme for the adjoint turbulence
   *        equations (centered or upwind).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the adjoint turbulence equations.
   */
  unsigned short GetKind_ConvNumScheme_AdjTurb(void);
  
  /*!
   * \brief Get the kind of convective numerical scheme for the heat equation.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of convective numerical scheme for the heat equation.
   */
  unsigned short GetKind_ConvNumScheme_Heat(void);

  /*!
   * \brief Get the kind of center convective numerical scheme for the adjoint turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of center convective numerical scheme for the adjoint turbulence equations.
   */
  unsigned short GetKind_Centered_AdjTurb(void);
  
  /*!
   * \brief Get the kind of upwind convective numerical scheme for the adjoint turbulence equations.
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of upwind convective numerical scheme for the adjoint turbulence equations.
   */
  unsigned short GetKind_Upwind_AdjTurb(void);
  
  /*!
   * \brief Provides information about the way in which the turbulence will be treated by the
   *        cont. adjoint method.
   * \return <code>FALSE</code> means that the adjoint turbulence equations will be used.
   */
  bool GetFrozen_Visc_Cont(void);
  
  /*!
   * \brief Provides information about the way in which the turbulence will be treated by the
   *        disc. adjoint method.
   * \return <code>FALSE</code> means that the adjoint turbulence equations will be used.
   */
  bool GetFrozen_Visc_Disc(void);
  
  /*!
   * \brief Provides information about using an inconsistent (primal/dual) discrete adjoint formulation
   * \return <code>FALSE</code> means that the adjoint use the same numerical methods than the primal problem.
   */
  bool GetInconsistent_Disc(void);

  /*!
   * \brief Provides information about the way in which the limiter will be treated by the
   *        disc. adjoint method.
   * \return <code>FALSE</code> means that the limiter computation is included.
   */
  bool GetFrozen_Limiter_Disc(void);
  
  /*!
   * \brief Write convergence file for FSI problems
   * \return <code>FALSE</code> means no file is written.
   */
  bool GetWrite_Conv_FSI(void);
  
  /*!
   * \brief Provides information about if the sharp edges are going to be removed from the sensitivity.
   * \return <code>FALSE</code> means that the sharp edges will be removed from the sensitivity.
   */
  bool GetSens_Remove_Sharp(void);
  
  /*!
   * \brief Get the kind of inlet boundary condition treatment (total conditions or mass flow).
   * \return Kind of inlet boundary condition.
   */
  unsigned short GetKind_Inlet(void);
  
  /*!
   * \brief Check if the inlet profile(s) are specified in an input file
   * \return True if an input file is to be used for the inlet profile(s)
   */
  bool GetInlet_Profile_From_File(void);

  /*!
   * \brief Get name of the input file for the specified inlet profile.
   * \return Name of the input file for the specified inlet profile.
   */
  string GetInlet_FileName(void);

  /*!
   * \brief Get the tolerance used for matching two points on a specified inlet
   * \return Tolerance used for matching a point to a specified inlet
   */
  su2double GetInlet_Profile_Matching_Tolerance(void);
  
  /*!
   * \brief Get the type of incompressible inlet from the list.
   * \return Kind of the incompressible inlet.
   */
  unsigned short GetKind_Inc_Inlet(string val_marker);

  /*!
   * \brief Get the total number of types in Kind_Inc_Inlet list
   * \return Total number of types in Kind_Inc_Inlet list
   */
  unsigned short GetnInc_Inlet(void);

  /*!
   * \brief Flag for whether the local boundary normal is used as the flow direction for an incompressible pressure inlet.
   * \return <code>FALSE</code> means the prescribed flow direction is used.
   */
  bool GetInc_Inlet_UseNormal(void);

  /*!
   * \brief Get the type of incompressible outlet from the list.
   * \return Kind of the incompressible outlet.
   */
  unsigned short GetKind_Inc_Outlet(string val_marker);
  
  /*!
   * \brief Get the damping factor applied to velocity updates at incompressible pressure inlets.
   * \return Damping factor applied to velocity updates at incompressible pressure inlets.
   */
  su2double GetInc_Inlet_Damping(void);
  
  /*!
   * \brief Get the damping factor applied to pressure updates at incompressible mass flow outlet.
   * \return Damping factor applied to pressure updates at incompressible mass flow outlet.
   */
  su2double GetInc_Outlet_Damping(void);
  
  /*!
   * \brief Get the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  unsigned short GetKind_AverageProcess(void);

  /*!
   * \brief Get the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  unsigned short GetKind_PerformanceAverageProcess(void);

  /*!
   * \brief Set the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  void SetKind_AverageProcess(unsigned short new_AverageProcess);

  /*!
   * \brief Set the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  void SetKind_PerformanceAverageProcess(unsigned short new_AverageProcess);

  /*!
   * \brief Get coeff for Rotating Frame Ramp.
   * \return coeff Ramp Rotating Frame.
   */
  su2double GetRampRotatingFrame_Coeff(unsigned short iCoeff);

  /*!
   * \brief Get Rotating Frame Ramp option.
   * \return Ramp Rotating Frame option.
   */
  bool GetRampRotatingFrame(void);

  /*!
   * \brief Get coeff for Outlet Pressure Ramp.
   * \return coeff Ramp Outlet Pressure.
   */
  su2double GetRampOutletPressure_Coeff(unsigned short iCoeff);

  /*!
   * \brief Get final Outlet Pressure value for the ramp.
   * \return final Outlet Pressure value.
   */
  su2double GetFinalOutletPressure(void);

  /*!
   * \brief Get final Outlet Pressure value for the ramp.
   * \return Monitor Outlet Pressure value.
   */
  su2double GetMonitorOutletPressure(void);

  /*!
   * \brief Set Monitor Outlet Pressure value for the ramp.
   */
  void SetMonitotOutletPressure(su2double newMonPres);

  /*!
   * \brief Get Outlet Pressure Ramp option.
   * \return Ramp Outlet pressure option.
   */
  bool GetRampOutletPressure(void);

  /*!
   * \brief Get mixedout coefficients.
   * \return mixedout coefficient.
   */
  su2double GetMixedout_Coeff(unsigned short iCoeff);

  /*!
   * \brief Get extra relaxation factor coefficients for the Giels BC.
   * \return mixedout coefficient.
   */
  su2double GetExtraRelFacGiles(unsigned short iCoeff);

  /*!
   * \brief Get mach limit for average massflow-based procedure .
   * \return mach limit.
   */
  su2double GetAverageMachLimit(void);

  /*!
   * \brief Get the kind of mixing process for averaging quantities at the boundaries.
   * \return Kind of mixing process.
   */
  unsigned short GetKind_MixingPlaneInterface(void);

  /*!
   * \brief Get the kind of turbomachinery architecture.
   * \return Kind of turbomachinery architecture.
   */
  unsigned short GetKind_TurboMachinery(unsigned short val_iZone);

  /*!
   * \brief Get the kind of turbomachinery architecture.
   * \return Kind of turbomachinery architecture.
   */
  unsigned short GetKind_SpanWise(void);
  
  /*!
   * \brief Verify if there is mixing plane interface specified from config file.
   * \return boolean.
   */
  bool GetBoolMixingPlaneInterface(void);

  /*!
   * \brief Verify if there is mixing plane interface specified from config file.
   * \return boolean.
   */
  bool GetBoolTurbMixingPlane(void);

  /*!
   * \brief Verify if there is mixing plane interface specified from config file.
   * \return boolean.
   */
  bool GetSpatialFourier(void);

  /*!
   * \brief number mixing plane interface specified from config file.
   * \return number of bound.
   */
  unsigned short GetnMarker_MixingPlaneInterface(void);
  
  /*!
   * \brief Verify if there is Turbomachinery performance option specified from config file.
   * \return boolean.
   */
  bool GetBoolTurbomachinery(void);

  /*!
   * \brief Verify if there are zone specific solvers entered in the config file.
   * \return boolean.
   */
  bool GetBoolZoneSpecific(void);
  
  /*!
   * \brief number Turbomachinery blades computed using the pitch information.
   * \return nBlades.
   */
  su2double GetnBlades(unsigned short val_iZone);

  /*!
   * \brief number Turbomachinery blades computed using the pitch information.
   * \return nBlades.
   */
  void SetnBlades(unsigned short val_iZone, su2double nblades);

  /*!
   * \brief Verify if there is any Giles Boundary Condition option specified from config file.
   * \return boolean.
   */
  bool GetBoolGiles(void);
  
  /*!
   * \brief Verify if there is any Riemann Boundary Condition option specified from config file.
   * \return boolean.
   */
  bool GetBoolRiemann(void);

  /*!
   * \brief number Turbomachinery performance option specified from config file.
   * \return number of bound.
   */
  unsigned short GetnMarker_Turbomachinery(void);

  /*!
   * \brief Get number of shroud markers.
   * \return number of marker shroud.
   */
  unsigned short GetnMarker_Shroud(void);

  /*!
   * \brief Get the marker shroud.
   * \return marker shroud.
   */
  string GetMarker_Shroud(unsigned short val_marker);

  /*!
   * \brief number Turbomachinery performance option specified from config file.
   * \return number of bound.
   */
  unsigned short GetnMarker_TurboPerformance(void);

  /*!
   * \brief number span-wise sections to compute 3D BC and performance for turbomachinery specified by the user.
   * \return number of span-wise sections.
   */
  unsigned short Get_nSpanWiseSections_User(void);

  /*!
   * \brief number span-wise sections to compute 3D BC and performance for turbomachinery.
   * \return number of span-wise sections.
   */
  unsigned short GetnSpanWiseSections(void);

  /*!
   * \brief set number of maximum span-wise sections among all zones .
   */
  void SetnSpanMaxAllZones(unsigned short val_nSpna_max);

  /*!
   * \brief number span-wise sections to compute performance for turbomachinery.
   * \return number of max span-wise sections.
   */
  unsigned short GetnSpanMaxAllZones(void);
	
  /*!
   * \brief set number span-wise sections to compute 3D BC and performance for turbomachinery.
   */
  void SetnSpanWiseSections(unsigned short nSpan);

  /*!
   * \brief set number span-wise sections to compute 3D BC and performance for turbomachinery.
   */
  unsigned short GetnSpan_iZones(unsigned short iZone);

  /*!
   * \brief set number span-wise sections to compute 3D BC and performance for turbomachinery.
   */
  void SetnSpan_iZones(unsigned short nSpan, unsigned short iZone);

  /*!
   * \brief get inlet bounds name for Turbomachinery performance calculation.
   * \return name of the bound.
   */
  string GetMarker_TurboPerf_BoundIn(unsigned short index);
  
  /*!
   * \brief get outlet bounds name for Turbomachinery performance calculation.
   * \return name of the bound.
   */
  string GetMarker_TurboPerf_BoundOut(unsigned short index);
  
  /*!
   * \brief get marker kind for Turbomachinery performance calculation.
   * \return kind index.
   */
  unsigned short GetKind_TurboPerf(unsigned short index);
  
  /*!
   * \brief get outlet bounds name for Turbomachinery performance calculation.
   * \return name of the bound.
   */
  string GetMarker_PerBound(unsigned short val_marker);
  
  /*!
   * \brief Get the kind of inlet boundary condition treatment (total conditions or mass flow).
   * \return Kind of inlet boundary condition.
   */
  unsigned short GetKind_Engine_Inflow(void);
  
  /*!
   * \brief Get the kind of inlet boundary condition treatment (total conditions or mass flow).
   * \return Kind of inlet boundary condition.
   */
  unsigned short GetKind_ActDisk(void);
  
  /*!
   * \brief Get the number of sections.
   * \return Number of sections
   */
  unsigned short GetnLocationStations(void);
  
  /*!
   * \brief Get the number of sections for computing internal volume.
   * \return Number of sections for computing internal volume.
   */
  unsigned short GetnWingStations(void);
  
  /*!
   * \brief Get the location of the waterline.
   * \return Z location of the waterline.
   */
  su2double GetGeo_Waterline_Location(void);
  
  /*!
   * \brief Provides information about the the nodes that are going to be moved on a deformation
   *        volumetric grid deformation.
   * \return <code>TRUE</code> means that only the points on the FFD box will be moved.
   */
  bool GetHold_GridFixed(void);
  
  /*!
   * \brief Get the kind of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \return Kind of objective function.
   */
  unsigned short GetKind_ObjFunc(void);
  
  /*!
   * \author H. Kline
   * \brief Get the kind of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \return Kind of objective function.
   */
  unsigned short GetKind_ObjFunc(unsigned short val_obj);
  
  /*!
   * \author H. Kline
   * \brief Get the weight of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \return Weight of objective function.
   */
  su2double GetWeight_ObjFunc(unsigned short val_obj);
  
  /*!
   * \author H. Kline
   * \brief Set the weight of objective function. There are several options: Drag coefficient,
   *        Lift coefficient, efficiency, etc.
   * \note The objective function will determine the boundary condition of the adjoint problem.
   * \return Weight of objective function.
   */
  void SetWeight_ObjFunc(unsigned short val_obj, su2double val);
  
  /*!
   * \author H. Kline
   * \brief Get the coefficients of the objective defined by the chain rule with primitive variables.
   * \note This objective is only applicable to gradient calculations. Objective value must be
   * calculated using the area averaged outlet values of density, velocity, and pressure.
   * Gradients are w.r.t density, velocity[3], and pressure. when 2D gradient w.r.t. 3rd component of velocity set to 0.
   */
  su2double GetCoeff_ObjChainRule(unsigned short iVar);
  
  /*!
   * \author H. Kline
   * \brief Get the flag indicating whether to comput a combined objective.
   */
  bool GetComboObj(void);
  
  /*!
   * \brief Get the kind of sensitivity smoothing technique.
   * \return Kind of sensitivity smoothing technique.
   */
  unsigned short GetKind_SensSmooth(void);
  
  /*!
   * \brief Provides information about the time integration, and change the write in the output
   *        files information about the iteration.
   * \return The kind of time integration: Steady state, time stepping method (unsteady) or
   *         dual time stepping method (unsteady).
   */
  unsigned short GetUnsteady_Simulation(void);
  
  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  unsigned short GetnReactions(void);
  
  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  su2double GetArrheniusCoeff(unsigned short iReaction);
  
  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  su2double GetArrheniusEta(unsigned short iReaction);
  
  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  su2double GetArrheniusTheta(unsigned short iReaction);
  
  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  su2double* GetRxnTcf_a(void);
  
  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  su2double* GetRxnTcf_b(void);
  
  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  su2double* GetRxnTcb_a(void);
  
  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  su2double* GetRxnTcb_b(void);
  
  /*!
   * \brief Dissociation potential of species.
   * \return: Dissociation potential.
   */
  su2double* GetDissociationPot(void);
  
  /*!
   * \brief Provides the number of rotational modes of energy storage
   * \return: Vector of rotational mode count
   */
  su2double* GetRotationModes(void);
  
  /*!
   * \brief Provides the characteristic vibrational temperature for calculating e_vib
   * \return: Vector of characteristic vibrational temperatures [K]
   */
  su2double* GetCharVibTemp(void);
  
  /*!
   * \brief Provides the characteristic electronic temperature for calculating e_el
   * \return: Vector of characteristic vibrational temperatures [K]
   */
  su2double** GetCharElTemp(void);
  
  /*!
   * \brief Provides the degeneracy of electron states for calculating e_el
   * \return: Vector of characteristic vibrational temperatures [K]
   */
  su2double** GetElDegeneracy(void);
  
  /*!
   * \brief Provides number electron states for calculating e_el
   * \return: Vector of number of electron states for each species
   */
  unsigned short* GetnElStates(void);
  
  
  /*!
   * \brief Provides the thermodynamic reference temperatures from the JANAF tables
   * \return: Vector of reference temperatures [K]
   */
  su2double* GetRefTemperature(void);
  
  /*!
   * \brief Provides the characteristic vibrational temperature for calculating e_vib
   * \return: The number of chemical reactions, read from input file
   */
  su2double GetCharVibTemp(unsigned short iSpecies);
  
  /*!
   * \brief Provides the molar mass of each species present in multi species fluid
   * \return: Vector of molar mass of each species in kg/kmol
   */
  su2double* GetMolar_Mass(void);
  
  /*!
   * \brief Provides the molar mass of each species present in multi species fluid
   * \return: Mass of each species in Kg
   */
  su2double GetMolar_Mass(unsigned short iSpecies);
  
  /*!
   * \brief Retrieves the number of monatomic species in the multicomponent gas.
   * \return: Number of monatomic species.
   */
  unsigned short GetnMonatomics(void);
  
  /*!
   * \brief Retrieves the number of monatomic species in the multicomponent gas.
   * \return: Number of monatomic species.
   */
  unsigned short GetnDiatomics(void);
  
  /*!
   * \brief Provides the molar mass of each species present in multi species fluid
   * \return: Molar mass of the specified gas consituent [kg/kmol]
   */
  su2double GetInitial_Gas_Composition(unsigned short iSpecies);
  
  /*!
   * \brief Provides the formation enthalpy of the specified species at standard conditions
   * \return: Enthalpy of formation
   */
  su2double* GetEnthalpy_Formation(void);
  
  /*!
   * \brief Provides the formation enthalpy of the specified species at standard conditions
   * \return: Enthalpy of formation
   */
  su2double GetEnthalpy_Formation(unsigned short iSpecies);
  
  /*!
   * \brief Provides the restart information.
   * \return Restart information, if <code>TRUE</code> then the code will use the solution as restart.
   */
  bool GetRestart(void);

  /*!
   * \brief Flag for whether binary SU2 native restart files are written.
   * \return Flag for whether binary SU2 native restart files are written, if <code>TRUE</code> then the code will output binary restart files.
   */
  bool GetWrt_Binary_Restart(void);

  /*!
   * \brief Flag for whether binary SU2 native restart files are read.
   * \return Flag for whether binary SU2 native restart files are read, if <code>TRUE</code> then the code will load binary restart files.
   */
  bool GetRead_Binary_Restart(void);

  /*!
   * \brief Provides the number of varaibles.
   * \return Number of variables.
   */
  unsigned short GetnVar(void);
  
  /*!
   * \brief Provides the number of varaibles.
   * \return Number of variables.
   */
  unsigned short GetnZone(void);
  
  /*!
   * \brief Provides the number of varaibles.
   * \return Number of variables.
   */
  unsigned short GetiZone(void);
  
  /*!
   * \brief For some problems like adjoint or the linearized equations it
   *		  is necessary to restart the flow solution.
   * \return Flow restart information, if <code>TRUE</code> then the code will restart the flow solution.
   */
  
  bool GetRestart_Flow(void);
  
  /*!
   * \brief Indicates whether electron gas is present in the gas mixture.
   */
  bool GetIonization(void);
  
  /*!
   * \brief Information about computing and plotting the equivalent area distribution.
   * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are computing the equivalent area.
   */
  bool GetEquivArea(void);
  
  /*!
   * \brief Information about computing and plotting the equivalent area distribution.
   * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are computing the equivalent area.
   */
  bool GetInvDesign_Cp(void);
  
  /*!
   * \brief Information about computing and plotting the equivalent area distribution.
   * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are computing the equivalent area.
   */
  bool GetInvDesign_HeatFlux(void);
  
  /*!
   * \brief Get name of the input grid.
   * \return File name of the input grid.
   */
  string GetMesh_FileName(void);
  
  /*!
   * \brief Get name of the output grid, this parameter is important for grid
   *        adaptation and deformation.
   * \return File name of the output grid.
   */
  string GetMesh_Out_FileName(void);
  
  /*!
   * \brief Get the name of the file with the solution of the flow problem.
   * \return Name of the file with the solution of the flow problem.
   */
  string GetSolution_FlowFileName(void);
  
  /*!
   * \brief Get the name of the file with the solution of the adjoint flow problem
   *		  with drag objective function.
   * \return Name of the file with the solution of the adjoint flow problem with
   *         drag objective function.
   */
  string GetSolution_AdjFileName(void);
  
  /*!
   * \brief Get the name of the file with the solution of the structural problem.
   * \return Name of the file with the solution of the structural problem.
   */
  string GetSolution_FEMFileName(void);
  
  /*!
   * \brief Get the name of the file with the solution of the adjoint structural problem.
   * \return Name of the file with the solution of the structural problem.
   */
  string GetSolution_AdjFEMFileName(void);
  
  /*!
   * \brief Get the name of the file with the residual of the problem.
   * \return Name of the file with the residual of the problem.
   */
  string GetResidual_FileName(void);
  
  /*!
   * \brief Get the format of the input/output grid.
   * \return Format of the input/output grid.
   */
  unsigned short GetMesh_FileFormat(void);
  
  /*!
   * \brief Get the format of the output solution.
   * \return Format of the output solution.
   */
  unsigned short GetOutput_FileFormat(void);
  
  /*!
   * \brief Get the format of the output solution.
   * \return Format of the output solution.
   */
  unsigned short GetActDisk_Jump(void);
  
  /*!
   * \brief Get the name of the file with the convergence history of the problem.
   * \return Name of the file with convergence history of the problem.
   */
  string GetConv_FileName(void);
  
  /*!
   * \brief Get the name of the file with the convergence history of the problem for FSI applications.
   * \return Name of the file with convergence history of the problem.
   */
  string GetConv_FileName_FSI(void);
  
  /*!
   * \brief Get the name of the file with the forces breakdown of the problem.
   * \return Name of the file with forces breakdown of the problem.
   */
  string GetBreakdown_FileName(void);
  
  /*!
   * \brief Get the name of the file with the flow variables.
   * \return Name of the file with the primitive variables.
   */
  string GetFlow_FileName(void);
  
  /*!
   * \brief Get the name of the file with the structure variables.
   * \return Name of the file with the structure variables.
   */
  string GetStructure_FileName(void);
  
  /*!
   * \brief Get the name of the file with the structure variables.
   * \return Name of the file with the structure variables.
   */
  string GetSurfStructure_FileName(void);
  
  /*!
   * \brief Get the name of the file with the adjoint structure variables.
   * \return Name of the file with the adjoint structure variables.
   */
  string GetAdjStructure_FileName(void);
  
  /*!
   * \brief Get the name of the file with the adjoint structure variables.
   * \return Name of the file with the adjoint structure variables.
   */
  string GetAdjSurfStructure_FileName(void);
  
  /*!
   * \brief Get the name of the file with the structure variables.
   * \return Name of the file with the structure variables.
   */
  string GetSurfHeat_FileName(void);
  
  /*!
   * \brief Get the name of the file with the wave variables.
   * \return Name of the file with the wave variables.
   */
  string GetHeat_FileName(void);
  
  /*!
   * \brief Get the name of the restart file for the heat variables.
   * \return Name of the restart file for the flow variables.
   */
  string GetRestart_HeatFileName(void);
  
  /*!
   * \brief Append the zone index to the restart or the solution files.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultizone_FileName(string val_filename, int val_iZone);

  /*!
   * \brief Append the zone index to the restart or the solution files.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultizone_HistoryFileName(string val_filename, int val_iZone);
  
  /*!
   * \brief Append the instance index to the restart or the solution files.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultiInstance_FileName(string val_filename, int val_iInst);

  /*!
   * \brief Append the instance index to the restart or the solution files.
   * \return Name of the restart file for the flow variables.
   */
  string GetMultiInstance_HistoryFileName(string val_filename, int val_iInst);

  /*!
   * \brief Get the name of the restart file for the flow variables.
   * \return Name of the restart file for the flow variables.
   */
  string GetRestart_FlowFileName(void);
  
  /*!
   * \brief Get the name of the restart file for the adjoint variables (drag objective function).
   * \return Name of the restart file for the adjoint variables (drag objective function).
   */
  string GetRestart_AdjFileName(void);
  
  /*!
   * \brief Get the name of the restart file for the structural variables.
   * \return Name of the restart file for the structural variables.
   */
  string GetRestart_FEMFileName(void);
  
  /*!
   * \brief Get the name of the restart file for the structural adjoint variables.
   * \return Name of the restart file for the structural adjoint variables.
   */
  string GetRestart_AdjFEMFileName(void);
  
  /*!
   * \brief Get the name of the file with the adjoint variables.
   * \return Name of the file with the adjoint variables.
   */
  string GetAdj_FileName(void);
  
  /*!
   * \brief Get the name of the file with the gradient of the objective function.
   * \return Name of the file with the gradient of the objective function.
   */
  string GetObjFunc_Grad_FileName(void);
  
  /*!
   * \brief Get the name of the file with the gradient of the objective function.
   * \return Name of the file with the gradient of the objective function.
   */
  string GetObjFunc_Value_FileName(void);
  
  /*!
   * \brief Get the name of the file with the surface information for the flow problem.
   * \return Name of the file with the surface information for the flow problem.
   */
  string GetSurfFlowCoeff_FileName(void);
  
  /*!
   * \brief Get the name of the file with the surface information for the adjoint problem.
   * \return Name of the file with the surface information for the adjoint problem.
   */
  string GetSurfAdjCoeff_FileName(void);
  
  /*!
   * \brief Get the name of the file with the surface sensitivity (discrete adjoint).
   * \return Name of the file with the surface sensitivity (discrete adjoint).
   */
  string GetSurfSens_FileName(void);
  
  /*!
   * \brief Get the name of the file with the volume sensitivity (discrete adjoint).
   * \return Name of the file with the volume sensitivity (discrete adjoint).
   */
  string GetVolSens_FileName(void);
  
  /*!
   * \brief Augment the input filename with the iteration number for an unsteady file.
   * \param[in] val_filename - String value of the base filename.
   * \param[in] val_iter - Unsteady iteration number or time instance.
   * \return Name of the file with the iteration number for an unsteady solution file.
   */
  string GetUnsteady_FileName(string val_filename, int val_iter);
  
  /*!
   * \brief Append the input filename string with the appropriate objective function extension.
   * \param[in] val_filename - String value of the base filename.
   * \return Name of the file with the appropriate objective function extension.
   */
  string GetObjFunc_Extension(string val_filename);
  
  /*!
   * \brief Get the criteria for structural residual (relative/absolute).
   * \return Relative/Absolute criteria for structural convergence.
   */
  unsigned short GetResidual_Criteria_FEM(void);
  
  /*!
   * \brief Get functional that is going to be used to evaluate the residual flow convergence.
   * \return Functional that is going to be used to evaluate the residual flow convergence.
   */
  unsigned short GetResidual_Func_Flow(void);
  
  /*!
   * \brief Get functional that is going to be used to evaluate the flow convergence.
   * \return Functional that is going to be used to evaluate the flow convergence.
   */
  unsigned short GetCauchy_Func_Flow(void);
  
  /*!
   * \brief Get functional that is going to be used to evaluate the adjoint flow convergence.
   * \return Functional that is going to be used to evaluate the adjoint flow convergence.
   */
  unsigned short GetCauchy_Func_AdjFlow(void);
  
  /*!
   * \brief Get the number of iterations that are considered in the Cauchy convergence criteria.
   * \return Number of elements in the Cauchy criteria.
   */
  unsigned short GetCauchy_Elems(void);
  
  /*!
   * \brief Get the number of iterations that are not considered in the convergence criteria.
   * \return Number of iterations before starting with the convergence criteria.
   */
  unsigned long GetStartConv_Iter(void);
  
  /*!
   * \brief Get the value of convergence criteria for the Cauchy method in the direct,
   *        adjoint or linearized problem.
   * \return Value of the convergence criteria.
   */
  su2double GetCauchy_Eps(void);
  
  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetDelta_UnstTimeND(void);
  
  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetTotal_UnstTimeND(void);
  
  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation.
   */
  su2double GetDelta_UnstTime(void);
  
  /*!
   * \brief Set the value of the unsteadty time step using the CFL number.
   * \param[in] val_delta_unsttimend - Value of the unsteady time step using CFL number.
   */
  void SetDelta_UnstTimeND(su2double val_delta_unsttimend);
  
  /*!
   * \brief If we are performing an unsteady simulation, this is the
   * 	value of max physical time for which we run the simulation
   * \return Value of the physical time in an unsteady simulation.
   */
  su2double GetTotal_UnstTime(void);
  
  /*!
   * \brief If we are performing an unsteady simulation, this is the
   * 	value of current time.
   * \return Value of the physical time in an unsteady simulation.
   */
  su2double GetCurrent_UnstTime(void);
  
  /*!
   * \brief Divide the rectbles and hexahedron.
   * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
   */
  bool GetSubsonicEngine(void);
  
  /*!
   * \brief Actuator disk defined with a double surface.
   * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
   */
  bool GetActDisk_DoubleSurface(void);
  
  /*!
   * \brief Only halg of the engine is in the compputational grid.
   * \return <code>TRUE</code> if the engine is complete; otherwise <code>FALSE</code>.
   */
  bool GetEngine_HalfModel(void);
  
  /*!
   * \brief Actuator disk defined with a double surface.
   * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
   */
  bool GetActDisk_SU2_DEF(void);
  
  /*!
   * \brief Value of the design variable step, we use this value in design problems.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \param[in] val_value - Value of the design variable that we want to read.
   * \return Design variable step.
   */
  su2double GetDV_Value(unsigned short val_dv, unsigned short val_val = 0);
  
  /*!
   * \brief Set the value of the design variable step, we use this value in design problems.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \param[in] val    - Value of the design variable.
   */
  void SetDV_Value(unsigned short val_dv, unsigned short val_ind, su2double val);
  
  /*!
   * \brief Get information about the grid movement.
   * \return <code>TRUE</code> if there is a grid movement; otherwise <code>FALSE</code>.
   */
  bool GetGrid_Movement(void);
  
  /*!
   * \brief Get the type of dynamic mesh motion.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Type of dynamic mesh motion.
   */
  unsigned short GetKind_GridMovement(unsigned short val_iZone);
  
  /*!
   * \brief Get the type of dynamic mesh motion. Each zone gets a config file.
   * \return Type of dynamic mesh motion.
   */
  unsigned short GetKind_GridMovement();

  /*!
   * \brief Set the type of dynamic mesh motion.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \param[in] motion_Type - Specify motion type.
   */
  void SetKind_GridMovement(unsigned short val_iZone, unsigned short motion_Type);
  
  /*!
   * \brief Get the mach number based on the mesh velocity and freestream quantities.
   * \return Mach number based on the mesh velocity and freestream quantities.
   */
  su2double GetMach_Motion(void);
  
  /*!
   * \brief Get x-coordinate of the mesh motion origin.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return X-coordinate of the mesh motion origin.
   */
  su2double GetMotion_Origin_X(unsigned short val_iZone);
  
  /*!
   * \brief Get y-coordinate of the mesh motion origin
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Y-coordinate of the mesh motion origin.
   */
  su2double GetMotion_Origin_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get z-coordinate of the mesh motion origin
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Z-coordinate of the mesh motion origin.
   */
  su2double GetMotion_Origin_Z(unsigned short val_iZone);
  
  /*!
   * \brief Set x-coordinate of the mesh motion origin.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \param[in] val_origin - New x-coordinate of the mesh motion origin.
   */
  void SetMotion_Origin_X(unsigned short val_iZone, su2double val_origin);
  
  /*!
   * \brief Set y-coordinate of the mesh motion origin
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \param[in] val_origin - New y-coordinate of the mesh motion origin.
   */
  void SetMotion_Origin_Y(unsigned short val_iZone, su2double val_origin);
  
  /*!
   * \brief Set z-coordinate of the mesh motion origin
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \param[in] val_origin - New y-coordinate of the mesh motion origin.
   */
  void SetMotion_Origin_Z(unsigned short val_iZone, su2double val_origin);
  
  /*!
   * \brief Get the translational velocity of the mesh in the x-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Translational velocity of the mesh in the x-direction.
   */
  su2double GetTranslation_Rate_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the translational velocity of the mesh in the y-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Translational velocity of the mesh in the y-direction.
   */
  su2double GetTranslation_Rate_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the translational velocity of the mesh in the z-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Translational velocity of the mesh in the z-direction.
   */
  su2double GetTranslation_Rate_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular velocity of the mesh about the x-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular velocity of the mesh about the x-axis.
   */
  su2double GetRotation_Rate_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular velocity of the mesh about the y-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular velocity of the mesh about the y-axis.
   */
  su2double GetRotation_Rate_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular velocity of the mesh about the z-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular velocity of the mesh about the z-axis.
   */
  su2double GetRotation_Rate_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular velocity of the mesh about the z-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular velocity of the mesh about the z-axis.
   */
  su2double GetFinalRotation_Rate_Z(unsigned short val_iZone);

  /*!
   * \brief Set the angular velocity of the mesh about the z-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \param[in] newRotation_Rate_Z - new rotation rate after computing the ramp value.
   */
  void SetRotation_Rate_Z(su2double newRotation_Rate_Z, unsigned short val_iZone);

  /*!
   * \brief Get the angular frequency of a mesh pitching about the x-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular frequency of a mesh pitching about the x-axis.
   */
  su2double GetPitching_Omega_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular frequency of a mesh pitching about the y-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular frequency of a mesh pitching about the y-axis.
   */
  su2double GetPitching_Omega_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular frequency of a mesh pitching about the z-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular frequency of a mesh pitching about the z-axis.
   */
  su2double GetPitching_Omega_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the pitching amplitude about the x-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Pitching amplitude about the x-axis.
   */
  su2double GetPitching_Ampl_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the pitching amplitude about the y-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Pitching amplitude about the y-axis.
   */
  su2double GetPitching_Ampl_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the pitching amplitude about the z-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Pitching amplitude about the z-axis.
   */
  su2double GetPitching_Ampl_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the pitching phase offset about the x-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Pitching phase offset about the x-axis.
   */
  su2double GetPitching_Phase_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the pitching phase offset about the y-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Pitching phase offset about the y-axis.
   */
  su2double GetPitching_Phase_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the pitching phase offset about the z-axis.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Pitching phase offset about the z-axis.
   */
  su2double GetPitching_Phase_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular frequency of a mesh plunging in the x-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular frequency of a mesh plunging in the x-direction.
   */
  su2double GetPlunging_Omega_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular frequency of a mesh plunging in the y-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular frequency of a mesh plunging in the y-direction.
   */
  su2double GetPlunging_Omega_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the angular frequency of a mesh plunging in the z-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Angular frequency of a mesh plunging in the z-direction.
   */
  su2double GetPlunging_Omega_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the plunging amplitude in the x-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Plunging amplitude in the x-direction.
   */
  su2double GetPlunging_Ampl_X(unsigned short val_iZone);
  
  /*!
   * \brief Get the plunging amplitude in the y-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Plunging amplitude in the y-direction.
   */
  su2double GetPlunging_Ampl_Y(unsigned short val_iZone);
  
  /*!
   * \brief Get the plunging amplitude in the z-direction.
   * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
   * \return Plunging amplitude in the z-direction.
   */
  su2double GetPlunging_Ampl_Z(unsigned short val_iZone);
  
  /*!
   * \brief Get the Harmonic Balance frequency pointer.
   * \return Harmonic Balance Frequency pointer.
   */
  su2double* GetOmega_HB(void);
	
  /*!
   * \brief Get if harmonic balance source term is to be preconditioned
   * \return yes or no to harmonic balance preconditioning
   */
  bool GetHB_Precondition(void);
  
  /*!
   * \brief Get if we should update the motion origin.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return yes or no to update motion origin.
   */
  unsigned short GetMoveMotion_Origin(unsigned short val_marker);
  
  /*!
   * \brief Get the minimum value of Beta for Roe-Turkel preconditioner
   * \return the minimum value of Beta for Roe-Turkel preconditioner
   */
  su2double GetminTurkelBeta();
  
  /*!
   * \brief Get the minimum value of Beta for Roe-Turkel preconditioner
   * \return the minimum value of Beta for Roe-Turkel preconditioner
   */
  su2double GetmaxTurkelBeta();
  
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
  bool Low_Mach_Preconditioning(void);
  
  /*!
   * \brief Get information about the Low Mach Correction
   * \return <code>TRUE</code> if we are using low Mach correction; otherwise <code>FALSE</code>.
   */
  bool Low_Mach_Correction(void);
  
  /*!
   * \brief Get information about the poisson solver condition
   * \return <code>TRUE</code> if it is a poisson solver condition; otherwise <code>FALSE</code>.
   */
  bool GetPoissonSolver(void);
  
  /*!
   * \brief Get information about the gravity force.
   * \return <code>TRUE</code> if it uses the gravity force; otherwise <code>FALSE</code>.
   */
  bool GetGravityForce(void);

  /*!
   * \brief Get information about the body force.
   * \return <code>TRUE</code> if it uses a body force; otherwise <code>FALSE</code>.
   */
  bool GetBody_Force(void);

  /*!
   * \brief Get a pointer to the body force vector.
   * \return A pointer to the body force vector.
   */
  su2double* GetBody_Force_Vector(void);

  /*!
   * \brief Get information about the rotational frame.
   * \return <code>TRUE</code> if there is a rotational frame; otherwise <code>FALSE</code>.
   */
  bool GetRotating_Frame(void);
  
  /*!
   * \brief Get information about the axisymmetric frame.
   * \return <code>TRUE</code> if there is a rotational frame; otherwise <code>FALSE</code>.
   */
  bool GetAxisymmetric(void);
    
  /*!
   * \brief Get information about the axisymmetric frame.
   * \return <code>TRUE</code> if there is a rotational frame; otherwise <code>FALSE</code>.
   */
  bool GetDebugMode(void);
  
  /*!
   * \brief Get information about there is a smoothing of the grid coordinates.
   * \return <code>TRUE</code> if there is smoothing of the grid coordinates; otherwise <code>FALSE</code>.
   */
  bool GetAdaptBoundary(void);
  
  /*!
   * \brief Get information about there is a smoothing of the grid coordinates.
   * \return <code>TRUE</code> if there is smoothing of the grid coordinates; otherwise <code>FALSE</code>.
   */
  bool GetSmoothNumGrid(void);
  
  /*!
   * \brief Set information about there is a smoothing of the grid coordinates.
   * \param[in] val_smoothnumgrid - <code>TRUE</code> if there is smoothing of the grid coordinates; otherwise <code>FALSE</code>.
   */
  void SetSmoothNumGrid(bool val_smoothnumgrid);
  
  /*!
   * \brief Subtract one to the index of the finest grid (full multigrid strategy).
   * \return Change the index of the finest grid.
   */
  void SubtractFinestMesh(void);
  
  /*!
   * \brief Obtain the kind of design variable.
   * \param[in] val_dv - Number of the design variable that we want to read.
   * \return Design variable identification.
   */
  unsigned short GetDesign_Variable(unsigned short val_dv);
    
  /*!
   * \brief Provides the buffet monitoring information.
   * \return Buffet monitoring information, if <code>TRUE</code> then the code will compute the buffet sensor.
   */
  bool GetBuffet_Monitoring(void);
    
  /*!
   * \brief Get the buffet sensor sharpness coefficient.
   * \return Sharpness coefficient for buffet sensor.
   */
  su2double GetBuffet_k(void);
    
  /*!
   * \brief Get the buffet sensor offset parameter.
   * \return Offset parameter for buffet sensor.
   */
  su2double GetBuffet_lambda(void);
  
  /*!
   * \brief Obtain the kind of convergence criteria to establish the convergence of the CFD code.
   * \return Kind of convergence criteria.
   */
  unsigned short GetConvCriteria(void);
  
  /*!
   * \brief Get the index in the config information of the marker <i>val_marker</i>.
   * \note When we read the config file, it stores the markers in a particular vector.
   * \return Index in the config information of the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_TagBound(string val_marker);
  
  /*!
   * \brief Get the name in the config information of the marker number <i>val_marker</i>.
   * \note When we read the config file, it stores the markers in a particular vector.
   * \return Name of the marker in the config information of the marker <i>val_marker</i>.
   */
  string GetMarker_CfgFile_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the boundary information (kind of boundary) in the config information of the marker <i>val_marker</i>.
   * \return Kind of boundary in the config information of the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_KindBC(string val_marker);
  
  /*!
   * \brief Get the monitoring information from the config definition for the marker <i>val_marker</i>.
   * \return Monitoring information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Monitoring(string val_marker);
  
  /*!
   * \brief Get the monitoring information from the config definition for the marker <i>val_marker</i>.
   * \return Monitoring information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_GeoEval(string val_marker);
  
  /*!
   * \brief Get the monitoring information from the config definition for the marker <i>val_marker</i>.
   * \return Monitoring information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Designing(string val_marker);
  
  /*!
   * \brief Get the plotting information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Plotting(string val_marker);
  
  /*!
   * \brief Get the plotting information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Analyze(string val_marker);
  
  /*!
   * \brief Get the FSI interface information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_ZoneInterface(string val_marker);
  
  /*!
   * \brief Get the TurboPerformance information from the config definition for the marker <i>val_marker</i>.
   * \return TurboPerformance information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Turbomachinery(string val_marker);

  /*!
   * \brief Get the TurboPerformance flag information from the config definition for the marker <i>val_marker</i>.
   * \return TurboPerformance flag information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_TurbomachineryFlag(string val_marker);

  /*!
   * \brief Get the MixingPlane interface information from the config definition for the marker <i>val_marker</i>.
   * \return Plotting information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_MixingPlaneInterface(string val_marker);
  
  /*!
   * \brief Get the DV information from the config definition for the marker <i>val_marker</i>.
   * \return DV information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_DV(string val_marker);
  
  /*!
   * \brief Get the motion information from the config definition for the marker <i>val_marker</i>.
   * \return Motion information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Moving(string val_marker);

  /*!
   * \brief Get the Python customization information from the config definition for the marker <i>val_marker</i>.
   * \return Python customization information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_PyCustom(string val_marker);
  
  /*!
   * \brief Get the periodic information from the config definition of the marker <i>val_marker</i>.
   * \return Periodic information of the boundary in the config information of the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_PerBound(string val_marker);
  
  /*!
   * \brief  Get the name of the marker <i>val_marker</i>.
   * \return The interface which owns that marker <i>val_marker</i>.
   */
  int GetMarker_ZoneInterface(string val_marker);
  
  /*!
   * \brief  Get the name of the marker <i>val_iMarker</i>.
   * \return The name of the marker in the interface
   */
  string GetMarkerTag_ZoneInterface(unsigned short val_iMarker);

  /*!
   * \brief  Get the number of markers in the multizone interface.
   * \return The number markers in the multizone interface
   */
  unsigned short GetnMarker_ZoneInterface(void);

  /*!
   * \brief Determines if problem is adjoint
   * \return true if Adjoint
   */
  bool GetContinuous_Adjoint(void);
  
  /*!
   * \brief Determines if problem is viscous
   * \return true if Viscous
   */
  bool GetViscous(void);
  
  /*!
   * \brief Provides the index of the solution in the container.
   * \param[in] val_eqsystem - Equation that is being solved.
   * \return Index on the solution container.
   */
  unsigned short GetContainerPosition(unsigned short val_eqsystem);
  
  /*!
   * \brief Value of the order of magnitude reduction of the residual.
   * \return Value of the order of magnitude reduction of the residual.
   */
  su2double GetOrderMagResidual(void);
  
  /*!
   * \brief Value of the minimum residual value (log10 scale).
   * \return Value of the minimum residual value (log10 scale).
   */
  su2double GetMinLogResidual(void);
  
  /*!
   * \brief Value of the order of magnitude reduction of the residual for FSI applications.
   * \return Value of the order of magnitude reduction of the residual.
   */
  su2double GetOrderMagResidualFSI(void);
  
  /*!
   * \brief Value of the minimum residual value for FSI applications (log10 scale).
   * \return Value of the minimum residual value (log10 scale).
   */
  su2double GetMinLogResidualFSI(void);
  
  /*!
   * \brief Value of the order of magnitude reduction of the flow residual for BGS applications.
   * \return Value of the order of magnitude reduction of the residual.
   */
  su2double GetOrderMagResidual_BGS_F(void);
  
  /*!
   * \brief Value of the minimum flow residual value for BGS applications (log10 scale).
   * \return Value of the minimum residual value (log10 scale).
   */
  su2double GetMinLogResidual_BGS_F(void);
  
  /*!
   * \brief Value of the order of magnitude reduction of the flow residual for BGS applications.
   * \return Value of the order of magnitude reduction of the residual.
   */
  su2double GetOrderMagResidual_BGS_S(void);
  
  /*!
   * \brief Value of the minimum flow residual value for BGS applications (log10 scale).
   * \return Value of the minimum residual value (log10 scale).
   */
  su2double GetMinLogResidual_BGS_S(void);
  
  /*!
   * \brief Value of the displacement tolerance UTOL for FEM structural analysis (log10 scale).
   * \return Value of Res_FEM_UTOL (log10 scale).
   */
  su2double GetResidual_FEM_UTOL(void);
  
  /*!
   * \brief Value of the displacement tolerance UTOL for FEM structural analysis (log10 scale).
   * \return Value of Res_FEM_UTOL (log10 scale).
   */
  su2double GetResidual_FEM_RTOL(void);
  
  /*!
   * \brief Value of the displacement tolerance UTOL for FEM structural analysis (log10 scale).
   * \return Value of Res_FEM_UTOL (log10 scale).
   */
  su2double GetResidual_FEM_ETOL(void);
  
  /*!
   * \brief Value of the maximum objective function for FEM elasticity adjoint (log10 scale).
   * \return Value of Res_FEM_ADJ (log10 scale).
   */
  su2double GetCriteria_FEM_ADJ(void);
  
  /*!
   * \brief Value of the damping factor for the engine inlet bc.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Engine_Inflow(void);
  
  /*!
   * \brief Value of the damping factor for the engine exhaust inlet bc.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Engine_Exhaust(void);
  
  /*!
   * \brief Value of the damping factor for the residual restriction.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Res_Restric(void);
  
  /*!
   * \brief Value of the damping factor for the correction prolongation.
   * \return Value of the damping factor.
   */
  su2double GetDamp_Correc_Prolong(void);
  
  /*!
   * \brief Value of the position of the Near Field (y coordinate for 2D, and z coordinate for 3D).
   * \return Value of the Near Field position.
   */
  su2double GetPosition_Plane(void);
  
  /*!
   * \brief Value of the weight of the drag coefficient in the Sonic Boom optimization.
   * \return Value of the weight of the drag coefficient in the Sonic Boom optimization.
   */
  su2double GetWeightCd(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdNetThrust_dBCThrust(su2double val_dnetthrust_dbcthrust);
  
  /*!
   * \brief Value of the azimuthal line to fix due to a misalignments of the nearfield.
   * \return Azimuthal line to fix due to a misalignments of the nearfield.
   */
  su2double GetFixAzimuthalLine(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCD_dCMy(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetCM_Target(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCD_dCL(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCD_dCL(su2double val_dcd_dcl);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCMx_dCL(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCMx_dCL(su2double val_dcmx_dcl);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCMy_dCL(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCMy_dCL(su2double val_dcmy_dcl);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetdCMz_dCL(void);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCMz_dCL(su2double val_dcmz_dcl);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCL_dAlpha(su2double val_dcl_dalpha);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCM_diH(su2double val_dcm_dhi);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  void SetdCD_dCMy(su2double val_dcd_dcmy);
  
  /*!
   * \brief Value of the weight of the CD, CL, CM optimization.
   * \return Value of the weight of the CD, CL, CM optimization.
   */
  su2double GetCL_Target(void);
  
  /*!
   * \brief Set the global parameters of each simulation for each runtime system.
   * \param[in] val_solver - Solver of the simulation.
   * \param[in] val_system - Runtime system that we are solving.
   */
  void SetGlobalParam(unsigned short val_solver, unsigned short val_system, unsigned long val_extiter);
  
  /*!
   * \brief Center of rotation for a rotational periodic boundary.
   */
  su2double *GetPeriodicRotCenter(string val_marker);
  
  /*!
   * \brief Angles of rotation for a rotational periodic boundary.
   */
  su2double *GetPeriodicRotAngles(string val_marker);
  
  /*!
   * \brief Translation vector for a rotational periodic boundary.
   */
  su2double *GetPeriodicTranslation(string val_marker);
  
  /*!
   * \brief Get the rotationally periodic donor marker for boundary <i>val_marker</i>.
   * \return Periodic donor marker from the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_Periodic_Donor(string val_marker);
  
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_NetThrust(string val_marker);
  
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_Power(string val_marker);
  
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_MassFlow(string val_marker);
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_Mach(string val_marker);
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_Force(string val_marker);
  
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_BCThrust(string val_marker);
  
  /*!
   * \brief Get the origin of the actuator disk.
   */
  su2double GetActDisk_BCThrust_Old(string val_marker);
  
  /*!
   * \brief Get the tip radius of th actuator disk.
   */
  su2double GetActDisk_Area(string val_marker);
  
  /*!
   * \brief Get the tip radius of th actuator disk.
   */
  su2double GetActDisk_ReverseMassFlow(string val_marker);
  
  /*!
   * \brief Get the thrust corffient of the actuator disk.
   */
  su2double GetActDisk_PressJump(string val_marker, unsigned short val_index);
  
  /*!
   * \brief Get the thrust corffient of the actuator disk.
   */
  su2double GetActDisk_TempJump(string val_marker, unsigned short val_index);
  
  /*!
   * \brief Get the rev / min of the actuator disk.
   */
  su2double GetActDisk_Omega(string val_marker, unsigned short val_index);
  
  /*!
   * \brief Get Actuator Disk Outlet for boundary <i>val_marker</i> (actuator disk inlet).
   * \return Actuator Disk Outlet from the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_ActDiskOutlet(string val_marker);
  
  /*!
   * \brief Get Actuator Disk Outlet for boundary <i>val_marker</i> (actuator disk inlet).
   * \return Actuator Disk Outlet from the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_EngineExhaust(string val_marker);
  
  /*!
   * \brief Get the internal index for a moving boundary <i>val_marker</i>.
   * \return Internal index for a moving boundary <i>val_marker</i>.
   */
  unsigned short GetMarker_Moving(string val_marker);
  
  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Moving_TagBound(unsigned short val_marker);

  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_PyCustom_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the name of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Name that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_Analyze_TagBound(unsigned short val_marker);
  
  /*!
   * \brief Get the total temperature at a nacelle boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total temperature.
   */
  su2double GetExhaust_Temperature_Target(string val_index);
  
  /*!
   * \brief Get the total temperature at an inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total temperature.
   */
  su2double GetInlet_Ttotal(string val_index);
  
  /*!
   * \brief Get the temperature at a supersonic inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet density.
   */
  su2double GetInlet_Temperature(string val_index);
  
  /*!
   * \brief Get the pressure at a supersonic inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet pressure.
   */
  su2double GetInlet_Pressure(string val_index);
  
  /*!
   * \brief Get the velocity vector at a supersonic inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The inlet velocity vector.
   */
  su2double* GetInlet_Velocity(string val_index);
  
  /*!
   * \brief Get the fixed value at the Dirichlet boundary.
   * \param[in] val_index - Index corresponding to the Dirichlet boundary.
   * \return The total temperature.
   */
  su2double GetDirichlet_Value(string val_index);
  
  /*!
   * \brief Get whether this is a Dirichlet or a Neumann boundary.
   * \param[in] val_index - Index corresponding to the Dirichlet boundary.
   * \return Yes or No.
   */
  bool GetDirichlet_Boundary(string val_index);
  
  /*!
   * \brief Get the total pressure at an inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total pressure.
   */
  su2double GetInlet_Ptotal(string val_index);

  /*!
   * \brief Set the total pressure at an inlet boundary.
   * \param[in] val_pressure - Pressure value at the inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   */
  void SetInlet_Ptotal(su2double val_pressure, string val_marker);

  /*!
   * \brief Get the total pressure at an nacelle boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The total pressure.
   */
  su2double GetExhaust_Pressure_Target(string val_index);
  
  /*!
   * \brief Value of the CFL reduction in LevelSet problems.
   * \return Value of the CFL reduction in LevelSet problems.
   */
  su2double GetCFLRedCoeff_Turb(void);
  
  /*!
   * \brief Get the flow direction unit vector at an inlet boundary.
   * \param[in] val_index - Index corresponding to the inlet boundary.
   * \return The flow direction vector.
   */
  su2double* GetInlet_FlowDir(string val_index);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_Pressure(string val_index);

  /*!
   * \brief Set the back pressure (static) at an outlet boundary.
   * \param[in] val_pressure - Pressure value at the outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   */
  void SetOutlet_Pressure(su2double val_pressure, string val_marker);

  /*!
   * \brief Get the var 1 at Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return The var1
   */
  su2double GetRiemann_Var1(string val_marker);
  
  /*!
   * \brief Get the var 2 at Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return The var2
   */
  su2double GetRiemann_Var2(string val_marker);
  
  /*!
   * \brief Get the Flowdir at Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return The Flowdir
   */
  su2double* GetRiemann_FlowDir(string val_marker);
  
  /*!
   * \brief Get Kind Data of Riemann boundary.
   * \param[in] val_marker - Index corresponding to the Riemann boundary.
   * \return Kind data
   */
  unsigned short GetKind_Data_Riemann(string val_marker);
  
  /*!
   * \brief Get the var 1 for the Giels BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The var1
   */
  su2double GetGiles_Var1(string val_marker);
  
  /*!
   * \brief Get the var 2 for the Giles boundary.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The var2
   */
  su2double GetGiles_Var2(string val_marker);
  
  /*!
   * \brief Get the Flowdir for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The Flowdir
   */
  su2double* GetGiles_FlowDir(string val_marker);
  
  /*!
   * \brief Get Kind Data for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return Kind data
   */
  unsigned short GetKind_Data_Giles(string val_marker);
  
  /*!
   * \brief Set the var 1 for Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   */
  void SetGiles_Var1(su2double newVar1, string val_marker);

  /*!
   * \brief Get the relax factor for the average component for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The relax factor for the average component
   */
  su2double GetGiles_RelaxFactorAverage(string val_marker);

  /*!
   * \brief Get the relax factor for the fourier component for the Giles BC.
   * \param[in] val_marker - Index corresponding to the Giles BC.
   * \return The relax factor for the fourier component
   */
  su2double GetGiles_RelaxFactorFourier(string val_marker);

  /*!
   * \brief Get the outlet pressure imposed as BC for internal flow.
   * \return outlet pressure
   */
  su2double GetPressureOut_BC();

  /*!
   * \brief Set the outlet pressure imposed as BC for internal flow.
   * \param[in] val_temp - New value of the outlet pressure.
   */
  void SetPressureOut_BC(su2double val_press);

  /*!
   * \brief Get the inlet velocity or pressure imposed for incompressible flow.
   * \return inlet velocity or pressure
   */
  su2double GetIncInlet_BC();

  /*!
   * \brief Set the inlet velocity or pressure imposed as BC for incompressible flow.
   * \param[in] val_in - New value of the inlet velocity or pressure.
   */
  void SetIncInlet_BC(su2double val_in);

  /*!
   * \brief Get the inlet temperature imposed as BC for incompressible flow.
   * \return inlet temperature
   */
  su2double GetIncTemperature_BC();

  /*!
   * \brief Set the inlet temperature imposed as BC for incompressible flow.
   * \param[in] val_temperature - New value of the inlet temperature.
   */
  void SetIncTemperature_BC(su2double val_temperature);

  /*!
   * \brief Get the outlet pressure imposed as BC for incompressible flow.
   * \return outlet pressure
   */
  su2double GetIncPressureOut_BC();

  /*!
   * \brief Set the outlet pressure imposed as BC for incompressible flow.
   * \param[in] val_pressure - New value of the outlet pressure.
   */
  void SetIncPressureOut_BC(su2double val_pressure);

  /*!
   * \brief Get the inlet total pressure imposed as BC for internal flow.
   * \return inlet total pressure
   */
  su2double GetTotalPressureIn_BC();

  /*!
   * \brief Get the inlet total temperature imposed as BC for internal flow.
   * \return inlet total temperature
   */
  su2double GetTotalTemperatureIn_BC();

  /*!
   * \brief Set the inlet total temperature imposed as BC for internal flow.
   * \param[in] val_temp - New value of the total temperature.
   */
  void SetTotalTemperatureIn_BC(su2double val_temp);

  /*!
   * \brief Get the inlet flow angle imposed as BC for internal flow.
   * \return inlet flow angle
   */
  su2double GetFlowAngleIn_BC();

  /*!
   * \brief Get the wall temperature (static) at an isothermal boundary.
   * \param[in] val_index - Index corresponding to the isothermal boundary.
   * \return The wall temperature.
   */
  su2double GetIsothermal_Temperature(string val_index);
  
  /*!
   * \brief Get the wall heat flux on a constant heat flux boundary.
   * \param[in] val_index - Index corresponding to the constant heat flux boundary.
   * \return The heat flux.
   */
  su2double GetWall_HeatFlux(string val_index);

  /*!
   * \brief Get the wall function treatment for the given boundary marker.
   * \param[in] val_marker - String of the viscous wall marker.
   * \return The type of wall function treatment.
   */
  unsigned short GetWallFunction_Treatment(string val_marker);

  /*!
   * \brief Get the additional integer info for the wall function treatment
            for the given boundary marker.
   * \param[in] val_marker - String of the viscous wall marker.
   * \return Pointer to the integer info for the given marker.
   */
  unsigned short* GetWallFunction_IntInfo(string val_marker);

  /*!
   * \brief Get the additional double info for the wall function treatment
            for the given boundary marker.
   * \param[in] val_marker - String of the viscous wall marker.
   * \return Pointer to the double info for the given marker.
   */
  su2double* GetWallFunction_DoubleInfo(string val_marker);
  
  /*!
   * \brief Get the target (pressure, massflow, etc) at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \return Target (pressure, massflow, etc) .
   */
  su2double GetEngineInflow_Target(string val_marker);
  
  /*!
   * \brief Get the fan face Mach number at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The fan face Mach number.
   */
  su2double GetInflow_Mach(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow pressure.
   */
  su2double GetInflow_Pressure(string val_marker);
  
  /*!
   * \brief Get the mass flow rate at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine mass flow rate.
   */
  su2double GetInflow_MassFlow(string val_marker);
  
  /*!
   * \brief Get the percentage of reverse flow at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The percentage of reverse flow.
   */
  su2double GetInflow_ReverseMassFlow(string val_marker);
  
  /*!
   * \brief Get the percentage of reverse flow at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \return The percentage of reverse flow.
   */
  su2double GetInflow_ReverseMassFlow(unsigned short val_marker);
  
  /*!
   * \brief Get the total pressure at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The total pressure.
   */
  su2double GetInflow_TotalPressure(string val_marker);
  
  /*!
   * \brief Get the temperature (static) at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow temperature.
   */
  su2double GetInflow_Temperature(string val_marker);
  
  /*!
   * \brief Get the total temperature at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow total temperature.
   */
  su2double GetInflow_TotalTemperature(string val_marker);
  
  /*!
   * \brief Get the ram drag at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow ram drag.
   */
  su2double GetInflow_RamDrag(string val_marker);
  
  /*!
   * \brief Get the force balance at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow force balance.
   */
  su2double GetInflow_Force(string val_marker);
  
  /*!
   * \brief Get the power at an engine inflow boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine inflow power.
   */
  su2double GetInflow_Power(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust pressure.
   */
  su2double GetExhaust_Pressure(string val_marker);
  
  /*!
   * \brief Get the temperature (static) at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust temperature.
   */
  su2double GetExhaust_Temperature(string val_marker);
  
  /*!
   * \brief Get the massflow at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust massflow.
   */
  su2double GetExhaust_MassFlow(string val_marker);
  
  /*!
   * \brief Get the total pressure at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The engine exhaust total pressure.
   */
  su2double GetExhaust_TotalPressure(string val_marker);
  
  /*!
   * \brief Get the total temperature at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return The total temperature.
   */
  su2double GetExhaust_TotalTemperature(string val_marker);
  
  /*!
   * \brief Get the gross thrust at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return Gross thrust.
   */
  su2double GetExhaust_GrossThrust(string val_marker);
  
  /*!
   * \brief Get the force balance at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return Force balance.
   */
  su2double GetExhaust_Force(string val_marker);
  
  /*!
   * \brief Get the power at an engine exhaust boundary.
   * \param[in] val_marker - Name of the boundary.
   * \return Power.
   */
  su2double GetExhaust_Power(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetInflow_Mach(unsigned short val_imarker, su2double val_fanface_mach);
  
  /*!
   * \brief Set the fan face static pressure at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_pressure - Fan face static pressure.
   */
  void SetInflow_Pressure(unsigned short val_imarker, su2double val_fanface_pressure);
  
  /*!
   * \brief Set the massflow at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_massflow - Massflow.
   */
  void SetInflow_MassFlow(unsigned short val_imarker, su2double val_fanface_massflow);
  
  /*!
   * \brief Set the reverse flow at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_reversemassflow - reverse flow.
   */
  void SetInflow_ReverseMassFlow(unsigned short val_imarker, su2double val_fanface_reversemassflow);
  
  /*!
   * \brief Set the fan face total pressure at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_totalpressure - Fan face total pressure.
   */
  void SetInflow_TotalPressure(unsigned short val_imarker, su2double val_fanface_totalpressure);
  
  /*!
   * \brief Set the fan face static temperature at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_pressure - Fan face static temperature.
   */
  void SetInflow_Temperature(unsigned short val_imarker, su2double val_fanface_temperature);
  
  /*!
   * \brief Set the fan face total temperature at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_totaltemperature - Fan face total temperature.
   */
  void SetInflow_TotalTemperature(unsigned short val_imarker, su2double val_fanface_totaltemperature);
  
  /*!
   * \brief Set the ram drag temperature at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_ramdrag - Ram drag value.
   */
  void SetInflow_RamDrag(unsigned short val_imarker, su2double val_fanface_ramdrag);
  
  /*!
   * \brief Set the force balance at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_force - Fan face force.
   */
  void SetInflow_Force(unsigned short val_imarker, su2double val_fanface_force);
  
  /*!
   * \brief Set the power at an engine inflow boundary.
   * \param[in] val_index - Index corresponding to the engine inflow boundary.
   * \param[in] val_fanface_force - Power.
   */
  void SetInflow_Power(unsigned short val_imarker, su2double val_fanface_power);
  
  /*!
   * \brief Set the back pressure (static) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_pressure - Exhaust static pressure.
   */
  void SetExhaust_Pressure(unsigned short val_imarker, su2double val_exhaust_pressure);
  
  /*!
   * \brief Set the temperature (static) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_temp - Exhaust static temperature.
   */
  void SetExhaust_Temperature(unsigned short val_imarker, su2double val_exhaust_temp);
  
  /*!
   * \brief Set the back pressure (static) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_temp - Exhaust static temperature.
   */
  void SetExhaust_MassFlow(unsigned short val_imarker, su2double val_exhaust_massflow);
  
  /*!
   * \brief Set the back pressure (total) at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_totalpressure - Exhaust total pressure.
   */
  void SetExhaust_TotalPressure(unsigned short val_imarker, su2double val_exhaust_totalpressure);
  
  /*!
   * \brief Set the total temperature at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_totaltemp - Exhaust total temperature.
   */
  void SetExhaust_TotalTemperature(unsigned short val_imarker, su2double val_exhaust_totaltemp);
  
  /*!
   * \brief Set the gross thrust at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_grossthrust - Exhaust gross thrust temperature.
   */
  void SetExhaust_GrossThrust(unsigned short val_imarker, su2double val_exhaust_grossthrust);
  
  /*!
   * \brief Set the force balance at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_force - Exhaust force balance.
   */
  void SetExhaust_Force(unsigned short val_imarker, su2double val_exhaust_force);
  
  /*!
   * \brief Set the power at an engine exhaust boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \param[in] val_exhaust_power - Exhaust power.
   */
  void SetExhaust_Power(unsigned short val_imarker, su2double val_exhaust_power);
  
  /*!
   * \brief Set the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_mach - Exhaust power.
   */
  void SetEngine_Mach(unsigned short val_imarker, su2double val_engine_mach);
  
  /*!
   * \brief Set the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_force - Exhaust power.
   */
  void SetEngine_Force(unsigned short val_imarker, su2double val_engine_force);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_power - Exhaust power.
   */
  void SetEngine_Power(unsigned short val_imarker, su2double val_engine_power);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_netthrust - Exhaust power.
   */
  void SetEngine_NetThrust(unsigned short val_imarker, su2double val_engine_netthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_grossthrust - Exhaust power.
   */
  void SetEngine_GrossThrust(unsigned short val_imarker, su2double val_engine_grossthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \param[in] val_engine_area - Exhaust power.
   */
  void SetEngine_Area(unsigned short val_imarker, su2double val_engine_area);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Mach(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Force(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Power(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  
  su2double GetEngine_NetThrust(unsigned short val_imarker);
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  
  su2double GetEngine_GrossThrust(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_imarker - Index corresponding to a particular engine boundary.
   * \return The outlet pressure.
   */
  su2double GetEngine_Area(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Temperature(unsigned short val_imarker, su2double val_actdisk_temp);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_TotalTemperature(unsigned short val_imarker, su2double val_actdisk_totaltemp);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Temperature(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_TotalTemperature(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Temperature(unsigned short val_imarker, su2double val_actdisk_temp);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_TotalTemperature(unsigned short val_imarker, su2double val_actdisk_totaltemp);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Temperature(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_TotalTemperature(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_MassFlow(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_MassFlow(unsigned short val_imarker, su2double val_actdisk_massflow);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_MassFlow(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_MassFlow(unsigned short val_imarker, su2double val_actdisk_massflow);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Pressure(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_TotalPressure(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_DeltaPress(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_DeltaTemp(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_TotalPressRatio(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_TotalTempRatio(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_StaticPressRatio(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_StaticTempRatio(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_NetThrust(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_BCThrust(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_BCThrust_Old(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_GrossThrust(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Area(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_ReverseMassFlow(unsigned short val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_RamDrag(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Force(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskInlet_Power(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Pressure(unsigned short val_imarker, su2double val_actdisk_pressure);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_TotalPressure(unsigned short val_imarker, su2double val_actdisk_totalpressure);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_DeltaPress(unsigned short val_imarker, su2double val_actdisk_deltapress);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Power(unsigned short val_imarker, su2double val_actdisk_power);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_MassFlow(unsigned short val_imarker, su2double val_actdisk_massflow);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Mach(unsigned short val_imarker, su2double val_actdisk_mach);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Force(unsigned short val_imarker, su2double val_actdisk_force);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_MassFlow(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetOutlet_MassFlow(unsigned short val_imarker, su2double val_massflow);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_Density(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetOutlet_Density(unsigned short val_imarker, su2double val_density);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetOutlet_Area(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetOutlet_Area(unsigned short val_imarker, su2double val_area);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_DC60(unsigned short val_imarker, su2double val_surface_distortion);
  
  /*!
   * \brief Set the massflow at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the mass flow.
   */
  void SetSurface_MassFlow(unsigned short val_imarker, su2double val_surface_massflow);
  
  /*!
   * \brief Set the mach number at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the mach number.
   */
  void SetSurface_Mach(unsigned short val_imarker, su2double val_surface_mach);

  /*!
   * \brief Set the temperature at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the temperature.
   */
  void SetSurface_Temperature(unsigned short val_imarker, su2double val_surface_temperature);

  /*!
   * \brief Set the pressure at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_massflow - Value of the pressure.
   */
  void SetSurface_Pressure(unsigned short val_imarker, su2double val_surface_pressure);

  /*!
   * \brief Set the density at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_density - Value of the density.
   */
  void SetSurface_Density(unsigned short val_imarker, su2double val_surface_density);

  /*!
   * \brief Set the enthalpy at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_density - Value of the density.
   */
  void SetSurface_Enthalpy(unsigned short val_imarker, su2double val_surface_enthalpy);

  /*!
   * \brief Set the normal velocity at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_normalvelocity - Value of the normal velocity.
   */
  void SetSurface_NormalVelocity(unsigned short val_imarker, su2double val_surface_normalvelocity);

  /*!
   * \brief Set the streamwise flow uniformity at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_streamwiseuniformity - Value of the streamwise flow uniformity.
   */
  void SetSurface_Uniformity(unsigned short val_imarker, su2double val_surface_streamwiseuniformity);

  /*!
   * \brief Set the secondary flow strength at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_secondarystrength - Value of the secondary flow strength.
   */
  void SetSurface_SecondaryStrength(unsigned short val_imarker, su2double val_surface_secondarystrength);

  /*!
   * \brief Set the relative secondary flow strength at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_secondaryoverstream - Value of the relative seondary flow strength.
   */
  void SetSurface_SecondOverUniform(unsigned short val_imarker, su2double val_surface_secondaryoverstream);

  /*!
   * \brief Set the momentum distortion at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_momentumdistortion - Value of the momentum distortion.
   */
  void SetSurface_MomentumDistortion(unsigned short val_imarker, su2double val_surface_momentumdistortion);

  /*!
   * \brief Set the total temperature at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_totaltemperature - Value of the total temperature.
   */
  void SetSurface_TotalTemperature(unsigned short val_imarker, su2double val_surface_totaltemperature);

  /*!
   * \brief Set the total pressure at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_totalpressure - Value of the total pressure.
   */
  void SetSurface_TotalPressure(unsigned short val_imarker, su2double val_surface_totalpressure);

  /*!
   * \brief Set the pressure drop between two surfaces.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \param[in] val_surface_pressuredrop - Value of the pressure drop.
   */
  void SetSurface_PressureDrop(unsigned short val_imarker, su2double val_surface_pressuredrop);

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_IDC(unsigned short val_imarker, su2double val_surface_distortion);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_IDC_Mach(unsigned short val_imarker, su2double val_surface_distortion);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetSurface_IDR(unsigned short val_imarker, su2double val_surface_distortion);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_DeltaTemp(unsigned short val_imarker, su2double val_actdisk_deltatemp);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_TotalPressRatio(unsigned short val_imarker, su2double val_actdisk_pressratio);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_TotalTempRatio(unsigned short val_imarker, su2double val_actdisk_tempratio);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_StaticPressRatio(unsigned short val_imarker, su2double val_actdisk_pressratio);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_StaticTempRatio(unsigned short val_imarker, su2double val_actdisk_tempratio);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_NetThrust(unsigned short val_imarker, su2double val_actdisk_netthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust(string val_marker, su2double val_actdisk_bcthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust(unsigned short val_imarker, su2double val_actdisk_bcthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust_Old(string val_marker, su2double val_actdisk_bcthrust_old);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_BCThrust_Old(unsigned short val_imarker, su2double val_actdisk_bcthrust_old);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_GrossThrust(unsigned short val_imarker, su2double val_actdisk_grossthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDisk_Area(unsigned short val_imarker, su2double val_actdisk_area);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_ReverseMassFlow(unsigned short val_imarker, su2double val_actdisk_area);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_RamDrag(unsigned short val_imarker, su2double val_actdisk_ramdrag);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Force(unsigned short val_imarker, su2double val_actdisk_force);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskInlet_Power(unsigned short val_imarker, su2double val_actdisk_power);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Power(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_MassFlow(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Mach(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDisk_Force(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_DC60(unsigned short val_imarker);
  
  /*!
   * \brief Get the massflow at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The massflow.
   */
  su2double GetSurface_MassFlow(unsigned short val_imarker);
  
  /*!
   * \brief Get the mach number at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The mach number.
   */
  su2double GetSurface_Mach(unsigned short val_imarker);

  /*!
   * \brief Get the temperature at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The temperature.
   */
  su2double GetSurface_Temperature(unsigned short val_imarker);

  /*!
   * \brief Get the pressure at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The pressure.
   */
  su2double GetSurface_Pressure(unsigned short val_imarker);

  /*!
   * \brief Get the density at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The density.
   */
  su2double GetSurface_Density(unsigned short val_imarker);

  /*!
   * \brief Get the enthalpy at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The density.
   */
  su2double GetSurface_Enthalpy(unsigned short val_imarker);

  /*!
   * \brief Get the normal velocity at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The normal velocity.
   */
  su2double GetSurface_NormalVelocity(unsigned short val_imarker);

  /*!
   * \brief Get the streamwise flow uniformity at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \return The streamwise flow uniformity.
   */
  su2double GetSurface_Uniformity(unsigned short val_imarker);

  /*!
   * \brief Get the secondary flow strength at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \return The secondary flow strength.
   */
  su2double GetSurface_SecondaryStrength(unsigned short val_imarker);

  /*!
   * \brief Get the relative secondary flow strength at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \return The relative seondary flow strength.
   */
  su2double GetSurface_SecondOverUniform(unsigned short val_imarker);

  /*!
   * \brief Get the momentum distortion at the surface.
   * \param[in] val_imarker - Index corresponding to the outlet boundary.
   * \return The momentum distortion.
   */
  su2double GetSurface_MomentumDistortion(unsigned short val_imarker);

  /*!
   * \brief Get the total temperature at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The total temperature.
   */
  su2double GetSurface_TotalTemperature(unsigned short val_imarker);

  /*!
   * \brief Get the total pressure at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The total pressure.
   */
  su2double GetSurface_TotalPressure(unsigned short val_imarker);

  /*!
   * \brief Get the pressure drop between two surfaces.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The pressure drop.
   */
  su2double GetSurface_PressureDrop(unsigned short val_imarker);

  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_IDC(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_IDC_Mach(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetSurface_IDR(unsigned short val_imarker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Pressure(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_TotalPressure(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_GrossThrust(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Force(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  su2double GetActDiskOutlet_Power(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Pressure(unsigned short val_imarker, su2double val_actdisk_pressure);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_TotalPressure(unsigned short val_imarker, su2double val_actdisk_totalpressure);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_GrossThrust(unsigned short val_imarker, su2double val_actdisk_grossthrust);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Force(unsigned short val_imarker, su2double val_actdisk_force);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetActDiskOutlet_Power(unsigned short val_imarker, su2double val_actdisk_power);
  
  /*!
   * \brief Get the displacement value at an displacement boundary.
   * \param[in] val_index - Index corresponding to the displacement boundary.
   * \return The displacement value.
   */
  su2double GetDispl_Value(string val_index);
  
  /*!
   * \brief Get the force value at an load boundary.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetLoad_Value(string val_index);
  
  /*!
   * \brief Get the constant value at a damper boundary.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The damper constant.
   */
  su2double GetDamper_Constant(string val_index);
  
  /*!
   * \brief Get the force value at a load boundary defined in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetLoad_Dir_Value(string val_index);
  
  /*!
   * \brief Get the force multiplier at a load boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load multiplier.
   */
  su2double GetLoad_Dir_Multiplier(string val_index);
  
  /*!
   * \brief Get the force value at a load boundary defined in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetDisp_Dir_Value(string val_index);
  
  /*!
   * \brief Get the force multiplier at a load boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load multiplier.
   */
  su2double GetDisp_Dir_Multiplier(string val_index);
  
  /*!
   * \brief Get the force direction at a loaded boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load direction.
   */
  su2double* GetLoad_Dir(string val_index);
  
  /*!
   * \brief Get the force direction at a loaded boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load direction.
   */
  su2double* GetDisp_Dir(string val_index);
  
  /*!
   * \brief Get the amplitude of the sine-wave at a load boundary defined in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetLoad_Sine_Amplitude(string val_index);
  
  /*!
   * \brief Get the frequency of the sine-wave at a load boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load frequency.
   */
  su2double GetLoad_Sine_Frequency(string val_index);
  
  /*!
   * \brief Get the force direction at a sine-wave loaded boundary in cartesian coordinates.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load direction.
   */
  su2double* GetLoad_Sine_Dir(string val_index);
  
  /*!
   * \brief Get the force value at an load boundary.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The load value.
   */
  su2double GetFlowLoad_Value(string val_index);
  
  /*!
   * \brief Cyclic pitch amplitude for rotor blades.
   * \return The specified cyclic pitch amplitude.
   */
  su2double GetCyclic_Pitch(void);
  
  /*!
   * \brief Collective pitch setting for rotor blades.
   * \return The specified collective pitch setting.
   */
  su2double GetCollective_Pitch(void);
  
  /*!
   * \brief Get name of the arbitrary mesh motion input file.
   * \return File name of the arbitrary mesh motion input file.
   */
  string GetDV_Filename(void);
  
  /*!
   * \brief Get name of the unordered ASCII volume sensitivity file.
   * \return File name of the unordered ASCII volume sensitivity file.
   */
  string GetDV_Unordered_Sens_Filename(void);
  
  /*!
   * \brief Get name of the unordered ASCII surface sensitivity file.
   * \return File name of the unordered ASCII surface sensitivity file.
   */
  string GetDV_Sens_Filename(void);
  
  /*!
   * \brief Set the config options.
   */
  void SetConfig_Options(unsigned short val_iZone, unsigned short val_nZone);
  
  /*!
   * \brief Set the config options.
   */
  void SetRunTime_Options(void);
  
  /*!
   * \brief Set the config file parsing.
   */
  void SetConfig_Parsing(char case_filename[MAX_STRING_SIZE]);
  
  /*!
   * \brief Set the config file parsing.
   */
  bool SetRunTime_Parsing(char case_filename[MAX_STRING_SIZE]);
  
  /*!
   * \brief Config file postprocessing.
   */
  void SetPostprocessing(unsigned short val_software, unsigned short val_izone, unsigned short val_nDim);
  
  /*!
   * \brief Config file markers processing.
   */
  void SetMarkers(unsigned short val_software);
  
  /*!
   * \brief Config file output.
   */
  void SetOutput(unsigned short val_software, unsigned short val_izone);
  
  /*!
   * \brief Value of Aeroelastic solution coordinate at time n+1.
   */
  vector<vector<su2double> > GetAeroelastic_np1(unsigned short iMarker);
  
  /*!
   * \brief Value of Aeroelastic solution coordinate at time n.
   */
  vector<vector<su2double> > GetAeroelastic_n(unsigned short iMarker);
  
  /*!
   * \brief Value of Aeroelastic solution coordinate at time n-1.
   */
  vector<vector<su2double> > GetAeroelastic_n1(unsigned short iMarker);
  
  /*!
   * \brief Value of Aeroelastic solution coordinate at time n+1.
   */
  void SetAeroelastic_np1(unsigned short iMarker, vector<vector<su2double> > solution);
  
  /*!
   * \brief Value of Aeroelastic solution coordinate at time n from time n+1.
   */
  void SetAeroelastic_n(void);
  
  /*!
   * \brief Value of Aeroelastic solution coordinate at time n-1 from time n.
   */
  void SetAeroelastic_n1(void);
  
  /*!
   * \brief Aeroelastic Flutter Speed Index.
   */
  su2double GetAeroelastic_Flutter_Speed_Index(void);
  
  /*!
   * \brief Uncoupled Aeroelastic Frequency Plunge.
   */
  su2double GetAeroelastic_Frequency_Plunge(void);
  
  /*!
   * \brief Uncoupled Aeroelastic Frequency Pitch.
   */
  su2double GetAeroelastic_Frequency_Pitch(void);
  
  /*!
   * \brief Aeroelastic Airfoil Mass Ratio.
   */
  su2double GetAeroelastic_Airfoil_Mass_Ratio(void);
  
  /*!
   * \brief Aeroelastic center of gravity location.
   */
  su2double GetAeroelastic_CG_Location(void);
  
  /*!
   * \brief Aeroelastic radius of gyration squared.
   */
  su2double GetAeroelastic_Radius_Gyration_Squared(void);
  
  /*!
   * \brief Aeroelastic solve every x inner iteration.
   */
  unsigned short GetAeroelasticIter(void);
  
  /*!
   * \brief Value of plunging coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Value of plunging coordinate.
   */
  su2double GetAeroelastic_plunge(unsigned short val_marker);
  
  /*!
   * \brief Value of pitching coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \return Value of pitching coordinate.
   */
  su2double GetAeroelastic_pitch(unsigned short val_marker);
  
  /*!
   * \brief Value of plunging coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val - value of plunging coordinate.
   */
  void SetAeroelastic_plunge(unsigned short val_marker, su2double val);
  
  /*!
   * \brief Value of pitching coordinate.
   * \param[in] val_marker - the marker we are monitoring.
   * \param[in] val - value of pitching coordinate.
   */
  void SetAeroelastic_pitch(unsigned short val_marker, su2double val);
  
  /*!
   * \brief Get information about the aeroelastic simulation.
   * \return <code>TRUE</code> if it is an aeroelastic case; otherwise <code>FALSE</code>.
   */
  bool GetAeroelastic_Simulation(void);
  
  /*!
   * \brief Get information about the wind gust.
   * \return <code>TRUE</code> if there is a wind gust; otherwise <code>FALSE</code>.
   */
  bool GetWind_Gust(void);
  
  /*!
   * \brief Get the type of gust to simulate.
   * \return type of gust to use for the simulation.
   */
  unsigned short GetGust_Type(void);
  
  /*!
   * \brief Get the gust direction.
   * \return the gust direction.
   */
  unsigned short GetGust_Dir(void);
  
  /*!
   * \brief Value of the gust wavelength.
   */
  su2double GetGust_WaveLength(void);
  
  /*!
   * \brief Value of the number of gust periods.
   */
  su2double GetGust_Periods(void);
  
  /*!
   * \brief Value of the gust amplitude.
   */
  su2double GetGust_Ampl(void);
  
  /*!
   * \brief Value of the time at which to begin the gust.
   */
  su2double GetGust_Begin_Time(void);
  
  /*!
   * \brief Value of the location ath which the gust begins.
   */
  su2double GetGust_Begin_Loc(void);
  
  /*!
   * \brief Get the number of iterations to evaluate the parametric coordinates.
   * \return Number of iterations to evaluate the parametric coordinates.
   */
  unsigned short GetnFFD_Iter(void);
  
  /*!
   * \brief Get the tolerance of the point inversion algorithm.
   * \return Tolerance of the point inversion algorithm.
   */
  su2double GetFFD_Tol(void);
  
  /*!
   * \brief Get the scale factor for the line search.
   * \return Scale factor for the line search.
   */
  su2double GetOpt_RelaxFactor(void);

  /*!
   * \brief Get the bound for the line search.
   * \return Bound for the line search.
   */
  su2double GetOpt_LineSearch_Bound(void);
  
  /*!
   * \brief Set the scale factor for the line search.
   * \param[in] val_scale - scale of the deformation.
   */
  void SetOpt_RelaxFactor(su2double val_scale);
  
  /*!
   * \brief Get the node number of the CV to visualize.
   * \return Node number of the CV to visualize.
   */
  long GetVisualize_CV(void);
  
  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetFixed_CL_Mode(void);
  
  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetFixed_CM_Mode(void);
  
  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetEval_dOF_dCX(void);
  
  /*!
   * \brief Get information about whether to use fixed CL mode.
   * \return <code>TRUE</code> if fixed CL mode is active; otherwise <code>FALSE</code>.
   */
  bool GetDiscard_InFiles(void);
  
  /*!
   * \brief Get the value specified for the target CL.
   * \return Value of the target CL.
   */
  su2double GetTarget_CL(void);
  
  /*!
   * \brief Get the value for the lift curve slope for fixed CL mode.
   * \return Lift curve slope for fixed CL mode.
   */
  su2double GetdCL_dAlpha(void);
  
  /*!
   * \brief Get the value of iterations to re-evaluate the angle of attack.
   * \return Number of iterations.
   */
  unsigned long GetUpdate_Alpha(void);
  
  /*!
   * \brief Number of iterations to evaluate dCL_dAlpha.
   * \return Number of iterations.
   */
  unsigned long GetIter_dCL_dAlpha(void);
  
  /*!
   * \brief Get the value of the damping coefficient for fixed CL mode.
   * \return Damping coefficient for fixed CL mode.
   */
  su2double GetdCM_diH(void);
  
  /*!
   * \brief Get the value of iterations to re-evaluate the angle of attack.
   * \return Number of iterations.
   */
  unsigned long GetIter_Fixed_CL(void);
  
  /*!
   * \brief Get the value of iterations to re-evaluate the angle of attack.
   * \return Number of iterations.
   */
  unsigned long GetIter_Fixed_NetThrust(void);
  
  /*!
   * \brief Get the value of the damping coefficient for fixed CL mode.
   * \return Damping coefficient for fixed CL mode.
   */
  su2double GetdNetThrust_dBCThrust(void);
  
  /*!
   * \brief Get the value of iterations to re-evaluate the angle of attack.
   * \return Number of iterations.
   */
  unsigned long GetUpdate_BCThrust(void);
  
  /*!
   * \brief Set the value of the boolean for updating AoA in fixed lift mode.
   * \param[in] val_update - the bool for whether to update the AoA.
   */
  void SetUpdate_BCThrust_Bool(bool val_update);
  
  /*!
   * \brief Set the value of the boolean for updating AoA in fixed lift mode.
   * \param[in] val_update - the bool for whether to update the AoA.
   */
  void SetUpdate_AoA(bool val_update);
  
  /*!
   * \brief Get information about whether to update the AoA for fixed lift mode.
   * \return <code>TRUE</code> if we should update the AoA for fixed lift mode; otherwise <code>FALSE</code>.
   */
  bool GetUpdate_BCThrust_Bool(void);
  
  /*!
   * \brief Get information about whether to update the AoA for fixed lift mode.
   * \return <code>TRUE</code> if we should update the AoA for fixed lift mode; otherwise <code>FALSE</code>.
   */
  bool GetUpdate_AoA(void);
  
  /*!
   * \brief Set the current number of non-physical nodes in the solution.
   * \param[in] val_nonphys_points - current number of non-physical points.
   */
  void SetNonphysical_Points(unsigned long val_nonphys_points);
  
  /*!
   * \brief Get the current number of non-physical nodes in the solution.
   * \return Current number of non-physical points.
   */
  unsigned long GetNonphysical_Points(void);
  
  /*!
   * \brief Set the current number of non-physical reconstructions for 2nd-order upwinding.
   * \param[in] val_nonphys_reconstr - current number of non-physical reconstructions for 2nd-order upwinding.
   */
  void SetNonphysical_Reconstr(unsigned long val_nonphys_reconstr);
  
  /*!
   * \brief Get the current number of non-physical reconstructions for 2nd-order upwinding.
   * \return Current number of non-physical reconstructions for 2nd-order upwinding.
   */
  unsigned long GetNonphysical_Reconstr(void);
  
  /*!
   * \brief Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
   x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
   function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
   the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
   ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
   condition for a natural spline, with zero second derivative on that boundary.
   Numerical Recipes: The Art of Scientific Computing, Third Edition in C++.
   */
  void SetSpline(vector<su2double> &x, vector<su2double> &y, unsigned long n, su2double yp1, su2double ypn, vector<su2double> &y2);
  
  /*!
   * \brief Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
   and given the array y2a[1..n], which is the output from spline above, and given a value of
   x, this routine returns a cubic-spline interpolated value y.
   Numerical Recipes: The Art of Scientific Computing, Third Edition in C++.
   * \returns The interpolated value of for x.
   */
  su2double GetSpline(vector<su2double> &xa, vector<su2double> &ya, vector<su2double> &y2a, unsigned long n, su2double x);

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
  void Tock(double val_start_time, string val_function_name, int val_group_id);

  /*!
   * \brief Write a CSV file containing the results of the profiling.
   */
  void SetProfilingCSV(void);

  /*!
   * \brief Start the timer for profiling subroutines.
   * \param[in] val_start_time - the value of the start time.
   */
  void GEMM_Tick(double *val_start_time);

  /*!
   * \brief Stop the timer for the GEMM profiling and store results.
   * \param[in] val_start_time - The value of the start time.
   * \param[in] M, N, K        - Matrix size of the GEMM call.
   */
  void GEMM_Tock(double val_start_time, int M, int N, int K);

  /*!
   * \brief Write a CSV file containing the results of the profiling.
   */
  void GEMMProfilingCSV(void);

  /*!
   *
   * \brief Set freestream turbonormal for initializing solution.
   */
  void SetFreeStreamTurboNormal(su2double* turboNormal);

  /*!
   *
   * \brief Set freestream turbonormal for initializing solution.
   */
  su2double* GetFreeStreamTurboNormal(void);

  /*!
   *
   * \brief Set multizone properties.
   */
  void SetMultizone(CConfig *driver_config, CConfig **config_container);

  /*!
   * \brief Get the verbosity level of the console output.
   * \return Verbosity level for the console output.
   */
  unsigned short GetConsole_Output_Verb(void);
  
  /*!
   * \brief Get the kind of marker analyze marker (area-averaged, mass flux averaged, etc).
   * \return Kind of average.
   */
  unsigned short GetKind_Average(void);

  /*!
   *
   * \brief Get the direct differentation method.
   * \return direct differentiation method.
   */
  unsigned short GetDirectDiff();
  
  /*!
   * \brief Get the indicator whether we are solving an discrete adjoint problem.
   * \return the discrete adjoint indicator.
   */
  bool GetDiscrete_Adjoint(void);
  
  /*!
   * \brief Get the indicator whether we want to benchmark the MPI performance of FSI problems
   * \return The value for checking
   */
  bool CheckFSI_MPI(void);
  
  /*!
   * \brief Get the number of fluid subiterations roblems.
   * \return Number of FSI subiters.
   */
  unsigned short GetnIterFSI(void);
  
  /*!
   * \brief Get the number of subiterations while a ramp is applied.
   * \return Number of FSI subiters.
   */
  unsigned short GetnIterFSI_Ramp(void);
  
  /*!
   * \brief Get Aitken's relaxation parameter for static relaxation cases.
   * \return Aitken's relaxation parameters.
   */
  su2double GetAitkenStatRelax(void);
  
  /*!
   * \brief Get Aitken's maximum relaxation parameter for dynamic relaxation cases and first iteration.
   * \return Aitken's relaxation parameters.
   */
  su2double GetAitkenDynMaxInit(void);
  
  /*!
   * \brief Get Aitken's maximum relaxation parameter for dynamic relaxation cases and first iteration.
   * \return Aitken's relaxation parameters.
   */
  su2double GetAitkenDynMinInit(void);
  
  
  /*!
   * \brief Decide whether to apply dead loads to the model.
   * \return <code>TRUE</code> if the dead loads are to be applied, <code>FALSE</code> otherwise.
   */
  
  bool GetDeadLoad(void);
  
  /*!
   * \brief Identifies if the mesh is matching or not (temporary, while implementing interpolation procedures).
   * \return <code>TRUE</code> if the mesh is matching, <code>FALSE</code> otherwise.
   */
  
  bool GetPseudoStatic(void);
  
  /*!
   * \brief Identifies if we want to restart from a steady or an unsteady solution.
   * \return <code>TRUE</code> if we restart from steady state solution, <code>FALSE</code> otherwise.
   */
  
  bool GetSteadyRestart(void);
  
  
  /*!
   * \brief Provides information about the time integration of the structural analysis, and change the write in the output
   *        files information about the iteration.
   * \return The kind of time integration: Static or dynamic analysis
   */
  unsigned short GetDynamic_Analysis(void);
  
  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetDelta_DynTime(void);
  
  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetTotal_DynTime(void);
  
  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation (non dimensional).
   */
  su2double GetCurrent_DynTime(void);
  
  /*!
   * \brief Get the current instance.
   * \return Current instance identifier.
   */
  unsigned short GetiInst(void);

  /*!
   * \brief Set the current instance.
   * \param[in] iInst - current instance identifier.
   */
  void SetiInst(unsigned short val_iInst);

  /*!
   * \brief Get information about writing dynamic structural analysis headers and file extensions.
   * \return 	<code>TRUE</code> means that dynamic structural analysis solution files will be written.
   */
  bool GetWrt_Dynamic(void);
  
  /*!
   * \brief Get Newmark alpha parameter.
   * \return Value of the Newmark alpha parameter.
   */
  su2double GetNewmark_beta(void);
  
  /*!
   * \brief Get Newmark delta parameter.
   * \return Value of the Newmark delta parameter.
   */
  su2double GetNewmark_gamma(void);
  
  /*!
   * \brief Get the number of integration coefficients provided by the user.
   * \return Number of integration coefficients.
   */
  unsigned short GetnIntCoeffs(void);
  
  /*!
   * \brief Get the number of different values for the elasticity modulus.
   * \return Number of different values for the elasticity modulus.
   */
  unsigned short GetnElasticityMod(void);
  
  /*!
   * \brief Get the number of different values for the Poisson ratio.
   * \return Number of different values for the Poisson ratio.
   */
  unsigned short GetnPoissonRatio(void);
  
  /*!
   * \brief Get the number of different values for the Material density.
   * \return Number of different values for the Material density.
   */
  unsigned short GetnMaterialDensity(void);
  
  /*!
   * \brief Get the integration coefficients for the Generalized Alpha - Newmark integration integration scheme.
   * \param[in] val_coeff - Index of the coefficient.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double Get_Int_Coeffs(unsigned short val_coeff);
  
  /*!
   * \brief Get the number of different values for the modulus of the electric field.
   * \return Number of different values for the modulus of the electric field.
   */
  unsigned short GetnElectric_Field(void);
  
  /*!
   * \brief Get the dimensionality of the electric field.
   * \return Number of integration coefficients.
   */
  unsigned short GetnDim_Electric_Field(void);
  
  /*!
   * \brief Get the values for the electric field modulus.
   * \param[in] val_coeff - Index of the coefficient.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double Get_Electric_Field_Mod(unsigned short val_coeff);
  
  /*!
   * \brief Set the values for the electric field modulus.
   * \param[in] val_coeff - Index of the electric field.
   * \param[in] val_el_field - Value of the electric field.
   */
  void Set_Electric_Field_Mod(unsigned short val_coeff, su2double val_el_field);
  
  /*!
   * \brief Get the direction of the electric field in reference configuration.
   * \param[in] val_coeff - Index of the coefficient.
   * \return Alpha coefficient for the Runge-Kutta integration scheme.
   */
  su2double* Get_Electric_Field_Dir(void);
  
  
  /*!
   * \brief Check if the user wants to apply the load as a ramp.
   * \return 	<code>TRUE</code> means that the load is to be applied as a ramp.
   */
  bool GetRamp_Load(void);
  
  /*!
   * \brief Get the maximum time of the ramp.
   * \return 	Value of the max time while the load is linearly increased
   */
  su2double GetRamp_Time(void);
  
  /*!
   * \brief Check if the user wants to apply the load as a ramp.
   * \return  <code>TRUE</code> means that the load is to be applied as a ramp.
   */
  bool GetRampAndRelease_Load(void);

  /*!
   * \brief Check if the user wants to apply the load as a ramp.
   * \return  <code>TRUE</code> means that the load is to be applied as a ramp.
   */
  bool GetSine_Load(void);

  /*!
   * \brief Get the sine load properties.
   * \param[in] val_index - Index corresponding to the load boundary.
   * \return The pointer to the sine load values.
   */
  su2double* GetLoad_Sine(void);
   /*!
    * \brief Get the kind of load transfer method we want to use for dynamic problems
    * \note This value is obtained from the config file, and it is constant
    *       during the computation.
    * \return Kind of transfer method for multiphysics problems
    */
   unsigned short GetDynamic_LoadTransfer(void);
  
   /*!
    * \brief Get the penalty weight value for the objective function.
    * \return  Penalty weight value for the reference geometry objective function.
    */
   su2double GetRefGeom_Penalty(void);
  
   /*!
    * \brief Get the penalty weight value for the objective function.
    * \return  Penalty weight value for the reference geometry objective function.
    */
   su2double GetTotalDV_Penalty(void);
  
   /*!
    * \brief Get whether a predictor is used for FSI applications.
    * \return  Bool: determines if predictor is used or not
    */
   bool GetPredictor(void);

  /*!
   * \brief Get the order of the predictor for FSI applications.
   * \return 	Order of predictor
   */
  unsigned short GetPredictorOrder(void);

  /*!
   * \brief Get boolean for using Persson's shock capturing method in Euler flow DG-FEM
   * \return Boolean for using Persson's shock capturing method in Euler flow DG-FEM
   */
  bool GetEulerPersson(void);

  /*!
   * \brief Set boolean for using Persson's shock capturing method in Euler flow DG-FEM
   * \param[in] val_EulerPersson - Boolean for using Persson's shock capturing method in Euler flow DG-FEM
   */
  void SetEulerPersson(bool val_EulerPersson);

  /*!
   * \brief Get whether a relaxation parameter is used for FSI applications.
   * \return Bool: determines if relaxation parameter  is used or not
   */
  bool GetRelaxation(void);

  /*!
   * \brief Check if the simulation we are running is a FSI simulation
   * \return Value of the physical time in an unsteady simulation.
   */
  bool GetFSI_Simulation(void);
  
  /*!
   * \brief Set that the simulation we are running is a FSI simulation
   * \param[in] FSI_sim - boolean that determines is FSI_Problem is true/false.
   */
  void SetFSI_Simulation(bool FSI_sim);

  /*!
   * \brief Set that the simulation we are running is a multizone simulation
   * \param[in] MZ_problem - boolean that determines is Multizone_Problem is true/false.
   */
  void SetMultizone_Problem(bool MZ_problem);

  /*!
   * \brief Get whether the simulation we are running is a multizone simulation
   * \return Multizone_Problem - boolean that determines is Multizone_Problem is true/false.
   */
  bool GetMultizone_Problem(void);

   /*!
    * \brief Get the ID for the FEA region that we want to compute the gradient for using direct differentiation
    * \return ID
    */
   unsigned short GetnID_DV(void);
  
  /*!
   * \brief Check if we want to apply an incremental load to the nonlinear structural simulation
   * \return <code>TRUE</code> means that the load is to be applied in increments.
   */
  bool GetIncrementalLoad(void);
  
  /*!
   * \brief Get the number of increments for an incremental load.
   * \return 	Number of increments.
   */
  unsigned long GetNumberIncrements(void);

  /*!
   * \brief Get the value of the criteria for applying incremental loading.
   * \return Value of the log10 of the residual.
   */
  su2double GetIncLoad_Criteria(unsigned short val_var);
  
  /*!
   * \brief Get the relaxation method chosen for the simulation
   * \return Value of the relaxation method
   */
  unsigned short GetRelaxation_Method_FSI(void);

  /*!
   * \brief Get the kind of Riemann solver for the DG method (FEM flow solver).
   * \note This value is obtained from the config file, and it is constant
   *       during the computation.
   * \return Kind of Riemann solver for the DG method (FEM flow solver).
   */
  unsigned short GetRiemann_Solver_FEM(void);

  /*!
   * \brief Get the factor applied during quadrature of straight elements.
   * \return The specified straight element quadrature factor.
   */
  su2double GetQuadrature_Factor_Straight(void);

  /*!
   * \brief Get the factor applied during quadrature of curved elements.
   * \return The specified curved element quadrature factor.
   */
  su2double GetQuadrature_Factor_Curved(void);

  /*!
   * \brief Get the factor applied during time quadrature for ADER-DG.
   * \return The specified ADER-DG time quadrature factor.
   */
  su2double GetQuadrature_Factor_Time_ADER_DG(void);

  /*!
   * \brief Function to make available the multiplication factor theta of the
            symmetrizing terms in the DG discretization of the viscous terms.
   * \return The specified factor for the DG discretization.
   */
  su2double GetTheta_Interior_Penalty_DGFEM(void);

  /*!
   * \brief Function to make available the matrix size in vectorization in
            order to optimize the gemm performance.
   * \return The matrix size in this direction.
   */
  unsigned short GetSizeMatMulPadding(void);

  /*!
   * \brief Function to make available whether or not the entropy must be computed.
   * \return The boolean whether or not the entropy must be computed.
   */
  bool GetCompute_Entropy(void);

  /*!
   * \brief Function to make available whether or not the lumped mass matrix
            must be used for steady computations.
   * \return The boolean whether or not to use the lumped mass matrix.
   */
  bool GetUse_Lumped_MassMatrix_DGFEM(void);

  /*!
   * \brief Function to make available whether or not only the exact Jacobian
            of the spatial discretization must be computed.
   * \return The boolean whether or not the Jacobian must be computed.
   */
  bool GetJacobian_Spatial_Discretization_Only(void);

  /*!
   * \brief Get the interpolation method used for matching between zones.
   */
  inline unsigned short GetKindInterpolation(void);
	
	/*!
	 * \brief Get option of whether to use conservative interpolation between zones.
	 */
  inline bool GetConservativeInterpolation(void);
	
  /*!
   * \brief Get the basis function to use for radial basis function interpolation for FSI.
   */
  inline unsigned short GetKindRadialBasisFunction(void);
	
  /*!
   * \brief Get option of whether to use polynomial terms in Radial Basis Function interpolation.
   */
  inline bool GetRadialBasisFunctionPolynomialOption(void);
	
  /*!
   * \brief Get the basis function radius to use for radial basis function interpolation for FSI.
   */
  inline su2double GetRadialBasisFunctionParameter(void);

  /*!
   * \brief Get information about using UQ methodology
   * \return <code>TRUE</code> means that UQ methodology of eigenspace perturbation will be used
   */
  bool GetUsing_UQ(void);

  /*!
   * \brief Get the amount of eigenvalue perturbation to be done
   * \return Value of the uq_delta_b parameter
   */
  su2double GetUQ_Delta_B(void);

  /*!
   * \brief Get the kind of eigenspace perturbation to be done
   * \return Value of the eig_val_comp
   */
  unsigned short GetEig_Val_Comp(void);

  /*!
   * \brief Get the underelaxation factor
   * \return Value of the uq_urlx parameter
   */
  su2double GetUQ_URLX(void);

  /*!
   * \brief Get information about eigenspace perturbation
   * \return <code>TRUE</code> means eigenspace perterturbation will be used
   */
  bool GetUQ_Permute(void);
  
  /*!
   * \brief Get information about whether to use wall functions.
   * \return <code>TRUE</code> if wall functions are on; otherwise <code>FALSE</code>.
   */
  bool GetWall_Functions(void);
  /*!
   * \brief Get the AD support.
   */
  bool GetAD_Mode(void);

  /*!
   * \brief Set the maximum velocity^2 in the domain for the incompressible preconditioner.
   * \param[in] Value of the maximum velocity^2 in the domain for the incompressible preconditioner.
   */
  void SetMax_Vel2(su2double val_max_vel2);

  /*!
   * \brief Get the maximum velocity^2 in the domain for the incompressible preconditioner.
   * \return Value of the maximum velocity^2 in the domain for the incompressible preconditioner.
   */
  su2double GetMax_Vel2(void);
  
  /*!
   * \brief Set the sum of the bandwidth for writing binary restarts (to be averaged later).
   * \param[in] Sum of the bandwidth for writing binary restarts.
   */
  void SetRestart_Bandwidth_Agg(su2double val_restart_bandwidth_sum);
  
  /*!
   * \brief Set the sum of the bandwidth for writing binary restarts (to be averaged later).
   * \return Sum of the bandwidth for writing binary restarts.
   */
  su2double GetRestart_Bandwidth_Agg(void);

  /*!
   * \brief Get the frequency for writing the surface solution file in Dual Time.
   * \return It writes the surface solution file with this frequency.
   */
  unsigned long GetWrt_Surf_Freq_DualTime(void);
    
  /*!
   * \brief Get the Kind of Hybrid RANS/LES.
   * \return Value of Hybrid RANS/LES method.
   */
  unsigned short GetKind_HybridRANSLES(void);

  /*!
   * \brief Get the Kind of Roe Low Dissipation Scheme for Unsteady flows.
   * \return Value of Low dissipation approach.
   */
   unsigned short GetKind_RoeLowDiss(void);
    
  /*!
   * \brief Get the DES Constant.
   * \return Value of DES constant.
   */
   su2double GetConst_DES(void);

  /*!
   * \brief Get QCR (SA-QCR2000).
   */
  bool GetQCR(void);

  /*!
   * \brief Get if AD preaccumulation should be performed.
   */
  bool GetAD_Preaccumulation(void);

  /*!
   * \brief Get the heat equation.
   * \return YES if weakly coupled heat equation for inc. flow is enabled.
   */
  bool GetWeakly_Coupled_Heat(void);

  /*!
   * \brief Check if values passed to the BC_HeatFlux-Routine are already integrated.
   * \return YES if the passed values is the integrated heat flux over the marker's surface.
   */
  bool GetIntegrated_HeatFlux(void);
  
  /*!
   * \brief Get Compute Average.
   * \return YES if start computing averages
   */
  bool GetCompute_Average(void);

  /*!
   * \brief Get the verification solution.
   * \return The verification solution to be used.
   */
  unsigned short GetVerification_Solution(void);
  
  /*!
   * \brief Get topology optimization.
   */
  bool GetTopology_Optimization(void) const;

  /*!
   * \brief Get name of output file for topology optimization derivatives.
   */
  string GetTopology_Optim_FileName(void) const;

  /*!
   * \brief Get exponent for density-based stiffness penalization.
   */
  su2double GetSIMP_Exponent(void) const;

  /*!
   * \brief Get lower bound for density-based stiffness penalization.
   */
  su2double GetSIMP_MinStiffness(void) const;
  
  /*!
   * \brief Number of kernels to use in filtering the design density field.
   */
  unsigned short GetTopology_Optim_Num_Kernels(void) const;
  
  /*!
   * \brief Get the i'th kernel to use, its parameter, and the radius.
   */
  void GetTopology_Optim_Kernel(const unsigned short iKernel, unsigned short &type,
                                su2double &param, su2double &radius) const;
  /*!
   * \brief Get the type and parameter for the projection function used in topology optimization
   */
  void GetTopology_Optim_Projection(unsigned short &type, su2double &param) const;

  /*!
   * \brief Retrieve the ofstream of the history file for the current zone.
   */
  ofstream* GetHistFile(void);

  /*!
   * \brief Set the ofstream of the history file for the current zone.
   */
  void SetHistFile(ofstream *HistFile);

  /*!
   * \brief Get the filenames of the individual config files
   * \return File name of the config file for zone "index"
   */
  string GetConfigFilename(unsigned short index);

  /*!
   * \brief Get the number of config files
   * \return Number of config filenames in CONFIG_LIST
   */
  unsigned short GetnConfigFiles(void);

  /*!
   * \brief Check if the multizone problem is solved for time domain.
   * \return YES if time-domain is considered.
   */
  bool GetTime_Domain(void);

  /*!
   * \brief Get the number of inner iterations
   * \return Number of inner iterations on each multizone block
   */
  unsigned long GetnInner_Iter(void);

  /*!
   * \brief Get the number of outer iterations
   * \return Number of outer iterations for the multizone problem
   */
  unsigned long GetnOuter_Iter(void);

  /*!
   * \brief Get the number of time iterations
   * \return Number of time steps run for the multizone problem
   */
  unsigned long GetnTime_Iter(void);

  /*!
   * \brief Get the number of pseudo-time iterations
   * \return Number of pseudo-time steps run for the single-zone problem
   */
  unsigned long GetnIter(void);

  /*!
   * \brief Get the restart iteration
   * \return Iteration for the restart of multizone problems
   */
  unsigned long GetRestart_Iter(void);

  /*!
   * \brief Get the time step for multizone problems
   * \return Time step for multizone problems, it is set on all the zones
   */
  su2double GetTime_Step(void);

  /*!
   * \brief Get the maximum simulation time for time-domain problems
   * \return Simulation time for multizone problems, it is set on all the zones
   */
  su2double GetMax_Time(void);

  /*!
   * \brief Get the level of MPI communications to be performed.
   * \return Level of MPI communications.
   */
  unsigned short GetComm_Level(void);
  
  /*
   * \brief Check if the mesh read supports multiple zones.
   * \return YES if multiple zones can be contained in the mesh file.
   */
  bool GetMultizone_Mesh(void);

  /*!
   * \brief Check if the mesh read supports multiple zones.
   * \return YES if multiple zones can be contained in the mesh file.
   */
  bool GetMultizone_Residual(void);

  /*!
   * \brief Check if the (new) single-zone driver is to be used (temporary)
   * \return YES if the (new) single-zone driver is to be used.
   */
  bool GetSinglezone_Driver(void);

  /*!
   * \brief Check if the special output is written
   * \return YES if the special output is written.
   */
  bool GetSpecial_Output(void);

  /*!
   * \brief Check if the forces breakdown file is written
   * \return YES if the forces breakdown file is written.
   */
  bool GetWrt_ForcesBreakdown(void);

};

#include "config_structure.inl"
