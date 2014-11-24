/*!
 * \file config_structure.hpp
 * \brief All the information about the definition of the physical problem.
 *        The subroutines and functions are in the <i>config_structure.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.4 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <assert.h>

#include "./option_structure.hpp"

using namespace std;

/*!
 * \class CConfig
 * \brief Main class for defining the problem; basically this class reads the configuration file, and
 *        stores all the information.
 * \author F. Palacios.
 * \version 3.2.4 "eagle"
 */
class CConfig {
private:
	unsigned short Kind_SU2; /*!< \brief Kind of SU2 software component. */
	unsigned short iZone, nZone; /*!< \brief Number of zones in the mesh. */
	double OrderMagResidual; /*!< \brief Order of magnitude reduction. */
	double MinLogResidual; /*!< \brief Minimum value of the log residual. */
	double EA_ScaleFactor; /*!< \brief Equivalent Area scaling factor */
	double* EA_IntLimit; /*!< \brief Integration limits of the Equivalent Area computation */
  double AdjointLimit; /*!< \brief Adjoint variable limit */
  double* Subsonic_Engine_Box; /*!< \brief Coordinates of the box subsonic region */
	double* Hold_GridFixed_Coord; /*!< \brief Coordinates of the box to hold fixed the nbumerical grid */
	unsigned short ConvCriteria;	/*!< \brief Kind of convergence criteria. */
  unsigned short nFFD_Iter; 	/*!< \brief Iteration for the point inversion problem. */
  double FFD_Tol;  	/*!< \brief Tolerance in the point inversion problem. */
	bool Adjoint,			/*!< \brief Flag to know if the code is solving an adjoint problem. */
    Viscous,                /*!< \brief Flag to know if the code is solving a viscous problem. */
	EquivArea,				/*!< \brief Flag to know if the code is going to compute and plot the equivalent area. */
	InvDesign_Cp,				/*!< \brief Flag to know if the code is going to compute and plot the inverse design. */
  InvDesign_HeatFlux,				/*!< \brief Flag to know if the code is going to compute and plot the inverse design. */
	OneShot,				/*!< \brief Flag to know if the code is solving a one shot problem. */
	Linearized,				/*!< \brief Flag to know if the code is solving a linearized problem. */
	Grid_Movement,			/*!< \brief Flag to know if there is grid movement. */
    Wind_Gust,              /*!< \brief Flag to know if there is a wind gust. */
    Aeroelastic_Simulation, /*!< \brief Flag to know if there is an aeroelastic simulation. */
	Rotating_Frame,			/*!< \brief Flag to know if there is a rotating frame. */
	PoissonSolver,			/*!< \brief Flag to know if we are solving  poisson forces  in plasma solver. */
	Low_Mach_Precon,		/*!< \brief Flag to know if we are using a low Mach number preconditioner. */
	GravityForce,			/*!< \brief Flag to know if the gravity force is incuded in the formulation. */
	SmoothNumGrid,			/*!< \brief Smooth the numerical grid. */
	AdaptBoundary,			/*!< \brief Adapt the elements on the boundary. */
	FullMG,					/*!< \brief Full multigrid strategy. */
	Divide_Element,			/*!< \brief Divide rectables and hexahedrom. */
	Engine_Intake,			/*!< \brief Engine intake subsonic region. */
	Frozen_Visc,			/*!< \brief Flag for adjoint problem with/without frozen viscosity. */
	Sens_Remove_Sharp,			/*!< \brief Flag for removing or not the sharp edges from the sensitivity computation. */
	Hold_GridFixed,	/*!< \brief Flag hold fixed some part of the mesh during the deformation. */
	Axisymmetric, /*!< \brief Flag for axisymmetric calculations */
	DebugMode, /*!< \brief Flag for debug mode */
	Show_Adj_Sens, /*!< \brief Flag for outputting sensitivities on exit */
  ionization;  /*!< \brief Flag for determining if free electron gas is in the mixture */
	bool Visualize_Partition;	/*!< \brief Flag to visualize each partition in the DDM. */
  double Damp_Engine_Inflow;	/*!< \brief Damping factor for the engine inlet. */
	double Damp_Res_Restric,	/*!< \brief Damping factor for the residual restriction. */
	Damp_Correc_Prolong; /*!< \brief Damping factor for the correction prolongation. */
	double Position_Plane; /*!< \brief Position of the Near-Field (y coordinate 2D, and z coordinate 3D). */
	double WeightCd; /*!< \brief Weight of the drag coefficient. */
	unsigned short Unsteady_Simulation;	/*!< \brief Steady or unsteady (time stepping or dual time stepping) computation. */
	unsigned short nStartUpIter;	/*!< \brief Start up iterations using the fine grid. */
  double FixAzimuthalLine; /*!< \brief Fix an azimuthal line due to misalignments of the nearfield. */
	double *DV_Value;		/*!< \brief Previous value of the design variable. */
	double LimiterCoeff;				/*!< \brief Limiter coefficient */
  unsigned long LimiterIter;	/*!< \brief Freeze the value of the limiter after a number of iterations */
	double SharpEdgesCoeff;				/*!< \brief Coefficient to identify the limit of a sharp edge. */
  unsigned short SystemMeasurements; /*!< \brief System of measurements. */
  unsigned short Kind_Regime;  /*!< \brief Kind of adjoint function. */
  unsigned short Kind_ObjFunc;  /*!< \brief Kind of objective function. */
  unsigned short Kind_SensSmooth; /*!< \brief Kind of sensitivity smoothing technique. */
  unsigned short Continuous_Eqns; /*!< \brief Which equations to treat continuously (Hybrid adjoint)*/
  unsigned short Discrete_Eqns; /*!< \brief Which equations to treat discretely (Hybrid adjoint). */
	unsigned short *Design_Variable; /*!< \brief Kind of design variable. */
	double RatioDensity,				/*!< \brief Ratio of density for a free surface problem. */
	RatioViscosity,				/*!< \brief Ratio of viscosity for a free surface problem. */
	FreeSurface_Thickness,  /*!< \brief Thickness of the interfase for a free surface problem. */
	FreeSurface_Outlet,  /*!< \brief Outlet of the interfase for a free surface problem. */
	FreeSurface_Damping_Coeff,  /*!< \brief Damping coefficient of the free surface for a free surface problem. */
	FreeSurface_Damping_Length;  /*!< \brief Damping length of the free surface for a free surface problem. */
	unsigned short Kind_Adaptation;	/*!< \brief Kind of numerical grid adaptation. */
	unsigned short nTimeInstances;  /*!< \brief Number of periodic time instances for Time Spectral integration. */
	double TimeSpectral_Period;		/*!< \brief Period of oscillation to be used with time-spectral computations. */
	double New_Elem_Adapt;			/*!< \brief Elements to adapt in the numerical grid adaptation process. */
	double Delta_UnstTime,			/*!< \brief Time step for unsteady computations. */
	Delta_UnstTimeND;						/*!< \brief Time step for unsteady computations (non dimensional). */
	double Total_UnstTime,						/*!< \brief Total time for unsteady computations. */
	Total_UnstTimeND;								/*!< \brief Total time for unsteady computations (non dimensional). */
	double Current_UnstTime,									/*!< \brief Global time of the unsteady simulation. */
	Current_UnstTimeND;									/*!< \brief Global time of the unsteady simulation. */
	unsigned short nMarker_Euler,	/*!< \brief Number of Euler wall markers. */
	nMarker_FarField,				/*!< \brief Number of far-field markers. */
	nMarker_Custom,
	nMarker_SymWall,				/*!< \brief Number of symmetry wall markers. */
  nMarker_Pressure,				/*!< \brief Number of pressure wall markers. */
	nMarker_PerBound,				/*!< \brief Number of periodic boundary markers. */
	nMarker_NearFieldBound,				/*!< \brief Number of near field boundary markers. */
  nMarker_ActDisk_Inlet, nMarker_ActDisk_Outlet,
	nMarker_InterfaceBound,				/*!< \brief Number of interface boundary markers. */
	nMarker_Dirichlet,				/*!< \brief Number of interface boundary markers. */
	nMarker_Dirichlet_Elec,				/*!< \brief Number of interface boundary markers. */
	nMarker_Inlet,					/*!< \brief Number of inlet flow markers. */
	nMarker_Riemann,					/*!< \brief Number of Riemann flow markers. */
	nMarker_Supersonic_Inlet,					/*!< \brief Number of supersonic inlet flow markers. */
	nMarker_Outlet,					/*!< \brief Number of outlet flow markers. */
	nMarker_Out_1D,         /*!< \brief Number of outlet flow markers over which to calculate 1D outputs */
	nMarker_Isothermal,     /*!< \brief Number of isothermal wall boundaries. */
  nMarker_IsothermalNonCatalytic, /*!< \brief Number of constant temperature wall boundaries. */
  nMarker_IsothermalCatalytic, /*!< \brief Number of constant temperature wall boundaries. */
	nMarker_HeatFlux,       /*!< \brief Number of constant heat flux wall boundaries. */
  nMarker_HeatFluxNonCatalytic, /*!< \brief Number of constant heat flux wall boundaries. */
  nMarker_HeatFluxCatalytic, /*!< \brief Number of constant heat flux wall boundaries. */
	nMarker_EngineExhaust,					/*!< \brief Number of nacelle exhaust flow markers. */
	nMarker_EngineInflow,					/*!< \brief Number of nacelle inflow flow markers. */
  nMarker_EngineBleed,					/*!< \brief Number of nacelle inflow flow markers. */
  nMarker_Displacement,					/*!< \brief Number of displacement surface markers. */
	nMarker_Load,					/*!< \brief Number of load surface markers. */
	nMarker_FlowLoad,					/*!< \brief Number of load surface markers. */
	nMarker_Neumann,				/*!< \brief Number of Neumann flow markers. */
	nMarker_Neumann_Elec,				/*!< \brief Number of Neumann flow markers. */
	nMarker_All,					/*!< \brief Total number of markers using the grid information. */
	nMarker_Config;					/*!< \brief Total number of markers using the config file
									(note that using parallel computation this number can be different
									from nMarker_All). */
	string *Marker_Euler,			/*!< \brief Euler wall markers. */
	*Marker_FarField,				/*!< \brief Far field markers. */
	*Marker_Custom,
	*Marker_SymWall,				/*!< \brief Symmetry wall markers. */
  *Marker_Pressure,				/*!< \brief Pressure boundary markers. */
	*Marker_PerBound,				/*!< \brief Periodic boundary markers. */
	*Marker_PerDonor,				/*!< \brief Rotationally periodic boundary donor markers. */
	*Marker_NearFieldBound,				/*!< \brief Near Field boundaries markers. */
	*Marker_InterfaceBound,				/*!< \brief Interface boundaries markers. */
  *Marker_ActDisk_Inlet,
  *Marker_ActDisk_Outlet,
	*Marker_Dirichlet,				/*!< \brief Interface boundaries markers. */
	*Marker_Dirichlet_Elec,				/*!< \brief Interface boundaries markers. */
	*Marker_Inlet,					/*!< \brief Inlet flow markers. */
	*Marker_Riemann,					/*!< \brief Riemann markers. */
	*Marker_Supersonic_Inlet,					/*!< \brief Supersonic inlet flow markers. */
	*Marker_Outlet,					/*!< \brief Outlet flow markers. */
	*Marker_Out_1D,         /*!< \brief Outlet flow markers over which to calculate 1D output. */
	*Marker_Isothermal,     /*!< \brief Isothermal wall markers. */
  *Marker_IsothermalNonCatalytic,     /*!< \brief Isothermal wall markers. */
  *Marker_IsothermalCatalytic,     /*!< \brief Isothermal wall markers. */
	*Marker_HeatFlux,       /*!< \brief Constant heat flux wall markers. */
  *Marker_HeatFluxNonCatalytic,       /*!< \brief Constant heat flux wall markers. */
  *Marker_HeatFluxCatalytic,       /*!< \brief Constant heat flux wall markers. */
	*Marker_EngineInflow,					/*!< \brief Engine Inflow flow markers. */
  *Marker_EngineBleed,					/*!< \brief Engine Inflow flow markers. */
  *Marker_EngineExhaust,					/*!< \brief Engine Exhaust flow markers. */
	*Marker_Displacement,					/*!< \brief Displacement markers. */
	*Marker_Load,					/*!< \brief Load markers. */
	*Marker_FlowLoad,					/*!< \brief Flow Load markers. */
	*Marker_Neumann,					/*!< \brief Neumann flow markers. */
	*Marker_Neumann_Elec,					/*!< \brief Neumann flow markers. */
	*Marker_All_TagBound;				/*!< \brief Global index for markers using grid information. */
	double *Dirichlet_Value;    /*!< \brief Specified Dirichlet value at the boundaries. */
	double *Nozzle_Ttotal;    /*!< \brief Specified total temperatures for nacelle boundaries. */
	double *Nozzle_Ptotal;    /*!< \brief Specified total pressures for nacelle boundaries. */
	double *Inlet_Ttotal;    /*!< \brief Specified total temperatures for inlet boundaries. */
	double *Riemann_Var1, *Riemann_Var2;    /*!< \brief Specified values for Riemann boundary. */
	double **Riemann_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for Riemann boundaries. */
	double *Inlet_Ptotal;    /*!< \brief Specified total pressures for inlet boundaries. */
	double **Inlet_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
	double *Inlet_Temperature;    /*!< \brief Specified temperatures for a supersonic inlet boundaries. */
	double *Inlet_Pressure;    /*!< \brief Specified static pressures for supersonic inlet boundaries. */
	double **Inlet_Velocity;  /*!< \brief Specified flow velocity vectors for supersonic inlet boundaries. */
	double *Inflow_Mach_Target;    /*!< \brief Specified fan face mach for nacelle boundaries. */
	double *Inflow_Mach;    /*!< \brief Specified fan face mach for nacelle boundaries. */
	double *Inflow_Pressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  double *Bleed_MassFlow_Target;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  double *Bleed_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  double *Bleed_Temperature_Target;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  double *Bleed_Temperature;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  double *Bleed_Pressure;    /*!< \brief Specified fan face mach for nacelle boundaries. */
  double *Outlet_Pressure;    /*!< \brief Specified back pressures (static) for outlet boundaries. */
	double *Isothermal_Temperature; /*!< \brief Specified isothermal wall temperatures (static). */
  double *Wall_Catalycity; /*!< \brief Specified wall species mass-fractions for catalytic boundaries. */
	double *Heat_Flux;  /*!< \brief Specified wall heat fluxes. */
  double *Heat_FluxNonCatalytic;  /*!< \brief Specified wall heat fluxes. */
  double *Heat_FluxCatalytic;  /*!< \brief Specified wall heat fluxes. */
	double *Displ_Value;    /*!< \brief Specified displacement for displacement boundaries. */
	double *Load_Value;    /*!< \brief Specified force for load boundaries. */
	double *FlowLoad_Value;    /*!< \brief Specified force for flow load boundaries. */
  double **ActDisk_Origin;
  double *ActDisk_RootRadius;
  double *ActDisk_TipRadius;
  double *ActDisk_CT;
  double *ActDisk_Omega;
  double **Periodic_RotCenter;  /*!< \brief Rotational center for each periodic boundary. */
	double **Periodic_RotAngles;      /*!< \brief Rotation angles for each periodic boundary. */
	double **Periodic_Translation;      /*!< \brief Translation vector for each periodic boundary. */
	unsigned short nPeriodic_Index;     /*!< \brief Number of SEND_RECEIVE periodic transformations. */
	double **Periodic_Center;         /*!< \brief Rotational center for each SEND_RECEIVE boundary. */
	double **Periodic_Rotation;      /*!< \brief Rotation angles for each SEND_RECEIVE boundary. */
	double **Periodic_Translate;      /*!< \brief Translation vector for each SEND_RECEIVE boundary. */
	string *Marker_CfgFile_TagBound;			/*!< \brief Global index for markers using config file. */
	unsigned short *Marker_All_KindBC,			/*!< \brief Global index for boundaries using grid information. */
	*Marker_CfgFile_KindBC;		/*!< \brief Global index for boundaries using config file. */
	short *Marker_All_SendRecv;		/*!< \brief Information about if the boundary is sended (+), received (-). */
	short *Marker_All_PerBound;	/*!< \brief Global index for periodic bc using the grid information. */
	unsigned long nExtIter;			/*!< \brief Number of external iterations. */
	unsigned long ExtIter;			/*!< \brief Current external iteration number. */
	unsigned long IntIter;			/*!< \brief Current internal iteration number. */
	unsigned long Unst_nIntIter;			/*!< \brief Number of internal iterations (Dual time Method). */
  long Unst_RestartIter;			/*!< \brief Iteration number to restart an unsteady simulation (Dual time Method). */
  long Unst_AdjointIter;			/*!< \brief Iteration number to begin the reverse time integration in the direct solver for the unsteady adjoint. */
	unsigned short nRKStep;			/*!< \brief Number of steps of the explicit Runge-Kutta method. */
	double *RK_Alpha_Step;			/*!< \brief Runge-Kutta beta coefficients. */
	unsigned short nMultiLevel;		/*!< \brief Number of multigrid levels (coarse levels). */
	unsigned short nCFL;			/*!< \brief Number of CFL, one for each multigrid level. */
	double
	CFLRedCoeff_Turb,		/*!< \brief CFL reduction coefficient on the LevelSet problem. */
	CFLRedCoeff_AdjFlow,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
	CFLRedCoeff_AdjTurb,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
	CFLFineGrid,		/*!< \brief CFL of the finest grid. */
  Max_DeltaTime,  		/*!< \brief Max delta time. */
	Unst_CFL;		/*!< \brief Unsteady CFL number. */
	bool AddIndNeighbor;			/*!< \brief Include indirect neighbor in the agglomeration process. */
	unsigned short nDV;		/*!< \brief Number of design variables. */
  unsigned short nGridMovement;		/*!< \brief Number of grid movement types specified. */
	unsigned short nParamDV;		/*!< \brief Number of parameters of the design variable. */
	double **ParamDV;				/*!< \brief Parameters of the design variable. */
  string *FFDTag;				/*!< \brief Parameters of the design variable. */
	unsigned short GeometryMode;			/*!< \brief Gemoetry mode (analysis or gradient computation). */
	unsigned short MGCycle;			/*!< \brief Kind of multigrid cycle. */
	unsigned short FinestMesh;		/*!< \brief Finest mesh for the full multigrid approach. */
	unsigned short nMG_PreSmooth,                 /*!< \brief Number of MG pre-smooth parameters found in config file. */
	nMG_PostSmooth,                             /*!< \brief Number of MG post-smooth parameters found in config file. */
	nMG_CorrecSmooth;                           /*!< \brief Number of MG correct-smooth parameters found in config file. */
	unsigned short *MG_PreSmooth,	/*!< \brief Multigrid Pre smoothing. */
	*MG_PostSmooth,					/*!< \brief Multigrid Post smoothing. */
	*MG_CorrecSmooth;					/*!< \brief Multigrid Jacobi implicit smoothing of the correction. */
	unsigned short Kind_Solver,	/*!< \brief Kind of solver Euler, NS, Continuous adjoint, etc.  */
	Kind_FluidModel,			/*!< \brief Kind of the Fluid Model: Ideal or Van der Walls, ... . */
	Kind_ViscosityModel,			/*!< \brief Kind of the Viscosity Model*/
	Kind_ConductivityModel,			/*!< \brief Kind of the Thermal Conductivity Model*/
	Kind_FreeStreamOption,			/*!< \brief Kind of free stream option to choose if initializing with density or temperature  */
	Kind_InitOption,			/*!< \brief Kind of Init option to choose if initializing with Reynolds number or with thermodynamic conditions   */
	Kind_GasModel,				/*!< \brief Kind of the Gas Model. */
	*Kind_GridMovement,    /*!< \brief Kind of the unsteady mesh movement. */
	Kind_Gradient_Method,		/*!< \brief Numerical method for computation of spatial gradients. */
	Kind_Linear_Solver,		/*!< \brief Numerical solver for the implicit scheme. */
	Kind_Linear_Solver_Prec,		/*!< \brief Preconditioner of the linear solver. */
	Kind_AdjTurb_Linear_Solver,		/*!< \brief Numerical solver for the turbulent adjoint implicit scheme. */
	Kind_AdjTurb_Linear_Prec,		/*!< \brief Preconditioner of the turbulent adjoint linear solver. */
	Kind_SlopeLimit,				/*!< \brief Global slope limiter. */
	Kind_SlopeLimit_Flow,		/*!< \brief Slope limiter for flow equations.*/
	Kind_SlopeLimit_TNE2,		/*!< \brief Slope limiter for flow equations.*/
  Kind_SlopeLimit_AdjTNE2,		/*!< \brief Slope limiter for flow equations.*/
	Kind_SlopeLimit_Turb,		/*!< \brief Slope limiter for the turbulence equation.*/
	Kind_SlopeLimit_AdjLevelSet,		/*!< \brief Slope limiter for the adjoint level set equation.*/
	Kind_SlopeLimit_AdjTurb,	/*!< \brief Slope limiter for the adjoint turbulent equation.*/
	Kind_SlopeLimit_AdjFlow,	/*!< \brief Slope limiter for the adjoint equation.*/
	Kind_TimeNumScheme,			/*!< \brief Global explicit or implicit time integration. */
	Kind_TimeIntScheme_Flow,	/*!< \brief Time integration for the flow equations. */
	Kind_TimeIntScheme_AdjFlow,		/*!< \brief Time integration for the adjoint flow equations. */
  Kind_TimeIntScheme_TNE2,	/*!< \brief Time integration for the flow equations. */
  Kind_TimeIntScheme_AdjTNE2, /*!< \brief Time integration for the flow equations. */
	Kind_TimeIntScheme_LinFlow,		/*!< \brief Time integration for the linearized flow equations. */
	Kind_TimeIntScheme_Turb,	/*!< \brief Time integration for the turbulence model. */
	Kind_TimeIntScheme_AdjLevelSet,	/*!< \brief Time integration for the adjoint level set model. */
	Kind_TimeIntScheme_AdjTurb,	/*!< \brief Time integration for the adjoint turbulence model. */
	Kind_TimeIntScheme_Wave,	/*!< \brief Time integration for the wave equations. */
	Kind_TimeIntScheme_Heat,	/*!< \brief Time integration for the wave equations. */
	Kind_TimeIntScheme_Poisson,	/*!< \brief Time integration for the wave equations. */
	Kind_TimeIntScheme_FEA,	/*!< \brief Time integration for the FEA equations. */
	Kind_ConvNumScheme,			/*!< \brief Global definition of the convective term. */
	Kind_ConvNumScheme_Flow,	/*!< \brief Centered or upwind scheme for the flow equations. */
	Kind_ConvNumScheme_Heat,	/*!< \brief Centered or upwind scheme for the flow equations. */
	Kind_ConvNumScheme_TNE2,	/*!< \brief Centered or upwind scheme for the flow equations. */
	Kind_ConvNumScheme_AdjFlow,		/*!< \brief Centered or upwind scheme for the adjoint flow equations. */
  Kind_ConvNumScheme_AdjTNE2,		/*!< \brief Centered or upwind scheme for the adjoint TNE2 equations. */
	Kind_ConvNumScheme_LinFlow,		/*!< \brief Centered or upwind scheme for the linearized flow equations. */
	Kind_ConvNumScheme_Turb,	/*!< \brief Centered or upwind scheme for the turbulence model. */
	Kind_ConvNumScheme_AdjTurb,	/*!< \brief Centered or upwind scheme for the adjoint turbulence model. */
	Kind_ConvNumScheme_AdjLevelSet,	/*!< \brief Centered or upwind scheme for the adjoint level set equation. */
	Kind_ConvNumScheme_Template,	/*!< \brief Centered or upwind scheme for the level set equation. */
	Kind_Centered,				/*!< \brief Centered scheme. */
	Kind_Centered_Flow,			/*!< \brief Centered scheme for the flow equations. */
	Kind_Centered_TNE2,			/*!< \brief Centered scheme for the flow equations. */
	Kind_Centered_AdjLevelSet,			/*!< \brief Centered scheme for the level set equation. */
	Kind_Centered_AdjFlow,			/*!< \brief Centered scheme for the adjoint flow equations. */
  Kind_Centered_AdjTNE2,			/*!< \brief Centered scheme for the adjoint TNE2 equations. */
	Kind_Centered_LinFlow,			/*!< \brief Centered scheme for the linearized flow equations. */
	Kind_Centered_Turb,			/*!< \brief Centered scheme for the turbulence model. */
	Kind_Centered_AdjTurb,		/*!< \brief Centered scheme for the adjoint turbulence model. */
	Kind_Centered_Template,		/*!< \brief Centered scheme for the template model. */
	Kind_Upwind,				/*!< \brief Upwind scheme. */
	Kind_Upwind_Flow,			/*!< \brief Upwind scheme for the flow equations. */
	Kind_Upwind_TNE2,			/*!< \brief Upwind scheme for the flow equations. */
	Kind_Upwind_AdjLevelSet,			/*!< \brief Upwind scheme for the level set equations. */
	Kind_Upwind_AdjFlow,			/*!< \brief Upwind scheme for the adjoint flow equations. */
  Kind_Upwind_AdjTNE2,			/*!< \brief Upwind scheme for the adjoint TNE2 equations. */
	Kind_Upwind_LinFlow,			/*!< \brief Upwind scheme for the linearized flow equations. */
	Kind_Upwind_Turb,			/*!< \brief Upwind scheme for the turbulence model. */
	Kind_Upwind_AdjTurb,		/*!< \brief Upwind scheme for the adjoint turbulence model. */
	Kind_Upwind_Template,			/*!< \brief Upwind scheme for the template model. */
  SpatialOrder,		/*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_Flow,		/*!< \brief Order of the spatial numerical integration.*/
	SpatialOrder_Turb,		/*!< \brief Order of the spatial numerical integration.*/
	SpatialOrder_TNE2,		/*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_AdjFlow,		/*!< \brief Order of the spatial numerical integration.*/
	SpatialOrder_AdjTurb,		/*!< \brief Order of the spatial numerical integration.*/
	SpatialOrder_AdjTNE2,     /*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_AdjLevelSet;		/*!< \brief Order of the spatial numerical integration.*/

  unsigned short Kind_Turb_Model;			/*!< \brief Turbulent model definition. */
  string ML_Turb_Model_File;  /*!< \brief File containing turbulence model. */
  string ML_Turb_Model_FeatureSet; /*! <\brief What are the input and ouput features > */
  string *ML_Turb_Model_Extra; /*! <\brief Store for extra variables coming from ML turb model */
  unsigned short nML_Turb_Model_Extra; /*!<\brief number of strings there */

  unsigned short Kind_Trans_Model,			/*!< \brief Transition model definition. */
	Kind_Inlet, *Kind_Data_Riemann;           /*!< \brief Kind of inlet boundary treatment. */
	double Linear_Solver_Error;		/*!< \brief Min error of the linear solver for the implicit formulation. */
	unsigned long Linear_Solver_Iter;		/*!< \brief Max iterations of the linear solver for the implicit formulation. */
	unsigned long Linear_Solver_Restart_Frequency;   /*!< \brief Restart frequency of the linear solver for the implicit formulation. */
	double Linear_Solver_Relax;		/*!< \brief Relaxation coefficient of the linear solver. */
	double AdjTurb_Linear_Error;		/*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  double EntropyFix_Coeff;              /*!< \brief Entropy fix coefficient. */
	unsigned short AdjTurb_Linear_Iter;		/*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
	double *Section_Location;                  /*!< \brief Airfoil section limit. */
  unsigned short nSections,      /*!< \brief Number of section cuts to make when calculating internal volume. */
  nVolSections;               /*!< \brief Number of sections. */
	double* Kappa_Flow,           /*!< \brief Numerical dissipation coefficients for the flow equations. */
	*Kappa_AdjFlow,                  /*!< \brief Numerical dissipation coefficients for the adjoint equations. */
  *Kappa_TNE2,             /*!< \brief Numerical dissipation coefficients for the TNE2 equations. */
  *Kappa_AdjTNE2,          /*!< \brief Numerical dissipation coefficients for the adjoint TNE2 equations. */
	*Kappa_LinFlow;                  /*!< \brief Numerical dissipation coefficients for the linearized equations. */
	double Kappa_1st_AdjFlow,	/*!< \brief JST 1st order dissipation coefficient for adjoint flow equations (coarse multigrid levels). */
	Kappa_2nd_AdjFlow,			/*!< \brief JST 2nd order dissipation coefficient for adjoint flow equations. */
	Kappa_4th_AdjFlow,			/*!< \brief JST 4th order dissipation coefficient for adjoint flow equations. */
	Kappa_1st_LinFlow,			/*!< \brief JST 1st order dissipation coefficient for linearized flow equations (coarse multigrid levels). */
	Kappa_4th_LinFlow,			/*!< \brief JST 4th order dissipation coefficient for linearized flow equations. */
	Kappa_1st_Flow,			/*!< \brief JST 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
	Kappa_2nd_Flow,			/*!< \brief JST 2nd order dissipation coefficient for flow equations. */
	Kappa_4th_Flow,			/*!< \brief JST 4th order dissipation coefficient for flow equations. */
	Kappa_1st_TNE2,			/*!< \brief JST 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
	Kappa_2nd_TNE2,			/*!< \brief JST 2nd order dissipation coefficient for flow equations. */
	Kappa_4th_TNE2,			/*!< \brief JST 4th order dissipation coefficient for flow equations. */
  Kappa_1st_AdjTNE2,			/*!< \brief JST 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
	Kappa_2nd_AdjTNE2,			/*!< \brief JST 2nd order dissipation coefficient for flow equations. */
	Kappa_4th_AdjTNE2;			/*!< \brief JST 4th order dissipation coefficient for flow equations. */

	double Min_Beta_RoeTurkel,		/*!< \brief Minimum value of Beta for the Roe-Turkel low Mach preconditioner. */
	Max_Beta_RoeTurkel;		/*!< \brief Maximum value of Beta for the Roe-Turkel low Mach preconditioner. */
  unsigned long GridDef_Nonlinear_Iter, /*!< \brief Number of nonlinear increments for grid deformation. */
  GridDef_Linear_Iter; /*!< \brief Number of linear smoothing iterations for grid deformation. */
  unsigned short Deform_Stiffness_Type; /*!< \brief Type of element stiffness imposed for FEA mesh deformation. */
  bool Deform_Output;  /*!< \brief Print the residuals during mesh deformation to the console. */
  double Deform_Tol_Factor; /*!< Factor to multiply smallest volume for deform tolerance (0.001 default) */
  double Deform_ElasticityMod, Deform_PoissonRatio; /*!< young's modulus and poisson ratio for volume deformation stiffness model */
  bool Visualize_Deformation;	/*!< \brief Flag to visualize the deformation in MDC. */
	double Mach;		/*!< \brief Mach number. */
	double Reynolds;	/*!< \brief Reynolds number. */
	double Froude;	/*!< \brief Froude number. */
	double Length_Reynolds;	/*!< \brief Reynolds length (dimensional). */
	double AoA,			/*!< \brief Angle of attack (just external flow). */
	AoS;				/*!< \brief Angle of sideSlip (just external flow). */
  bool Fixed_CL_Mode;			/*!< \brief Activate fixed CL mode (external flow only). */
  double Target_CL;			/*!< \brief Specify a target CL instead of AoA (external flow only). */
  double Damp_Fixed_CL;			/*!< \brief Damping coefficient for fixed CL mode (external flow only). */
  bool Update_AoA;			/*!< \brief Boolean flag for whether to update the AoA for fixed lift mode on a given iteration. */
	double ChargeCoeff;		/*!< \brief Charge coefficient (just for poisson problems). */
	double *U_FreeStreamND;			/*!< \brief Reference variables at the infinity, free stream values. */
	unsigned short Cauchy_Func_Flow,	/*!< \brief Function where to apply the convergence criteria in the flow problem. */
	Cauchy_Func_AdjFlow,				/*!< \brief Function where to apply the convergence criteria in the adjoint problem. */
	Cauchy_Func_LinFlow,				/*!< \brief Function where to apply the convergence criteria in the linearized problem. */
	Cauchy_Elems;						/*!< \brief Number of elements to evaluate. */
	unsigned short Residual_Func_Flow;	/*!< \brief Equation to apply residual convergence to. */
	unsigned long StartConv_Iter;	/*!< \brief Start convergence criteria at iteration. */
	double Cauchy_Eps,	/*!< \brief Epsilon used for the convergence. */
	Cauchy_Eps_OneShot,	/*!< \brief Epsilon used for the one shot method convergence. */
	Cauchy_Eps_FullMG;	/*!< \brief Epsilon used for the full multigrid method convergence. */
	unsigned long Wrt_Sol_Freq,	/*!< \brief Writing solution frequency. */
	Wrt_Sol_Freq_DualTime,	/*!< \brief Writing solution frequency for Dual Time. */
	Wrt_Con_Freq,				/*!< \brief Writing convergence history frequency. */
	Wrt_Con_Freq_DualTime;				/*!< \brief Writing convergence history frequency. */
	bool Wrt_Unsteady;  /*!< \brief Write unsteady data adding header and prefix. */
	bool LowFidelitySim;  /*!< \brief Compute a low fidelity simulation. */
	bool Restart,	/*!< \brief Restart solution (for direct, adjoint, and linearized problems).*/
	Restart_Flow;	/*!< \brief Restart flow solution for adjoint and linearized problems. */
	unsigned short nMarker_Monitoring,	/*!< \brief Number of markers to monitor. */
	nMarker_Designing,					/*!< \brief Number of markers for the objective function. */
	nMarker_GeoEval,					/*!< \brief Number of markers for the objective function. */
	nMarker_Plotting,					/*!< \brief Number of markers to plot. */
  nMarker_Moving,               /*!< \brief Number of markers in motion (DEFORMING, MOVING_WALL, or FLUID_STRUCTURE). */
	nMarker_DV;               /*!< \brief Number of markers affected by the design variables. */
  string *Marker_Monitoring,     /*!< \brief Markers to monitor. */
  *Marker_Designing,         /*!< \brief Markers to plot. */
  *Marker_GeoEval,         /*!< \brief Markers to plot. */
  *Marker_Plotting,          /*!< \brief Markers to plot. */
  *Marker_Moving,            /*!< \brief Markers in motion (DEFORMING, MOVING_WALL, or FLUID_STRUCTURE). */
  *Marker_DV;            /*!< \brief Markers affected by the design variables. */
  unsigned short  *Marker_All_Monitoring,        /*!< \brief Global index for monitoring using the grid information. */
  *Marker_All_GeoEval,       /*!< \brief Global index for geometrical evaluation. */
  *Marker_All_Plotting,        /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_DV,          /*!< \brief Global index for design variable markers using the grid information. */
  *Marker_All_Moving,          /*!< \brief Global index for moving surfaces using the grid information. */
  *Marker_All_Designing,         /*!< \brief Global index for moving using the grid information. */
  *Marker_All_Out_1D,      /*!< \brief Global index for moving using 1D integrated output. */
  *Marker_CfgFile_Monitoring,     /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Designing,      /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_GeoEval,      /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Plotting,     /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_Out_1D,      /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_Moving,       /*!< \brief Global index for moving surfaces using the config information. */
  *Marker_CfgFile_DV,       /*!< \brief Global index for design variable markers using the config information. */
  *Marker_CfgFile_PerBound;     /*!< \brief Global index for periodic boundaries using the config information. */
  string *PlaneTag;      /*!< \brief Global index for the plane adaptation (upper, lower). */
	unsigned short nDomain;			/*!< \brief Number of domains in the MPI parallelization. */
	double DualVol_Power;			/*!< \brief Power for the dual volume in the grid adaptation sensor. */
	unsigned short Analytical_Surface;	/*!< \brief Information about the analytical definition of the surface for grid adaptation. */
	unsigned short Axis_Orientation;	/*!< \brief Axis orientation. */
	unsigned short Mesh_FileFormat;	/*!< \brief Mesh input format. */
	unsigned short Output_FileFormat;	/*!< \brief Format of the output files. */
	double RefAreaCoeff,		/*!< \brief Reference area for coefficient computation. */
	RefElemLength,				/*!< \brief Reference element length for computing the slope limiting epsilon. */
	RefSharpEdges,				/*!< \brief Reference coefficient for detecting sharp edges. */
	RefLengthMoment,			/*!< \brief Reference length for moment computation. */
  *RefOriginMoment,           /*!< \brief Origin for moment computation. */
  *RefOriginMoment_X,      /*!< \brief X Origin for moment computation. */
  *RefOriginMoment_Y,      /*!< \brief Y Origin for moment computation. */
  *RefOriginMoment_Z,      /*!< \brief Z Origin for moment computation. */
  *CFLRamp,      /*!< \brief Information about the CFL ramp. */
  *CFL,
	DomainVolume;		/*!< \brief Volume of the computational grid. */
  unsigned short nRefOriginMoment_X,    /*!< \brief Number of X-coordinate moment computation origins. */
	nRefOriginMoment_Y,           /*!< \brief Number of Y-coordinate moment computation origins. */
	nRefOriginMoment_Z;           /*!< \brief Number of Z-coordinate moment computation origins. */
	string Mesh_FileName,			/*!< \brief Mesh input file. */
	Mesh_Out_FileName,				/*!< \brief Mesh output file. */
	Solution_FlowFileName,			/*!< \brief Flow solution input file. */
	Solution_LinFileName,			/*!< \brief Linearized flow solution input file. */
	Solution_AdjFileName,			/*!< \brief Adjoint solution input file for drag functional. */
	Flow_FileName,					/*!< \brief Flow variables output file. */
	Structure_FileName,					/*!< \brief Structure variables output file. */
	SurfStructure_FileName,					/*!< \brief Surface structure variables output file. */
  SurfWave_FileName,					/*!< \brief Surface structure variables output file. */
	SurfHeat_FileName,					/*!< \brief Surface structure variables output file. */
	Wave_FileName,					/*!< \brief Wave variables output file. */
	Heat_FileName,					/*!< \brief Heat variables output file. */
	AdjWave_FileName,					/*!< \brief Adjoint wave variables output file. */
	Residual_FileName,				/*!< \brief Residual variables output file. */
	Conv_FileName,					/*!< \brief Convergence history output file. */
	Restart_FlowFileName,			/*!< \brief Restart file for flow variables. */
	Restart_WaveFileName,			/*!< \brief Restart file for wave variables. */
	Restart_HeatFileName,			/*!< \brief Restart file for heat variables. */
	Restart_LinFileName,			/*!< \brief Restart file for linearized flow variables. */
	Restart_AdjFileName,			/*!< \brief Restart file for adjoint variables, drag functional. */
	Adj_FileName,					/*!< \brief Output file with the adjoint variables. */
	Lin_FileName,					/*!< \brief Output file with the linearized variables. */
	ObjFunc_Grad_FileName,			/*!< \brief Gradient of the objective function. */
	ObjFunc_Value_FileName,			/*!< \brief Objective function. */
	SurfFlowCoeff_FileName,			/*!< \brief Output file with the flow variables on the surface. */
	SurfAdjCoeff_FileName,			/*!< \brief Output file with the adjoint variables on the surface. */
	SurfLinCoeff_FileName,			/*!< \brief Output file with the linearized variables on the surface. */
	New_SU2_FileName;        		/*!< \brief Output SU2 mesh file converted from CGNS format. */
	bool CGNS_To_SU2;      		 	/*!< \brief Flag to specify whether a CGNS mesh is converted to SU2 format. */
	unsigned short nSpecies, 		/*!< \brief No of species present in plasma */
	nReactions;									/*!< \brief Number of reactions in chemical model. */
	bool Low_MemoryOutput,      /*!< \brief Write a volume solution file */
  Wrt_Vol_Sol,                /*!< \brief Write a volume solution file */
	Wrt_Srf_Sol,                /*!< \brief Write a surface solution file */
	Wrt_Restart,                /*!< \brief Write a restart solution file */
	Wrt_Csv_Sol,                /*!< \brief Write a surface comma-separated values solution file */
	Wrt_Residuals,              /*!< \brief Write residuals to solution file */
  Wrt_Limiters,              /*!< \brief Write residuals to solution file */
	Wrt_SharpEdges,              /*!< \brief Write residuals to solution file */
  Wrt_Halo,                   /*!< \brief Write rind layers in solution files */
  Plot_Section_Forces,       /*!< \brief Write sectional forces for specified markers. */
	Wrt_1D_Output;                /*!< \brief Write average stagnation pressure specified markers. */
	double *ArrheniusCoefficient,					/*!< \brief Arrhenius reaction coefficient */
	*ArrheniusEta,								/*!< \brief Arrhenius reaction temperature exponent */
	*ArrheniusTheta,							/*!< \brief Arrhenius reaction characteristic temperature */
	*CharVibTemp,									/*!< \brief Characteristic vibrational temperature for e_vib */
  *RotationModes,				/*!< \brief Rotational modes of energy storage */
  *Ref_Temperature,   			/*!< \brief Reference temperature for thermodynamic relations */
  *Tcf_a,   /*!< \brief Rate controlling temperature exponent (fwd) */
  *Tcf_b,   /*!< \brief Rate controlling temperature exponent (fwd) */
  *Tcb_a,   /*!< \brief Rate controlling temperature exponent (bkw) */
  *Tcb_b,   /*!< \brief Rate controlling temperature exponent (bkw) */
  *Diss;                /*!< \brief Dissociation potential. */
	unsigned short nMass,                 /*!< \brief No of particle masses */
	nTemp;						/*!< \brief No of freestream temperatures specified */
	bool Inlet_Outlet_Defined; /*!< \brief  that inlet and outlet conditions are defined for each species*/
  double *Particle_Mass,         /*!< \brief Mass of all particles present in the plasma */
  *Molar_Mass,               /*!< \brief Molar mass of species in the plasma [kg/kmol] */
  Mixture_Molar_mass,       /*!< \brief Molar mass of the multi-species fluid [kg/kmol] */
  *Gas_Composition,          /*!< \brief Initial mass fractions of flow [dimensionless] */
  *Enthalpy_Formation,     /*!< \brief Enthalpy of formation */
  **Blottner,               /*!< \brief Blottner viscosity coefficients */
  *Species_Ref_Temperature,  /*!< \brief Reference Temperature for viscosity of all particles present in the plasma */
  *Species_Ref_Viscosity;    /*!< \brief Reference viscosity  of all particles present in the plasma */
  unsigned short *nElStates; /*!< \brief Number of electron states. */
  double **CharElTemp, /*!< \brief Characteristic temperature of electron states. */
  **degen; /*!< \brief Degeneracy of electron states. */
	double Gamma,			/*!< \brief Ratio of specific heats of the gas. */
	Bulk_Modulus,			/*!< \brief Value of the bulk modulus for incompressible flows. */
	ArtComp_Factor,			/*!< \brief Value of the artificial compresibility factor for incompressible flows. */
	Gas_Constant,     /*!< \brief Specific gas constant. */
	Gas_ConstantND,     /*!< \brief Non-dimensional specific gas constant. */
	Gas_Constant_Ref, /*!< \brief Reference specific gas constant. */
	Temperature_Critical,   /*!< \brief Critical Temperature for real fluid model.  */
	Pressure_Critical,   /*!< \brief Critical Pressure for real fluid model.  */
	Density_Critical,   /*!< \brief Critical Density for real fluid model.  */
	Acentric_Factor,   /*!< \brief Acentric Factor for real fluid model.  */
	Mu_ConstantND,   /*!< \brief Constant Viscosity for CostantViscosity model.  */
	Kt_ConstantND,   /*!< \brief Constant Thermal Conductivity for CostantConductivity model.  */
	Mu_RefND,   /*!< \brief reference viscosity for Sutherland model.  */
	Mu_Temperature_RefND,   /*!< \brief reference Temperature for Sutherland model.  */
	Mu_SND,   /*!< \brief reference S for Sutherland model.  */
	FreeSurface_Zero,	/*!< \brief Coordinate of the level set zero. */
	FreeSurface_Depth,	/*!< \brief Coordinate of the level set zero. */
	*Velocity_FreeStream,     /*!< \brief Total velocity of the fluid.  */
	Energy_FreeStream,     /*!< \brief Total energy of the fluid.  */
	ModVel_FreeStream,     /*!< \brief Total density of the fluid.  */
	ModVel_FreeStreamND,     /*!< \brief Total density of the fluid.  */
	Density_FreeStream,     /*!< \brief Total density of the fluid. */
	Viscosity_FreeStream,     /*!< \brief Total density of the fluid.  */
	Tke_FreeStream,     /*!< \brief Total turbulent kinetic energy of the fluid.  */
	Intermittency_FreeStream,     /*!< \brief Freestream intermittency (for sagt transition model) of the fluid.  */
	TurbulenceIntensity_FreeStream,     /*!< \brief Freestream turbulent intensity (for sagt transition model) of the fluid.  */
	Turb2LamViscRatio_FreeStream,          /*!< \brief Ratio of turbulent to laminar viscosity. */
	NuFactor_FreeStream,  /*!< \brief Ratio of turbulent to laminar viscosity. */
	Pressure_FreeStream,     /*!< \brief Total pressure of the fluid. */
	Temperature_FreeStream,  /*!< \brief Total temperature of the fluid.  */
  Temperature_ve_FreeStream,  /*!< \brief Total vibrational-electronic temperature of the fluid.  */
  *MassFrac_FreeStream, /*!< \brief Mixture mass fractions of the fluid. */
	Prandtl_Lam,      /*!< \brief Laminar Prandtl number for the gas.  */
	Prandtl_Turb,     /*!< \brief Turbulent Prandtl number for the gas.  */
	Length_Ref,       /*!< \brief Reference length for non-dimensionalization. */
	Mesh_Scale_Change,       /*!< \brief Conversion factor from grid units to meters. */
	Pressure_Ref,     /*!< \brief Reference pressure for non-dimensionalization.  */
	Temperature_Ref,  /*!< \brief Reference temperature for non-dimensionalization.*/
	Density_Ref,      /*!< \brief Reference density for non-dimensionalization.*/
	Velocity_Ref,     /*!< \brief Reference velocity for non-dimensionalization.*/
	Time_Ref,         /*!< \brief Reference time for non-dimensionalization. */
	Viscosity_Ref,    /*!< \brief Reference viscosity for non-dimensionalization. */
	Conductivity_Ref,    /*!< \brief Reference conductivity for non-dimensionalization. */
	Energy_Ref,    /*!< \brief Reference viscosity for non-dimensionalization. */
	Wall_Temperature,    /*!< \brief Temperature at an isotropic wall in Kelvin. */
	Omega_Ref,        /*!< \brief Reference angular velocity for non-dimensionalization. */
	Force_Ref,        /*!< \brief Reference body force for non-dimensionalization. */
	Pressure_FreeStreamND,     /*!< \brief Farfield pressure value (external flow). */
	Temperature_FreeStreamND,  /*!< \brief Farfield temperature value (external flow). */
	Density_FreeStreamND,      /*!< \brief Farfield density value (external flow). */
  Velocity_FreeStreamND[3],    /*!< \brief Farfield velocity values (external flow). */
	Energy_FreeStreamND,       /*!< \brief Farfield energy value (external flow). */
	Viscosity_FreeStreamND,    /*!< \brief Farfield viscosity value (external flow). */
	Tke_FreeStreamND,    /*!< \brief Farfield kinetic energy (external flow). */
  Omega_FreeStreamND, /*!< \brief Specific dissipation (external flow). */
  Omega_FreeStream, /*!< \brief Specific dissipation (external flow). */
  pnorm_heat;           /*!< \brief pnorm for heat-flux objective functions. */
	int ***Reactions;					/*!< \brief Reaction map for chemically reacting, multi-species flows. */
  double ***Omega00,        /*!< \brief Collision integrals (Omega(0,0)) */
  ***Omega11;                  /*!< \brief Collision integrals (Omega(1,1)) */
	bool Mesh_Output; /*!< \brief Flag to specify whether a new mesh should be written in the converted units. */
	double ElasticyMod,			/*!< \brief Young's modulus of elasticity. */
	PoissonRatio,						/*!< \brief Poisson's ratio. */
	MaterialDensity;								/*!< \brief Material density. */
	double Wave_Speed;			/*!< \brief Wave speed used in the wave solver. */
	double Thermal_Diffusivity;			/*!< \brief Thermal diffusivity used in the heat solver. */
	double Cyclic_Pitch,          /*!< \brief Cyclic pitch for rotorcraft simulations. */
	Collective_Pitch;             /*!< \brief Collective pitch for rotorcraft simulations. */
	string Motion_Filename;				/*!< \brief Arbitrary mesh motion input base filename. */
	double Mach_Motion;			/*!< \brief Mach number based on mesh velocity and freestream quantities. */
  double *Motion_Origin_X,    /*!< \brief X-coordinate of the mesh motion origin. */
  *Motion_Origin_Y,           /*!< \brief Y-coordinate of the mesh motion origin. */
  *Motion_Origin_Z,           /*!< \brief Z-coordinate of the mesh motion origin. */
  *Translation_Rate_X,           /*!< \brief Translational velocity of the mesh in the x-direction. */
  *Translation_Rate_Y,           /*!< \brief Translational velocity of the mesh in the y-direction. */
  *Translation_Rate_Z,           /*!< \brief Translational velocity of the mesh in the z-direction. */
  *Rotation_Rate_X,           /*!< \brief Angular velocity of the mesh about the x-axis. */
  *Rotation_Rate_Y,           /*!< \brief Angular velocity of the mesh about the y-axis. */
  *Rotation_Rate_Z,           /*!< \brief Angular velocity of the mesh about the z-axis. */
  *Pitching_Omega_X,           /*!< \brief Angular frequency of the mesh pitching about the x-axis. */
  *Pitching_Omega_Y,           /*!< \brief Angular frequency of the mesh pitching about the y-axis. */
  *Pitching_Omega_Z,           /*!< \brief Angular frequency of the mesh pitching about the z-axis. */
  *Pitching_Ampl_X,           /*!< \brief Pitching amplitude about the x-axis. */
  *Pitching_Ampl_Y,           /*!< \brief Pitching amplitude about the y-axis. */
  *Pitching_Ampl_Z,           /*!< \brief Pitching amplitude about the z-axis. */
  *Pitching_Phase_X,           /*!< \brief Pitching phase offset about the x-axis. */
  *Pitching_Phase_Y,           /*!< \brief Pitching phase offset about the y-axis. */
  *Pitching_Phase_Z,           /*!< \brief Pitching phase offset about the z-axis. */
  *Plunging_Omega_X,           /*!< \brief Angular frequency of the mesh plunging in the x-direction. */
  *Plunging_Omega_Y,           /*!< \brief Angular frequency of the mesh plunging in the y-direction. */
  *Plunging_Omega_Z,           /*!< \brief Angular frequency of the mesh plunging in the z-direction. */
  *Plunging_Ampl_X,           /*!< \brief Plunging amplitude in the x-direction. */
  *Plunging_Ampl_Y,           /*!< \brief Plunging amplitude in the y-direction. */
  *Plunging_Ampl_Z;           /*!< \brief Plunging amplitude in the z-direction. */
  unsigned short nMotion_Origin_X,    /*!< \brief Number of X-coordinate mesh motion origins. */
	nMotion_Origin_Y,           /*!< \brief Number of Y-coordinate mesh motion origins. */
	nMotion_Origin_Z,           /*!< \brief Number of Z-coordinate mesh motion origins. */
	nTranslation_Rate_X,           /*!< \brief Number of Translational x-velocities for mesh motion. */
	nTranslation_Rate_Y,           /*!< \brief Number of Translational y-velocities for mesh motion. */
	nTranslation_Rate_Z,           /*!< \brief Number of Translational z-velocities for mesh motion. */
	nRotation_Rate_X,           /*!< \brief Number of Angular velocities about the x-axis for mesh motion. */
	nRotation_Rate_Y,           /*!< \brief Number of Angular velocities about the y-axis for mesh motion. */
	nRotation_Rate_Z,           /*!< \brief Number of Angular velocities about the z-axis for mesh motion. */
	nPitching_Omega_X,           /*!< \brief Number of Angular frequencies about the x-axis for pitching. */
	nPitching_Omega_Y,           /*!< \brief Number of Angular frequencies about the y-axis for pitching. */
	nPitching_Omega_Z,           /*!< \brief Number of Angular frequencies about the z-axis for pitching. */
	nPitching_Ampl_X,           /*!< \brief Number of Pitching amplitudes about the x-axis. */
	nPitching_Ampl_Y,           /*!< \brief Number of Pitching amplitudes about the y-axis. */
	nPitching_Ampl_Z,           /*!< \brief Number of Pitching amplitudes about the z-axis. */
	nPitching_Phase_X,           /*!< \brief Number of Pitching phase offsets about the x-axis. */
	nPitching_Phase_Y,           /*!< \brief Number of Pitching phase offsets about the y-axis. */
	nPitching_Phase_Z,           /*!< \brief Number of Pitching phase offsets about the z-axis. */
	nPlunging_Omega_X,           /*!< \brief Number of Angular frequencies in the x-direction for plunging. */
	nPlunging_Omega_Y,           /*!< \brief Number of Angular frequencies in the y-direction for plunging. */
	nPlunging_Omega_Z,           /*!< \brief Number of Angular frequencies in the z-direction for plunging. */
	nPlunging_Ampl_X,           /*!< \brief Number of Plunging amplitudes in the x-direction. */
	nPlunging_Ampl_Y,           /*!< \brief Number of Plunging amplitudes in the y-direction. */
	nPlunging_Ampl_Z,           /*!< \brief Number of Plunging amplitudes in the z-direction. */
  nMoveMotion_Origin,         /*!< \brief Number of motion origins. */
  *MoveMotion_Origin;         /*!< \brief Keeps track if we should move moment origin. */
  vector<vector<vector<double> > > Aeroelastic_np1, /*!< \brief Aeroelastic solution at time level n+1. */
  Aeroelastic_n, /*!< \brief Aeroelastic solution at time level n. */
	Aeroelastic_n1; /*!< \brief Aeroelastic solution at time level n-1. */
  double FreqPlungeAeroelastic, /*!< \brief Plunging natural frequency for Aeroelastic. */
	FreqPitchAeroelastic; /*!< \brief Pitch natural frequency for Aeroelastic. */
  double *Aeroelastic_plunge, /*!< \brief Value of plunging coordinate at the end of an external iteration. */
	*Aeroelastic_pitch; /*!< \brief Value of pitching coordinate at the end of an external iteration. */
  unsigned short Gust_Type,	/*!< \brief Type of Gust. */
  Gust_Dir;   /*!< \brief Direction of the gust */
  double Gust_WaveLength,     /*!< \brief The gust wavelength. */
  Gust_Periods,              /*!< \brief Number of gust periods. */
  Gust_Ampl,                  /*!< \brief Gust amplitude. */
  Gust_Begin_Time,            /*!< \brief Time at which to begin the gust. */
  Gust_Begin_Loc;             /*!< \brief Location at which the gust begins. */
  long Visualize_CV; /*!< \brief Node number for the CV to be visualized */
  bool ExtraOutput;
  unsigned long Nonphys_Points, /*!< \brief Current number of non-physical points in the solution. */
  Nonphys_Reconstr;      /*!< \brief Current number of non-physical reconstructions for 2nd-order upwinding. */
  
  /*!< \brief param is a map from the option name (config file string) to a pointer to an option child class */
//	map<string, CAnyOptionRef*> param;

  /*!<brief all_options is a map containing all of the options. This is used during config file parsing
  to track the options which have not been set (so the default values can be used). Without this map
   there would be no list of all the config file options. > */
  map<string, bool> all_options;

  /*<brief param is a map from the option name (config file string) to its decoder (the specific child
   class of COptionBase that turns the string into a value) */
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
   value can be represented by a double.*/
  void addDoubleOption(const string name, double & option_field, double default_value){
    // Check if the key is already in the map. If this fails, it is coder error
    // and not user error, so throw.
    assert(option_map.find(name) == option_map.end());

    // Add this option to the list of all the options
    all_options.insert(pair<string,bool>(name,true));

    // Create the parser for a double option with a reference to the option_field and the desired
    // default value. This will take the string in the config file, convert it to a double, and
    // place that double in the memory location specified by the reference.
    COptionBase* val = new COptionDouble(name, option_field, default_value);

    // Create an association between the option name ("CFL") and the parser generated above.
    // During configuration, the parsing script will get the option name, and use this map
    // to find how to parse that option.
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addStringOption(const string name, string & option_field, string default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionString(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addIntegerOption(const string name, int & option_field, int default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionInt(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addUnsignedLongOption(const string name, unsigned long & option_field, unsigned long default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionULong(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addUnsignedShortOption(const string name, unsigned short & option_field, unsigned short default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionUShort(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addLongOption(const string name, long & option_field, long default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionLong(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addBoolOption(const string name, bool & option_field, bool default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionBool(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  // enum types work differently than all of the others because there are a small number of valid
  // string entries for the type. One must also provide a list of all the valid strings of that type.
  template <class Tenum>
  void addEnumOption(const string name, unsigned short & option_field, const map<string, Tenum> & enum_map, Tenum default_value){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionEnum<Tenum>(name, enum_map, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
    return;
  }


  // input_size is the number of options read in from the config file
  template <class Tenum>
	void addEnumListOption(const string name, unsigned short & input_size, unsigned short * & option_field, const map<string, Tenum> & enum_map) {
    input_size = 0;
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
		COptionBase* val = new COptionEnumList<Tenum>(name, enum_map, option_field, input_size);
    option_map.insert( pair<string, COptionBase*>(name, val) );
	}

  void addDoubleArrayOption(const string name, const int size, double * & option_field, double * default_value){

    double * def = new double [size];
    for (int i = 0; i < size; i++){
      def[i] = default_value[i];
    }
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionDoubleArray(name, size, option_field, def);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addDoubleListOption(const string name, unsigned short & size, double * & option_field){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionDoubleList(name, size, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addUShortListOption(const string name, unsigned short & size, unsigned short * & option_field){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionUShortList(name, size, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addStringListOption(const string name, unsigned short & num_marker, string* & option_field){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionStringList(name, num_marker, option_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addConvectOption(const string name, unsigned short & space_field, unsigned short & centered_field, unsigned short & upwind_field){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionConvect(name, space_field, centered_field, upwind_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addMathProblemOption(const string name, bool & Adjoint, const bool & Adjoint_default,
                      bool & OneShot, const bool & OneShot_default,
                      bool & Linearized, const bool & Linearized_default,
                            bool & Restart_Flow, const bool & Restart_Flow_default){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionMathProblem(name, Adjoint, Adjoint_default, OneShot, OneShot_default, Linearized, Linearized_default, Restart_Flow, Restart_Flow_default);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addDVParamOption(const string name, unsigned short & nDV_field, double** & paramDV, string* & FFDTag,
                        unsigned short* & design_variable){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionDVParam(name, nDV_field, paramDV, FFDTag, design_variable);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addStringDoubleListOption(const string name, unsigned short & list_size, string * & string_field,
                        double* & double_field){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionStringDoubleList(name, list_size, string_field, double_field);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addInletOption(const string name, unsigned short & nMarker_Inlet, string * & Marker_Inlet,
                                 double* & Ttotal, double* & Ptotal, double** & FlowDir){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionInlet(name, nMarker_Inlet, Marker_Inlet, Ttotal, Ptotal, FlowDir);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  template <class Tenum>
  
  void addRiemannOption(const string name, unsigned short & nMarker_Riemann, string * & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                                 double* & var1, double* & var2, double** & FlowDir){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionRiemann<Tenum>(name, nMarker_Riemann, Marker_Riemann, option_field, enum_map, var1, var2, FlowDir);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addExhaustOption(const string name, unsigned short & nMarker_Exhaust, string * & Marker_Exhaust,
                      double* & Ttotal, double* & Ptotal){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionExhaust(name, nMarker_Exhaust, Marker_Exhaust, Ttotal, Ptotal);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addBleedOption(const string name, unsigned short & nMarker_Bleed, string * & Marker_Bleed,
                        double* & MassFlow_Target, double* & Temp_Target){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionBleed(name, nMarker_Bleed, Marker_Bleed, MassFlow_Target, Temp_Target);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addPeriodicOption(const string & name, unsigned short & nMarker_PerBound,
                    string* & Marker_PerBound, string* & Marker_PerDonor,
                         double** & RotCenter, double** & RotAngles, double** & Translation){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionPeriodic(name, nMarker_PerBound, Marker_PerBound, Marker_PerDonor, RotCenter, RotAngles, Translation);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addActuatorDiskOption(const string & name, unsigned short & nMarker_ActDisk_Inlet, unsigned short & nMarker_ActDisk_Outlet,
                                      string* & Marker_ActDisk_Inlet, string* & Marker_ActDisk_Outlet,
                                      double** & ActDisk_Origin, double* & ActDisk_RootRadius, double* & ActDisk_TipRadius,
                                      double* & ActDisk_CT, double* & ActDisk_Omega) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionActuatorDisk(name, nMarker_ActDisk_Inlet, nMarker_ActDisk_Outlet, Marker_ActDisk_Inlet, Marker_ActDisk_Outlet, ActDisk_Origin, ActDisk_RootRadius, ActDisk_TipRadius, ActDisk_CT, ActDisk_Omega);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

  void addPythonOption(const string name){
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string,bool>(name,true));
    COptionBase* val = new COptionPython(name);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }

public:

	vector<string> fields; /*!< \brief Tags for the different fields in a restart file. */

	/*!
	 * \brief Constructor of the class which reads the input file.
	 */
	CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_nDim, unsigned short verb_level);

	/*!
	 * \brief Constructor of the class which reads the input file.
	 */
	CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software);

	/*!
	 * \brief Destructor of the class.
	 */
	~CConfig(void);

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
	 * \brief Get information about whether this is a Python config option for design.
	 * \return <code>TRUE</code> if this is a Python config option for design; otherwise <code>FALSE</code>.
	 */
//	bool GetPython_Option(string & option_name);

	/*!
	 * \brief Get reference origin for moment computation.
     * \param[in] val_marker - the marker we are monitoring.
	 * \return Reference origin (in cartesians coordinates) for moment computation.
	 */
	double *GetRefOriginMoment(unsigned short val_marker);

  /*!
	 * \brief Get reference origin x-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
	 * \return Reference origin x-coordinate (in cartesians coordinates) for moment computation.
	 */
	double GetRefOriginMoment_X(unsigned short val_marker);

  /*!
	 * \brief Get reference origin y-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
	 * \return Reference origin y-coordinate (in cartesians coordinates) for moment computation.
	 */
	double GetRefOriginMoment_Y(unsigned short val_marker);

  /*!
	 * \brief Get reference origin z-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
	 * \return Reference origin z-coordinate (in cartesians coordinates) for moment computation.
	 */
	double GetRefOriginMoment_Z(unsigned short val_marker);

  /*!
	 * \brief Set reference origin x-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
	 * \param[in] val_origin - New x-coordinate of the mesh motion origin.
	 */
	void SetRefOriginMoment_X(unsigned short val_marker, double val_origin);

  /*!
	 * \brief Set reference origin y-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
	 * \param[in] val_origin - New y-coordinate of the mesh motion origin.
	 */
	void SetRefOriginMoment_Y(unsigned short val_marker, double val_origin);

  /*!
	 * \brief Set reference origin z-coordinate for moment computation.
   * \param[in] val_marker - the marker we are monitoring.
	 * \param[in] val_origin - New z-coordinate of the mesh motion origin.
	 */
	void SetRefOriginMoment_Z(unsigned short val_marker, double val_origin);

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
	double GetEA_IntLimit(unsigned short index);
  
  /*!
	 * \brief Get the integration limits for the equivalent area computation.
	 * \param[in] index - 0 means x_min, and 1 means x_max.
	 * \return Integration limits for the equivalent area computation.
	 */
	double GetEA_ScaleFactor(void);

  /*!
	 * \brief Get the limit value for the adjoint variables.
	 * \return Limit value for the adjoint variables.
	 */
	double GetAdjointLimit(void);

	/*!
	 * \brief Get the the coordinates where of the box where the grid is going to be deformed.
	 * \return Coordinates where of the box where the grid is going to be deformed.
	 */
	double *GetHold_GridFixed_Coord(void);

  /*!
   * \brief Get the the coordinates where of the box where a subsonic region is imposed.
   * \return Coordinates where of the box where the grid is going to be a subsonic region.
   */
  double *GetSubsonic_Engine_Box(void);
  
	/*!
	 * \brief Get the power of the dual volume in the grid adaptation sensor.
	 * \return Power of the dual volume in the grid adaptation sensor.
	 */
	double GetDualVol_Power(void);

	/*!
	 * \brief Get Information about if there is an analytical definition of the surface for doing the
	 *        grid adaptation.
	 * \return Definition of the surfaces. NONE implies that there isn't any analytical definition
	 *         and it will use and interpolation.
	 */
	unsigned short GetAnalytical_Surface(void);

  /*!
	 * \brief Get Information about if there is an analytical definition of the surface for doing the
	 *        grid adaptation.
	 * \return Definition of the surfaces. NONE implies that there isn't any analytical definition
	 *         and it will use and interpolation.
	 */
	unsigned short GetAxis_Orientation(void);

	/*!
	 * \brief Get the ratio of density for a free surface problem.
	 * \return Ratio of density for a free surface problem.
	 */
	double GetRatioDensity(void);

	/*!
	 * \brief Get the ratio of viscosity for a free surface problem.
	 * \return Ratio of viscosity for a free surface problem.
	 */
	double GetRatioViscosity(void);

	/*!
	 * \brief Get the thickness of the interfase for a free surface problem.
	 * \return Thickness of the interfase for a free surface problem.
	 */
	double GetFreeSurface_Thickness(void);

	/*!
	 * \brief Get the damping of the free surface for a free surface problem.
	 * \return Damping of the interfase for a free surface problem.
	 */
	double GetFreeSurface_Damping_Coeff(void);

	/*!
	 * \brief Get the damping of the free surface for a free surface problem.
	 * \return Damping of the interfase for a free surface problem.
	 */
	double GetFreeSurface_Damping_Length(void);

	/*!
	 * \brief Get the outlet position of the free surface for a free surface problem.
	 * \return Outlet position of the interfase for a free surface problem.
	 */
	double GetFreeSurface_Outlet(void);

	/*!
	 * \brief Creates a tecplot file to visualize the partition made by the DDC software.
	 * \return <code>TRUE</code> if the partition is going to be plotted; otherwise <code>FALSE</code>.
	 */
	bool GetVisualize_Partition(void);

  /*!
	 * \brief Creates a tecplot file to visualize the partition made by the DDC software.
	 * \return <code>TRUE</code> if the partition is going to be plotted; otherwise <code>FALSE</code>.
	 */
  bool GetExtraOutput(void);

	/*!
	 * \brief Get the value of the Mach number (velocity divided by speed of sound).
	 * \return Value of the Mach number.
	 */
	double GetMach(void);

	/*!
	 * \brief Get the value of the Gamma of fluid (ratio of specific heats).
	 * \return Value of the constant: Gamma
	 */
	double GetGamma(void);

	/*!
	 * \brief Get the value of the Gamma of fluid (ratio of specific heats) for a particular species.
	 * \param[in] - val_Species: Index of desired species specific heat ratio.
	 * \return Value of the constant: Species_Gamma[iSpecies]
	 */
	double GetSpecies_Gamma(unsigned short val_Species);

	/*!
	 * \brief Get the value of the charge number for a particular species (1 for ions, -1 for electrons, 0 for neutral).
	 * \param[in] - val_Species: Index of desired species charge number.
	 * \return Value of the constant: Charge_Number[val_Species]
	 */
	int GetCharge_Number(unsigned short val_Species);

  /*!
	 * \brief Get the value of the limits for the sections.
	 * \return Value of the limits for the sections.
	 */
	double GetSection_Location(unsigned short val_var);

	/*!
	 * \brief Get the array that maps chemical consituents to each chemical reaction.
	 * \return Memory location of the triple pointer to the 3-D reaction map array.
	 */
	int ***GetReaction_Map(void);

	/*!
	 * \brief Get the array containing the curve fit coefficients for the Omega(0,0) collision integrals.
	 * \return Memory location of the triple pointer to the 3-D collision integral array.
	 */
	double ***GetCollisionIntegral00(void);

	/*!
	 * \brief Get the array containing the curve fit coefficients for the Omega(1,1) collision integrals.
	 * \return Memory location of the triple pointer to the 3-D collision integral array.
	 */
	double ***GetCollisionIntegral11(void);

	/*!
	 * \brief Get the value of the bulk modulus.
	 * \return Value of the bulk modulus.
	 */
	double GetBulk_Modulus(void);

	/*!
	 * \brief Get the value of the Gamma of fluid (ratio of specific heats) for monatomic species.
	 * \return Value of the constant: GammaMonatomic
	 */
	double GetGammaMonatomic(void);

	/*!
	 * \brief Get the value of the Gamma of fluid (ratio of specific heats) for diatomic species.
	 * \return Value of the constant: Gamma
	 */
	double GetGammaDiatomic(void);

	/*!
	 * \brief Get the artificial compresibility factor.
	 * \return Value of the artificial compresibility factor.
	 */
	double GetArtComp_Factor(void);

	/*!
	 * \brief Get the Level set zero for free surface .
	 * \return Value of the level set zero coordinate
	 */
	double GetFreeSurface_Zero(void);

	/*!
	 * \brief Get the Level set zero for free surface .
	 * \return Value of the level set zero coordinate
	 */
	double GetFreeSurface_Depth(void);

	/*!
	 * \brief Get the value of specific gas constant.
	 * \return Value of the constant: Gamma
	 */
	double GetGas_Constant(void);

  /*!
	 * \brief Get the value of specific gas constant.
	 * \return Value of the constant: Gamma
	 */
	double GetGas_ConstantND(void);

	/*!
	 * \brief Get the value of specific gas constant for a particular species.
	 * \param[in] val_Species - Index of desired species gas constant.
	 * \return Value of the constant: R
	 */
	double GetSpecies_Gas_Constant(unsigned short val_Species);

	/*!
	 * \brief Get the coefficients of the Blottner viscosity model
	 * \param[in] val_Species - Index of the species
	 * \param[in] val_Coeff - Index of the coefficient (As, Bs, Cs)
	 * \return Value of the Blottner coefficient
	 */
	double GetBlottnerCoeff(unsigned short val_Species, unsigned short val_Coeff);

  /*!
	 * \brief Get the p-norm for heat-flux objective functions (adjoint problem).
	 * \return Value of the heat flux p-norm
	 */
	double GetPnormHeat(void);

	/*!
	 * \brief Get the value of wall temperature.
	 * \return Value of the constant: Temperature
	 */
	double GetWallTemperature(void);

	/*!
	 * \brief Get the reference value for the specific gas constant.
	 * \return Reference value for the specific gas constant.
	 */
	double GetGas_Constant_Ref(void);

	/*!
	 * \brief Get the value of the frestream temperature.
	 * \return Freestream temperature.
	 */
	double GetTemperature_FreeStream(void);

  /*!
	 * \brief Get the value of the frestream temperature.
	 * \return Freestream temperature.
	 */
	double GetEnergy_FreeStream(void);

  /*!
	 * \brief Get the value of the frestream temperature.
	 * \return Freestream temperature.
	 */
	double GetViscosity_FreeStream(void);

  /*!
	 * \brief Get the value of the frestream temperature.
	 * \return Freestream temperature.
	 */
	double GetDensity_FreeStream(void);

  /*!
	 * \brief Get the value of the frestream temperature.
	 * \return Freestream temperature.
	 */
	double GetModVel_FreeStream(void);

  /*!
	 * \brief Get the value of the frestream temperature.
	 * \return Freestream temperature.
	 */
	double GetModVel_FreeStreamND(void);

  /*!
	 * \brief Get the value of the frestream vibrational-electronic temperature.
	 * \return Freestream temperature.
	 */
	double GetTemperature_ve_FreeStream(void);

	/*!
	 * \brief Get the value of the laminar Prandtl number.
	 * \return Laminar Prandtl number.
	 */
	double GetPrandtl_Lam(void);

	/*!
	 * \brief Get the value of the turbulent Prandtl number.
	 * \return Turbulent Prandtl number.
	 */
	double GetPrandtl_Turb(void);

	/*!
	 * \brief Get the value of the reference length for non-dimensionalization.
	 *        This value should always be 1 internally, and is not user-specified.
	 * \return Reference length for non-dimensionalization.
	 */
	double GetLength_Ref(void);

	/*!
	 * \brief Get the value of the reference pressure for non-dimensionalization.
	 * \return Reference pressure for non-dimensionalization.
	 */
	double GetPressure_Ref(void);

  /*!
	 * \brief Get the value of the reference pressure for non-dimensionalization.
	 * \return Reference pressure for non-dimensionalization.
	 */
	double GetEnergy_Ref(void);

	/*!
	 * \brief Get the value of the reference temperature for non-dimensionalization.
	 * \return Reference temperature for non-dimensionalization.
	 */
	double GetTemperature_Ref(void);

	/*!
	 * \brief Get the value of the reference density for non-dimensionalization.
	 * \return Reference density for non-dimensionalization.
	 */
	double GetDensity_Ref(void);

	/*!
	 * \brief Get the value of the reference velocity for non-dimensionalization.
	 * \return Reference velocity for non-dimensionalization.
	 */
	double GetVelocity_Ref(void);

	/*!
	 * \brief Get the value of the reference time for non-dimensionalization.
	 * \return Reference time for non-dimensionalization.
	 */
	double GetTime_Ref(void);

	/*!
	 * \brief Get the value of the reference viscosity for non-dimensionalization.
	 * \return Reference viscosity for non-dimensionalization.
	 */
	double GetViscosity_Ref(void);

	/*!
	 * \brief Get the value of the reference conductivity for non-dimensionalization.
	 * \return Reference conductivity for non-dimensionalization.
	 */
	double GetConductivity_Ref(void);

	/*!
	 * \brief Get the value of the reference angular velocity for non-dimensionalization.
	 * \return Reference angular velocity for non-dimensionalization.
	 */
	double GetOmega_Ref(void);

	/*!
	 * \brief Get the value of the reference force for non-dimensionalization.
	 * \return Reference force for non-dimensionalization.
	 */
	double GetForce_Ref(void);

  /*!
	 * \brief Get the value of the non-dimensionalized freestream pressure.
	 * \return Non-dimensionalized freestream pressure.
	 */
	double GetPressure_FreeStream(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream pressure.
	 * \return Non-dimensionalized freestream pressure.
	 */
	double GetPressure_FreeStreamND(void);

	/*!
	 * \brief Get the vector of the dimensionalized freestream velocity.
	 * \return Dimensionalized freestream velocity vector.
	 */
	double* GetVelocity_FreeStream(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream temperature.
	 * \return Non-dimensionalized freestream temperature.
	 */
	double GetTemperature_FreeStreamND(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream density.
	 * \return Non-dimensionalized freestream density.
	 */
	double GetDensity_FreeStreamND(void);

	/*!
	 * \brief Get the vector of the non-dimensionalized freestream velocity.
	 * \return Non-dimensionalized freestream velocity vector.
	 */
	double* GetVelocity_FreeStreamND(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream energy.
	 * \return Non-dimensionalized freestream energy.
	 */
	double GetEnergy_FreeStreamND(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream viscosity.
	 * \return Non-dimensionalized freestream viscosity.
	 */
	double GetViscosity_FreeStreamND(void);

  /*!
	 * \brief Get the value of the non-dimensionalized freestream viscosity.
	 * \return Non-dimensionalized freestream viscosity.
	 */
	double GetTke_FreeStreamND(void);

  /*!
	 * \brief Get the value of the non-dimensionalized freestream viscosity.
	 * \return Non-dimensionalized freestream viscosity.
	 */
	double GetOmega_FreeStreamND(void);

  /*!
	 * \brief Get the value of the non-dimensionalized freestream viscosity.
	 * \return Non-dimensionalized freestream viscosity.
	 */
	double GetTke_FreeStream(void);

  /*!
	 * \brief Get the value of the non-dimensionalized freestream viscosity.
	 * \return Non-dimensionalized freestream viscosity.
	 */
	double GetOmega_FreeStream(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream intermittency.
	 * \return Non-dimensionalized freestream intermittency.
	 */
	double GetIntermittency_FreeStream(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream turbulence intensity.
	 * \return Non-dimensionalized freestream intensity.
	 */
	double GetTurbulenceIntensity_FreeStream(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream turbulence intensity.
	 * \return Non-dimensionalized freestream intensity.
	 */
	double GetNuFactor_FreeStream(void);

	/*!
	 * \brief Get the value of the turbulent to laminar viscosity ratio.
	 * \return Ratio of turbulent to laminar viscosity ratio.
	 */
	double GetTurb2LamViscRatio_FreeStream(void);

  /*!
	 * \brief Get the vector of free stream mass fraction values.
	 * \return Ratio of species mass to mixture mass.
	 */
	double* GetMassFrac_FreeStream(void);

	/*!
	 * \brief Get the value of the Reynolds length.
	 * \return Reynolds length.
	 */
	double GetLength_Reynolds(void);

	/*!
	 * \brief Get the conversion factor for converting the grid to meters.
	 * \return Conversion factor for converting the grid to meters.
	 */
	double GetMesh_Scale_Change(void);

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
	double GetRefAreaCoeff(void);

	/*!
	 * \brief Get the wave speed.
	 * \return Value of the wave speed.
	 */
	double GetWaveSpeed(void);

	/*!
	 * \brief Get the wave speed.
	 * \return Value of the wave speed.
	 */
	double GetThermalDiffusivity(void);

	/*!
	 * \brief Get the Young's modulus of elasticity.
	 * \return Value of the Young's modulus of elasticity.
	 */
	double GetElasticyMod(void);

	/*!
	 * \brief Get the Poisson's ratio.
	 * \return Value of the Poisson's ratio.
	 */
	double GetPoissonRatio(void);

	/*!
	 * \brief Get the Material Density.
	 * \return Value of the Material Density.
	 */
	double GetMaterialDensity(void);

	/*!
	 * \brief Get the reference length for computing moment (the default value is 1).
	 * \return Reference length for moment computation.
	 */
	double GetRefLengthMoment(void);

	/*!
	 * \brief Get the reference element length for computing the slope limiting epsilon.
	 * \return Reference element length for slope limiting epsilon.
	 */
	double GetRefElemLength(void);

  /*!
	 * \brief Get the reference coefficient for detecting sharp edges.
	 * \return Reference coefficient for detecting sharp edges.
	 */
	double GetRefSharpEdges(void);

	/*!
	 * \brief Get the volume of the whole domain using the fine grid, this value is common for all the grids
	 *        in the multigrid method.
	 * \return Volume of the whole domain.
	 */
	double GetDomainVolume(void);

	/*!
	 * \brief In case the <i>RefAreaCoeff</i> is equal to 0 then, it is necessary to compute a reference area,
	 *        with this function we set the value of the reference area.
	 * \param[in] val_area - Value of the reference area for non dimensional coefficient computation.
	 */
	void SetRefAreaCoeff(double val_area);

	/*!
	 * \brief Set the value of the domain volume computed on the finest grid.
	 * \note This volume do not include the volume of the body that is being simulated.
	 * \param[in] val_volume - Value of the domain volume computed on the finest grid.
	 */
	void SetDomainVolume(double val_volume);

	/*!
	 * \brief Get number of domains in a parallel computation.
	 * \note The number of domains is read from the configuration file or from mpirun -np option.
	 * \return Number of domains in the parallel computation.
	 */
	unsigned short GetnDomain(void);

	/*!
	 * \brief In case we are running the CFD software, the number of domains is read
	 *        from <i>MPI::COMM_WORLD.Get_size()</i>.
	 * \param[in] val_ndomain - Number of domains for the MPI parallelization.
	 */
	void SetnDomain(unsigned short val_ndomain);

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
	 */
	void SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme, unsigned short val_kind_centered,
			unsigned short val_kind_upwind, unsigned short val_kind_slopelimit, unsigned short val_order_spatial_int);

	/*!
	 * \brief Get the value of limiter coefficient.
	 * \return Value of the limiter coefficient.
	 */
	double GetLimiterCoeff(void);
  
  /*!
	 * \brief Freeze the value of the limiter after a number of iterations.
	 * \return Number of iterations.
	 */
	unsigned long GetLimiterIter(void);

  /*!
	 * \brief Get the value of sharp edge limiter.
	 * \return Value of the sharp edge limiter coefficient.
	 */
	double GetSharpEdgesCoeff(void);

	/*!
	 * \brief Get the Reynolds number. Dimensionless number that gives a measure of the ratio of inertial forces
	 *        to viscous forces and consequently quantifies the relative importance of these two types of forces
	 *        for given flow condition.
	 * \return Value of the Reynolds number.
	 */
	double GetReynolds(void);

	/*!
	 * \brief Get the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	double GetFroude(void);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetFroude(double val_froude);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetMach(double val_mach);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetReynolds(double val_reynolds);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetLength_Ref(double val_length_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetVelocity_Ref(double val_velocity_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetPressure_Ref(double val_pressure_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetDensity_Ref(double val_density_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetTime_Ref(double val_time_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetEnergy_Ref(double val_energy_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetOmega_Ref(double val_omega_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetForce_Ref(double val_force_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetGas_Constant_Ref(double val_gas_constant_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetGas_Constant(double val_gas_constant);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetViscosity_Ref(double val_viscosity_ref);

   /*!
    * \brief Set the Froude number for free surface problems.
	* \return Value of the Froude number.
	*/
	void SetConductivity_Ref(double val_conductivity_ref);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetPressure_FreeStreamND(double val_pressure_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetPressure_FreeStream(double val_pressure_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetDensity_FreeStreamND(double val_density_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetDensity_FreeStream(double val_density_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetViscosity_FreeStream(double val_viscosity_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetModVel_FreeStream(double val_modvel_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetModVel_FreeStreamND(double val_modvel_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetTemperature_FreeStream(double val_temperature_freestream);
  
  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetTemperature_FreeStreamND(double val_temperature_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetGas_ConstantND(double val_gas_constantnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetVelocity_FreeStreamND(double val_velocity_freestreamnd, unsigned short val_dim);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetViscosity_FreeStreamND(double val_viscosity_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetTke_FreeStreamND(double val_tke_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetOmega_FreeStreamND(double val_omega_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetTke_FreeStream(double val_tke_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetOmega_FreeStream(double val_omega_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetEnergy_FreeStreamND(double val_energy_freestreamnd);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetEnergy_FreeStream(double val_energy_freestream);

  /*!
	 * \brief Set the Froude number for free surface problems.
	 * \return Value of the Froude number.
	 */
	void SetTotal_UnstTimeND(double val_total_unsttimend);

	/*!
	 * \brief Get the angle of attack of the body. This is the angle between a reference line on a lifting body
	 *        (often the chord line of an airfoil) and the vector representing the relative motion between the
	 *        lifting body and the fluid through which it is moving.
	 * \return Value of the angle of attack.
	 */
	double GetAoA(void);

	/*!
	 * \brief Set the angle of attack.
	 * \param[in] val_AoA - Value of the angle of attack.
	 */
	void SetAoA(double val_AoA);

  /*!
	 * \brief Set the angle of attack.
	 * \param[in] val_AoA - Value of the angle of attack.
	 */
	void SetAoS(double val_AoS);

	/*!
	 * \brief Get the angle of sideslip of the body. It relates to the rotation of the aircraft centerline from
	 *        the relative wind.
	 * \return Value of the angle of sideslip.
	 */
	double GetAoS(void);

	/*!
	 * \brief Get the charge coefficient that is used in the poissonal potential simulation.
	 * \return Value of the charge coefficient.
	 */
	double GetChargeCoeff(void);

	/*!
	 * \brief Get the number of multigrid levels.
	 * \return Number of multigrid levels (without including the original grid).
	 */
	unsigned short GetMGLevels(void);

	/*!
	 * \brief Set the number of multigrid levels.
	 * \param[in] val_nMultiLevel - Index of the mesh were the CFL is applied
	 */
	void SetMGLevels(unsigned short val_nMultiLevel);

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
	double GetCFL(unsigned short val_mesh);
  
  /*!
	 * \brief Get the Courant Friedrich Levi number for each grid.
	 * \param[in] val_mesh - Index of the mesh were the CFL is applied.
	 * \return CFL number for each grid.
	 */
	void SetCFL(unsigned short val_mesh, double val_cfl);

	/*!
	 * \brief Get the Courant Friedrich Levi number for each grid, for each species
	 * \param[in] val_mesh - Index of the mesh were the CFL is applied.
	 * \param[in] val_Species - Index of the chemical species
	 * \return CFL number for each grid.
	 */
	double GetCFL(unsigned short val_mesh, unsigned short val_Species);

	/*!
	 * \brief Get the Courant Friedrich Levi number for unsteady simulations.
	 * \return CFL number for unsteady simulations.
	 */
	double GetUnst_CFL(void);
  
  /*!
   * \brief Get the Courant Friedrich Levi number for unsteady simulations.
   * \return CFL number for unsteady simulations.
   */
  double GetMax_DeltaTime(void);
  
	/*!
	 * \brief Get a parameter of the particular design variable.
	 * \param[in] val_dv - Number of the design variable that we want to read.
	 * \param[in] val_param - Index of the parameter that we want to read.
	 * \return Design variable parameter.
	 */
	double GetParamDV(unsigned short val_dv, unsigned short val_param);

  /*!
	 * \brief Get the FFD Tag of a particular design variable.
	 * \param[in] val_dv - Number of the design variable that we want to read.
	 * \return Design variable parameter.
	 */
	string GetFFDTag(unsigned short val_dv);

	/*!
	 * \brief Get the number of design variables.
	 * \return Number of the design variables.
	 */
	unsigned short GetnDV(void);

	/*!
	 * \brief Get the number of Runge-Kutta steps.
	 * \return Number of Runge-Kutta steps.
	 */
	unsigned short GetnRKStep(void);

	/*!
	 * \brief Get the total number of boundary markers.
	 * \return Total number of boundary markers.
	 */
	unsigned short GetnMarker_All(void);

    /*!
	 * \brief Get the total number of boundary markers.
	 * \return Total number of boundary markers.
	 */
	unsigned short GetnMarker_EngineInflow(void);
  
  /*!
   * \brief Get the total number of boundary markers.
   * \return Total number of boundary markers.
   */
  unsigned short GetnMarker_EngineBleed(void);

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
	unsigned short GetnMarker_ActDisk_Inlet(void);

  /*!
	 * \brief Get the total number of boundary markers.
	 * \return Total number of boundary markers.
	 */
	unsigned short GetnMarker_ActDisk_Outlet(void);

  /*!
   * \brief Get the total number of 1D output markers.
   * \return Total number of monitoring markers.
   */
  unsigned short GetnMarker_Out_1D(void);


    /*!
	 * \brief Get the total number of monitoring markers.
	 * \return Total number of monitoring markers.
	 */
	unsigned short GetnMarker_Monitoring(void);

  /*!
	 * \brief Get the total number of moving markers.
	 * \return Total number of moving markers.
	 */
	unsigned short GetnMarker_Moving(void);

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
	 * \brief Retrieves the number of periodic time instances for Time Spectral.
	 * \return: Number of periodic time instances for Time Spectral.
	 */
	unsigned short GetnTimeInstances(void);

	/*!
	 * \brief Retrieves the period of oscillations to be used with Time Spectral.
	 * \return: Period for Time Spectral.
	 */
	double GetTimeSpectral_Period(void);

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
	 * \brief Set the current internal iteration number.
	 * \param[in] val_iter - Current external iteration number.
	 */
	void SetIntIter(unsigned long val_iter);

	/*!
	 * \brief Get the current internal iteration number.
	 * \return Current external iteration.
	 */
	unsigned long GetExtIter(void);

	/*!
	 * \brief Get the current external iteration number.
	 * \return Current external iteration.
	 */
	unsigned long GetIntIter(void);

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
	 * \brief Get information about performing a low fidelity simulation.
	 * \return 	<code>TRUE</code> means that a low fidelity simulation will be performed.
	 */
	bool GetLowFidelitySim(void);

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
	 * \brief Get information about writing a restart solution file.
	 * \return <code>TRUE</code> means that a restart solution file will be written.
	 */
	bool GetWrt_Restart(void);

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
	 * \brief Get information about writing sectional force files.
	 * \return <code>TRUE</code> means that sectional force files will be written for specified markers.
	 */
	bool GetPlot_Section_Forces(void);

  /*!
   * \brief Get information about writing average stagnation pressure
   * \return <code>TRUE</code> means that the average stagnation pressure will be output for specified markers.
   */
  bool GetWrt_1D_Output(void);

	/*!
	 * \brief Get the alpha (convective) coefficients for the Runge-Kutta integration scheme.
	 * \param[in] val_step - Index of the step.
	 * \return Alpha coefficient for the Runge-Kutta integration scheme.
	 */
	double Get_Alpha_RKStep(unsigned short val_step);

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
	string GetMarker_EngineInflow(unsigned short val_marker);
  
  /*!
   * \brief Get the index of the surface defined in the geometry file.
   * \param[in] val_marker - Value of the marker in which we are interested.
   * \return Value of the index that is in the geometry file for the surface that
   *         has the marker <i>val_marker</i>.
   */
  string GetMarker_EngineBleed(unsigned short val_marker);

	/*!
	 * \brief Get the index of the surface defined in the geometry file.
	 * \param[in] val_marker - Value of the marker in which we are interested.
	 * \return Value of the index that is in the geometry file for the surface that
	 *         has the marker <i>val_marker</i>.
	 */
	string GetMarker_EngineExhaust(unsigned short val_marker);

    /*!
	 * \brief Get the name of the surface defined in the geometry file.
	 * \param[in] val_marker - Value of the marker in which we are interested.
	 * \return Name that is in the geometry file for the surface that
	 *         has the marker <i>val_marker</i>.
	 */
	string GetMarker_Monitoring(unsigned short val_marker);

	/*!
	 * \brief Get the tag if the iMarker defined in the geometry file.
	 * \param[in] val_tag - Value of the tag in which we are interested.
	 * \return Value of the marker <i>val_marker</i> that is in the geometry file
	 *         for the surface that has the tag.
	 */
  short GetTagBound_Marker_All(string val_tag);

	/*!
	 * \brief Get the kind of boundary for each marker.
	 * \param[in] val_marker - Index of the marker in which we are interested.
	 * \return Kind of boundary for the marker <i>val_marker</i>.
	 */
	unsigned short GetMarker_All_KindBC(unsigned short val_marker);

  /*!
   * \brief Get the kind of boundary for each marker.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \return Kind of boundary for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_All_Out_1D(unsigned short val_marker);

  /*!
   * \brief Set the value of the boundary <i>val_boundary</i> (read from the config file)
   *        for the marker <i>val_marker</i>.
   * \param[in] val_marker - Index of the marker in which we are interested.
   * \param[in] val_boundary - Kind of boundary read from config file.
   */
  void SetMarker_All_Out_1D(unsigned short val_marker, unsigned short val_boundary);


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
	double GetPressure_Critical(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	double GetTemperature_Critical(void);

	/*!
	 * \brief Get the value of the critical pressure.
	 * \return Critical pressure.
	 */
	double GetAcentric_Factor(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	unsigned short GetKind_ViscosityModel(void);

	/*!
	 * \brief Get the value of the thermal conductivity .
	 * \return Critical temperature.
	 */
	unsigned short GetKind_ConductivityModel(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	double GetMu_ConstantND(void);

	/*!
	 * \brief Get the value of the non-dimensional thermal conductivity.
	 * \return Critical temperature.
	 */
	double GetKt_ConstantND(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	double GetMu_RefND(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	double GetMu_Temperature_RefND(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	double GetMu_SND(void);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	void SetMu_ConstantND(double mu_const);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	void SetKt_ConstantND(double kt_const);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	void SetMu_RefND(double mu_ref);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	void SetMu_Temperature_RefND(double mu_Tref);

	/*!
	 * \brief Get the value of the critical temperature.
	 * \return Critical temperature.
	 */
	void SetMu_SND(double mu_s);

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
	 * \brief Set the kind of preconditioner for the implicit solver.
	 * \return Numerical preconditioner for implicit formulation (solving the linear system).
	 */
	void SetKind_Linear_Solver_Prec(unsigned short val_kind_prec);

	/*!
	 * \brief Get min error of the linear solver for the implicit formulation.
	 * \return Min error of the linear solver for the implicit formulation.
	 */
	double GetLinear_Solver_Error(void);

	/*!
	 * \brief Get max number of iterations of the linear solver for the implicit formulation.
	 * \return Max number of iterations of the linear solver for the implicit formulation.
	 */
	unsigned long GetLinear_Solver_Iter(void);

  /*!
   * \brief Get restart frequency of the linear solver for the implicit formulation.
   * \return Restart frequency of the linear solver for the implicit formulation.
   */
  unsigned long GetLinear_Solver_Restart_Frequency(void);

	/*!
	 * \brief Get the relaxation coefficient of the linear solver for the implicit formulation.
	 * \return relaxation coefficient of the linear solver for the implicit formulation.
	 */
	double GetLinear_Solver_Relax(void);

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
	 * \brief Set the kind of preconditioner for the implicit solver.
	 * \return Numerical preconditioner for implicit formulation (solving the linear system).
	 */
	void SetKind_AdjTurb_Linear_Prec(unsigned short val_kind_prec);

	/*!
	 * \brief Get min error of the linear solver for the implicit formulation.
	 * \return Min error of the linear solver for the implicit formulation.
	 */
	double GetAdjTurb_Linear_Error(void);
  
  /*!
	 * \brief Get the entropy fix.
	 * \return Vaule of the entropy fix.
	 */
	double GetEntropyFix_Coeff(void);

	/*!
	 * \brief Get max number of iterations of the linear solver for the implicit formulation.
	 * \return Max number of iterations of the linear solver for the implicit formulation.
	 */
	unsigned short GetAdjTurb_Linear_Iter(void);

	/*!
	 * \brief Get CFL reduction factor for adjoint turbulence model.
	 * \return CFL reduction factor.
	 */
	double GetCFLRedCoeff_AdjTurb(void);

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
	double GetDeform_Tol_Factor(void);

  /*!
   * \brief Get Young's modulus for deformation (constant stiffness deformation)
   */
  double GetDeform_ElasticityMod(void);

  /*!
   * \brief Get Poisson's ratio for deformation (constant stiffness deformation)
   * \
   */
  double GetDeform_PoissonRatio(void);

  /*!
	 * \brief Get the type of stiffness to impose for FEA mesh deformation.
	 * \return type of stiffness to impose for FEA mesh deformation.
	 */
	unsigned short GetDeform_Stiffness_Type(void);

	/*!
	 * \brief Creates a teot file to visualize the deformation made by the MDC software.
	 * \return <code>TRUE</code> if the deformation is going to be plotted; otherwise <code>FALSE</code>.
	 */
	bool GetVisualize_Deformation(void);

	/*!
	 * \brief Get the kind of SU2 software component.
	 * \return Kind of the SU2 software component.
	 */
	unsigned short GetKind_SU2(void);

	/*!
	 * \brief Get the kind of the turbulence model.
	 * \return Kind of the turbulence model.
	 */
	unsigned short GetKind_Turb_Model(void);

  /*!
	 * \brief Get the file containing the ML model
	 */
	string GetML_Turb_Model_File(void);

  /*!
	 * \brief File containing a check for the proper creation of the turb model
	 * \return Temporary ml->SU2 file name.
	 */
  string GetML_Turb_Model_FeatureSet(void);

  /*!
	 * \brief File containing a check for the proper creation of the turb model
	 * \return Temporary ml->SU2 file name.
	 */
  string* GetML_Turb_Model_Extra(void);

  /*!
	 * \brief File containing a check for the proper creation of the turb model
	 * \return Temporary ml->SU2 file name.
	 */
  unsigned short GetNumML_Turb_Model_Extra(void);

	/*!
	 * \brief Get the kind of the transition model.
	 * \return Kind of the transion model.
	 */
	unsigned short GetKind_Trans_Model(void);

	/*!
	 * \brief Get the kind of adaptation technique.
	 * \return Kind of adaptation technique.
	 */
	unsigned short GetKind_Adaptation(void);

	/*!
	 * \brief Get the number of new elements added in the adaptation process.
	 * \return percentage of new elements that are going to be added in the adaptation.
	 */
	double GetNew_Elem_Adapt(void);

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
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder(void);

  /*!
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder_Flow(void);

  /*!
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder_Turb(void);

  /*!
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder_TNE2(void);

  /*!
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder_AdjLevelSet(void);

  /*!
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder_AdjFlow(void);

  /*!
	 * \brief Get the order of the spatial integration.
	 * \note This is the information that the code will use, the method will
	 *       change in runtime depending of the specific equation (direct, adjoint,
	 *       linearized) that is being solved.
	 * \return Kind of upwind scheme for the convective terms.
	 */
	unsigned short GetSpatialOrder_AdjTNE2(void);

	/*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_Flow(void);

  /*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_TNE2(void);

	/*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the plasma equations.
	 */
	unsigned short GetKind_TimeIntScheme_Wave(void);

  /*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the plasma equations.
	 */
	unsigned short GetKind_TimeIntScheme_Heat(void);

  /*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the plasma equations.
	 */
	unsigned short GetKind_TimeIntScheme_Poisson(void);

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
	 * \brief Get the kind of convective numerical scheme for the flow
	 *        equations (centered or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_Flow(void);

  /*!
	 * \brief Get the kind of convective numerical scheme for the flow
	 *        equations (centered or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_TNE2(void);

  /*!
	 * \brief Get the kind of convective numerical scheme for the flow
	 *        equations (centered or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_AdjTNE2(void);

	/*!
	 * \brief Get the kind of convective numerical scheme for the template
	 *        equations (centered or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_Template(void);

	/*!
	 * \brief Get the kind of convective numerical scheme for the adjoint level set
	 *        equations (centered or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the level set equation.
	 */
	unsigned short GetKind_ConvNumScheme_AdjLevelSet(void);

	/*!
	 * \brief Get the kind of center convective numerical scheme for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Centered_Flow(void);

  /*!
	 * \brief Get the kind of center convective numerical scheme for the two-temperature model.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Centered_TNE2(void);

  /*!
	 * \brief Get the kind of center convective numerical scheme for the two-temperature model.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Centered_AdjTNE2(void);

	/*!
	 * \brief Get the kind of center convective numerical scheme for the adjoint level set equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the level set equations.
	 */
	unsigned short GetKind_Centered_AdjLevelSet(void);

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
	 * \brief Get the kind of upwind convective numerical scheme for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_TNE2(void);

  /*!
	 * \brief Get the kind of upwind convective numerical scheme for the flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_AdjTNE2(void);

	/*!
	 * \brief Get the kind of upwind convective numerical scheme for the adjoint level set equation.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_AdjLevelSet(void);

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
	 * \return Method for limiting the spatial gradients solving the flow equations.
	 */
	unsigned short GetKind_SlopeLimit_TNE2(void);

  /*!
	 * \brief Get the method for limiting the spatial gradients.
	 * \return Method for limiting the spatial gradients solving the flow equations.
	 */
	unsigned short GetKind_SlopeLimit_AdjTNE2(void);

	/*!
	 * \brief Get the method for limiting the spatial gradients.
	 * \return Method for limiting the spatial gradients solving the turbulent equation.
	 */
	unsigned short GetKind_SlopeLimit_Turb(void);

	/*!
	 * \brief Get the method for limiting the spatial gradients.
	 * \return Method for limiting the spatial gradients solving the level set equation.
	 */
	unsigned short GetKind_SlopeLimit_AdjLevelSet(void);

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
	double GetKappa_1st_Flow(void);

	/*!
	 * \brief Value of the calibrated constant for the JST method (center scheme).
	 * \return Calibrated constant for the JST method for the flow equations.
	 */
	double GetKappa_2nd_Flow(void);

	/*!
	 * \brief Value of the calibrated constant for the JST method (center scheme).
	 * \return Calibrated constant for the JST method for the flow equations.
	 */
	double GetKappa_4th_Flow(void);

  /*!
	 * \brief Value of the calibrated constant for the Lax method (center scheme).
	 * \note This constant is used in coarse levels and with first order methods.
	 * \return Calibrated constant for the Lax method.
	 */
	double GetKappa_1st_TNE2(void);

	/*!
	 * \brief Value of the calibrated constant for the JST method (center scheme).
	 * \return Calibrated constant for the JST method for the flow equations.
	 */
	double GetKappa_2nd_TNE2(void);

	/*!
	 * \brief Value of the calibrated constant for the JST method (center scheme).
	 * \return Calibrated constant for the JST method for the flow equations.
	 */
	double GetKappa_4th_TNE2(void);

	/*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the adjoint flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_AdjFlow(void);

  /*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the adjoint flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_AdjTNE2(void);

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
	double GetKappa_2nd_AdjFlow(void);

	/*!
	 * \brief Value of the calibrated constant for the high order method (center scheme).
	 * \return Calibrated constant for the high order center method for the adjoint flow equations.
	 */
	double GetKappa_4th_AdjFlow(void);

	/*!
	 * \brief Value of the calibrated constant for the low order method (center scheme).
	 * \return Calibrated constant for the low order center method for the adjoint flow equations.
	 */
	double GetKappa_1st_AdjFlow(void);

  /*!
	 * \brief Value of the calibrated constant for the high order method (center scheme).
	 * \return Calibrated constant for the high order center method for the adjoint flow equations.
	 */
	double GetKappa_2nd_AdjTNE2(void);

	/*!
	 * \brief Value of the calibrated constant for the high order method (center scheme).
	 * \return Calibrated constant for the high order center method for the adjoint flow equations.
	 */
	double GetKappa_4th_AdjTNE2(void);

	/*!
	 * \brief Value of the calibrated constant for the low order method (center scheme).
	 * \return Calibrated constant for the low order center method for the adjoint flow equations.
	 */
	double GetKappa_1st_AdjTNE2(void);

	/*!
	 * \brief Get the kind of integration scheme (explicit or implicit)
	 *        for the linearized flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the linearized flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_LinFlow(void);

	/*!
	 * \brief Get the kind of convective numerical scheme for the linearized flow
	 *        equations (centered or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the linearized flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_LinFlow(void);

	/*!
	 * \brief Get the kind of center convective numerical scheme for the linearized flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the linearized flow equations.
	 */
	unsigned short GetKind_Centered_LinFlow(void);

	/*!
	 * \brief Get the kind of upwind convective numerical scheme for the linearized flow equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the linearized flow equations.
	 */
	unsigned short GetKind_Upwind_LinFlow(void);

	/*!
	 * \brief Value of the calibrated constant for the high order method (center scheme).
	 * \return Calibrated constant for the high order center method for the linearized flow equations.
	 */
	double GetKappa_4th_LinFlow(void);

	/*!
	 * \brief Value of the calibrated constant for the low order method (center scheme).
	 * \return Calibrated constant for the low order center method for the linearized flow equations.
	 */
	double GetKappa_1st_LinFlow(void);

	/*!
	 * \brief Get the kind of integration scheme (implicit)
	 *        for the turbulence equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the turbulence equations.
	 */
	unsigned short GetKind_TimeIntScheme_Turb(void);

	/*!
	 * \brief Get the kind of integration scheme (implicit)
	 *        for the level set equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the level set equations.
	 */
	unsigned short GetKind_TimeIntScheme_AdjLevelSet(void);

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
	 *        adjoint method.
	 * \return <code>FALSE</code> means that the adjoint turbulence equations will be used.
	 */
	bool GetFrozen_Visc(void);

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
	 * \brief Get the number of sections.
	 * \return Number of sections
	 */
	unsigned short GetnSections(void);

  /*!
	 * \brief Get the number of sections for computing internal volume.
	 * \return Number of sections for computing internal volume.
	 */
	unsigned short GetnVolSections(void);

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
	 * \brief Get the kind of sensitivity smoothing technique.
	 * \return Kind of sensitivity smoothing technique.
	 */
	unsigned short GetKind_SensSmooth(void);

	/*!
	 * \brief Get equations to be treated continuously. There are several options: Euler, Navier Stokes
	 * \return Continuous equations.
	 */
	unsigned short GetContinuous_Eqns(void);

	/*!
	 * \brief Get equations to be treated discretely. There are several options: None, SA, SST
	 * \return Discrete equations.
	 */
	unsigned short GetDiscrete_Eqns(void);

	/*!
	 * \brief Provides information about the time integration, and change the write in the output
	 *        files information about the iteration.
	 * \return The kind of time integration: Steady state, time stepping method (unsteady) or
	 *         dual time stepping method (unsteady).
	 */
	unsigned short GetUnsteady_Simulation(void);

	/*!
	 * \brief Provides the number of species present in the plasma
	 * \return: The number of species present in the plasma, read from input file
	 */
	unsigned short GetnSpecies(void);

	/*!
	 * \brief Provides the number of chemical reactions in the chemistry model
	 * \return: The number of chemical reactions, read from input file
	 */
	unsigned short GetnReactions(void);

	/*!
	 * \brief Provides the number of chemical reactions in the chemistry model
	 * \return: The number of chemical reactions, read from input file
	 */
	double GetArrheniusCoeff(unsigned short iReaction);

	/*!
	 * \brief Provides the number of chemical reactions in the chemistry model
	 * \return: The number of chemical reactions, read from input file
	 */
	double GetArrheniusEta(unsigned short iReaction);

	/*!
	 * \brief Provides the number of chemical reactions in the chemistry model
	 * \return: The number of chemical reactions, read from input file
	 */
	double GetArrheniusTheta(unsigned short iReaction);

  /*!
	 * \brief Provides the rate controlling temperature exponents for chemistry.
	 * \return: Rate controlling temperature exponents.
	 */
  double* GetRxnTcf_a(void);

  /*!
	 * \brief Provides the rate controlling temperature exponents for chemistry.
	 * \return: Rate controlling temperature exponents.
	 */
  double* GetRxnTcf_b(void);

  /*!
	 * \brief Provides the rate controlling temperature exponents for chemistry.
	 * \return: Rate controlling temperature exponents.
	 */
  double* GetRxnTcb_a(void);

  /*!
	 * \brief Provides the rate controlling temperature exponents for chemistry.
	 * \return: Rate controlling temperature exponents.
	 */
  double* GetRxnTcb_b(void);

  /*!
	 * \brief Dissociation potential of species.
	 * \return: Dissociation potential.
	 */
	double* GetDissociationPot(void);

	/*!
	 * \brief Provides the number of rotational modes of energy storage
	 * \return: Vector of rotational mode count
	 */
  double* GetRotationModes(void);

	/*!
	 * \brief Provides the characteristic vibrational temperature for calculating e_vib
	 * \return: Vector of characteristic vibrational temperatures [K]
	 */
	double* GetCharVibTemp(void);

  /*!
	 * \brief Provides the characteristic electronic temperature for calculating e_el
	 * \return: Vector of characteristic vibrational temperatures [K]
	 */
	double** GetCharElTemp(void);

  /*!
	 * \brief Provides the degeneracy of electron states for calculating e_el
	 * \return: Vector of characteristic vibrational temperatures [K]
	 */
	double** GetElDegeneracy(void);

  /*!
	 * \brief Provides number electron states for calculating e_el
	 * \return: Vector of number of electron states for each species
	 */
	unsigned short* GetnElStates(void);


  /*!
	 * \brief Provides the thermodynamic reference temperatures from the JANAF tables
	 * \return: Vector of reference temperatures [K]
	 */
  double* GetRefTemperature(void);

  /*!
	 * \brief Provides the characteristic vibrational temperature for calculating e_vib
	 * \return: The number of chemical reactions, read from input file
	 */
	double GetCharVibTemp(unsigned short iSpecies);

	/*!
	 * \brief Provides a table of equilibrium constants for a particular chemical reaction for a supplied gas model.
	 * \return: Matrix of reaction constants
	 */
	void GetChemistryEquilConstants(double **RxnConstantTable, unsigned short iReaction);

	/*!
	 * \brief Provides the molar mass of each species present in multi species fluid
	 * \return: Vector of molar mass of each species in kg/kmol
	 */
	double* GetMolar_Mass(void);

  /*!
	 * \brief Provides the molar mass of each species present in multi species fluid
	 * \return: Mass of each species in Kg
	 */
	double GetMolar_Mass(unsigned short iSpecies);

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
	double GetInitial_Gas_Composition(unsigned short iSpecies);

	/*!
	 * \brief Retrieves the multi-species fluid mixture molar mass.
	 * \return: Molar mass of the fluid mixture
	 */
	double GetMixtureMolar_Mass();

  /*!
	 * \brief Provides the formation enthalpy of the specified species at standard conditions
	 * \return: Enthalpy of formation
	 */
	double* GetEnthalpy_Formation(void);

	/*!
	 * \brief Provides the formation enthalpy of the specified species at standard conditions
	 * \return: Enthalpy of formation
	 */
	double GetEnthalpy_Formation(unsigned short iSpecies);

	/*!
	 * \brief Provides the restart information.
	 * \return Restart information, if <code>TRUE</code> then the code will use the solution as restart.
	 */
	bool GetRestart(void);

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
	 * \brief Information about doing a full multigrid strategy (start in the coarse level).
	 * \return <code>TRUE</code> or <code>FALSE</code>  depending if we are performing a full multigrid strategy.
	 */
	bool GetFullMG(void);

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
	 * \brief Get the name of the file with the solution of the linearized flow problem.
	 * \return Name of the file with the solution of the linearized flow problem.
	 */
	string GetSolution_LinFileName(void);

	/*!
	 * \brief Get the name of the file with the solution of the adjoint flow problem
	 *		  with drag objective function.
	 * \return Name of the file with the solution of the adjoint flow problem with
	 *         drag objective function.
	 */
	string GetSolution_AdjFileName(void);

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
	 * \brief Get the name of the file with the convergence history of the problem.
	 * \return Name of the file with convergence history of the problem.
	 */
	string GetConv_FileName(void);

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
	 * \brief Get the name of the file with the structure variables.
	 * \return Name of the file with the structure variables.
	 */
	string GetSurfWave_FileName(void);

  /*!
	 * \brief Get the name of the file with the structure variables.
	 * \return Name of the file with the structure variables.
	 */
	string GetSurfHeat_FileName(void);

	/*!
	 * \brief Get the name of the file with the wave variables.
	 * \return Name of the file with the wave variables.
	 */
	string GetWave_FileName(void);

  /*!
	 * \brief Get the name of the file with the wave variables.
	 * \return Name of the file with the wave variables.
	 */
	string GetHeat_FileName(void);

	/*!
	 * \brief Get the name of the file with the adjoint wave variables.
	 * \return Name of the file with the adjoint wave variables.
	 */
	string GetAdjWave_FileName(void);

	/*!
	 * \brief Get the name of the restart file for the wave variables.
	 * \return Name of the restart file for the flow variables.
	 */
	string GetRestart_WaveFileName(void);

	/*!
	 * \brief Get the name of the restart file for the heat variables.
	 * \return Name of the restart file for the flow variables.
	 */
	string GetRestart_HeatFileName(void);

	/*!
	 * \brief Get the name of the restart file for the flow variables.
	 * \return Name of the restart file for the flow variables.
	 */
	string GetRestart_FlowFileName(void);

	/*!
	 * \brief Get the name of the restart file for the linearized flow variables.
	 * \return Name of the restart file for the linearized flow variables.
	 */
	string GetRestart_LinFileName(void);

	/*!
	 * \brief Get the name of the restart file for the adjoint variables (drag objective function).
	 * \return Name of the restart file for the adjoint variables (drag objective function).
	 */
	string GetRestart_AdjFileName(void);

	/*!
	 * \brief Get the name of the file with the adjoint variables.
	 * \return Name of the file with the adjoint variables.
	 */
	string GetAdj_FileName(void);

	/*!
	 * \brief Get the name of the file with the linearized flow variables.
	 * \return Name of the file with the linearized flow variables.
	 */
	string GetLin_FileName(void);

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
	 * \brief Get the name of the file with the surface information for the linearized flow problem.
	 * \return Name of the file with the surface information for the linearized flow problem.
	 */
	string GetSurfLinCoeff_FileName(void);

  /*!
	 * \brief Augment the input filename with the iteration number for an unsteady file.
   * \param[in] val_filename - String value of the base filename.
   * \param[in] val_iter - Unsteady iteration number or time spectral instance.
	 * \return Name of the file with the iteration numer for an unsteady solution file.
	 */
  string GetUnsteady_FileName(string val_filename, int val_iter);

  /*!
	 * \brief Append the input filename string with the appropriate objective function extension.
   * \param[in] val_filename - String value of the base filename.
	 * \return Name of the file with the appropriate objective function extension.
	 */
  string GetObjFunc_Extension(string val_filename);
  
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
	 * \brief Get functional that is going to be used to evaluate the linearized flow convergence.
	 * \return Functional that is going to be used to evaluate the linearized flow convergence.
	 */
	unsigned short GetCauchy_Func_LinFlow(void);

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
	double GetCauchy_Eps(void);

	/*!
	 * \brief Get the value of convergence criteria for the one-shot problem.
	 * \return Value of the convergence criteria.
	 */
	double GetCauchy_Eps_OneShot(void);

	/*!
	 * \brief Get the value of convergence criteria for the full multigrid method.
	 * \return Value of the convergence criteria.
	 */
	double GetCauchy_Eps_FullMG(void);

	/*!
	 * \brief If we are prforming an unsteady simulation, there is only
	 *        one value of the time step for the complete simulation.
	 * \return Value of the time step in an unsteady simulation (non dimensional).
	 */
	double GetDelta_UnstTimeND(void);

  /*!
	 * \brief If we are prforming an unsteady simulation, there is only
	 *        one value of the time step for the complete simulation.
	 * \return Value of the time step in an unsteady simulation (non dimensional).
	 */
	double GetTotal_UnstTimeND(void);

	/*!
	 * \brief If we are prforming an unsteady simulation, there is only
	 *        one value of the time step for the complete simulation.
	 * \return Value of the time step in an unsteady simulation.
	 */
	double GetDelta_UnstTime(void);

	/*!
	 * \brief Set the value of the unsteadty time step using the CFL number.
	 * \param[in] val_delta_unsttimend - Value of the unsteady time step using CFL number.
	 */
	void SetDelta_UnstTimeND(double val_delta_unsttimend);

	/*!
	 * \brief If we are performing an unsteady simulation, this is the
	 * 	value of max physical time for which we run the simulation
	 * \return Value of the physical time in an unsteady simulation.
	 */
	double GetTotal_UnstTime(void);

	/*!
	 * \brief If we are performing an unsteady simulation, this is the
	 * 	value of current time.
	 * \return Value of the physical time in an unsteady simulation.
	 */
	double GetCurrent_UnstTime(void);

	/*!
	 * \brief Divide the rectbles and hexahedron.
	 * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
	 */
	bool GetDivide_Element(void);

  /*!
	 * \brief Divide the rectbles and hexahedron.
	 * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
	 */
	bool GetEngine_Intake(void);

	/*!
	 * \brief Value of the design variable step, we use this value in design problems.
	 * \param[in] val_dv - Number of the design variable that we want to read.
	 * \return Design variable step.
	 */
	double GetDV_Value(unsigned short val_dv);

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
	 * \brief Set the type of dynamic mesh motion.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \param[in] motion_Type - Specify motion type.
	 */
	void SetKind_GridMovement(unsigned short val_iZone, unsigned short motion_Type);

	/*!
	 * \brief Get the mach number based on the mesh velocity and freestream quantities.
	 * \return Mach number based on the mesh velocity and freestream quantities.
	 */
	double GetMach_Motion(void);

	/*!
	 * \brief Get x-coordinate of the mesh motion origin.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return X-coordinate of the mesh motion origin.
	 */
	double GetMotion_Origin_X(unsigned short val_iZone);

	/*!
	 * \brief Get y-coordinate of the mesh motion origin
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Y-coordinate of the mesh motion origin.
	 */
	double GetMotion_Origin_Y(unsigned short val_iZone);

	/*!
	 * \brief Get z-coordinate of the mesh motion origin
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Z-coordinate of the mesh motion origin.
	 */
	double GetMotion_Origin_Z(unsigned short val_iZone);

	/*!
	 * \brief Set x-coordinate of the mesh motion origin.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \param[in] val_origin - New x-coordinate of the mesh motion origin.
	 */
	void SetMotion_Origin_X(unsigned short val_iZone, double val_origin);

	/*!
	 * \brief Set y-coordinate of the mesh motion origin
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \param[in] val_origin - New y-coordinate of the mesh motion origin.
	 */
	void SetMotion_Origin_Y(unsigned short val_iZone, double val_origin);

	/*!
	 * \brief Set z-coordinate of the mesh motion origin
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \param[in] val_origin - New y-coordinate of the mesh motion origin.
	 */
	void SetMotion_Origin_Z(unsigned short val_iZone, double val_origin);

	/*!
	 * \brief Get the translational velocity of the mesh in the x-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Translational velocity of the mesh in the x-direction.
	 */
	double GetTranslation_Rate_X(unsigned short val_iZone);

	/*!
	 * \brief Get the translational velocity of the mesh in the y-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Translational velocity of the mesh in the y-direction.
	 */
	double GetTranslation_Rate_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the translational velocity of the mesh in the z-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Translational velocity of the mesh in the z-direction.
	 */
	double GetTranslation_Rate_Z(unsigned short val_iZone);

	/*!
	 * \brief Get the angular velocity of the mesh about the x-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular velocity of the mesh about the x-axis.
	 */
	double GetRotation_Rate_X(unsigned short val_iZone);

	/*!
	 * \brief Get the angular velocity of the mesh about the y-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular velocity of the mesh about the y-axis.
	 */
	double GetRotation_Rate_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the angular velocity of the mesh about the z-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular velocity of the mesh about the z-axis.
	 */
	double GetRotation_Rate_Z(unsigned short val_iZone);

	/*!
	 * \brief Get the angular frequency of a mesh pitching about the x-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular frequency of a mesh pitching about the x-axis.
	 */
	double GetPitching_Omega_X(unsigned short val_iZone);

	/*!
	 * \brief Get the angular frequency of a mesh pitching about the y-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular frequency of a mesh pitching about the y-axis.
	 */
	double GetPitching_Omega_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the angular frequency of a mesh pitching about the z-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular frequency of a mesh pitching about the z-axis.
	 */
	double GetPitching_Omega_Z(unsigned short val_iZone);

	/*!
	 * \brief Get the pitching amplitude about the x-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Pitching amplitude about the x-axis.
	 */
	double GetPitching_Ampl_X(unsigned short val_iZone);

	/*!
	 * \brief Get the pitching amplitude about the y-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Pitching amplitude about the y-axis.
	 */
	double GetPitching_Ampl_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the pitching amplitude about the z-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Pitching amplitude about the z-axis.
	 */
	double GetPitching_Ampl_Z(unsigned short val_iZone);

	/*!
	 * \brief Get the pitching phase offset about the x-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Pitching phase offset about the x-axis.
	 */
	double GetPitching_Phase_X(unsigned short val_iZone);

	/*!
	 * \brief Get the pitching phase offset about the y-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Pitching phase offset about the y-axis.
	 */
	double GetPitching_Phase_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the pitching phase offset about the z-axis.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Pitching phase offset about the z-axis.
	 */
	double GetPitching_Phase_Z(unsigned short val_iZone);

	/*!
	 * \brief Get the angular frequency of a mesh plunging in the x-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular frequency of a mesh plunging in the x-direction.
	 */
	double GetPlunging_Omega_X(unsigned short val_iZone);

	/*!
	 * \brief Get the angular frequency of a mesh plunging in the y-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular frequency of a mesh plunging in the y-direction.
	 */
	double GetPlunging_Omega_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the angular frequency of a mesh plunging in the z-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Angular frequency of a mesh plunging in the z-direction.
	 */
	double GetPlunging_Omega_Z(unsigned short val_iZone);

	/*!
	 * \brief Get the plunging amplitude in the x-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Plunging amplitude in the x-direction.
	 */
	double GetPlunging_Ampl_X(unsigned short val_iZone);

	/*!
	 * \brief Get the plunging amplitude in the y-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Plunging amplitude in the y-direction.
	 */
	double GetPlunging_Ampl_Y(unsigned short val_iZone);

	/*!
	 * \brief Get the plunging amplitude in the z-direction.
	 * \param[in] val_iZone - Number for the current zone in the mesh (each zone has independent motion).
	 * \return Plunging amplitude in the z-direction.
	 */
	double GetPlunging_Ampl_Z(unsigned short val_iZone);

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
	double GetminTurkelBeta();

	/*!
	 * \brief Get the minimum value of Beta for Roe-Turkel preconditioner
	 * \return the minimum value of Beta for Roe-Turkel preconditioner
	 */
	double GetmaxTurkelBeta();

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
	 * \brief Get information about the catalytic wall condition
	 * \return <code>TRUE</code> if it is a catalytic wall condition; otherwise <code>FALSE</code>.
	 */
	bool GetCatalyticWall(void);

	/*!
	 * \brief Get information about the Low Mach Preconditioning
	 * \return <code>TRUE</code> if we are using low Mach preconditioner; otherwise <code>FALSE</code>.
	 */
	bool Low_Mach_Preconditioning(void);

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
	 * \brief Obtain the kind of convergence criteria to establish the convergence of the CFD code.
	 * \return Kind of convergence criteria.
	 */
	unsigned short GetConvCriteria(void);

	/*!
	 * \brief This subroutine adds the domain index to the name of some input-output file names.
	 * \param[in] val_domain - Index of the domain.
	 */
	void SetFileNameDomain(unsigned short val_domain);

	/*!
	 * \brief Get the index in the config information of the marker <i>val_marker</i>.
	 * \note When we read the config file, it stores the markers in a particular vector.
	 * \return Index in the config information of the marker <i>val_marker</i>.
	 */
	unsigned short GetMarker_CfgFile_TagBound(string val_marker);

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
   * \brief Get the 1-D output (ie, averaged pressure) information from the config definition for the marker <i>val_marker</i>.
   * \return 1D output information of the boundary in the config information for the marker <i>val_marker</i>.
   */
  unsigned short GetMarker_CfgFile_Out_1D(string val_marker);

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
	 * \brief Get the periodic information from the config definition of the marker <i>val_marker</i>.
	 * \return Periodic information of the boundary in the config information of the marker <i>val_marker</i>.
	 */
	unsigned short GetMarker_CfgFile_PerBound(string val_marker);

	/*!
	 * \brief Determines if problem is adjoint
	 * \return true if Adjoint
	 */
	bool GetAdjoint(void);

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
	double GetOrderMagResidual(void);

	/*!
	 * \brief Value of the minimum residual value (log10 scale).
	 * \return Value of the minimum residual value (log10 scale).
	 */
	double GetMinLogResidual(void);

    /*!
	 * \brief Value of the damping factor for the engine inlet bc.
	 * \return Value of the damping factor.
	 */
	double GetDamp_Engine_Inflow(void);

	/*!
	 * \brief Value of the damping factor for the residual restriction.
	 * \return Value of the damping factor.
	 */
	double GetDamp_Res_Restric(void);

	/*!
	 * \brief Value of the damping factor for the correction prolongation.
	 * \return Value of the damping factor.
	 */
	double GetDamp_Correc_Prolong(void);

	/*!
	 * \brief Value of the position of the Near Field (y coordinate for 2D, and z coordinate for 3D).
	 * \return Value of the Near Field position.
	 */
	double GetPosition_Plane(void);

	/*!
	 * \brief Value of the weight of the drag coefficient in the Sonic Boom optimization.
	 * \return Value of the weight of the drag coefficient in the Sonic Boom optimization.
	 */
	double GetWeightCd(void);

  /*!
	 * \brief Value of the azimuthal line to fix due to a misalignments of the nearfield.
	 * \return Azimuthal line to fix due to a misalignments of the nearfield.
	 */
	double GetFixAzimuthalLine(void);

	/*!
	 * \brief Update the CFL number using the ramp information.
	 * \param[in] val_iter - Current solver iteration.
	 */
	void UpdateCFL(unsigned long val_iter);

	/*!
	 * \brief Set the global parameters of each simulation for each runtime system.
	 * \param[in] val_solver - Solver of the simulation.
	 * \param[in] val_system - Runtime system that we are solving.
	 */
	void SetGlobalParam(unsigned short val_solver, unsigned short val_system, unsigned long val_extiter);

	/*!
	 * \brief Center of rotation for a rotational periodic boundary.
	 */
	double *GetPeriodicRotCenter(string val_marker);

	/*!
	 * \brief Angles of rotation for a rotational periodic boundary.
	 */
	double *GetPeriodicRotAngles(string val_marker);

	/*!
	 * \brief Translation vector for a rotational periodic boundary.
	 */
	double *GetPeriodicTranslation(string val_marker);

	/*!
	 * \brief Get the rotationally periodic donor marker for boundary <i>val_marker</i>.
	 * \return Periodic donor marker from the config information for the marker <i>val_marker</i>.
	 */
	unsigned short GetMarker_Periodic_Donor(string val_marker);

  /*!
	 * \brief Get the origin of the actuator disk.
	 */
  double* GetActDisk_Origin(string val_marker);

  /*!
	 * \brief Get the root radius of the actuator disk.
	 */
  double GetActDisk_RootRadius(string val_marker);

  /*!
	 * \brief Get the tip radius of th actuator disk.
	 */
  double GetActDisk_TipRadius(string val_marker);

  /*!
	 * \brief Get the thurst corffient of the actuator disk.
	 */
  double GetActDisk_CT(string val_marker);

  /*!
	 * \brief Get the rev / min of the actuator disk.
	 */
  double GetActDisk_Omega(string val_marker);

  /*!
	 * \brief Get Actuator Disk Outlet for boundary <i>val_marker</i> (actuator disk inlet).
	 * \return Actuator Disk Outlet from the config information for the marker <i>val_marker</i>.
	 */
	unsigned short GetMarker_ActDisk_Outlet(string val_marker);

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
	string GetMarker_Moving(unsigned short val_marker);

	/*!
	 * \brief Get information about converting a mesh from CGNS to SU2 format.
	 * \return <code>TRUE</code> if a conversion is requested; otherwise <code>FALSE</code>.
	 */
	bool GetCGNS_To_SU2(void);

	/*!
	 * \brief Get information about whether a converted mesh should be written.
	 * \return <code>TRUE</code> if the converted mesh should be written; otherwise <code>FALSE</code>.
	 */
	bool GetMesh_Output(void);

	/*!
	 * \brief Set the total number of SEND_RECEIVE periodic transformations.
	 * \param[in] val_index - Total number of transformations.
	 */
	void SetnPeriodicIndex(unsigned short val_index);

	/*!
	 * \brief Get the total number of SEND_RECEIVE periodic transformations.
	 * \return Total number of transformations.
	 */
	unsigned short GetnPeriodicIndex(void);

	/*!
	 * \brief Set the rotation center for a periodic transformation.
	 * \param[in] val_index - Index corresponding to the periodic transformation.
	 * \param[in] center - Pointer to a vector containing the coordinate of the center.
	 */
	void SetPeriodicCenter(unsigned short val_index, double* center);

	/*!
	 * \brief Get the rotation center for a periodic transformation.
	 * \param[in] val_index - Index corresponding to the periodic transformation.
	 * \return A vector containing coordinates of the center point.
	 */
	double* GetPeriodicCenter(unsigned short val_index);

	/*!
	 * \brief Set the rotation angles for a periodic transformation.
	 * \param[in] val_index - Index corresponding to the periodic transformation.
	 * \param[in] rotation - Pointer to a vector containing the rotation angles.
	 */
	void SetPeriodicRotation(unsigned short val_index, double* rotation);

	/*!
	 * \brief Get the rotation angles for a periodic transformation.
	 * \param[in] val_index - Index corresponding to the periodic transformation.
	 * \return A vector containing the angles of rotation.
	 */
	double* GetPeriodicRotation(unsigned short val_index);

	/*!
	 * \brief Set the translation vector for a periodic transformation.
	 * \param[in] val_index - Index corresponding to the periodic transformation.
	 * \param[in] translate - Pointer to a vector containing the coordinate of the center.
	 */
	void SetPeriodicTranslate(unsigned short val_index, double* translate);

	/*!
	 * \brief Get the translation vector for a periodic transformation.
	 * \param[in] val_index - Index corresponding to the periodic transformation.
	 * \return The translation vector.
	 */
	double* GetPeriodicTranslate(unsigned short val_index);

	/*!
	 * \brief Get the total temperature at a nacelle boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The total temperature.
	 */
	double GetNozzle_Ttotal(string val_index);

	/*!
	 * \brief Get the total temperature at an inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The total temperature.
	 */
	double GetInlet_Ttotal(string val_index);

	/*!
	 * \brief Get the temperature at a supersonic inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The inlet density.
	 */
	double GetInlet_Temperature(string val_index);

	/*!
	 * \brief Get the pressure at a supersonic inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The inlet pressure.
	 */
	double GetInlet_Pressure(string val_index);

	/*!
	 * \brief Get the velocity vector at a supersonic inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The inlet velocity vector.
	 */
	double* GetInlet_Velocity(string val_index);

	/*!
	 * \brief Get the fixed value at the Dirichlet boundary.
	 * \param[in] val_index - Index corresponding to the Dirichlet boundary.
	 * \return The total temperature.
	 */
	double GetDirichlet_Value(string val_index);

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
	double GetInlet_Ptotal(string val_index);

	/*!
	 * \brief Get the total pressure at an nacelle boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The total pressure.
	 */
	double GetNozzle_Ptotal(string val_index);

	/*!
	 * \brief If inlet and outlet conditions are defined for multi species
	 * \return true/false
	 */
	bool GetInletConditionsDefined();

	/*!
	 * \brief Get the temperature at an inlet boundary.
	 * \param[in] iSpecies - Index of the species
	 * \return The total temperature.
	 */
	double GetInlet_Species_Temperature(unsigned short iSpecies);

	/*!
	 * \brief Get the temperature at an outlet boundary.
	 * \param[in] iSpecies - Index of the species
	 * \return The total temperature.
	 */
	double GetOutlet_Species_Temperature(unsigned short iSpecies);

	/*!
	 * \brief Get the pressure at an inlet boundary.
	 * \param[in] iSpecies - Index of the species
	 * \return The total temperature.
	 */
	double GetInlet_Species_Pressure(unsigned short iSpecies);

	/*!
	 * \brief Get the pressure at an outlet boundary.
	 * \param[in] iSpecies - Index of the species
	 * \return The total temperature.
	 */
	double GetOutlet_Species_Pressure(unsigned short iSpecies);

	/*!
	 * \brief Get the velocity at an inlet boundary.
	 * \param[in] iSpecies - Index of the species
	 * \return The total temperature.
	 */
	double GetInlet_Species_Velocity(unsigned short iSpecies);

	/*!
	 * \brief Get the velocity at an outlet boundary.
	 * \param[in] iSpecies - Index of the species
	 * \return The total temperature.
	 */
	double GetOutlet_Species_Velocity(unsigned short iSpecies);

  /*!
	 * \brief Value of the CFL reduction in LevelSet problems.
	 * \return Value of the CFL reduction in LevelSet problems.
	 */
	double GetCFLRedCoeff_Turb(void);

	/*!
	 * \brief Get the flow direction unit vector at an inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The flow direction vector.
	 */
	double* GetInlet_FlowDir(string val_index);

	/*!
	 * \brief Get the back pressure (static) at an outlet boundary.
	 * \param[in] val_index - Index corresponding to the outlet boundary.
	 * \return The outlet pressure.
	 */
	double GetOutlet_Pressure(string val_index);

	/*!
	 * \brief Get the var 1 at Riemann boundary.
	 * \param[in] val_marker - Index corresponding to the Riemann boundary.
	 * \return The var1
	 */
	double GetRiemann_Var1(string val_marker);

	/*!
	 * \brief Get the var 2 at Riemann boundary.
	 * \param[in] val_marker - Index corresponding to the Riemann boundary.
	 * \return The var2
	 */

	double GetRiemann_Var2(string val_marker);

	/*!
	 * \brief Get the Flowdir at Riemann boundary.
	 * \param[in] val_marker - Index corresponding to the Riemann boundary.
	 * \return The Flowdir
	 */
	double* GetRiemann_FlowDir(string val_marker);

	/*!
	 * \brief Get Kind Data of Riemann boundary.
	 * \param[in] val_marker - Index corresponding to the Riemann boundary.
	 * \return Kind data
	 */
	unsigned short GetKind_Data_Riemann(string val_marker);

	/*!
	 * \brief Get the wall temperature (static) at an isothermal boundary.
	 * \param[in] val_index - Index corresponding to the isothermal boundary.
	 * \return The wall temperature.
	 */
	double GetIsothermal_Temperature(string val_index);

	/*!
	 * \brief Get the wall heat flux on a constant heat flux boundary.
	 * \param[in] val_index - Index corresponding to the constant heat flux boundary.
	 * \return The heat flux.
	 */
	double GetWall_HeatFlux(string val_index);

	/*!
	 * \brief Get the wall heat flux on a constant heat flux boundary.
	 * \return The heat flux.
	 */
	double *GetWall_Catalycity(void);

	/*!
	 * \brief Get the back pressure (static) at an outlet boundary.
	 * \param[in] val_index - Index corresponding to the outlet boundary.
	 * \return The outlet pressure.
	 */
	double GetInflow_Mach_Target(string val_marker);

    /*!
	 * \brief Get the back pressure (static) at an outlet boundary.
	 * \param[in] val_index - Index corresponding to the outlet boundary.
	 * \return The outlet pressure.
	 */
	double GetInflow_Mach(string val_marker);

    /*!
	 * \brief Get the back pressure (static) at an outlet boundary.
	 * \param[in] val_index - Index corresponding to the outlet boundary.
	 * \return The outlet pressure.
	 */
	void SetInflow_Mach(unsigned short val_imarker, double val_fanface_mach);

    /*!
	 * \brief Get the back pressure (static) at an outlet boundary.
	 * \param[in] val_index - Index corresponding to the outlet boundary.
	 * \return The outlet pressure.
	 */
	double GetInflow_Pressure(string val_marker);

    /*!
	 * \brief Get the back pressure (static) at an outlet boundary.
	 * \param[in] val_index - Index corresponding to the outlet boundary.
	 * \return The outlet pressure.
	 */
	void SetInflow_Pressure(unsigned short val_imarker, double val_fanface_pressure);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  double GetBleed_Temperature_Target(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  double GetBleed_Temperature(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetBleed_Temperature(unsigned short val_imarker, double val_bleed_temp);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  double GetBleed_MassFlow_Target(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  double GetBleed_MassFlow(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetBleed_MassFlow(unsigned short val_imarker, double val_bleed_massflow);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  double GetBleed_Pressure(string val_marker);
  
  /*!
   * \brief Get the back pressure (static) at an outlet boundary.
   * \param[in] val_index - Index corresponding to the outlet boundary.
   * \return The outlet pressure.
   */
  void SetBleed_Pressure(unsigned short val_imarker, double val_fanface_pressure);
  
	/*!
	 * \brief Get the displacement value at an displacement boundary.
	 * \param[in] val_index - Index corresponding to the displacement boundary.
	 * \return The displacement value.
	 */
	double GetDispl_Value(string val_index);

	/*!
	 * \brief Get the force value at an load boundary.
	 * \param[in] val_index - Index corresponding to the load boundary.
	 * \return The load value.
	 */
	double GetLoad_Value(string val_index);

	/*!
	 * \brief Get the force value at an load boundary.
	 * \param[in] val_index - Index corresponding to the load boundary.
	 * \return The load value.
	 */
	double GetFlowLoad_Value(string val_index);

	/*!
	 * \brief Cyclic pitch amplitude for rotor blades.
	 * \return The specified cyclic pitch amplitude.
	 */
	double GetCyclic_Pitch(void);

	/*!
	 * \brief Collective pitch setting for rotor blades.
	 * \return The specified collective pitch setting.
	 */
	double GetCollective_Pitch(void);

	/*!
	 * \brief Get name of the arbitrary mesh motion input file.
	 * \return File name of the arbitrary mesh motion input file.
	 */
	string GetMotion_FileName(void);

  /*!
	 * \brief Set the config options.
	 */
	void SetConfig_Options(unsigned short val_iZone, unsigned short val_nZone);

  /*!
	 * \brief Set the config file parsing.
	 */
  void SetParsing(char case_filename[MAX_STRING_SIZE]);

	/*!
	 * \brief Config file postprocessing.
	 */
	void SetPostprocessing(unsigned short val_software, unsigned short val_izone, unsigned short val_nDim);

	/*!
	 * \brief Config file markers processing.
	 */
	void SetMarkers(unsigned short val_software, unsigned short val_izone);

	/*!
	 * \brief Config file output.
	 */
	void SetOutput(unsigned short val_software, unsigned short val_izone);

	/*!
	 * \brief Value of Aeroelastic solution coordinate at time n+1.
	 */
	vector<vector<double> > GetAeroelastic_np1(unsigned short iMarker);

	/*!
	 * \brief Value of Aeroelastic solution coordinate at time n.
	 */
	vector<vector<double> > GetAeroelastic_n(unsigned short iMarker);

	/*!
	 * \brief Value of Aeroelastic solution coordinate at time n-1.
	 */
	vector<vector<double> > GetAeroelastic_n1(unsigned short iMarker);

	/*!
	 * \brief Value of Aeroelastic solution coordinate at time n+1.
	 */
	void SetAeroelastic_np1(unsigned short iMarker, vector<vector<double> > solution);

	/*!
	 * \brief Value of Aeroelastic solution coordinate at time n from time n+1.
	 */
	void SetAeroelastic_n(void);

	/*!
	 * \brief Value of Aeroelastic solution coordinate at time n-1 from time n.
	 */
	void SetAeroelastic_n1(void);

	/*!
	 * \brief Uncoupled Aeroelastic Frequency Plunge.
	 */
	double GetAeroelastic_Frequency_Plunge(void);

	/*!
	 * \brief Uncoupled Aeroelastic Frequency Pitch.
	 */
	double GetAeroelastic_Frequency_Pitch(void);

	/*!
	 * \brief Value of plunging coordinate.
     * \param[in] val_marker - the marker we are monitoring.
	 * \return Value of plunging coordinate.
	 */
	double GetAeroelastic_plunge(unsigned short val_marker);

    /*!
	 * \brief Value of pitching coordinate.
     * \param[in] val_marker - the marker we are monitoring.
	 * \return Value of pitching coordinate.
	 */
	double GetAeroelastic_pitch(unsigned short val_marker);

	/*!
	 * \brief Value of plunging coordinate.
     * \param[in] val_marker - the marker we are monitoring.
     * \param[in] val - value of plunging coordinate.
	 */
	void SetAeroelastic_plunge(unsigned short val_marker, double val);

	/*!
	 * \brief Value of pitching coordinate.
     * \param[in] val_marker - the marker we are monitoring.
     * \param[in] val - value of pitching coordinate.
	 */
	void SetAeroelastic_pitch(unsigned short val_marker, double val);

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
	double GetGust_WaveLength(void);

    /*!
	 * \brief Value of the number of gust periods.
	 */
	double GetGust_Periods(void);

    /*!
	 * \brief Value of the gust amplitude.
	 */
	double GetGust_Ampl(void);

    /*!
	 * \brief Value of the time at which to begin the gust.
	 */
	double GetGust_Begin_Time(void);

    /*!
	 * \brief Value of the location ath which the gust begins.
	 */
	double GetGust_Begin_Loc(void);

  /*!
	 * \brief Value of the time at which to begin the gust.
	 */
	unsigned short GetnFFD_Iter(void);

  /*!
	 * \brief Value of the location ath which the gust begins.
	 */
	double GetFFD_Tol(void);

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
	 * \brief Get the value specified for the target CL.
	 * \return Value of the target CL.
	 */
	double GetTarget_CL(void);

  /*!
	 * \brief Get the value of the damping coefficient for fixed CL mode.
	 * \return Damping coefficient for fixed CL mode.
	 */
	double GetDamp_Fixed_CL(void);

  /*!
	 * \brief Set the value of the boolean for updating AoA in fixed lift mode.
   * \param[in] val_update - the bool for whether to update the AoA.
	 */
	void SetUpdate_AoA(bool val_update);

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
	void SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2);

	/*!
	 * \brief Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
	          and given the array y2a[1..n], which is the output from spline above, and given a value of
	          x, this routine returns a cubic-spline interpolated value y.
         	  Numerical Recipes: The Art of Scientific Computing, Third Edition in C++.
	 * \returns The interpolated value of for x.
	 */
	double GetSpline(vector<double> &xa, vector<double> &ya, vector<double> &y2a, unsigned long n, double x);

};

#include "config_structure.inl"
