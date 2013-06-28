/*!
 * \file config_structure.hpp
 * \brief All the information about the definition of the physical problem.
 *        The subroutines and functions are in the <i>config_structure.cpp</i> file.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h> 
#include <cmath>
#include <map>
#include "./option_structure.hpp"

#ifndef NO_MPI
#include <mpi.h>
#endif

using namespace std;

/*! 
 * \class CConfig
 * \brief Main class for defining the problem; basically this class reads the configuration file, and
 *        stores all the information.
 * \author F. Palacios.
 * \version 1.1.
 */
class CConfig {
private:
	unsigned short Kind_SU2; /*!< \brief Kind of SU2 software component. */
	double OrderMagResidual; /*!< \brief Order of magnitude reduction. */
	double *RotAxisOrigin,	 /*!< \brief Axis of rotation (origin) for rotational frame problem. */
	*Omega,						/*!< \brief Angular velocity vector for rotational frame problem. */
	Omega_Mag,						/*!< \brief Angular velocity magnitude for rotational frame problem. */
	Rot_Radius;						/*!< \brief Reference length (e.g. rotor radius) for computing force coefficients in a rotating frame. */
	double MinLogResidual; /*!< \brief Minimum value of the log residual. */
	double* EA_IntLimit; /*!< \brief Integration limits of the Equivalent Area computation */
	unsigned short ConvCriteria;	/*!< \brief Kind of convergence criteria. */
	bool ContAdj,			/*!< \brief Flag to know if the code is solving a continuous adjoint problem. */
	EquivArea,				/*!< \brief Flag to know if the code is going to compute and plot the equivalent area. */
	OneShot,				/*!< \brief Flag to know if the code is solving a one shot problem. */
	Linerized,				/*!< \brief Flag to know if the code is solving a linearized problem. */
	BiGrid,					/*!< \brief BiGrid filtering strategy for sensitivity computation. */
	Grid_Movement,			/*!< \brief Flag to know if there is grid movement. */
	Rotating_Frame,			/*!< \brief Flag to know if there is a rotating frame. */
	Incompressible,			/*!< \brief Flag to know if we are using the incompressible formulation. */
	AdiabaticWall,			/*!< \brief Flag to know if we are using the Adiabatic Wall. */
	IsothermalWall,			/*!< \brief Flag to know if we are using the Isothermal Wall. */
	CatalyticWall,			/*!< \brief Flag to know if we are using the Catalytic Wall. */
	PlasmaMultiTimeSteps,	/*!< \brief Flag to know if we are using multiple time steps for different species in plasma. */
	ElectricSolver,		/*!< \brief Flag to know if we are solving  electric forces  in plasma solver. */
	FreeSurface,			/*!< \brief Flag to know if we are solving a free surface problem. */
	GravityForce,			/*!< \brief Flag to know if the gravity force is incuded in the formulation. */
	SmoothNumGrid,			/*!< \brief Smooth the numerical grid. */
	FullMG,					/*!< \brief Full multigrid strategy. */
	Divide_Element,			/*!< \brief Divide rectables and hexahedrom. */
	Frozen_Visc,			/*!< \brief Flag for adjoint problem with/without frozen viscosity. */
	Move_FFDBox,	/*!< \brief Flag for rotational frame. */
	Axisymmetric; /*!< \brief Flag for axisymmetric calculations */
	bool ReStartMGCycle;			/*!< \brief Flag to know if MG cycle must be restarted with the interpolated solution. */
	bool ShockTube; /*!< \brief Flag to define a shock tube problem. */
	bool Visualize_Partition;	/*!< \brief Flag to visualize each partition in the DDM. */
	bool Visualize_Deformation;	/*!< \brief Flag to visualize the deformation in the MDC. */
	double Damp_Res_Restric,	/*!< \brief Damping factor for the residual restriction. */
	Damp_Correc_Prolong; /*!< \brief Damping factor for the correction prolongation. */
	double Position_Plane; /*!< \brief Position of the Near-Field (y coordinate 2D, and z coordinate 3D). */
	double WeightCd; /*!< \brief Weight of the drag coefficient. */
	unsigned short Unsteady_Simulation;	/*!< \brief Steady or unsteady (time stepping or dual time stepping) computation. */
	unsigned short nStartUpIter;	/*!< \brief Start up iterations using the fine grid. */
	double CteViscDrag;		/*!< \brief Constant value of the viscous drag. */
	double *DV_Value_New,		/*!< \brief Finite difference step for gradient computation. */
	*DV_Value_Old;		/*!< \brief Previous value of the design variable. */
	unsigned short nSmooth;			/*!< \brief Number of iterations of the residual smoothing technique. */
	double SmoothCoeff;				/*!< \brief Relaxation factor of the residual smoothing technique */ 
	double LimiterCoeff;				/*!< \brief Limiter coefficient */ 
	unsigned short Kind_ObjFunc;	/*!< \brief Kind of objective function. */
	unsigned short *Design_Variable; /*!< \brief Kind of design variable. */
	double RatioDensity,				/*!< \brief Ratio of density for a free surface problem. */
	RatioViscosity,				/*!< \brief Ratio of viscosity for a free surface problem. */
	FreeSurface_Thickness,  /*!< \brief Thickness of the interfase for a free surface problem. */
	FreeSurface_Outlet,  /*!< \brief Outlet of the interfase for a free surface problem. */
	FreeSurface_Inlet,  /*!< \brief Inlet of the interfase for a free surface problem. */
	FreeSurface_Damping_Coeff,  /*!< \brief Damping coefficient of the free surface for a free surface problem. */
	FreeSurface_Damping_Length;  /*!< \brief Damping length of the free surface for a free surface problem. */
	unsigned short Kind_Adaptation;	/*!< \brief Kind of numerical grid adaptation. */
	double New_Elem_Adapt;			/*!< \brief Elements to adapt in the numerical grid adaptation process. */
	double Delta_UnstTime,			/*!< \brief Time step for unsteady computations. */
	Delta_UnstTimeND;						/*!< \brief Time step for unsteady computations (non dimensional). */
	double Unst_Time,						/*!< \brief Total time for unsteady computations. */
	Unst_TimeND;								/*!< \brief Total time for unsteady computations (non dimensional). */
	bool Unst_Mesh_Motion;			/*!< \brief Flag for including mesh motion in unsteady computations. */
	double Reduced_Frequency,		/*!< \brief Reduced frequency for airfoil movement. */
	Pitching_Amplitude;				/*!< \brief Pitching amplitude for airfoil movement. */
	unsigned short nMarker_Euler,	/*!< \brief Number of Euler wall markers. */
	nMarker_NS,						/*!< \brief Number of no slip wall markers. */
	nMarker_FarField,				/*!< \brief Number of far-field markers. */
	nMarker_Custom,
	nMarker_SymWall,				/*!< \brief Number of symmetry wall markers. */
	nMarker_PerBound,				/*!< \brief Number of periodic boundary markers. */
	nMarker_NearFieldBound,				/*!< \brief Number of near field boundary markers. */
	nMarker_InterfaceBound,				/*!< \brief Number of interface boundary markers. */
	nMarker_Dirichlet,				/*!< \brief Number of interface boundary markers. */
	nMarker_Inlet,					/*!< \brief Number of inlet flow markers. */
	nMarker_Outlet,					/*!< \brief Number of outlet flow markers. */
	nMarker_Neumann,				/*!< \brief Number of Neumann flow markers. */
	nMarker_All,					/*!< \brief Total number of markers using the grid information. */
	nMarker_Config;					/*!< \brief Total number of markers using the config file 
									(note that using parallel computation this number can be different 
									from nMarker_All). */
	string *Marker_Euler,			/*!< \brief Euler wall markers. */
	*Marker_NS,						/*!< \brief No slip wall markers. */
	*Marker_FarField,				/*!< \brief Far field markers. */
	*Marker_Custom,
	*Marker_SymWall,				/*!< \brief Symmetry wall markers. */
	*Marker_PerBound,				/*!< \brief Periodic boundaries markers. */
	*Marker_PerDonor,		/*!< \brief Rotationally periodic boundary donor markers. */
	*Marker_NearFieldBound,				/*!< \brief Near Field boundaries markers. */
	*Marker_InterfaceBound,				/*!< \brief Interface boundaries markers. */
	*Marker_Dirichlet,				/*!< \brief Interface boundaries markers. */
	*Marker_Inlet,					/*!< \brief Inlet flow markers. */
	*Marker_Outlet,					/*!< \brief Outlet flow markers. */
	*Marker_Neumann,					/*!< \brief Neumann flow markers. */
	*Marker_Electrode,				/*!< \brief Electrode flow markers. */
	*Marker_Dielectric,				/*!< \brief Dielectric flow markers. */
	*Marker_All_Tag;				/*!< \brief Global index for markers using grid information. */
	double *Inlet_Ttotal;    /*!< \brief Specified total temperatures for inlet boundaries. */
	double *Inlet_Ptotal;    /*!< \brief Specified total pressures for inlet boundaries. */
	double **Inlet_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
	double *Outlet_Pressure;    /*!< \brief Specified back pressures (static) for outlet boundaries. */
	double **Periodic_RotCenter;  /*!< \brief Rotational center for each periodic boundary. */
	double **Periodic_RotAngles;      /*!< \brief Rotation angles for each periodic boundary. */
	double **Periodic_Translation;      /*!< \brief Translation vector for each periodic boundary. */
	unsigned short nPeriodic_Index;     /*!< \brief Number of SEND_RECEIVE periodic transformations. */
	double **Periodic_Center;         /*!< \brief Rotational center for each SEND_RECEIVE boundary. */
	double **Periodic_Rotation;      /*!< \brief Rotation angles for each SEND_RECEIVE boundary. */
	double **Periodic_Translate;      /*!< \brief Translation vector for each SEND_RECEIVE boundary. */
	string *Marker_Config_Tag;			/*!< \brief Global index for markers using config file. */
	unsigned short *Marker_All_Boundary,			/*!< \brief Global index for boundaries using grid information. */
	*Marker_Config_Boundary;		/*!< \brief Global index for boundaries using config file. */
	short *Marker_All_SendRecv;		/*!< \brief Information about if the boundary is sended (+), received (-). */
	short *Marker_All_PerBound;	/*!< \brief Global index for periodic bc using the grid information. */
	unsigned long nExtIter;			/*!< \brief Number of external iterations. */
	unsigned long Unst_nIntIter;			/*!< \brief Number of internal iterations (Dual time Method). */
	unsigned short nRKStep;			/*!< \brief Number of steps of the explicit Runge-Kutta method. */
	double *RK_Alpha_Step,	/*!< \brief Runge-Kutta alfa coefficients. */
	*RK_Beta_Step;			/*!< \brief Runge-Kutta beta coefficients. */
	unsigned short nMultiLevel;		/*!< \brief Number of multigrid levels (coarse levels). */
	unsigned short nCFL;			/*!< \brief Number of CFL, one for each multigrid level. */
	double *CFL,		/*!< \brief CFL number for each multigrid level. */
	MG_CFLRedCoeff,		/*!< \brief CFL reduction coefficient on the MG coarse level. */
	LevelSet_CFLRedCoeff,		/*!< \brief CFL reduction coefficient on the LevelSet problem. */
	Adj_CFLRedCoeff,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
	CFLFineGrid,		/*!< \brief CFL of the finest grid. */
	Unst_CFL;		/*!< \brief Unsteady CFL number. */
	unsigned short MaxChildren;		/*!< \brief Maximum number of children. */
	double MaxDimension;			/*!< \brief Maximum dimension of the aglomerated element compared with the whole domain. */
	bool AddIndNeighbor;			/*!< \brief Include indirect neighbor in the agglomeration process. */
	unsigned short nDV;		/*!< \brief Number of design variables. */
	unsigned short nParamDV;		/*!< \brief Number of parameters of the design variable. */
	double **ParamDV;				/*!< \brief Parameters of the design variable. */
	double PrimGrad_Threshold;		/*!< \brief Primitive variables gradient threshold for the adjoint problem. */
	unsigned short MGCycle;			/*!< \brief Kind of multigrid cycle. */
	unsigned short FinestMesh;		/*!< \brief Finest mesh for the full multigrid approach. */
	unsigned short nMG_PreSmooth,                 /*!< \brief Number of MG pre-smooth parameters found in config file. */
	nMG_PostSmooth,                             /*!< \brief Number of MG post-smooth parameters found in config file. */
	nMG_CorrecSmooth;                           /*!< \brief Number of MG correct-smooth parameters found in config file. */
	unsigned short *MG_PreSmooth,	/*!< \brief Multigrid Pre smoothing. */
	*MG_PostSmooth,					/*!< \brief Multigrid Post smoothing. */
	*MG_CorrecSmooth;					/*!< \brief Multigrid Jacobi implicit smoothing of the correction. */
	unsigned short Kind_Solver,	/*!< \brief Kind of solver Euler, NS, Continuous adjoint, etc. */
	Kind_GasModel,				/*!< \brief Kind of the Gas Model. */
	Kind_Solver_MPI,			/*!< \brief Copy of Kind_Solver. */
	Kind_Gradient_Method,		/*!< \brief Numerical method for computation of spatial gradients. */
	Kind_GridDef_Method,		/*!< \brief Numerical method for the grid deformation. */
	Kind_GridDef_Solver,		/*!< \brief Numerical solver for the grid deformation (linear system solving). */
	Kind_Linear_Solver,		/*!< \brief Numerical solver for the implicit scheme. */
	Kind_SlopeLimit,				/*!< \brief Global slope limiter. */
	Kind_SlopeLimit_Flow,		/*!< \brief Slope limiter for flow equations.*/
	Kind_SlopeLimit_Turb,		/*!< \brief Slope limiter for the turbulence equation.*/
	Kind_SlopeLimit_LevelSet,		/*!< \brief Slope limiter for the level set equation.*/
	Kind_SlopeLimit_AdjLevelSet,		/*!< \brief Slope limiter for the adjoint level set equation.*/
	Kind_SlopeLimit_Plasma,		/*!< \brief Slope limiter for the plasma equation.*/
	Kind_SlopeLimit_AdjTurb,	/*!< \brief Slope limiter for the adjoint turbulent equation.*/
	Kind_SlopeLimit_AdjFlow,	/*!< \brief Slope limiter for the adjoint equation.*/
	Kind_TimeNumScheme,			/*!< \brief Global explicit or implicit time integration. */
	Kind_TimeIntScheme_Flow,	/*!< \brief Time integration for the flow equations. */
	Kind_TimeIntScheme_AdjFlow,		/*!< \brief Time integration for the adjoint flow equations. */
	Kind_TimeIntScheme_LinFlow,		/*!< \brief Time integration for the linearized flow equations. */
	Kind_TimeIntScheme_Turb,	/*!< \brief Time integration for the turbulence model. */
	Kind_TimeIntScheme_LevelSet,	/*!< \brief Time integration for the level set model. */
	Kind_TimeIntScheme_AdjLevelSet,	/*!< \brief Time integration for the adjoint level set model. */
	Kind_TimeIntScheme_Combustion,	/*!< \brief Time integration for the turbulence model. */
	Kind_TimeIntScheme_AdjTurb,	/*!< \brief Time integration for the adjoint turbulence model. */
	Kind_TimeIntScheme_Plasma,	/*!< \brief Time integration for the plasma equations. */
	Kind_TimeIntScheme_Wave,	/*!< \brief Time integration for the wave equations. */
	Kind_ConvNumScheme,			/*!< \brief Global definition of the convective term. */
	Kind_ConvNumScheme_Flow,	/*!< \brief Centered or upwind scheme for the flow equations. */
	Kind_ConvNumScheme_AdjFlow,		/*!< \brief Centered or upwind scheme for the adjoint flow equations. */
	Kind_ConvNumScheme_LinFlow,		/*!< \brief Centered or upwind scheme for the linearized flow equations. */
	Kind_ConvNumScheme_Turb,	/*!< \brief Centered or upwind scheme for the turbulence model. */
	Kind_ConvNumScheme_AdjTurb,	/*!< \brief Centered or upwind scheme for the adjoint turbulence model. */
	Kind_ConvNumScheme_Plasma,	/*!< \brief Centered or upwind scheme for the plasma equations. */
	Kind_ConvNumScheme_LevelSet,	/*!< \brief Centered or upwind scheme for the level set equation. */
	Kind_ConvNumScheme_AdjLevelSet,	/*!< \brief Centered or upwind scheme for the adjoint level set equation. */
	Kind_ConvNumScheme_Combustion,	/*!< \brief Centered or upwind scheme for the level set equation. */
	Kind_ConvNumScheme_Template,	/*!< \brief Centered or upwind scheme for the level set equation. */
	Kind_ViscNumScheme,			/*!< \brief Global definition of the viscous term. */
	Kind_ViscNumScheme_Flow,	/*!< \brief Viscous scheme for the flow equations. */
	Kind_ViscNumScheme_AdjFlow,		/*!< \brief Viscous scheme for the adjoint flow equations. */
	Kind_ViscNumScheme_LinFlow,		/*!< \brief Viscous scheme for the linearized flow equations. */
	Kind_ViscNumScheme_Turb,	/*!< \brief Viscous scheme for the turbulence model. */
	Kind_ViscNumScheme_Elec,	/*!< \brief Viscous scheme for the electric potential. */
	Kind_ViscNumScheme_Wave,	/*!< \brief Viscous scheme for the wave equation. */
	Kind_ViscNumScheme_AdjTurb,	/*!< \brief Viscous scheme for the adjoint turbulence model. */
	Kind_ViscNumScheme_Plasma,	/*!< \brief Viscous scheme for the plasma equations. */
	Kind_ViscNumScheme_AdjLevelSet,	/*!< \brief Viscous scheme for the adjoint level set equation. */
	Kind_ViscNumScheme_Template,	/*!< \brief Viscous scheme for the template. */
	Kind_SourNumScheme,			/*!< \brief Global definition of the source term. */
	Kind_SourNumScheme_Flow,	/*!< \brief Source numerical scheme for the flow equations. */
	Kind_SourNumScheme_AdjFlow,		/*!< \brief Source numerical scheme for the adjoint flow equations. */
	Kind_SourNumScheme_LinFlow,		/*!< \brief Source numerical scheme for the linearized flow equations. */
	Kind_SourNumScheme_Turb,	/*!< \brief Source numerical scheme for the turbulence model. */
	Kind_SourNumScheme_Elec,	/*!< \brief Source numerical scheme for the electric potential. */
	Kind_SourNumScheme_AdjTurb,	/*!< \brief Source numerical scheme for the adjoint turbulence model. */
	Kind_SourNumScheme_AdjLevelSet,	/*!< \brief Source numerical scheme for the adjoint level set model. */
	Kind_SourNumScheme_Plasma,	/*!< \brief Source numerical scheme for the plasma equations. */
	Kind_SourNumScheme_Combustion,	/*!< \brief Source numerical scheme for the plasma equations. */
	Kind_SourNumScheme_LevelSet,	/*!< \brief Source scheme for the level set equation. */
	Kind_SourNumScheme_Wave,	/*!< \brief Source scheme for the wave equation. */
	Kind_SourNumScheme_Template,	/*!< \brief Source numerical scheme for the template. */
	Kind_Centred,				/*!< \brief Centred scheme. */
	Kind_Centred_Flow,			/*!< \brief Centred scheme for the flow equations. */
	Kind_Centred_LevelSet,			/*!< \brief Centred scheme for the level set equation. */
	Kind_Centred_AdjLevelSet,			/*!< \brief Centred scheme for the level set equation. */
	Kind_Centred_Combustion,			/*!< \brief Centred scheme for the level set equation. */
	Kind_Centred_AdjFlow,			/*!< \brief Centred scheme for the adjoint flow equations. */
	Kind_Centred_LinFlow,			/*!< \brief Centred scheme for the linearized flow equations. */
	Kind_Centred_Turb,			/*!< \brief Centred scheme for the turbulence model. */
	Kind_Centred_AdjTurb,		/*!< \brief Centred scheme for the adjoint turbulence model. */
	Kind_Centred_Plasma,		/*!< \brief Centred scheme for the adjoint plasma model. */
	Kind_Centred_Template,		/*!< \brief Centred scheme for the template model. */
	Kind_Upwind,				/*!< \brief Upwind scheme. */
	Kind_Upwind_Flow,			/*!< \brief Upwind scheme for the flow equations. */
	Kind_Upwind_LevelSet,			/*!< \brief Upwind scheme for the level set equations. */
	Kind_Upwind_AdjLevelSet,			/*!< \brief Upwind scheme for the level set equations. */
	Kind_Upwind_Combustion,			/*!< \brief Upwind scheme for the level set equations. */
	Kind_Upwind_AdjFlow,			/*!< \brief Upwind scheme for the adjoint flow equations. */
	Kind_Upwind_LinFlow,			/*!< \brief Upwind scheme for the linearized flow equations. */
	Kind_Upwind_Turb,			/*!< \brief Upwind scheme for the turbulence model. */
	Kind_Upwind_AdjTurb,		/*!< \brief Upwind scheme for the adjoint turbulence model. */
	Kind_Upwind_Template,			/*!< \brief Upwind scheme for the template model. */
	Kind_Upwind_Plasma,			/*!< \brief Upwind scheme for the plasma model. */
	Kind_Turb_Model;			/*!< \brief Turbulent model definition. */
	double Linear_Solver_Error;		/*!< \brief Min error of the linear solver for the implicit formulation. */
	unsigned short Linear_Solver_Iter;		/*!< \brief Min error of the linear solver for the implicit formulation. */
	double* Kappa_Flow,           /*!< \brief Numerical dissipation coefficients for the flow equations. */
	*Kappa_AdjFlow,                  /*!< \brief Numerical dissipation coefficients for the adjoint equations. */
	*Kappa_LinFlow;                  /*!< \brief Numerical dissipation coefficients for the linearized equations. */
	double Kappa_1st_AdjFlow,	/*!< \brief JST 1st order dissipation coefficient for adjoint flow equations (coarse multigrid levels). */
	Kappa_2nd_AdjFlow,			/*!< \brief JST 2nd order dissipation coefficient for adjoint flow equations. */
	Kappa_4th_AdjFlow,			/*!< \brief JST 4th order dissipation coefficient for adjoint flow equations. */
	Kappa_1st_LinFlow,			/*!< \brief JST 1st order dissipation coefficient for linearized flow equations (coarse multigrid levels). */
	Kappa_4th_LinFlow,			/*!< \brief JST 4th order dissipation coefficient for linearized flow equations. */
	Kappa_1st_Flow,			/*!< \brief JST 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
	Kappa_2nd_Flow,			/*!< \brief JST 2nd order dissipation coefficient for flow equations. */
	Kappa_4th_Flow;			/*!< \brief JST 4th order dissipation coefficient for flow equations. */
	double GridDef_Error;	/*!< \brief Error of the numerical solver for the linear system of the grid deformation. */
	double Mach;		/*!< \brief Mach number. */
	double Reynolds;	/*!< \brief Reynolds number. */
	double Froude;	/*!< \brief Froude number. */	
	double Length_Reynolds;	/*!< \brief Reynolds length (dimensional). */
	double AoA,			/*!< \brief Angle of attack (just external flow). */
	AoS;				/*!< \brief Angle of sideSlip (just external flow). */
	double ChargeCoeff;		/*!< \brief Charge coefficient (just for electric problems). */
	double *U_FreeStreamND;			/*!< \brief Reference variables at the infinity, free stream values. */ 
	unsigned short Cauchy_Func_Flow,	/*!< \brief Function where to apply the convergence criteria in the flow problem. */
	Cauchy_Func_AdjFlow,				/*!< \brief Function where to apply the convergence criteria in the adjoint problem. */
	Cauchy_Func_LinFlow,				/*!< \brief Function where to apply the convergence criteria in the linearized problem. */
	Cauchy_Elems;						/*!< \brief Number of elements to evaluate. */
	unsigned long StartConv_Iter;	/*!< \brief Start convergence criteria at iteration. */
	double Cauchy_Eps,	/*!< \brief Epsilon used for the convergence. */
	Cauchy_Eps_OneShot,	/*!< \brief Epsilon used for the one shot method convergence. */
	Cauchy_Eps_FullMG;	/*!< \brief Epsilon used for the full multigrid method convergence. */
	unsigned long Wrt_Sol_Freq,	/*!< \brief Writing solution frequency. */
	Wrt_Con_Freq;				/*!< \brief Writing convergence history frequency. */
	bool Wrt_Unsteady;  /*!< \brief Write unsteady data adding header and prefix. */
	bool Restart,	/*!< \brief Restart solution (for direct, adjoint, and linearized problems). */
	Restart_Flow;	/*!< \brief Restart flow solution for adjoint and linearized problems. */
	bool Block_Diagonal_Jacobian; /*!< \brief Block Diagonal Jacobian of multi species flow solution. */
	unsigned short nMarker_Monitoring,	/*!< \brief Number of markers to monitor. */
	nMarker_Plotting,					/*!< \brief Number of markers to plot. */
	nMarker_Moving;						/*!< \brief Number of markers to move. */
	string *Marker_Monitoring,			/*!< \brief Markers to monitor. */
	*Marker_Plotting,					/*!< \brief Markers to plot. */
	*Marker_Moving;						/*!< \brief Markers to move. */
	unsigned short  *Marker_All_Monitoring,				/*!< \brief Global index for monitoring using the grid information. */
	*Marker_All_Plotting,				/*!< \brief Global index for plotting using the grid information. */
	*Marker_All_Moving,					/*!< \brief Global index for moving using the grid information. */
	*Marker_Config_Monitoring,			/*!< \brief Global index for monitoring using the config information. */
	*Marker_Config_Plotting,			/*!< \brief Global index for plotting using the config information. */
	*Marker_Config_Moving,				/*!< \brief Global index for moving using the config information. */
	*Marker_Config_PerBound;			/*!< \brief Global index for periodic boundaries using the config information. */
	string *PlaneTag;			/*!< \brief Global index for the plane adaptation (upper, lower). */
	unsigned short nDomain;			/*!< \brief Number of domains in the MPI parallelization. */
	double DualVol_Power;			/*!< \brief Power for the dual volume in the grid adaptation sensor. */
	unsigned short Analytical_Surface;	/*!< \brief Information about the analytical definition of the surface for grid adaptation. */
	unsigned short Mesh_FileFormat;	/*!< \brief Mesh input format. */
	unsigned short Output_FileFormat;	/*!< \brief Format of the output files. */
	double RefAreaCoeff,		/*!< \brief Reference area for coefficient computation. */
	RefElemLength,				/*!< \brief Reference element length for computing the slope limiting epsilon. */
	RefLengthMoment,			/*!< \brief Reference length for moment computation. */
	*RefOriginMoment,			/*!< \brief Origin for moment computation. */
	*CFLRamp,			/*!< \brief Information about the CFL ramp. */
	DomainVolume;		/*!< \brief Volume of the computational grid. */	
	string Mesh_FileName,			/*!< \brief Mesh input file. */
	Mesh_Out_FileName,				/*!< \brief Mesh output file. */
	Solution_FlowFileName,			/*!< \brief Flow solution input file. */
	Solution_LinFileName,			/*!< \brief Linearized flow solution input file. */
	Solution_AdjFileName,			/*!< \brief Adjoint solution input file for drag functional. */
	Flow_FileName,					/*!< \brief Flow variables output file. */
	Residual_FileName,				/*!< \brief Residual variables output file. */
	Conv_FileName,					/*!< \brief Convergence history output file. */
	ReStart_FlowFileName,			/*!< \brief Restart file for flow variables. */
	ReStart_LinFileName,			/*!< \brief Restart file for linearized flow variables. */
	ReStart_AdjFileName,			/*!< \brief Restart file for adjoint variables, drag functional. */
	Adj_FileName,					/*!< \brief Output file with the adjoint variables. */
	Lin_FileName,					/*!< \brief Output file with the linearized variables. */
	ObjFunc_Grad_FileName,			/*!< \brief Gradient of the objective function. */
	SurfFlowCoeff_FileName,			/*!< \brief Output file with the flow variables on the surface. */
	SurfFlowCSV_FileName,			/*!< \brief Output file with the flow variables on the surface. */
	SurfAdjCoeff_FileName,			/*!< \brief Output file with the adjoint variables on the surface. */
	SurfAdjCSV_FileName,			/*!< \brief Output file with the adjoint variables on the surface. */
	SurfLinCoeff_FileName,			/*!< \brief Output file with the linearized variables on the surface. */
	New_SU2_FileName;        		/*!< \brief Output SU2 mesh file converted from CGNS format. */
	bool CGNS_To_SU2;      		 	/*!< \brief Flag to specify whether a CGNS mesh is converted to SU2 format. */
	unsigned short nSpecies, 		/*!< \brief No of species present in plasma */
	nReactions,									/*!< \brief Number of reactions in chemical model. */
	nFluids;						/*!< \brief No of fluids modeled in plasma */						/*!< \brief No of fluids modeled in plasma */
	bool Write_Mean_Solution;		/*!< \brief Write a mean solution file */
	double *ArrheniusCoefficient,					/*!< \brief Arrhenius reaction coefficient */
	*ArrheniusEta,								/*!< \brief Arrhenius reaction temperature exponent */
	*ArrheniusTheta;							/*!< \brief Arrhenius reaction characteristic temperature */
	unsigned short nMass,                 /*!< \brief No of particle masses */
	nRef_Temperature,   			/*!< \brief No of particle Reference Temperature */
	nRef_Viscosity,   				/*!< \brief No of particle Reference Viscosity */
	nMagnet;							/*!< \brief This value must always be 3 for a magnet */
	int *Charge_Number;			/*!< \brief Charge number of all species present (+1/0/-1) */
	double *Particle_Mass,					/*!< \brief Mass of all particles present in the plasma */
	*Molar_Mass,								/*!< \brief Molar mass of species in the plasma [kg/kmol] */
	*Molecular_Diameter,			/*!< \brief Molecular diameter of species [m] */
	*Gas_Composition,					/*!< \brief Initial mass fractions of flow [dimensionless] */
	*Enthalpy_Formation,			/*!< \brief Enthalpy of formation */
	*Species_Ref_Temperature,	/*!< \brief Reference Temperature for viscosity of all particles present in the plasma */
	*Species_Ref_Viscosity,		/*!< \brief Reference viscosity  of all particles present in the plasma */
	*MagneticDipole;			/*!< \brief Magnetic dipole present in the plasma */
	double Gamma,			/*!< \brief Ratio of specific heats of the gas. */
	GammaDiatomic,			/*!< \brief Ratio of specific heats of the diatomic gas. */ 
	GammaMonatomic,			/*!< \brief Ratio of specific heats of the monatomic gas. */ 
	Bulk_Modulus,			/*!< \brief Value of the bulk modulus for incompressible flows. */ 
	ArtComp_Factor,			/*!< \brief Value of the artificial compresibility factor for incompressible flows. */
	Gas_Constant,     /*!< \brief Specific gas constant. */
	Gas_Constant_Ref, /*!< \brief Reference specific gas constant. */
	FreeSurface_Zero,	/*!< \brief Coordinate of the level set zero. */
	FreeSurface_Depth,	/*!< \brief Coordinate of the level set zero. */
	*Velocity_FreeStream,     /*!< \brief Total velocity of the fluid.  */
	Density_FreeStream,     /*!< \brief Total density of the fluid.  */
	Viscosity_FreeStream,     /*!< \brief Total density of the fluid.  */
	Pressure_FreeStream,     /*!< \brief Total pressure of the fluid.  */
	Temperature_FreeStream,  /*!< \brief Total temperature of the fluid.  */
	Prandtl_Lam,      /*!< \brief Laminar Prandtl number for the gas.  */
	Prandtl_Turb,     /*!< \brief Turbulent Prandtl number for the gas.  */
	Length_Ref,       /*!< \brief Reference length for non-dimensionalization. */
	Conversion_Factor,       /*!< \brief Conversion factor from grid units to meters. */
	Pressure_Ref,     /*!< \brief Reference pressure for non-dimensionalization. */
	Temperature_Ref,  /*!< \brief Reference temperature for non-dimensionalization. */
	Density_Ref,      /*!< \brief Reference density for non-dimensionalization. */
	Velocity_Ref,     /*!< \brief Reference velocity for non-dimensionalization. */
	Time_Ref,         /*!< \brief Reference time for non-dimensionalization. */
	Viscosity_Ref,    /*!< \brief Reference viscosity for non-dimensionalization. */
	Wall_Temperature,    /*!< \brief Temperature at an isotropic wall in Kelvin. */
	Omega_Ref,        /*!< \brief Reference angular velocity for non-dimensionalization. */
	Force_Ref,        /*!< \brief Reference body force for non-dimensionalization. */
	Pressure_FreeStreamND,     /*!< \brief Farfield pressure value (external flow). */
	Temperature_FreeStreamND,  /*!< \brief Farfield temperature value (external flow). */
	Density_FreeStreamND,      /*!< \brief Farfield density value (external flow). */
	*Velocity_FreeStreamND,    /*!< \brief Farfield velocity values (external flow). */
	*Omega_FreeStreamND,       /*!< \brief Farfield angular velocity values (external flow). */
	Energy_FreeStreamND,       /*!< \brief Farfield energy value (external flow). */
	Viscosity_FreeStreamND;    /*!< \brief Farfield viscosity value (external flow). */
	bool Write_Converted_Mesh; /*!< \brief Flag to specify whether a new mesh should be written in the converted units. */
	double Scale_GradOF;			 /*!< \brief Scale of the gradient of the objective function. */

	map<string, CAnyOptionRef*> param; /*!< \brief associates option names (strings) with options */

public:

	/*! 
	 * \brief Constructor of the class which reads the input file.
	 */
	CConfig(char case_filename[200], unsigned short val_software);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CConfig(void);

	/*!
	 * \brief add a scalar option to the param map and set its default value
	 * \param[in] name - name of the option as it appears in the .cfg file
	 * \param[in,out] option - the option to associate with name
	 * \param[in] default_value - option is set to default_value
	 * \tparam T - an arbitary type (int, double, enum, etc)
	 */
	template <class T, class T_default>
	void AddScalarOption(const string & name, T & option, const T_default & default_value) {
		//cout << "Adding Scalar option " << name << endl;
		option = static_cast<T>(default_value);
		CAnyOptionRef* option_ref = new COptionRef<T>(option);
		param.insert( pair<string, CAnyOptionRef*>(string(name), option_ref) );
	}

	/*!
	 * \brief add an enum-based option to the param map and set its default value
	 * \param[in] name - name of the option as it appears in the .cfg file
	 * \param[in,out] option - the option to associate with name
	 * \param[in] Tmap - a map from strings to type Tenum
	 * \tparam T - the option type (usually unsigned short)
	 * \tparam Tenum - an enumeration assocaited with T
	 * \param[in] default_value - option is set to default_value
	 */
	template <class T, class Tenum>
	void AddEnumOption(const string & name, T & option, const map<string, Tenum> & Tmap,
			const string & default_value) {
		//cout << "Adding Enum option " << name << endl;
		typename map<string,Tenum>::const_iterator it;
		it = Tmap.find(default_value);
		if (it == Tmap.end()) {
			cerr << "Error in CConfig::AddEnumOption(string, T&, const map<string, Tenum> &, const string): "
					<< "cannot find " << default_value << " in given map."
					<< endl;
			throw(-1);
		}
		option = it->second;
		CAnyOptionRef* option_ref = new CEnumOptionRef<T,Tenum>(option, Tmap);
		param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
	}

	/*!
	 * \brief add an enum-based array option to the param map
	 * \param[in] name - name of the option as it appears in the .cfg file
	 * \param[in,out] size - a reference to the size of option
	 * \param[in,out] option - the option to associate with name
	 * \param[in] Tmap - a map from strings to type Tenum
	 * \param[in] update - set to true if the option has already been initialized
	 * \tparam T - the option type (usually unsigned short)
	 * \tparam Tenum - an enumeration assocaited with T
	 */
	template <class T, class Tenum>
	void AddEnumListOption(const string & name, unsigned short & size, T* & option,
			const map<string, Tenum> & Tmap,
			const bool & update = false) {
		//cout << "Adding Enum-List option " << name << endl;
		size = 0;
		if (update && option != NULL) delete [] option;
		CAnyOptionRef* option_ref = new CEnumOptionRef<T,Tenum>(size, option, Tmap);
		param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
	}

	/*!
	 * \brief add an array option to the param map and set its default value
	 * \param[in] name - name of the option as it appears in the .cfg file
	 * \param[in] size - length of option and default_value arrays
	 * \param[in,out] option - the option to associate with name
	 * \param[in] default_value - option is set to default_value
	 * \param[in] update - set to true if the option has already been initialized
	 * \tparam T - an arbitary type (int, double, enum, etc)
	 *
	 * The option array is allocated in this function.  If memory for
	 * the option array is already allocated, this memory is first released.
	 */
	template <class T>
	void AddArrayOption(const string & name, const int & size, T* & option,
			const T* default_value, const bool & update = false) {
		//cout << "Adding Array option " << name << endl;
		if (update && option != NULL) delete [] option;
		option = new T[size];
		for (int i = 0; i < size; i++)
			option[i] = default_value[i];
		CAnyOptionRef* option_ref = new COptionRef<T>(option, size);
		param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
	}

	/*!
	 * \brief add a list option to the param map
	 * \param[in] name - name of the option as it appears in the .cfg file
	 * \param[in] size - length of option (will be set by data in .cfg file)
	 * \param[in] option - the option to associate with name
	 * \param[in] update - set to true if the option has already been initialized
	 * \tparam T - an arbitary type (int, double, enum, etc)
	 *
	 * This is similar to AddArrayOption, but is useful for options of
	 * variable size.  Also, this routine does not allow the default
	 * value to be set.
	 */
	template <class T>
	void AddListOption(const string & name, unsigned short & size, T* & option,
			const bool & update = false) {
		//cout << "Adding List option " << name << endl;
		if (update && option != NULL) delete [] option;
		CAnyOptionRef* option_ref = new CListOptionRef<T>(size, option);
		param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
	}

	/*!
	 * \brief add a special option to the param map and set its default value
	 * \param[in] name - string name of the option as it appears in the .cfg file
	 * \param[in,out] option - the option to associate with name
	 * \param[in] set_value - function to parse string vectors and set option
	 * \param[in] default_value - option is set to default_value
	 * \tparam T - an arbitrary type (int, double, bool, enum, etc)
	 */
	template <class T>
	void AddSpecialOption(const string & name, T & option,
			void (*set_value)(T*, const vector<string>&),
			const T & default_value) {
		//cout << "Adding Special option " << name << endl;
		option = default_value;
		CAnyOptionRef* option_ref = new COptionRef<T>(option, set_value);
		param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
	}

	/*!
	 * \brief add a marker-type option to the param map
	 * \param[in] name - name of the marker option as it appears in the .cfg file
	 * \param[in,out] num_marker - number of boundary markers
	 * \param[in,out] marker - an array of boundary marker names
	 */
	void AddMarkerOption(const string & name, unsigned short & num_marker, string* & marker);

	/*!
	 * \brief add a convection-discretization type option to the param map
	 * \param[in] name - name of the convection option as it appears in the .cfg file
	 * \param[in] space - the spatial discretization type associataed with name
	 * \param[in] centred - the centred spatial discretization type of name
	 * \param[in] upwind - the upwind spatial discretization type of name
	 */
	void AddConvectOption(const string & name, unsigned short & space, unsigned short & centred,
			unsigned short & upwind);

	/*!
	 * \brief adds the math problem option to the param map
	 * \param[in] name - name of the math problem as it appears in the .cfg file
	 * \param[in] ContAdj - is the continuous adjoint solved?
	 * \param[in] ContAdj_default - the default value for ContAdj
	 * \param[in] OneShot - is a one-shot problem solved?
	 * \param[in] OneShot_default - the default value for OneShot
	 * \param[in] Linearized - is a linearized problem solved?
	 * \param[in] Linearized_default - the default value for Linearized
	 * \param[in] Restart_Flow - is the flow restarted for adjoint and linearized problems?
	 * \param[in] Restart_Flow_default - the default value for Restart_Flow
	 */
	void AddMathProblem(const string & name, bool & ContAdj, const bool & ContAdj_default,
			bool & OneShot, const bool & OneShot_default,
			bool & Linearized, const bool & Linearized_default,
			bool & Restart_Flow, const bool & Restart_Flow_default);

	/*!
	 * \brief adds the design variable parameters option to the param map
	 * \param[in] name - name of the design-variable parameters option in the config file
	 * \param[in] nDV - the number of design variables
	 * \param[in] ParamDV - the parameter values of each design variable
	 * \param[in] Design_Variable - the type of each design variable
	 */
	void AddDVParamOption(const string & name, unsigned short & nDV, double** & ParamDV,
			unsigned short* & Design_Variable);

	/*!
	 * \brief adds a periodic marker option to the param map
	 * \param[in] name - name of the periodic marker option in the config file
	 * \param[in] nMarker_PerBound - the number of periodic marker boundaries
	 * \param[in] Marker_PerBound - string names of periodic boundaries
	 * \param[in] Marker_PerDonor - names of boundaries that supply data to periodic boundaries
	 * \param[in] RotCenter - rotational center for each periodic boundary
	 * \param[in] RotAngles - rotation angles for each periodic boundary
	 * \param[in] Translation - translation vector for each periodic boundary.
	 */
	void AddMarkerPeriodic(const string & name, unsigned short & nMarker_PerBound,
			string* & Marker_PerBound, string* & Marker_PerDonor,
			double** & RotCenter, double** & RotAngles, double** & Translation);

	/*!
	 * \brief adds an inlet marker option to the param map
	 * \param[in] name - name of the inlet marker option in the config file
	 * \param[in] nMarker_Inlet - the number of inlet marker boundaries
	 * \param[in] Marker_Inlet - string names of inlet boundaries
	 * \param[in] Ttotal - specified total temperatures for inlet boundaries
	 * \param[in] Ptotal - specified total pressures for inlet boundaries
	 * \param[in] FlowDir - specified flow direction vector (unit vector) for inlet boundaries
	 */
	void AddMarkerInlet(const string & name, unsigned short & nMarker_Inlet,
			string* & Marker_Inlet, double* & Ttotal, double* & Ptotal,
			double** & FlowDir);

	/*!
	 * \brief adds an outlet marker option to the param map
	 * \param[in] name - name of the outlet marker option in the config file
	 * \param[in] nMarker_Outlet - the number of outlet marker boundaries
	 * \param[in] Marker_Outlet - string names of outlet boundaries
	 * \param[in] Pressure - Specified back pressures (static) for outlet boundaries
	 */
	void AddMarkerOutlet(const string & name, unsigned short & nMarker_Outlet,
			string* & Marker_Outlet, double* & Pressure);

	/*!
	 * \brief used to set Boolean values based on strings "YES" and "NO"
	 * \param[in] ref - a pointer to the boolean value being assigned
	 * \param[in] value - value[0] is "YES" or "NO" and determines the value of ref
	 */
	static void SetBoolOption(bool* ref, const vector<string> & value);

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
	bool GetPython_Option(string & option_name);

	/*! 
	 * \brief Get reference origin for moment computation.
	 * \return Reference origin (in cartesians coordinates) for moment computation.
	 */
	double *GetRefOriginMoment(void);

	/*! 
	 * \brief Get maximum number of children in the agglomeration process.
	 * \return Maximum number of children.
	 */
	unsigned short GetMaxChildren(void);

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
	 * \brief Get the maximum dimension of the agglomerated element compared with the whole domain.
	 * \return Maximum dimension of the agglomerated element.
	 */
	double GetMaxDimension(void);

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
	 * \brief Get info about a free surface problem.
	 * \return Free surface problem.
	 */
	bool GetFreeSurface(void);

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
	 * \brief Get the inlet position of the free surface for a free surface problem.
	 * \return Inlet position of the interfase for a free surface problem.
	 */
	double GetFreeSurface_Inlet(void);

	/*!
	 * \brief Get the outlet position of the free surface for a free surface problem.
	 * \return Outlet position of the interfase for a free surface problem.
	 */
	double GetFreeSurface_Outlet(void);


	/*! 
	 * \brief Get information about the the restart cycling in the MG.
	 * \return <code>TRUE</code> if the MG cycle solution is restarted at each cycle; otherwise <code>FALSE</code>.
	 */
	bool GetReStartMGCycle(void);

	/*! 
	 * \brief Get information about the kind of problem that we are solving.
	 * \return <code>TRUE</code> if a shock tube problem is defined; otherwise <code>FALSE</code>.
	 */
	bool GetShockTube(void);

	/*! 
	 * \brief Creates a paraview file to visualize the partition made by the DDM software.
	 * \return <code>TRUE</code> if the partition is going to be plotted; otherwise <code>FALSE</code>.
	 */
	bool GetVisualize_Partition(void);

	/*! 
	 * \brief Creates a paraview file to visualize the deformation made by the MDC software.
	 * \return <code>TRUE</code> if the deformation is going to be plotted; otherwise <code>FALSE</code>.
	 */
	bool GetVisualize_Deformation(void);

	/*! 
	 * \brief Get the value of the Mach number (velocity divided by speed of sound).
	 * \return Value of the Mach number.
	 */
	double GetMach_FreeStreamND(void);

	/*! 
	 * \brief Get the value of the Gamma of fluid (ratio of specific heats).
	 * \return Value of the constant: Gamma
	 */
	double GetGamma(void);

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
	double GetPressure_FreeStreamND(void);

	/*!
	 * \brief Get the value of the non-dimensionalized freestream pressure.
	 * \return Non-dimensionalized freestream pressure.
	 */
	double GetPressure_FreeStream(void);

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
	 * \brief Get the vector of the non-dimensionalized freestream angular velocity (rotating frame).
	 * \return Non-dimensionalized freestream angular velocity vector (rotating frame).
	 */
	double* GetOmega_FreeStreamND(void);

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
	 * \brief Get the value of the Reynolds length.
	 * \return Reynolds length.
	 */
	double GetLength_Reynolds(void);

	/*!
	 * \brief Get the conversion factor for converting the grid to meters.
	 * \return Conversion factor for converting the grid to meters.
	 */
	double GetConversion_Factor(void);

	/*! 
	 * \brief Get the value of the Primitive variables gradient threshold (just for the adjoint problem).
	 * \return Value of the Primitive variables gradient threshold.
	 */
	double GetPrimGrad_Threshold(void);

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
	 * \param[in] val_kind_centred - If centred scheme, kind of centred scheme (JST, etc.).
	 * \param[in] val_kind_upwind - If upwind scheme, kind of upwind scheme (Roe, etc.).
	 * \param[in] val_kind_slopelimit - If upwind scheme, kind of slope limit.
	 */		
	void SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme, unsigned short val_kind_centred, 
			unsigned short val_kind_upwind, unsigned short val_kind_slopelimit);

	/*! 
	 * \brief Set the parameters of the viscous numerical scheme.
	 * \note The parameters will change because we are solving different kind of equations.
	 * \param[in] val_kind_viscnumscheme - Kind of viscous scheme.
	 */
	void SetKind_ViscNumScheme(unsigned short val_kind_viscnumscheme);

	/*! 
	 * \brief Set the parameters of the source term.
	 * \note The parameters will change because we are solving different kind of equations.
	 * \param[in] val_kind_sournumscheme - Kind of source term.
	 */
	void SetKind_SourNumScheme(unsigned short val_kind_sournumscheme);

	/*! 
	 * \brief Get the number of smoothing iterations that the implicit residual smoothing 
	 *            technique is performing.
	 * \note If <i>nSmooth</i> is 0 then, it doesn't make an implicit smoothing.
	 * \return Number of smoothing iteration.
	 */
	unsigned short GetnSmooth(void);

	/*! 
	 * \brief Get the value of the relaxation factor of the implicit residual smoothing technique.
	 * \return Value of the relaxation factor.
	 */
	double GetSmoothCoeff(void);

	/*! 
	 * \brief Get the value of limiter coefficient.
	 * \return Value of the limiter coefficient.
	 */
	double GetLimiterCoeff(void);

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
	 * \brief Get the angle of attack of the body. This is the angle between a reference line on a lifting body 
	 *        (often the chord line of an airfoil) and the vector representing the relative motion between the 
	 *        lifting body and the fluid through which it is moving.
	 * \return Value of the angle of attack.
	 */		
	double GetAoA(void);

	/*! 
	 * \brief Get the angle of sideslip of the body. It relates to the rotation of the aircraft centerline from 
	 *        the relative wind.
	 * \return Value of the angle of sideslip.
	 */		
	double GetAoS(void);

	/*! 
	 * \brief Get the charge coefficient that is used in the electrical potential simulation.
	 * \return Value of the charge coefficient.
	 */		
	double GetChargeCoeff(void);	

	/*! 
	 * \brief Get the number of multigrid levels.
	 * \return Number of multigrid levels (without including the original grid).
	 */
	unsigned short GetMGLevels(void);

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
	 * \brief Get the Courant Friedrich Levi number for each grid.
	 * \param[in] val_mesh - Index of the mesh were the CFL is applied.
	 * \return CFL number for each grid.
	 */		
	double GetCFL(unsigned short val_mesh);

	/*! 
	 * \brief Get the Courant Friedrich Levi number for unsteady simulations.
	 * \return CFL number for unsteady simulations.
	 */			
	double GetUnst_CFL(void);

	/*! 
	 * \brief Get a parameter of the particular design variable.
	 * \param[in] val_dv - Number of the design variable that we want to read.
	 * \param[in] val_param - Index of the parameter that we want to read.
	 * \return Design variable parameter.
	 */		
	double GetParamDV(unsigned short val_dv, unsigned short val_param);

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
	 * \brief Set the number of external iterations.
	 * \note This is important in no time depending methods, where only 
	 *       one external iteration is needed.
	 * \param[in] val_niter - Set the number of external iterations.
	 */
	void SetnExtIter(unsigned long val_niter);

	/*! 
	 * \brief Get the frequency for writing the solution file.
	 * \return It writes the solution file with this frequency.
	 */		
	unsigned long GetWrt_Sol_Freq(void);

	/*! 
	 * \brief Get the frequency for writing the convergence file.
	 * \return It writes the convergence file with this frequency.
	 */		
	unsigned long GetWrt_Con_Freq(void);

	/*!
	 * \brief Get information about writing headers and prefix on the unsteady data.
	 * \return 	<code>TRUE</code> means that headers and prefix on the unsteady data will be written.
	 */
	bool GetWrt_Unsteady(void);

	/*! 
	 * \brief Get the alpha (convective) coefficients for the Runge-Kutta integration scheme.
	 * \param[in] val_step - Index of the step.
	 * \return Alpha coefficient for the Runge-Kutta integration scheme.
	 */		
	double Get_Alpha_RKStep(unsigned short val_step);

	/*! 
	 * \brief Get the beta (viscous) coefficients for the Runge-Kutta integration scheme.
	 * \param[in] val_step - Index of the step.
	 * \return Beta coefficient for the Runge-Kutta integration scheme.
	 */		
	double Get_Beta_RKStep(unsigned short val_step);

	/*! 
	 * \brief Get the index of the surface defined in the geometry file.
	 * \param[in] val_marker - Value of the marker in which we are interested.
	 * \return Value of the index that is in the geometry file for the surface that 
	 *         has the marker <i>val_marker</i>.
	 */		
	string GetMarker_All_Tag(unsigned short val_marker);
	
	/*! 
	 * \brief Get the tag if the iMarker defined in the geometry file.
	 * \param[in] val_tag - Value of the tag in which we are interested.
	 * \return Value of the marker <i>val_marker</i> that is in the geometry file 
	 *         for the surface that has the tag.
	 */		
	unsigned short GetTag_Marker_All(string val_tag);

	/*! 
	 * \brief Get the kind of boundary for each marker.
	 * \param[in] val_marker - Index of the marker in which we are interested.
	 * \return Kind of boundary for the marker <i>val_marker</i>.
	 */		
	unsigned short GetMarker_All_Boundary(unsigned short val_marker);

	/*! 
	 * \brief Set the value of the boundary <i>val_boundary</i> (read from the config file) 
	 *        for the marker <i>val_marker</i>.
	 * \param[in] val_marker - Index of the marker in which we are interested.
	 * \param[in] val_boundary - Kind of boundary read from config file.
	 */		
	void SetMarker_All_Boundary(unsigned short val_marker, unsigned short val_boundary);

	/*! 
	 * \brief Set the value of the index <i>val_index</i> (read from the geometry file) for 
	 *        the marker <i>val_marker</i>.
	 * \param[in] val_marker - Index of the marker in which we are interested.
	 * \param[in] val_index - Index of the surface read from geometry file.
	 */	
	void SetMarker_All_Tag(unsigned short val_marker, string val_index);

	/*! 
	 * \brief Set if a marker <i>val_marker</i> is going to be monitored <i>val_monitoring</i> 
	 *        (read from the config file).
	 * \note This is important for non dimensional coefficient computation.
	 * \param[in] val_marker - Index of the marker in which we are interested.
	 * \param[in] val_monitoring - 0 or 1 depending if the the marker is going to be monitored.
	 */	
	void SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring);

	/*! 
	 * \brief Set if a marker <i>val_marker</i> is going to be plot <i>val_plotting</i> 
	 *        (read from the config file).
	 * \param[in] val_marker - Index of the marker in which we are interested.
	 * \param[in] val_plotting - 0 or 1 depending if the the marker is going to be plot.
	 */	
	void SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting);

	/*! 
	 * \brief Set if a marker <i>val_marker</i> is going to be move <i>val_moving</i> 
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
	 * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
	 * \return 0 or 1 depending if the marker is going to be monitorized.
	 */		
	unsigned short GetMarker_All_Monitoring(unsigned short val_marker);

	/*! 
	 * \brief Get the plotting information for a marker <i>val_marker</i>.
	 * \param[in] val_marker - 0 or 1 depending if the the marker is going to be moved.
	 * \return 0 or 1 depending if the marker is going to be plotted.
	 */		
	unsigned short GetMarker_All_Plotting(unsigned short val_marker);

	/*! 
	 * \brief Get the moving information for a marker <i>val_marker</i>.
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
	 * \return Governing equation that we are solving.
	 */		
	unsigned short GetKind_Solver(void);

	/*! 
	 * \brief Gas model that we are using.
	 * \return Gas model that we are using.
	 */		
	unsigned short GetKind_GasModel(void);

	/*! 
	 * \brief Set the kind of solver that is going to be used.
	 * \note Only in case of the master node in the MPI parallelization it 
	 *       is necessary to change the original solver.
	 * \param[in] val_kind_solver - Kind of solver.
	 */		
	void SetKind_Solver(unsigned short val_kind_solver);

	/*! 
	 * \brief Get the kind of method for computation of spatial gradients.
	 * \return Numerical method for computation of spatial gradients.
	 */		
	unsigned short GetKind_Gradient_Method(void);

	/*! 
	 * \brief Get the kind of method for deforming the numerical grid.
	 * \return Numerical method for deforming the numerical grid.
	 */
	unsigned short GetKind_GridDef_Method(void);

	/*!
	 * \brief Get the kind of solver for deforming the numerical grid (linear system solving).
	 * \return Numerical solver for deforming the numerical grid (solving the linear system).
	 */
	unsigned short GetKind_GridDef_Solver(void);

	/*!
	 * \brief Get the kind of solver for the implicit solver.
	 * \return Numerical solver for implicit formulation (solving the linear system).
	 */
	unsigned short GetKind_Linear_Solver(void);

	/*!
	 * \brief Get min error of the linear solver for the implicit formulation.
	 * \return Min error of the linear solver for the implicit formulation.
	 */
	double GetLinear_Solver_Error(void);

	/*! 
	 * \brief Get max number of iterations of the linear solver for the implicit formulation.
	 * \return Max number of iterations of the linear solver for the implicit formulation.
	 */
	unsigned short GetLinear_Solver_Iter(void);

	/*!
	 * \brief Get the error of the solver for deforming the numerical grid (linear system solving).
	 * \return Numerical error of the solver for deforming the numerical grid (solving the linear system).
	 */
	double GetGridDef_Error(void);

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
	 * \brief Get the kind of viscous numerical scheme.
	 * \note This is the information that the code will use, the method will 
	 *       change in runtime depending of the specific equation (direct, adjoint, 
	 *       linearized) that is being solved.	 
	 * \return Kind of the viscous scheme.
	 */		
	unsigned short GetKind_ViscNumScheme(void);

	/*! 
	 * \brief Get the kind scheme for the source term integration.
	 * \note This is the information that the code will use, the method will 
	 *       change in runtime depending of the specific equation (direct, adjoint, 
	 *       linearized) that is being solved.	 
	 * \return Scheme for the source term integration.
	 */		
	unsigned short GetKind_SourNumScheme(void);

	/*! 
	 * \brief Get kind of center scheme for the convective terms.
	 * \note This is the information that the code will use, the method will 
	 *       change in runtime depending of the specific equation (direct, adjoint, 
	 *       linearized) that is being solved.	 
	 * \return Kind of center scheme for the convective terms.
	 */		
	unsigned short GetKind_Centred(void);

	/*! 
	 * \brief Get kind of upwind scheme for the convective terms.
	 * \note This is the information that the code will use, the method will 
	 *       change in runtime depending of the specific equation (direct, adjoint, 
	 *       linearized) that is being solved.	 
	 * \return Kind of upwind scheme for the convective terms.
	 */		
	unsigned short GetKind_Upwind(void);

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
	 * \return Kind of integration scheme for the plasma equations.
	 */
	unsigned short GetKind_TimeIntScheme_Plasma(void);
	
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
	 *        for the template equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of integration scheme for the plasma equations.
	 */
	unsigned short GetKind_TimeIntScheme_Template(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the flow 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */		
	unsigned short GetKind_ConvNumScheme_Flow(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the plasma 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */		
	unsigned short GetKind_ConvNumScheme_Plasma(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the template 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the flow equations.
	 */		
	unsigned short GetKind_ConvNumScheme_Template(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the level set 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the level set equation.
	 */		
	unsigned short GetKind_ConvNumScheme_LevelSet(void);

	/*!
	 * \brief Get the kind of convective numerical scheme for the adjoint level set
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the level set equation.
	 */
	unsigned short GetKind_ConvNumScheme_AdjLevelSet(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the level set 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the level set equation.
	 */		
	unsigned short GetKind_ConvNumScheme_Combustion(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the flow
	 *        equations (Galerkin, Average of gradients, Average of gradients
	 *        with correction).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the flow equations.
	 */		
	unsigned short GetKind_ViscNumScheme_Flow(void);

	/*!
	 * \brief Get the kind of viscous numerical scheme for the level set
	 (        equation.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the levelset equations.
	 */
	unsigned short GetKind_SourNumScheme_LevelSet(void);

	/*!
	 * \brief Get the kind of viscous numerical scheme for the wave
	 (        equation.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the levelset equations.
	 */		
	unsigned short GetKind_SourNumScheme_Wave(void);
	
	/*! 
	 * \brief Get the kind of viscous numerical scheme for the adjoint level set
	 (        equation.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the levelset equations.
	 */
	unsigned short GetKind_ViscNumScheme_AdjLevelSet(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the plasma
	 *        equations (Galerkin, Average of gradients)
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the flow equations.
	 */		
	unsigned short GetKind_ViscNumScheme_Plasma(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the plasma
	 *        equations (Galerkin, Average of gradients)
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the flow equations.
	 */		
	unsigned short GetKind_ViscNumScheme_Template(void);

	/*! 
	 * \brief Get the kind of source term for the flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the flow equations.
	 */			
	unsigned short GetKind_SourNumScheme_Flow(void);

	/*! 
	 * \brief Get the kind of source term for the plasma equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the flow equations.
	 */			
	unsigned short GetKind_SourNumScheme_Plasma(void);

	/*! 
	 * \brief Get the kind of source term for the combustion equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the flow equations.
	 */			
	unsigned short GetKind_SourNumScheme_Combustion(void);

	/*! 
	 * \brief Get the kind of source term for the plasma equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the flow equations.
	 */			
	unsigned short GetKind_SourNumScheme_Template(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Centred_Flow(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the level set equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the level set equations.
	 */
	unsigned short GetKind_Centred_LevelSet(void);

	/*!
	 * \brief Get the kind of center convective numerical scheme for the adjoint level set equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the level set equations.
	 */
	unsigned short GetKind_Centred_AdjLevelSet(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the level set equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the level set equations.
	 */
	unsigned short GetKind_Centred_Combustion(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the plasma equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Centred_Plasma(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the plasma equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Centred_Template(void);

	/*! 
	 * \brief Get the kind of upwind convective numerical scheme for the flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_Flow(void);

	/*! 
	 * \brief Get the kind of upwind convective numerical scheme for the level set equation.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_LevelSet(void);

	/*!
	 * \brief Get the kind of upwind convective numerical scheme for the adjoint level set equation.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_AdjLevelSet(void);

	/*! 
	 * \brief Get the kind of upwind convective numerical scheme for the level set equation.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_Combustion(void);

	/*! 
	 * \brief Get the kind of upwind convective numerical scheme for the plasma equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the flow equations.
	 */
	unsigned short GetKind_Upwind_Plasma(void);

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
	 * \return Method for limiting the spatial gradients solving the level set equation.
	 */		
	unsigned short GetKind_SlopeLimit_LevelSet(void);

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
	 * \brief Get the method for limiting the spatial gradients.
	 * \return Method for limiting the spatial gradients solving the Plasma equations.
	 */
	unsigned short GetKind_SlopeLimit_Plasma(void);

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
	 * \brief Get the kind of integration scheme (explicit or implicit) 
	 *        for the adjoint flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of integration scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_AdjFlow(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the adjoint flow 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_AdjFlow(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the adjoint flow
	 *        equations (Galerkin, Average of gradients, Average of gradients
	 *        with correction).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_ViscNumScheme_AdjFlow(void);
	
	/*! 
	 * \brief Get the kind of viscous numerical scheme for the wave
	 *        equations (Galerkin, Average of gradients, Average of gradients
	 *        with correction).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_ViscNumScheme_Wave(void);

	/*! 
	 * \brief Get the kind of source term for the adjoint flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the adjoint flow equations.
	 */
	unsigned short GetKind_SourNumScheme_AdjFlow(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the adjoint flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the adjoint flow equations.
	 */
	unsigned short GetKind_Centred_AdjFlow(void);

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
	 * \brief Get the kind of integration scheme (explicit or implicit) 
	 *        for the linearized flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of integration scheme for the linearized flow equations.
	 */
	unsigned short GetKind_TimeIntScheme_LinFlow(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the linearized flow 
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the linearized flow equations.
	 */
	unsigned short GetKind_ConvNumScheme_LinFlow(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the linearized flow
	 *        equations (Galerkin, Divergence theorem or Weiss correction).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the linearized flow equations.
	 */
	unsigned short GetKind_ViscNumScheme_LinFlow(void);

	/*! 
	 * \brief Get the kind of source term for the linearized flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the linearized flow equations.
	 */
	unsigned short GetKind_SourNumScheme_LinFlow(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the linearized flow equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the linearized flow equations.
	 */
	unsigned short GetKind_Centred_LinFlow(void);

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
	unsigned short GetKind_TimeIntScheme_LevelSet(void);

	/*!
	 * \brief Get the kind of integration scheme (implicit)
	 *        for the level set equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of integration scheme for the level set equations.
	 */
	unsigned short GetKind_TimeIntScheme_AdjLevelSet(void);

	/*! 
	 * \brief Get the kind of integration scheme (implicit) 
	 *        for the level set equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of integration scheme for the level set equations.
	 */
	unsigned short GetKind_TimeIntScheme_Combustion(void);

	/*! 
	 * \brief Get the kind of convective numerical scheme for the turbulence 
	 *        equations (upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the turbulence equations.
	 */
	unsigned short GetKind_ConvNumScheme_Turb(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the turbulence
	 *        equations (Galerkin, Average of gradients, Average of gradients
	 *        with correction).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the turbulence equations.
	 */	
	unsigned short GetKind_ViscNumScheme_Turb(void);

	/*! 
	 * \brief Get the kind of source term for the turbulence equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the turbulence equations.
	 */
	unsigned short GetKind_SourNumScheme_Turb(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the turbulence equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the turbulence equations.
	 */
	unsigned short GetKind_Centred_Turb(void);

	/*! 
	 * \brief Get the kind of upwind convective numerical scheme for the turbulence equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of upwind convective numerical scheme for the turbulence equations.
	 */	
	unsigned short GetKind_Upwind_Turb(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the electric potential
	 *        equation (Galerkin).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the electric potential equation.
	 */
	unsigned short GetKind_ViscNumScheme_Elec(void);

	/*! 
	 * \brief Get the kind of source term for the electric potential equation.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the electric potential equation.
	 */	
	unsigned short GetKind_SourNumScheme_Elec(void);

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
	 *        equations (centred or upwind).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of convective numerical scheme for the adjoint turbulence equations.
	 */
	unsigned short GetKind_ConvNumScheme_AdjTurb(void);

	/*! 
	 * \brief Get the kind of viscous numerical scheme for the adjoint turbulence
	 *        equations (Galerkin, Average of gradients, Average of gradients
	 *        with correction).
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of viscous numerical scheme for the adjoint turbulence equations.
	 */
	unsigned short GetKind_ViscNumScheme_AdjTurb(void);

	/*! 
	 * \brief Get the kind of source term for the adjoint turbulence equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of source term for the adjoint turbulence equations.
	 */
	unsigned short GetKind_SourNumScheme_AdjTurb(void);

	/*!
	 * \brief Get the kind of source term for the adjoint levelset equations.
	 * \note This value is obtained from the config file, and it is constant
	 *       during the computation.
	 * \return Kind of source term for the adjoint levelset equations.
	 */
	unsigned short GetKind_SourNumScheme_AdjLevelSet(void);

	/*! 
	 * \brief Get the kind of center convective numerical scheme for the adjoint turbulence equations.
	 * \note This value is obtained from the config file, and it is constant 
	 *       during the computation.
	 * \return Kind of center convective numerical scheme for the adjoint turbulence equations.
	 */
	unsigned short GetKind_Centred_AdjTurb(void);

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
	 * \brief Provides information about the the nodes that are going to be moved on a deformation 
	 *        volumetric grid deformation.
	 * \return <code>TRUE</code> means that only the points on the FFD box will be moved.
	 */
	bool GetMove_FFDBox(void);

	/*! 
	 * \brief Get the kind of objective function. There are several options: Drag coefficient, 
	 *        Lift coefficient, efficiency, etc.
	 * \note The objective function will determine the boundary condition of the adjoint problem.
	 * \return Kind of objective function.
	 */
	unsigned short GetKind_ObjFunc(void);

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
	 * \brief Provides the number of fluids to imagine the plasma is made up of. 
	 * \return: The number of fluids present in the plasma, read from input file
	 */
	unsigned short GetnFluids(void);

	/*! 
	 * \brief Provides a table of equilibrium constants for a particular chemical reaction for a supplied gas model.
	 * \return: Matrix of reaction constants
	 */
	void GetChemistryEquilConstants(double **RxnConstantTable, unsigned short iReaction);

	/*!
	 * \brief Provides the nMass of each species present in multi species fluid 
	 * \return: Mass of each species in Kg
	 */
	double GetParticle_Mass(unsigned short iSpecies);

	/*!
	 * \brief Provides the molar mass of each species present in multi species fluid
	 * \return: Mass of each species in Kg
	 */
	double GetMolar_Mass(unsigned short iSpecies);

	/*!
	 * \brief Provides the molar mass of each species present in multi species fluid
	 * \return: Mass of each species in Kg
	 */
	double GetMolecular_Diameter(unsigned short iSpecies);

	/*!
	 * \brief Provides the molar mass of each species present in multi species fluid
	 * \return: Molar mass of the specified gas consituent [kg/kmol]
	 */
	double GetInitial_Gas_Composition(unsigned short iSpecies);

	/*!
	 * \brief Provides the formation enthalpy of the specified species at standard conditions
	 * \return: Enthalpy of formation
	 */
	double GetEnthalpy_Formation(unsigned short iSpecies);

	/*!
	 * \brief Provides the charge number for each species present in the multi-species fluid.
	 * Charge numbers (Z) are either +1/0/-1 for positively, neutrally, or negatively ionized
	 * particles.
	 * \return: Charge number of each species (+1/0/-1)
	 */
	int GetParticle_ChargeNumber(unsigned short iSpecies);

	/*!
	 * \brief Determines whether a mean-solution restart file should be written for the plasma solution.
	 * \return: True/false boolean
	 */
	bool GetWrite_Mean_Solution(void);

	/*! 
	 * \brief Provides the Ref Temperature of each species present in multi species fluid
	 * \return: Reference Temperature for viscosity of each species in K
	 */
	double GetTemperature_Ref(unsigned short iSpecies);

	/*!
	 * \brief Provides the ref viscosity of each species present in multi species fluid
	 * \return: Reference viscosity of each species
	 */
	double GetViscosity_Ref(unsigned short iSpecies);

	/*!
	 * \brief Provides the Dipole moment of a magnet
	 * \return: dipole moment value
	 */
	double GetMagneticDipole(unsigned short iDim);

	/*!
	 * \brief Provides the restart information.
	 * \return Restart information, if <code>TRUE</code> then the code will use the solution as restart.
	 */		
	bool GetRestart(void);

	unsigned short GetnVar(void);

	/*! 
	 * \brief For some problems like adjoint or the linearized equations it 
	 *		  is necessary to restart the flow solution.
	 * \return Flow restart information, if <code>TRUE</code> then the code will restart the flow solution.
	 */

	bool GetRestart_Flow(void);

	/*!
	 * \brief For some problems like multi species  flow equations , we have block diagonal jacobians which
	 * can solved more effectively than inverting the whole big jacobian.
	 * \return block diagonal jacobian information, if <code>TRUE</code> then the code will use a modified version of Gauss Elimination
	 * for matrix inversion
	 */

	bool GetBlockDiagonalJacobian(void);

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
	 * \brief The bigrid technique is a way to filter the sensitivity obtained from the adjoint problem.
	 * \return If <code>TRUE</code> then the code will do a bigrid filtering of the sensitivity.
	 */		
	bool GetBiGrid(void);

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
	 * \brief Get the name of the restart file for the flow variables.
	 * \return Name of the restart file for the flow variables.
	 */
	string GetReStart_FlowFileName(void);

	/*! 
	 * \brief Get the name of the restart file for the linearized flow variables.
	 * \return Name of the restart file for the linearized flow variables.
	 */
	string GetReStart_LinFileName(void);

	/*! 
	 * \brief Get the name of the restart file for the adjoint variables (drag objective function).
	 * \return Name of the restart file for the adjoint variables (drag objective function).
	 */
	string GetReStart_AdjFileName(void);

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
	 * \brief Get the name of the file with the surface information for the flow problem.
	 * \return Name of the file with the surface information for the flow problem.
	 */
	string GetSurfFlowCoeff_FileName(void);

	/*! 
	 * \brief Get the name of the file with the surface information for the flow problem.
	 * \return Name of the file with the surface information for the flow problem.
	 */
	string GetSurfFlowCSV_FileName(void);

	/*! 
	 * \brief Get the name of the file with the surface information for the adjoint problem.
	 * \return Name of the file with the surface information for the adjoint problem.
	 */
	string GetSurfAdjCoeff_FileName(void);

	/*! 
	 * \brief Get the name of the file with the surface information for the adjoint problem.
	 * \return Name of the file with the surface information for the adjoint problem.
	 */
	string GetSurfAdjCSV_FileName(void);

	/*! 
	 * \brief Get the name of the file with the surface information for the linearized flow problem.
	 * \return Name of the file with the surface information for the linearized flow problem.
	 */
	string GetSurfLinCoeff_FileName(void);

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
	 * \brief Get the value of the reduced frequency.
	 * \return Value of the reduced frequency in a non-steady problem.
	 */
	double GetReduced_Frequency(void);

	/*! 
	 * \brief Get the value of the pitching amplitude.
	 * \return Value of the pitching amplitude in a non-steady problem.
	 */
	double GetPitching_Amplitude(void);

	/*! 
	 * \brief If we are prforming an unsteady simulation, there is only 
	 *        one value of the time step for the complete simulation.
	 * \return Value of the time step in an unsteady simulation (non dimensional). 
	 */
	double GetDelta_UnstTimeND(void);

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
	double GetUnst_Time(void);

	/*!
	 * \brief Flag for including mesh motion in an unsteady simulation.
	 * \return <code>TRUE</code> or <code>FALSE</code>  depending if there is unsteady mesh motion.
	 */
	bool GetUnst_Mesh_Motion(void);

	/*! 
	 * \brief Divide the rectbles and hexahedron.
	 * \return <code>TRUE</code> if the elements must be divided; otherwise <code>FALSE</code>.
	 */
	bool GetDivide_Element(void);

	/*! 
	 * \brief Value of the design variable step, we use this value in design problems.
	 * \param[in] val_dv - Number of the design variable that we want to read.
	 * \return Design variable step.
	 */	
	double GetDV_Value_New(unsigned short val_dv);

	/*! 
	 * \brief If we are doing and incremental deformation, this is the origin value.
	 * \param[in] val_dv - Number of the design variable that we want to read.
	 * \return Origin value for incremental deformations.
	 */	
	double GetDV_Value_Old(unsigned short val_dv);

	/*! 
	 * \brief Scale the value of gradient of the objective function.
	 * \return The gradient of the objective function is scaled by this value.
	 */
	double GetScale_GradOF(void);

	/*! 
	 * \brief Get information about the grid movement.
	 * \return <code>TRUE</code> if there is a grid movement; otherwise <code>FALSE</code>.
	 */
	bool GetGrid_Movement(void);

	/*! 
	 * \brief Get information about the compressible or imcompressible solver.
	 * \return <code>TRUE</code> if it is a incompressible formulation; otherwise <code>FALSE</code>.
	 */
	bool GetIncompressible(void);

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
	 * \brief Get information for multiple time stepping for plasma
	 * \return <code>TRUE</code> if it is a catalytic wall condition; otherwise <code>FALSE</code>.
	 */
	bool MultipleTimeSteps(void);

	/*!
	 * \brief Get information about the electric solver condition
	 * \return <code>TRUE</code> if it is a electric solver condition; otherwise <code>FALSE</code>.
	 */
	bool GetElectricSolver(void);
	/*! 
	 * \brief Get information about the the gravity force.
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
	unsigned short GetMarker_Config_Tag(string val_marker);

	/*! 
	 * \brief Get the boundary information (kind of boundary) in the config information of the marker <i>val_marker</i>.
	 * \return Kind of boundary in the config information of the marker <i>val_marker</i>.
	 */	
	unsigned short GetMarker_Config_Boundary(string val_marker);

	/*! 
	 * \brief Get the monitoring information in the config information of the marker <i>val_marker</i>.
	 * \return Monitoring information of the boundary in the config information of the marker <i>val_marker</i>.
	 */	
	unsigned short GetMarker_Config_Monitoring(string val_marker);

	/*! 
	 * \brief Get the plotting information in the config information of the marker <i>val_marker</i>.
	 * \return Plotting information of the boundary in the config information of the marker <i>val_marker</i>.
	 */	
	unsigned short GetMarker_Config_Plotting(string val_marker);

	/*! 
	 * \brief Get the moving information in the config information of the marker <i>val_marker</i>.
	 * \return Moving information of the boundary in the config information of the marker <i>val_marker</i>.
	 */	
	unsigned short GetMarker_Config_Moving(string val_marker);

	/*! 
	 * \brief Get the periodic information in the config information of the marker <i>val_marker</i>.
	 * \return Periodic information of the boundary in the config information of the marker <i>val_marker</i>.
	 */	
	unsigned short GetMarker_Config_PerBound(string val_marker);

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
	 * \brief Value of ther constant viscous drag for Cl/Cd computation.
	 * \return Value of ther constant viscous drag for Cl/Cd computation.
	 */
	double GetCteViscDrag(void);

	/*! 
	 * \brief Value of the origin of the rotation axis for a rotating frame problem.
	 * \return Value of the rotation axis origin.
	 */	
	double *GetRotAxisOrigin(void);

	/*! 
	 * \brief Angular velocity vector for a rotating frame problem.
	 * \return The specified angular velocity vector.
	 */	
	double *GetOmega(void);

	/*!
	 * \brief Angular velocity magnitude for a rotating frame problem.
	 * \return The specified angular velocity magnitude.
	 */	
	double GetOmegaMag(void);

	/*!
	 * \brief Reference length (m, e.g. rotor radius) for computing force coefficients in a rotating frame.
	 * \return The specified reference length (m, e.g. rotor radius).
	 */	
	double GetRotRadius(void);

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
	void SetGlobalParam(unsigned short val_solver, unsigned short val_system);

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
	 * \brief Get information about converting a mesh from CGNS to SU2 format.
	 * \return <code>TRUE</code> if a conversion is requested; otherwise <code>FALSE</code>.
	 */
	bool GetCGNS_To_SU2(void);

	/*! 
	 * \brief Get information about whether a converted mesh should be written.
	 * \return <code>TRUE</code> if the converted mesh should be written; otherwise <code>FALSE</code>.
	 */
	bool GetWrite_Converted_Mesh(void);

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
	 * \brief Get the total temperature at an inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The total temperature.
	 */	

	double GetInlet_Ttotal(string val_index);

	/*!
	 * \brief Get the total pressure at an inlet boundary.
	 * \param[in] val_index - Index corresponding to the inlet boundary.
	 * \return The total pressure.
	 */
	double GetInlet_Ptotal(string val_index);

	/*!
	 * \brief Value of the CFL reduction in LevelSet problems.
	 * \return Value of the CFL reduction in LevelSet problems.
	 */
	double GetLevelSet_CFLRedCoeff(void);

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
	 * \brief Set the non-dimensionalization for SU2_CFD.
	 * \param[in] val_nDim - Number of dimensions for this particular problem.
	 * \param[in] val_rank - Processor rank.
	 * \param[in] val_iDomain - Current grid domain number.
	 */	
	void SetNondimensionalization(unsigned short val_nDim, int val_rank, unsigned short val_iDomain);

};

#include "config_structure.inl"
