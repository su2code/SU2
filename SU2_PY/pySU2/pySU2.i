/*
################################################################################
#
# \file pySU2.i
# \brief Configuration file for the Swig compilation of the Python wrapper.
# \author D. Thomas
# \version 5.0.0 "Raven"
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
*/

%feature("autodoc","1");

%module(docstring=
"'pysu2' module",
directors="1",
threads="1"
) pysu2
%{

#include "../../SU2_CFD/include/driver_structure.hpp"

%}

// ----------- USED MODULES ------------
%import "../../Common/include/datatypes/primitive_structure.hpp"
%import "../../Common/include/mpi_structure.hpp"
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "typemaps.i"
//%include "numpy.i"
#ifdef HAVE_MPI                    //Need mpi4py only for a parallel build of the wrapper.
  %include "mpi4py/mpi4py.i"
  %mpi4py_typemap(Comm, MPI_Comm)
#endif

namespace std {
   %template() vector<int>;
   %template() vector<double>;
   %template() vector<string>;
   %template() map<string, int>;
   %template() map<string, string>;
}

// ----------- API CLASSES ----------------

//Constants definitions
/*!
 * \brief different software components of SU2
 */
enum SU2_COMPONENT {
  SU2_CFD = 1,	/*!< \brief Running the SU2_CFD software. */
  SU2_DEF = 2,	/*!< \brief Running the SU2_DEF software. */
  SU2_DOT = 3,	/*!< \brief Running the SU2_DOT software. */
  SU2_MSH = 4,	/*!< \brief Running the SU2_MSH software. */
  SU2_GEO = 5,	/*!< \brief Running the SU2_GEO software. */
  SU2_SOL = 6 	/*!< \brief Running the SU2_SOL software. */
};

const unsigned int MESH_0 = 0; /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1; /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0; /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1; /*!< \brief Definition of the first grid domain. */

// CDriver class
%include "../../SU2_CFD/include/driver_structure.hpp"

// CConfig class
class CConfig {
private:
  SU2_Comm SU2_Communicator; /*!< \brief MPI communicator of SU2.*/
  int rank;
  unsigned short Kind_SU2; /*!< \brief Kind of SU2 software component.*/
  unsigned short Ref_NonDim; /*!< \brief Kind of non dimensionalization.*/
  unsigned short Kind_MixingProcess; /*!< \brief Kind of mixing process.*/
  unsigned short *Kind_TurboPerformance; /*!< \brief Kind of Turbomachinery performance calculation.*/
  unsigned short iZone, nZone; /*!< \brief Number of zones in the mesh. */
  su2double Highlite_Area; /*!< \brief Highlite area. */
  su2double Fan_Poly_Eff; /*!< \brief Highlite area. */
  su2double OrderMagResidual; /*!< \brief Order of magnitude reduction. */
  su2double MinLogResidual; /*!< \brief Minimum value of the log residual. */
  su2double OrderMagResidualFSI; /*!< \brief Order of magnitude reduction. */
  su2double MinLogResidualFSI; /*!< \brief Minimum value of the log residual. */
  su2double Res_FEM_UTOL; 		/*!< \brief UTOL criteria for structural FEM. */
  su2double Res_FEM_RTOL; 		/*!< \brief RTOL criteria for structural FEM. */
  su2double Res_FEM_ETOL; 		/*!< \brief ETOL criteria for structural FEM. */
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
  su2double FFD_Scale;  	/*!< \brief Scale factor between the design variable value and the control point movement. */
  bool Viscous_Limiter_Flow, Viscous_Limiter_Turb;			/*!< \brief Viscous limiters. */
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
  Rotating_Frame,			/*!< \brief Flag to know if there is a rotating frame. */
  PoissonSolver,			/*!< \brief Flag to know if we are solving  poisson forces  in plasma solver. */
  Low_Mach_Precon,		/*!< \brief Flag to know if we are using a low Mach number preconditioner. */
  Low_Mach_Corr,			/*!< \brief Flag to know if we are using a low Mach number correction. */
  GravityForce,			/*!< \brief Flag to know if the gravity force is incuded in the formulation. */
  SmoothNumGrid,			/*!< \brief Smooth the numerical grid. */
  AdaptBoundary,			/*!< \brief Adapt the elements on the boundary. */
  SubsonicEngine,			/*!< \brief Engine intake subsonic region. */
  Frozen_Visc,			/*!< \brief Flag for adjoint problem with/without frozen viscosity. */
  Sens_Remove_Sharp,			/*!< \brief Flag for removing or not the sharp edges from the sensitivity computation. */
  Hold_GridFixed,	/*!< \brief Flag hold fixed some part of the mesh during the deformation. */
  Axisymmetric; /*!< \brief Flag for axisymmetric calculations */
  su2double Damp_Engine_Inflow;	/*!< \brief Damping factor for the engine inlet. */
  su2double Damp_Engine_Exhaust;	/*!< \brief Damping factor for the engine exhaust. */
  su2double Damp_Res_Restric,	/*!< \brief Damping factor for the residual restriction. */
  Damp_Correc_Prolong; /*!< \brief Damping factor for the correction prolongation. */
  su2double Position_Plane; /*!< \brief Position of the Near-Field (y coordinate 2D, and z coordinate 3D). */
  su2double WeightCd; /*!< \brief Weight of the drag coefficient. */
  su2double dCD_dCL; /*!< \brief Weight of the drag coefficient. */
  su2double dCD_dCM; /*!< \brief Weight of the drag coefficient. */
  su2double CL_Target; /*!< \brief Weight of the drag coefficient. */
  su2double CM_Target; /*!< \brief Weight of the drag coefficient. */
  su2double *HTP_Min_XCoord, *HTP_Min_YCoord; /*!< \brief Identification of the HTP. */
  unsigned short Unsteady_Simulation;	/*!< \brief Steady or unsteady (time stepping or dual time stepping) computation. */
  unsigned short Dynamic_Analysis;	/*!< \brief Static or dynamic structural analysis. */
  unsigned short nStartUpIter;	/*!< \brief Start up iterations using the fine grid. */
  su2double FixAzimuthalLine; /*!< \brief Fix an azimuthal line due to misalignments of the nearfield. */
  su2double **DV_Value;		/*!< \brief Previous value of the design variable. */
  su2double DVBound_Upper;		/*!< \brief Previous value of the design variable. */
  su2double DVBound_Lower;		/*!< \brief Previous value of the design variable. */
  su2double LimiterCoeff;				/*!< \brief Limiter coefficient */
  unsigned long LimiterIter;	/*!< \brief Freeze the value of the limiter after a number of iterations */
  su2double SharpEdgesCoeff;				/*!< \brief Coefficient to identify the limit of a sharp edge. */
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
  nMarker_Pressure,				/*!< \brief Number of pressure wall markers. */
  nMarker_PerBound,				/*!< \brief Number of periodic boundary markers. */
  nMarker_MixBound,				/*!< \brief Number of mixing boundary markers. */
  nMarker_TurboPerf,				/*!< \brief Number of mixing boundary markers. */
  nMarker_NearFieldBound,				/*!< \brief Number of near field boundary markers. */
  nMarker_ActDiskInlet, nMarker_ActDiskOutlet,
  nMarker_InterfaceBound,				/*!< \brief Number of interface boundary markers. */
  nMarker_Fluid_InterfaceBound,				/*!< \brief Number of fluid interface markers. */
  nMarker_Dirichlet,				/*!< \brief Number of interface boundary markers. */
  nMarker_Inlet,					/*!< \brief Number of inlet flow markers. */
  nMarker_Riemann,					/*!< \brief Number of Riemann flow markers. */
  nMarker_NRBC,					/*!< \brief Number of NRBC flow markers. */
  nMarker_Supersonic_Inlet,					/*!< \brief Number of supersonic inlet flow markers. */
  nMarker_Supersonic_Outlet,					/*!< \brief Number of supersonic outlet flow markers. */
  nMarker_Outlet,					/*!< \brief Number of outlet flow markers. */
  nMarker_Out_1D,         /*!< \brief Number of outlet flow markers over which to calculate 1D outputs */
  nMarker_Isothermal,     /*!< \brief Number of isothermal wall boundaries. */
  nMarker_HeatFlux,       /*!< \brief Number of constant heat flux wall boundaries. */
  nMarker_EngineExhaust,					/*!< \brief Number of nacelle exhaust flow markers. */
  nMarker_EngineInflow,					/*!< \brief Number of nacelle inflow flow markers. */
  nMarker_Clamped,						/*!< \brief Number of clamped markers in the FEM. */
  nMarker_Displacement,					/*!< \brief Number of displacement surface markers. */
  nMarker_Load,					/*!< \brief Number of load surface markers. */
  nMarker_Load_Dir,					/*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_Load_Sine,					/*!< \brief Number of load surface markers defined by magnitude and direction. */
  nMarker_FlowLoad,					/*!< \brief Number of load surface markers. */
  nMarker_Neumann,				/*!< \brief Number of Neumann flow markers. */
  nMarker_Internal,				/*!< \brief Number of Neumann flow markers. */
  nMarker_All,					/*!< \brief Total number of markers using the grid information. */
  nMarker_Max,					/*!< \brief Max number of number of markers using the grid information. */
  nMarker_CfgFile;					/*!< \brief Total number of markers using the config file
                             (note that using parallel computation this number can be different
                             from nMarker_All). */
  string *Marker_Euler,			/*!< \brief Euler wall markers. */
  *Marker_FarField,				/*!< \brief Far field markers. */
  *Marker_Custom,
  *Marker_SymWall,				/*!< \brief Symmetry wall markers. */
  *Marker_Pressure,				/*!< \brief Pressure boundary markers. */
  *Marker_PerBound,				/*!< \brief Periodic boundary markers. */
  *Marker_PerDonor,				/*!< \brief Rotationally periodic boundary donor markers. */
  *Marker_MixBound,				/*!< \brief MixingPlane boundary markers. */
  *Marker_MixDonor,				/*!< \brief MixingPlane boundary donor markers. */
  *Marker_TurboBoundIn,				/*!< \brief Turbomachinery performance boundary markers. */
  *Marker_TurboBoundOut,				/*!< \brief Turbomachinery performance boundary donor markers. */
  *Marker_NearFieldBound,				/*!< \brief Near Field boundaries markers. */
  *Marker_InterfaceBound,				/*!< \brief Interface boundaries markers. */
  *Marker_Fluid_InterfaceBound,				/*!< \brief Fluid interface markers. */
  *Marker_ActDiskInlet,
  *Marker_ActDiskOutlet,
  *Marker_Dirichlet,				/*!< \brief Interface boundaries markers. */
  *Marker_Inlet,					/*!< \brief Inlet flow markers. */
  *Marker_Riemann,					/*!< \brief Riemann markers. */
  *Marker_NRBC,					/*!< \brief NRBC markers. */
  *Marker_Supersonic_Inlet,					/*!< \brief Supersonic inlet flow markers. */
  *Marker_Supersonic_Outlet,					/*!< \brief Supersonic outlet flow markers. */
  *Marker_Outlet,					/*!< \brief Outlet flow markers. */
  *Marker_Out_1D,         /*!< \brief Outlet flow markers over which to calculate 1D output. */
  *Marker_Isothermal,     /*!< \brief Isothermal wall markers. */
  *Marker_HeatFlux,       /*!< \brief Constant heat flux wall markers. */
  *Marker_EngineInflow,					/*!< \brief Engine Inflow flow markers. */
  *Marker_EngineExhaust,					/*!< \brief Engine Exhaust flow markers. */
  *Marker_Clamped,						/*!< \brief Clamped markers. */
  *Marker_Displacement,					/*!< \brief Displacement markers. */
  *Marker_Load,					/*!< \brief Load markers. */
  *Marker_Load_Dir,					/*!< \brief Load markers defined in cartesian coordinates. */
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
  su2double *NRBC_Var1, *NRBC_Var2;    /*!< \brief Specified values for NRBC boundary. */
  su2double **NRBC_FlowDir;  /*!< \brief Specified flow direction vector (unit vector) for NRBC boundaries. */
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
  su2double *Load_Dir_Value;    /*!< \brief Specified force for load boundaries defined in cartesian coordinates. */
  su2double *Load_Dir_Multiplier;    /*!< \brief Specified multiplier for load boundaries defined in cartesian coordinates. */
  su2double **Load_Dir;  /*!< \brief Specified flow direction vector (unit vector) for inlet boundaries. */
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
  su2double *Surface_MassFlow;    /*!< \brief Specified fan face mach for nacelle boundaries. */
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
  unsigned long FSIIter;			/*!< \brief Current Fluid Structure Interaction sub-iteration number. */
  unsigned long Unst_nIntIter;			/*!< \brief Number of internal iterations (Dual time Method). */
  unsigned long Dyn_nIntIter;			/*!< \brief Number of internal iterations (Newton-Raphson Method for nonlinear structural analysis). */
  long Unst_RestartIter;			/*!< \brief Iteration number to restart an unsteady simulation (Dual time Method). */
  long Unst_AdjointIter;			/*!< \brief Iteration number to begin the reverse time integration in the direct solver for the unsteady adjoint. */
  long Iter_Avg_Objective;			/*!< \brief Iteration the number of time steps to be averaged, counting from the back */
  long Dyn_RestartIter;			/*!< \brief Iteration number to restart a dynamic structural analysis. */
  unsigned short nRKStep;			/*!< \brief Number of steps of the explicit Runge-Kutta method. */
  su2double *RK_Alpha_Step;			/*!< \brief Runge-Kutta beta coefficients. */
  unsigned short nMGLevels;		/*!< \brief Number of multigrid levels (coarse levels). */
  unsigned short nCFL;			/*!< \brief Number of CFL, one for each multigrid level. */
  su2double
  CFLRedCoeff_Turb,		/*!< \brief CFL reduction coefficient on the LevelSet problem. */
  CFLRedCoeff_AdjFlow,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
  CFLRedCoeff_AdjTurb,	/*!< \brief CFL reduction coefficient for the adjoint problem. */
  CFLFineGrid,		/*!< \brief CFL of the finest grid. */
  Max_DeltaTime,  		/*!< \brief Max delta time. */
  Unst_CFL;		/*!< \brief Unsteady CFL number. */
  bool AddIndNeighbor;			/*!< \brief Include indirect neighbor in the agglomeration process. */
  unsigned short nDV,		/*!< \brief Number of design variables. */
  nObj, nObjW;              /*! \brief Number of objective functions. */
  unsigned short* nDV_Value;		/*!< \brief Number of values for each design variable (might be different than 1 if we allow arbitrary movement). */
  unsigned short nFFDBox;		/*!< \brief Number of ffd boxes. */
  unsigned short nGridMovement;		/*!< \brief Number of grid movement types specified. */
  unsigned short nParamDV;		/*!< \brief Number of parameters of the design variable. */
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
  unsigned short Kind_Solver,	/*!< \brief Kind of solver Euler, NS, Continuous adjoint, etc.  */
  Kind_FluidModel,			/*!< \brief Kind of the Fluid Model: Ideal or Van der Walls, ... . */
  Kind_ViscosityModel,			/*!< \brief Kind of the Viscosity Model*/
  Kind_ConductivityModel,			/*!< \brief Kind of the Thermal Conductivity Model*/
  Kind_FreeStreamOption,			/*!< \brief Kind of free stream option to choose if initializing with density or temperature  */
  Kind_InitOption,			/*!< \brief Kind of Init option to choose if initializing with Reynolds number or with thermodynamic conditions   */
  Kind_GasModel,				/*!< \brief Kind of the Gas Model. */
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
  Kind_SlopeLimit,				/*!< \brief Global slope limiter. */
  Kind_SlopeLimit_Flow,		/*!< \brief Slope limiter for flow equations.*/
  Kind_SlopeLimit_Turb,		/*!< \brief Slope limiter for the turbulence equation.*/
  Kind_SlopeLimit_AdjTurb,	/*!< \brief Slope limiter for the adjoint turbulent equation.*/
  Kind_SlopeLimit_AdjFlow,	/*!< \brief Slope limiter for the adjoint equation.*/
  Kind_TimeNumScheme,			/*!< \brief Global explicit or implicit time integration. */
  Kind_TimeIntScheme_Flow,	/*!< \brief Time integration for the flow equations. */
  Kind_TimeIntScheme_AdjFlow,		/*!< \brief Time integration for the adjoint flow equations. */
  Kind_TimeIntScheme_Turb,	/*!< \brief Time integration for the turbulence model. */
  Kind_TimeIntScheme_AdjTurb,	/*!< \brief Time integration for the adjoint turbulence model. */
  Kind_TimeIntScheme_Wave,	/*!< \brief Time integration for the wave equations. */
  Kind_TimeIntScheme_Heat,	/*!< \brief Time integration for the wave equations. */
  Kind_TimeIntScheme_Poisson,	/*!< \brief Time integration for the wave equations. */
  Kind_TimeIntScheme_FEA,	/*!< \brief Time integration for the FEA equations. */
  Kind_SpaceIteScheme_FEA,	/*!< \brief Iterative scheme for nonlinear structural analysis. */
  Kind_ConvNumScheme,			/*!< \brief Global definition of the convective term. */
  Kind_ConvNumScheme_Flow,	/*!< \brief Centered or upwind scheme for the flow equations. */
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
  Kind_Solver_Fluid_FSI,		/*!< \brief Kind of solver for the fluid in FSI applications. */
  Kind_Solver_Struc_FSI,		/*!< \brief Kind of solver for the structure in FSI applications. */
  Kind_BGS_RelaxMethod,				/*!< \brief Kind of relaxation method for Block Gauss Seidel method in FSI problems. */
  Kind_TransferMethod,	/*!< \brief Iterative scheme for nonlinear structural analysis. */
  SpatialOrder,		/*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_Flow,		/*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_Turb,		/*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_AdjFlow,		/*!< \brief Order of the spatial numerical integration.*/
  SpatialOrder_AdjTurb;		/*!< \brief Order of the spatial numerical integration.*/
  bool FSI_Problem;			/*!< \brief Boolean to determine whether the simulation is FSI or not. */
  bool AD_Mode;         /*!< \brief Algorithmic Differentiation support. */
  unsigned short Kind_Material_Compress,	/*!< \brief Determines if the material is compressible or incompressible (structural analysis). */
  Kind_Material,			/*!< \brief Determines the material model to be used (structural analysis). */
  Kind_Struct_Solver;		/*!< \brief Determines the geometric condition (small or large deformations) for structural analysis. */
  unsigned short Kind_Turb_Model;			/*!< \brief Turbulent model definition. */
  unsigned short Kind_Trans_Model,			/*!< \brief Transition model definition. */
  Kind_ActDisk, Kind_Engine_Inflow, Kind_Inlet, *Kind_Data_Riemann, *Kind_Data_NRBC;           /*!< \brief Kind of inlet boundary treatment. */
  su2double Linear_Solver_Error;		/*!< \brief Min error of the linear solver for the implicit formulation. */
  su2double Linear_Solver_Error_FSI_Struc;		/*!< \brief Min error of the linear solver for the implicit formulation in the structural side for FSI problems . */
  unsigned long Linear_Solver_Iter;		/*!< \brief Max iterations of the linear solver for the implicit formulation. */
  unsigned long Linear_Solver_Iter_FSI_Struc;		/*!< \brief Max iterations of the linear solver for FSI applications and structural solver. */
  unsigned long Linear_Solver_Restart_Frequency;   /*!< \brief Restart frequency of the linear solver for the implicit formulation. */
  su2double SemiSpan;		/*!< \brief Wing Semi span. */
  su2double Roe_Kappa;		/*!< \brief Relaxation of the Roe scheme. */
  su2double Relaxation_Factor_Flow;		/*!< \brief Relaxation coefficient of the linear solver mean flow. */
  su2double Relaxation_Factor_Turb;		/*!< \brief Relaxation coefficient of the linear solver turbulence. */
  su2double Relaxation_Factor_AdjFlow;		/*!< \brief Relaxation coefficient of the linear solver adjoint mean flow. */
  su2double AdjTurb_Linear_Error;		/*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  su2double EntropyFix_Coeff;              /*!< \brief Entropy fix coefficient. */
  unsigned short AdjTurb_Linear_Iter;		/*!< \brief Min error of the turbulent adjoint linear solver for the implicit formulation. */
  su2double *Section_WingBounds;                  /*!< \brief Airfoil section limit. */
  unsigned short nLocationStations,      /*!< \brief Number of section cuts to make when outputting mesh and cp . */
  nWingStations;               /*!< \brief Number of section cuts to make when calculating internal volume. */
  su2double* Kappa_Flow,           /*!< \brief Numerical dissipation coefficients for the flow equations. */
  *Kappa_AdjFlow;                  /*!< \brief Numerical dissipation coefficients for the linearized equations. */
  su2double* FFD_Axis;       /*!< \brief Numerical dissipation coefficients for the adjoint equations. */
  su2double Kappa_1st_AdjFlow,	/*!< \brief JST 1st order dissipation coefficient for adjoint flow equations (coarse multigrid levels). */
  Kappa_2nd_AdjFlow,			/*!< \brief JST 2nd order dissipation coefficient for adjoint flow equations. */
  Kappa_4th_AdjFlow,			/*!< \brief JST 4th order dissipation coefficient for adjoint flow equations. */
  Kappa_1st_Flow,			/*!< \brief JST 1st order dissipation coefficient for flow equations (coarse multigrid levels). */
  Kappa_2nd_Flow,			/*!< \brief JST 2nd order dissipation coefficient for flow equations. */
  Kappa_4th_Flow;			/*!< \brief JST 4th order dissipation coefficient for flow equations. */

  su2double Min_Beta_RoeTurkel,		/*!< \brief Minimum value of Beta for the Roe-Turkel low Mach preconditioner. */
  Max_Beta_RoeTurkel;		/*!< \brief Maximum value of Beta for the Roe-Turkel low Mach preconditioner. */
  unsigned long GridDef_Nonlinear_Iter, /*!< \brief Number of nonlinear increments for grid deformation. */
  GridDef_Linear_Iter; /*!< \brief Number of linear smoothing iterations for grid deformation. */
  unsigned short Deform_Stiffness_Type; /*!< \brief Type of element stiffness imposed for FEA mesh deformation. */
  bool Deform_Output;  /*!< \brief Print the residuals during mesh deformation to the console. */
  su2double Deform_Tol_Factor; /*!< Factor to multiply smallest volume for deform tolerance (0.001 default) */
  su2double Deform_Coeff; /*!< Deform coeffienct */
  unsigned short FFD_Continuity; /*!< Surface continuity at the intersection with the FFD */
  unsigned short FFD_CoordSystem; /*!< Define the coordinates system */
  su2double Deform_ElasticityMod, Deform_PoissonRatio; /*!< young's modulus and poisson ratio for volume deformation stiffness model */
  bool Visualize_Deformation;	/*!< \brief Flag to visualize the deformation in MDC. */
  bool FFD_Symmetry_Plane;	/*!< \brief FFD symmetry plane. */
  su2double Mach;		/*!< \brief Mach number. */
  su2double Reynolds;	/*!< \brief Reynolds number. */
  su2double Froude;	/*!< \brief Froude number. */
  su2double Length_Reynolds;	/*!< \brief Reynolds length (dimensional). */
  su2double AoA,			/*!< \brief Angle of attack (just external flow). */
  iH, AoS, AoA_Offset, AoS_Offset, AoA_Sens;		/*!< \brief Angle of sideSlip (just external flow). */
  bool Fixed_CL_Mode;			/*!< \brief Activate fixed CL mode (external flow only). */
  bool Fixed_CM_Mode;			/*!< \brief Activate fixed CL mode (external flow only). */
  bool Eval_dCD_dCX;			/*!< \brief Activate fixed CL mode (external flow only). */
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
  unsigned long StartConv_Iter;	/*!< \brief Start convergence criteria at iteration. */
  su2double Cauchy_Eps;	/*!< \brief Epsilon used for the convergence. */
  unsigned long Wrt_Sol_Freq,	/*!< \brief Writing solution frequency. */
  Wrt_Sol_Freq_DualTime,	/*!< \brief Writing solution frequency for Dual Time. */
  Wrt_Con_Freq,				/*!< \brief Writing convergence history frequency. */
  Wrt_Con_Freq_DualTime;				/*!< \brief Writing convergence history frequency. */
  bool Wrt_Unsteady;  /*!< \brief Write unsteady data adding header and prefix. */
  bool Wrt_Dynamic;  		/*!< \brief Write dynamic data adding header and prefix. */
  bool LowFidelitySim;  /*!< \brief Compute a low fidelity simulation. */
  bool Restart,	/*!< \brief Restart solution (for direct, adjoint, and linearized problems).*/
  Restart_Flow;	/*!< \brief Restart flow solution for adjoint and linearized problems. */
  unsigned short nMarker_Monitoring,	/*!< \brief Number of markers to monitor. */
  nMarker_Designing,					/*!< \brief Number of markers for the objective function. */
  nMarker_GeoEval,					/*!< \brief Number of markers for the objective function. */
  nMarker_Plotting,					/*!< \brief Number of markers to plot. */
  nMarker_Analyze,					/*!< \brief Number of markers to plot. */
  nMarker_FSIinterface,					/*!< \brief Number of markers in the FSI interface. */
  nMarker_Moving,               /*!< \brief Number of markers in motion (DEFORMING, MOVING_WALL, or FLUID_STRUCTURE). */
  nMarker_DV;               /*!< \brief Number of markers affected by the design variables. */
  string *Marker_Monitoring,     /*!< \brief Markers to monitor. */
  *Marker_Designing,         /*!< \brief Markers to plot. */
  *Marker_GeoEval,         /*!< \brief Markers to plot. */
  *Marker_Plotting,          /*!< \brief Markers to plot. */
  *Marker_Analyze,          /*!< \brief Markers to plot. */
  *Marker_FSIinterface,          /*!< \brief Markers in the FSI interface. */
  *Marker_Moving,            /*!< \brief Markers in motion (DEFORMING, MOVING_WALL, or FLUID_STRUCTURE). */
  *Marker_DV;            /*!< \brief Markers affected by the design variables. */
  unsigned short  *Marker_All_Monitoring,        /*!< \brief Global index for monitoring using the grid information. */
  *Marker_All_GeoEval,       /*!< \brief Global index for geometrical evaluation. */
  *Marker_All_Plotting,        /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_Analyze,        /*!< \brief Global index for plotting using the grid information. */
  *Marker_All_FSIinterface,        /*!< \brief Global index for FSI interface markers using the grid information. */
  *Marker_All_DV,          /*!< \brief Global index for design variable markers using the grid information. */
  *Marker_All_Moving,          /*!< \brief Global index for moving surfaces using the grid information. */
  *Marker_All_Designing,         /*!< \brief Global index for moving using the grid information. */
  *Marker_All_Out_1D,      /*!< \brief Global index for moving using 1D integrated output. */
  *Marker_CfgFile_Monitoring,     /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Designing,      /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_GeoEval,      /*!< \brief Global index for monitoring using the config information. */
  *Marker_CfgFile_Plotting,     /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_Analyze,     /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_FSIinterface,     /*!< \brief Global index for FSI interface using the config information. */
  *Marker_CfgFile_Out_1D,      /*!< \brief Global index for plotting using the config information. */
  *Marker_CfgFile_Moving,       /*!< \brief Global index for moving surfaces using the config information. */
  *Marker_CfgFile_DV,       /*!< \brief Global index for design variable markers using the config information. */
  *Marker_CfgFile_PerBound;     /*!< \brief Global index for periodic boundaries using the config information. */
  string *PlaneTag;      /*!< \brief Global index for the plane adaptation (upper, lower). */
  su2double DualVol_Power;			/*!< \brief Power for the dual volume in the grid adaptation sensor. */
  unsigned short Analytical_Surface;	/*!< \brief Information about the analytical definition of the surface for grid adaptation. */
  unsigned short Axis_Stations;	/*!< \brief Axis orientation. */
  unsigned short Mesh_FileFormat;	/*!< \brief Mesh input format. */
  unsigned short Output_FileFormat;	/*!< \brief Format of the output files. */
  unsigned short ActDisk_Jump;	/*!< \brief Format of the output files. */
  bool CFL_Adapt;      /*!< \brief Adaptive CFL number. */
  su2double RefAreaCoeff,		/*!< \brief Reference area for coefficient computation. */
  RefElemLength,				/*!< \brief Reference element length for computing the slope limiting epsilon. */
  RefSharpEdges,				/*!< \brief Reference coefficient for detecting sharp edges. */
  RefLengthMoment,			/*!< \brief Reference length for moment computation. */
  *RefOriginMoment,           /*!< \brief Origin for moment computation. */
  *RefOriginMoment_X,      /*!< \brief X Origin for moment computation. */
  *RefOriginMoment_Y,      /*!< \brief Y Origin for moment computation. */
  *RefOriginMoment_Z,      /*!< \brief Z Origin for moment computation. */
  *CFL_AdaptParam,      /*!< \brief Information about the CFL ramp. */
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
  Solution_FEMFileName,			/*!< \brief Adjoint solution input file for drag functional. */
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
  Breakdown_FileName,			    /*!< \brief Breakdown output file. */
  Conv_FileName_FSI,					/*!< \brief Convergence history output file. */
  Restart_FlowFileName,			/*!< \brief Restart file for flow variables. */
  Restart_WaveFileName,			/*!< \brief Restart file for wave variables. */
  Restart_HeatFileName,			/*!< \brief Restart file for heat variables. */
  Restart_AdjFileName,			/*!< \brief Restart file for adjoint variables, drag functional. */
  Restart_FEMFileName,			/*!< \brief Restart file for FEM elasticity. */
  Adj_FileName,					/*!< \brief Output file with the adjoint variables. */
  ObjFunc_Grad_FileName,			/*!< \brief Gradient of the objective function. */
  ObjFunc_Value_FileName,			/*!< \brief Objective function. */
  SurfFlowCoeff_FileName,			/*!< \brief Output file with the flow variables on the surface. */
  SurfAdjCoeff_FileName,			/*!< \brief Output file with the adjoint variables on the surface. */
  New_SU2_FileName,       		/*!< \brief Output SU2 mesh file converted from CGNS format. */
  SurfSens_FileName,			/*!< \brief Output file for the sensitivity on the surface (discrete adjoint). */
  VolSens_FileName;			/*!< \brief Output file for the sensitivity in the volume (discrete adjoint). */
  bool Low_MemoryOutput,      /*!< \brief Write a volume solution file */
  Wrt_Vol_Sol,                /*!< \brief Write a volume solution file */
  Wrt_Srf_Sol,                /*!< \brief Write a surface solution file */
  Wrt_Csv_Sol,                /*!< \brief Write a surface comma-separated values solution file */
  Wrt_Residuals,              /*!< \brief Write residuals to solution file */
  Wrt_Limiters,              /*!< \brief Write residuals to solution file */
  Wrt_SharpEdges,              /*!< \brief Write residuals to solution file */
  Wrt_Halo,                   /*!< \brief Write rind layers in solution files */
  Plot_Section_Forces,       /*!< \brief Write sectional forces for specified markers. */
  Wrt_1D_Output;                /*!< \brief Write average stagnation pressure specified markers. */
  unsigned short Console_Output_Verb;  /*!< \brief Level of verbosity for console output */
  su2double Gamma,			/*!< \brief Ratio of specific heats of the gas. */
  Bulk_Modulus,			/*!< \brief Value of the bulk modulus for incompressible flows. */
  ArtComp_Factor,			/*!< \brief Value of the artificial compresibility factor for incompressible flows. */
  Gas_Constant,     /*!< \brief Specific gas constant. */
  Gas_ConstantND,     /*!< \brief Non-dimensional specific gas constant. */
  Gas_Constant_Ref, /*!< \brief Reference specific gas constant. */
  Temperature_Critical,   /*!< \brief Critical Temperature for real fluid model.  */
  Pressure_Critical,   /*!< \brief Critical Pressure for real fluid model.  */
  Density_Critical,   /*!< \brief Critical Density for real fluid model.  */
  Acentric_Factor,   /*!< \brief Acentric Factor for real fluid model.  */
  Mu_ConstantND,   /*!< \brief Constant Viscosity for ConstantViscosity model.  */
  Kt_ConstantND,   /*!< \brief Constant Thermal Conductivity for ConstantConductivity model.  */
  Mu_RefND,   /*!< \brief reference viscosity for Sutherland model.  */
  Mu_Temperature_RefND,   /*!< \brief reference Temperature for Sutherland model.  */
  Mu_SND,   /*!< \brief reference S for Sutherland model.  */
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
  NuFactor_Engine,  /*!< \brief Ratio of turbulent to laminar viscosity at the engine. */
  SecondaryFlow_ActDisk,  /*!< \brief Ratio of turbulent to laminar viscosity at the actuator disk. */
  Initial_BCThrust,  /*!< \brief Ratio of turbulent to laminar viscosity at the actuator disk. */
  Pressure_FreeStream,     /*!< \brief Total pressure of the fluid. */
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
  Temperature_FreeStreamND,   /*!< \brief Farfield temperature value (external flow). */
  Density_FreeStreamND,       /*!< \brief Farfield density value (external flow). */
  Velocity_FreeStreamND[3],   /*!< \brief Farfield velocity values (external flow). */
  Energy_FreeStreamND,        /*!< \brief Farfield energy value (external flow). */
  Viscosity_FreeStreamND,     /*!< \brief Farfield viscosity value (external flow). */
  Tke_FreeStreamND,           /*!< \brief Farfield kinetic energy (external flow). */
  Omega_FreeStreamND,         /*!< \brief Specific dissipation (external flow). */
  Omega_FreeStream;           /*!< \brief Specific dissipation (external flow). */
  su2double ElasticyMod,			/*!< \brief Young's modulus of elasticity. */
  PoissonRatio,						    /*!< \brief Poisson's ratio. */
  MaterialDensity,					  /*!< \brief Material density. */
  Bulk_Modulus_Struct;				/*!< \brief Bulk modulus (on the structural side). */
  unsigned short Kind_2DElasForm;			/*!< \brief Kind of bidimensional elasticity solver. */
  unsigned short nIterFSI;	  /*!< \brief Number of maximum number of subiterations in a FSI problem. */
  su2double AitkenStatRelax;	/*!< \brief Aitken's relaxation factor (if set as static) */
  su2double AitkenDynMaxInit;	/*!< \brief Aitken's maximum dynamic relaxation factor for the first iteration */
  su2double AitkenDynMinInit;	/*!< \brief Aitken's minimum dynamic relaxation factor for the first iteration */
  su2double Wave_Speed;			  /*!< \brief Wave speed used in the wave solver. */
  su2double Thermal_Diffusivity;			/*!< \brief Thermal diffusivity used in the heat solver. */
  su2double Cyclic_Pitch,     /*!< \brief Cyclic pitch for rotorcraft simulations. */
  Collective_Pitch;           /*!< \brief Collective pitch for rotorcraft simulations. */
  string Motion_Filename;			/*!< \brief Arbitrary mesh motion input base filename. */
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
  bool DeadLoad; 	          	/*!< Application of dead loads to the FE analysis */
  bool MatchingMesh; 	        /*!< Matching mesh (while implementing interpolation procedures). */
  bool SteadyRestart; 	      /*!< Restart from a steady state for FSI problems. */
  su2double Newmark_alpha,		/*!< \brief Parameter alpha for Newmark method. */
  Newmark_delta;				      /*!< \brief Parameter delta for Newmark method. */
  unsigned short nIntCoeffs;	/*!< \brief Number of integration coeffs for structural calculations. */
  su2double *Int_Coeffs;		  /*!< \brief Time integration coefficients for structural method. */
  bool Sigmoid_Load,		      /*!< \brief Apply the load using a sigmoid. */
  Ramp_Load;				          /*!< \brief Apply the load with linear increases. */
  bool IncrementalLoad;		    /*!< \brief Apply the load in increments (for nonlinear structural analysis). */
  unsigned long IncLoad_Nincrements; /*!< \brief Number of increments. */
  su2double *IncLoad_Criteria;/*!< \brief Criteria for the application of incremental loading. */
  su2double Ramp_Time;			  /*!< \brief Time until the maximum load is applied. */
  su2double Sigmoid_Time;			/*!< \brief Time until the maximum load is applied, using a sigmoid. */
  su2double Sigmoid_K;			  /*!< \brief Sigmoid parameter determining its steepness. */
  su2double Static_Time;			/*!< \brief Time while the structure is not loaded in FSI applications. */
  unsigned short Pred_Order;  /*!< \brief Order of the predictor for FSI applications. */
  unsigned short Kind_Interpolation; /*!\brief type of interpolation to use for FSI applications. */
  bool Prestretch;            /*!< Read a reference geometry for optimization purposes. */
  string Prestretch_FEMFileName;         /*!< \brief File name for reference geometry. */
  unsigned long Nonphys_Points, /*!< \brief Current number of non-physical points in the solution. */
  Nonphys_Reconstr;           /*!< \brief Current number of non-physical reconstructions for 2nd-order upwinding. */
  bool ParMETIS;              /*!< \brief Boolean for activating ParMETIS mode (while testing). */
  unsigned short DirectDiff;  /*!< \brief Direct Differentation mode. */
  bool DiscreteAdjoint;       /*!< \brief AD-based discrete adjoint mode. */
  su2double *default_vel_inf, /*!< \brief Default freestream velocity array for the COption class. */
  *default_eng_cyl,           /*!< \brief Default engine box array for the COption class. */
  *default_eng_val,           /*!< \brief Default engine box array values for the COption class. */
  *default_cfl_adapt,         /*!< \brief Default CFL adapt param array for the COption class. */
  *default_ad_coeff_flow,     /*!< \brief Default artificial dissipation (flow) array for the COption class. */
  *default_ad_coeff_adj,      /*!< \brief Default artificial dissipation (adjoint) array for the COption class. */
  *default_obj_coeff,         /*!< \brief Default objective array for the COption class. */
  *default_geo_loc,           /*!< \brief Default SU2_GEO section locations array for the COption class. */
  *default_distortion,        /*!< \brief Default SU2_GEO section locations array for the COption class. */
  *default_ea_lim,            /*!< \brief Default equivalent area limit array for the COption class. */
  *default_grid_fix,          /*!< \brief Default fixed grid (non-deforming region) array for the COption class. */
  *default_htp_axis,          /*!< \brief Default HTP axis for the COption class. */
  *default_ffd_axis,          /*!< \brief Default FFD axis for the COption class. */
  *default_inc_crit;          /*!< \brief Default incremental criteria array for the COption class. */

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
  void addNRBCOption(const string name, unsigned short & nMarker_NRBC, string * & Marker_NRBC, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                     su2double* & var1, su2double* & var2, su2double** & FlowDir) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionNRBC<Tenum>(name, nMarker_NRBC, Marker_NRBC, option_field, enum_map, var1, var2, FlowDir);
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

  void addMixingPlaneOption(const string & name, unsigned short & nMarker_MixBound,
                            string* & Marker_MixBound, string* & Marker_MixDonor) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionMixingPlane(name, nMarker_MixBound, Marker_MixBound, Marker_MixDonor);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  template <class Tenum>
  void addTurboPerfOption(const string & name, unsigned short & nMarker_TurboPerf,
                          string* & Marker_TurboBoundIn, string* & Marker_TurboBoundOut,  unsigned short* & Kind_TurboPerformance, const map<string, Tenum> & TurboPerformance_Map) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionTurboPerformance<Tenum>(name, nMarker_TurboPerf, Marker_TurboBoundIn, Marker_TurboBoundOut, Kind_TurboPerformance, TurboPerformance_Map );
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
  CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_nDim, unsigned short verb_level);

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
  SU2_Comm GetMPICommunicator();

  /*!
   * \brief Set the MPI communicator for SU2.
   * \param[in] Communicator - MPI communicator for SU2.
   */
  void SetMPICommunicator(SU2_Comm Communicator);

  /*!
   * \brief Get the number of external iterations.
   * \return Number of external iterations.
   */
  unsigned long GetnExtIter(void);

  /*!
   * \brief If we are prforming an unsteady simulation, there is only
   *        one value of the time step for the complete simulation.
   * \return Value of the time step in an unsteady simulation.
   */
  su2double GetDelta_UnstTime(void);

  void SetMarker_All_BCCustom(unsigned short val_marker, unsigned short val_custom);

};

//CGeometry class
class CGeometry {
protected:
        unsigned long nPoint,	/*!< \brief Number of points of the mesh. */
        nPointDomain,						/*!< \brief Number of real points of the mesh. */
        nPointGhost,					/*!< \brief Number of ghost points of the mesh. */
        nPointNode,					/*!< \brief Size of the node array allocated to hold CPoint objects. */
  Global_nPoint,	/*!< \brief Total number of nodes in a simulation across all processors (including halos). */
        Global_nPointDomain,	/*!< \brief Total number of nodes in a simulation across all processors (excluding halos). */
        nElem,					/*!< \brief Number of elements of the mesh. */
  Global_nElem,	/*!< \brief Total number of elements in a simulation across all processors (all types). */
  Global_nElemDomain,  /*!< \brief Total number of elements in a simulation across all processors (excluding halos). */
        nEdge,					/*!< \brief Number of edges of the mesh. */
        nFace,					/*!< \brief Number of faces of the mesh. */
  nelem_edge,             /*!< \brief Number of edges in the mesh. */
  Global_nelem_edge,      /*!< \brief Total number of edges in the mesh across all processors. */
  nelem_triangle,       /*!< \brief Number of triangles in the mesh. */
  Global_nelem_triangle,       /*!< \brief Total number of triangles in the mesh across all processors. */
  nelem_quad,           /*!< \brief Number of quadrangles in the mesh. */
  Global_nelem_quad,           /*!< \brief Total number of quadrangles in the mesh across all processors. */
  nelem_tetra,          /*!< \brief Number of tetrahedra in the mesh. */
  Global_nelem_tetra,          /*!< \brief Total number of tetrahedra in the mesh across all processors. */
  nelem_hexa,           /*!< \brief Number of hexahedra in the mesh. */
  Global_nelem_hexa,           /*!< \brief Total number of hexahedra in the mesh across all processors. */
  nelem_prism,          /*!< \brief Number of prisms in the mesh. */
  Global_nelem_prism,          /*!< \brief Total number of prisms in the mesh across all processors. */
  nelem_pyramid,        /*!< \brief Number of pyramids in the mesh. */
  Global_nelem_pyramid,        /*!< \brief Total number of pyramids in the mesh across all processors. */
  nelem_edge_bound,           /*!< \brief Number of edges on the mesh boundaries. */
  Global_nelem_edge_bound,           /*!< \brief Total number of edges on the mesh boundaries across all processors. */
  nelem_triangle_bound,          /*!< \brief Number of triangles on the mesh boundaries. */
  Global_nelem_triangle_bound,          /*!< \brief Total number of triangles on the mesh boundaries across all processors. */
  nelem_quad_bound,        /*!< \brief Number of quads on the mesh boundaries. */
  Global_nelem_quad_bound;        /*!< \brief Total number of quads on the mesh boundaries across all processors. */
        unsigned short nDim,	/*!< \brief Number of dimension of the problem. */
        nZone,								/*!< \brief Number of zones in the problem. */
        nMarker;				/*!< \brief Number of different markers of the mesh. */
  unsigned long Max_GlobalPoint;  /*!< \brief Greater global point in the domain local structure. */

public:
        unsigned long *nElem_Bound;			/*!< \brief Number of elements of the boundary. */
        string *Tag_to_Marker;	/*!< \brief If you know the index of the boundary (depend of the
                                                         grid definition), it gives you the maker (where the boundary
                                                         is stored from 0 to boundaries). */
        CPrimalGrid** elem;	/*!< \brief Element vector (primal grid information). */
        CPrimalGrid** face;			/*!< \brief Face vector (primal grid information). */
        CPrimalGrid*** bound;	/*!< \brief Boundary vector (primal grid information). */
        CPoint** node;			/*!< \brief Node vector (dual grid information). */
        CEdge** edge;			/*!< \brief Edge vector (dual grid information). */
        CVertex*** vertex;		/*!< \brief Boundary Vertex vector (dual grid information). */
        unsigned long *nVertex;	/*!< \brief Number of vertex for each marker. */
        unsigned short nCommLevel;		/*!< \brief Number of non-blocking communication levels. */
        vector<unsigned long> PeriodicPoint[MAX_NUMBER_PERIODIC][2];			/*!< \brief PeriodicPoint[Periodic bc] and return the point that
                                                                                                                                                         must be sent [0], and the image point in the periodic bc[1]. */
        vector<unsigned long> PeriodicElem[MAX_NUMBER_PERIODIC];				/*!< \brief PeriodicElem[Periodic bc] and return the elements that
                                                                                                                                                         must be sent. */

  short *Marker_All_SendRecv;

        /*--- Create vectors and distribute the values among the different planes queues ---*/
        vector<vector<su2double> > Xcoord_plane; /*!< \brief Vector containing x coordinates of new points appearing on a single plane */
        vector<vector<su2double> > Ycoord_plane; /*!< \brief Vector containing y coordinates of  new points appearing on a single plane */
        vector<vector<su2double> > Zcoord_plane; 	/*!< \brief Vector containing z coordinates of  new points appearing on a single plane */
        vector<vector<su2double> > FaceArea_plane; /*!< \brief Vector containing area/volume associated with  new points appearing on a single plane */
        vector<vector<unsigned long> > Plane_points; /*!< \brief Vector containing points appearing on a single plane */

        vector<su2double> XCoordList;	/*!< \brief Vector containing points appearing on a single plane */
        CPrimalGrid*** newBound;            /*!< \brief Boundary vector for new periodic elements (primal grid information). */
        unsigned long *nNewElem_Bound;			/*!< \brief Number of new periodic elements of the boundary. */


  /*--- Partitioning-specific variables ---*/
  map<unsigned long,unsigned long> Global_to_Local_Elem;
  unsigned long xadj_size;
  unsigned long adjacency_size;
  unsigned long *starting_node;
  unsigned long *ending_node;
  unsigned long *npoint_procs;
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  idx_t * adjacency;
  idx_t * xadj;
#endif
#endif

        /*!
         * \brief Constructor of the class.
         */
        CGeometry(void);

        /*!
         * \brief Destructor of the class.
         */
        virtual ~CGeometry(void);

        /*!
         * \brief Get number of coordinates.
         * \return Number of coordinates.
         */
        unsigned short GetnDim(void);

        /*!
         * \brief Get number of zones.
         * \return Number of zones.
         */
        unsigned short GetnZone(void);

        /*!
         * \brief Get number of points.
         * \return Number of points.
         */
        unsigned long GetnPoint(void);

        /*!
         * \brief Get number of real points (that belong to the domain).
         * \return Number of real points.
         */
        unsigned long GetnPointDomain(void);
        /*!
         * \brief Get number of elements.
         * \return Number of elements.
         */
        unsigned long GetnElem(void);

        /*!
         * \brief Get number of edges.
         * \return Number of edges.
         */
        unsigned long GetnEdge(void);

        /*!
         * \brief Get number of markers.
         * \return Number of markers.
         */
        unsigned short GetnMarker(void);

        /*!
         * \brief Get number of vertices.
         * \param[in] val_marker - Marker of the boundary.
         * \return Number of vertices.
         */
        unsigned long GetnVertex(unsigned short val_marker);

        void SetImposedVelocity(unsigned short val_marker, unsigned long iVertex, su2double val_vx, su2double val_vy, su2double val_vz);
};
