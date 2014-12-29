/*!
 * \file config_structure.cpp
 * \brief Main file for managing the config file
 * \author F. Palacios, T. Economon, B. Tracey
 * \version 3.2.7 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 Core Developers.
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

#include "../include/config_structure.hpp"

CConfig::CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software, unsigned short val_iZone, unsigned short val_nZone, unsigned short val_nDim, unsigned short verb_level) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Initialize pointers to Null---*/
  
  SetPointersNull();

  /*--- Reading config options  ---*/
  
  SetConfig_Options(val_iZone, val_nZone);

  /*--- Parsing the config file  ---*/
  
  SetParsing(case_filename);

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, val_iZone, val_nDim);

  /*--- Configuration file boundaries/markers setting ---*/
  
  SetMarkers(val_software, val_iZone);

  /*--- Configuration file output ---*/

  if ((rank == MASTER_NODE) && (verb_level == VERB_HIGH) && (val_iZone != 1))
    SetOutput(val_software, val_iZone);

}

CConfig::CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software) {
  
  /*--- Initialize pointers to Null---*/
  
  SetPointersNull();

  /*--- Reading config options  ---*/
  
  SetConfig_Options(1, 1);

  /*--- Parsing the config file  ---*/
  
  SetParsing(case_filename);

  /*--- Configuration file postprocessing ---*/
  
  SetPostprocessing(val_software, 1, 1);
  
}

void CConfig::SetPointersNull(void){

  /*--- Marker Pointers ---*/

  Marker_Euler=NULL;          Marker_FarField=NULL;         Marker_Custom=NULL;
  Marker_SymWall=NULL;        Marker_Pressure=NULL;         Marker_PerBound=NULL;
  Marker_PerDonor=NULL;       Marker_NearFieldBound=NULL;   Marker_InterfaceBound=NULL;
  Marker_Dirichlet=NULL;      Marker_Dirichlet_Elec=NULL;   Marker_Inlet=NULL;
  Marker_Supersonic_Inlet=NULL;    Marker_Outlet=NULL;      Marker_Out_1D=NULL;
  Marker_Isothermal=NULL;     Marker_HeatFlux=NULL;         Marker_EngineInflow=NULL;
  Marker_EngineBleed=NULL;
  Marker_IsothermalCatalytic=NULL; Marker_IsothermalNonCatalytic=NULL;
  Marker_HeatFluxNonCatalytic=NULL; Marker_HeatFluxCatalytic=NULL;
  Marker_EngineExhaust=NULL; Marker_Displacement=NULL;     Marker_Load=NULL;
  Marker_FlowLoad=NULL;       Marker_Neumann=NULL;          Marker_Neumann_Elec=NULL;
  Marker_All_TagBound=NULL;        Marker_CfgFile_TagBound=NULL;       Marker_All_KindBC=NULL;
  Marker_CfgFile_KindBC=NULL;    Marker_All_SendRecv=NULL; Marker_All_PerBound=NULL;

  /*--- Boundary Condition settings ---*/

  Dirichlet_Value=NULL;       Nozzle_Ttotal=NULL;
  Nozzle_Ptotal=NULL;         Inlet_Ttotal=NULL;            Inlet_Ptotal=NULL;
  Inlet_FlowDir=NULL;         Inlet_Temperature=NULL;       Inlet_Pressure=NULL;
  Inlet_Velocity=NULL;        Inflow_Mach_Target=NULL;     Inflow_Mach=NULL;
  Inflow_Pressure=NULL;      Bleed_Temperature_Target=NULL;       Bleed_Temperature=NULL;
  Bleed_MassFlow_Target=NULL; Bleed_MassFlow=NULL;
  Bleed_Pressure=NULL;        Outlet_Pressure=NULL;         Isothermal_Temperature=NULL;
  Heat_Flux=NULL;             Displ_Value=NULL;             Load_Value=NULL;
  FlowLoad_Value=NULL;        Periodic_RotCenter=NULL;      Periodic_RotAngles=NULL;
  Periodic_Translation=NULL;  Periodic_Center=NULL;         Periodic_Rotation=NULL;
  Periodic_Translate=NULL;    Wall_Catalycity=NULL;

  /*--- Miscellaneous/unsorted ---*/

  Aeroelastic_plunge=NULL;    Aeroelastic_pitch=NULL;
  MassFrac_FreeStream=NULL;
  Velocity_FreeStream=NULL;
  RefOriginMoment=NULL;     RefOriginMoment_X=NULL;  RefOriginMoment_Y=NULL;
  RefOriginMoment_Z=NULL;   CFLRamp=NULL;            CFL=NULL;
  PlaneTag=NULL;
  Kappa_Flow=NULL;    Kappa_AdjFlow=NULL;  Kappa_TNE2=NULL;
  Kappa_AdjTNE2=NULL;  Kappa_LinFlow=NULL;
  Section_Location=NULL;
  U_FreeStreamND=NULL;

  /*--- Moving mesh pointers ---*/

  Kind_GridMovement = NULL;
  Motion_Origin_X = NULL;     Motion_Origin_Y = NULL;     Motion_Origin_Z = NULL;
  Translation_Rate_X = NULL;  Translation_Rate_Y = NULL;  Translation_Rate_Z = NULL;
  Rotation_Rate_X = NULL;     Rotation_Rate_Y = NULL;     Rotation_Rate_Z = NULL;
  Pitching_Omega_X = NULL;    Pitching_Omega_Y = NULL;    Pitching_Omega_Z = NULL;
  Pitching_Ampl_X = NULL;     Pitching_Ampl_Y = NULL;     Pitching_Ampl_Z = NULL;
  Pitching_Phase_X = NULL;    Pitching_Phase_Y = NULL;    Pitching_Phase_Z = NULL;
  Plunging_Omega_X = NULL;    Plunging_Omega_Y = NULL;    Plunging_Omega_Z = NULL;
  Plunging_Ampl_X = NULL;     Plunging_Ampl_Y = NULL;     Plunging_Ampl_Z = NULL;
  RefOriginMoment_X = NULL;   RefOriginMoment_Y = NULL;   RefOriginMoment_Z = NULL;
  MoveMotion_Origin = NULL;

  /*--- Variable initialization ---*/
  
  ExtIter = 0;
  IntIter = 0;
  
  /*--- Reacting chemistry, collisions, plasma ---*/

  Reactions=NULL;                 Omega00=NULL;               Omega11=NULL;
  Gas_Composition=NULL;           Enthalpy_Formation=NULL;    Blottner=NULL;
  Species_Ref_Temperature=NULL;   Species_Ref_Viscosity=NULL; nElStates=NULL;
  CharElTemp=NULL;                degen=NULL;
  Molar_Mass=NULL;                Particle_Mass=NULL;
  ArrheniusCoefficient=NULL;    ArrheniusEta=NULL;    ArrheniusTheta=NULL;
  CharVibTemp=NULL;             RotationModes=NULL;   Ref_Temperature=NULL;
  Tcf_a=NULL;    Tcf_b=NULL;    Tcb_a=NULL;    Tcb_b=NULL;
  Diss=NULL;
  
}

void CConfig::SetConfig_Options(unsigned short val_iZone, unsigned short val_nZone) {
  
  double default_vec_3d[3];
  double default_vec_2d[2];
  double default_vec_6d[6];
  
  nZone = val_nZone;
  iZone = val_iZone;


  // This config file is parsed by a number of programs to make it easy to write SU2
  // wrapper scripts (in python, go, etc.) so please do
  // the best you can to follow the established format. It's very hard to parse c++ code
  // and none of us that write the parsers want to write a full c++ interpreter. Please
  // play nice with the existing format so that you don't break the existing scripts.


  /* BEGIN_CONFIG_OPTIONS */

  /* CONFIG_CATEGORY: Problem Definition */
  /*--- Options related to problem definition and partitioning ---*/

  /*!\par REGIME_TYPE
   *  DESCRIPTION: Adjoint type \ingroup Config*/
  addEnumOption("REGIME_TYPE", Kind_Regime, Regime_Map, COMPRESSIBLE);
  
  /* DESCRIPTION: Engine subsonic intake region */
  addBoolOption("DEBUG_MODE",DebugMode, false);
  
  /*!\par PHYSICAL_PROBLEM
   *  DESCRIPTION: Physical governing equations \n
   *  Options: NONE (default),EULER, NAVIER_STOKES, RANS, POISSON_EQUATION, ADJ_EULER, ADJ_NAVIER_STOKES, ADJ_RANS, LIN_EULER,
   LIN_NAVIER_STOKES, TNE2_EULER, TNE2_NAVIER_STOKES, ADJ_TNE2_EULER, ADJ_TNE2_NAVIER_STOKES, WAVE_EQUATION, HEAT_EQUATION ,LINEAR_ELASTICITY,
   FLUID_STRUCTURE_EULER, FLUID_STRUCTURE_NAVIER_STOKES, FLUID_STRUCTURE_RANS, TEMPLATE_SOLVER
   *  \ingroup Config*/
  addEnumOption("PHYSICAL_PROBLEM", Kind_Solver, Solver_Map, NO_SOLVER);
  /*!\par MATH_PROBLEM
   *  DESCRIPTION: Mathematical problem \n  Options: DIRECT, ADJOINT \ingroup Config*/
  addMathProblemOption("MATH_PROBLEM" , Adjoint, false , Linearized, false, Restart_Flow, false);
  /*!\par KIND_TURB_MODEL
   *  DESCRIPTION: Specify turbulence model \n Options: SA, SST, ML, NONE (default) \ingroup Config*/
  addEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, NO_TURB_MODEL);

  /* DESCRIPTION: Specify transition model */
  addEnumOption("KIND_TRANS_MODEL", Kind_Trans_Model, Trans_Model_Map, NO_TRANS_MODEL);

  /* DESCRIPTION: Axisymmetric simulation */
  addBoolOption("AXISYMMETRIC", Axisymmetric, false);
  /* DESCRIPTION: Add the gravity force */
  addBoolOption("GRAVITY_FORCE", GravityForce, false);
  /* DESCRIPTION: Perform a low fidelity simulation */
  addBoolOption("LOW_FIDELITY_SIMULATION", LowFidelitySim, false);
  /*!\par RESTART_SOL
   *  DESCRIPTION: Restart solution from native solution file \n Options: NO, YES \ingroup Config */
  addBoolOption("RESTART_SOL", Restart, false);
  /* DESCRIPTION: System of measurements */
  addEnumOption("SYSTEM_MEASUREMENTS", SystemMeasurements, Measurements_Map, SI);

  /* CONFIG_CATEGORY: FluidModel */
  /* DESCRIPTION: Fluid model */
  addEnumOption("FLUID_MODEL", Kind_FluidModel, FluidModel_Map, STANDARD_AIR);



  /* CONFIG_CATEGORY: Freestream Conditions */
  /*--- Options related to freestream specification ---*/

  /*!\par GAS_CONSTANT
   *  DESCRIPTION: Specific gas constant (287.058 J/kg*K (air), only for compressible flows) \ingroup Config*/
  addDoubleOption("GAS_CONSTANT", Gas_Constant, 287.058);
  /*!\par GAMMA_VALUE
   *  DESCRIPTION: Ratio of specific heats (1.4 (air), only for compressible flows) \ingroup Config*/
  addDoubleOption("GAMMA_VALUE", Gamma, 1.4);


  /*--- Options related to VAN der WAALS MODEL and PENG ROBINSON ---*/

  /* DESCRIPTION: Critical Temperature, default value for AIR */
  addDoubleOption("CRITICAL_TEMPERATURE", Temperature_Critical, 131.00);
  /* DESCRIPTION: Critical Pressure, default value for MDM */
  addDoubleOption("CRITICAL_PRESSURE", Pressure_Critical, 3588550.0);
  /* DESCRIPTION: Critical Density, default value for MDM */
  addDoubleOption("CRITICAL_DENSITY", Density_Critical, 263.0);

  /*--- Options related to VAN der WAALS MODEL and PENG ROBINSON ---*/
  /* DESCRIPTION: Critical Density, default value for MDM */
   addDoubleOption("ACENTRIC_FACTOR", Acentric_Factor, 0.035);

   /*--- Options related to Viscosity Model ---*/
  /* DESCRIPTION: model of the viscosity */
  addEnumOption("VISCOSITY_MODEL", Kind_ViscosityModel, ViscosityModel_Map, SUTHERLAND);

  /*--- Options related to Constant Viscosity Model ---*/

  /* DESCRIPTION: default value for AIR */
  addDoubleOption("MU_CONSTANT", Mu_ConstantND , 1.716E-5);

  /*--- Options related to Sutherland Viscosity Model ---*/

  /* DESCRIPTION: Sutherland Viscosity Ref default value for AIR SI */
  addDoubleOption("MU_REF", Mu_RefND, 1.716E-5);
  /* DESCRIPTION: Sutherland Temperature Ref, default value for AIR SI */
  addDoubleOption("MU_T_REF", Mu_Temperature_RefND, 273.15);
  /* DESCRIPTION: Sutherland constant, default value for AIR SI */
  addDoubleOption("SUTHERLAND_CONSTANT", Mu_SND, 110.4);

  /*--- Options related to Thermal Conductivity Model ---*/

  addEnumOption("CONDUCTIVITY_MODEL", Kind_ConductivityModel, ConductivityModel_Map, CONSTANT_PRANDTL);

 /*--- Options related to Constant Thermal Conductivity Model ---*/

 /* DESCRIPTION: default value for AIR */
  addDoubleOption("KT_CONSTANT", Kt_ConstantND , 0.0257);

  /* DESCRIPTION: Reynolds number (non-dimensional, based on the free-stream values) */
  addDoubleOption("REYNOLDS_NUMBER", Reynolds, 0.0);
  /* DESCRIPTION: Reynolds length (1 m by default) */
  addDoubleOption("REYNOLDS_LENGTH", Length_Reynolds, 1.0);
  /* DESCRIPTION: Laminar Prandtl number (0.72 (air), only for compressible flows) */
  addDoubleOption("PRANDTL_LAM", Prandtl_Lam, 0.72);
  /* DESCRIPTION: Turbulent Prandtl number (0.9 (air), only for compressible flows) */
  addDoubleOption("PRANDTL_TURB", Prandtl_Turb, 0.90);
  /* DESCRIPTION: Value of the Bulk Modulus  */
  addDoubleOption("BULK_MODULUS", Bulk_Modulus, 1.42E5);
  /* DESCRIPTION: Artifical compressibility factor  */
  addDoubleOption("ARTCOMP_FACTOR", ArtComp_Factor, 1.0);
  /*!\par MACH_NUMBER
   *  DESCRIPTION:  Mach number (non-dimensional, based on the free-stream values). 0.0 by default \ingroup Config*/
  addDoubleOption("MACH_NUMBER", Mach, 0.0);
  /* DESCRIPTION: Init option to choose between Reynolds or thermodynamics quantities for initializing the solution */
  addEnumOption("INIT_OPTION", Kind_InitOption, InitOption_Map, REYNOLDS);
  /* DESCRIPTION: Free-stream option to choose between density and temperature for initializing the solution */
  addEnumOption("FREESTREAM_OPTION", Kind_FreeStreamOption, FreeStreamOption_Map, TEMPERATURE_FS);
  /*!\par FREESTREAM_PRESSURE
   * DESCRIPTION: Free-stream pressure (101325.0 N/m^2 by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_PRESSURE", Pressure_FreeStream, 101325.0);
  /*!\par FREESTREAM_DENSITY
   * DESCRIPTION: Free-stream density (1.2886 Kg/m^3 (air), 998.2 Kg/m^3 (water)) \ingroup Config*/
  addDoubleOption("FREESTREAM_DENSITY", Density_FreeStream, -1.0);
  /*!\par FREESTREAM_TEMPERATURE
   * DESCRIPTION: Free-stream temperature (288.15 K by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_TEMPERATURE", Temperature_FreeStream, 288.15);
  /*!\par FREESTREAM_TEMPERATURE_VE
   * DESCRIPTION: Free-stream vibrational-electronic temperature (288.15 K by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_TEMPERATURE_VE", Temperature_ve_FreeStream, 288.15);
  default_vec_3d[0] = 1.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
  /*!\par FREESTREAM_VELOCITY
   * DESCRIPTION: Free-stream velocity (m/s) */
  addDoubleArrayOption("FREESTREAM_VELOCITY", 3, Velocity_FreeStream, default_vec_3d);
  /* DESCRIPTION: Free-stream viscosity (1.853E-5 Ns/m^2 (air), 0.798E-3 Ns/m^2 (water)) */
  addDoubleOption("FREESTREAM_VISCOSITY", Viscosity_FreeStream, -1.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_INTERMITTENCY", Intermittency_FreeStream, 1.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_TURBULENCEINTENSITY", TurbulenceIntensity_FreeStream, 0.05);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_NU_FACTOR", NuFactor_FreeStream, 3.0);
  /* DESCRIPTION:  */
  addDoubleOption("ENGINE_NU_FACTOR", NuFactor_Engine, 30.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_TURB2LAMVISCRATIO", Turb2LamViscRatio_FreeStream, 10.0);
  /* DESCRIPTION: Side-slip angle (degrees, only for compressible flows) */
  addDoubleOption("SIDESLIP_ANGLE", AoS, 0.0);
  /*!\par AOA
   *  DESCRIPTION: Angle of attack (degrees, only for compressible flows) \ingroup Config*/
  addDoubleOption("AOA", AoA, 0.0);
  /* DESCRIPTION: Activate fixed CL mode (specify a CL instead of AoA). */
  addBoolOption("FIXED_CL_MODE", Fixed_CL_Mode, false);
  /* DESCRIPTION: Specify a fixed coefficient of lift instead of AoA (only for compressible flows) */
  addDoubleOption("TARGET_CL", Target_CL, 0.0);
  /* DESCRIPTION: Damping factor for fixed CL mode. */
  addDoubleOption("DAMP_FIXED_CL", Damp_Fixed_CL, 0.1);


  /* CONFIG_CATEGORY: Reference Conditions */
  /*--- Options related to reference values for nondimensionalization ---*/

  Length_Ref = 1.0; //<---- NOTE: this should be given an option or set as a const

  /*!\par REF_ORIGIN_MOMENT_X
   *  DESCRIPTION: X Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_X", nRefOriginMoment_X, RefOriginMoment_X);
  /*!\par REF_ORIGIN_MOMENT_Y
   * DESCRIPTION: Y Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_Y", nRefOriginMoment_Y, RefOriginMoment_Y);
  /*!\par REF_ORIGIN_MOMENT_Z
   * DESCRIPTION: Z Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_Z", nRefOriginMoment_Z, RefOriginMoment_Z);
  /*!\par REF_AREA
   * DESCRIPTION: Reference area for force coefficients (0 implies automatic calculation) \ingroup Config*/
  addDoubleOption("REF_AREA", RefAreaCoeff, 1.0);
  /*!\par REF_LENGTH_MOMENT
   * DESCRIPTION: Reference length for pitching, rolling, and yawing non-dimensional moment \ingroup Config*/
  addDoubleOption("REF_LENGTH_MOMENT", RefLengthMoment, 1.0);
  /*!\par REF_ELEM_LENGTH
   * DESCRIPTION: Reference element length for computing the slope limiter epsilon \ingroup Config*/
  addDoubleOption("REF_ELEM_LENGTH", RefElemLength, 0.1);
  /*!\par REF_SHARP_EDGES
   * DESCRIPTION: Reference coefficient for detecting sharp edges \ingroup Config*/
  addDoubleOption("REF_SHARP_EDGES", RefSharpEdges, 3.0);
	/*!\par REF_PRESSURE
	 *  DESCRIPTION: Reference pressure (1.0 N/m^2 by default, only for compressible flows)  \ingroup Config*/
  addDoubleOption("REF_PRESSURE", Pressure_Ref, 1.0);
	/* !\par REF_TEMPERATURE
   *  DESCRIPTION: Reference temperature (1.0 K by default, only for compressible flows)  \ingroup Config*/
  addDoubleOption("REF_TEMPERATURE", Temperature_Ref, 1.0);
	/* !\par REF_DENSITY
   *  DESCRIPTION: Reference density (1.0 Kg/m^3 by default, only for compressible flows)  \ingroup Config*/
  addDoubleOption("REF_DENSITY", Density_Ref, 1.0);
	/* !\par REF_VELOCITY
   *  DESCRIPTION: Reference velocity (incompressible only)  \ingroup Config*/
  addDoubleOption("REF_VELOCITY", Velocity_Ref, -1.0);
	/* !\par REF_VISCOSITY
   *  DESCRIPTION: Reference viscosity (incompressible only)  \ingroup Config*/
  addDoubleOption("REF_VISCOSITY", Viscosity_Ref, -1.0);


  /* CONFIG_CATEGORY: Boundary Markers */
  /*--- Options related to various boundary markers ---*/

  /*!\par MARKER_PLOTTING
   *  DESCRIPTION: Marker(s) of the surface in the surface flow solution file  \ingroup Config*/
  addStringListOption("MARKER_PLOTTING", nMarker_Plotting, Marker_Plotting);
  /*!\par MARKER_MONITORING
   *  DESCRIPTION: Marker(s) of the surface where evaluate the non-dimensional coefficients \ingroup Config*/
  addStringListOption("MARKER_MONITORING", nMarker_Monitoring, Marker_Monitoring);
  /*!\par MARKER_DESIGNING
   * DESCRIPTION: Marker(s) of the surface where objective function (design problem) will be evaluated \ingroup Config*/
  addStringListOption("MARKER_DESIGNING", nMarker_Designing, Marker_Designing);
  /*!\par GEO_MARKER
   * DESCRIPTION: Marker(s) of the surface where evaluate the geometrical functions \ingroup Config*/
  addStringListOption("GEO_MARKER", nMarker_GeoEval, Marker_GeoEval);
  /*!\par MARKER_EULER
   * DESCRIPTION: Euler wall boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_EULER", nMarker_Euler, Marker_Euler);
  /*!\par MARKER_FAR
   * DESCRIPTION: Far-field boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_FAR", nMarker_FarField, Marker_FarField);
  /*!\par MARKER_SYM
   * DESCRIPTION: Symmetry boundary condition \ingroup Config*/
  addStringListOption("MARKER_SYM", nMarker_SymWall, Marker_SymWall);
  /*!\par MARKER_PRESSURE
   * DESCRIPTION: Symmetry boundary condition \ingroup Config*/
  addStringListOption("MARKER_PRESSURE", nMarker_Pressure, Marker_Pressure);
  /*!\par MARKER_NEARFIELD
   * DESCRIPTION: Near-Field boundary condition \ingroup Config*/
  addStringListOption("MARKER_NEARFIELD", nMarker_NearFieldBound, Marker_NearFieldBound);
  /*!\par MARKER_INTERFACE
   * DESCRIPTION: Zone interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_INTERFACE", nMarker_InterfaceBound, Marker_InterfaceBound);
  /*!\par MARKER_DIRICHLET
   *  DESCRIPTION: Dirichlet boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_DIRICHLET", nMarker_Dirichlet, Marker_Dirichlet);
  /* DESCRIPTION: Neumann boundary marker(s) */
  addStringListOption("MARKER_NEUMANN", nMarker_Neumann, Marker_Neumann);
  /* DESCRIPTION: poisson dirichlet boundary marker(s) */
  addStringDoubleListOption("ELEC_DIRICHLET", nMarker_Dirichlet_Elec, Marker_Dirichlet_Elec, Dirichlet_Value );
  /* DESCRIPTION: poisson neumann boundary marker(s) */
  addStringListOption("ELEC_NEUMANN", nMarker_Neumann_Elec, Marker_Neumann_Elec);
  /* DESCRIPTION: Custom boundary marker(s) */
  addStringListOption("MARKER_CUSTOM", nMarker_Custom, Marker_Custom);
  /* DESCRIPTION: Periodic boundary marker(s) for use with SU2_MSH
   Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
   rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
   rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) */
  addPeriodicOption("MARKER_PERIODIC", nMarker_PerBound, Marker_PerBound, Marker_PerDonor,
                    Periodic_RotCenter, Periodic_RotAngles, Periodic_Translation);

  /*!\par MARKER_ACTDISK
   * DESCRIPTION: Periodic boundary marker(s) for use with SU2_MSH
   Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
   rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
   rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) \ingroup Config*/
  addActuatorDiskOption("MARKER_ACTDISK", nMarker_ActDisk_Inlet, nMarker_ActDisk_Outlet,
                        Marker_ActDisk_Inlet, Marker_ActDisk_Outlet,
                        ActDisk_Origin, ActDisk_RootRadius, ActDisk_TipRadius,
                        ActDisk_CT, ActDisk_Omega);

  /* DESCRIPTION: Inlet boundary type */
  addEnumOption("INLET_TYPE", Kind_Inlet, Inlet_Map, TOTAL_CONDITIONS);

  /*!\par MARKER_INLET
   *  DESCRIPTION: Inlet boundary marker(s) with the following formats,
   Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector.
   Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector. \ingroup Config*/
  addInletOption("MARKER_INLET", nMarker_Inlet, Marker_Inlet, Inlet_Ttotal, Inlet_Ptotal, Inlet_FlowDir);

  /* DESCRIPTION: Riemann boundary marker(s) with the following formats, a unit vector. */
  addRiemannOption("MARKER_RIEMANN", nMarker_Riemann, Marker_Riemann, Kind_Data_Riemann, Riemann_Map, Riemann_Var1, Riemann_Var2, Riemann_FlowDir);
  /*!\par MARKER_SUPERSONIC_INLET
   *  DESCRIPTION: Supersonic inlet boundary marker(s) \n   Format: (inlet marker, temperature, static pressure, velocity_x,   velocity_y, velocity_z, ... ), i.e. primitive variables specified. \ingroup Config*/
  addInletOption("MARKER_SUPERSONIC_INLET", nMarker_Supersonic_Inlet, Marker_Supersonic_Inlet,
                 Inlet_Temperature, Inlet_Pressure, Inlet_Velocity);
  /*!\par MARKER_OUTLET
   *  DESCRIPTION: Outlet boundary marker(s)\n
   Format: ( outlet marker, back pressure (static), ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_OUTLET", nMarker_Outlet, Marker_Outlet, Outlet_Pressure);
  /*!\par MARKER_ISOTHERMAL
   * DESCRIPTION: Isothermal wall boundary marker(s)\n
   * Format: ( isothermal marker, wall temperature (static), ... ) \ingroup Config  */
  addStringDoubleListOption("MARKER_ISOTHERMAL", nMarker_Isothermal, Marker_Isothermal, Isothermal_Temperature);
  /*!\par MARKER_ISOTHERMAL_NONCATALYTIC
   *  DESCRIPTION: Isothermal wall boundary marker(s)\n
   Format: ( isothermal marker, wall temperature (static), ... ) \ingroup Config */
  addStringDoubleListOption("MARKER_ISOTHERMAL_NONCATALYTIC", nMarker_IsothermalNonCatalytic, Marker_IsothermalNonCatalytic, Isothermal_Temperature);
  /*!\par MARKER_ISOTHERMAL_CATALYTIC
   *  DESCRIPTION: Isothermal wall boundary marker(s)
   Format: ( isothermal marker, wall temperature (static), ... ) \ingroup Config */
  addStringDoubleListOption("MARKER_ISOTHERMAL_CATALYTIC", nMarker_IsothermalCatalytic, Marker_IsothermalCatalytic, Isothermal_Temperature);
  /*!\par MARKER_HEATFLUX
   *  DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_HEATFLUX", nMarker_HeatFlux, Marker_HeatFlux, Heat_Flux);
  /* DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) */
  addStringDoubleListOption("MARKER_HEATFLUX_NONCATALYTIC", nMarker_HeatFluxNonCatalytic, Marker_HeatFluxNonCatalytic, Heat_Flux);
  /* DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) */
  addStringDoubleListOption("MARKER_HEATFLUX_CATALYTIC", nMarker_HeatFluxCatalytic, Marker_HeatFluxCatalytic, Heat_Flux);
  /*!\par MARKER_ENGINE_INFLOW
   *  DESCRIPTION: Engine inflow boundary marker(s)
   Format: ( nacelle inflow marker, fan face Mach, ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_ENGINE_INFLOW", nMarker_EngineInflow, Marker_EngineInflow, Inflow_Mach_Target);
  /*!\par MARKER_ENGINE_BLEED
   *  DESCRIPTION: Engine bleed boundary marker(s)
   Format: ( nacelle inflow marker, flux, ... ) \ingroup Config*/
  addBleedOption("MARKER_ENGINE_BLEED", nMarker_EngineBleed, Marker_EngineBleed, Bleed_MassFlow_Target, Bleed_Temperature_Target);
  /* DESCRIPTION: Engine subsonic intake region */
  addBoolOption("SUBSONIC_ENGINE", Engine_Intake, false);
  default_vec_6d[0] = -1E15; default_vec_6d[1] = -1E15; default_vec_6d[2] = -1E15;
  default_vec_6d[3] =  1E15; default_vec_6d[4] =  1E15; default_vec_6d[5] =  1E15;
  /* DESCRIPTION: Coordinates of the box to impose a subsonic nacellle (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax) */
  addDoubleArrayOption("SUBSONIC_ENGINE_BOX", 6, Subsonic_Engine_Box, default_vec_6d);
  /* DESCRIPTION: Engine exhaust boundary marker(s)
   Format: (nacelle exhaust marker, total nozzle temp, total nozzle pressure, ... )*/
  addExhaustOption("MARKER_ENGINE_EXHAUST", nMarker_EngineExhaust, Marker_EngineExhaust, Nozzle_Ttotal, Nozzle_Ptotal);
  /* DESCRIPTION: Displacement boundary marker(s) */
  addStringDoubleListOption("MARKER_NORMAL_DISPL", nMarker_Displacement, Marker_Displacement, Displ_Value);
  /* DESCRIPTION: Load boundary marker(s) */
  addStringDoubleListOption("MARKER_NORMAL_LOAD", nMarker_Load, Marker_Load, Load_Value);
  /* DESCRIPTION: Flow load boundary marker(s) */
  addStringDoubleListOption("MARKER_FLOWLOAD", nMarker_FlowLoad, Marker_FlowLoad, FlowLoad_Value);
  /* DESCRIPTION: Damping factor for engine inlet condition */
  addDoubleOption("DAMP_ENGINE_INFLOW", Damp_Engine_Inflow, 0.75);
  /* DESCRIPTION: Damping factor for engine bleed condition */
  addDoubleOption("DAMP_ENGINE_BLEED", Damp_Engine_Bleed, 0.01);
  /* DESCRIPTION: Outlet boundary marker(s) over which to calculate 1-D flow properties
   Format: ( outlet marker) */
  addStringListOption("MARKER_OUT_1D", nMarker_Out_1D, Marker_Out_1D);


  /* CONFIG_CATEGORY: Time-marching */
  /*--- Options related to time-marching ---*/

  /* DESCRIPTION: Unsteady simulation  */
  addEnumOption("UNSTEADY_SIMULATION", Unsteady_Simulation, Unsteady_Map, STEADY);
  /* DESCRIPTION:  Courant-Friedrichs-Lewy condition of the finest grid */
  addDoubleOption("CFL_NUMBER", CFLFineGrid, 1.25);
  /* DESCRIPTION:  Max time step in local time stepping simulations */
  addDoubleOption("MAX_DELTA_TIME", Max_DeltaTime, 1000000);
  default_vec_3d[0] = 1.0; default_vec_3d[1] = 100.0; default_vec_3d[2] = 1.0;
  /* DESCRIPTION: CFL ramp (factor, number of iterations, CFL limit) */
  addDoubleArrayOption("CFL_RAMP", 3, CFLRamp, default_vec_3d);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the adjoint problem */
  addDoubleOption("CFL_REDUCTION_ADJFLOW", CFLRedCoeff_AdjFlow, 0.8);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
  addDoubleOption("CFL_REDUCTION_TURB", CFLRedCoeff_Turb, 1.0);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the turbulent adjoint problem */
  addDoubleOption("CFL_REDUCTION_ADJTURB", CFLRedCoeff_AdjTurb, 1.0);
  /* DESCRIPTION: Number of total iterations */
  addUnsignedLongOption("EXT_ITER", nExtIter, 999999);
  // these options share nRKStep as their size, which is not a good idea in general
  /* DESCRIPTION: Runge-Kutta alpha coefficients */
  addDoubleListOption("RK_ALPHA_COEFF", nRKStep, RK_Alpha_Step);
  /* DESCRIPTION: Time Step for dual time stepping simulations (s) */
  addDoubleOption("UNST_TIMESTEP", Delta_UnstTime, 0.0);
  /* DESCRIPTION: Total Physical Time for dual time stepping simulations (s) */
  addDoubleOption("UNST_TIME", Total_UnstTime, 1.0);
  /* DESCRIPTION: Unsteady Courant-Friedrichs-Lewy number of the finest grid */
  addDoubleOption("UNST_CFL_NUMBER", Unst_CFL, 0.0);
  /* DESCRIPTION: Number of internal iterations (dual time method) */
  addUnsignedLongOption("UNST_INT_ITER", Unst_nIntIter, 100);
  /* DESCRIPTION: Integer number of periodic time instances for Time Spectral */
  addUnsignedShortOption("TIME_INSTANCES", nTimeInstances, 1);
  /* DESCRIPTION: Iteration number to begin unsteady restarts (dual time method) */
  addLongOption("UNST_RESTART_ITER", Unst_RestartIter, 0);
  /* DESCRIPTION: Starting direct solver iteration for the unsteady adjoint */
  addLongOption("UNST_ADJOINT_ITER", Unst_AdjointIter, 0);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FLOW", Kind_TimeIntScheme_Flow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_TNE2", Kind_TimeIntScheme_TNE2, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJTNE2", Kind_TimeIntScheme_AdjTNE2, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJLEVELSET", Kind_TimeIntScheme_AdjLevelSet, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJFLOW", Kind_TimeIntScheme_AdjFlow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_LIN", Kind_TimeIntScheme_LinFlow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_TURB", Kind_TimeIntScheme_Turb, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJTURB", Kind_TimeIntScheme_AdjTurb, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_WAVE", Kind_TimeIntScheme_Wave, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FEA", Kind_TimeIntScheme_FEA, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_HEAT", Kind_TimeIntScheme_Heat, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_POISSON", Kind_TimeIntScheme_Poisson, Time_Int_Map, EULER_IMPLICIT);

  /* CONFIG_CATEGORY: Linear solver definition */
  /*--- Options related to the linear solvers ---*/

  /* DESCRIPTION: Linear solver for the implicit, mesh deformation, or discrete adjoint systems */
  addEnumOption("LINEAR_SOLVER", Kind_Linear_Solver, Linear_Solver_Map, FGMRES);
  /* DESCRIPTION: Preconditioner for the Krylov linear solvers */
  addEnumOption("LINEAR_SOLVER_PREC", Kind_Linear_Solver_Prec, Linear_Solver_Prec_Map, LU_SGS);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("LINEAR_SOLVER_ERROR", Linear_Solver_Error, 1E-5);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("LINEAR_SOLVER_ITER", Linear_Solver_Iter, 10);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("LINEAR_SOLVER_RESTART_FREQUENCY", Linear_Solver_Restart_Frequency, 10);
  /* DESCRIPTION: Relaxation of the linear solver for the implicit formulation */
  addDoubleOption("LINEAR_SOLVER_RELAX", Linear_Solver_Relax, 1.0);
  /* DESCRIPTION: Roe-Turkel preconditioning for low Mach number flows */
  addBoolOption("ROE_TURKEL_PREC", Low_Mach_Precon, false);
  /* DESCRIPTION: Time Step for dual time stepping simulations (s) */
  addDoubleOption("MIN_ROE_TURKEL_PREC", Min_Beta_RoeTurkel, 0.01);
  /* DESCRIPTION: Time Step for dual time stepping simulations (s) */
  addDoubleOption("MAX_ROE_TURKEL_PREC", Max_Beta_RoeTurkel, 0.2);
  /* DESCRIPTION: Linear solver for the turbulent adjoint systems */
  addEnumOption("ADJTURB_LIN_SOLVER", Kind_AdjTurb_Linear_Solver, Linear_Solver_Map, FGMRES);
  /* DESCRIPTION: Preconditioner for the turbulent adjoint Krylov linear solvers */
  addEnumOption("ADJTURB_LIN_PREC", Kind_AdjTurb_Linear_Prec, Linear_Solver_Prec_Map, LU_SGS);
  /* DESCRIPTION: Minimum error threshold for the turbulent adjoint linear solver for the implicit formulation */
  addDoubleOption("ADJTURB_LIN_ERROR", AdjTurb_Linear_Error, 1E-5);
  /* DESCRIPTION: Maximum number of iterations of the turbulent adjoint linear solver for the implicit formulation */
  addUnsignedShortOption("ADJTURB_LIN_ITER", AdjTurb_Linear_Iter, 10);
  /* DESCRIPTION: Entropy fix factor */
  addDoubleOption("ENTROPY_FIX_COEFF", EntropyFix_Coeff, 0.0);
  
  /* CONFIG_CATEGORY: Convergence*/
  /*--- Options related to convergence ---*/
  
  /* DESCRIPTION: Convergence criteria */
  addEnumOption("CONV_CRITERIA", ConvCriteria, Converge_Crit_Map, RESIDUAL);
  /* DESCRIPTION: Residual reduction (order of magnitude with respect to the initial value) */
  addDoubleOption("RESIDUAL_REDUCTION", OrderMagResidual, 3.0);
  /* DESCRIPTION: Min value of the residual (log10 of the residual) */
  addDoubleOption("RESIDUAL_MINVAL", MinLogResidual, -8.0);
  /* DESCRIPTION: Flow functional for the Residual criteria */
  addEnumOption("RESIDUAL_FUNC_FLOW", Residual_Func_Flow, Residual_Map, RHO_RESIDUAL);
  /* DESCRIPTION: Iteration number to begin convergence monitoring */
  addUnsignedLongOption("STARTCONV_ITER", StartConv_Iter, 5);
  /* DESCRIPTION: Number of elements to apply the criteria */
  addUnsignedShortOption("CAUCHY_ELEMS", Cauchy_Elems, 100);
  /* DESCRIPTION: Epsilon to control the series convergence */
  addDoubleOption("CAUCHY_EPS", Cauchy_Eps, 1E-10);
  /* DESCRIPTION: Flow functional for the Cauchy criteria */
  addEnumOption("CAUCHY_FUNC_FLOW", Cauchy_Func_Flow, Objective_Map, DRAG_COEFFICIENT);
  /* DESCRIPTION: Adjoint functional for the Cauchy criteria */
  addEnumOption("CAUCHY_FUNC_ADJFLOW", Cauchy_Func_AdjFlow, Sens_Map, SENS_GEOMETRY);
  /* DESCRIPTION: Linearized functional for the Cauchy criteria */
  addEnumOption("CAUCHY_FUNC_LIN", Cauchy_Func_LinFlow, Linear_Obj_Map, DELTA_DRAG_COEFFICIENT);

  /* CONFIG_CATEGORY: Multi-grid */
  /*--- Options related to Multi-grid ---*/

  /* DESCRIPTION: Start up iterations using the fine grid only */
  addUnsignedShortOption("START_UP_ITER", nStartUpIter, 0);
  /* DESCRIPTION: Multi-grid Levels */
  addUnsignedShortOption("MGLEVEL", nMultiLevel, 3);
  /* DESCRIPTION: Multi-grid cycle */
  addEnumOption("MGCYCLE", MGCycle, MG_Cycle_Map, V_CYCLE);
  /* DESCRIPTION: Multi-grid pre-smoothing level */
  addUShortListOption("MG_PRE_SMOOTH", nMG_PreSmooth, MG_PreSmooth);
  /* DESCRIPTION: Multi-grid post-smoothing level */
  addUShortListOption("MG_POST_SMOOTH", nMG_PostSmooth, MG_PostSmooth);
  /* DESCRIPTION: Jacobi implicit smoothing of the correction */
  addUShortListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth, MG_CorrecSmooth);
  /* DESCRIPTION: Damping factor for the residual restriction */
  addDoubleOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.75);
  /* DESCRIPTION: Damping factor for the correction prolongation */
  addDoubleOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.75);

  /* CONFIG_CATEGORY: Spatial Discretization */
  /*--- Options related to the spatial discretization ---*/

  /* DESCRIPTION: Numerical method for spatial gradients */
  addEnumOption("NUM_METHOD_GRAD", Kind_Gradient_Method, Gradient_Map, WEIGHTED_LEAST_SQUARES);
  /* DESCRIPTION: Coefficient for the limiter */
  addDoubleOption("LIMITER_COEFF", LimiterCoeff, 0.5);
  /* DESCRIPTION: Freeze the value of the limiter after a number of iterations */
  addUnsignedLongOption("LIMITER_ITER", LimiterIter, 999999);
  /* DESCRIPTION: Coefficient for detecting the limit of the sharp edges */
  addDoubleOption("SHARP_EDGES_COEFF", SharpEdgesCoeff, 3.0);

  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centered_Flow, Kind_Upwind_Flow);
  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_FLOW", SpatialOrder_Flow, SpatialOrder_Map, SECOND_ORDER);
  /*!\par SLOPE_LIMITER_FLOW
   * DESCRIPTION: Slope limiter for the direct solution. \n OPTIONS: VENKATAKRISHNAAN (default), BARTH_JESPERSEN, SHARP_EDGES \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_FLOW", Kind_SlopeLimit_Flow, Limiter_Map, VENKATAKRISHNAN);
  default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
  /* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
  addDoubleArrayOption("AD_COEFF_FLOW", 3, Kappa_Flow, default_vec_3d);

  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_ADJFLOW", Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow, Kind_Upwind_AdjFlow);
  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_ADJFLOW", SpatialOrder_AdjFlow, SpatialOrder_Map, SECOND_ORDER);
  /*!\par SLOPE_LIMITER_ADJFLOW
     * DESCRIPTION: Slope limiter for the adjoint solution. \n OPTIONS: VENKATAKRISHNAAN (default), BARTH_JESPERSEN, SHARP_EDGES \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_ADJFLOW", Kind_SlopeLimit_AdjFlow, Limiter_Map, VENKATAKRISHNAN);
  default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
  /* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
  addDoubleArrayOption("AD_COEFF_ADJFLOW", 3, Kappa_AdjFlow, default_vec_3d);

  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_TURB", SpatialOrder_Turb, SpatialOrder_Map, FIRST_ORDER);
  /* DESCRIPTION: Slope limiter */
  addEnumOption("SLOPE_LIMITER_TURB", Kind_SlopeLimit_Turb, Limiter_Map, VENKATAKRISHNAN);
  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_TURB", Kind_ConvNumScheme_Turb, Kind_Centered_Turb, Kind_Upwind_Turb);
  
  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_ADJTURB", SpatialOrder_AdjTurb, SpatialOrder_Map, FIRST_ORDER);
  /* DESCRIPTION: Slope limiter */
  addEnumOption("SLOPE_LIMITER_ADJTURB", Kind_SlopeLimit_AdjTurb, Limiter_Map, VENKATAKRISHNAN);
  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_ADJTURB", Kind_ConvNumScheme_AdjTurb, Kind_Centered_AdjTurb, Kind_Upwind_AdjTurb);

  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_LIN", Kind_ConvNumScheme_LinFlow, Kind_Centered_LinFlow, Kind_Upwind_LinFlow);
  default_vec_2d[0] = 0.15; default_vec_2d[1] = 0.02;
  /* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
  addDoubleArrayOption("AD_COEFF_LIN", 2, Kappa_LinFlow, default_vec_2d);

  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_ADJLEVELSET", SpatialOrder_AdjLevelSet, SpatialOrder_Map, SECOND_ORDER);
  /* DESCRIPTION: Slope limiter */
  addEnumOption("SLOPE_LIMITER_ADJLEVELSET", Kind_SlopeLimit_AdjLevelSet, Limiter_Map, VENKATAKRISHNAN);
  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_ADJLEVELSET", Kind_ConvNumScheme_AdjLevelSet, Kind_Centered_AdjLevelSet, Kind_Upwind_AdjLevelSet);

  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_TNE2", Kind_ConvNumScheme_TNE2, Kind_Centered_TNE2, Kind_Upwind_TNE2);
  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_TNE2", SpatialOrder_TNE2, SpatialOrder_Map, SECOND_ORDER);
  /* DESCRIPTION: Slope limiter */
  addEnumOption("SLOPE_LIMITER_TNE2", Kind_SlopeLimit_TNE2, Limiter_Map, VENKATAKRISHNAN);
  default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
  /* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
  addDoubleArrayOption("AD_COEFF_TNE2", 3, Kappa_TNE2, default_vec_3d);

  /* DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_ADJTNE2", Kind_ConvNumScheme_AdjTNE2, Kind_Centered_AdjTNE2, Kind_Upwind_AdjTNE2);
  /* DESCRIPTION: Spatial numerical order integration */
  addEnumOption("SPATIAL_ORDER_ADJTNE2", SpatialOrder_AdjTNE2, SpatialOrder_Map, SECOND_ORDER);
  /* DESCRIPTION: Slope limiter */
  addEnumOption("SLOPE_LIMITER_ADJTNE2", Kind_SlopeLimit_AdjTNE2, Limiter_Map, VENKATAKRISHNAN);
  default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
  /* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
  addDoubleArrayOption("AD_COEFF_ADJTNE2", 3, Kappa_AdjTNE2, default_vec_3d);

  /* CONFIG_CATEGORY: Adjoint and Gradient */
  /*--- Options related to the adjoint and gradient ---*/

  /* DESCRIPTION: Limit value for the adjoint variable */
  addDoubleOption("LIMIT_ADJFLOW", AdjointLimit, 1E6);
  /* DESCRIPTION: Multigrid with the adjoint problem */
  addBoolOption("MG_ADJFLOW", MG_AdjointFlow, true);
  /* DESCRIPTION: Adjoint problem boundary condition */
  addEnumOption("OBJECTIVE_FUNCTION", Kind_ObjFunc, Objective_Map, DRAG_COEFFICIENT);
  default_vec_2d[0] = 0.0; default_vec_2d[1] = 1.0;
  /* DESCRIPTION: Definition of the airfoil section */
  addDoubleArrayOption("GEO_LOCATION_SECTIONS", 2, Section_Location, default_vec_2d);
  /* DESCRIPTION: Identify the axis of the section */
  addEnumOption("GEO_ORIENTATION_SECTIONS", Axis_Orientation, Axis_Orientation_Map, Y_AXIS);
  /* DESCRIPTION: Percentage of new elements (% of the original number of elements) */
  addUnsignedShortOption("GEO_NUMBER_SECTIONS", nSections, 5);
  /* DESCRIPTION: Number of section cuts to make when calculating internal volume */
  addUnsignedShortOption("GEO_VOLUME_SECTIONS", nVolSections, 101);
  /* DESCRIPTION: Output sectional forces for specified markers. */
  addBoolOption("GEO_PLOT_SECTIONS", Plot_Section_Forces, false);
  /* DESCRIPTION: Mode of the GDC code (analysis, or gradient) */
  addEnumOption("GEO_MODE", GeometryMode, GeometryMode_Map, FUNCTION);

  /* DESCRIPTION: Drag weight in sonic boom Objective Function (from 0.0 to 1.0) */
  addDoubleOption("DRAG_IN_SONICBOOM", WeightCd, 0.0);
  /* DESCRIPTION: Sensitivity smoothing  */
  addEnumOption("SENS_SMOOTHING", Kind_SensSmooth, Sens_Smoothing_Map, NO_SMOOTH);
  /* DESCRIPTION: Continuous governing equation set  */
  addEnumOption("CONTINUOUS_EQNS", Continuous_Eqns, ContinuousEqns_Map, EULER_EQNS);
  /* DESCRIPTION: Discrete governing equation set */
  addEnumOption("DISCRETE_EQNS", Discrete_Eqns, DiscreteEqns_Map, NONE_EQNS);
  /* DESCRIPTION: Adjoint frozen viscosity */
  addBoolOption("FROZEN_VISC", Frozen_Visc, true);
   /* DESCRIPTION:  */
  addDoubleOption("FIX_AZIMUTHAL_LINE", FixAzimuthalLine, 90.0);
  /*!\par SENS_REMOVE_SHARP
   * DESCRIPTION: Remove sharp edges from the sensitivity evaluation  \n Format: SENS_REMOVE_SHARP = YES \n Default: NO \ingroup Config*/
  addBoolOption("SENS_REMOVE_SHARP", Sens_Remove_Sharp, false);
  /* DESCRIPTION: P-norm for heat-flux based objective functions. */
	addDoubleOption("PNORM_HEATFLUX", pnorm_heat, 1.0);

  /* CONFIG_CATEGORY: Input/output files and formats */
  /*--- Options related to input/output files and formats ---*/

  /* DESCRIPTION: I/O */
  addEnumOption("OUTPUT_FORMAT", Output_FileFormat, Output_Map, TECPLOT);
  /* DESCRIPTION: Mesh input file format */
  addEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, SU2);
  /* DESCRIPTION: Convert a CGNS mesh to SU2 format */
  addBoolOption("CGNS_TO_SU2", CGNS_To_SU2, false);
  /* DESCRIPTION:  Mesh input file */
  addStringOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));

  /* DESCRIPTION: Factor for scaling the mesh */
  addDoubleOption("MESH_SCALE_CHANGE", Mesh_Scale_Change, 1.0);
  /* DESCRIPTION: Write a new mesh converted to meters */
  addBoolOption("MESH_OUTPUT", Mesh_Output, false);
  /* DESCRIPTION: Mesh output file */
  addStringOption("MESH_OUT_FILENAME", Mesh_Out_FileName, string("mesh_out.su2"));

  /* DESCRIPTION: Output file convergence history (w/o extension) */
  addStringOption("CONV_FILENAME", Conv_FileName, string("history"));
  /* DESCRIPTION: Output file forces breakdown */
  addStringOption("BREAKDOWN_FILENAME", Breakdown_FileName, string("forces_breakdown.dat"));
  /* DESCRIPTION: Restart flow input file */
  addStringOption("SOLUTION_FLOW_FILENAME", Solution_FlowFileName, string("solution_flow.dat"));
  /* DESCRIPTION: Restart linear flow input file */
  addStringOption("SOLUTION_LIN_FILENAME", Solution_LinFileName, string("solution_lin.dat"));
  /* DESCRIPTION: Restart adjoint input file */
  addStringOption("SOLUTION_ADJ_FILENAME", Solution_AdjFileName, string("solution_adj.dat"));
  /* DESCRIPTION: Output file restart flow */
  addStringOption("RESTART_FLOW_FILENAME", Restart_FlowFileName, string("restart_flow.dat"));
  /* DESCRIPTION: Output file linear flow */
  addStringOption("RESTART_LIN_FILENAME",Restart_LinFileName, string("restart_lin.dat"));
  /* DESCRIPTION: Output file restart adjoint */
  addStringOption("RESTART_ADJ_FILENAME", Restart_AdjFileName, string("restart_adj.dat"));
  /* DESCRIPTION: Output file restart wave */
  addStringOption("RESTART_WAVE_FILENAME", Restart_WaveFileName, string("restart_wave.dat"));
  /* DESCRIPTION: Output file flow (w/o extension) variables */
  addStringOption("VOLUME_FLOW_FILENAME", Flow_FileName, string("flow"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  addStringOption("VOLUME_STRUCTURE_FILENAME", Structure_FileName, string("structure"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  addStringOption("SURFACE_STRUCTURE_FILENAME", SurfStructure_FileName, string("surface_structure"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  addStringOption("SURFACE_WAVE_FILENAME", SurfWave_FileName, string("surface_wave"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  addStringOption("SURFACE_HEAT_FILENAME", SurfHeat_FileName, string("surface_heat"));
  /* DESCRIPTION: Output file wave (w/o extension) variables */
  addStringOption("VOLUME_WAVE_FILENAME", Wave_FileName, string("wave"));
  /* DESCRIPTION: Output file wave (w/o extension) variables */
  addStringOption("VOLUME_HEAT_FILENAME", Heat_FileName, string("heat"));
  /* DESCRIPTION: Output file adj. wave (w/o extension) variables */
  addStringOption("VOLUME_ADJWAVE_FILENAME", AdjWave_FileName, string("adjoint_wave"));
  /* DESCRIPTION: Output file adjoint (w/o extension) variables */
  addStringOption("VOLUME_ADJ_FILENAME", Adj_FileName, string("adjoint"));
  /* DESCRIPTION: Output file linear (w/o extension) variables */
  addStringOption("VOLUME_LIN_FILENAME", Lin_FileName, string("linearized"));
  /* DESCRIPTION: Output objective function gradient */
  addStringOption("GRAD_OBJFUNC_FILENAME", ObjFunc_Grad_FileName, string("of_grad.dat"));
  /* DESCRIPTION: Output objective function */
  addStringOption("VALUE_OBJFUNC_FILENAME", ObjFunc_Value_FileName, string("of_func.dat"));
  /* DESCRIPTION: Output file surface flow coefficient (w/o extension) */
  addStringOption("SURFACE_FLOW_FILENAME", SurfFlowCoeff_FileName, string("surface_flow"));
  /* DESCRIPTION: Output file surface adjoint coefficient (w/o extension) */
  addStringOption("SURFACE_ADJ_FILENAME", SurfAdjCoeff_FileName, string("surface_adjoint"));
  /* DESCRIPTION: Output file surface linear coefficient (w/o extension) */
  addStringOption("SURFACE_LIN_FILENAME", SurfLinCoeff_FileName, string("surface_linear"));
  /* DESCRIPTION: Writing solution file frequency */
  addUnsignedLongOption("WRT_SOL_FREQ", Wrt_Sol_Freq, 1000);
  /* DESCRIPTION: Writing solution file frequency */
  addUnsignedLongOption("WRT_SOL_FREQ_DUALTIME", Wrt_Sol_Freq_DualTime, 1);
  /* DESCRIPTION: Writing convergence history frequency */
  addUnsignedLongOption("WRT_CON_FREQ",  Wrt_Con_Freq, 1);
  /* DESCRIPTION: Writing convergence history frequency for the dual time */
  addUnsignedLongOption("WRT_CON_FREQ_DUALTIME",  Wrt_Con_Freq_DualTime, 10);
  /* DESCRIPTION: Write a volume solution file */
  addBoolOption("LOW_MEMORY_OUTPUT", Low_MemoryOutput, false);
  /* DESCRIPTION: Write a volume solution file */
  addBoolOption("WRT_VOL_SOL", Wrt_Vol_Sol, true);
  /* DESCRIPTION: Write a surface solution file */
  addBoolOption("WRT_SRF_SOL", Wrt_Srf_Sol, true);
  /* DESCRIPTION: Write a surface CSV solution file */
  addBoolOption("WRT_CSV_SOL", Wrt_Csv_Sol, true);
  /* DESCRIPTION: Output residual info to solution/restart file */
  addBoolOption("WRT_RESIDUALS", Wrt_Residuals, false);
  /* DESCRIPTION: Output residual info to solution/restart file */
  addBoolOption("WRT_LIMITERS", Wrt_Limiters, false);
  /* DESCRIPTION: Output residual info to solution/restart file */
  addBoolOption("WRT_SHARPEDGES", Wrt_SharpEdges, false);
  /* DESCRIPTION: Output the rind layers in the solution files */
  addBoolOption("WRT_HALO", Wrt_Halo, false);
  /* DESCRIPTION: Output averaged stagnation pressure on specified exit marker. */
  addBoolOption("ONE_D_OUTPUT", Wrt_1D_Output, false);

  /* CONFIG_CATEGORY: Dynamic mesh definition */
  /*--- Options related to dynamic meshes ---*/

  /* DESCRIPTION: Mesh motion for unsteady simulations */
  addBoolOption("GRID_MOVEMENT", Grid_Movement, false);
  /* DESCRIPTION: Type of mesh motion */
  addEnumListOption("GRID_MOVEMENT_KIND", nGridMovement, Kind_GridMovement, GridMovement_Map);
  /* DESCRIPTION: Marker(s) of moving surfaces (MOVING_WALL or DEFORMING grid motion). */
  addStringListOption("MARKER_MOVING", nMarker_Moving, Marker_Moving);
  /* DESCRIPTION: Mach number (non-dimensional, based on the mesh velocity and freestream vals.) */
  addDoubleOption("MACH_MOTION", Mach_Motion, 0.0);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  addDoubleListOption("MOTION_ORIGIN_X", nMotion_Origin_X, Motion_Origin_X);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  addDoubleListOption("MOTION_ORIGIN_Y", nMotion_Origin_Y, Motion_Origin_Y);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  addDoubleListOption("MOTION_ORIGIN_Z", nMotion_Origin_Z, Motion_Origin_Z);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("TRANSLATION_RATE_X", nTranslation_Rate_X, Translation_Rate_X);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("TRANSLATION_RATE_Y", nTranslation_Rate_Y, Translation_Rate_Y);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("TRANSLATION_RATE_Z", nTranslation_Rate_Z, Translation_Rate_Z);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("ROTATION_RATE_X", nRotation_Rate_X, Rotation_Rate_X);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("ROTATION_RATE_Y", nRotation_Rate_Y, Rotation_Rate_Y);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("ROTATION_RATE_Z", nRotation_Rate_Z, Rotation_Rate_Z);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_OMEGA_X", nPitching_Omega_X, Pitching_Omega_X);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_OMEGA_Y", nPitching_Omega_Y, Pitching_Omega_Y);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_OMEGA_Z", nPitching_Omega_Z, Pitching_Omega_Z);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_AMPL_X", nPitching_Ampl_X, Pitching_Ampl_X);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_AMPL_Y", nPitching_Ampl_Y, Pitching_Ampl_Y);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_AMPL_Z", nPitching_Ampl_Z, Pitching_Ampl_Z);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_PHASE_X", nPitching_Phase_X, Pitching_Phase_X);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_PHASE_Y", nPitching_Phase_Y, Pitching_Phase_Y);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleListOption("PITCHING_PHASE_Z", nPitching_Phase_Z, Pitching_Phase_Z);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("PLUNGING_OMEGA_X", nPlunging_Omega_X, Plunging_Omega_X);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("PLUNGING_OMEGA_Y", nPlunging_Omega_Y, Plunging_Omega_Y);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("PLUNGING_OMEGA_Z", nPlunging_Omega_Z, Plunging_Omega_Z);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("PLUNGING_AMPL_X", nPlunging_Ampl_X, Plunging_Ampl_X);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("PLUNGING_AMPL_Y", nPlunging_Ampl_Y, Plunging_Ampl_Y);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleListOption("PLUNGING_AMPL_Z", nPlunging_Ampl_Z, Plunging_Ampl_Z);
  /* DESCRIPTION: Value to move motion origins (1 or 0) */
  addUShortListOption("MOVE_MOTION_ORIGIN", nMoveMotion_Origin, MoveMotion_Origin);
  /* DESCRIPTION:  */
  addStringOption("MOTION_FILENAME", Motion_Filename, string("mesh_motion.dat"));
  /* DESCRIPTION: Uncoupled Aeroelastic Frequency Plunge. */
  addDoubleOption("FREQ_PLUNGE_AEROELASTIC", FreqPlungeAeroelastic, 100);
  /* DESCRIPTION: Uncoupled Aeroelastic Frequency Pitch. */
  addDoubleOption("FREQ_PITCH_AEROELASTIC", FreqPitchAeroelastic, 100);

  /* CONFIG_CATEGORY: Grid adaptation */
  /*--- Options related to grid adaptation ---*/

  /* DESCRIPTION: Kind of grid adaptation */
  addEnumOption("KIND_ADAPT", Kind_Adaptation, Adapt_Map, NO_ADAPT);
  /* DESCRIPTION: Percentage of new elements (% of the original number of elements) */
  addDoubleOption("NEW_ELEMS", New_Elem_Adapt, -1.0);
  /* DESCRIPTION: Scale factor for the dual volume */
  addDoubleOption("DUALVOL_POWER", DualVol_Power, 0.5);
  /* DESCRIPTION: Use analytical definition for surfaces */
  addEnumOption("ANALYTICAL_SURFDEF", Analytical_Surface, Geo_Analytic_Map, NO_GEO_ANALYTIC);
  /* DESCRIPTION: Before each computation, implicitly smooth the nodal coordinates */
  addBoolOption("SMOOTH_GEOMETRY", SmoothNumGrid, false);
  /* DESCRIPTION: Adapt the boundary elements */
  addBoolOption("ADAPT_BOUNDARY", AdaptBoundary, true);
  /* DESCRIPTION: Divide rectangles into triangles */
  addBoolOption("DIVIDE_ELEMENTS", Divide_Element, false);

  /* CONFIG_CATEGORY: Wind Gust */
  /*--- Options related to wind gust simulations ---*/

  /* DESCRIPTION: Apply a wind gust */
  addBoolOption("WIND_GUST", Wind_Gust, false);
  /* DESCRIPTION: Type of gust */
  addEnumOption("GUST_TYPE", Gust_Type, Gust_Type_Map, NO_GUST);
  /* DESCRIPTION: Gust wavelenght (meters) */
  addDoubleOption("GUST_WAVELENGTH", Gust_WaveLength, 0.0);
  /* DESCRIPTION: Number of gust periods */
  addDoubleOption("GUST_PERIODS", Gust_Periods, 1.0);
  /* DESCRIPTION: Gust amplitude (m/s) */
  addDoubleOption("GUST_AMPL", Gust_Ampl, 0.0);
  /* DESCRIPTION: Time at which to begin the gust (sec) */
  addDoubleOption("GUST_BEGIN_TIME", Gust_Begin_Time, 0.0);
  /* DESCRIPTION: Location at which the gust begins (meters) */
  addDoubleOption("GUST_BEGIN_LOC", Gust_Begin_Loc, 0.0);
  /* DESCRIPTION: Direction of the gust X or Y dir */
  addEnumOption("GUST_DIR", Gust_Dir, Gust_Dir_Map, Y_DIR);


  /* CONFIG_CATEGORY: Equivalent Area */
  /*--- Options related to the equivalent area ---*/

  /* DESCRIPTION: Evaluate equivalent area on the Near-Field  */
  addBoolOption("EQUIV_AREA", EquivArea, false);
  default_vec_3d[0] = 0.0; default_vec_3d[1] = 1.0; default_vec_3d[2] = 1.0;
  /* DESCRIPTION: Integration limits of the equivalent area ( xmin, xmax, Dist_NearField ) */
  addDoubleArrayOption("EA_INT_LIMIT", 3, EA_IntLimit, default_vec_3d);
  /* DESCRIPTION: Equivalent area scaling factor */
  addDoubleOption("EA_SCALE_FACTOR", EA_ScaleFactor, 1E-7);

	/* CONFIG_CATEGORY: Reacting Flow */
  /*--- Options related to the reacting gas mixtures ---*/

	/* DESCRIPTION: Specify chemical model for multi-species simulations */
	addEnumOption("GAS_MODEL", Kind_GasModel, GasModel_Map, ARGON);
	/* DESCRIPTION:  */
	addDoubleListOption("GAS_COMPOSITION", nTemp, Gas_Composition);

	/* CONFIG_CATEGORY: Free surface simulation */
	/*--- Options related to free surface simulation ---*/

	/* DESCRIPTION: Ratio of density for two phase problems */
  addDoubleOption("RATIO_DENSITY", RatioDensity, 0.1);
	/* DESCRIPTION: Ratio of viscosity for two phase problems */
  addDoubleOption("RATIO_VISCOSITY", RatioViscosity, 0.1);
	/* DESCRIPTION: Location of the freesurface (y or z coordinate) */
  addDoubleOption("FREESURFACE_ZERO", FreeSurface_Zero, 0.0);
	/* DESCRIPTION: Free surface depth surface (x or y coordinate) */
  addDoubleOption("FREESURFACE_DEPTH", FreeSurface_Depth, 1.0);
	/* DESCRIPTION: Thickness of the interface in a free surface problem */
  addDoubleOption("FREESURFACE_THICKNESS", FreeSurface_Thickness, 0.1);
	/* DESCRIPTION: Free surface damping coefficient */
  addDoubleOption("FREESURFACE_DAMPING_COEFF", FreeSurface_Damping_Coeff, 0.0);
	/* DESCRIPTION: Free surface damping length (times the baseline wave) */
  addDoubleOption("FREESURFACE_DAMPING_LENGTH", FreeSurface_Damping_Length, 1.0);
	/* DESCRIPTION: Location of the free surface outlet surface (x or y coordinate) */
  addDoubleOption("FREESURFACE_OUTLET", FreeSurface_Outlet, 0.0);

	// these options share nDV as their size in the option references; not a good idea
	/* CONFIG_CATEGORY: Grid deformation */
  /*--- Options related to the grid deformation ---*/

	/* DESCRIPTION: Kind of deformation */
	addEnumListOption("DV_KIND", nDV, Design_Variable, Param_Map);
	/* DESCRIPTION: Marker of the surface to which we are going apply the shape deformation */
  addStringListOption("DV_MARKER", nMarker_DV, Marker_DV);
	/* DESCRIPTION: New value of the shape deformation */
	addDoubleListOption("DV_VALUE", nDV, DV_Value);
	/* DESCRIPTION: Parameters of the shape deformation
   - FFD_CONTROL_POINT_2D ( FFDBox ID, i_Ind, j_Ind, x_Disp, y_Disp )
   - FFD_RADIUS_2D ( FFDBox ID )
   - FFD_CAMBER_2D ( FFDBox ID, i_Ind )
   - FFD_THICKNESS_2D ( FFDBox ID, i_Ind )
   - HICKS_HENNE ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc )
   - COSINE_BUMP ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc, Thickness )
   - FOURIER ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), index, cos(0)/sin(1) )
   - NACA_4DIGITS ( 1st digit, 2nd digit, 3rd and 4th digit )
   - PARABOLIC ( Center, Thickness )
   - DISPLACEMENT ( x_Disp, y_Disp, z_Disp )
   - ROTATION ( x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - OBSTACLE ( Center, Bump size )
   - SPHERICAL ( ControlPoint_Index, Theta_Disp, R_Disp )
   - FFD_CONTROL_POINT ( FFDBox ID, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )
   - FFD_DIHEDRAL_ANGLE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_TWIST_ANGLE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_ROTATION ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CONTROL_SURFACE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CAMBER ( FFDBox ID, i_Ind, j_Ind )
   - FFD_THICKNESS ( FFDBox ID, i_Ind, j_Ind ) */
	addDVParamOption("DV_PARAM", nDV, ParamDV, FFDTag, Design_Variable);
	/* DESCRIPTION: Hold the grid fixed in a region */
  addBoolOption("HOLD_GRID_FIXED", Hold_GridFixed, false);
	default_vec_6d[0] = -1E15; default_vec_6d[1] = -1E15; default_vec_6d[2] = -1E15;
	default_vec_6d[3] =  1E15; default_vec_6d[4] =  1E15; default_vec_6d[5] =  1E15;
	/* DESCRIPTION: Coordinates of the box where the grid will be deformed (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax) */
	addDoubleArrayOption("HOLD_GRID_FIXED_COORD", 6, Hold_GridFixed_Coord, default_vec_6d);
	/* DESCRIPTION: Visualize the deformation */
  addBoolOption("VISUALIZE_DEFORMATION", Visualize_Deformation, false);
  /* DESCRIPTION: Print the residuals during mesh deformation to the console */
  addBoolOption("DEFORM_CONSOLE_OUTPUT", Deform_Output, false);
  /* DESCRIPTION: Number of nonlinear deformation iterations (surface deformation increments) */
  addUnsignedLongOption("DEFORM_NONLINEAR_ITER", GridDef_Nonlinear_Iter, 5);
  /* DESCRIPTION: Number of smoothing iterations for FEA mesh deformation */
  addUnsignedLongOption("DEFORM_LINEAR_ITER", GridDef_Linear_Iter, 500);
  /* DESCRIPTION: Factor to multiply smallest volume for deform tolerance (0.001 default) */
  addDoubleOption("DEFORM_TOL_FACTOR", Deform_Tol_Factor, 0.001);
  /* DESCRIPTION: Type of element stiffness imposed for FEA mesh deformation (INVERSE_VOLUME, WALL_DISTANCE, CONSTANT_STIFFNESS) */
  addEnumOption("DEFORM_STIFFNESS_TYPE", Deform_Stiffness_Type, Deform_Stiffness_Map, INVERSE_VOLUME);
  /* DESCRIPTION: Poisson's ratio for constant stiffness FEA method of grid deformation*/
  addDoubleOption("DEFORM_ELASTICITY_MODULUS", Deform_ElasticityMod, 2E11);
  /* DESCRIPTION: Young's modulus and Poisson's ratio for constant stiffness FEA method of grid deformation*/
  addDoubleOption("DEFORM_POISSONS_RATIO", Deform_PoissonRatio, 0.3);

  /* CONFIG_CATEGORY: Rotorcraft problem */
  /*--- option related to rotorcraft problems ---*/

  /* DESCRIPTION: MISSING --- */
  addDoubleOption("CYCLIC_PITCH", Cyclic_Pitch, 0.0);
  /* DESCRIPTION: MISSING --- */
  addDoubleOption("COLLECTIVE_PITCH", Collective_Pitch, 0.0);


  /* CONFIG_CATEGORY: FEA solver */
  /*--- Options related to the FEA solver ---*/

  /* DESCRIPTION: Modulus of elasticity */
  addDoubleOption("ELASTICITY_MODULUS", ElasticyMod, 2E11);
  /* DESCRIPTION: Poisson ratio */
  addDoubleOption("POISSON_RATIO", PoissonRatio, 0.30);
  /* DESCRIPTION: Material density */
  addDoubleOption("MATERIAL_DENSITY", MaterialDensity, 7854);

  /* CONFIG_CATEGORY: Wave solver */
  /*--- options related to the wave solver ---*/

  /* DESCRIPTION: Constant wave speed */
  addDoubleOption("WAVE_SPEED", Wave_Speed, 331.79);

  /* CONFIG_CATEGORY: Heat solver */
  /*--- options related to the heat solver ---*/

  /* DESCRIPTION: Thermal diffusivity constant */
  addDoubleOption("THERMAL_DIFFUSIVITY", Thermal_Diffusivity, 1.172E-5);

  /* CONFIG_CATEGORY: Visualize Control Volumes */
  /*--- options related to visualizing control volumes ---*/

  /* DESCRIPTION: Node number for the CV to be visualized */
  addLongOption("VISUALIZE_CV", Visualize_CV, -1);

  /* CONFIG_CATEGORY: Inverse design problem */
  /*--- options related to inverse design problem ---*/

  /* DESCRIPTION: Evaluate inverse design on the surface  */
  addBoolOption("INV_DESIGN_CP", InvDesign_Cp, false);

  /* DESCRIPTION: Evaluate inverse design on the surface  */
  addBoolOption("INV_DESIGN_HEATFLUX", InvDesign_HeatFlux, false);

  /* CONFIG_CATEGORY: Unsupported options */
  /*--- Options that are experimental and not intended for general use ---*/

  /* DESCRIPTION: Write extra output */
  addBoolOption("EXTRA_OUTPUT", ExtraOutput, false);

  /* DESCRIPTION: Location of the turb model itself */
  addStringOption("ML_TURB_MODEL_FILE", ML_Turb_Model_File, string("model.json"));

  /* DESCRIPTION: what kind of input/output feature map is there */
  addStringOption("ML_TURB_MODEL_FEATURESET", ML_Turb_Model_FeatureSet, string("none"));

  /* DESCRIPTION: Extra values for ML Turb model */
  addStringListOption("ML_TURB_MODEL_EXTRA",nML_Turb_Model_Extra, ML_Turb_Model_Extra);

  /*--- options related to the FFD problem ---*/
  /* CONFIG_CATEGORY:FFD point inversion */

  /* DESCRIPTION: Number of total iterations in the FFD point inversion */
  addUnsignedShortOption("FFD_ITERATIONS", nFFD_Iter, 500);

  /* DESCRIPTION: Free surface damping coefficient */
	addDoubleOption("FFD_TOLERANCE", FFD_Tol, 1E-8);

  /*--- options that are used in the python optimization scripts. These have no effect on the c++ toolsuite ---*/
  /* CONFIG_CATEGORY:Python Options*/

  /* DESCRIPTION: Gradient method */
  addPythonOption("GRADIENT_METHOD");

  /* DESCRIPTION: Geometrical Parameter */
  addPythonOption("GEO_PARAM");

  /* DESCRIPTION: Setup for design variables */
  addPythonOption("DEFINITION_DV");

  /* DESCRIPTION: Maximum number of iterations */
  addPythonOption("OPT_ITERATIONS");
  
  /* DESCRIPTION: Requested accuracy */
  addPythonOption("OPT_ACCURACY");
  
  /* DESCRIPTION: Setup for design variables */
  addPythonOption("BOUND_DV");
  
  /* DESCRIPTION: Current value of the design variables */
  addPythonOption("DV_VALUE_NEW");

  /* DESCRIPTION: Previous value of the design variables */
  addPythonOption("DV_VALUE_OLD");

  /* DESCRIPTION: Number of partitions of the mesh */
  addPythonOption("NUMBER_PART");

  /* DESCRIPTION: Optimization objective function with optional scaling factor*/
  addPythonOption("OPT_OBJECTIVE");

  /* DESCRIPTION: Optimization constraint functions with optional scaling factor */
  addPythonOption("OPT_CONSTRAINT");

  /* DESCRIPTION: Finite different step for gradient estimation */
  addPythonOption("FIN_DIFF_STEP");

  /* DESCRIPTION: Verbosity of the python scripts to Stdout */
  addPythonOption("CONSOLE");

  /* DESCRIPTION: Flag specifying if the mesh was decomposed */
  addPythonOption("DECOMPOSED");

  /* END_CONFIG_OPTIONS */

}

void CConfig::SetParsing(char case_filename[MAX_STRING_SIZE]) {
  string text_line, option_name;
  ifstream case_file;
  vector<string> option_value;
  int rank = MASTER_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Read the configuration file ---*/
  case_file.open(case_filename, ios::in);

  if (case_file.fail()) {
    if (rank == MASTER_NODE) cout << endl << "The configuration file (.cfg) is missing!!" << endl << endl;
    exit(EXIT_FAILURE);
  }

  string errorString;

  int  err_count = 0;  // How many errors have we found in the config file
  int max_err_count = 30; // Maximum number of errors to print before stopping

  map<string,bool> included_options;

  /*--- Parse the configuration file and set the options ---*/
  while (getline (case_file,text_line)) {
    
    if (err_count >= max_err_count){
      errorString.append("too many errors. Stopping parse");

      //SU2MPI::PrintAndFinalize(errorString);
      cout << errorString << endl;
      throw(1);
    }
    
    if (TokenizeString(text_line, option_name, option_value)) {
      
      if (option_map.find(option_name) == option_map.end()){
        // See if it's a python option
          string newString;
          newString.append(option_name);
          newString.append(": invalid option name");
          newString.append(". Check current SU2 options in config_template.cfg.");
          newString.append("\n");
          errorString.append(newString);
          err_count++;
        continue;
      }

      // Option exists, check if the option has already been in the config file
      if (included_options.find(option_name) != included_options.end()){
        string newString;
        newString.append(option_name);
        newString.append(": option appears twice");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }


      // New found option. Add it to the map, and delete from all options
      included_options.insert(pair<string,bool>(option_name, true));
      all_options.erase(option_name);

      // Set the value and check error
      string out = option_map[option_name]->SetValue(option_value);
      if (out.compare("") != 0){
        errorString.append(out);
        errorString.append("\n");
        err_count++;
      }
    }
  }

  // See if there were any errors parsing the config file
  if (errorString.size() != 0){
//    SU2MPI::PrintAndFinalize(errorString);
    if (rank == MASTER_NODE) cout << errorString << endl;
    exit(EXIT_FAILURE);
  }

  // Set the default values for all of the options that weren't set
  for(map<string, bool>::iterator iter = all_options.begin(); iter != all_options.end(); ++iter){
    option_map[iter->first]->SetDefault();
  }

  case_file.close();
}

void CConfig::SetPostprocessing(unsigned short val_software, unsigned short val_izone, unsigned short val_nDim) {
    
    unsigned short iZone, iCFL;
    bool ideal_gas       = (Kind_FluidModel == STANDARD_AIR || Kind_FluidModel == IDEAL_GAS );
    bool standard_air       = (Kind_FluidModel == STANDARD_AIR);
    
    int size = SINGLE_NODE;
#ifdef HAVE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    
#ifndef HAVE_TECIO
    if (Output_FileFormat == TECPLOT_BINARY) {
        cout << "Tecplot binary file requested but SU2 was built without TecIO support." << "\n";
        Output_FileFormat = TECPLOT;
    }
#endif
    
    /*--- Store the SU2 module that we are executing. ---*/
    
    Kind_SU2 = val_software;
    
    /*--- Only SU2_CFD work with CGNS ---*/
    
    if ((Kind_SU2 != SU2_CFD) && (Kind_SU2 != SU2_SOL)) {
        if (Mesh_FileFormat == CGNS) {
            cout << "This software is not prepared for CGNS, please switch to SU2" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Deactivate the multigrid in the adjoint problem ---*/
  
    if (Adjoint && !MG_AdjointFlow) { nMultiLevel = 0; }
  
    /*--- Initialize non-physical points/reconstructions to zero ---*/
    
    Nonphys_Points   = 0;
    Nonphys_Reconstr = 0;
    
    /*--- Don't do any deformation if there is no Design variable information ---*/
    
    if (Design_Variable == NULL) {
        Design_Variable = new unsigned short [1];
        nDV = 1; Design_Variable[0] = NONE;
    }
    
    /*--- If multiple processors the grid should be always in native .su2 format ---*/
    
    if ((size > SINGLE_NODE) && ((Kind_SU2 == SU2_CFD) || (Kind_SU2 == SU2_SOL))) Mesh_FileFormat = SU2;
    
    /*--- Don't divide the numerical grid unless running SU2_DEF ---*/
    
    if (Kind_SU2 != SU2_DEF) Divide_Element = false;
    
    /*--- Identification of free-surface problem, this problems are always unsteady and incompressible. ---*/
    
    if (Kind_Regime == FREESURFACE) {
        if (Unsteady_Simulation != DT_STEPPING_2ND) Unsteady_Simulation = DT_STEPPING_1ST;
    }
    
    if (Kind_Solver == POISSON_EQUATION) {
        Unsteady_Simulation = STEADY;
    }
    
    /*--- Set the number of external iterations to 1 for the steady state problem ---*/
    
    if ((Kind_Solver == LINEAR_ELASTICITY) || (Kind_Solver == HEAT_EQUATION) ||
        (Kind_Solver == WAVE_EQUATION) || (Kind_Solver == POISSON_EQUATION)) {
        nMultiLevel = 0;
        if (Unsteady_Simulation == STEADY) nExtIter = 1;
        else Unst_nIntIter = 2;
    }
    
    /*--- Decide whether we should be writing unsteady solution files. ---*/
    
    if (Unsteady_Simulation == STEADY ||
        Unsteady_Simulation == TIME_STEPPING ||
        Unsteady_Simulation == TIME_SPECTRAL  ||
        Kind_Regime == FREESURFACE) { Wrt_Unsteady = false; }
    else { Wrt_Unsteady = true; }
    
    
    /*--- Check for Measurement System ---*/
    
    if (SystemMeasurements == US && !standard_air) {
        cout << "Only STANDARD_AIR fluid model can be used with US Measurement System" << endl;
        exit(EXIT_FAILURE);
    }
    
    /*--- Check for Convective scheme available for NICF ---*/
    
    if (!ideal_gas) {
        if (Kind_ConvNumScheme_Flow != SPACE_UPWIND){
            cout << "Only ROE Upwind scheme can be used for Not Ideal Compressible Fluids" << endl;
            exit(EXIT_FAILURE);
        }
        else {
            if (Kind_Upwind_Flow != ROE){
                cout << "Only ROE Upwind scheme can be used for Not Ideal Compressible Fluids" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /*--- Check for Boundary condition available for NICF ---*/
    
    if (!ideal_gas) {
        if (nMarker_Inlet != 0) {
            cout << "Riemann Boundary conditions must be used for inlet and outlet with Not Ideal Compressible Fluids " << endl;
            exit(EXIT_FAILURE);
        }
        if (nMarker_Outlet != 0) {
            cout << "Riemann Boundary conditions must be used outlet with Not Ideal Compressible Fluids " << endl;
            exit(EXIT_FAILURE);
        }
        
        if (nMarker_FarField != 0) {
            cout << "Far Field BC is not generalized for Not Ideal Compressible Fluids " << endl;
            exit(EXIT_FAILURE);
        }
        
    }
    
    /*--- Check for Boundary condition available for NICF ---*/
    
    if (ideal_gas) {
        if (SystemMeasurements == US && standard_air){
            if (Kind_ViscosityModel != SUTHERLAND){
                cout << "Only SUTHERLAND viscosity model can be used with US Measurement  " << endl;
                exit(EXIT_FAILURE);
            }
        }
        if (Kind_ConductivityModel != CONSTANT_PRANDTL ){
            cout << "Only CONSTANT_PRANDTL thermal conductivity model can be used with STANDARD_AIR and IDEAL_GAS" << endl;
            exit(EXIT_FAILURE);
        }
        
    }
    
    /*--- Set grid movement kind to NO_MOVEMENT if not specified, which means
     that we also set the Grid_Movement flag to false. We initialize to the
     number of zones here, because we are guaranteed to at least have one. ---*/
    
    if (Kind_GridMovement == NULL) {
        Kind_GridMovement = new unsigned short[nZone];
        for (unsigned short iZone = 0; iZone < nZone; iZone++ )
            Kind_GridMovement[iZone] = NO_MOVEMENT;
        if (Grid_Movement == true) {
            cout << "GRID_MOVEMENT = YES but no type provided in GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- If we're solving a purely steady problem with no prescribed grid
     movement (both rotating frame and moving walls can be steady), make sure that
     there is no grid motion ---*/
    
    if ((Kind_SU2 == SU2_CFD || Kind_SU2 == SU2_SOL) &&
        (Unsteady_Simulation == STEADY) &&
        ((Kind_GridMovement[ZONE_0] != MOVING_WALL) &&
         (Kind_GridMovement[ZONE_0] != ROTATING_FRAME)))
        Grid_Movement = false;
    
    /*--- If it is not specified, set the mesh motion mach number
     equal to the freestream value. ---*/
    
    if (Grid_Movement && Mach_Motion == 0.0)
        Mach_Motion = Mach;
    
    /*--- Set the boolean flag if we are in a rotating frame (source term). ---*/
    
    if (Grid_Movement && Kind_GridMovement[ZONE_0] == ROTATING_FRAME)
        Rotating_Frame = true;
    else
        Rotating_Frame = false;
    
    /*--- Check the number of moving markers against the number of grid movement
     types provided (should be equal, except that rigid motion and rotating frame
     do not depend on surface specification). ---*/
    
    if (Grid_Movement && (Kind_GridMovement[ZONE_0] != RIGID_MOTION) &&
        (Kind_GridMovement[ZONE_0] != ROTATING_FRAME) &&
        (nGridMovement != nMarker_Moving)) {
        cout << "Number of GRID_MOVEMENT_KIND must match number of MARKER_MOVING!!" << endl;
        exit(EXIT_FAILURE);
    }
    
    /*--- Make sure that there aren't more than one rigid motion or
     rotating frame specified in GRID_MOVEMENT_KIND. ---*/
    
    if (Grid_Movement && (Kind_GridMovement[ZONE_0] == RIGID_MOTION) &&
        (nGridMovement > 1)) {
        cout << "Can not support more than one type of rigid motion in GRID_MOVEMENT_KIND!!" << endl;
        exit(EXIT_FAILURE);
    }
    if (Grid_Movement && (Kind_GridMovement[ZONE_0] == ROTATING_FRAME) &&
        (nGridMovement > 1)) {
        cout << "Can not support more than one rotating frame in GRID_MOVEMENT_KIND!!" << endl;
        exit(EXIT_FAILURE);
    }
    
    /*--- In case the grid movement parameters have not been declared in the
     config file, set them equal to zero for safety. Also check to make sure
     that for each option, a value has been declared for each moving marker. ---*/
    
    unsigned short nMoving;
    if (nGridMovement > nZone) nMoving = nGridMovement;
    else nMoving = nZone;
    
    /*--- Motion Origin: ---*/
    if (Motion_Origin_X == NULL) {
        Motion_Origin_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Motion_Origin_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nMotion_Origin_X != nGridMovement)) {
            cout << "Length of MOTION_ORIGIN_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Motion_Origin_Y == NULL) {
        Motion_Origin_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Motion_Origin_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nMotion_Origin_Y != nGridMovement)) {
            cout << "Length of MOTION_ORIGIN_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Motion_Origin_Z == NULL) {
        Motion_Origin_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Motion_Origin_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nMotion_Origin_Z != nGridMovement)) {
            cout << "Length of MOTION_ORIGIN_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (MoveMotion_Origin == NULL) {
        MoveMotion_Origin = new unsigned short[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            MoveMotion_Origin[iZone] = 0;
    } else {
        if (Grid_Movement && (nMoveMotion_Origin != nGridMovement)) {
            cout << "Length of MOVE_MOTION_ORIGIN must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Translation: ---*/
    
    if (Translation_Rate_X == NULL) {
        Translation_Rate_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Translation_Rate_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nTranslation_Rate_X != nGridMovement)) {
            cout << "Length of TRANSLATION_RATE_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Translation_Rate_Y == NULL) {
        Translation_Rate_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Translation_Rate_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nTranslation_Rate_Y != nGridMovement)) {
            cout << "Length of TRANSLATION_RATE_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Translation_Rate_Z == NULL) {
        Translation_Rate_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Translation_Rate_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nTranslation_Rate_Z != nGridMovement)) {
            cout << "Length of TRANSLATION_RATE_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Rotation: ---*/
    
    if (Rotation_Rate_X == NULL) {
        Rotation_Rate_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Rotation_Rate_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nRotation_Rate_X != nGridMovement)) {
            cout << "Length of ROTATION_RATE_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Rotation_Rate_Y == NULL) {
        Rotation_Rate_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Rotation_Rate_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nRotation_Rate_Y != nGridMovement)) {
            cout << "Length of ROTATION_RATE_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Rotation_Rate_Z == NULL) {
        Rotation_Rate_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Rotation_Rate_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nRotation_Rate_Z != nGridMovement)) {
            cout << "Length of ROTATION_RATE_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Pitching: ---*/
    
    if (Pitching_Omega_X == NULL) {
        Pitching_Omega_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Omega_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Omega_X != nGridMovement)) {
            cout << "Length of PITCHING_OMEGA_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Pitching_Omega_Y == NULL) {
        Pitching_Omega_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Omega_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Omega_Y != nGridMovement)) {
            cout << "Length of PITCHING_OMEGA_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Pitching_Omega_Z == NULL) {
        Pitching_Omega_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Omega_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Omega_Z != nGridMovement)) {
            cout << "Length of PITCHING_OMEGA_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Pitching Amplitude: ---*/
    
    if (Pitching_Ampl_X == NULL) {
        Pitching_Ampl_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Ampl_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Ampl_X != nGridMovement)) {
            cout << "Length of PITCHING_AMPL_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Pitching_Ampl_Y == NULL) {
        Pitching_Ampl_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Ampl_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Ampl_Y != nGridMovement)) {
            cout << "Length of PITCHING_AMPL_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Pitching_Ampl_Z == NULL) {
        Pitching_Ampl_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Ampl_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Ampl_Z != nGridMovement)) {
            cout << "Length of PITCHING_AMPL_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Pitching Phase: ---*/
    
    if (Pitching_Phase_X == NULL) {
        Pitching_Phase_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Phase_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Phase_X != nGridMovement)) {
            cout << "Length of PITCHING_PHASE_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Pitching_Phase_Y == NULL) {
        Pitching_Phase_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Phase_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Phase_Y != nGridMovement)) {
            cout << "Length of PITCHING_PHASE_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Pitching_Phase_Z == NULL) {
        Pitching_Phase_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Pitching_Phase_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPitching_Phase_Z != nGridMovement)) {
            cout << "Length of PITCHING_PHASE_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Plunging: ---*/
    
    if (Plunging_Omega_X == NULL) {
        Plunging_Omega_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Plunging_Omega_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPlunging_Omega_X != nGridMovement)) {
            cout << "Length of PLUNGING_OMEGA_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Plunging_Omega_Y == NULL) {
        Plunging_Omega_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Plunging_Omega_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPlunging_Omega_Y != nGridMovement)) {
            cout << "Length of PLUNGING_OMEGA_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Plunging_Omega_Z == NULL) {
        Plunging_Omega_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Plunging_Omega_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPlunging_Omega_Z != nGridMovement)) {
            cout << "Length of PLUNGING_OMEGA_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Plunging Amplitude: ---*/
    
    if (Plunging_Ampl_X == NULL) {
        Plunging_Ampl_X = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Plunging_Ampl_X[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPlunging_Ampl_X != nGridMovement)) {
            cout << "Length of PLUNGING_AMPL_X must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Plunging_Ampl_Y == NULL) {
        Plunging_Ampl_Y = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Plunging_Ampl_Y[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPlunging_Ampl_Y != nGridMovement)) {
            cout << "Length of PLUNGING_AMPL_Y must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (Plunging_Ampl_Z == NULL) {
        Plunging_Ampl_Z = new double[nMoving];
        for (iZone = 0; iZone < nMoving; iZone++ )
            Plunging_Ampl_Z[iZone] = 0.0;
    } else {
        if (Grid_Movement && (nPlunging_Ampl_Z != nGridMovement)) {
            cout << "Length of PLUNGING_AMPL_Z must match GRID_MOVEMENT_KIND!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Use the various rigid-motion input frequencies to determine the period to be used with time-spectral cases.
     There are THREE types of motion to consider, namely: rotation, pitching, and plunging.
     The largest period of motion is the one to be used for time-spectral calculations. ---*/
    
    if (Unsteady_Simulation == TIME_SPECTRAL) {
        
        unsigned short N_MOTION_TYPES = 3;
        double *periods;
        periods = new double[N_MOTION_TYPES];
        
        /*--- rotation: ---*/
        
        double Omega_mag_rot = sqrt(pow(Rotation_Rate_X[ZONE_0],2)+pow(Rotation_Rate_Y[ZONE_0],2)+pow(Rotation_Rate_Z[ZONE_0],2));
        if (Omega_mag_rot > 0)
            periods[0] = 2*PI_NUMBER/Omega_mag_rot;
        else
            periods[0] = 0.0;
        
        /*--- pitching: ---*/
        
        double Omega_mag_pitch = sqrt(pow(Pitching_Omega_X[ZONE_0],2)+pow(Pitching_Omega_Y[ZONE_0],2)+pow(Pitching_Omega_Z[ZONE_0],2));
        if (Omega_mag_pitch > 0)
            periods[1] = 2*PI_NUMBER/Omega_mag_pitch;
        else
            periods[1] = 0.0;
        
        /*--- plunging: ---*/
        
        double Omega_mag_plunge = sqrt(pow(Plunging_Omega_X[ZONE_0],2)+pow(Plunging_Omega_Y[ZONE_0],2)+pow(Plunging_Omega_Z[ZONE_0],2));
        if (Omega_mag_plunge > 0)
            periods[2] = 2*PI_NUMBER/Omega_mag_plunge;
        else
            periods[2] = 0.0;
        
        /*--- determine which period is largest ---*/
        
        unsigned short iVar;
        TimeSpectral_Period = 0.0;
        for (iVar = 0; iVar < N_MOTION_TYPES; iVar++) {
            if (periods[iVar] > TimeSpectral_Period)
                TimeSpectral_Period = periods[iVar];
        }
        
        delete periods;
        
    }
    
    /*--- Initialize the RefOriginMoment Pointer ---*/
    
    RefOriginMoment = NULL;
    RefOriginMoment = new double[3];
    RefOriginMoment[0] = 0.0; RefOriginMoment[1] = 0.0; RefOriginMoment[2] = 0.0;
    
    /*--- In case the moment origin coordinates have not been declared in the
     config file, set them equal to zero for safety. Also check to make sure
     that for each marker, a value has been declared for the moment origin.
     Unless only one value was specified, then set this value for all the markers
     being monitored. ---*/
    
    unsigned short iMarker;
    
    if ((nRefOriginMoment_X != nRefOriginMoment_Y) || (nRefOriginMoment_X != nRefOriginMoment_Z) ) {
        cout << "ERROR: Length of REF_ORIGIN_MOMENT_X, REF_ORIGIN_MOMENT_Y and REF_ORIGIN_MOMENT_Z must be the same!!" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (RefOriginMoment_X == NULL) {
        RefOriginMoment_X = new double[nMarker_Monitoring];
        for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
            RefOriginMoment_X[iMarker] = 0.0;
    } else {
        if (nRefOriginMoment_X == 1) {
            
            double aux_RefOriginMoment_X = RefOriginMoment_X[0];
            delete [] RefOriginMoment_X;
            RefOriginMoment_X = new double[nMarker_Monitoring];
            nRefOriginMoment_X = nMarker_Monitoring;
            
            for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
                RefOriginMoment_X[iMarker] = aux_RefOriginMoment_X;
        }
        else if (nRefOriginMoment_X != nMarker_Monitoring) {
            cout << "ERROR: Length of REF_ORIGIN_MOMENT_X must match number of Monitoring Markers!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (RefOriginMoment_Y == NULL) {
        RefOriginMoment_Y = new double[nMarker_Monitoring];
        for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
            RefOriginMoment_Y[iMarker] = 0.0;
    } else {
        if (nRefOriginMoment_Y == 1) {
            
            double aux_RefOriginMoment_Y = RefOriginMoment_Y[0];
            delete [] RefOriginMoment_Y;
            RefOriginMoment_Y = new double[nMarker_Monitoring];
            nRefOriginMoment_Y = nMarker_Monitoring;
            
            for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
                RefOriginMoment_Y[iMarker] = aux_RefOriginMoment_Y;
        }
        else if (nRefOriginMoment_Y != nMarker_Monitoring) {
            cout << "ERROR: Length of REF_ORIGIN_MOMENT_Y must match number of Monitoring Markers!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (RefOriginMoment_Z == NULL) {
        RefOriginMoment_Z = new double[nMarker_Monitoring];
        for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
            RefOriginMoment_Z[iMarker] = 0.0;
    } else {
        if (nRefOriginMoment_Z == 1) {
            
            double aux_RefOriginMoment_Z = RefOriginMoment_Z[0];
            delete [] RefOriginMoment_Z;
            RefOriginMoment_Z = new double[nMarker_Monitoring];
            nRefOriginMoment_Z = nMarker_Monitoring;
            
            for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
                RefOriginMoment_Z[iMarker] = aux_RefOriginMoment_Z;
        }
        else if (nRefOriginMoment_Z != nMarker_Monitoring) {
            cout << "ERROR: Length of REF_ORIGIN_MOMENT_Z must match number of Monitoring Markers!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    /*--- Set the boolean flag if we are carrying out an aeroelastic simulation. ---*/
    
    if (Grid_Movement && (Kind_GridMovement[ZONE_0] == AEROELASTIC || Kind_GridMovement[ZONE_0] == AEROELASTIC_RIGID_MOTION))
        Aeroelastic_Simulation = true;
    else
        Aeroelastic_Simulation = false;
    
    /*--- Initializing the size for the solutions of the Aeroelastic problem. ---*/
    
    
    if (Grid_Movement && Aeroelastic_Simulation) {
        Aeroelastic_np1.resize(nMarker_Monitoring);
        Aeroelastic_n.resize(nMarker_Monitoring);
        Aeroelastic_n1.resize(nMarker_Monitoring);
        for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++) {
            Aeroelastic_np1[iMarker].resize(2);
            Aeroelastic_n[iMarker].resize(2);
            Aeroelastic_n1[iMarker].resize(2);
            for (int i =0; i<2; i++) {
                Aeroelastic_np1[iMarker][i].resize(2);
                Aeroelastic_n[iMarker][i].resize(2);
                Aeroelastic_n1[iMarker][i].resize(2);
                for (int j=0; j<2; j++) {
                    Aeroelastic_np1[iMarker][i][j] = 0.0;
                    Aeroelastic_n[iMarker][i][j] = 0.0;
                    Aeroelastic_n1[iMarker][i][j] = 0.0;
                }
            }
        }
    }
    
    /*--- Allocate memory for the plunge and pitch and initialized them to zero ---*/
    
    if (Grid_Movement && Aeroelastic_Simulation) {
        Aeroelastic_pitch = new double[nMarker_Monitoring];
        Aeroelastic_plunge = new double[nMarker_Monitoring];
        for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ ) {
            Aeroelastic_pitch[iMarker] = 0.0;
            Aeroelastic_plunge[iMarker] = 0.0;
        }
    }
    
    /*--- Fluid-Structure problems always have grid movement ---*/
    
    if (Kind_Solver == FLUID_STRUCTURE_EULER ||
        Kind_Solver == FLUID_STRUCTURE_NAVIER_STOKES) {
        if (Kind_Turb_Model != NONE)
            Kind_Solver = FLUID_STRUCTURE_RANS;
        Grid_Movement = true;
    }
    
    if (MGCycle == FULLMG_CYCLE) FinestMesh = nMultiLevel;
    else FinestMesh = MESH_0;
    
    if ((Kind_Solver == NAVIER_STOKES) &&
        (Kind_Turb_Model != NONE))
        Kind_Solver = RANS;
    
    if (Kind_Regime == FREESURFACE) GravityForce = true;
    
    Kappa_1st_Flow = Kappa_Flow[0];
    Kappa_2nd_Flow = Kappa_Flow[1];
    Kappa_4th_Flow = Kappa_Flow[2];
    Kappa_1st_AdjFlow = Kappa_AdjFlow[0];
    Kappa_2nd_AdjFlow = Kappa_AdjFlow[1];
    Kappa_4th_AdjFlow = Kappa_AdjFlow[2];
    Kappa_1st_LinFlow = Kappa_AdjFlow[0];
    Kappa_4th_LinFlow = Kappa_AdjFlow[1];
    
    // make the MG_PreSmooth, MG_PostSmooth, and MG_CorrecSmooth arrays consistent with nMultiLevel
    unsigned short * tmp_smooth = new unsigned short[nMultiLevel+1];
    
    if ((nMG_PreSmooth != nMultiLevel+1) && (nMG_PreSmooth != 0)) {
        if (nMG_PreSmooth > nMultiLevel+1) {
            
            // truncate by removing unnecessary elements at the end
            for (unsigned int i = 0; i <= nMultiLevel; i++)
                tmp_smooth[i] = MG_PreSmooth[i];
            delete [] MG_PreSmooth;
        } else {
            
            // add additional elements equal to last element
            for (unsigned int i = 0; i < nMG_PreSmooth; i++)
                tmp_smooth[i] = MG_PreSmooth[i];
            for (unsigned int i = nMG_PreSmooth; i <= nMultiLevel; i++)
                tmp_smooth[i] = MG_PreSmooth[nMG_PreSmooth-1];
            delete [] MG_PreSmooth;
        }
        
        nMG_PreSmooth = nMultiLevel+1;
        MG_PreSmooth = new unsigned short[nMG_PreSmooth];
        for (unsigned int i = 0; i < nMG_PreSmooth; i++)
            MG_PreSmooth[i] = tmp_smooth[i];
    }
    if ((nMultiLevel != 0) && (nMG_PreSmooth == 0)) {
        delete [] MG_PreSmooth;
        nMG_PreSmooth = nMultiLevel+1;
        MG_PreSmooth = new unsigned short[nMG_PreSmooth];
        for (unsigned int i = 0; i < nMG_PreSmooth; i++)
            MG_PreSmooth[i] = i+1;
    }
    
    if ((nMG_PostSmooth != nMultiLevel+1) && (nMG_PostSmooth != 0)) {
        if (nMG_PostSmooth > nMultiLevel+1) {
            // truncate by removing unnecessary elements at the end
            for (unsigned int i = 0; i <= nMultiLevel; i++)
                tmp_smooth[i] = MG_PostSmooth[i];
            delete [] MG_PostSmooth;
        } else {
            // add additional elements equal to last element
            for (unsigned int i = 0; i < nMG_PostSmooth; i++)
                tmp_smooth[i] = MG_PostSmooth[i];
            for (unsigned int i = nMG_PostSmooth; i <= nMultiLevel; i++)
                tmp_smooth[i] = MG_PostSmooth[nMG_PostSmooth-1];
            delete [] MG_PostSmooth;
        }
        nMG_PostSmooth = nMultiLevel+1;
        MG_PostSmooth = new unsigned short[nMG_PostSmooth];
        for (unsigned int i = 0; i < nMG_PostSmooth; i++)
            MG_PostSmooth[i] = tmp_smooth[i];
    }
    if ((nMultiLevel != 0) && (nMG_PostSmooth == 0)) {
        delete [] MG_PostSmooth;
        nMG_PostSmooth = nMultiLevel+1;
        MG_PostSmooth = new unsigned short[nMG_PostSmooth];
        for (unsigned int i = 0; i < nMG_PostSmooth; i++)
            MG_PostSmooth[i] = 0;
    }
    
    if ((nMG_CorrecSmooth != nMultiLevel+1) && (nMG_CorrecSmooth != 0)) {
        if (nMG_CorrecSmooth > nMultiLevel+1) {
            // truncate by removing unnecessary elements at the end
            for (unsigned int i = 0; i <= nMultiLevel; i++)
                tmp_smooth[i] = MG_CorrecSmooth[i];
            delete [] MG_CorrecSmooth;
        } else {
            // add additional elements equal to last element
            for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
                tmp_smooth[i] = MG_CorrecSmooth[i];
            for (unsigned int i = nMG_CorrecSmooth; i <= nMultiLevel; i++)
                tmp_smooth[i] = MG_CorrecSmooth[nMG_CorrecSmooth-1];
            delete [] MG_CorrecSmooth;
        }
        nMG_CorrecSmooth = nMultiLevel+1;
        MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
        for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
            MG_CorrecSmooth[i] = tmp_smooth[i];
    }
    if ((nMultiLevel != 0) && (nMG_CorrecSmooth == 0)) {
        delete [] MG_CorrecSmooth;
        nMG_CorrecSmooth = nMultiLevel+1;
        MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
        for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
            MG_CorrecSmooth[i] = 0;
    }
    
    // override MG Smooth parameters
    if (nMG_PreSmooth != 0)
        MG_PreSmooth[MESH_0] = 1;
    if (nMG_PostSmooth != 0) {
        MG_PostSmooth[MESH_0] = 0;
        MG_PostSmooth[nMultiLevel] = 0;
    }
    if (nMG_CorrecSmooth != 0)
        MG_CorrecSmooth[nMultiLevel] = 0;
    
    if (Restart) MGCycle = V_CYCLE;
    
    if (Adjoint) {
        if (Kind_Solver == EULER) Kind_Solver = ADJ_EULER;
        if (Kind_Solver == NAVIER_STOKES) Kind_Solver = ADJ_NAVIER_STOKES;
        if (Kind_Solver == RANS) Kind_Solver = ADJ_RANS;
        if (Kind_Solver == TNE2_EULER) Kind_Solver = ADJ_TNE2_EULER;
        if (Kind_Solver == TNE2_NAVIER_STOKES) Kind_Solver = ADJ_TNE2_NAVIER_STOKES;
    }
    
    if (Linearized) {
        if (Kind_Solver == EULER) Kind_Solver = LIN_EULER;
    }
    
    if (Unsteady_Simulation == TIME_STEPPING) {
        nMultiLevel = 0;
        MGCycle = V_CYCLE;
    }
    
    nCFL = nMultiLevel+1;
    CFL = new double[nCFL];
    CFL[0] = CFLFineGrid;
    if (Adjoint) CFL[0] = CFL[0] * CFLRedCoeff_AdjFlow;
    
    for (iCFL = 1; iCFL < nCFL; iCFL++)
        CFL[iCFL] = CFL[iCFL-1];
    
    if (nRKStep == 0) {
        nRKStep = 1;
        RK_Alpha_Step = new double[1]; RK_Alpha_Step[0] = 1.0;
    }
    
    if ((Kind_SU2 == SU2_CFD) && (Kind_Solver == NO_SOLVER)) {
        cout << "PHYSICAL_PROBLEM must be set in the configuration file" << endl;
        exit(EXIT_FAILURE);
    }
    
    /*--- Set a flag for viscous simulations ---*/
    
    Viscous = (( Kind_Solver == NAVIER_STOKES          ) ||
               ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
               ( Kind_Solver == RANS                   ) ||
               ( Kind_Solver == ADJ_RANS               ) ||
               ( Kind_Solver == TNE2_NAVIER_STOKES     ) ||
               ( Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   );
    
    
    /*--- Re-scale the length based parameters. The US system uses feet,
     but SU2 assumes that the grid is in inches ---*/
    
    if ((SystemMeasurements == US) && (Kind_SU2 == SU2_CFD)) {
        
        for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++) {
            RefOriginMoment_X[iMarker] = RefOriginMoment_X[iMarker]/12.0;
            RefOriginMoment_Y[iMarker] = RefOriginMoment_Y[iMarker]/12.0;
            RefOriginMoment_Z[iMarker] = RefOriginMoment_Z[iMarker]/12.0;
        }
        
        for (iMarker = 0; iMarker < nGridMovement; iMarker++) {
            Motion_Origin_X[iMarker] = Motion_Origin_X[iMarker]/12.0;
            Motion_Origin_Y[iMarker] = Motion_Origin_Y[iMarker]/12.0;
            Motion_Origin_Z[iMarker] = Motion_Origin_Z[iMarker]/12.0;
        }
        
        RefLengthMoment = RefLengthMoment/12.0;
        if (val_nDim == 2) RefAreaCoeff = RefAreaCoeff/12.0;
        else RefAreaCoeff = RefAreaCoeff/144.0;
        Length_Reynolds = Length_Reynolds/12.0;
        RefElemLength = RefElemLength/12.0;
        
        EA_IntLimit[0] = EA_IntLimit[0]/12.0;
        EA_IntLimit[1] = EA_IntLimit[1]/12.0;
        EA_IntLimit[2] = EA_IntLimit[2]/12.0;
        
        Section_Location[0] = Section_Location[0]/12.0;
        Section_Location[1] = Section_Location[1]/12.0;
        
    }
    
    /*--- Reacting flows initialization ---*/
    
    if (( Kind_Solver == TNE2_EULER             ) ||
        ( Kind_Solver == TNE2_NAVIER_STOKES     ) ||
        ( Kind_Solver == ADJ_TNE2_EULER         ) ||
        ( Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
        
        Kappa_1st_TNE2    = Kappa_TNE2[0];
        Kappa_2nd_TNE2    = Kappa_TNE2[1];
        Kappa_4th_TNE2    = Kappa_TNE2[2];
        Kappa_1st_AdjTNE2 = Kappa_AdjTNE2[0];
        Kappa_2nd_AdjTNE2 = Kappa_AdjTNE2[1];
        Kappa_4th_AdjTNE2 = Kappa_AdjTNE2[2];
        
        
        unsigned short maxEl = 0;
        unsigned short iSpecies, jSpecies, iEl;
        
        switch (Kind_GasModel) {
            case ONESPECIES:
                /*--- Define parameters of the gas model ---*/
                nSpecies    = 1;
                ionization  = false;
                
                /*--- Allocate vectors for gas properties ---*/
                Molar_Mass         = new double[nSpecies];
                CharVibTemp        = new double[nSpecies];
                RotationModes      = new double[nSpecies];
                Enthalpy_Formation = new double[nSpecies];
                Wall_Catalycity    = new double[nSpecies];
                Ref_Temperature    = new double[nSpecies];
                nElStates          = new unsigned short[nSpecies];
                
                
                
                MassFrac_FreeStream = new double[nSpecies];
                MassFrac_FreeStream[0] = 1.0;
                
                /*--- Assign gas properties ---*/
                // Rotational modes of energy storage
                RotationModes[0] = 2.0;
                // Molar mass [kg/kmol]
                Molar_Mass[0] = 14.0067+15.9994;
                // Characteristic vibrational temperatures for calculating e_vib [K]
                //CharVibTemp[0] = 3395.0;
                CharVibTemp[0] = 1000.0;
                // Formation enthalpy: (JANAF values, [KJ/Kmol])
                Enthalpy_Formation[0] = 0.0;					//N2
                // Reference temperature (JANAF values, [K])
                Ref_Temperature[0] = 0.0;
                
                /*        nElStates[0] = 0;
                 CharElTemp   = new double *[nSpecies];
                 degen        = new double *[nSpecies];
                 
                 OSPthetae    = new double[nElStates[0]];
                 OSPthetae[0] = 1.0;
                 OSPg         = new double[nElStates[0]];
                 OSPg[0]      = 1.0;
                 
                 CharElTemp[0] = OSPthetae;
                 degen[0] = OSPg;*/
                
                break;
                
            case N2:
                
                /*--- Define parameters of the gas model ---*/
                nSpecies    = 2;
                nReactions  = 2;
                ionization  = false;
                
                /*--- Allocate vectors for gas properties ---*/
                Wall_Catalycity      = new double[nSpecies];
                Molar_Mass           = new double[nSpecies];
                CharVibTemp          = new double[nSpecies];
                RotationModes        = new double[nSpecies];
                Enthalpy_Formation   = new double[nSpecies];
                Ref_Temperature      = new double[nSpecies];
                Diss                 = new double[nSpecies];
                ArrheniusCoefficient = new double[nReactions];
                ArrheniusEta         = new double[nReactions];
                ArrheniusTheta       = new double[nReactions];
                Tcf_a                = new double[nReactions];
                Tcf_b                = new double[nReactions];
                Tcb_a                = new double[nReactions];
                Tcb_b                = new double[nReactions];
                nElStates            = new unsigned short[nSpecies];
                Reactions = new int**[nReactions];
                for (unsigned short iRxn = 0; iRxn < nReactions; iRxn++) {
                    Reactions[iRxn] = new int*[2];
                    for (unsigned short ii = 0; ii < 2; ii++)
                        Reactions[iRxn][ii] = new int[6];
                }
                
                // Omega[iSpecies][jSpecies][iCoeff]
                Omega00 = new double**[nSpecies];
                Omega11 = new double**[nSpecies];
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    Omega00[iSpecies] = new double*[nSpecies];
                    Omega11[iSpecies] = new double*[nSpecies];
                    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
                        Omega00[iSpecies][jSpecies] = new double[4];
                        Omega11[iSpecies][jSpecies] = new double[4];
                    }
                }
                
                MassFrac_FreeStream = new double[nSpecies];
                MassFrac_FreeStream[0] = 0.99;
                MassFrac_FreeStream[1] = 0.01;
                
                /*--- Assign gas properties ---*/
                
                // Wall mass fractions for catalytic boundaries
                Wall_Catalycity[0] = 0.7;
                Wall_Catalycity[1] = 0.3;
                
                // Rotational modes of energy storage
                RotationModes[0] = 2.0;
                RotationModes[1] = 0.0;
                
                // Molar mass [kg/kmol]
                Molar_Mass[0] = 2.0*14.0067;
                Molar_Mass[1] = 14.0067;
                
                // Characteristic vibrational temperatures
                CharVibTemp[0] = 3395.0;
                CharVibTemp[1] = 0.0;
                
                // Formation enthalpy: (JANAF values [KJ/Kmol])
                Enthalpy_Formation[0] = 0.0;					//N2
                Enthalpy_Formation[1] = 472.683E3;		//N
                
                // Reference temperature (JANAF values, [K])
                Ref_Temperature[0] = 0.0;
                Ref_Temperature[1] = 298.15;
                
                // Number of electron states
                nElStates[0] = 15;                    // N2
                nElStates[1] = 3;                     // N
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
                    maxEl = max(maxEl, nElStates[iSpecies]);
                
                /*--- Allocate electron data arrays ---*/
                CharElTemp = new double*[nSpecies];
                degen      = new double*[nSpecies];
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    CharElTemp[iSpecies] = new double[maxEl];
                    degen[iSpecies]      = new double[maxEl];
                }
                
                /*--- Initialize the arrays ---*/
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    for (iEl = 0; iEl < maxEl; iEl++) {
                        CharElTemp[iSpecies][iEl] = 0.0;
                        degen[iSpecies][iEl] = 0.0;
                    }
                }
                
                /*--- Assign values to data structures ---*/
                // N2: 15 states
                CharElTemp[0][0]  = 0.000000000000000E+00;
                CharElTemp[0][1]  = 7.223156514095200E+04;
                CharElTemp[0][2]  = 8.577862640384000E+04;
                CharElTemp[0][3]  = 8.605026716160000E+04;
                CharElTemp[0][4]  = 9.535118627874400E+04;
                CharElTemp[0][5]  = 9.805635702203200E+04;
                CharElTemp[0][6]  = 9.968267656935200E+04;
                CharElTemp[0][7]  = 1.048976467715200E+05;
                CharElTemp[0][8]  = 1.116489555200000E+05;
                CharElTemp[0][9]  = 1.225836470400000E+05;
                CharElTemp[0][10] = 1.248856873600000E+05;
                CharElTemp[0][11] = 1.282476158188320E+05;
                CharElTemp[0][12] = 1.338060936000000E+05;
                CharElTemp[0][13] = 1.404296391107200E+05;
                CharElTemp[0][14] = 1.504958859200000E+05;
                degen[0][0]  = 1;
                degen[0][1]  = 3;
                degen[0][2]  = 6;
                degen[0][3]  = 6;
                degen[0][4]  = 3;
                degen[0][5]  = 1;
                degen[0][6]  = 2;
                degen[0][7]  = 2;
                degen[0][8]  = 5;
                degen[0][9]  = 1;
                degen[0][10] = 6;
                degen[0][11] = 6;
                degen[0][12] = 10;
                degen[0][13] = 6;
                degen[0][14] = 6;
                // N: 3 states
                CharElTemp[1][0] = 0.000000000000000E+00;
                CharElTemp[1][1] = 2.766469645581980E+04;
                CharElTemp[1][2] = 4.149309313560210E+04;
                degen[1][0] = 4;
                degen[1][1] = 10;
                degen[1][2] = 6;
                
                /*--- Set Arrhenius coefficients for chemical reactions ---*/
                // Pre-exponential factor
                ArrheniusCoefficient[0]  = 7.0E21;
                ArrheniusCoefficient[1]  = 3.0E22;
                // Rate-controlling temperature exponent
                ArrheniusEta[0]  = -1.60;
                ArrheniusEta[1]  = -1.60;
                // Characteristic temperature
                ArrheniusTheta[0] = 113200.0;
                ArrheniusTheta[1] = 113200.0;
                
                /*--- Set reaction maps ---*/
                // N2 + N2 -> 2N + N2
                Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;
                Reactions[0][1][0]=1;		Reactions[0][1][1]=1;		Reactions[0][1][2] =0;
                // N2 + N -> 2N + N
                Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;
                Reactions[1][1][0]=1;		Reactions[1][1][1]=1;		Reactions[1][1][2]=1;
                
                /*--- Set rate-controlling temperature exponents ---*/
                //  -----------  Tc = Ttr^a * Tve^b  -----------
                //
                // Forward Reactions
                //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
                //   Exchange:          a = 1,   b = 0
                //   Impact ionization: a = 0,   b = 1
                //
                // Backward Reactions
                //   Recomb ionization:      a = 0, b = 1
                //   Impact ionization:      a = 0, b = 1
                //   N2 impact dissociation: a = 0, b = 1
                //   Others:                 a = 1, b = 0
                Tcf_a[0] = 0.5; Tcf_b[0] = 0.5; Tcb_a[0] = 1;  Tcb_b[0] = 0;
                Tcf_a[1] = 0.5; Tcf_b[1] = 0.5; Tcb_a[1] = 1;  Tcb_b[1] = 0;
                
                /*--- Dissociation potential [KJ/kg] ---*/
                Diss[0] = 3.36E4;
                Diss[1] = 0.0;
                
                /*--- Collision integral data ---*/
                Omega00[0][0][0] = -6.0614558E-03;  Omega00[0][0][1] = 1.2689102E-01;   Omega00[0][0][2] = -1.0616948E+00;  Omega00[0][0][3] = 8.0955466E+02;
                Omega00[0][1][0] = -1.0796249E-02;  Omega00[0][1][1] = 2.2656509E-01;   Omega00[0][1][2] = -1.7910602E+00;  Omega00[0][1][3] = 4.0455218E+03;
                Omega00[1][0][0] = -1.0796249E-02;  Omega00[1][0][1] = 2.2656509E-01;   Omega00[1][0][2] = -1.7910602E+00;  Omega00[1][0][3] = 4.0455218E+03;
                Omega00[1][1][0] = -9.6083779E-03;  Omega00[1][1][1] = 2.0938971E-01;   Omega00[1][1][2] = -1.7386904E+00;  Omega00[1][1][3] = 3.3587983E+03;
                
                Omega11[0][0][0] = -7.6303990E-03;  Omega11[0][0][1] = 1.6878089E-01;   Omega11[0][0][2] = -1.4004234E+00;  Omega11[0][0][3] = 2.1427708E+03;
                Omega11[0][1][0] = -8.3493693E-03;  Omega11[0][1][1] = 1.7808911E-01;   Omega11[0][1][2] = -1.4466155E+00;  Omega11[0][1][3] = 1.9324210E+03;
                Omega11[1][0][0] = -8.3493693E-03;  Omega11[1][0][1] = 1.7808911E-01;   Omega11[1][0][2] = -1.4466155E+00;  Omega11[1][0][3] = 1.9324210E+03;
                Omega11[1][1][0] = -7.7439615E-03;  Omega11[1][1][1] = 1.7129007E-01;   Omega11[1][1][2] = -1.4809088E+00;  Omega11[1][1][3] = 2.1284951E+03;
                
                break;
                
            case AIR5:
                /*--- Define parameters of the gas model ---*/
                nSpecies    = 5;
                nReactions  = 17;
                ionization  = false;
                
                /*--- Allocate vectors for gas properties ---*/
                Wall_Catalycity      = new double[nSpecies];
                Molar_Mass           = new double[nSpecies];
                CharVibTemp          = new double[nSpecies];
                RotationModes        = new double[nSpecies];
                Enthalpy_Formation   = new double[nSpecies];
                Ref_Temperature      = new double[nSpecies];
                ArrheniusCoefficient = new double[nReactions];
                ArrheniusEta         = new double[nReactions];
                ArrheniusTheta       = new double[nReactions];
                Tcf_a                = new double[nReactions];
                Tcf_b                = new double[nReactions];
                Tcb_a                = new double[nReactions];
                Tcb_b                = new double[nReactions];
                nElStates            = new unsigned short[nSpecies];
                Reactions            = new int**[nReactions];
                for (unsigned short iRxn = 0; iRxn < nReactions; iRxn++) {
                    Reactions[iRxn] = new int*[2];
                    for (unsigned short ii = 0; ii < 2; ii++)
                        Reactions[iRxn][ii] = new int[6];
                }
                
                // Omega[iSpecies][jSpecies][iCoeff]
                Omega00 = new double**[nSpecies];
                Omega11 = new double**[nSpecies];
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    Omega00[iSpecies] = new double*[nSpecies];
                    Omega11[iSpecies] = new double*[nSpecies];
                    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
                        Omega00[iSpecies][jSpecies] = new double[4];
                        Omega11[iSpecies][jSpecies] = new double[4];
                    }
                }
                
                // Wall mass fractions for catalytic boundaries
                Wall_Catalycity[0] = 0.4;
                Wall_Catalycity[1] = 0.4;
                Wall_Catalycity[2] = 0.1;
                Wall_Catalycity[3] = 0.05;
                Wall_Catalycity[4] = 0.05;
                
                // Free stream mass fractions
                MassFrac_FreeStream = new double[nSpecies];
                MassFrac_FreeStream[0] = 0.78;
                MassFrac_FreeStream[1] = 0.19;
                MassFrac_FreeStream[2] = 0.01;
                MassFrac_FreeStream[3] = 0.01;
                MassFrac_FreeStream[4] = 0.01;
                
                /*--- Assign gas properties ---*/
                // Rotational modes of energy storage
                RotationModes[0] = 2.0;
                RotationModes[1] = 2.0;
                RotationModes[2] = 2.0;
                RotationModes[3] = 0.0;
                RotationModes[4] = 0.0;
                
                // Molar mass [kg/kmol]
                Molar_Mass[0] = 2.0*14.0067;
                Molar_Mass[1] = 2.0*15.9994;
                Molar_Mass[2] = 14.0067+15.9994;
                Molar_Mass[3] = 14.0067;
                Molar_Mass[4] = 15.9994;
                
                //Characteristic vibrational temperatures
                CharVibTemp[0] = 3395.0;
                CharVibTemp[1] = 2239.0;
                CharVibTemp[2] = 2817.0;
                CharVibTemp[3] = 0.0;
                CharVibTemp[4] = 0.0;
                
                // Formation enthalpy: (JANAF values [KJ/Kmol])
                Enthalpy_Formation[0] = 0.0;					//N2
                Enthalpy_Formation[1] = 0.0;					//O2
                Enthalpy_Formation[2] = 90.291E3;			//NO
                Enthalpy_Formation[3] = 472.683E3;		//N
                Enthalpy_Formation[4] = 249.173E3;		//O
                
                // Reference temperature (JANAF values, [K])
                Ref_Temperature[0] = 0.0;
                Ref_Temperature[1] = 0.0;
                Ref_Temperature[2] = 298.15;
                Ref_Temperature[3] = 298.15;
                Ref_Temperature[4] = 298.15;
                
                // Number of electron states
                nElStates[0] = 15;                    // N2
                nElStates[1] = 7;                     // O2
                nElStates[2] = 16;                    // NO
                nElStates[3] = 3;                     // N
                nElStates[4] = 5;                     // O
                /////
                
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
                    maxEl = max(maxEl, nElStates[iSpecies]);
                
                /*--- Allocate electron data arrays ---*/
                CharElTemp = new double*[nSpecies];
                degen      = new double*[nSpecies];
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    CharElTemp[iSpecies] = new double[maxEl];
                    degen[iSpecies]      = new double[maxEl];
                }
                
                /*--- Initialize the arrays ---*/
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    for (iEl = 0; iEl < maxEl; iEl++) {
                        CharElTemp[iSpecies][iEl] = 0.0;
                        degen[iSpecies][iEl] = 0.0;
                    }
                }
                
                //N2: 15 states
                CharElTemp[0][0]  = 0.000000000000000E+00;
                CharElTemp[0][1]  = 7.223156514095200E+04;
                CharElTemp[0][2]  = 8.577862640384000E+04;
                CharElTemp[0][3]  = 8.605026716160000E+04;
                CharElTemp[0][4]  = 9.535118627874400E+04;
                CharElTemp[0][5]  = 9.805635702203200E+04;
                CharElTemp[0][6]  = 9.968267656935200E+04;
                CharElTemp[0][7]  = 1.048976467715200E+05;
                CharElTemp[0][8]  = 1.116489555200000E+05;
                CharElTemp[0][9]  = 1.225836470400000E+05;
                CharElTemp[0][10] = 1.248856873600000E+05;
                CharElTemp[0][11] = 1.282476158188320E+05;
                CharElTemp[0][12] = 1.338060936000000E+05;
                CharElTemp[0][13] = 1.404296391107200E+05;
                CharElTemp[0][14] = 1.504958859200000E+05;
                degen[0][0]  = 1;
                degen[0][1]  = 3;
                degen[0][2]  = 6;
                degen[0][3]  = 6;
                degen[0][4]  = 3;
                degen[0][5]  = 1;
                degen[0][6]  = 2;
                degen[0][7]  = 2;
                degen[0][8]  = 5;
                degen[0][9]  = 1;
                degen[0][10] = 6;
                degen[0][11] = 6;
                degen[0][12] = 10;
                degen[0][13] = 6;
                degen[0][14] = 6;
                
                // O2: 7 states
                CharElTemp[1][0] = 0.000000000000000E+00;
                CharElTemp[1][1] = 1.139156019700800E+04;
                CharElTemp[1][2] = 1.898473947826400E+04;
                CharElTemp[1][3] = 4.755973576639200E+04;
                CharElTemp[1][4] = 4.991242097343200E+04;
                CharElTemp[1][5] = 5.092268575561600E+04;
                CharElTemp[1][6] = 7.189863255967200E+04;
                degen[1][0] = 3;
                degen[1][1] = 2;
                degen[1][2] = 1;
                degen[1][3] = 1;
                degen[1][4] = 6;
                degen[1][5] = 3;
                degen[1][6] = 3;
                
                // NO: 16 states
                CharElTemp[2][0]  = 0.000000000000000E+00;
                CharElTemp[2][1]  = 5.467345760000000E+04;
                CharElTemp[2][2]  = 6.317139627802400E+04;
                CharElTemp[2][3]  = 6.599450342445600E+04;
                CharElTemp[2][4]  = 6.906120960000000E+04;
                CharElTemp[2][5]  = 7.049998480000000E+04;
                CharElTemp[2][6]  = 7.491055017560000E+04;
                CharElTemp[2][7]  = 7.628875293968000E+04;
                CharElTemp[2][8]  = 8.676188537552000E+04;
                CharElTemp[2][9]  = 8.714431182368000E+04;
                CharElTemp[2][10] = 8.886077063728000E+04;
                CharElTemp[2][11] = 8.981755614528000E+04;
                CharElTemp[2][12] = 8.988445919208000E+04;
                CharElTemp[2][13] = 9.042702132000000E+04;
                CharElTemp[2][14] = 9.064283760000000E+04;
                CharElTemp[2][15] = 9.111763341600000E+04;
                degen[2][0]  = 4;
                degen[2][1]  = 8;
                degen[2][2]  = 2;
                degen[2][3]  = 4;
                degen[2][4]  = 4;
                degen[2][5]  = 4;
                degen[2][6]  = 4;
                degen[2][7]  = 2;
                degen[2][8]  = 4;
                degen[2][9]  = 2;
                degen[2][10] = 4;
                degen[2][11] = 4;
                degen[2][12] = 2;
                degen[2][13] = 2;
                degen[2][14] = 2;
                degen[2][15] = 4;
                
                // N: 3 states
                CharElTemp[3][0] = 0.000000000000000E+00;
                CharElTemp[3][1] = 2.766469645581980E+04;
                CharElTemp[3][2] = 4.149309313560210E+04;
                degen[3][0] = 4;
                degen[3][1] = 10;
                degen[3][2] = 6;
                
                // O: 5 states
                CharElTemp[4][0] = 0.000000000000000E+00;
                CharElTemp[4][1] = 2.277077570280000E+02;
                CharElTemp[4][2] = 3.265688785704000E+02;
                CharElTemp[4][3] = 2.283028632262240E+04;
                CharElTemp[4][4] = 4.861993036434160E+04;
                degen[4][0] = 5;
                degen[4][1] = 3;
                degen[4][2] = 1;
                degen[4][3] = 5;
                degen[4][4] = 1;
                
                /*--- Set reaction maps ---*/
                // N2 dissociation
                Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=3;		Reactions[0][1][1]=3;		Reactions[0][1][2] =0;
                Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=3;		Reactions[1][1][1]=3;		Reactions[1][1][2] =1;
                Reactions[2][0][0]=0;		Reactions[2][0][1]=2;		Reactions[2][0][2]=nSpecies;		Reactions[2][1][0]=3;		Reactions[2][1][1]=3;		Reactions[2][1][2] =2;
                Reactions[3][0][0]=0;		Reactions[3][0][1]=3;		Reactions[3][0][2]=nSpecies;		Reactions[3][1][0]=3;		Reactions[3][1][1]=3;		Reactions[3][1][2] =3;
                Reactions[4][0][0]=0;		Reactions[4][0][1]=4;		Reactions[4][0][2]=nSpecies;		Reactions[4][1][0]=3;		Reactions[4][1][1]=3;		Reactions[4][1][2] =4;
                // O2 dissociation
                Reactions[5][0][0]=1;		Reactions[5][0][1]=0;		Reactions[5][0][2]=nSpecies;		Reactions[5][1][0]=4;		Reactions[5][1][1]=4;		Reactions[5][1][2] =0;
                Reactions[6][0][0]=1;		Reactions[6][0][1]=1;		Reactions[6][0][2]=nSpecies;		Reactions[6][1][0]=4;		Reactions[6][1][1]=4;		Reactions[6][1][2] =1;
                Reactions[7][0][0]=1;		Reactions[7][0][1]=2;		Reactions[7][0][2]=nSpecies;		Reactions[7][1][0]=4;		Reactions[7][1][1]=4;		Reactions[7][1][2] =2;
                Reactions[8][0][0]=1;		Reactions[8][0][1]=3;		Reactions[8][0][2]=nSpecies;		Reactions[8][1][0]=4;		Reactions[8][1][1]=4;		Reactions[8][1][2] =3;
                Reactions[9][0][0]=1;		Reactions[9][0][1]=4;		Reactions[9][0][2]=nSpecies;		Reactions[9][1][0]=4;		Reactions[9][1][1]=4;		Reactions[9][1][2] =4;
                // NO dissociation
                Reactions[10][0][0]=2;		Reactions[10][0][1]=0;		Reactions[10][0][2]=nSpecies;		Reactions[10][1][0]=3;		Reactions[10][1][1]=4;		Reactions[10][1][2] =0;
                Reactions[11][0][0]=2;		Reactions[11][0][1]=1;		Reactions[11][0][2]=nSpecies;		Reactions[11][1][0]=3;		Reactions[11][1][1]=4;		Reactions[11][1][2] =1;
                Reactions[12][0][0]=2;		Reactions[12][0][1]=2;		Reactions[12][0][2]=nSpecies;		Reactions[12][1][0]=3;		Reactions[12][1][1]=4;		Reactions[12][1][2] =2;
                Reactions[13][0][0]=2;		Reactions[13][0][1]=3;		Reactions[13][0][2]=nSpecies;		Reactions[13][1][0]=3;		Reactions[13][1][1]=4;		Reactions[13][1][2] =3;
                Reactions[14][0][0]=2;		Reactions[14][0][1]=4;		Reactions[14][0][2]=nSpecies;		Reactions[14][1][0]=3;		Reactions[14][1][1]=4;		Reactions[14][1][2] =4;
                // N2 + O -> NO + N
                Reactions[15][0][0]=0;		Reactions[15][0][1]=4;		Reactions[15][0][2]=nSpecies;		Reactions[15][1][0]=2;		Reactions[15][1][1]=3;		Reactions[15][1][2]= nSpecies;
                // NO + O -> O2 + N
                Reactions[16][0][0]=2;		Reactions[16][0][1]=4;		Reactions[16][0][2]=nSpecies;		Reactions[16][1][0]=1;		Reactions[16][1][1]=3;		Reactions[16][1][2]= nSpecies;
                
                /*--- Set Arrhenius coefficients for reactions ---*/
                // Pre-exponential factor
                ArrheniusCoefficient[0]  = 7.0E21;
                ArrheniusCoefficient[1]  = 7.0E21;
                ArrheniusCoefficient[2]  = 7.0E21;
                ArrheniusCoefficient[3]  = 3.0E22;
                ArrheniusCoefficient[4]  = 3.0E22;
                ArrheniusCoefficient[5]  = 2.0E21;
                ArrheniusCoefficient[6]  = 2.0E21;
                ArrheniusCoefficient[7]  = 2.0E21;
                ArrheniusCoefficient[8]  = 1.0E22;
                ArrheniusCoefficient[9]  = 1.0E22;
                ArrheniusCoefficient[10] = 5.0E15;
                ArrheniusCoefficient[11] = 5.0E15;
                ArrheniusCoefficient[12] = 5.0E15;
                ArrheniusCoefficient[13] = 1.1E17;
                ArrheniusCoefficient[14] = 1.1E17;
                ArrheniusCoefficient[15] = 6.4E17;
                ArrheniusCoefficient[16] = 8.4E12;
                
                // Rate-controlling temperature exponent
                ArrheniusEta[0]  = -1.60;
                ArrheniusEta[1]  = -1.60;
                ArrheniusEta[2]  = -1.60;
                ArrheniusEta[3]  = -1.60;
                ArrheniusEta[4]  = -1.60;
                ArrheniusEta[5]  = -1.50;
                ArrheniusEta[6]  = -1.50;
                ArrheniusEta[7]  = -1.50;
                ArrheniusEta[8]  = -1.50;
                ArrheniusEta[9]  = -1.50;
                ArrheniusEta[10] = 0.0;
                ArrheniusEta[11] = 0.0;
                ArrheniusEta[12] = 0.0;
                ArrheniusEta[13] = 0.0;
                ArrheniusEta[14] = 0.0;
                ArrheniusEta[15] = -1.0;
                ArrheniusEta[16] = 0.0;
                
                // Characteristic temperature
                ArrheniusTheta[0]  = 113200.0;
                ArrheniusTheta[1]  = 113200.0;
                ArrheniusTheta[2]  = 113200.0;
                ArrheniusTheta[3]  = 113200.0;
                ArrheniusTheta[4]  = 113200.0;
                ArrheniusTheta[5]  = 59500.0;
                ArrheniusTheta[6]  = 59500.0;
                ArrheniusTheta[7]  = 59500.0;
                ArrheniusTheta[8]  = 59500.0;
                ArrheniusTheta[9]  = 59500.0;
                ArrheniusTheta[10] = 75500.0;
                ArrheniusTheta[11] = 75500.0;
                ArrheniusTheta[12] = 75500.0;
                ArrheniusTheta[13] = 75500.0;
                ArrheniusTheta[14] = 75500.0;
                ArrheniusTheta[15] = 38400.0;
                ArrheniusTheta[16] = 19450.0;
                
                /*--- Set rate-controlling temperature exponents ---*/
                //  -----------  Tc = Ttr^a * Tve^b  -----------
                //
                // Forward Reactions
                //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
                //   Exchange:          a = 1,   b = 0
                //   Impact ionization: a = 0,   b = 1
                //
                // Backward Reactions
                //   Recomb ionization:      a = 0, b = 1
                //   Impact ionization:      a = 0, b = 1
                //   N2 impact dissociation: a = 0, b = 1
                //   Others:                 a = 1, b = 0
                Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;  Tcb_b[0] = 0;
                Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;  Tcb_b[1] = 0;
                Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;  Tcb_b[2] = 0;
                Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;  Tcb_b[3] = 0;
                Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;  Tcb_b[4] = 0;
                
                Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;  Tcb_b[5] = 0;
                Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;  Tcb_b[6] = 0;
                Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;  Tcb_b[7] = 0;
                Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;  Tcb_b[8] = 0;
                Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;  Tcb_b[9] = 0;
                
                Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;  Tcb_b[10] = 0;
                Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;  Tcb_b[11] = 0;
                Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;  Tcb_b[12] = 0;
                Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;  Tcb_b[13] = 0;
                Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;  Tcb_b[14] = 0;
                
                Tcf_a[15] = 1.0; Tcf_b[15] = 0.0; Tcb_a[15] = 1;  Tcb_b[15] = 0;
                Tcf_a[16] = 1.0; Tcf_b[16] = 0.0; Tcb_a[16] = 1;  Tcb_b[16] = 0;
                
                /*--- Collision integral data ---*/
                // Omega(0,0) ----------------------
                //N2
                Omega00[0][0][0] = -6.0614558E-03;  Omega00[0][0][1] = 1.2689102E-01;   Omega00[0][0][2] = -1.0616948E+00;  Omega00[0][0][3] = 8.0955466E+02;
                Omega00[0][1][0] = -3.7959091E-03;  Omega00[0][1][1] = 9.5708295E-02;   Omega00[0][1][2] = -1.0070611E+00;  Omega00[0][1][3] = 8.9392313E+02;
                Omega00[0][2][0] = -1.9295666E-03;  Omega00[0][2][1] = 2.7995735E-02;   Omega00[0][2][2] = -3.1588514E-01;  Omega00[0][2][3] = 1.2880734E+02;
                Omega00[0][3][0] = -1.0796249E-02;  Omega00[0][3][1] = 2.2656509E-01;   Omega00[0][3][2] = -1.7910602E+00;  Omega00[0][3][3] = 4.0455218E+03;
                Omega00[0][4][0] = -2.7244269E-03;  Omega00[0][4][1] = 6.9587171E-02;   Omega00[0][4][2] = -7.9538667E-01;  Omega00[0][4][3] = 4.0673730E+02;
                //O2
                Omega00[1][0][0] = -3.7959091E-03;  Omega00[1][0][1] = 9.5708295E-02;   Omega00[1][0][2] = -1.0070611E+00;  Omega00[1][0][3] = 8.9392313E+02;
                Omega00[1][1][0] = -8.0682650E-04;  Omega00[1][1][1] = 1.6602480E-02;   Omega00[1][1][2] = -3.1472774E-01;  Omega00[1][1][3] = 1.4116458E+02;
                Omega00[1][2][0] = -6.4433840E-04;  Omega00[1][2][1] = 8.5378580E-03;   Omega00[1][2][2] = -2.3225102E-01;  Omega00[1][2][3] = 1.1371608E+02;
                Omega00[1][3][0] = -1.1453028E-03;  Omega00[1][3][1] = 1.2654140E-02;   Omega00[1][3][2] = -2.2435218E-01;  Omega00[1][3][3] = 7.7201588E+01;
                Omega00[1][4][0] = -4.8405803E-03;  Omega00[1][4][1] = 1.0297688E-01;   Omega00[1][4][2] = -9.6876576E-01;  Omega00[1][4][3] = 6.1629812E+02;
                //NO
                Omega00[2][0][0] = -1.9295666E-03;  Omega00[2][0][1] = 2.7995735E-02;   Omega00[2][0][2] = -3.1588514E-01;  Omega00[2][0][3] = 1.2880734E+02;
                Omega00[2][1][0] = -6.4433840E-04;  Omega00[2][1][1] = 8.5378580E-03;   Omega00[2][1][2] = -2.3225102E-01;  Omega00[2][1][3] = 1.1371608E+02;
                Omega00[2][2][0] = -0.0000000E+00;  Omega00[2][2][1] = -1.1056066E-02;  Omega00[2][2][2] = -5.9216250E-02;  Omega00[2][2][3] = 7.2542367E+01;
                Omega00[2][3][0] = -1.5770918E-03;  Omega00[2][3][1] = 1.9578381E-02;   Omega00[2][3][2] = -2.7873624E-01;  Omega00[2][3][3] = 9.9547944E+01;
                Omega00[2][4][0] = -1.0885815E-03;  Omega00[2][4][1] = 1.1883688E-02;   Omega00[2][4][2] = -2.1844909E-01;  Omega00[2][4][3] = 7.5512560E+01;
                //N
                Omega00[3][0][0] = -1.0796249E-02;  Omega00[3][0][1] = 2.2656509E-01;   Omega00[3][0][2] = -1.7910602E+00;  Omega00[3][0][3] = 4.0455218E+03;
                Omega00[3][1][0] = -1.1453028E-03;  Omega00[3][1][1] = 1.2654140E-02;   Omega00[3][1][2] = -2.2435218E-01;  Omega00[3][1][3] = 7.7201588E+01;
                Omega00[3][2][0] = -1.5770918E-03;  Omega00[3][2][1] = 1.9578381E-02;   Omega00[3][2][2] = -2.7873624E-01;  Omega00[3][2][3] = 9.9547944E+01;
                Omega00[3][3][0] = -9.6083779E-03;  Omega00[3][3][1] = 2.0938971E-01;   Omega00[3][3][2] = -1.7386904E+00;  Omega00[3][3][3] = 3.3587983E+03;
                Omega00[3][4][0] = -7.8147689E-03;  Omega00[3][4][1] = 1.6792705E-01;   Omega00[3][4][2] = -1.4308628E+00;  Omega00[3][4][3] = 1.6628859E+03;
                //O
                Omega00[4][0][0] = -2.7244269E-03;  Omega00[4][0][1] = 6.9587171E-02;   Omega00[4][0][2] = -7.9538667E-01;  Omega00[4][0][3] = 4.0673730E+02;
                Omega00[4][1][0] = -4.8405803E-03;  Omega00[4][1][1] = 1.0297688E-01;   Omega00[4][1][2] = -9.6876576E-01;  Omega00[4][1][3] = 6.1629812E+02;
                Omega00[4][2][0] = -1.0885815E-03;  Omega00[4][2][1] = 1.1883688E-02;   Omega00[4][2][2] = -2.1844909E-01;  Omega00[4][2][3] = 7.5512560E+01;
                Omega00[4][3][0] = -7.8147689E-03;  Omega00[4][3][1] = 1.6792705E-01;   Omega00[4][3][2] = -1.4308628E+00;  Omega00[4][3][3] = 1.6628859E+03;
                Omega00[4][4][0] = -6.4040535E-03;  Omega00[4][4][1] = 1.4629949E-01;   Omega00[4][4][2] = -1.3892121E+00;  Omega00[4][4][3] = 2.0903441E+03;
                
                // Omega(1,1) ----------------------
                //N2
                Omega11[0][0][0] = -7.6303990E-03;  Omega11[0][0][1] = 1.6878089E-01;   Omega11[0][0][2] = -1.4004234E+00;  Omega11[0][0][3] = 2.1427708E+03;
                Omega11[0][1][0] = -8.0457321E-03;  Omega11[0][1][1] = 1.9228905E-01;   Omega11[0][1][2] = -1.7102854E+00;  Omega11[0][1][3] = 5.2213857E+03;
                Omega11[0][2][0] = -6.8237776E-03;  Omega11[0][2][1] = 1.4360616E-01;   Omega11[0][2][2] = -1.1922240E+00;  Omega11[0][2][3] = 1.2433086E+03;
                Omega11[0][3][0] = -8.3493693E-03;  Omega11[0][3][1] = 1.7808911E-01;   Omega11[0][3][2] = -1.4466155E+00;  Omega11[0][3][3] = 1.9324210E+03;
                Omega11[0][4][0] = -8.3110691E-03;  Omega11[0][4][1] = 1.9617877E-01;   Omega11[0][4][2] = -1.7205427E+00;  Omega11[0][4][3] = 4.0812829E+03;
                //O2
                Omega11[1][0][0] = -8.0457321E-03;  Omega11[1][0][1] = 1.9228905E-01;   Omega11[1][0][2] = -1.7102854E+00;  Omega11[1][0][3] = 5.2213857E+03;
                Omega11[1][1][0] = -6.2931612E-03;  Omega11[1][1][1] = 1.4624645E-01;   Omega11[1][1][2] = -1.3006927E+00;  Omega11[1][1][3] = 1.8066892E+03;
                Omega11[1][2][0] = -6.8508672E-03;  Omega11[1][2][1] = 1.5524564E-01;   Omega11[1][2][2] = -1.3479583E+00;  Omega11[1][2][3] = 2.0037890E+03;
                Omega11[1][3][0] = -1.0608832E-03;  Omega11[1][3][1] = 1.1782595E-02;   Omega11[1][3][2] = -2.1246301E-01;  Omega11[1][3][3] = 8.4561598E+01;
                Omega11[1][4][0] = -3.7969686E-03;  Omega11[1][4][1] = 7.6789981E-02;   Omega11[1][4][2] = -7.3056809E-01;  Omega11[1][4][3] = 3.3958171E+02;
                //NO
                Omega11[2][0][0] = -6.8237776E-03;  Omega11[2][0][1] = 1.4360616E-01;   Omega11[2][0][2] = -1.1922240E+00;  Omega11[2][0][3] = 1.2433086E+03;
                Omega11[2][1][0] = -6.8508672E-03;  Omega11[2][1][1] = 1.5524564E-01;   Omega11[2][1][2] = -1.3479583E+00;  Omega11[2][1][3] = 2.0037890E+03;
                Omega11[2][2][0] = -7.4942466E-03;  Omega11[2][2][1] = 1.6626193E-01;   Omega11[2][2][2] = -1.4107027E+00;  Omega11[2][2][3] = 2.3097604E+03;
                Omega11[2][3][0] = -1.4719259E-03;  Omega11[2][3][1] = 1.8446968E-02;   Omega11[2][3][2] = -2.6460411E-01;  Omega11[2][3][3] = 1.0911124E+02;
                Omega11[2][4][0] = -1.0066279E-03;  Omega11[2][4][1] = 1.1029264E-02;   Omega11[2][4][2] = -2.0671266E-01;  Omega11[2][4][3] = 8.2644384E+01;
                //N
                Omega11[3][0][0] = -8.3493693E-03;  Omega11[3][0][1] = 1.7808911E-01;   Omega11[3][0][2] = -1.4466155E+00;  Omega11[3][0][3] = 1.9324210E+03;
                Omega11[3][1][0] = -1.0608832E-03;  Omega11[3][1][1] = 1.1782595E-02;   Omega11[3][1][2] = -2.1246301E-01;  Omega11[3][1][3] = 8.4561598E+01;
                Omega11[3][2][0] = -1.4719259E-03;  Omega11[3][2][1] = 1.8446968E-02;   Omega11[3][2][2] = -2.6460411E-01;  Omega11[3][2][3] = 1.0911124E+02;
                Omega11[3][3][0] = -7.7439615E-03;  Omega11[3][3][1] = 1.7129007E-01;   Omega11[3][3][2] = -1.4809088E+00;  Omega11[3][3][3] = 2.1284951E+03;
                Omega11[3][4][0] = -5.0478143E-03;  Omega11[3][4][1] = 1.0236186E-01;   Omega11[3][4][2] = -9.0058935E-01;  Omega11[3][4][3] = 4.4472565E+02;
                //O
                Omega11[4][0][0] = -8.3110691E-03;  Omega11[4][0][1] = 1.9617877E-01;   Omega11[4][0][2] = -1.7205427E+00;  Omega11[4][0][3] = 4.0812829E+03;
                Omega11[4][1][0] = -3.7969686E-03;  Omega11[4][1][1] = 7.6789981E-02;   Omega11[4][1][2] = -7.3056809E-01;  Omega11[4][1][3] = 3.3958171E+02;
                Omega11[4][2][0] = -1.0066279E-03;  Omega11[4][2][1] = 1.1029264E-02;   Omega11[4][2][2] = -2.0671266E-01;  Omega11[4][2][3] = 8.2644384E+01;
                Omega11[4][3][0] = -5.0478143E-03;  Omega11[4][3][1] = 1.0236186E-01;   Omega11[4][3][2] = -9.0058935E-01;  Omega11[4][3][3] = 4.4472565E+02;
                Omega11[4][4][0] = -4.2451096E-03;  Omega11[4][4][1] = 9.6820337E-02;   Omega11[4][4][2] = -9.9770795E-01;  Omega11[4][4][3] = 8.3320644E+02;
                
                break;
        }
    }
    
    /*--- Check for constant lift mode  ---*/
    if (Fixed_CL_Mode) {
        
        /*--- The initial AoA will be taken as the value input in the 
         config file (the default is zero). ---*/
        
        /*--- We will force the use of the cauchy convergence criteria for
         constant lift mode. ---*/
        
        ConvCriteria = CAUCHY;
        Cauchy_Func_Flow = LIFT_COEFFICIENT;
        
        /*--- Initialize the update flag for the AoA with each iteration to false ---*/
        
        Update_AoA = false;
        
    }
    
    /*--- Check for 2nd order w/ limiting for JST and correct ---*/
    
    if ((Kind_ConvNumScheme_Flow == SPACE_CENTERED) && (Kind_Centered_Flow == JST) && (SpatialOrder_Flow == SECOND_ORDER_LIMITER))
        SpatialOrder_Flow = SECOND_ORDER;
    
    if ((Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED) && (Kind_Centered_AdjFlow == JST) && (SpatialOrder_AdjFlow == SECOND_ORDER_LIMITER))
        SpatialOrder_AdjFlow = SECOND_ORDER;
    
    delete [] tmp_smooth;
    
}

void CConfig::SetMarkers(unsigned short val_software, unsigned short val_izone) {

  /*--- Identify the solvers that work in serial ---*/

#ifndef HAVE_MPI
  nDomain = SINGLE_NODE;
#else
  if (val_software != SU2_MSH)
    MPI_Comm_size(MPI_COMM_WORLD, (int*)&nDomain);   // any issue with type conversion here? MC
  else
    nDomain = SINGLE_NODE;
#endif

  /*--- Boundary (marker) treatment ---*/
  
  nMarker_All = nMarker_Euler + nMarker_FarField + nMarker_SymWall +
  nMarker_PerBound + nMarker_NearFieldBound + nMarker_Supersonic_Inlet +
  nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Riemann+ nMarker_Inlet +
  nMarker_Outlet + nMarker_Isothermal + nMarker_IsothermalCatalytic +
  nMarker_IsothermalNonCatalytic + nMarker_HeatFlux + nMarker_HeatFluxCatalytic +
  nMarker_HeatFluxNonCatalytic + nMarker_EngineInflow + nMarker_EngineBleed + nMarker_EngineExhaust +
  nMarker_Dirichlet_Elec + nMarker_Displacement + nMarker_Load +
  nMarker_FlowLoad + nMarker_Pressure + nMarker_Custom +
  nMarker_ActDisk_Inlet + nMarker_ActDisk_Outlet + nMarker_Out_1D + 2*nDomain;

  Marker_All_TagBound   = new string[nMarker_All+2];			    // Store the tag that correspond with each marker.
  Marker_All_SendRecv   = new short[nMarker_All+2];						// +#domain (send), -#domain (receive) or 0 (neither send nor receive).
  Marker_All_KindBC     = new unsigned short[nMarker_All+2];	// Store the kind of boundary condition.
  Marker_All_Monitoring = new unsigned short[nMarker_All+2];	// Store whether the boundary should be monitored.
  Marker_All_Designing  = new unsigned short[nMarker_All+2];  // Store whether the boundary should be designed.
  Marker_All_Plotting   = new unsigned short[nMarker_All+2];	// Store whether the boundary should be plotted.
  Marker_All_GeoEval    = new unsigned short[nMarker_All+2];	// Store whether the boundary should be geometry evaluation.
  Marker_All_DV         = new unsigned short[nMarker_All+2];	// Store whether the boundary should be affected by design variables.
  Marker_All_Moving     = new unsigned short[nMarker_All+2];	// Store whether the boundary should be in motion.
  Marker_All_PerBound   = new short[nMarker_All+2];						// Store whether the boundary belongs to a periodic boundary.
  Marker_All_Out_1D     = new unsigned short[nMarker_All+2];           // Store whether the boundary belongs to a 1-d output boundary.

  unsigned short iMarker_All, iMarker_Config, iMarker_Euler, iMarker_Custom,
  iMarker_FarField, iMarker_SymWall, iMarker_Pressure, iMarker_PerBound,
  iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Dirichlet,
  iMarker_Inlet, iMarker_Riemann, iMarker_Outlet, iMarker_Isothermal, iMarker_IsothermalCatalytic,
  iMarker_IsothermalNonCatalytic, iMarker_HeatFlux, iMarker_HeatFluxNoncatalytic,
  iMarker_HeatFluxCatalytic, iMarker_EngineInflow, iMarker_EngineBleed, iMarker_EngineExhaust,
  iMarker_Displacement, iMarker_Load, iMarker_FlowLoad, iMarker_Neumann,
  iMarker_Monitoring, iMarker_Designing, iMarker_GeoEval, iMarker_Plotting,
  iMarker_DV, iMarker_Moving, iMarker_Supersonic_Inlet,
  iMarker_ActDisk_Inlet, iMarker_ActDisk_Outlet, iMarker_Out_1D;

  for (iMarker_All = 0; iMarker_All < nMarker_All; iMarker_All++) {
    Marker_All_TagBound[iMarker_All]   = "SEND_RECEIVE";
    Marker_All_SendRecv[iMarker_All]   = 0;
    Marker_All_KindBC[iMarker_All]     = 0;
    Marker_All_Monitoring[iMarker_All] = 0;
    Marker_All_GeoEval[iMarker_All]    = 0;
    Marker_All_Designing[iMarker_All]  = 0;
    Marker_All_Plotting[iMarker_All]   = 0;
    Marker_All_DV[iMarker_All]         = 0;
    Marker_All_Moving[iMarker_All]     = 0;
    Marker_All_PerBound[iMarker_All]   = 0;
    Marker_All_Out_1D[iMarker_All]     = 0;
  }

  nMarker_Config = nMarker_Euler + nMarker_FarField + nMarker_SymWall +
  nMarker_Pressure + nMarker_PerBound + nMarker_NearFieldBound +
  nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Inlet + nMarker_Riemann +
  nMarker_Outlet + nMarker_Isothermal + nMarker_IsothermalNonCatalytic +
  nMarker_IsothermalCatalytic + nMarker_HeatFlux + nMarker_HeatFluxNonCatalytic +
  nMarker_HeatFluxCatalytic + nMarker_EngineInflow + nMarker_EngineBleed + nMarker_EngineExhaust +
  nMarker_Supersonic_Inlet + nMarker_Displacement + nMarker_Load +
  nMarker_FlowLoad + nMarker_Custom +
  nMarker_ActDisk_Inlet + nMarker_ActDisk_Outlet + nMarker_Out_1D;

  Marker_CfgFile_TagBound        = new string[nMarker_Config];
  Marker_CfgFile_KindBC   = new unsigned short[nMarker_Config];
  Marker_CfgFile_Monitoring = new unsigned short[nMarker_Config];
  Marker_CfgFile_GeoEval    = new unsigned short[nMarker_Config];
  Marker_CfgFile_Plotting   = new unsigned short[nMarker_Config];
  Marker_CfgFile_DV         = new unsigned short[nMarker_Config];
  Marker_CfgFile_Moving     = new unsigned short[nMarker_Config];
  Marker_CfgFile_Designing  = new unsigned short[nMarker_Config];
  Marker_CfgFile_PerBound   = new unsigned short[nMarker_Config];
  Marker_CfgFile_Out_1D   = new unsigned short[nMarker_Config];

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_TagBound[iMarker_Config] = "SEND_RECEIVE";
    Marker_CfgFile_KindBC[iMarker_Config]   = 0;
    Marker_CfgFile_Monitoring[iMarker_Config] = 0;
    Marker_CfgFile_GeoEval[iMarker_Config]    = 0;
    Marker_CfgFile_Designing[iMarker_Config]  = 0;
    Marker_CfgFile_Plotting[iMarker_Config]   = 0;
    Marker_CfgFile_DV[iMarker_Config]         = 0;
    Marker_CfgFile_Moving[iMarker_Config]     = 0;
    Marker_CfgFile_PerBound[iMarker_Config]   = 0;
    Marker_CfgFile_Out_1D[iMarker_Config]   = 0;
  }

  iMarker_Config = 0;
  for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Euler[iMarker_Euler];
    Marker_CfgFile_KindBC[iMarker_Config] = EULER_WALL;
    iMarker_Config++;
  }

  for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_FarField[iMarker_FarField];
    Marker_CfgFile_KindBC[iMarker_Config] = FAR_FIELD;
    iMarker_Config++;
  }

  for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_SymWall[iMarker_SymWall];
    Marker_CfgFile_KindBC[iMarker_Config] = SYMMETRY_PLANE;
    iMarker_Config++;
  }

  for (iMarker_Pressure = 0; iMarker_Pressure < nMarker_Pressure; iMarker_Pressure++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Pressure[iMarker_Pressure];
    Marker_CfgFile_KindBC[iMarker_Config] = PRESSURE_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_PerBound[iMarker_PerBound];
    Marker_CfgFile_KindBC[iMarker_Config] = PERIODIC_BOUNDARY;
    Marker_CfgFile_PerBound[iMarker_Config] = iMarker_PerBound + 1;
    iMarker_Config++;
  }

  for (iMarker_ActDisk_Inlet = 0; iMarker_ActDisk_Inlet < nMarker_ActDisk_Inlet; iMarker_ActDisk_Inlet++) {
		Marker_CfgFile_TagBound[iMarker_Config] = Marker_ActDisk_Inlet[iMarker_ActDisk_Inlet];
		Marker_CfgFile_KindBC[iMarker_Config] = ACTDISK_INLET;
		iMarker_Config++;
	}

  for (iMarker_ActDisk_Outlet = 0; iMarker_ActDisk_Outlet < nMarker_ActDisk_Outlet; iMarker_ActDisk_Outlet++) {
		Marker_CfgFile_TagBound[iMarker_Config] = Marker_ActDisk_Outlet[iMarker_ActDisk_Outlet];
		Marker_CfgFile_KindBC[iMarker_Config] = ACTDISK_OUTLET;
		iMarker_Config++;
	}

  for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_NearFieldBound[iMarker_NearFieldBound];
    Marker_CfgFile_KindBC[iMarker_Config] = NEARFIELD_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_InterfaceBound[iMarker_InterfaceBound];
    Marker_CfgFile_KindBC[iMarker_Config] = INTERFACE_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Dirichlet[iMarker_Dirichlet];
    Marker_CfgFile_KindBC[iMarker_Config] = DIRICHLET;
    iMarker_Config++;
  }

  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Inlet[iMarker_Inlet];
    Marker_CfgFile_KindBC[iMarker_Config] = INLET_FLOW;
    iMarker_Config++;
  }

  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Riemann[iMarker_Riemann];
    Marker_CfgFile_KindBC[iMarker_Config] = RIEMANN_BOUNDARY;
    iMarker_Config++;
  }

  Inflow_Mach = new double[nMarker_EngineInflow];
  Inflow_Pressure = new double[nMarker_EngineInflow];

  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_EngineInflow[iMarker_EngineInflow];
    Marker_CfgFile_KindBC[iMarker_Config] = ENGINE_INFLOW;
    Inflow_Mach[iMarker_EngineInflow] = 0.0;
    Inflow_Pressure[iMarker_EngineInflow] = 0.0;
    iMarker_Config++;
  }
  
  Bleed_MassFlow = new double[nMarker_EngineBleed];
  Bleed_Temperature = new double[nMarker_EngineBleed];
  Bleed_Pressure = new double[nMarker_EngineInflow];

  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_EngineBleed[iMarker_EngineBleed];
    Marker_CfgFile_KindBC[iMarker_Config] = ENGINE_BLEED;
    Bleed_MassFlow[iMarker_EngineBleed] = 0.0;
    Bleed_Temperature[iMarker_EngineBleed] = 0.0;
    Bleed_Pressure[iMarker_EngineBleed] = 0.0;
    iMarker_Config++;
  }

  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_EngineExhaust[iMarker_EngineExhaust];
    Marker_CfgFile_KindBC[iMarker_Config] = ENGINE_EXHAUST;
    iMarker_Config++;
  }

  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
    Marker_CfgFile_KindBC[iMarker_Config] = SUPERSONIC_INLET;
    iMarker_Config++;
  }

  for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Neumann[iMarker_Neumann];
    Marker_CfgFile_KindBC[iMarker_Config] = NEUMANN;
    iMarker_Config++;
  }

  for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Custom[iMarker_Custom];
    Marker_CfgFile_KindBC[iMarker_Config] = CUSTOM_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Outlet[iMarker_Outlet];
    Marker_CfgFile_KindBC[iMarker_Config] = OUTLET_FLOW;
    iMarker_Config++;
  }

  for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Isothermal[iMarker_Isothermal];
    Marker_CfgFile_KindBC[iMarker_Config] = ISOTHERMAL;
    iMarker_Config++;
  }

  for (iMarker_IsothermalCatalytic = 0; iMarker_IsothermalCatalytic < nMarker_IsothermalCatalytic; iMarker_IsothermalCatalytic++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_IsothermalCatalytic[iMarker_IsothermalCatalytic];
    Marker_CfgFile_KindBC[iMarker_Config] = ISOTHERMAL_CATALYTIC;
    iMarker_Config++;
  }

  for (iMarker_IsothermalNonCatalytic = 0; iMarker_IsothermalNonCatalytic < nMarker_IsothermalNonCatalytic; iMarker_IsothermalNonCatalytic++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_IsothermalNonCatalytic[iMarker_IsothermalNonCatalytic];
    Marker_CfgFile_KindBC[iMarker_Config] = ISOTHERMAL_NONCATALYTIC;
    iMarker_Config++;
  }

  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_HeatFlux[iMarker_HeatFlux];
    Marker_CfgFile_KindBC[iMarker_Config] = HEAT_FLUX;
    iMarker_Config++;
  }

  for (iMarker_HeatFluxCatalytic = 0; iMarker_HeatFluxCatalytic < nMarker_HeatFluxCatalytic; iMarker_HeatFluxCatalytic++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_HeatFluxCatalytic[iMarker_HeatFluxCatalytic];
    Marker_CfgFile_KindBC[iMarker_Config] = HEAT_FLUX_CATALYTIC;
    iMarker_Config++;
  }

  for (iMarker_HeatFluxNoncatalytic = 0; iMarker_HeatFluxNoncatalytic < nMarker_HeatFluxNonCatalytic; iMarker_HeatFluxNoncatalytic++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_HeatFluxNonCatalytic[iMarker_HeatFluxNoncatalytic];
    Marker_CfgFile_KindBC[iMarker_Config] = HEAT_FLUX_NONCATALYTIC;
    iMarker_Config++;
  }

  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Displacement[iMarker_Displacement];
    Marker_CfgFile_KindBC[iMarker_Config] = DISPLACEMENT_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_Load[iMarker_Load];
    Marker_CfgFile_KindBC[iMarker_Config] = LOAD_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
    Marker_CfgFile_TagBound[iMarker_Config] = Marker_FlowLoad[iMarker_FlowLoad];
    Marker_CfgFile_KindBC[iMarker_Config] = FLOWLOAD_BOUNDARY;
    iMarker_Config++;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_Monitoring[iMarker_Config] = NO;
    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_Monitoring[iMarker_Monitoring])
        Marker_CfgFile_Monitoring[iMarker_Config] = YES;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_GeoEval[iMarker_Config] = NO;
    for (iMarker_GeoEval = 0; iMarker_GeoEval < nMarker_GeoEval; iMarker_GeoEval++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_GeoEval[iMarker_GeoEval])
        Marker_CfgFile_GeoEval[iMarker_Config] = YES;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_Designing[iMarker_Config] = NO;
    for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_Designing[iMarker_Designing])
        Marker_CfgFile_Designing[iMarker_Config] = YES;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_Plotting[iMarker_Config] = NO;
    for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_Plotting[iMarker_Plotting])
        Marker_CfgFile_Plotting[iMarker_Config] = YES;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_DV[iMarker_Config] = NO;
    for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_DV[iMarker_DV])
        Marker_CfgFile_DV[iMarker_Config] = YES;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_Moving[iMarker_Config] = NO;
    for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_Moving[iMarker_Moving])
        Marker_CfgFile_Moving[iMarker_Config] = YES;
  }

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_CfgFile_Out_1D[iMarker_Config] = NO;
    for (iMarker_Out_1D = 0; iMarker_Out_1D < nMarker_Out_1D; iMarker_Out_1D++)
      if (Marker_CfgFile_TagBound[iMarker_Config] == Marker_Out_1D[iMarker_Out_1D])
        Marker_CfgFile_Out_1D[iMarker_Config] = YES;
  }


}

void CConfig::SetOutput(unsigned short val_software, unsigned short val_izone) {

  unsigned short iMarker_Euler, iMarker_Custom, iMarker_FarField,
  iMarker_SymWall, iMarker_PerBound, iMarker_Pressure, iMarker_NearFieldBound,
  iMarker_InterfaceBound, iMarker_Dirichlet, iMarker_Inlet,iMarker_Riemann, iMarker_Outlet,
  iMarker_Isothermal, iMarker_IsothermalNonCatalytic, iMarker_IsothermalCatalytic,
  iMarker_HeatFlux, iMarker_HeatFluxNonCatalytic, iMarker_HeatFluxCatalytic,
  iMarker_EngineInflow, iMarker_EngineBleed, iMarker_EngineExhaust, iMarker_Displacement,
  iMarker_Load, iMarker_FlowLoad,  iMarker_Neumann, iMarker_Monitoring,
  iMarker_Designing, iMarker_GeoEval, iMarker_Plotting, iMarker_DV,
  iMarker_Moving, iMarker_Supersonic_Inlet, iMarker_ActDisk_Inlet,
  iMarker_ActDisk_Outlet;

  cout << endl <<"-------------------------------------------------------------------------" << endl;
  cout <<"|    _____   _    _   ___                                               |" << endl;
  cout <<"|   / ____| | |  | | |__ \\    Release 3.2.7 \"eagle\"                     |" << endl;
  cout <<"|  | (___   | |  | |    ) |                                             |" << endl;
  cout <<"|   \\___ \\  | |  | |   / /                                              |" << endl;
  cout <<"|   ____) | | |__| |  / /_                                              |" << endl;
  switch (val_software) {
    case SU2_CFD: cout << "|  |_____/   \\____/  |____|   Suite (Computational Fluid Dynamics Code) |" << endl; break;
    case SU2_DEF: cout << "|  |_____/   \\____/  |____|   Suite (Mesh Deformation Code)             |" << endl; break;
    case SU2_DOT: cout << "|  |_____/   \\____/  |____|   Suite (Gradient Projection Code)          |" << endl; break;
    case SU2_MSH: cout << "|  |_____/   \\____/  |____|   Suite (Mesh Adaptation Code)              |" << endl; break;
    case SU2_GEO: cout << "|  |_____/   \\____/  |____|   Suite (Geometry Definition Code)          |" << endl; break;
    case SU2_SOL: cout << "|  |_____/   \\____/  |____|   Suite (Solution Exporting Code)           |" << endl; break;
  }

  cout << "|                                                                       |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << "| Copyright (C) 2012-2014 SU2 Core Developers.                          |" << endl;
  cout << "| SU2 is distributed in the hope that it will be useful,                |" << endl;
  cout << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |" << endl;
  cout << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |" << endl;
  cout << "| Lesser General Public License (version 2.1) for more details.         |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;

  cout << endl <<"------------------------ Physical case definition -----------------------" << endl;
  if (val_software == SU2_CFD) {
    switch (Kind_Solver) {
      case EULER:
        if (Kind_Regime == COMPRESSIBLE) cout << "Compressible Euler equations." << endl;
        if (Kind_Regime == INCOMPRESSIBLE) cout << "Incompressible Euler equations." << endl;
        if (Kind_Regime == FREESURFACE) {
          cout << "Incompressible Euler equations with FreeSurface." << endl;
          cout << "Free surface flow equation. Density ratio: " << RatioDensity << "." << endl;
          cout << "The free surface is located at: " << FreeSurface_Zero <<", and its thickness is: " << FreeSurface_Thickness << "." << endl;
        }
        break;
      case NAVIER_STOKES:
        if (Kind_Regime == COMPRESSIBLE) cout << "Compressible Laminar Navier-Stokes' equations." << endl;
        if (Kind_Regime == INCOMPRESSIBLE) cout << "Incompressible Laminar Navier-Stokes' equations." << endl;
        if (Kind_Regime == FREESURFACE) {
          cout << "Incompressible Laminar Navier-Stokes' equations with FreeSurface." << endl;
          cout << "Free surface flow equation. Density ratio: " << RatioDensity <<". Viscosity ratio: "<< RatioViscosity << "." << endl;
          cout << "The free surface is located at: " << FreeSurface_Zero <<", and its thickness is: " << FreeSurface_Thickness << "." << endl;
        }
        break;
      case RANS:
        if (Kind_Regime == COMPRESSIBLE) cout << "Compressible RANS equations." << endl;
        if (Kind_Regime == INCOMPRESSIBLE) cout << "Incompressible RANS equations." << endl;
        if (Kind_Regime == FREESURFACE) {
          cout << "Incompressible RANS equations with FreeSurface." << endl;
          cout << "Free surface flow equation. Density ratio: " << RatioDensity <<". Viscosity ratio: "<< RatioViscosity << "." << endl;
          cout << "The free surface is located at: " << FreeSurface_Zero <<", and its thickness is: " << FreeSurface_Thickness << "." << endl;
        }
        cout << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case SA:  cout << "Spalart Allmaras" << endl; break;
          case SST: cout << "Menter's SST"     << endl; break;
          case ML: cout << "Machine Learning" << endl;break;
        }
        break;
      case TNE2_EULER:
        cout << "Compressible TNE2 Euler equations." << endl;
        break;
      case TNE2_NAVIER_STOKES:
        cout << "Compressible TNE2 Laminar Navier-Stokes' equations." << endl;
        break;
      case POISSON_EQUATION: cout << "Poisson equation." << endl; break;
      case WAVE_EQUATION: cout << "Wave equation." << endl; break;
      case HEAT_EQUATION: cout << "Heat equation." << endl; break;
      case LINEAR_ELASTICITY: cout << "Linear elasticity solver." << endl; break;
      case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: cout << "Fluid-structure interaction." << endl; break;
      case ADJ_EULER: cout << "Continuous Euler adjoint equations." << endl; break;
      case ADJ_NAVIER_STOKES:
        if (Frozen_Visc)
          cout << "Continuous Navier-Stokes adjoint equations with frozen (laminar) viscosity." << endl;
        else
          cout << "Continuous Navier-Stokes adjoint equations." << endl;
        break;
      case ADJ_TNE2_EULER: cout << "Continuous TNE2 Euler adjoint equations." << endl; break;
      case ADJ_TNE2_NAVIER_STOKES:
        cout << "Continuous TNE2 Navier-Stokes adjoint equations with frozen (laminar) viscosity." << endl;
        break;
      case ADJ_RANS:
        if (Frozen_Visc)
          cout << "Continuous RANS adjoint equations with frozen (laminar and eddy) viscosity." << endl;
        else
          cout << "Continuous RANS adjoint equations." << endl;

        break;
      case LIN_EULER: cout << "Linearized Euler equations." << endl; break;

    }

    if ((Kind_Regime == COMPRESSIBLE) && (Kind_Solver != LINEAR_ELASTICITY) &&
        (Kind_Solver != HEAT_EQUATION) && (Kind_Solver != WAVE_EQUATION)) {
      cout << "Mach number: " << Mach <<"."<< endl;
      cout << "Angle of attack (AoA): " << AoA <<" deg, and angle of sideslip (AoS): " << AoS <<" deg."<< endl;
      if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
          (Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS))
        cout << "Reynolds number: " << Reynolds <<"."<< endl;
    }

    if (EquivArea) {
      cout <<"The equivalent area is going to be evaluated on the near-field."<< endl;
      cout <<"The lower integration limit is "<<EA_IntLimit[0]<<", and the upper is "<<EA_IntLimit[1]<<"."<< endl;
      cout <<"The near-field is situated at "<<EA_IntLimit[2]<<"."<< endl;
    }

    if (Grid_Movement) {
      cout << "Performing a dynamic mesh simulation: ";
      switch (Kind_GridMovement[ZONE_0]) {
        case NO_MOVEMENT:     cout << "no movement." << endl; break;
        case DEFORMING:       cout << "deforming mesh motion." << endl; break;
        case RIGID_MOTION:    cout << "rigid mesh motion." << endl; break;
        case MOVING_WALL:     cout << "moving walls." << endl; break;
        case ROTATING_FRAME:  cout << "rotating reference frame." << endl; break;
        case AEROELASTIC:     cout << "aeroelastic motion." << endl; break;
        case FLUID_STRUCTURE: cout << "fluid-structure motion." << endl; break;
        case EXTERNAL:        cout << "externally prescribed motion." << endl; break;
        case AEROELASTIC_RIGID_MOTION:  cout << "rigid mesh motion plus aeroelastic motion." << endl; break;
      }
    }

    if (Restart) {
      if (!Adjoint && !Linearized) cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
      if (Adjoint) cout << "Read adjoint solution from: " << Solution_AdjFileName << "." << endl;
      if (Linearized) cout << "Read linearized solution from: " << Solution_LinFileName << "." << endl;
    }
    else {
      cout << "No restart solution, use the values at infinity (freestream)." << endl;
    }

    if (Adjoint || Linearized)
      cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;

    if (RefAreaCoeff == 0) cout << "The reference length/area will be computed using y(2D) or z(3D) projection." << endl;
    else cout << "The reference length/area (force coefficient) is " << RefAreaCoeff << "." << endl;
    cout << "The reference length (moment computation) is " << RefLengthMoment << "." << endl;

    if ((nRefOriginMoment_X > 1) || (nRefOriginMoment_Y > 1) || (nRefOriginMoment_Z > 1)) {
      cout << "Surface(s) where the force coefficients are evaluated and their reference origin for moment computation: ";
      for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
        cout << Marker_Monitoring[iMarker_Monitoring] << " (" << RefOriginMoment_X[iMarker_Monitoring] <<", "<<RefOriginMoment_Y[iMarker_Monitoring] <<", "<< RefOriginMoment_Z[iMarker_Monitoring] << ")";
        if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
        else cout <<"."<< endl;
      }
    }
    else {
      cout << "Reference origin (moment computation) is (" << RefOriginMoment_X[0] << ", " << RefOriginMoment_Y[0] << ", " << RefOriginMoment_Z[0] << ")." << endl;
      cout << "Surface(s) where the force coefficients are evaluated: ";
      for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
        cout << Marker_Monitoring[iMarker_Monitoring];
        if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
        else cout <<"."<< endl;
      }
    }

    if (nMarker_Designing != 0) {
      cout << "Surface(s) where the objective function is evaluated: ";
      for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++) {
        cout << Marker_Designing[iMarker_Designing];
        if (iMarker_Designing < nMarker_Designing-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

    cout << "Surface(s) plotted in the output file: ";
    for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++) {
      cout << Marker_Plotting[iMarker_Plotting];
      if (iMarker_Plotting < nMarker_Plotting-1) cout << ", ";
      else cout <<".";
    }
    cout<< endl;

    if (nMarker_DV != 0) {
      cout << "Surface(s) affected by the design variables: ";
      for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
        cout << Marker_DV[iMarker_DV];
        if (iMarker_DV < nMarker_DV-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

    if ((Kind_GridMovement[ZONE_0] == DEFORMING) || (Kind_GridMovement[ZONE_0] == MOVING_WALL)) {
      cout << "Surface(s) in motion: ";
      for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
        cout << Marker_Moving[iMarker_Moving];
        if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

  }

  if (val_software == SU2_GEO) {
    if (nMarker_GeoEval != 0) {
      cout << "Surface(s) where the geometrical based functions is evaluated: ";
      for (iMarker_GeoEval = 0; iMarker_GeoEval < nMarker_GeoEval; iMarker_GeoEval++) {
        cout << Marker_GeoEval[iMarker_GeoEval];
        if (iMarker_GeoEval < nMarker_GeoEval-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }
  }

  cout << "Input mesh file name: " << Mesh_FileName << endl;

	if (val_software == SU2_DOT) {
		cout << "Input sensitivity file name: " << SurfAdjCoeff_FileName << "." << endl;
	}

	if (val_software == SU2_MSH) {
		switch (Kind_Adaptation) {
		case FULL: case WAKE: case FULL_FLOW: case FULL_ADJOINT: case FULL_LINEAR: case SMOOTHING: case SUPERSONIC_SHOCK:
			break;
		case GRAD_FLOW:
			cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
			break;
		case GRAD_ADJOINT:
			cout << "Read adjoint flow solution from: " << Solution_AdjFileName << "." << endl;
			break;
		case GRAD_FLOW_ADJ: case ROBUST: case COMPUTABLE_ROBUST: case COMPUTABLE: case REMAINING:
			cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
			cout << "Read adjoint flow solution from: " << Solution_AdjFileName << "." << endl;
			break;
		}
	}

	if (val_software == SU2_DEF) {
		cout << endl <<"---------------------- Grid deformation parameters ----------------------" << endl;
		cout << "Grid deformation using a linear elasticity method." << endl;

    if (Hold_GridFixed == YES) cout << "Hold some regions of the mesh fixed (hardcode implementation)." << endl;
  }

  if (val_software == SU2_DOT) {
  cout << endl <<"-------------------- Surface deformation parameters ---------------------" << endl;
  }

  if ((val_software == SU2_DEF) || (val_software == SU2_DOT)) {

    cout << "Design variables definition (markers <-> value <-> param):" << endl;

    for (unsigned short iDV = 0; iDV < nDV; iDV++) {

      switch (Design_Variable[iDV]) {
        case FFD_SETTING:           cout << "Setting the FFD box structure." ; break;
        case FFD_CONTROL_POINT_2D:  cout << "FFD 2D (control point) <-> "; break;
        case FFD_RADIUS_2D:        cout << "FFD 2D (radious)"; break;
        case FFD_CAMBER_2D:         cout << "FFD 2D (camber) <-> "; break;
        case FFD_THICKNESS_2D:      cout << "FFD 2D (thickness) <-> "; break;
        case HICKS_HENNE:           cout << "Hicks Henne <-> " ; break;
        case COSINE_BUMP:           cout << "Cosine bump <-> " ; break;
        case FOURIER:               cout << "Fourier <-> " ; break;
        case SPHERICAL:             cout << "Spherical design <-> " ; break;
        case DISPLACEMENT:          cout << "Displacement design variable."; break;
        case NACA_4DIGITS:          cout << "NACA four digits <-> "; break;
        case PARABOLIC:             cout << "Parabolic <-> "; break;
        case OBSTACLE:              cout << "Obstacle <-> "; break;
        case AIRFOIL:               cout << "Airfoil <-> "; break;
        case STRETCH:               cout << "Stretch <-> "; break;
        case ROTATION:              cout << "Rotation <-> "; break;
        case FFD_CONTROL_POINT:     cout << "FFD (control point) <-> "; break;
        case FFD_DIHEDRAL_ANGLE:    cout << "FFD (dihedral angle) <-> "; break;
        case FFD_TWIST_ANGLE:       cout << "FFD (twist angle) <-> "; break;
        case FFD_ROTATION:          cout << "FFD (rotation) <-> "; break;
        case FFD_CONTROL_SURFACE:   cout << "FFD (control surface) <-> "; break;
        case FFD_CAMBER:            cout << "FFD (camber) <-> "; break;
        case FFD_THICKNESS:         cout << "FFD (thickness) <-> "; break;
        case SURFACE_FILE:          cout << "Surface file based deformation." ; break;
      }

      if ((Design_Variable[iDV] != FFD_SETTING) &&
          (Design_Variable[iDV] != SURFACE_FILE)) {

        for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
          cout << Marker_DV[iMarker_DV];
          if (iMarker_DV < nMarker_DV-1) cout << ", ";
          else cout << " <-> ";
        }
        cout << DV_Value[iDV] << " <-> ";

        if (Design_Variable[iDV] == FFD_SETTING) nParamDV = 0;
        if (Design_Variable[iDV] == SURFACE_FILE) nParamDV = 0;
        if (Design_Variable[iDV] == FFD_CONTROL_POINT_2D) nParamDV = 5;
        if (Design_Variable[iDV] == FFD_RADIUS_2D) nParamDV = 1;
        if (Design_Variable[iDV] == FFD_CAMBER_2D) nParamDV = 2;
        if (Design_Variable[iDV] == FFD_THICKNESS_2D) nParamDV = 2;
        if (Design_Variable[iDV] == HICKS_HENNE) nParamDV = 2;
        if (Design_Variable[iDV] == SPHERICAL) nParamDV = 3;
        if (Design_Variable[iDV] == COSINE_BUMP) nParamDV = 3;
        if (Design_Variable[iDV] == FOURIER) nParamDV = 3;
        if (Design_Variable[iDV] == DISPLACEMENT) nParamDV = 3;
        if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
        if (Design_Variable[iDV] == NACA_4DIGITS) nParamDV = 3;
        if (Design_Variable[iDV] == PARABOLIC) nParamDV = 2;
        if (Design_Variable[iDV] == OBSTACLE) nParamDV = 2;
        if (Design_Variable[iDV] == AIRFOIL) nParamDV = 2;
        if (Design_Variable[iDV] == STRETCH) nParamDV = 2;
        if (Design_Variable[iDV] == FFD_CONTROL_POINT) nParamDV = 7;
        if (Design_Variable[iDV] == FFD_DIHEDRAL_ANGLE) nParamDV = 7;
        if (Design_Variable[iDV] == FFD_TWIST_ANGLE) nParamDV = 7;
        if (Design_Variable[iDV] == FFD_ROTATION) nParamDV = 7;
        if (Design_Variable[iDV] == FFD_CONTROL_SURFACE) nParamDV = 7;
        if (Design_Variable[iDV] == FFD_CAMBER) nParamDV = 3;
        if (Design_Variable[iDV] == FFD_THICKNESS) nParamDV = 3;

        for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {

          if (iParamDV == 0) cout << "( ";

          if ((iParamDV == 0) &&
              ((Design_Variable[iDV] == FFD_SETTING) ||
               (Design_Variable[iDV] == FFD_CONTROL_POINT_2D) ||
               (Design_Variable[iDV] == FFD_CAMBER_2D) ||
               (Design_Variable[iDV] == FFD_THICKNESS_2D) ||
               (Design_Variable[iDV] == FFD_RADIUS_2D) ||
               (Design_Variable[iDV] == FFD_CONTROL_POINT) ||
               (Design_Variable[iDV] == FFD_DIHEDRAL_ANGLE) ||
               (Design_Variable[iDV] == FFD_TWIST_ANGLE) ||
               (Design_Variable[iDV] == FFD_ROTATION) ||
               (Design_Variable[iDV] == FFD_CONTROL_SURFACE) ||
               (Design_Variable[iDV] == FFD_CAMBER) ||
               (Design_Variable[iDV] == FFD_THICKNESS))) cout << FFDTag[iDV];
          else cout << ParamDV[iDV][iParamDV];

          if (iParamDV < nParamDV-1) cout << ", ";
          else cout <<" )"<< endl;;
        }

      } else cout << endl;

		}
	}

	if (((val_software == SU2_CFD) && ( Adjoint )) || (val_software == SU2_DOT)) {

		cout << endl <<"----------------------- Design problem definition -----------------------" << endl;
		switch (Kind_ObjFunc) {
      case DRAG_COEFFICIENT: cout << "Drag objective function." << endl; break;
      case LIFT_COEFFICIENT: cout << "Lift objective function." << endl; break;
      case INVERSE_DESIGN_PRESSURE: cout << "Inverse design (Cp) objective function." << endl; break;
      case INVERSE_DESIGN_HEATFLUX: cout << "Inverse design (Heat Flux) objective function." << endl; break;
      case SIDEFORCE_COEFFICIENT: cout << "Side force objective function." << endl; break;
      case MOMENT_X_COEFFICIENT: cout << "Mx objective function." << endl; break;
      case MOMENT_Y_COEFFICIENT: cout << "My objective function." << endl; break;
      case MOMENT_Z_COEFFICIENT: cout << "Mz objective function." << endl; break;
      case EFFICIENCY: cout << "Efficiency objective function." << endl; break;
      case EQUIVALENT_AREA:
        cout << "Equivalent area objective function." << endl;
        cout << "Drag coefficient weight in the objective function: " << WeightCd <<"."<< endl;  break;
      case NEARFIELD_PRESSURE:
        cout << "Nearfield pressure objective function." << endl;
        cout << "Drag coefficient weight in the objective function: " << WeightCd <<"."<< endl;  break;
        break;
      case FORCE_X_COEFFICIENT: cout << "X-force objective function." << endl; break;
      case FORCE_Y_COEFFICIENT: cout << "Y-force moment objective function." << endl; break;
      case FORCE_Z_COEFFICIENT: cout << "Z-force moment objective function." << endl; break;
      case THRUST_COEFFICIENT: cout << "Thrust objective function." << endl; break;
      case TORQUE_COEFFICIENT: cout << "Torque efficiency objective function." << endl; break;
      case TOTAL_HEATFLUX: cout << "Total heat flux objective function." << endl; break;
      case MAXIMUM_HEATFLUX: cout << "Maximum heat flux objective function." << endl; break;
      case FIGURE_OF_MERIT: cout << "Rotor Figure of Merit objective function." << endl; break;
      case FREE_SURFACE: cout << "Free-Surface objective function." << endl; break;
      case AVG_TOTAL_PRESSURE: cout << "Average total objective pressure." << endl; break;
      case MASS_FLOW_RATE: cout << "Mass flow rate objective function." << endl; break;
		}

	}

	if (val_software == SU2_CFD) {
		cout << endl <<"---------------------- Space numerical integration ----------------------" << endl;

		if (SmoothNumGrid) cout << "There are some smoothing iterations on the grid coordinates." << endl;

		if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {

      if (Kind_ConvNumScheme_Flow == SPACE_CENTERED) {
        if (Kind_Centered_Flow == JST) {
          cout << "Jameson-Schmidt-Turkel scheme for the flow inviscid terms."<< endl;
          cout << "JST viscous coefficients (1st, 2nd & 4th): " << Kappa_1st_Flow
          << ", " << Kappa_2nd_Flow << ", " << Kappa_4th_Flow <<"."<< endl;
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
          cout << "Second order integration." << endl;
        }
        if (Kind_Centered_Flow == JST_KE) {
          cout << "Jameson-Schmidt-Turkel scheme for the flow inviscid terms."<< endl;
          cout << "JST viscous coefficients (1st, 2nd): " << Kappa_1st_Flow
          << ", " << Kappa_2nd_Flow << "."<< endl;
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
          cout << "Second order integration." << endl;
        }
        if (Kind_Centered_Flow == LAX) {
          cout << "Lax-Friedrich scheme for the flow inviscid terms."<< endl;
          cout << "First order integration." << endl;
        }
      }

			if (Kind_ConvNumScheme_Flow == SPACE_UPWIND) {
				if (Kind_Upwind_Flow == ROE) cout << "Roe solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == TURKEL) cout << "Roe-Turkel solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == AUSM)	cout << "AUSM solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == HLLC)	cout << "HLLC solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == SW)	cout << "Steger-Warming solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == MSW)	cout << "Modified Steger-Warming solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == CUSP)	cout << "CUSP solver for the flow inviscid terms."<< endl;
        switch (SpatialOrder_Flow) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (Kind_SlopeLimit_Flow) {
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
			}

		}

    if (Kind_Solver == RANS) {
      if (Kind_ConvNumScheme_Turb == SPACE_UPWIND) {
        if (Kind_Upwind_Turb == SCALAR_UPWIND) cout << "Scalar upwind solver (first order) for the turbulence model."<< endl;
        switch (SpatialOrder_Turb) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (Kind_SlopeLimit_Turb) {
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
      }
    }

    if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {

      if (Kind_ConvNumScheme_TNE2 == SPACE_CENTERED) {
        if (Kind_ConvNumScheme_TNE2 == LAX) {
          cout << "Lax-Friedrich convective scheme for the inviscid terms of the two-temperature model." << endl;
          cout << "First order integration." << endl;
        }
      }
      if (Kind_ConvNumScheme_TNE2 == SPACE_UPWIND) {
        if (Kind_Upwind_TNE2 == ROE) cout << "Roe solver for the inviscid terms of the two-temperature model."<< endl;
        if (Kind_Upwind_TNE2 == TURKEL) cout << "Roe-Turkel solver for the inviscid terms of the two-temperature model."<< endl;
        if (Kind_Upwind_TNE2 == AUSM)	cout << "AUSM solver for the inviscid terms of the two-temperature model."<< endl;
        if (Kind_Upwind_TNE2 == HLLC)	cout << "HLLC solver for the inviscid terms of the two-temperature model."<< endl;
        if (Kind_Upwind_TNE2 == SW)	cout << "Steger-Warming solver for the inviscid terms of the two-temperature model."<< endl;
        if (Kind_Upwind_TNE2 == MSW)	cout << "Modified Steger-Warming solver for the inviscid terms of the two-temperature model."<< endl;
        switch (SpatialOrder_TNE2) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (Kind_SlopeLimit_TNE2) {
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
      }

    }

    if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {

      if (Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED) {
        if (Kind_Centered_AdjFlow == JST) {
          cout << "Jameson-Schmidt-Turkel scheme for the adjoint inviscid terms."<< endl;
          cout << "JST viscous coefficients (1st, 2nd, & 4th): " << Kappa_1st_AdjFlow
          << ", " << Kappa_2nd_AdjFlow << ", " << Kappa_4th_AdjFlow <<"."<< endl;
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
          cout << "Second order integration." << endl;
        }
        if (Kind_Centered_AdjFlow == LAX) {
          cout << "Lax-Friedrich scheme for the adjoint inviscid terms."<< endl;
          cout << "First order integration." << endl;
        }
      }

      if (Kind_ConvNumScheme_AdjFlow == SPACE_UPWIND) {
        if (Kind_Upwind_AdjFlow == ROE) cout << "Roe solver for the adjoint inviscid terms."<< endl;
        switch (SpatialOrder_AdjFlow) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (Kind_SlopeLimit_AdjFlow) {
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case SHARP_EDGES:
                cout << "Sharp edges slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                cout << "The reference sharp edge distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
      }
      
      cout << "The reference sharp edge distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;

    }

    if ((Kind_Solver == ADJ_RANS) && (!Frozen_Visc)) {
      if (Kind_ConvNumScheme_AdjTurb == SPACE_UPWIND) {
        if (Kind_Upwind_Turb == SCALAR_UPWIND) cout << "Scalar upwind solver (first order) for the adjoint turbulence model."<< endl;
        switch (SpatialOrder_AdjTurb) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (Kind_SlopeLimit_AdjTurb) {
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case SHARP_EDGES:
                cout << "Sharp edges slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                cout << "The reference sharp edge distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
      }
    }

    if (Kind_Solver == ADJ_TNE2_EULER || Kind_Solver == ADJ_TNE2_NAVIER_STOKES) {

      if (Kind_ConvNumScheme_AdjTNE2 == SPACE_CENTERED) {
        if (Kind_Centered_AdjTNE2 == JST) {
          cout << "Jameson-Schmidt-Turkel scheme for the adjoint inviscid terms."<< endl;
          cout << "JST viscous coefficients (1st, 2nd, & 4th): " << Kappa_1st_AdjTNE2
          << ", " << Kappa_2nd_AdjTNE2 << ", " << Kappa_4th_AdjTNE2 <<"."<< endl;
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
          cout << "The reference sharp edge distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;
          cout << "Second order integration." << endl;
        }
        if (Kind_Centered_AdjTNE2 == LAX) {
          cout << "Lax-Friedrich scheme for the adjoint inviscid terms."<< endl;
          cout << "First order integration." << endl;
        }
      }

      if (Kind_ConvNumScheme_AdjTNE2 == SPACE_UPWIND) {
        if (Kind_Upwind_AdjTNE2 == ROE) cout << "Roe solver for the adjoint inviscid terms."<< endl;

        switch (SpatialOrder_AdjTNE2) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (SpatialOrder_AdjTNE2) {
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case SHARP_EDGES:
                cout << "Sharp edges slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                cout << "The reference sharp edge distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
      }

    }

    if (Kind_Solver == LIN_EULER) {

      if (Kind_ConvNumScheme_LinFlow == SPACE_CENTERED) {
        if (Kind_Centered_LinFlow == JST) {
          cout << "Jameson-Schmidt-Turkel scheme for the linearized inviscid terms."<< endl;
          cout << "JST viscous coefficients (1st, & 4th): " << Kappa_1st_LinFlow
          << ", " << Kappa_4th_LinFlow <<"."<< endl;
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
          cout << "Second order integration." << endl;
        }
        if (Kind_Centered_LinFlow == LAX) {
          cout << "Lax-Friedrich scheme for the linearized inviscid terms."<< endl;
          cout << "First order integration." << endl;
        }
      }

    }

    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        cout << "Average of gradients with correction (viscous flow terms)." << endl;
    }

    if (Kind_Solver == TNE2_NAVIER_STOKES) {
        cout << "Average of gradients with correction (viscous flow terms)." << endl;
    }

    if ((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
      cout << "Average of gradients with correction (viscous adjoint terms)." << endl;
    }

    if (Kind_Solver == RANS) {
      cout << "Average of gradients with correction (viscous turbulence terms)." << endl;
    }

    if (Kind_Solver == POISSON_EQUATION) {
      cout << "Galerkin method for viscous terms computation of the poisson potential equation." << endl;
    }

    if ((Kind_Solver == ADJ_RANS) && (!Frozen_Visc)) {
      cout << "Average of gradients with correction (2nd order) for computation of adjoint viscous turbulence terms." << endl;
      if (Kind_TimeIntScheme_AdjTurb == EULER_IMPLICIT) cout << "Euler implicit method for the turbulent adjoint equation." << endl;
    }

    switch (Kind_Gradient_Method) {
      case GREEN_GAUSS: cout << "Gradient computation using Green-Gauss theorem." << endl; break;
      case WEIGHTED_LEAST_SQUARES: cout << "Gradient Computation using weighted Least-Squares method." << endl; break;
    }

    if ((Kind_Regime == INCOMPRESSIBLE) || (Kind_Regime == FREESURFACE)) {
      cout << "Artificial compressibility factor: " << ArtComp_Factor << "." << endl;
    }

    cout << endl <<"---------------------- Time numerical integration -----------------------" << endl;
    switch (Unsteady_Simulation) {
      case NO:
        cout << "Local time stepping (steady state simulation)." << endl; break;
      case TIME_STEPPING:
        cout << "Unsteady simulation using a time stepping strategy."<< endl;
        if (Unst_CFL != 0.0) cout << "Time step computed by the code. Unsteady CFL number: " << Unst_CFL <<"."<< endl;
        else cout << "Unsteady time step provided by the user (s): "<< Delta_UnstTime << "." << endl;
        break;
      case DT_STEPPING_1ST: case DT_STEPPING_2ND:
        if (Unsteady_Simulation == DT_STEPPING_1ST) cout << "Unsteady simulation, dual time stepping strategy (first order in time)."<< endl;
        if (Unsteady_Simulation == DT_STEPPING_2ND) cout << "Unsteady simulation, dual time stepping strategy (second order in time)."<< endl;
        if (Unst_CFL != 0.0) cout << "Time step computed by the code. Unsteady CFL number: " << Unst_CFL <<"."<< endl;
        else cout << "Unsteady time step provided by the user (s): "<< Delta_UnstTime << "." << endl;
        cout << "Total number of internal Dual Time iterations: "<< Unst_nIntIter <<"." << endl;
        break;
    }

    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      switch (Kind_TimeIntScheme_Flow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the flow equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT: cout << "Euler explicit method for the flow equations." << endl; break;
        case EULER_IMPLICIT:
          cout << "Euler implicit method for the flow equations." << endl;
          switch (Kind_Linear_Solver) {
            case BCGSTAB:
              cout << "BCGSTAB is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
            case FGMRES || RFGMRES:
              cout << "FGMRES is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
            case SMOOTHER_JACOBI:
              cout << "A Jacobi method is used for smoothing the linear system." << endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
            case SMOOTHER_ILU:
              cout << "A ILU0 method is used for smoothing the linear system." << endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
            case SMOOTHER_LUSGS:
              cout << "A LU-SGS method is used for smoothing the linear system." << endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
            case SMOOTHER_LINELET:
              cout << "A Linelet method is used for smoothing the linear system." << endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
          }
          break;
      }
    }

    if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
      switch (Kind_TimeIntScheme_TNE2) {
        case EULER_IMPLICIT:
          cout << "Euler implicit method for the flow equations." << endl;
          switch (Kind_Linear_Solver) {
            case BCGSTAB:
              cout << "BCGSTAB is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
            case FGMRES || RFGMRES:
              cout << "FGMRES is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<< endl;
              break;
          }
          break;
      }
    }

    if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
      switch (Kind_TimeIntScheme_AdjFlow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the adjoint equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT: cout << "Euler explicit method for the adjoint equations." << endl; break;
        case EULER_IMPLICIT: cout << "Euler implicit method for the adjoint equations." << endl; break;
      }
    }

    if (Kind_Solver == LIN_EULER) {
      switch (Kind_TimeIntScheme_LinFlow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit explicit method for the linearized equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT: cout << "Euler explicit method for the linearized equations." << endl; break;
      }
    }

    if ((Kind_Solver == ADJ_TNE2_EULER) || (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)) {
      switch (Kind_TimeIntScheme_AdjTNE2) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the adjoint equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT: cout << "Euler explicit method for the adjoint equations." << endl; break;
        case EULER_IMPLICIT: cout << "Euler implicit method for the adjoint equations." << endl; break;
      }
    }

    if (nMultiLevel !=0) {
      
      if (nStartUpIter != 0) cout << "A total of " << nStartUpIter << " start up iterations on the fine grid."<< endl;
      if (MGCycle == V_CYCLE) cout << "V Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;
      if (MGCycle == W_CYCLE) cout << "W Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;
      if (MGCycle == FULLMG_CYCLE) cout << "Full Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;

      cout << "Damping factor for the residual restriction: " << Damp_Res_Restric <<"."<< endl;
      cout << "Damping factor for the correction prolongation: " << Damp_Correc_Prolong <<"."<< endl;
    }

    if ((Kind_Solver != LINEAR_ELASTICITY) && (Kind_Solver != HEAT_EQUATION) && (Kind_Solver != WAVE_EQUATION)) {

      if (CFLRamp[0] == 1.0) cout << "No CFL ramp." << endl;
      else cout << "CFL ramp definition. factor: "<< CFLRamp[0] <<", every "<< int(CFLRamp[1]) <<" iterations, with a limit of "<< CFLRamp[2] <<"." << endl;

      if (nMultiLevel !=0) {
        cout << "Multigrid Level:                  ";
        for (unsigned short iLevel = 0; iLevel < nMultiLevel+1; iLevel++) {
          cout.width(6); cout << iLevel;
        }
        cout << endl;
      }

      cout << "Courant-Friedrichs-Lewy number:   ";
      cout.precision(3);
      cout.width(6); cout << CFL[0];
      cout << endl;

      if (nMultiLevel !=0) {
        cout.precision(3);
        cout << "MG PreSmooth coefficients:        ";
        for (unsigned short iMG_PreSmooth = 0; iMG_PreSmooth < nMultiLevel+1; iMG_PreSmooth++) {
          cout.width(6); cout << MG_PreSmooth[iMG_PreSmooth];
        }
        cout << endl;
      }

      if (nMultiLevel !=0) {
        cout.precision(3);
        cout << "MG PostSmooth coefficients:       ";
        for (unsigned short iMG_PostSmooth = 0; iMG_PostSmooth < nMultiLevel+1; iMG_PostSmooth++) {
          cout.width(6); cout << MG_PostSmooth[iMG_PostSmooth];
        }
        cout << endl;
      }

      if (nMultiLevel !=0) {
        cout.precision(3);
        cout << "MG CorrecSmooth coefficients:     ";
        for (unsigned short iMG_CorrecSmooth = 0; iMG_CorrecSmooth < nMultiLevel+1; iMG_CorrecSmooth++) {
          cout.width(6); cout << MG_CorrecSmooth[iMG_CorrecSmooth];
        }
        cout << endl;
      }

    }

    if (Kind_Solver == RANS)
      if (Kind_TimeIntScheme_Turb == EULER_IMPLICIT)
        cout << "Euler implicit time integration for the turbulence model." << endl;
  }

  if (val_software == SU2_CFD) {

    cout << endl <<"------------------------- Convergence criteria --------------------------" << endl;

    cout << "Maximum number of iterations: " << nExtIter <<"."<< endl;

    if (ConvCriteria == CAUCHY) {
      if (!Adjoint && !Linearized)
        switch (Cauchy_Func_Flow) {
          case LIFT_COEFFICIENT: cout << "Cauchy criteria for Lift using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case DRAG_COEFFICIENT: cout << "Cauchy criteria for Drag using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }

      if (Adjoint)
        switch (Cauchy_Func_AdjFlow) {
          case SENS_GEOMETRY: cout << "Cauchy criteria for geo. sensitivity using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case SENS_MACH: cout << "Cauchy criteria for Mach number sensitivity using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }

      if (Linearized)
        switch (Cauchy_Func_LinFlow) {
          case DELTA_LIFT_COEFFICIENT: cout << "Cauchy criteria for linearized Lift using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case DELTA_DRAG_COEFFICIENT: cout << "Cauchy criteria for linearized Drag using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }

      cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
      
    }


    if (ConvCriteria == RESIDUAL) {
      if (!Adjoint && !Linearized) {
        cout << "Reduce the density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum bound for the density residual is 10^(" << MinLogResidual<< ")."<< endl;
        cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
      }

      if (Adjoint) {
        cout << "Reduce the adjoint density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum value for the adjoint density residual is 10^(" << MinLogResidual<< ")."<< endl;
      }

      if (Linearized) {
        cout << "Reduce the linearized density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum value for the linearized density residual is 10^(" << MinLogResidual<< ")."<< endl;
      }

    }

  }

  if (val_software == SU2_MSH) {
    cout << endl <<"----------------------- Grid adaptation strategy ------------------------" << endl;

    switch (Kind_Adaptation) {
      case NONE: break;
      case PERIODIC: cout << "Grid modification to run periodic bc problems." << endl; break;
      case FULL: cout << "Grid adaptation using a complete refinement." << endl; break;
      case WAKE: cout << "Grid adaptation of the wake." << endl; break;
      case FULL_FLOW: cout << "Flow grid adaptation using a complete refinement." << endl; break;
      case FULL_ADJOINT: cout << "Adjoint grid adaptation using a complete refinement." << endl; break;
      case FULL_LINEAR: cout << "Linear grid adaptation using a complete refinement." << endl; break;
      case GRAD_FLOW: cout << "Grid adaptation using gradient based strategy (density)." << endl; break;
      case GRAD_ADJOINT: cout << "Grid adaptation using gradient based strategy (adjoint density)." << endl; break;
      case GRAD_FLOW_ADJ: cout << "Grid adaptation using gradient based strategy (density and adjoint density)." << endl; break;
      case ROBUST: cout << "Grid adaptation using robust adaptation."<< endl; break;
      case COMPUTABLE: cout << "Grid adaptation using computable correction."<< endl; break;
      case COMPUTABLE_ROBUST: cout << "Grid adaptation using computable correction."<< endl; break;
      case REMAINING: cout << "Grid adaptation using remaining error."<< endl; break;
      case SMOOTHING: cout << "Grid smoothing using an implicit method."<< endl; break;
      case SUPERSONIC_SHOCK: cout << "Grid adaptation for a supersonic shock at Mach: " << Mach <<"."<< endl; break;
    }

    switch (Kind_Adaptation) {
      case GRAD_FLOW: case GRAD_ADJOINT: case GRAD_FLOW_ADJ: case ROBUST: case COMPUTABLE: case COMPUTABLE_ROBUST: case REMAINING:
        cout << "Power of the dual volume in the adaptation sensor: " << DualVol_Power << endl;
        cout << "Percentage of new elements in the adaptation process: " << New_Elem_Adapt << "."<< endl;
        break;
    }

    if (Analytical_Surface != NONE)
      cout << "Use analytical definition for including points in the surfaces." << endl;

  }

  cout << endl <<"-------------------------- Output information ---------------------------" << endl;

  if (val_software == SU2_CFD) {

    if (Low_MemoryOutput) cout << "Writing output files with low memory RAM requirements."<< endl;
    cout << "Writing a flow solution every " << Wrt_Sol_Freq <<" iterations."<< endl;
    cout << "Writing the convergence history every " << Wrt_Con_Freq <<" iterations."<< endl;
    if ((Unsteady_Simulation == DT_STEPPING_1ST) || (Unsteady_Simulation == DT_STEPPING_2ND))  {
      cout << "Writing the dual time flow solution every " << Wrt_Sol_Freq_DualTime <<" iterations."<< endl;
      cout << "Writing the dual time convergence history every " << Wrt_Con_Freq_DualTime <<" iterations."<< endl;
    }

    switch (Output_FileFormat) {
      case PARAVIEW: cout << "The output file format is Paraview ASCII (.vtk)." << endl; break;
      case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
      case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
      case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
    }

    cout << "Convergence history file name: " << Conv_FileName << "." << endl;

    cout << "Forces breakdown file name: " << Breakdown_FileName << "." << endl;

    if ((Kind_Solver != LINEAR_ELASTICITY) && (Kind_Solver != HEAT_EQUATION) && (Kind_Solver != WAVE_EQUATION)) {
      if (!Linearized && !Adjoint) {
        cout << "Surface flow coefficients file name: " << SurfFlowCoeff_FileName << "." << endl;
        cout << "Flow variables file name: " << Flow_FileName << "." << endl;
        cout << "Restart flow file name: " << Restart_FlowFileName << "." << endl;
      }

      if (Linearized) {
        cout << "Linearized flow solution file name: " << Solution_LinFileName << "." << endl;
        cout << "Restart linearized flow file name: " << Restart_LinFileName << "." << endl;
        cout << "Linearized variables file name: " << Lin_FileName << "." << endl;
        cout << "Surface linearized coefficients file name: " << SurfLinCoeff_FileName << "." << endl;
      }

      if (Adjoint) {
        cout << "Adjoint solution file name: " << Solution_AdjFileName << "." << endl;
        cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
        cout << "Adjoint variables file name: " << Adj_FileName << "." << endl;
        cout << "Surface adjoint coefficients file name: " << SurfAdjCoeff_FileName << "." << endl;
      }
    }
    else {
      cout << "Surface structure coefficients file name: " << SurfStructure_FileName << "." << endl;
      cout << "Structure variables file name: " << Structure_FileName << "." << endl;
      cout << "Restart structure file name: " << Restart_FlowFileName << "." << endl;
    }

  }

  if (val_software == SU2_SOL) {
    if (Low_MemoryOutput) cout << "Writing output files with low memory RAM requirements."<< endl;
    switch (Output_FileFormat) {
      case PARAVIEW: cout << "The output file format is Paraview ASCII (.dat)." << endl; break;
      case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
      case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
      case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
    }
    cout << "Flow variables file name: " << Flow_FileName << "." << endl;
  }

  if (val_software == SU2_DEF) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
    if (Visualize_Deformation) cout << "A file will be created to visualize the deformation." << endl;
    else cout << "No file for visualizing the deformation." << endl;
    switch (GetDeform_Stiffness_Type()) {
      case INVERSE_VOLUME:
        cout << "Cell stiffness scaled by inverse of the cell volume." << endl;
        break;
      case WALL_DISTANCE:
        cout << "Cell stiffness scaled by distance from the deforming surface." << endl;
        break;
      case CONSTANT_STIFFNESS:
        cout << "Imposing constant cell stiffness (steel)." << endl;
        break;
    }
  }

  if (val_software == SU2_MSH) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
  }

  if (val_software == SU2_DOT) {
    cout << "Output gradient file name: " << ObjFunc_Grad_FileName << ". " << endl;
  }

  if (val_software == SU2_MSH) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
    cout << "Restart flow file name: " << Restart_FlowFileName << "." << endl;
    if ((Kind_Adaptation == FULL_ADJOINT) || (Kind_Adaptation == GRAD_ADJOINT) || (Kind_Adaptation == GRAD_FLOW_ADJ) ||
        (Kind_Adaptation == ROBUST) || (Kind_Adaptation == COMPUTABLE_ROBUST) || (Kind_Adaptation == COMPUTABLE) ||
        (Kind_Adaptation == REMAINING)) {
      if (Kind_ObjFunc == DRAG_COEFFICIENT) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc == EQUIVALENT_AREA) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc == NEARFIELD_PRESSURE) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc == LIFT_COEFFICIENT) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
    }
  }

  cout << endl <<"------------------- Config file boundary information --------------------" << endl;

  if (nMarker_Euler != 0) {
    cout << "Euler wall boundary marker(s): ";
    for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
      cout << Marker_Euler[iMarker_Euler];
      if (iMarker_Euler < nMarker_Euler-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_FarField != 0) {
    cout << "Far-field boundary marker(s): ";
    for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
      cout << Marker_FarField[iMarker_FarField];
      if (iMarker_FarField < nMarker_FarField-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_SymWall != 0) {
    cout << "Symmetry plane boundary marker(s): ";
    for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
      cout << Marker_SymWall[iMarker_SymWall];
      if (iMarker_SymWall < nMarker_SymWall-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Pressure != 0) {
    cout << "Pressure boundary marker(s): ";
    for (iMarker_Pressure = 0; iMarker_Pressure < nMarker_Pressure; iMarker_Pressure++) {
      cout << Marker_Pressure[iMarker_Pressure];
      if (iMarker_Pressure < nMarker_Pressure-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_PerBound != 0) {
    cout << "Periodic boundary marker(s): ";
    for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
      cout << Marker_PerBound[iMarker_PerBound];
      if (iMarker_PerBound < nMarker_PerBound-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_NearFieldBound != 0) {
    cout << "Near-field boundary marker(s): ";
    for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
      cout << Marker_NearFieldBound[iMarker_NearFieldBound];
      if (iMarker_NearFieldBound < nMarker_NearFieldBound-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_InterfaceBound != 0) {
    cout << "Interface boundary marker(s): ";
    for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
      cout << Marker_InterfaceBound[iMarker_InterfaceBound];
      if (iMarker_InterfaceBound < nMarker_InterfaceBound-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Dirichlet != 0) {
    cout << "Dirichlet boundary marker(s): ";
    for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
      cout << Marker_Dirichlet[iMarker_Dirichlet];
      if (iMarker_Dirichlet < nMarker_Dirichlet-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_FlowLoad != 0) {
    cout << "Flow Load boundary marker(s): ";
    for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
      cout << Marker_FlowLoad[iMarker_FlowLoad];
      if (iMarker_FlowLoad < nMarker_FlowLoad-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Neumann != 0) {
    cout << "Neumann boundary marker(s): ";
    for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
      cout << Marker_Neumann[iMarker_Neumann];
      if (iMarker_Neumann < nMarker_Neumann-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Inlet != 0) {
    cout << "Inlet boundary marker(s): ";
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      cout << Marker_Inlet[iMarker_Inlet];
      if (iMarker_Inlet < nMarker_Inlet-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Riemann != 0) {
      cout << "Riemann boundary marker(s): ";
      for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++) {
        cout << Marker_Riemann[iMarker_Riemann];
        if (iMarker_Riemann < nMarker_Riemann-1) cout << ", ";
        else cout <<"."<< endl;
    }
  }
  
  if (nMarker_EngineInflow != 0) {
    cout << "Engine inflow boundary marker(s): ";
    for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
      cout << Marker_EngineInflow[iMarker_EngineInflow];
      if (iMarker_EngineInflow < nMarker_EngineInflow-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }
  
  if (nMarker_EngineBleed != 0) {
    cout << "Engine bleed boundary marker(s): ";
    for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
      cout << Marker_EngineBleed[iMarker_EngineBleed];
      if (iMarker_EngineBleed < nMarker_EngineBleed-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_EngineExhaust != 0) {
    cout << "Engine exhaust boundary marker(s): ";
    for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
      cout << Marker_EngineExhaust[iMarker_EngineExhaust];
      if (iMarker_EngineExhaust < nMarker_EngineExhaust-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Supersonic_Inlet != 0) {
    cout << "Supersonic inlet boundary marker(s): ";
    for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
      cout << Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
      if (iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Outlet != 0) {
    cout << "Outlet boundary marker(s): ";
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      cout << Marker_Outlet[iMarker_Outlet];
      if (iMarker_Outlet < nMarker_Outlet-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Isothermal != 0) {
    cout << "Isothermal wall boundary marker(s): ";
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
      cout << Marker_Isothermal[iMarker_Isothermal];
      if (iMarker_Isothermal < nMarker_Isothermal-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_IsothermalNonCatalytic != 0) {
    cout << "Non-catalytic isothermal wall boundary marker(s): ";
    for (iMarker_IsothermalNonCatalytic = 0; iMarker_IsothermalNonCatalytic < nMarker_IsothermalNonCatalytic; iMarker_IsothermalNonCatalytic++) {
      cout << Marker_IsothermalNonCatalytic[iMarker_IsothermalNonCatalytic];
      if (iMarker_IsothermalNonCatalytic < nMarker_IsothermalNonCatalytic-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_IsothermalCatalytic != 0) {
    cout << "Catalytic isothermal wall boundary marker(s): ";
    for (iMarker_IsothermalCatalytic = 0; iMarker_IsothermalCatalytic < nMarker_IsothermalCatalytic; iMarker_IsothermalCatalytic++) {
      cout << Marker_IsothermalCatalytic[iMarker_IsothermalCatalytic];
      if (iMarker_IsothermalCatalytic < nMarker_IsothermalCatalytic-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_HeatFlux != 0) {
    cout << "Constant heat flux wall boundary marker(s): ";
    for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
      cout << Marker_HeatFlux[iMarker_HeatFlux];
      if (iMarker_HeatFlux < nMarker_HeatFlux-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_HeatFluxNonCatalytic != 0) {
    cout << "Non-catalytic constant heat flux wall boundary marker(s): ";
    for (iMarker_HeatFluxNonCatalytic = 0; iMarker_HeatFluxNonCatalytic < nMarker_HeatFluxNonCatalytic; iMarker_HeatFluxNonCatalytic++) {
      cout << Marker_HeatFluxNonCatalytic[iMarker_HeatFluxNonCatalytic];
      if (iMarker_HeatFluxNonCatalytic < nMarker_HeatFluxNonCatalytic-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_HeatFluxCatalytic != 0) {
    cout << "Catalytic constant heat flux wall boundary marker(s): ";
    for (iMarker_HeatFluxCatalytic = 0; iMarker_HeatFluxCatalytic < nMarker_HeatFluxCatalytic; iMarker_HeatFluxCatalytic++) {
      cout << Marker_HeatFluxCatalytic[iMarker_HeatFluxCatalytic];
      if (iMarker_HeatFluxCatalytic < nMarker_HeatFluxCatalytic-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Displacement != 0) {
    cout << "Displacement boundary marker(s): ";
    for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
      cout << Marker_Displacement[iMarker_Displacement];
      if (iMarker_Displacement < nMarker_Displacement-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Load != 0) {
    cout << "Load boundary marker(s): ";
    for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
      cout << Marker_Load[iMarker_Load];
      if (iMarker_Load < nMarker_Load-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Neumann != 0) {
    cout << "Neumann boundary marker(s): ";
    for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
      cout << Marker_Neumann[iMarker_Neumann];
      if (iMarker_Neumann < nMarker_Neumann-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Custom != 0) {
    cout << "Custom boundary marker(s): ";
    for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
      cout << Marker_Custom[iMarker_Custom];
      if (iMarker_Custom < nMarker_Custom-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_ActDisk_Inlet != 0) {
		cout << "Actuator disk (inlet) boundary marker(s): ";
		for (iMarker_ActDisk_Inlet = 0; iMarker_ActDisk_Inlet < nMarker_ActDisk_Inlet; iMarker_ActDisk_Inlet++) {
			cout << Marker_ActDisk_Inlet[iMarker_ActDisk_Inlet];
			if (iMarker_ActDisk_Inlet < nMarker_ActDisk_Inlet-1) cout << ", ";
			else cout <<"."<< endl;
		}
	}

  if (nMarker_ActDisk_Outlet != 0) {
		cout << "Actuator disk (outlet) boundary marker(s): ";
		for (iMarker_ActDisk_Outlet = 0; iMarker_ActDisk_Outlet < nMarker_ActDisk_Outlet; iMarker_ActDisk_Outlet++) {
			cout << Marker_ActDisk_Outlet[iMarker_ActDisk_Outlet];
			if (iMarker_ActDisk_Outlet < nMarker_ActDisk_Outlet-1) cout << ", ";
			else cout <<"."<< endl;
		}
	}

}


void CConfig::GetChemistryEquilConstants(double **RxnConstantTable, unsigned short iReaction) {
  switch (Kind_GasModel) {
    case O2:
      //O2 + M -> 2O + M
      RxnConstantTable[0][0] = 1.8103;	RxnConstantTable[0][1] = 1.9607;	RxnConstantTable[0][2] = 3.5716;	RxnConstantTable[0][3] = -7.3623;		RxnConstantTable[0][4] = 0.083861;
      RxnConstantTable[1][0] = 0.91354;	RxnConstantTable[1][1] = 2.3160;	RxnConstantTable[1][2] = 2.2885;	RxnConstantTable[1][3] = -6.7969;		RxnConstantTable[1][4] = 0.046338;
      RxnConstantTable[2][0] = 0.64183;	RxnConstantTable[2][1] = 2.4253;	RxnConstantTable[2][2] = 1.9026;	RxnConstantTable[2][3] = -6.6277;		RxnConstantTable[2][4] = 0.035151;
      RxnConstantTable[3][0] = 0.55388;	RxnConstantTable[3][1] = 2.4600;	RxnConstantTable[3][2] = 1.7763;	RxnConstantTable[3][3] = -6.5720;		RxnConstantTable[3][4] = 0.031445;
      RxnConstantTable[4][0] = 0.52455;	RxnConstantTable[4][1] = 2.4715;	RxnConstantTable[4][2] = 1.7342;	RxnConstantTable[4][3] = -6.55534;	RxnConstantTable[4][4] = 0.030209;
      RxnConstantTable[5][0] = 0.50989;	RxnConstantTable[5][1] = 2.4773;	RxnConstantTable[5][2] = 1.7132;	RxnConstantTable[5][3] = -6.5441;		RxnConstantTable[5][4] = 0.029591;
      break;
    case N2:
      //N2 + M -> 2N + M
      RxnConstantTable[0][0] = 3.4907;	RxnConstantTable[0][1] = 0.83133;	RxnConstantTable[0][2] = 4.0978;	RxnConstantTable[0][3] = -12.728;	RxnConstantTable[0][4] = 0.07487;		//n = 1E14
      RxnConstantTable[1][0] = 2.0723;	RxnConstantTable[1][1] = 1.38970;	RxnConstantTable[1][2] = 2.0617;	RxnConstantTable[1][3] = -11.828;	RxnConstantTable[1][4] = 0.015105;	//n = 1E15
      RxnConstantTable[2][0] = 1.6060;	RxnConstantTable[2][1] = 1.57320;	RxnConstantTable[2][2] = 1.3923;	RxnConstantTable[2][3] = -11.533;	RxnConstantTable[2][4] = -0.004543;	//n = 1E16
      RxnConstantTable[3][0] = 1.5351;	RxnConstantTable[3][1] = 1.60610;	RxnConstantTable[3][2] = 1.2993;	RxnConstantTable[3][3] = -11.494;	RxnConstantTable[3][4] = -0.00698;	//n = 1E17
      RxnConstantTable[4][0] = 1.4766;	RxnConstantTable[4][1] = 1.62910;	RxnConstantTable[4][2] = 1.2153;	RxnConstantTable[4][3] = -11.457;	RxnConstantTable[4][4] = -0.00944;	//n = 1E18
      RxnConstantTable[5][0] = 1.4766;	RxnConstantTable[5][1] = 1.62910;	RxnConstantTable[5][2] = 1.2153;	RxnConstantTable[5][3] = -11.457;	RxnConstantTable[5][4] = -0.00944;	//n = 1E19
      break;

    case ARGON_SID:
      //N2 + M -> 2N + M
      RxnConstantTable[0][0] = 3.4907;	RxnConstantTable[0][1] = 0.83133;	RxnConstantTable[0][2] = 4.0978;	RxnConstantTable[0][3] = -12.728;	RxnConstantTable[0][4] = 0.07487;		//n = 1E14
      RxnConstantTable[1][0] = 2.0723;	RxnConstantTable[1][1] = 1.38970;	RxnConstantTable[1][2] = 2.0617;	RxnConstantTable[1][3] = -11.828;	RxnConstantTable[1][4] = 0.015105;	//n = 1E15
      RxnConstantTable[2][0] = 1.6060;	RxnConstantTable[2][1] = 1.57320;	RxnConstantTable[2][2] = 1.3923;	RxnConstantTable[2][3] = -11.533;	RxnConstantTable[2][4] = -0.004543;	//n = 1E16
      RxnConstantTable[3][0] = 1.5351;	RxnConstantTable[3][1] = 1.60610;	RxnConstantTable[3][2] = 1.2993;	RxnConstantTable[3][3] = -11.494;	RxnConstantTable[3][4] = -0.00698;	//n = 1E17
      RxnConstantTable[4][0] = 1.4766;	RxnConstantTable[4][1] = 1.62910;	RxnConstantTable[4][2] = 1.2153;	RxnConstantTable[4][3] = -11.457;	RxnConstantTable[4][4] = -0.00944;	//n = 1E18
      RxnConstantTable[5][0] = 1.4766;	RxnConstantTable[5][1] = 1.62910;	RxnConstantTable[5][2] = 1.2153;	RxnConstantTable[5][3] = -11.457;	RxnConstantTable[5][4] = -0.00944;	//n = 1E19
      break;

    case AIR5:
      if (iReaction <= 4) {
        //N2 + M -> 2N + M
        RxnConstantTable[0][0] = 3.4907;	RxnConstantTable[0][1] = 0.83133;	RxnConstantTable[0][2] = 4.0978;	RxnConstantTable[0][3] = -12.728;	RxnConstantTable[0][4] = 0.07487;		//n = 1E14
        RxnConstantTable[1][0] = 2.0723;	RxnConstantTable[1][1] = 1.38970;	RxnConstantTable[1][2] = 2.0617;	RxnConstantTable[1][3] = -11.828;	RxnConstantTable[1][4] = 0.015105;	//n = 1E15
        RxnConstantTable[2][0] = 1.6060;	RxnConstantTable[2][1] = 1.57320;	RxnConstantTable[2][2] = 1.3923;	RxnConstantTable[2][3] = -11.533;	RxnConstantTable[2][4] = -0.004543;	//n = 1E16
        RxnConstantTable[3][0] = 1.5351;	RxnConstantTable[3][1] = 1.60610;	RxnConstantTable[3][2] = 1.2993;	RxnConstantTable[3][3] = -11.494;	RxnConstantTable[3][4] = -0.00698;	//n = 1E17
        RxnConstantTable[4][0] = 1.4766;	RxnConstantTable[4][1] = 1.62910;	RxnConstantTable[4][2] = 1.2153;	RxnConstantTable[4][3] = -11.457;	RxnConstantTable[4][4] = -0.00944;	//n = 1E18
        RxnConstantTable[5][0] = 1.4766;	RxnConstantTable[5][1] = 1.62910;	RxnConstantTable[5][2] = 1.2153;	RxnConstantTable[5][3] = -11.457;	RxnConstantTable[5][4] = -0.00944;	//n = 1E19
      } else if (iReaction > 4 && iReaction <= 9) {
        //O2 + M -> 2O + M
        RxnConstantTable[0][0] = 1.8103;	RxnConstantTable[0][1] = 1.9607;	RxnConstantTable[0][2] = 3.5716;	RxnConstantTable[0][3] = -7.3623;		RxnConstantTable[0][4] = 0.083861;
        RxnConstantTable[1][0] = 0.91354;	RxnConstantTable[1][1] = 2.3160;	RxnConstantTable[1][2] = 2.2885;	RxnConstantTable[1][3] = -6.7969;		RxnConstantTable[1][4] = 0.046338;
        RxnConstantTable[2][0] = 0.64183;	RxnConstantTable[2][1] = 2.4253;	RxnConstantTable[2][2] = 1.9026;	RxnConstantTable[2][3] = -6.6277;		RxnConstantTable[2][4] = 0.035151;
        RxnConstantTable[3][0] = 0.55388;	RxnConstantTable[3][1] = 2.4600;	RxnConstantTable[3][2] = 1.7763;	RxnConstantTable[3][3] = -6.5720;		RxnConstantTable[3][4] = 0.031445;
        RxnConstantTable[4][0] = 0.52455;	RxnConstantTable[4][1] = 2.4715;	RxnConstantTable[4][2] = 1.7342;	RxnConstantTable[4][3] = -6.55534;	RxnConstantTable[4][4] = 0.030209;
        RxnConstantTable[5][0] = 0.50989;	RxnConstantTable[5][1] = 2.4773;	RxnConstantTable[5][2] = 1.7132;	RxnConstantTable[5][3] = -6.5441;		RxnConstantTable[5][4] = 0.029591;
      } else if (iReaction > 9 && iReaction <= 14) {
        //NO + M -> N + O + M
        RxnConstantTable[0][0] = 2.1649;	RxnConstantTable[0][1] = 0.078577;	RxnConstantTable[0][2] = 2.8508;	RxnConstantTable[0][3] = -8.5422;	RxnConstantTable[0][4] = 0.053043;
        RxnConstantTable[1][0] = 1.0072;	RxnConstantTable[1][1] = 0.53545;		RxnConstantTable[1][2] = 1.1911;	RxnConstantTable[1][3] = -7.8098;	RxnConstantTable[1][4] = 0.004394;
        RxnConstantTable[2][0] = 0.63817;	RxnConstantTable[2][1] = 0.68189;		RxnConstantTable[2][2] = 0.66336;	RxnConstantTable[2][3] = -7.5773;	RxnConstantTable[2][4] = -0.011025;
        RxnConstantTable[3][0] = 0.55889;	RxnConstantTable[3][1] = 0.71558;		RxnConstantTable[3][2] = 0.55396;	RxnConstantTable[3][3] = -7.5304;	RxnConstantTable[3][4] = -0.014089;
        RxnConstantTable[4][0] = 0.5150;	RxnConstantTable[4][1] = 0.73286;		RxnConstantTable[4][2] = 0.49096;	RxnConstantTable[4][3] = -7.5025;	RxnConstantTable[4][4] = -0.015938;
        RxnConstantTable[5][0] = 0.50765;	RxnConstantTable[5][1] = 0.73575;		RxnConstantTable[5][2] = 0.48042;	RxnConstantTable[5][3] = -7.4979;	RxnConstantTable[5][4] = -0.016247;
      } else if (iReaction == 15) {
        //N2 + O -> NO + N
        RxnConstantTable[0][0] = 1.3261;	RxnConstantTable[0][1] = 0.75268;	RxnConstantTable[0][2] = 1.2474;	RxnConstantTable[0][3] = -4.1857;	RxnConstantTable[0][4] = 0.02184;
        RxnConstantTable[1][0] = 1.0653;	RxnConstantTable[1][1] = 0.85417;	RxnConstantTable[1][2] = 0.87093;	RxnConstantTable[1][3] = -4.0188;	RxnConstantTable[1][4] = 0.010721;
        RxnConstantTable[2][0] = 0.96794;	RxnConstantTable[2][1] = 0.89131;	RxnConstantTable[2][2] = 0.7291;	RxnConstantTable[2][3] = -3.9555;	RxnConstantTable[2][4] = 0.006488;
        RxnConstantTable[3][0] = 0.97646;	RxnConstantTable[3][1] = 0.89043;	RxnConstantTable[3][2] = 0.74572;	RxnConstantTable[3][3] = -3.9642;	RxnConstantTable[3][4] = 0.007123;
        RxnConstantTable[4][0] = 0.96188;	RxnConstantTable[4][1] = 0.89617;	RxnConstantTable[4][2] = 0.72479;	RxnConstantTable[4][3] = -3.955;	RxnConstantTable[4][4] = 0.006509;
        RxnConstantTable[5][0] = 0.96921;	RxnConstantTable[5][1] = 0.89329;	RxnConstantTable[5][2] = 0.73531;	RxnConstantTable[5][3] = -3.9596;	RxnConstantTable[5][4] = 0.006818;
      } else if (iReaction == 16) {
        //NO + O -> O2 + N
        RxnConstantTable[0][0] = 0.35438;		RxnConstantTable[0][1] = -1.8821;	RxnConstantTable[0][2] = -0.72111;	RxnConstantTable[0][3] = -1.1797;		RxnConstantTable[0][4] = -0.030831;
        RxnConstantTable[1][0] = 0.093613;	RxnConstantTable[1][1] = -1.7806;	RxnConstantTable[1][2] = -1.0975;		RxnConstantTable[1][3] = -1.0128;		RxnConstantTable[1][4] = -0.041949;
        RxnConstantTable[2][0] = -0.003732;	RxnConstantTable[2][1] = -1.7434;	RxnConstantTable[2][2] = -1.2394;		RxnConstantTable[2][3] = -0.94952;	RxnConstantTable[2][4] = -0.046182;
        RxnConstantTable[3][0] = 0.004815;	RxnConstantTable[3][1] = -1.7443;	RxnConstantTable[3][2] = -1.2227;		RxnConstantTable[3][3] = -0.95824;	RxnConstantTable[3][4] = -0.045545;
        RxnConstantTable[4][0] = -0.009758;	RxnConstantTable[4][1] = -1.7386;	RxnConstantTable[4][2] = -1.2436;		RxnConstantTable[4][3] = -0.949;		RxnConstantTable[4][4] = -0.046159;
        RxnConstantTable[5][0] = -0.002428;	RxnConstantTable[5][1] = -1.7415;	RxnConstantTable[5][2] = -1.2331;		RxnConstantTable[5][3] = -0.95365;	RxnConstantTable[5][4] = -0.04585;
      }
      break;
    case AIR7:
      if (iReaction <= 6) {
        //N2 + M -> 2N + M
        RxnConstantTable[0][0] = 3.4907;	RxnConstantTable[0][1] = 0.83133;	RxnConstantTable[0][2] = 4.0978;	RxnConstantTable[0][3] = -12.728;	RxnConstantTable[0][4] = 0.07487;		//n = 1E14
        RxnConstantTable[1][0] = 2.0723;	RxnConstantTable[1][1] = 1.38970;	RxnConstantTable[1][2] = 2.0617;	RxnConstantTable[1][3] = -11.828;	RxnConstantTable[1][4] = 0.015105;	//n = 1E15
        RxnConstantTable[2][0] = 1.6060;	RxnConstantTable[2][1] = 1.57320;	RxnConstantTable[2][2] = 1.3923;	RxnConstantTable[2][3] = -11.533;	RxnConstantTable[2][4] = -0.004543;	//n = 1E16
        RxnConstantTable[3][0] = 1.5351;	RxnConstantTable[3][1] = 1.60610;	RxnConstantTable[3][2] = 1.2993;	RxnConstantTable[3][3] = -11.494;	RxnConstantTable[3][4] = -0.00698;	//n = 1E17
        RxnConstantTable[4][0] = 1.4766;	RxnConstantTable[4][1] = 1.62910;	RxnConstantTable[4][2] = 1.2153;	RxnConstantTable[4][3] = -11.457;	RxnConstantTable[4][4] = -0.00944;	//n = 1E18
        RxnConstantTable[5][0] = 1.4766;	RxnConstantTable[5][1] = 1.62910;	RxnConstantTable[5][2] = 1.2153;	RxnConstantTable[5][3] = -11.457;	RxnConstantTable[5][4] = -0.00944;	//n = 1E19
      } else if (iReaction > 6 && iReaction <= 13) {
        //O2 + M -> 2O + M
        RxnConstantTable[0][0] = 1.8103;	RxnConstantTable[0][1] = 1.9607;	RxnConstantTable[0][2] = 3.5716;	RxnConstantTable[0][3] = -7.3623;		RxnConstantTable[0][4] = 0.083861;
        RxnConstantTable[1][0] = 0.91354;	RxnConstantTable[1][1] = 2.3160;	RxnConstantTable[1][2] = 2.2885;	RxnConstantTable[1][3] = -6.7969;		RxnConstantTable[1][4] = 0.046338;
        RxnConstantTable[2][0] = 0.64183;	RxnConstantTable[2][1] = 2.4253;	RxnConstantTable[2][2] = 1.9026;	RxnConstantTable[2][3] = -6.6277;		RxnConstantTable[2][4] = 0.035151;
        RxnConstantTable[3][0] = 0.55388;	RxnConstantTable[3][1] = 2.4600;	RxnConstantTable[3][2] = 1.7763;	RxnConstantTable[3][3] = -6.5720;		RxnConstantTable[3][4] = 0.031445;
        RxnConstantTable[4][0] = 0.52455;	RxnConstantTable[4][1] = 2.4715;	RxnConstantTable[4][2] = 1.7342;	RxnConstantTable[4][3] = -6.55534;	RxnConstantTable[4][4] = 0.030209;
        RxnConstantTable[5][0] = 0.50989;	RxnConstantTable[5][1] = 2.4773;	RxnConstantTable[5][2] = 1.7132;	RxnConstantTable[5][3] = -6.5441;		RxnConstantTable[5][4] = 0.029591;
      } else if (iReaction > 13 && iReaction <= 20) {
        //NO + M -> N + O + M
        RxnConstantTable[0][0] = 2.1649;	RxnConstantTable[0][1] = 0.078577;	RxnConstantTable[0][2] = 2.8508;	RxnConstantTable[0][3] = -8.5422;	RxnConstantTable[0][4] = 0.053043;
        RxnConstantTable[1][0] = 1.0072;	RxnConstantTable[1][1] = 0.53545;		RxnConstantTable[1][2] = 1.1911;	RxnConstantTable[1][3] = -7.8098;	RxnConstantTable[1][4] = 0.004394;
        RxnConstantTable[2][0] = 0.63817;	RxnConstantTable[2][1] = 0.68189;		RxnConstantTable[2][2] = 0.66336;	RxnConstantTable[2][3] = -7.5773;	RxnConstantTable[2][4] = -0.011025;
        RxnConstantTable[3][0] = 0.55889;	RxnConstantTable[3][1] = 0.71558;		RxnConstantTable[3][2] = 0.55396;	RxnConstantTable[3][3] = -7.5304;	RxnConstantTable[3][4] = -0.014089;
        RxnConstantTable[4][0] = 0.5150;	RxnConstantTable[4][1] = 0.73286;		RxnConstantTable[4][2] = 0.49096;	RxnConstantTable[4][3] = -7.5025;	RxnConstantTable[4][4] = -0.015938;
        RxnConstantTable[5][0] = 0.50765;	RxnConstantTable[5][1] = 0.73575;		RxnConstantTable[5][2] = 0.48042;	RxnConstantTable[5][3] = -7.4979;	RxnConstantTable[5][4] = -0.016247;
      } else if (iReaction == 21) {
        //N2 + O -> NO + N
        RxnConstantTable[0][0] = 1.3261;	RxnConstantTable[0][1] = 0.75268;	RxnConstantTable[0][2] = 1.2474;	RxnConstantTable[0][3] = -4.1857;	RxnConstantTable[0][4] = 0.02184;
        RxnConstantTable[1][0] = 1.0653;	RxnConstantTable[1][1] = 0.85417;	RxnConstantTable[1][2] = 0.87093;	RxnConstantTable[1][3] = -4.0188;	RxnConstantTable[1][4] = 0.010721;
        RxnConstantTable[2][0] = 0.96794;	RxnConstantTable[2][1] = 0.89131;	RxnConstantTable[2][2] = 0.7291;	RxnConstantTable[2][3] = -3.9555;	RxnConstantTable[2][4] = 0.006488;
        RxnConstantTable[3][0] = 0.97646;	RxnConstantTable[3][1] = 0.89043;	RxnConstantTable[3][2] = 0.74572;	RxnConstantTable[3][3] = -3.9642;	RxnConstantTable[3][4] = 0.007123;
        RxnConstantTable[4][0] = 0.96188;	RxnConstantTable[4][1] = 0.89617;	RxnConstantTable[4][2] = 0.72479;	RxnConstantTable[4][3] = -3.955;	RxnConstantTable[4][4] = 0.006509;
        RxnConstantTable[5][0] = 0.96921;	RxnConstantTable[5][1] = 0.89329;	RxnConstantTable[5][2] = 0.73531;	RxnConstantTable[5][3] = -3.9596;	RxnConstantTable[5][4] = 0.006818;
      } else if (iReaction == 22) {
        //NO + O -> O2 + N
        RxnConstantTable[0][0] = 0.35438;		RxnConstantTable[0][1] = -1.8821;	RxnConstantTable[0][2] = -0.72111;	RxnConstantTable[0][3] = -1.1797;		RxnConstantTable[0][4] = -0.030831;
        RxnConstantTable[1][0] = 0.093613;	RxnConstantTable[1][1] = -1.7806;	RxnConstantTable[1][2] = -1.0975;		RxnConstantTable[1][3] = -1.0128;		RxnConstantTable[1][4] = -0.041949;
        RxnConstantTable[2][0] = -0.003732;	RxnConstantTable[2][1] = -1.7434;	RxnConstantTable[2][2] = -1.2394;		RxnConstantTable[2][3] = -0.94952;	RxnConstantTable[2][4] = -0.046182;
        RxnConstantTable[3][0] = 0.004815;	RxnConstantTable[3][1] = -1.7443;	RxnConstantTable[3][2] = -1.2227;		RxnConstantTable[3][3] = -0.95824;	RxnConstantTable[3][4] = -0.045545;
        RxnConstantTable[4][0] = -0.009758;	RxnConstantTable[4][1] = -1.7386;	RxnConstantTable[4][2] = -1.2436;		RxnConstantTable[4][3] = -0.949;		RxnConstantTable[4][4] = -0.046159;
        RxnConstantTable[5][0] = -0.002428;	RxnConstantTable[5][1] = -1.7415;	RxnConstantTable[5][2] = -1.2331;		RxnConstantTable[5][3] = -0.95365;	RxnConstantTable[5][4] = -0.04585;
      } else if (iReaction == 23) {
        //N + O -> NO+ + e-
        RxnConstantTable[0][0] = -2.1852;		RxnConstantTable[0][1] = -6.6709;	RxnConstantTable[0][2] = -4.2968;	RxnConstantTable[0][3] = -2.2175;	RxnConstantTable[0][4] = -0.050748;
        RxnConstantTable[1][0] = -1.0276;		RxnConstantTable[1][1] = -7.1278;	RxnConstantTable[1][2] = -2.637;	RxnConstantTable[1][3] = -2.95;		RxnConstantTable[1][4] = -0.0021;
        RxnConstantTable[2][0] = -0.65871;	RxnConstantTable[2][1] = -7.2742;	RxnConstantTable[2][2] = -2.1096;	RxnConstantTable[2][3] = -3.1823;	RxnConstantTable[2][4] = 0.01331;
        RxnConstantTable[3][0] = -0.57924;	RxnConstantTable[3][1] = -7.3079;	RxnConstantTable[3][2] = -1.9999;	RxnConstantTable[3][3] = -3.2294;	RxnConstantTable[3][4] = 0.016382;
        RxnConstantTable[4][0] = -0.53538;	RxnConstantTable[4][1] = -7.3252;	RxnConstantTable[4][2] = -1.937;	RxnConstantTable[4][3] = -3.2572;	RxnConstantTable[4][4] = 0.01823;
        RxnConstantTable[5][0] = -0.52801;	RxnConstantTable[5][1] = -7.3281;	RxnConstantTable[5][2] = -1.9264;	RxnConstantTable[5][3] = -3.2618;	RxnConstantTable[5][4] = 0.01854;
      }
      break;
  }
}

bool CConfig::TokenizeString(string & str, string & option_name,
                             vector<string> & option_value) {
  const string delimiters(" ()[]{}:,\t\n\v\f\r");
  // check for comments or empty string
  string::size_type pos, last_pos;
  pos = str.find_first_of("%");
  if ( (str.length() == 0) || (pos == 0) ) {
    // str is empty or a comment line, so no option here
    return false;
  }
  if (pos != string::npos) {
    // remove comment at end if necessary
    str.erase(pos);
  }

  // look for line composed on only delimiters (usually whitespace)
  pos = str.find_first_not_of(delimiters);
  if (pos == string::npos) {
    return false;
  }

  // find the equals sign and split string
  string name_part, value_part;
  pos = str.find("=");
  if (pos == string::npos) {
    cerr << "Error in TokenizeString(): "
    << "line in the configuration file with no \"=\" sign."
    << endl;
    cout << "Look for: " << str << endl;
    cout << "str.length() = " << str.length() << endl;
    throw(-1);
  }
  name_part = str.substr(0, pos);
  value_part = str.substr(pos+1,string::npos);
  //cout << "name_part  = |" << name_part  << "|" << endl;
  //cout << "value_part = |" << value_part << "|" << endl;

  // the first_part should consist of one string with no interior delimiters
  last_pos = name_part.find_first_not_of(delimiters, 0);
  pos = name_part.find_first_of(delimiters, last_pos);
  if ( (name_part.length() == 0) || (last_pos == string::npos) ) {
    cerr << "Error in CConfig::TokenizeString(): "
    << "line in the configuration file with no name before the \"=\" sign."
    << endl;
    throw(-1);
  }
  if (pos == string::npos) pos = name_part.length();
  option_name = name_part.substr(last_pos, pos - last_pos);
  last_pos = name_part.find_first_not_of(delimiters, pos);
  if (last_pos != string::npos) {
    cerr << "Error in TokenizeString(): "
    << "two or more options before an \"=\" sign in the configuration file."
    << endl;
    throw(-1);
  }
  StringToUpperCase(option_name);

  //cout << "option_name = |" << option_name << "|" << endl;
  //cout << "pos = " << pos << ": last_pos = " << last_pos << endl;

  // now fill the option value vector
  option_value.clear();
  last_pos = value_part.find_first_not_of(delimiters, 0);
  pos = value_part.find_first_of(delimiters, last_pos);
  while (string::npos != pos || string::npos != last_pos) {
    // add token to the vector<string>
    option_value.push_back(value_part.substr(last_pos, pos - last_pos));
    // skip delimiters
    last_pos = value_part.find_first_not_of(delimiters, pos);
    // find next "non-delimiter"
    pos = value_part.find_first_of(delimiters, last_pos);
  }
  if (option_value.size() == 0) {
    cerr << "Error inT okenizeString(): "
    << "option " << option_name << " in configuration file with no value assigned."
    << endl;
    throw(-1);
  }

#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif

  // look for ';' DV delimiters attached to values
  vector<string>::iterator it;
  it = option_value.begin();
  while (it != option_value.end()) {
    if (it->compare(";") == 0) {
      it++;
      continue;
    }

    pos = it->find(';');
    if (pos != string::npos) {
      string before_semi = it->substr(0, pos);
      string after_semi= it->substr(pos+1,string::npos);
      if (before_semi.empty()) {
        *it = ";";
        it++;
        option_value.insert(it, after_semi);
      } else {
        *it = before_semi;
        it++;
        vector<string> to_insert;
        to_insert.push_back(";");
        if (!after_semi.empty())
          to_insert.push_back(after_semi);
        option_value.insert(it, to_insert.begin(), to_insert.end());
      }
      it = option_value.begin(); // go back to beginning; not efficient
      continue;
    } else {
      it++;
    }
  }
#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  // remove any consecutive ";"
  it = option_value.begin();
  bool semi_at_prev = false;
  while (it != option_value.end()) {
    if (semi_at_prev) {
      if (it->compare(";") == 0) {
        option_value.erase(it);
        it = option_value.begin();
        semi_at_prev = false;
        continue;
      }
    }
    if (it->compare(";") == 0) {
      semi_at_prev = true;
    } else {
      semi_at_prev = false;
    }
    it++;
  }

#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  return true;
}

unsigned short CConfig::GetMarker_CfgFile_TagBound(string val_marker) {

  unsigned short iMarker_Config;

  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker)
      return iMarker_Config;

  cout <<"The configuration file doesn't have any definition for marker "<< val_marker <<"!!" << endl;
  exit(EXIT_FAILURE);
  
}

string CConfig::GetMarker_CfgFile_TagBound(unsigned short val_marker) {
  return Marker_CfgFile_TagBound[val_marker];
}

unsigned short CConfig::GetMarker_CfgFile_KindBC(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_KindBC[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_Monitoring(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_Monitoring[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_GeoEval(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_GeoEval[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_Designing(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_Designing[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_Plotting(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_Plotting[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_Out_1D(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_Out_1D[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_DV(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_DV[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_Moving(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_Moving[iMarker_Config];
}

unsigned short CConfig::GetMarker_CfgFile_PerBound(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_CfgFile_TagBound[iMarker_Config] == val_marker) break;
  return Marker_CfgFile_PerBound[iMarker_Config];
}

CConfig::~CConfig(void) {
  
  if (RK_Alpha_Step!=NULL) delete [] RK_Alpha_Step;
  if (MG_PreSmooth!=NULL) delete [] MG_PreSmooth;
  if (MG_PostSmooth!=NULL) delete [] MG_PostSmooth;
  if (U_FreeStreamND!=NULL) delete [] U_FreeStreamND;

  /*--- If allocated, delete arrays for Plasma solver ---*/
  
  if (Molar_Mass           != NULL) delete [] Molar_Mass;
  if (Gas_Composition      != NULL) delete [] Gas_Composition;
  if (Enthalpy_Formation   != NULL) delete [] Enthalpy_Formation;
  if (ArrheniusCoefficient != NULL) delete [] ArrheniusCoefficient;
  if (ArrheniusEta         != NULL) delete [] ArrheniusEta;
  if (ArrheniusTheta       != NULL) delete [] ArrheniusTheta;
  if (CharVibTemp          != NULL) delete [] CharVibTemp;
  if (CharElTemp           != NULL) {
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      delete[] CharElTemp[iSpecies];
    delete [] CharElTemp;
  }
  if (degen != NULL) {
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      delete[] degen[iSpecies];
    delete [] degen;
  }
  unsigned short ii, iReaction;
  if (Reactions != NULL) {
    for (iReaction = 0; iReaction < nReactions; iReaction++) {
      for (ii = 0; ii < 2; ii++) {
        delete [] Reactions[iReaction][ii];
      }
      delete[] Reactions[iReaction];
    }
    delete [] Reactions;
  }

  /*--- Free memory for Aeroelastic problems. ---*/
  
  if (Grid_Movement && Aeroelastic_Simulation) {

    delete[] Aeroelastic_pitch;
    delete[] Aeroelastic_plunge;
  }

  /*--- Free memory for unspecified grid motion parameters ---*/

  if (Kind_GridMovement != NULL)
    delete [] Kind_GridMovement;

  /*--- motion origin: ---*/
  
  if (Motion_Origin_X != NULL)
    delete [] Motion_Origin_X;
  if (Motion_Origin_Y != NULL)
    delete [] Motion_Origin_Y;
  if (Motion_Origin_Z != NULL)
    delete [] Motion_Origin_Z;
  if (MoveMotion_Origin != NULL)
    delete [] MoveMotion_Origin;

  /*--- rotation: ---*/
  
  if (Rotation_Rate_X != NULL)
    delete [] Rotation_Rate_X;
  if (Rotation_Rate_Y != NULL)
    delete [] Rotation_Rate_Y;
  if (Rotation_Rate_Z != NULL)
    delete [] Rotation_Rate_Z;

  /*--- pitching: ---*/
  
  if (Pitching_Omega_X != NULL)
    delete [] Pitching_Omega_X;
  if (Pitching_Omega_Y != NULL)
    delete [] Pitching_Omega_Y;
  if (Pitching_Omega_Z != NULL)
    delete [] Pitching_Omega_Z;

  /*--- pitching amplitude: ---*/
  
  if (Pitching_Ampl_X != NULL)
    delete [] Pitching_Ampl_X;
  if (Pitching_Ampl_Y != NULL)
    delete [] Pitching_Ampl_Y;
  if (Pitching_Ampl_Z != NULL)
    delete [] Pitching_Ampl_Z;

  /*--- pitching phase: ---*/
  
  if (Pitching_Phase_X != NULL)
    delete [] Pitching_Phase_X;
  if (Pitching_Phase_Y != NULL)
    delete [] Pitching_Phase_Y;
  if (Pitching_Phase_Z != NULL)
    delete [] Pitching_Phase_Z;

  /*--- plunging: ---*/
  
  if (Plunging_Omega_X != NULL)
    delete [] Plunging_Omega_X;
  if (Plunging_Omega_Y != NULL)
    delete [] Plunging_Omega_Y;
  if (Plunging_Omega_Z != NULL)
    delete [] Plunging_Omega_Z;

  /*--- plunging amplitude: ---*/
  
  if (Plunging_Ampl_X != NULL)
    delete [] Plunging_Ampl_X;
  if (Plunging_Ampl_Y != NULL)
    delete [] Plunging_Ampl_Y;
  if (Plunging_Ampl_Z != NULL)
    delete [] Plunging_Ampl_Z;

  if (RefOriginMoment != NULL)
    delete [] RefOriginMoment;
  if (RefOriginMoment_X != NULL)
    delete [] RefOriginMoment_X;
  if (RefOriginMoment_Y != NULL)
    delete [] RefOriginMoment_Y;
  if (RefOriginMoment_Z != NULL)
    delete [] RefOriginMoment_Z;

  /*--- Marker pointers ---*/
  
  if (Marker_CfgFile_Out_1D!=NULL)   delete[] Marker_CfgFile_Out_1D;
  if (Marker_All_Out_1D!=NULL)      delete[] Marker_All_Out_1D;
  if (Marker_CfgFile_GeoEval!=NULL)  delete[] Marker_CfgFile_GeoEval;
  if (Marker_All_GeoEval!=NULL)     delete[] Marker_All_GeoEval;
  if (Marker_CfgFile_TagBound!=NULL)      delete[] Marker_CfgFile_TagBound;
  if (Marker_All_TagBound!=NULL)         delete[] Marker_All_TagBound;
  if (Marker_CfgFile_KindBC!=NULL) delete[] Marker_CfgFile_KindBC;
  if (Marker_All_KindBC!=NULL)    delete[] Marker_All_KindBC;
  if (Marker_CfgFile_Monitoring!=NULL)    delete[] Marker_CfgFile_Monitoring;
  if (Marker_All_Monitoring!=NULL)   delete[] Marker_All_Monitoring;
  if (Marker_CfgFile_Designing!=NULL) delete[] Marker_CfgFile_Designing;
  if (Marker_All_Designing!=NULL)    delete[] Marker_All_Designing;
  if (Marker_CfgFile_Plotting!=NULL)  delete[] Marker_CfgFile_Plotting;
  if (Marker_All_Plotting!=NULL)     delete[] Marker_All_Plotting;
  if (Marker_CfgFile_DV!=NULL)        delete[] Marker_CfgFile_DV;
  if (Marker_All_DV!=NULL)           delete[] Marker_All_DV;
  if (Marker_CfgFile_Moving!=NULL)   delete[] Marker_CfgFile_Moving;
  if (Marker_All_Moving!=NULL)      delete[] Marker_All_Moving;
  if (Marker_CfgFile_PerBound!=NULL) delete[] Marker_CfgFile_PerBound;
  if (Marker_All_PerBound!=NULL)    delete[] Marker_All_PerBound;

  if (Marker_DV!=NULL)               delete[] Marker_DV;
  if (Marker_Moving!=NULL)           delete[] Marker_Moving;
  if (Marker_Monitoring!=NULL)      delete[] Marker_Monitoring;
  if (Marker_Designing!=NULL)       delete[] Marker_Designing;
  if (Marker_GeoEval!=NULL)         delete[] Marker_GeoEval;
  if (Marker_Plotting!=NULL)        delete[] Marker_Plotting;
  if (Marker_All_SendRecv!=NULL)    delete[] Marker_All_SendRecv;

  if (EA_IntLimit!=NULL)    delete[] EA_IntLimit;
  if (Hold_GridFixed_Coord!=NULL)    delete[] Hold_GridFixed_Coord ;
  if (Subsonic_Engine_Box!=NULL)    delete[] Subsonic_Engine_Box ;
  if (DV_Value!=NULL)    delete[] DV_Value;
  if (Design_Variable!=NULL)    delete[] Design_Variable;
  if (Dirichlet_Value!=NULL)    delete[] Dirichlet_Value;
  if (Nozzle_Ttotal!=NULL)    delete[]  Nozzle_Ttotal;
  if (Nozzle_Ptotal!=NULL)    delete[]  Nozzle_Ptotal;
  if (Inlet_Ttotal!=NULL)    delete[]  Inlet_Ttotal;
  if (Inlet_Ptotal!=NULL)    delete[]  Inlet_Ptotal;
  if (Inlet_FlowDir!=NULL)    delete[] Inlet_FlowDir;
  if (Inlet_Temperature!=NULL)    delete[] Inlet_Temperature;
  if (Inlet_Pressure!=NULL)    delete[] Inlet_Pressure;
  if (Inlet_Velocity!=NULL)    delete[] Inlet_Velocity ;
  if (Inflow_Mach_Target!=NULL)    delete[] Inflow_Mach_Target;
  if (Inflow_Mach!=NULL)    delete[]  Inflow_Mach;
  if (Inflow_Pressure!=NULL)    delete[] Inflow_Pressure;
  if (Bleed_MassFlow_Target!=NULL)    delete[] Bleed_MassFlow_Target;
  if (Bleed_MassFlow!=NULL)    delete[]  Bleed_MassFlow;
  if (Bleed_Temperature_Target!=NULL)    delete[] Bleed_Temperature_Target;
  if (Bleed_Temperature!=NULL)    delete[]  Bleed_Temperature;
  if (Bleed_Pressure!=NULL)    delete[] Bleed_Pressure;
  if (Outlet_Pressure!=NULL)    delete[] Outlet_Pressure;
  if (Isothermal_Temperature!=NULL)    delete[] Isothermal_Temperature;
  if (Heat_Flux!=NULL)    delete[] Heat_Flux;
  if (Heat_FluxNonCatalytic!=NULL)    delete[] Heat_FluxNonCatalytic;
  if (Heat_FluxCatalytic!=NULL)    delete[] Heat_FluxCatalytic;
  if (Displ_Value!=NULL)    delete[] Displ_Value;
  if (Load_Value!=NULL)    delete[] Load_Value;
  if (FlowLoad_Value!=NULL)    delete[] FlowLoad_Value;
  if (Periodic_RotCenter!=NULL)    delete[] Periodic_RotCenter;
  if (Periodic_RotAngles!=NULL)    delete[] Periodic_RotAngles;
  if (Periodic_Translation!=NULL)    delete[] Periodic_Translation;
  if (Periodic_Center!=NULL)    delete[] Periodic_Center;
  if (Periodic_Rotation!=NULL)    delete[] Periodic_Rotation;
  if (Periodic_Translate!=NULL)    delete[] Periodic_Translate;

  if (ParamDV!=NULL  )    delete[] ParamDV;
  if (MG_CorrecSmooth!=NULL    )    delete[] MG_CorrecSmooth;
  if (Section_Location!=NULL)    delete[] Section_Location;
  if (Kappa_Flow!=NULL      )    delete[] Kappa_Flow;
  if (Kappa_AdjFlow!=NULL             )    delete[] Kappa_AdjFlow;
  if (Kappa_TNE2!=NULL   )    delete[] Kappa_TNE2;
  if (Kappa_AdjTNE2!=NULL        )    delete[] Kappa_AdjTNE2;
  if (Kappa_LinFlow!=NULL  )    delete[] Kappa_LinFlow;
  if (PlaneTag!=NULL)    delete[] PlaneTag;
  if (CFLRamp!=NULL)    delete[] CFLRamp;
  if (CFL!=NULL)    delete[] CFL;
  
  /*--- String markers ---*/
  if (Marker_Euler!=NULL )              delete[] Marker_Euler;
  if (Marker_FarField!=NULL )           delete[] Marker_FarField;
  if (Marker_Custom!=NULL )             delete[] Marker_Custom;
  if (Marker_SymWall!=NULL )            delete[] Marker_SymWall;
  if (Marker_Pressure!=NULL )           delete[] Marker_Pressure;
  if (Marker_PerBound!=NULL )           delete[] Marker_PerBound;
  if (Marker_PerDonor!=NULL )           delete[] Marker_PerDonor;
  if (Marker_NearFieldBound!=NULL )     delete[] Marker_NearFieldBound;
  if (Marker_InterfaceBound!=NULL )     delete[] Marker_InterfaceBound;
  if (Marker_Dirichlet!=NULL )          delete[] Marker_Dirichlet;
  if (Marker_Dirichlet_Elec!=NULL )     delete[] Marker_Dirichlet_Elec;
  if (Marker_Inlet!=NULL )              delete[] Marker_Inlet;
  if (Marker_Supersonic_Inlet!=NULL )   delete[] Marker_Supersonic_Inlet;
  if (Marker_Outlet!=NULL )             delete[] Marker_Outlet;
  if (Marker_Out_1D!=NULL )             delete[] Marker_Out_1D;
  if (Marker_Isothermal!=NULL )         delete[] Marker_Isothermal;
  if (Marker_EngineInflow!=NULL )      delete[] Marker_EngineInflow;
  if (Marker_EngineBleed!=NULL )      delete[] Marker_EngineBleed;
  if (Marker_EngineExhaust!=NULL )     delete[] Marker_EngineExhaust;
  if (Marker_Displacement!=NULL )       delete[] Marker_Displacement;
  if (Marker_Load!=NULL )               delete[] Marker_Load;
  if (Marker_FlowLoad!=NULL )           delete[] Marker_FlowLoad;
  if (Marker_Neumann!=NULL )            delete[] Marker_Neumann;
  if (Marker_Neumann_Elec!=NULL )       delete[] Marker_Neumann_Elec;
  if (Marker_IsothermalNonCatalytic!=NULL ) delete[] Marker_IsothermalNonCatalytic;
  if (Marker_IsothermalCatalytic!=NULL )    delete[] Marker_IsothermalCatalytic;
  if (Marker_HeatFlux!=NULL )               delete[] Marker_HeatFlux;
  if (Marker_HeatFluxNonCatalytic!=NULL )   delete[] Marker_HeatFluxNonCatalytic;
  if (Marker_HeatFluxCatalytic!=NULL )      delete[] Marker_HeatFluxCatalytic;

}

string CConfig::GetUnsteady_FileName(string val_filename, int val_iter) {

  string UnstExt, UnstFilename = val_filename;
  char buffer[50];

  /*--- Check that a positive value iteration is requested (for now). ---*/
  if (val_iter < 0) {
    cout << "Requesting a negative iteration number for the restart file!!" << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Append iteration number for unsteady cases ---*/
  if ((Wrt_Unsteady) || (Unsteady_Simulation == TIME_SPECTRAL)) {
    unsigned short lastindex = UnstFilename.find_last_of(".");
    UnstFilename = UnstFilename.substr(0, lastindex);
    if ((val_iter >= 0)    && (val_iter < 10))    sprintf (buffer, "_0000%d.dat", val_iter);
    if ((val_iter >= 10)   && (val_iter < 100))   sprintf (buffer, "_000%d.dat",  val_iter);
    if ((val_iter >= 100)  && (val_iter < 1000))  sprintf (buffer, "_00%d.dat",   val_iter);
    if ((val_iter >= 1000) && (val_iter < 10000)) sprintf (buffer, "_0%d.dat",    val_iter);
    if (val_iter >= 10000) sprintf (buffer, "_%d.dat", val_iter);
    string UnstExt = string(buffer);
    UnstFilename.append(UnstExt);
  }

  return UnstFilename;
}

string CConfig::GetObjFunc_Extension(string val_filename) {

  string AdjExt, Filename = val_filename;

  if (Adjoint) {

    /*--- Remove filename extension (.dat) ---*/
    unsigned short lastindex = Filename.find_last_of(".");
    Filename = Filename.substr(0, lastindex);

    switch (Kind_ObjFunc) {
      case DRAG_COEFFICIENT:        AdjExt = "_cd";       break;
      case LIFT_COEFFICIENT:        AdjExt = "_cl";       break;
      case SIDEFORCE_COEFFICIENT:   AdjExt = "_csf";      break;
      case INVERSE_DESIGN_PRESSURE: AdjExt = "_invpress"; break;
      case INVERSE_DESIGN_HEATFLUX: AdjExt = "_invheat";  break;
      case MOMENT_X_COEFFICIENT:    AdjExt = "_cmx";      break;
      case MOMENT_Y_COEFFICIENT:    AdjExt = "_cmy";      break;
      case MOMENT_Z_COEFFICIENT:    AdjExt = "_cmz";      break;
      case EFFICIENCY:              AdjExt = "_eff";      break;
      case EQUIVALENT_AREA:         AdjExt = "_ea";       break;
      case NEARFIELD_PRESSURE:      AdjExt = "_nfp";      break;
      case FORCE_X_COEFFICIENT:     AdjExt = "_cfx";      break;
      case FORCE_Y_COEFFICIENT:     AdjExt = "_cfy";      break;
      case FORCE_Z_COEFFICIENT:     AdjExt = "_cfz";      break;
      case THRUST_COEFFICIENT:      AdjExt = "_ct";       break;
      case TORQUE_COEFFICIENT:      AdjExt = "_cq";       break;
      case TOTAL_HEATFLUX:          AdjExt = "_totheat";  break;
      case MAXIMUM_HEATFLUX:        AdjExt = "_maxheat";  break;
      case FIGURE_OF_MERIT:         AdjExt = "_merit";    break;
      case FREE_SURFACE:            AdjExt = "_fs";       break;
      case AVG_TOTAL_PRESSURE:      AdjExt = "_pt";       break;
      case MASS_FLOW_RATE:          AdjExt = "_mfr";       break;
    }
    Filename.append(AdjExt);

    /*--- Lastly, add the .dat extension ---*/
    Filename.append(".dat");

  }

  return Filename;
}

unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {

  switch (val_eqsystem) {
    case RUNTIME_FLOW_SYS:      return FLOW_SOL;
    case RUNTIME_TURB_SYS:      return TURB_SOL;
    case RUNTIME_TNE2_SYS:      return TNE2_SOL;
    case RUNTIME_TRANS_SYS:     return TRANS_SOL;
    case RUNTIME_POISSON_SYS:   return POISSON_SOL;
    case RUNTIME_WAVE_SYS:      return WAVE_SOL;
    case RUNTIME_HEAT_SYS:      return HEAT_SOL;
    case RUNTIME_FEA_SYS:       return FEA_SOL;
    case RUNTIME_ADJPOT_SYS:    return ADJFLOW_SOL;
    case RUNTIME_ADJFLOW_SYS:   return ADJFLOW_SOL;
    case RUNTIME_ADJTURB_SYS:   return ADJTURB_SOL;
    case RUNTIME_ADJTNE2_SYS:   return ADJTNE2_SOL;
    case RUNTIME_LINPOT_SYS:    return LINFLOW_SOL;
    case RUNTIME_LINFLOW_SYS:   return LINFLOW_SOL;
    case RUNTIME_MULTIGRID_SYS: return 0;
  }
  return 0;
}

void CConfig::SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme,
                                    unsigned short val_kind_centered, unsigned short val_kind_upwind,
                                    unsigned short val_kind_slopelimit, unsigned short val_order_spatial_int) {

  Kind_ConvNumScheme = val_kind_convnumscheme;
  Kind_Centered = val_kind_centered;
  Kind_Upwind = val_kind_upwind;
  Kind_SlopeLimit = val_kind_slopelimit;
  SpatialOrder = val_order_spatial_int;

}

void CConfig::UpdateCFL(unsigned long val_iter) {
  double coeff;
  unsigned short iCFL;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (Adjoint) coeff = CFLRedCoeff_AdjFlow;
  else coeff = 1.0;

  if ((CFLRamp[0] != 1.0) && (val_iter % int(CFLRamp[1]) == 0 ) && (val_iter != 0) && (CFL[0] < CFLRamp[2]*coeff)) {

    for (iCFL = 0; iCFL <= nMultiLevel; iCFL++)
        CFL[iCFL] *= CFLRamp[0];

    if (rank == MASTER_NODE) {
      cout <<"\n New value of the CFL number: ";
      for (iCFL = 0; iCFL < nMultiLevel; iCFL++)
        cout << CFL[iCFL] <<", ";
      cout << CFL[nMultiLevel] <<".\n"<< endl;
    }

  }
  
}

void CConfig::SetGlobalParam(unsigned short val_solver,
                             unsigned short val_system,
                             unsigned long val_extiter) {

  /*--- Set the simulation global time ---*/
  Current_UnstTime = static_cast<double>(val_extiter)*Delta_UnstTime;
  Current_UnstTimeND = static_cast<double>(val_extiter)*Delta_UnstTimeND;

  /*--- Set the solver methods ---*/
  switch (val_solver) {
    case EULER:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      break;
    case NAVIER_STOKES:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      break;
    case RANS:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_TURB_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Turb, Kind_Centered_Turb,
                              Kind_Upwind_Turb, Kind_SlopeLimit_Turb,
                              SpatialOrder_Turb);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Turb);
      }
      if (val_system == RUNTIME_TRANS_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Turb, Kind_Centered_Turb,
                              Kind_Upwind_Turb, Kind_SlopeLimit_Turb,
                              SpatialOrder_Turb);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Turb);
      }
      break;
    case TNE2_EULER:
      if (val_system == RUNTIME_TNE2_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_TNE2, Kind_Centered_TNE2,
                              Kind_Upwind_TNE2, Kind_SlopeLimit_TNE2,
                              SpatialOrder_TNE2);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_TNE2);
      }
      break;
    case TNE2_NAVIER_STOKES:
      if (val_system == RUNTIME_TNE2_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_TNE2, Kind_Centered_TNE2,
                              Kind_Upwind_TNE2, Kind_SlopeLimit_TNE2,
                              SpatialOrder_TNE2);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_TNE2);
      }
      break;
    case ADJ_EULER:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_ADJFLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow,
                              Kind_Upwind_AdjFlow, Kind_SlopeLimit_AdjFlow,
                              SpatialOrder_AdjFlow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjFlow);
      }
      break;
    case ADJ_NAVIER_STOKES:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_ADJFLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow,
                              Kind_Upwind_AdjFlow, Kind_SlopeLimit_AdjFlow,
                              SpatialOrder_AdjFlow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjFlow);
      }
      break;
    case ADJ_RANS:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_ADJFLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow,
                              Kind_Upwind_AdjFlow, Kind_SlopeLimit_AdjFlow,
                              SpatialOrder_AdjFlow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjFlow);
      }
      if (val_system == RUNTIME_TURB_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Turb, Kind_Centered_Turb,
                              Kind_Upwind_Turb, Kind_SlopeLimit_Turb,
                              SpatialOrder_Turb);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Turb);
      }
      if (val_system == RUNTIME_ADJTURB_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjTurb, Kind_Centered_AdjTurb,
                              Kind_Upwind_AdjTurb, Kind_SlopeLimit_AdjTurb,
                              SpatialOrder_AdjTurb);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjTurb);
      }
      break;
    case ADJ_TNE2_EULER:
      if (val_system == RUNTIME_TNE2_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_TNE2, Kind_Centered_TNE2,
                              Kind_Upwind_TNE2, Kind_SlopeLimit_TNE2,
                              SpatialOrder_TNE2);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_TNE2);
      }
      if (val_system == RUNTIME_ADJTNE2_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjTNE2, Kind_Centered_AdjTNE2,
                              Kind_Upwind_AdjTNE2, Kind_SlopeLimit_AdjTNE2,
                              SpatialOrder_AdjTNE2);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjTNE2);
      }
      break;
    case ADJ_TNE2_NAVIER_STOKES:
      if (val_system == RUNTIME_TNE2_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_TNE2, Kind_Centered_TNE2,
                              Kind_Upwind_TNE2, Kind_SlopeLimit_TNE2,
                              SpatialOrder_TNE2);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_TNE2);
      }
      if (val_system == RUNTIME_ADJTNE2_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjTNE2, Kind_Centered_AdjTNE2,
                              Kind_Upwind_AdjTNE2, Kind_SlopeLimit_AdjTNE2,
                              SpatialOrder_AdjTNE2);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjTNE2);
      }
      break;
    case LIN_EULER:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_LINFLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_LinFlow, Kind_Centered_LinFlow,
                              Kind_Upwind_LinFlow, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_LinFlow);
      }
      break;
    case POISSON_EQUATION:
      if (val_system == RUNTIME_POISSON_SYS) {
        SetKind_ConvNumScheme(NONE, NONE, NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Poisson);
      }
      break;
    case WAVE_EQUATION:
      if (val_system == RUNTIME_WAVE_SYS) {
        SetKind_ConvNumScheme(NONE, NONE, NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Wave);
      }
      break;
    case HEAT_EQUATION:
      if (val_system == RUNTIME_HEAT_SYS) {
        SetKind_ConvNumScheme(NONE, NONE, NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Heat);
      }
      break;
    case LINEAR_ELASTICITY:
      if (val_system == RUNTIME_FEA_SYS) {
        SetKind_ConvNumScheme(NONE, NONE, NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_FEA);
      }
      break;
    case FLUID_STRUCTURE_EULER:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_FEA_SYS) {
        SetKind_ConvNumScheme(NONE, NONE, NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_FEA);
      }
      break;
  }
}

double* CConfig::GetPeriodicRotCenter(string val_marker) {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_RotCenter[iMarker_PerBound];
}

double* CConfig::GetPeriodicRotAngles(string val_marker) {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_RotAngles[iMarker_PerBound];
}

double* CConfig::GetPeriodicTranslation(string val_marker) {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_Translation[iMarker_PerBound];
}

unsigned short CConfig::GetMarker_Periodic_Donor(string val_marker) {
  unsigned short iMarker_PerBound, jMarker_PerBound, kMarker_All;

  /*--- Find the marker for this periodic boundary. ---*/
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;

  /*--- Find corresponding donor. ---*/
  for (jMarker_PerBound = 0; jMarker_PerBound < nMarker_PerBound; jMarker_PerBound++)
    if (Marker_PerBound[jMarker_PerBound] == Marker_PerDonor[iMarker_PerBound]) break;

  /*--- Find and return global marker index for donor boundary. ---*/
  for (kMarker_All = 0; kMarker_All < nMarker_Config; kMarker_All++)
    if (Marker_PerBound[jMarker_PerBound] == Marker_All_TagBound[kMarker_All]) break;

  return kMarker_All;
}

double* CConfig::GetActDisk_Origin(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDisk_Inlet; iMarker_ActDisk++)
    if ((Marker_ActDisk_Inlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDisk_Outlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Origin[iMarker_ActDisk];
}

double CConfig::GetActDisk_RootRadius(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDisk_Inlet; iMarker_ActDisk++)
    if ((Marker_ActDisk_Inlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDisk_Outlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_RootRadius[iMarker_ActDisk];
}

double CConfig::GetActDisk_TipRadius(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDisk_Inlet; iMarker_ActDisk++)
    if ((Marker_ActDisk_Inlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDisk_Outlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_TipRadius[iMarker_ActDisk];
}

double CConfig::GetActDisk_CT(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDisk_Inlet; iMarker_ActDisk++)
    if ((Marker_ActDisk_Inlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDisk_Outlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_CT[iMarker_ActDisk];
}

double CConfig::GetActDisk_Omega(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDisk_Inlet; iMarker_ActDisk++)
    if ((Marker_ActDisk_Inlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDisk_Outlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Omega[iMarker_ActDisk];
}

unsigned short CConfig::GetMarker_ActDisk_Outlet(string val_marker) {
  unsigned short iMarker_ActDisk, kMarker_All;

  /*--- Find the marker for this actuator disk inlet. ---*/

  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDisk_Inlet; iMarker_ActDisk++)
    if (Marker_ActDisk_Inlet[iMarker_ActDisk] == val_marker) break;

  /*--- Find and return global marker index for the actuator disk outlet. ---*/

  for (kMarker_All = 0; kMarker_All < nMarker_Config; kMarker_All++)
    if (Marker_ActDisk_Outlet[iMarker_ActDisk] == Marker_All_TagBound[kMarker_All]) break;

  return kMarker_All;
}

void CConfig::SetnPeriodicIndex(unsigned short val_index) {

  /*--- Store total number of transformations. ---*/
  nPeriodic_Index = val_index;

  /*--- Allocate memory for centers, angles, translations. ---*/
  Periodic_Center    = new double*[nPeriodic_Index];
  Periodic_Rotation  = new double*[nPeriodic_Index];
  Periodic_Translate = new double*[nPeriodic_Index];

}

unsigned short CConfig::GetMarker_Moving(string val_marker) {
  unsigned short iMarker_Moving;

  /*--- Find the marker for this moving boundary. ---*/
  for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
    if (Marker_Moving[iMarker_Moving] == val_marker) break;

  return iMarker_Moving;
}

double CConfig::GetDirichlet_Value(string val_marker) {
  unsigned short iMarker_Dirichlet;
  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet_Elec; iMarker_Dirichlet++)
    if (Marker_Dirichlet_Elec[iMarker_Dirichlet] == val_marker) break;
  return Dirichlet_Value[iMarker_Dirichlet];
}

bool CConfig::GetDirichlet_Boundary(string val_marker) {
  unsigned short iMarker_Dirichlet;
  bool Dirichlet = false;
  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet_Elec; iMarker_Dirichlet++)
    if (Marker_Dirichlet_Elec[iMarker_Dirichlet] == val_marker) {
      Dirichlet = true;
      break;
    }
  return Dirichlet;
}

double CConfig::GetNozzle_Ttotal(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Nozzle_Ttotal[iMarker_EngineExhaust];
}

double CConfig::GetNozzle_Ptotal(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Nozzle_Ptotal[iMarker_EngineExhaust];
}

double CConfig::GetInlet_Ttotal(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ttotal[iMarker_Inlet];
}

double CConfig::GetInlet_Ptotal(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ptotal[iMarker_Inlet];
}

double* CConfig::GetInlet_FlowDir(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_FlowDir[iMarker_Inlet];
}

double CConfig::GetInlet_Temperature(string val_marker) {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Temperature[iMarker_Supersonic_Inlet];
}

double CConfig::GetInlet_Pressure(string val_marker) {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Pressure[iMarker_Supersonic_Inlet];
}

double* CConfig::GetInlet_Velocity(string val_marker) {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Velocity[iMarker_Supersonic_Inlet];
}

double CConfig::GetOutlet_Pressure(string val_marker) {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
  return Outlet_Pressure[iMarker_Outlet];
}

double CConfig::GetRiemann_Var1(string val_marker) {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_Var1[iMarker_Riemann];
}

double CConfig::GetRiemann_Var2(string val_marker) {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_Var2[iMarker_Riemann];
}

double* CConfig::GetRiemann_FlowDir(string val_marker) {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_FlowDir[iMarker_Riemann];
}

unsigned short CConfig::GetKind_Data_Riemann(string val_marker) {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Kind_Data_Riemann[iMarker_Riemann];
}

double CConfig::GetIsothermal_Temperature(string val_marker) {
  unsigned short iMarker_Isothermal = 0;

  if ((nMarker_Isothermal*nMarker_IsothermalCatalytic             != 0) ||
      (nMarker_Isothermal*nMarker_IsothermalNonCatalytic          != 0) ||
      (nMarker_IsothermalCatalytic*nMarker_IsothermalNonCatalytic != 0)   ) {
    cout << "Mixed catalytic boundaries not supported!!" << endl;
    exit(EXIT_FAILURE);
  }

  if (nMarker_Isothermal > 0) {
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++)
      if (Marker_Isothermal[iMarker_Isothermal] == val_marker) break;
  }

  if (nMarker_IsothermalCatalytic > 0) {
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_IsothermalCatalytic; iMarker_Isothermal++)
      if (Marker_IsothermalCatalytic[iMarker_Isothermal] == val_marker) break;
  }

  if (nMarker_IsothermalNonCatalytic > 0) {
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_IsothermalNonCatalytic; iMarker_Isothermal++)
      if (Marker_IsothermalNonCatalytic[iMarker_Isothermal] == val_marker) break;
  }

  return Isothermal_Temperature[iMarker_Isothermal];
}

double CConfig::GetWall_HeatFlux(string val_marker) {
  unsigned short iMarker_HeatFlux = 0;

  if ((nMarker_HeatFlux*nMarker_HeatFluxCatalytic             != 0) ||
      (nMarker_HeatFlux*nMarker_HeatFluxNonCatalytic          != 0) ||
      (nMarker_HeatFluxCatalytic*nMarker_HeatFluxNonCatalytic != 0)) {
    cout << "Mixed catalytic boundaries not supported!!" << endl;
    exit(EXIT_FAILURE);
  }

  if (nMarker_HeatFlux > 0) {
  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++)
    if (Marker_HeatFlux[iMarker_HeatFlux] == val_marker) break;
  }

  if (nMarker_HeatFluxCatalytic > 0) {
    for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFluxCatalytic; iMarker_HeatFlux++)
      if (Marker_HeatFluxCatalytic[iMarker_HeatFlux] == val_marker) break;
  }

  if (nMarker_HeatFluxNonCatalytic > 0) {
    for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFluxNonCatalytic; iMarker_HeatFlux++)
      if (Marker_HeatFluxNonCatalytic[iMarker_HeatFlux] == val_marker) break;
  }

  return Heat_Flux[iMarker_HeatFlux];
}

double CConfig::GetInflow_Mach_Target(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Mach_Target[iMarker_EngineInflow];
}

double CConfig::GetBleed_MassFlow_Target(string val_marker) {
  unsigned short iMarker_EngineBleed;
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++)
    if (Marker_EngineBleed[iMarker_EngineBleed] == val_marker) break;
  return Bleed_MassFlow_Target[iMarker_EngineBleed];
}

double CConfig::GetBleed_Temperature_Target(string val_marker) {
  unsigned short iMarker_EngineBleed;
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++)
    if (Marker_EngineBleed[iMarker_EngineBleed] == val_marker) break;
  return Bleed_Temperature_Target[iMarker_EngineBleed];
}

double CConfig::GetInflow_Pressure(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Pressure[iMarker_EngineInflow];
}

double CConfig::GetBleed_Pressure(string val_marker) {
  unsigned short iMarker_EngineBleed;
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++)
    if (Marker_EngineBleed[iMarker_EngineBleed] == val_marker) break;
  return Bleed_Pressure[iMarker_EngineBleed];
}

double CConfig::GetInflow_Mach(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Mach[iMarker_EngineInflow];
}

double CConfig::GetBleed_MassFlow(string val_marker) {
  unsigned short iMarker_EngineBleed;
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++)
    if (Marker_EngineBleed[iMarker_EngineBleed] == val_marker) break;
  return Bleed_MassFlow[iMarker_EngineBleed];
}

double CConfig::GetBleed_Temperature(string val_marker) {
  unsigned short iMarker_EngineBleed;
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++)
    if (Marker_EngineBleed[iMarker_EngineBleed] == val_marker) break;
  return Bleed_Temperature[iMarker_EngineBleed];
}

double CConfig::GetDispl_Value(string val_marker) {
  unsigned short iMarker_Displacement;
  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++)
    if (Marker_Displacement[iMarker_Displacement] == val_marker) break;
  return Displ_Value[iMarker_Displacement];
}

double CConfig::GetLoad_Value(string val_marker) {
  unsigned short iMarker_Load;
  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++)
    if (Marker_Load[iMarker_Load] == val_marker) break;
  return Load_Value[iMarker_Load];
}

double CConfig::GetFlowLoad_Value(string val_marker) {
  unsigned short iMarker_FlowLoad;
  for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++)
    if (Marker_FlowLoad[iMarker_FlowLoad] == val_marker) break;
  return FlowLoad_Value[iMarker_FlowLoad];
}

void CConfig::SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2) {
  unsigned long i, k;
  double p, qn, sig, un, *u;

  u = new double [n];

  if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
    y2[0]=u[0]=0.0;			  // -ural"
  else {									// or else to have a specified first derivative.
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for (i=2; i<=n-1; i++) {									//  This is the decomposition loop of the tridiagonal al-
    sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	gorithm. y2 and u are used for tem-
    p=sig*y2[i-2]+2.0;										//	porary storage of the decomposed
    y2[i-1]=(sig-1.0)/p;										//	factors.
    u[i-1]=(y[i]-y[i-1])/(x[i]-x[i-1]) - (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
    u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
  }

  if (ypn > 0.99e30)						// The upper boundary condition is set either to be
    qn=un=0.0;									// "natural"
  else {												// or else to have a specified first derivative.
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-1; k>=1; k--)					// This is the backsubstitution loop of the tridiagonal
    y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.

  delete[] u;

}

double CConfig::GetSpline(vector<double>&xa, vector<double>&ya, vector<double>&y2a, unsigned long n, double x) {
  unsigned long klo, khi, k;
  double h, b, a, y;

  klo=1;										// We will find the right place in the table by means of
  khi=n;										// bisection. This is optimal if sequential calls to this
  while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
    k=(khi+klo) >> 1;				// are in order, and closely spaced, one would do better
    if (xa[k-1] > x) khi=k;		// to store previous values of klo and khi and test if
    else klo=k;							// they remain appropriate on the next call.
  }								// klo and khi now bracket the input value of x
  h=xa[khi-1]-xa[klo-1];
  if (h == 0.0) cout << "Bad xa input to routine splint" << endl;	// The xa’s must be dis-
  a=(xa[khi-1]-x)/h;																					      // tinct.
  b=(x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
  y=a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;

  return y;
}
