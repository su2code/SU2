/*!
 * \file config_structure.cpp
 * \brief Main file for reading the config file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#include "../include/config_structure.hpp"

CConfig::CConfig(char case_filename[200], unsigned short val_software, unsigned short val_iZone, unsigned short val_nZone, unsigned short verb_level) {
	string text_line, text2find, option_name, keyword;
	ifstream case_file;
	vector<string> option_value;
	double default_vec_3d[3];
	double default_vec_6d[6];
	nZone = val_nZone;

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Intialize motion pointers to NULL. If we don't find these values
   in the config file, they will all be set to zero. ---*/
	Motion_Origin_X = NULL;
	Motion_Origin_Y = NULL;
	Motion_Origin_Z = NULL;
	Translation_Rate_X = NULL;
	Translation_Rate_Y = NULL;
	Translation_Rate_Z = NULL;
	Rotation_Rate_X = NULL;
	Rotation_Rate_Y = NULL;
	Rotation_Rate_Z = NULL;
	Pitching_Omega_X = NULL;
	Pitching_Omega_Y = NULL;
	Pitching_Omega_Z = NULL;
	Pitching_Ampl_X = NULL;
	Pitching_Ampl_Y = NULL;
	Pitching_Ampl_Z = NULL;
	Pitching_Phase_X = NULL;
	Pitching_Phase_Y = NULL;
	Pitching_Phase_Z = NULL;
	Plunging_Omega_X = NULL;
	Plunging_Omega_Y = NULL;
	Plunging_Omega_Z = NULL;
	Plunging_Ampl_X = NULL;
	Plunging_Ampl_Y = NULL;
	Plunging_Ampl_Z = NULL;

	/* BEGIN_CONFIG_OPTIONS */

	/*--- Options related to problem definition and partitioning ---*/
	/* CONFIG_CATEGORY: Problem Definition and Partitioning */

	/* DESCRIPTION: Physical governing equations */
	AddEnumOption("PHYSICAL_PROBLEM", Kind_Solver, Solver_Map, "NONE");
	/* DESCRIPTION: Mathematical problem */
	AddMathProblem("MATH_PROBLEM" , Adjoint, false , OneShot, false, 
			Linearized, false, Restart_Flow, false);
	/* DESCRIPTION: Free-surface problem */
	AddSpecialOption("FREE_SURFACE", FreeSurface, SetBoolOption, false);
	/* DESCRIPTION: Incompressible flow using artificial compressibility */
	AddSpecialOption("INCOMPRESSIBLE_FORMULATION", Incompressible, SetBoolOption, false);
	/* DESCRIPTION: Axisymmetric simulation */
	AddSpecialOption("AXISYMMETRIC", Axisymmetric, SetBoolOption, false);
	/* DESCRIPTION: Magnetic simulation */
	AddSpecialOption("MAGNET", MagneticForce, SetBoolOption, false);
	/* DESCRIPTION: Joule heating simulation */
	AddSpecialOption("JOULE_HEAT", JouleHeating, SetBoolOption, false);
	/* DESCRIPTION: Flag for running the electric potential solver as part of the plasma solver */
	AddSpecialOption("ELECTRIC_SOLVER", ElectricSolver, SetBoolOption, false);
	/* DESCRIPTION:  */
	AddSpecialOption("MACCORMACK_RELAXATION", MacCormackRelaxation, SetBoolOption, false);
	/* DESCRIPTION: Time stepping of the various species in a steady plasma solution */
	AddSpecialOption("PLASMA_MULTI_TIME_STEP", PlasmaMultiTimeSteps, SetBoolOption, false);
	/* DESCRIPTION: Perform a low fidelity simulation */
	AddSpecialOption("LOW_FIDELITY_SIMULATION", LowFidelitySim, SetBoolOption, false);

	/* DESCRIPTION: Roe-Turkel preconditioning for low Mach number flows */
	AddSpecialOption("ROE_TURKEL_PREC", Low_Mach_Precon, SetBoolOption, false);
	/* DESCRIPTION: Time Step for dual time stepping simulations (s) */
	AddScalarOption("MIN_ROE_TURKEL_PREC", Min_Beta_RoeTurkel, 0.01);
	/* DESCRIPTION: Time Step for dual time stepping simulations (s) */
	AddScalarOption("MAX_ROE_TURKEL_PREC", Max_Beta_RoeTurkel, 0.2);
	/* DESCRIPTION: Time Step for dual time stepping simulations (s) */
	AddScalarOption("STAGNATION_BFIELD", Stagnation_B, 0.2);
	/* DESCRIPTION: Time Step for dual time stepping simulations (s) */
	AddScalarOption("ELECTRICAL_CONDUCTIVITY", Electric_Cond, 2000.0);
	/* DESCRIPTION: Time Step for dual time stepping simulations (s) */
	AddScalarOption("DIPOLE_DIST", DipoleDist, 1E-6);
	/* DESCRIPTION: Restart solution from native solution file */
	AddSpecialOption("RESTART_SOL", Restart, SetBoolOption, false);
	/* DESCRIPTION: Restart a Plasma solution from an Euler native solution file */
	AddSpecialOption("RESTART_PLASMA_FROM_EULER", Restart_Euler2Plasma, SetBoolOption, false);
	/* DESCRIPTION: Write a tecplot file for each partition */
	AddSpecialOption("VISUALIZE_PART", Visualize_Partition, SetBoolOption, false);
	/* DESCRIPTION: Divide rectangles into triangles */
	AddSpecialOption("DIVIDE_ELEMENTS", Divide_Element, SetBoolOption, false);
	/* DESCRIPTION: Specify turbulence model */
	AddEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, "NONE");
	/* DESCRIPTION: Specify transition model */
	AddEnumOption("KIND_TRANS_MODEL", Kind_Trans_Model, Trans_Model_Map, "NONE");
	/* DESCRIPTION: Engine subsonic intake region */
	AddSpecialOption("SUBSONIC_NACELLE_INFLOW", Engine_Intake, SetBoolOption, false);
  
	/*--- options related to various boundary markers ---*/
	/* CONFIG_CATEGORY: Boundary Markers */

	/* DESCRIPTION: Marker(s) of the surface in the surface flow solution file */
	AddMarkerOption("MARKER_PLOTTING", nMarker_Plotting, Marker_Plotting);
	/* DESCRIPTION: Marker(s) of the surface where evaluate the non-dimensional coefficients */
	AddMarkerOption("MARKER_MONITORING", nMarker_Monitoring, Marker_Monitoring);
  /* DESCRIPTION: Marker(s) of the surface where objective function (design problem) will be evaluated */
	AddMarkerOption("MARKER_DESIGNING", nMarker_Designing, Marker_Designing);
  
	/* DESCRIPTION: Euler wall boundary marker(s) */
	AddMarkerOption("MARKER_EULER", nMarker_Euler, Marker_Euler);
	/* DESCRIPTION: Adiabatic wall boundary condition */
	AddSpecialOption("ADIABATIC_WALL", AdiabaticWall, SetBoolOption, true);
	/* DESCRIPTION: Isothermal wall boundary condition */
	AddSpecialOption("ISOTHERMAL_WALL", IsothermalWall, SetBoolOption, false);
	/* DESCRIPTION: Catalytic wall boundary condition */
	AddSpecialOption("CATALYTIC_WALL", CatalyticWall, SetBoolOption, false);
	/* DESCRIPTION: Far-field boundary marker(s) */
	AddMarkerOption("MARKER_FAR", nMarker_FarField, Marker_FarField);
	/* DESCRIPTION: Symmetry boundary marker(s) */
	AddMarkerOption("MARKER_SYM", nMarker_SymWall, Marker_SymWall);
	/* DESCRIPTION: Near-Field boundary marker(s) */
	AddMarkerOption("MARKER_NEARFIELD", nMarker_NearFieldBound, Marker_NearFieldBound);
	/* DESCRIPTION: Zone interface boundary marker(s) */
	AddMarkerOption("MARKER_INTERFACE", nMarker_InterfaceBound, Marker_InterfaceBound);
	/* DESCRIPTION: Dirichlet boundary marker(s) */
	AddMarkerOption("MARKER_DIRICHLET", nMarker_Dirichlet, Marker_Dirichlet);
	/* DESCRIPTION: Custom boundary marker(s) */
	AddMarkerOption("MARKER_CUSTOM", nMarker_Custom, Marker_Custom);
	/* DESCRIPTION: Neumann boundary marker(s) */
	AddMarkerOption("MARKER_NEUMANN", nMarker_Neumann, Marker_Neumann);


	/* DESCRIPTION: Electric dirichlet boundary marker(s) */
	AddMarkerDirichlet("ELEC_DIRICHLET", nMarker_Dirichlet_Elec, Marker_Dirichlet_Elec, Dirichlet_Value );
	/* DESCRIPTION: Electric neumann boundary marker(s) */
	AddMarkerOption("ELEC_NEUMANN", nMarker_Neumann_Elec, Marker_Neumann_Elec);


	/* DESCRIPTION: Periodic boundary marker(s) for use with SU2_PBC
     Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
     rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
     rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) */
	AddMarkerPeriodic("MARKER_PERIODIC", nMarker_PerBound, Marker_PerBound, Marker_PerDonor,
			Periodic_RotCenter, Periodic_RotAngles, Periodic_Translation);
	/* DESCRIPTION: Sliding mesh interface boundary marker(s) for use with SU2_SMC
     Format: ( sliding marker, zone # of sliding marker, donor marker, zone # of donor, ... ) */
	AddMarkerSliding("MARKER_SLIDING", nMarker_Sliding, Marker_SlideBound, Marker_SlideDonor,
			SlideBound_Zone, SlideDonor_Zone);
	/* DESCRIPTION: Inlet boundary type */
	AddEnumOption("INLET_TYPE", Kind_Inlet, Inlet_Map, "TOTAL_CONDITIONS");
	/* DESCRIPTION: Inlet boundary marker(s) with the following formats,
     Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x,
               flow_direction_y, flow_direction_z, ... ) where flow_direction is
               a unit vector.
     Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x,
               flow_direction_y, flow_direction_z, ... ) where flow_direction is
               a unit vector. */
	AddMarkerInlet("MARKER_INLET", nMarker_Inlet, Marker_Inlet, Inlet_Ttotal, Inlet_Ptotal, Inlet_FlowDir);
	/* DESCRIPTION: % Supersonic inlet boundary marker(s)
   Format: (inlet marker, temperature, static pressure, velocity_x,
   velocity_y, velocity_z, ... ), i.e. primitive variables specified. */
	AddMarkerInlet("MARKER_SUPERSONIC_INLET", nMarker_Supersonic_Inlet, Marker_Supersonic_Inlet,
			Inlet_Temperature, Inlet_Pressure, Inlet_Velocity);
	/* DESCRIPTION: Outlet boundary marker(s)
     Format: ( outlet marker, back pressure (static), ... ) */
	AddMarkerOutlet("MARKER_OUTLET", nMarker_Outlet, Marker_Outlet, Outlet_Pressure);
	/* DESCRIPTION: Isothermal wall boundary marker(s)
   Format: ( isothermal marker, wall temperature (static), ... ) */
	AddMarkerOutlet("MARKER_ISOTHERMAL", nMarker_Isothermal, Marker_Isothermal, Isothermal_Temperature);
	/* DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) */
	AddMarkerOutlet("MARKER_HEATFLUX", nMarker_HeatFlux, Marker_HeatFlux, Heat_Flux);
	/* DESCRIPTION: Nacelle inflow boundary marker(s)
     Format: ( nacelle inflow marker, fan face Mach, ... ) */
	AddMarkerOutlet("MARKER_NACELLE_INFLOW", nMarker_NacelleInflow, Marker_NacelleInflow, FanFace_Mach_Target);
	/* DESCRIPTION: Nacelle exhaust boundary marker(s)
     Format: (nacelle exhaust marker, total nozzle temp, total nozzle pressure, ... )*/
	AddMarkerInlet("MARKER_NACELLE_EXHAUST", nMarker_NacelleExhaust, Marker_NacelleExhaust, Nozzle_Ttotal, Nozzle_Ptotal);
	/* DESCRIPTION: Displacement boundary marker(s) */
	AddMarkerDisplacement("MARKER_DISPL", nMarker_Displacement, Marker_Displacement, Displ_Value);
	/* DESCRIPTION: Load boundary marker(s) */
	AddMarkerLoad("MARKER_LOAD", nMarker_Load, Marker_Load, Load_Value);
	/* DESCRIPTION: Flow load boundary marker(s) */
	AddMarkerFlowLoad("MARKER_FLOWLOAD", nMarker_FlowLoad, Marker_FlowLoad, FlowLoad_Value);
	/* DESCRIPTION: FW-H boundary marker(s) */
	AddMarkerOption("MARKER_FWH", nMarker_FWH, Marker_FWH);
	/* DESCRIPTION: Observer boundary marker(s) */
	AddMarkerOption("MARKER_OBSERVER", nMarker_Observer, Marker_Observer);
	/* DESCRIPTION: Damping factor for engine inlet condition */
	AddScalarOption("DAMP_NACELLE_INFLOW", Damp_Nacelle_Inflow, 0.1);
    
	/*--- options related to grid adaptation ---*/
	/* CONFIG_CATEGORY: Grid adaptation */

	/* DESCRIPTION: Kind of grid adaptation */
	AddEnumOption("KIND_ADAPT", Kind_Adaptation, Adapt_Map, "NONE");
	/* DESCRIPTION: Percentage of new elements (% of the original number of elements) */
	AddScalarOption("NEW_ELEMS", New_Elem_Adapt, -1.0);
	/* DESCRIPTION: Scale factor for the dual volume */
	AddScalarOption("DUALVOL_POWER", DualVol_Power, 0.5);
	/* DESCRIPTION: Use analytical definition for surfaces */
	AddEnumOption("ANALYTICAL_SURFDEF", Analytical_Surface, Geo_Analytic_Map, "NONE");
	/* DESCRIPTION: Before each computation, implicitly smooth the nodal coordinates */
	AddSpecialOption("SMOOTH_GEOMETRY", SmoothNumGrid, SetBoolOption, false);
	/* DESCRIPTION: Adapt the boundary elements */
	AddSpecialOption("ADAPT_BOUNDARY", AdaptBoundary, SetBoolOption, true);


	/*--- options related to time-marching ---*/
	/* CONFIG_CATEGORY: Time-marching */

	/* DESCRIPTION: Unsteady simulation  */
	AddEnumOption("UNSTEADY_SIMULATION", Unsteady_Simulation, Unsteady_Map, "NO");
	/* DESCRIPTION:  Unsteady farfield boundaries  */
	AddSpecialOption("UNSTEADY_FARFIELD", Unsteady_Farfield, SetBoolOption, false);
	/* DESCRIPTION:  Courant-Friedrichs-Lewy condition of the finest grid */
	AddScalarOption("CFL_NUMBER", CFLFineGrid, 1.25);
	/* DESCRIPTION: CFL ramp (factor, number of iterations, CFL limit) */
	default_vec_3d[0] = 1.0; default_vec_3d[1] = 100.0; default_vec_3d[2] = 1.0;
	AddArrayOption("CFL_RAMP", 3, CFLRamp, default_vec_3d);
	/* DESCRIPTION: Number of total iterations */
	AddScalarOption("EXT_ITER", nExtIter, 999999);
	// these options share nRKStep as their size, which is not a good idea in general
	/* DESCRIPTION: Runge-Kutta alpha coefficients */
	AddListOption("RK_ALPHA_COEFF", nRKStep, RK_Alpha_Step);
	/* DESCRIPTION: Time Step for dual time stepping simulations (s) */
	AddScalarOption("UNST_TIMESTEP", Delta_UnstTime, 0.0);
	/* DESCRIPTION: Total Physical Time for dual time stepping simulations (s) */
	AddScalarOption("UNST_TIME", Total_UnstTime, 1.0);
	/* DESCRIPTION: Unsteady Courant-Friedrichs-Lewy number of the finest grid */
	AddScalarOption("UNST_CFL_NUMBER", Unst_CFL, 0.0);
	/* DESCRIPTION: Number of internal iterations (dual time method) */
	AddScalarOption("UNST_INT_ITER", Unst_nIntIter, 100);
	/* DESCRIPTION: Integer number of periodic time instances for Time Spectral */
	AddScalarOption("TIME_INSTANCES", nTimeInstances, 1);
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_FLOW", Kind_TimeIntScheme_Flow, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_LEVELSET", Kind_TimeIntScheme_LevelSet, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_ADJLEVELSET", Kind_TimeIntScheme_AdjLevelSet, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_PLASMA", Kind_TimeIntScheme_Plasma, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_ADJPLASMA", Kind_TimeIntScheme_AdjPlasma, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_ADJ", Kind_TimeIntScheme_AdjFlow, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_LIN", Kind_TimeIntScheme_LinFlow, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_TURB", Kind_TimeIntScheme_Turb, Time_Int_Map, "EULER_IMPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_ADJTURB", Kind_TimeIntScheme_AdjTurb, Time_Int_Map, "EULER_IMPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_WAVE", Kind_TimeIntScheme_Wave, Time_Int_Map, "EULER_IMPLICIT");
	/* DESCRIPTION: Time discretization */
	AddEnumOption("TIME_DISCRE_FEA", Kind_TimeIntScheme_FEA, Time_Int_Map, "EULER_IMPLICIT");
	/* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
	AddScalarOption("LEVELSET_CFL_REDUCTION", LevelSet_CFLRedCoeff, 1E-2);
	/* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
	AddScalarOption("TURB_CFL_REDUCTION", Turb_CFLRedCoeff, 1.0);
  
	/*--- options related to the linear solvers ---*/
	/* CONFIG_CATEGORY: Linear solver definition */

	/* DESCRIPTION: Linear solver for the implicit, mesh deformation, or discrete adjoint systems */
	AddEnumOption("LINEAR_SOLVER", Kind_Linear_Solver, Linear_Solver_Map, "FGMRES");
	/* DESCRIPTION: Preconditioner for the Krylov linear solvers */
	AddEnumOption("LINEAR_SOLVER_PREC", Kind_Linear_Solver_Prec, Linear_Solver_Prec_Map, "LU_SGS");
	/* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
	AddScalarOption("LINEAR_SOLVER_ERROR", Linear_Solver_Error, 1E-5);
	/* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
	AddScalarOption("LINEAR_SOLVER_ITER", Linear_Solver_Iter, 10);
	/* DESCRIPTION: Relaxation of the linear solver for the implicit formulation */
	AddScalarOption("LINEAR_SOLVER_RELAX", Linear_Solver_Relax, 1.0);

	/* DESCRIPTION: Linear solver for the turbulent adjoint systems */
	AddEnumOption("ADJTURB_LIN_SOLVER", Kind_AdjTurb_Linear_Solver, Linear_Solver_Map, "FGMRES");
	/* DESCRIPTION: Preconditioner for the turbulent adjoint Krylov linear solvers */
	AddEnumOption("ADJTURB_LIN_PREC", Kind_AdjTurb_Linear_Prec, Linear_Solver_Prec_Map, "LU_SGS");
	/* DESCRIPTION: Minimum error threshold for the turbulent adjoint linear solver for the implicit formulation */
	AddScalarOption("ADJTURB_LIN_ERROR", AdjTurb_Linear_Error, 1E-5);
	/* DESCRIPTION: Maximum number of iterations of the turbulent adjoint linear solver for the implicit formulation */
	AddScalarOption("ADJTURB_LIN_ITER", AdjTurb_Linear_Iter, 10);

	/*--- options related to dynamic meshes ---*/
	/* CONFIG_CATEGORY: Dynamic mesh definition */

	/* DESCRIPTION: Mesh motion for unsteady simulations */
	AddSpecialOption("GRID_MOVEMENT", Grid_Movement, SetBoolOption, false);
  /* DESCRIPTION: Add the gravity force */
	AddSpecialOption("GRAVITY_FORCE", GravityForce, SetBoolOption, false);
	/* DESCRIPTION: Type of mesh motion */
	AddEnumListOption("GRID_MOVEMENT_KIND", nZone, Kind_GridMovement, GridMovement_Map);
	default_vec_3d[0] = 0; default_vec_3d[1] = 0;
	/* DESCRIPTION: Mesh motion ramp (starting iteration, ramp iterations (zero to max)) */
	AddArrayOption("MOTION_RAMP", 2, Motion_Ramp, default_vec_3d);
	/* DESCRIPTION: % Mach number (non-dimensional, based on the mesh velocity and freestream vals.) */
	AddScalarOption("MACH_MOTION", Mach_Motion, 0.0);
	/* DESCRIPTION: Coordinates of the rigid motion origin */
	AddListOption("MOTION_ORIGIN_X", nZone, Motion_Origin_X);
	/* DESCRIPTION: Coordinates of the rigid motion origin */
	AddListOption("MOTION_ORIGIN_Y", nZone, Motion_Origin_Y);
	/* DESCRIPTION: Coordinates of the rigid motion origin */
	AddListOption("MOTION_ORIGIN_Z", nZone, Motion_Origin_Z);
	/* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
	AddListOption("TRANSLATION_RATE_X", nZone, Translation_Rate_X);
	/* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
	AddListOption("TRANSLATION_RATE_Y", nZone, Translation_Rate_Y);
	/* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
	AddListOption("TRANSLATION_RATE_Z", nZone, Translation_Rate_Z);
	/* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("ROTATION_RATE_X", nZone, Rotation_Rate_X);
	/* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("ROTATION_RATE_Y", nZone, Rotation_Rate_Y);
	/* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("ROTATION_RATE_Z", nZone, Rotation_Rate_Z);
	/* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_OMEGA_X", nZone, Pitching_Omega_X);
	/* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_OMEGA_Y", nZone, Pitching_Omega_Y);
	/* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_OMEGA_Z", nZone, Pitching_Omega_Z);
	/* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_AMPL_X", nZone, Pitching_Ampl_X);
	/* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_AMPL_Y", nZone, Pitching_Ampl_Y);
	/* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_AMPL_Z", nZone, Pitching_Ampl_Z);
	/* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_PHASE_X", nZone, Pitching_Phase_X);
	/* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_PHASE_Y", nZone, Pitching_Phase_Y);
	/* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
	AddListOption("PITCHING_PHASE_Z", nZone, Pitching_Phase_Z);
	/* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
	AddListOption("PLUNGING_OMEGA_X", nZone, Plunging_Omega_X);
	/* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
	AddListOption("PLUNGING_OMEGA_Y", nZone, Plunging_Omega_Y);
	/* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
	AddListOption("PLUNGING_OMEGA_Z", nZone, Plunging_Omega_Z);
	/* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
	AddListOption("PLUNGING_AMPL_X", nZone, Plunging_Ampl_X);
	/* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
	AddListOption("PLUNGING_AMPL_Y", nZone, Plunging_Ampl_Y);
	/* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
	AddListOption("PLUNGING_AMPL_Z", nZone, Plunging_Ampl_Z);
	/* DESCRIPTION:  */
	AddScalarOption("MOTION_FILENAME", Motion_Filename, string("mesh_motion.dat"));
	/* DESCRIPTION: Reduced frequency for flutter pitching motion */
	AddScalarOption("RED_FREC", Reduced_Frequency, 0.0);
	/* DESCRIPTION: Pitching amplitude for flutter (degrees) */
	AddScalarOption("MAX_PITCH", Pitching_Amplitude, 0.0);
	/* DESCRIPTION: Uncoupled Aeroelastic Frequency Plunge. */
	AddScalarOption("FREQ_PLUNGE_AEROELASTIC", FreqPlungeAeroelastic, 100);
	/* DESCRIPTION: Uncoupled Aeroelastic Frequency Pitch. */
	AddScalarOption("FREQ_PITCH_AEROELASTIC", FreqPitchAeroelastic, 100);
	/* DESCRIPTION: Type of aeroelastic grid movement */
	AddEnumOption("TYPE_AEROELASTIC_GRID_MOVEMENT", Aeroelastic_Grid_Movement, Aeroelastic_Movement_Map, "RIGID");
	/* DESCRIPTION: Type of grid velocities for aeroelastic motion */
	AddEnumOption("TYPE_AEROELASTIC_GRID_VELOCITY", Aeroelastic_Grid_Velocity, Aeroelastic_Velocity_Map, "FD");

	/*--- options related to rotating frame problems ---*/
	/* CONFIG_CATEGORY: Rotating frame */

	/* DESCRIPTION: Rotating frame problem */
	AddSpecialOption("ROTATING_FRAME", Rotating_Frame, SetBoolOption, false);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	/* DESCRIPTION: Origin of the axis of rotation */
	AddArrayOption("ROTATIONAL_ORIGIN", 3, RotAxisOrigin, default_vec_3d);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	/* DESCRIPTION: Angular velocity vector (rad/s) */
	AddArrayOption("ROTATION_RATE", 3, Omega, default_vec_3d);
	/*--- Initializing this here because it might be needed in the geometry classes. ---*/
	Omega_FreeStreamND = new double[3];

	/*--- options related to convergence ---*/
	/* CONFIG_CATEGORY: Convergence*/

	/* DESCRIPTION: Convergence criteria */
	AddEnumOption("CONV_CRITERIA", ConvCriteria, Converge_Crit_Map, "RESIDUAL");

	/* DESCRIPTION: Residual reduction (order of magnitude with respect to the initial value) */
	AddScalarOption("RESIDUAL_REDUCTION", OrderMagResidual, 3.0);
	/* DESCRIPTION: Min value of the residual (log10 of the residual) */
	AddScalarOption("RESIDUAL_MINVAL", MinLogResidual, -8.0);
	/* DESCRIPTION: Iteration number to begin convergence monitoring */
	AddScalarOption("STARTCONV_ITER", StartConv_Iter, 5);
	/* DESCRIPTION: Number of elements to apply the criteria */
	AddScalarOption("CAUCHY_ELEMS", Cauchy_Elems, 100);
	/* DESCRIPTION: Epsilon to control the series convergence */
	AddScalarOption("CAUCHY_EPS", Cauchy_Eps, 1E-10);
	/* DESCRIPTION: Flow functional for the Cauchy criteria */
	AddEnumOption("CAUCHY_FUNC_FLOW", Cauchy_Func_Flow, Objective_Map, "DRAG");
	/* DESCRIPTION: Adjoint functional for the Cauchy criteria */
	AddEnumOption("CAUCHY_FUNC_ADJ", Cauchy_Func_AdjFlow, Sens_Map, "SENS_GEOMETRY");
	/* DESCRIPTION: Linearized functional for the Cauchy criteria */
	AddEnumOption("CAUCHY_FUNC_LIN", Cauchy_Func_LinFlow, Linear_Obj_Map, "DELTA_DRAG");
	/* DESCRIPTION: Epsilon for a full multigrid method evaluation */
	AddScalarOption("FULLMG_CAUCHY_EPS", Cauchy_Eps_FullMG, 1E-4);

	/*--- options related to Multi-grid ---*/
	/* CONFIG_CATEGORY: Multi-grid */

	/* DESCRIPTION: Full multi-grid  */
	AddSpecialOption("FULLMG", FullMG, SetBoolOption, false);
	/* DESCRIPTION: Start up iterations using the fine grid only */
	AddScalarOption("START_UP_ITER", nStartUpIter, 0);
	/* DESCRIPTION: Multi-grid Levels */
	AddScalarOption("MGLEVEL", nMultiLevel, 0);
	/* DESCRIPTION: Multi-grid Cycle (0 = V cycle, 1 = W Cycle) */
	AddScalarOption("MGCYCLE", MGCycle, 0);
	/* DESCRIPTION: Multi-grid pre-smoothing level */
	AddListOption("MG_PRE_SMOOTH", nMG_PreSmooth, MG_PreSmooth);
	/* DESCRIPTION: Multi-grid post-smoothing level */
	AddListOption("MG_POST_SMOOTH", nMG_PostSmooth, MG_PostSmooth);
	/* DESCRIPTION: Jacobi implicit smoothing of the correction */
	AddListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth, MG_CorrecSmooth);
	/* DESCRIPTION: Damping factor for the residual restriction */
	AddScalarOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.75);
	/* DESCRIPTION: Damping factor for the correction prolongation */
	AddScalarOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.75);
	/* DESCRIPTION: CFL reduction factor on the coarse levels */
	AddScalarOption("MG_CFL_REDUCTION", MG_CFLRedCoeff, 0.75);
	/* DESCRIPTION: Maximum number of children in the agglomeration stage */
	AddScalarOption("MAX_CHILDREN", MaxChildren, 500);
	/* DESCRIPTION: Maximum length of an agglomerated element (relative to the domain) */
	AddScalarOption("MAX_DIMENSION", MaxDimension, 0.1);

	/*--- options related to the discretization ---*/
	/* CONFIG_CATEGORY: Discretization */

	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_FLOW", Kind_SlopeLimit_Flow, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_ADJFLOW", Kind_SlopeLimit_AdjFlow, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_TURB", Kind_SlopeLimit_Turb, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_ADJTURB", Kind_SlopeLimit_AdjTurb, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_LEVELSET", Kind_SlopeLimit_LevelSet, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_ADJLEVELSET", Kind_SlopeLimit_AdjLevelSet, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_PLASMA", Kind_SlopeLimit_Plasma, Limiter_Map, "NONE");
	/* DESCRIPTION: Slope limiter */
	AddEnumOption("SLOPE_LIMITER_ADJPLASMA", Kind_SlopeLimit_AdjPlasma, Limiter_Map, "NONE");
	/* DESCRIPTION: Coefficient for the limiter */
	AddScalarOption("LIMITER_COEFF", LimiterCoeff, 0.3);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centered_Flow, Kind_Upwind_Flow);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_LEVELSET", Kind_ConvNumScheme_LevelSet, Kind_Centered_LevelSet, Kind_Upwind_LevelSet);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_ADJLEVELSET", Kind_ConvNumScheme_AdjLevelSet, Kind_Centered_AdjLevelSet, Kind_Upwind_AdjLevelSet);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_PLASMA", Kind_ConvNumScheme_Plasma, Kind_Centered_Plasma, Kind_Upwind_Plasma);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_ADJPLASMA", Kind_ConvNumScheme_AdjPlasma, Kind_Centered_AdjPlasma, Kind_Upwind_AdjPlasma);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_ADJ", Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow, Kind_Upwind_AdjFlow);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_LIN", Kind_ConvNumScheme_LinFlow, Kind_Centered_LinFlow, Kind_Upwind_LinFlow);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_TURB", Kind_ConvNumScheme_Turb, Kind_Centered_Turb, Kind_Upwind_Turb);
	/* DESCRIPTION: Convective numerical method */
	AddConvectOption("CONV_NUM_METHOD_ADJTURB", Kind_ConvNumScheme_AdjTurb, Kind_Centered_AdjTurb, Kind_Upwind_AdjTurb);
	/* DESCRIPTION: Numerical method for spatial gradients */
	AddEnumOption("NUM_METHOD_GRAD", Kind_Gradient_Method, Gradient_Map, "GREEN_GAUSS");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_FLOW", Kind_ViscNumScheme_Flow, Viscous_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_WAVE", Kind_ViscNumScheme_Wave, Viscous_Map, "GALERKIN");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_FEA", Kind_ViscNumScheme_FEA, Viscous_Map, "GALERKIN");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_WAVE", Kind_SourNumScheme_Wave, Source_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_FEA", Kind_SourNumScheme_FEA, Source_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_FLOW", Kind_SourNumScheme_Flow, Source_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_LEVELSET", Kind_SourNumScheme_LevelSet, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_ADJLEVELSET", Kind_ViscNumScheme_AdjLevelSet, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_ADJLEVELSET", Kind_SourNumScheme_AdjLevelSet, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_ADJ", Kind_ViscNumScheme_AdjFlow, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_ADJ", Kind_SourNumScheme_AdjFlow, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_LIN", Kind_ViscNumScheme_LinFlow, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_LIN", Kind_SourNumScheme_LinFlow, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_PLASMA", Kind_ViscNumScheme_Plasma, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_PLASMA", Kind_SourNumScheme_Plasma, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_ADJPLASMA", Kind_ViscNumScheme_AdjPlasma, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_ADJPLASMA", Kind_SourNumScheme_AdjPlasma, Source_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_JAC_METHOD_PLASMA", Kind_SourJac_Plasma, SourceJac_Map, "FINITE_DIFF");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_TEMPLATE", Kind_SourNumScheme_Template, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_TURB", Kind_ViscNumScheme_Turb, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_TURB", Kind_SourNumScheme_Turb, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_ELEC", Kind_ViscNumScheme_Elec, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_ELEC", Kind_SourNumScheme_Elec, Source_Map, "NONE");
	/* DESCRIPTION: Viscous numerical method */
	AddEnumOption("VISC_NUM_METHOD_ADJTURB", Kind_ViscNumScheme_AdjTurb, Viscous_Map, "NONE");
	/* DESCRIPTION: Source term numerical method */
	AddEnumOption("SOUR_NUM_METHOD_ADJTURB", Kind_SourNumScheme_AdjTurb, Source_Map, "NONE");

	/*--- options related to the artificial dissipation ---*/
	/* CONFIG_CATEGORY: Artificial dissipation */

	default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
	/* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
	AddArrayOption("AD_COEFF_FLOW", 3, Kappa_Flow, default_vec_3d);
	default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
	/* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
	AddArrayOption("AD_COEFF_ADJ", 3, Kappa_AdjFlow, default_vec_3d);
	default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
	/* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
	AddArrayOption("AD_COEFF_PLASMA", 3, Kappa_Plasma, default_vec_3d);
	default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
	/* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
	AddArrayOption("AD_COEFF_ADJPLASMA", 3, Kappa_AdjPlasma, default_vec_3d);
	default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.02;
	/* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
	AddArrayOption("AD_COEFF_LIN", 2, Kappa_LinFlow, default_vec_3d);


	/*--- options related to the adjoint and gradient ---*/
	/* CONFIG_CATEGORY: Adjoint and gradient */

	/* DESCRIPTION: Adjoint type */
	AddEnumOption("ADJOINT_TYPE", Kind_Adjoint, Adjoint_Map,"CONTINUOUS");
	/* DESCRIPTION: Reduction factor of the CFL coefficient in the adjoint problem */
	AddScalarOption("ADJ_CFL_REDUCTION", Adj_CFLRedCoeff, 0.8);
  /* DESCRIPTION: Limit value for the adjoint variable */
	AddScalarOption("ADJ_LIMIT", AdjointLimit, 1E6);
	/* DESCRIPTION: Adjoint problem boundary condition */
	AddEnumOption("ADJ_OBJFUNC", Kind_ObjFunc, Objective_Map, "DRAG");
	/* DESCRIPTION: Geometrical objective function */
	AddEnumOption("GEO_PARAM", Kind_GeoObjFunc, Objective_Map, "MAX_THICKNESS");
	/* DESCRIPTION: Mode of the GDC code (analysis, or gradient) */
	AddEnumOption("GEO_MODE", GeometryMode, GeometryMode_Map, "ANALYSIS");
	/* DESCRIPTION: Drag weight in sonic boom Objective Function (from 0.0 to 1.0) */
	AddScalarOption("DRAG_IN_SONICBOOM", WeightCd, 0.0);
	/* DESCRIPTION: Sensitivity smoothing  */
	AddEnumOption("SENS_SMOOTHING", Kind_SensSmooth, Sens_Smoothing_Map, "NONE");
	/* DESCRIPTION: Continuous governing equation set  */
	AddEnumOption("CONTINUOUS_EQNS", Continuous_Eqns, ContinuousEqns_Map, "EULER");
	/* DESCRIPTION: Discrete governing equation set */
	AddEnumOption("DISCRETE_EQNS", Discrete_Eqns, DiscreteEqns_Map, "NONE");
	/* DESCRIPTION: Adjoint frozen viscosity */
	AddSpecialOption("FROZEN_VISC", Frozen_Visc, SetBoolOption, true);
	/* DESCRIPTION:  */
	AddScalarOption("CTE_VISCOUS_DRAG", CteViscDrag, 0.0);
	/* DESCRIPTION: Print sensitivities to screen on exit */
	AddSpecialOption("SHOW_ADJ_SENS", Show_Adj_Sens, SetBoolOption, false);

	/* DESCRIPTION: Reduction factor of the CFL coefficient in the turbulent adjoint problem */
	AddScalarOption("ADJTURB_CFL_REDUCTION", AdjTurb_CFLRedCoeff, 0.1);

	/*--- options related to input/output files and formats ---*/
	/* CONFIG_CATEGORY: Input/output files and formats */

	/* DESCRIPTION: Output file format */
	AddEnumOption("OUTPUT_FORMAT", Output_FileFormat, Output_Map, "TECPLOT");
	/* DESCRIPTION: Mesh input file format */
	AddEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, "SU2");
	/* DESCRIPTION: Convert a CGNS mesh to SU2 format */
	AddSpecialOption("CGNS_TO_SU2", CGNS_To_SU2, SetBoolOption, false);
	/* DESCRIPTION:  Mesh input file */
	AddScalarOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));
	/* DESCRIPTION: Mesh output file */
	AddScalarOption("MESH_OUT_FILENAME", Mesh_Out_FileName, string("mesh_out.su2"));
	/* DESCRIPTION: Output file convergence history (w/o extension) */
	AddScalarOption("CONV_FILENAME", Conv_FileName, string("history"));
	/* DESCRIPTION: Restart flow input file */
	AddScalarOption("SOLUTION_FLOW_FILENAME", Solution_FlowFileName, string("solution_flow.dat"));
	/* DESCRIPTION: Restart flow input file */
	AddScalarOption("FARFIELD_FILENAME", Farfield_FileName, string("farfield.dat"));
	/* DESCRIPTION: Restart linear flow input file */
	AddScalarOption("SOLUTION_LIN_FILENAME", Solution_LinFileName, string("solution_lin.dat"));
	/* DESCRIPTION: Restart adjoint input file */
	AddScalarOption("SOLUTION_ADJ_FILENAME", Solution_AdjFileName, string("solution_adj.dat"));
	/* DESCRIPTION: Output file restart flow */
	AddScalarOption("RESTART_FLOW_FILENAME", Restart_FlowFileName, string("restart_flow.dat"));
	/* DESCRIPTION: Output file linear flow */
	AddScalarOption("RESTART_LIN_FILENAME",Restart_LinFileName, string("restart_lin.dat"));
	/* DESCRIPTION: Output file restart adjoint */
	AddScalarOption("RESTART_ADJ_FILENAME", Restart_AdjFileName, string("restart_adj.dat"));
	/* DESCRIPTION: Output file restart wave */
	AddScalarOption("RESTART_WAVE_FILENAME", Restart_WaveFileName, string("restart_wave.dat"));
	/* DESCRIPTION: Output file flow (w/o extension) variables */
	AddScalarOption("VOLUME_FLOW_FILENAME", Flow_FileName, string("flow"));
	/* DESCRIPTION: Output file structure (w/o extension) variables */
	AddScalarOption("VOLUME_STRUCTURE_FILENAME", Structure_FileName, string("structure"));
	/* DESCRIPTION: Output file wave (w/o extension) variables */
	AddScalarOption("VOLUME_WAVE_FILENAME", Wave_FileName, string("wave"));
	/* DESCRIPTION: Output file adj. wave (w/o extension) variables */
	AddScalarOption("VOLUME_ADJWAVE_FILENAME", AdjWave_FileName, string("adjoint_wave"));
	/* DESCRIPTION: Output file adjoint (w/o extension) variables */
	AddScalarOption("VOLUME_ADJ_FILENAME", Adj_FileName, string("adjoint"));
	/* DESCRIPTION: Output file linear (w/o extension) variables */
	AddScalarOption("VOLUME_LIN_FILENAME", Lin_FileName, string("linearized"));
	/* DESCRIPTION: Output objective function gradient */
	AddScalarOption("GRAD_OBJFUNC_FILENAME", ObjFunc_Grad_FileName, string("of_grad.dat"));
	/* DESCRIPTION: Output objective function */
	AddScalarOption("OBJFUNC_FILENAME", ObjFunc_Eval_FileName, string("of_eval.dat"));
	/* DESCRIPTION: Output file surface flow coefficient (w/o extension) */
	AddScalarOption("SURFACE_FLOW_FILENAME", SurfFlowCoeff_FileName, string("surface_flow"));
	/* DESCRIPTION: Output file surface adjoint coefficient (w/o extension) */
	AddScalarOption("SURFACE_ADJ_FILENAME", SurfAdjCoeff_FileName, string("surface_adjoint"));
	/* DESCRIPTION: Output file surface linear coefficient (w/o extension) */
	AddScalarOption("SURFACE_LIN_FILENAME", SurfLinCoeff_FileName, string("surface_linear"));
	/* DESCRIPTION: Writing solution file frequency */
	AddScalarOption("WRT_SOL_FREQ", Wrt_Sol_Freq, 1000);
	/* DESCRIPTION: Writing solution file frequency */
	AddScalarOption("WRT_SOL_FREQ_DUALTIME", Wrt_Sol_Freq_DualTime, 1);
	/* DESCRIPTION: Writing convergence history frequency */
	AddScalarOption("WRT_CON_FREQ",  Wrt_Con_Freq, 1);
	/* DESCRIPTION: Writing convergence history frequency for the dual time */
	AddScalarOption("WRT_CON_FREQ_DUALTIME",  Wrt_Con_Freq_DualTime, 10);
	/* DESCRIPTION: Write a volume solution file */
	AddSpecialOption("WRT_VOL_SOL", Wrt_Vol_Sol, SetBoolOption, true);
	/* DESCRIPTION: Write a surface solution file */
	AddSpecialOption("WRT_SRF_SOL", Wrt_Srf_Sol, SetBoolOption, true);
	/* DESCRIPTION: Write a surface CSV solution file */
	AddSpecialOption("WRT_CSV_SOL", Wrt_Csv_Sol, SetBoolOption, true);
	/* DESCRIPTION: Write a restart solution file */
	AddSpecialOption("WRT_RESTART", Wrt_Restart, SetBoolOption, true);
	/* DESCRIPTION: Write a CGNS solution file */
	AddSpecialOption("WRT_SOL_CGNS", Wrt_Sol_CGNS, SetBoolOption, false);
	/* DESCRIPTION: Write a Tecplot ASCII volume solution file */
	AddSpecialOption("WRT_SOL_TEC_ASCII", Wrt_Sol_Tec_ASCII, SetBoolOption, true);
	/* DESCRIPTION: Write a Tecplot binary volume solution file */
	AddSpecialOption("WRT_SOL_TEC_BINARY", Wrt_Sol_Tec_Binary, SetBoolOption, false);
	/* DESCRIPTION: Output residual info to solution/restart file */
	AddSpecialOption("WRT_RESIDUALS", Wrt_Residuals, SetBoolOption, false);
  /* DESCRIPTION: Output the rind layers in the solution files */
	AddSpecialOption("WRT_HALO", Wrt_Halo, SetBoolOption, false);

	/*--- options related to the equivalent area ---*/
	/* CONFIG_CATEGORY: Equivalent Area */

	/* DESCRIPTION: Evaluate equivalent area on the Near-Field  */
	AddSpecialOption("EQUIV_AREA", EquivArea, SetBoolOption, false);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 1.0; default_vec_3d[2] = 1.0;
	/* DESCRIPTION: Integration limits of the equivalent area ( xmin, xmax, Dist_NearField ) */
	AddArrayOption("EA_INT_LIMIT", 3, EA_IntLimit, default_vec_3d);


	/*--- options related to flow rate ---*/
	/* CONFIG_CATEGORY: Flow rate */

	/* DESCRIPTION:  */
	AddSpecialOption("FLOW_RATE", FlowRate, SetBoolOption, false);

	/*--- options related to the chemical system ---*/
	/* CONFIG_CATEGORY: Chemical system */

	/* DESCRIPTION: Specify chemical model for multi-species simulations */
	AddEnumOption("GAS_MODEL", Kind_GasModel, GasModel_Map, "ARGON");        

	/*--- options related to reference and freestream quantities ---*/
	/* CONFIG_CATEGORY: Reference and freestream quantities */

	Length_Ref = 1.0; //<---- NOTE: this should be given an option or set as a const

	/* DESCRIPTION: Ratio of specific heats (1.4 (air), only for compressible flows) */
	AddScalarOption("GAMMA_VALUE", Gamma, 1.4);
	/* DESCRIPTION: Ratio of specific heats for monatomic gas */
	AddScalarOption("GAMMA_MONATOMIC_VALUE", GammaMonatomic, 5.0/3.0);
	/* DESCRIPTION: Ratio of specific heats for diatomic gas */
	AddScalarOption("GAMMA_DIATOMIC_VALUE", GammaDiatomic, 7.0/5.0);
	/* DESCRIPTION:  */
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	/* DESCRIPTION: Free-stream velocity (m/s) */
	AddArrayOption("FREESTREAM_VELOCITY", 3, Velocity_FreeStream, default_vec_3d);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	/* DESCRIPTION: Reference origin for moment computation */
	AddArrayOption("REF_ORIGIN_MOMENT", 3, RefOriginMoment, default_vec_3d);
	/* DESCRIPTION: Reference area for force coefficients (0 implies automatic calculation) */
	AddScalarOption("REF_AREA", RefAreaCoeff, 1.0);
	/* DESCRIPTION: Reference length for pitching, rolling, and yawing non-dimensional moment */
	AddScalarOption("REF_LENGTH_MOMENT", RefLengthMoment, 1.0);
	/* DESCRIPTION: Reference element length for computing the slope limiter epsilon */
	AddScalarOption("REF_ELEM_LENGTH", RefElemLength, 0.1);
	/* DESCRIPTION: Side-slip angle (degrees, only for compressible flows) */
	AddScalarOption("SIDESLIP_ANGLE", AoS, 0.0);
	/* DESCRIPTION: Angle of attack (degrees, only for compressible flows) */
	AddScalarOption("AOA", AoA, 0.0);
	/* DESCRIPTION: Specific gas constant (287.87 J/kg*K (air), only for compressible flows) */
	AddScalarOption("GAS_CONSTANT", Gas_Constant, 287.87);
	//	AddScalarOption("MIXTURE_MOLAR_MASS", Mixture_Molar_mass, 28.97);
	/* DESCRIPTION: Free-stream pressure (101325.0 N/m^2 by default) */
	AddScalarOption("FREESTREAM_PRESSURE", Pressure_FreeStream, 101325.0);
	/* DESCRIPTION: Free-stream density (1.2886 Kg/m^3 (air), 998.2 Kg/m^3 (water)) */
	AddScalarOption("FREESTREAM_DENSITY", Density_FreeStream, -1.0);
	/* DESCRIPTION: Free-stream temperature (273.15 K by default) */
	AddScalarOption("FREESTREAM_TEMPERATURE", Temperature_FreeStream, 273.15);
	/* DESCRIPTION: Free-stream viscosity (1.853E-5 Ns/m^2 (air), 0.798E-3 Ns/m^2 (water)) */
	AddScalarOption("FREESTREAM_VISCOSITY", Viscosity_FreeStream, -1.0);
	/* DESCRIPTION:  */
	AddScalarOption("FREESTREAM_INTERMITTENCY", Intermittency_FreeStream, 1.0);
	/* DESCRIPTION:  */
	AddScalarOption("FREESTREAM_TURBULENCEINTENSITY", TurbulenceIntensity_FreeStream, 0.05);
	/* DESCRIPTION:  */
	AddScalarOption("FREESTREAM_NU_FACTOR", NuFactor_FreeStream, 3.0);
	/* DESCRIPTION:  */
	AddScalarOption("FREESTREAM_TURB2LAMVISCRATIO", Turb2LamViscRatio_FreeStream, 10.0);
	/* DESCRIPTION: Laminar Prandtl number (0.72 (air), only for compressible flows) */
	AddScalarOption("PRANDTL_LAM", Prandtl_Lam, 0.72);
	/* DESCRIPTION: Turbulent Prandtl number (0.9 (air), only for compressible flows) */
	AddScalarOption("PRANDTL_TURB", Prandtl_Turb, 0.90);
	/* DESCRIPTION: Reference pressure (1.0 N/m^2 by default, only for compressible flows)  */
	AddScalarOption("REF_PRESSURE", Pressure_Ref, 1.0);
	/* DESCRIPTION: Reference temperature (1.0 K by default, only for compressible flows) */
	AddScalarOption("REF_TEMPERATURE", Temperature_Ref, 1.0);
	/* DESCRIPTION: Reference density (1.0 Kg/m^3 by default, only for compressible flows) */
	AddScalarOption("REF_DENSITY", Density_Ref, 1.0);
	/* DESCRIPTION: Reference velocity (incompressible only) */
	AddScalarOption("REF_VELOCITY", Velocity_Ref, -1.0);
	/* DESCRIPTION: Reference viscosity (incompressible only) */
	AddScalarOption("REF_VISCOSITY", Viscosity_Ref, -1.0);
	/* DESCRIPTION:  */
	AddScalarOption("WALL_TEMPERATURE", Wall_Temperature, 300.0);
	/* DESCRIPTION:  Mach number (non-dimensional, based on the free-stream values) */
	AddScalarOption("MACH_NUMBER", Mach, 0.0);
	/* DESCRIPTION: Reynolds number (non-dimensional, based on the free-stream values) */
	AddScalarOption("REYNOLDS_NUMBER", Reynolds, 0.0);
	/* DESCRIPTION: Reynolds length (1 m by default) */
	AddScalarOption("REYNOLDS_LENGTH", Length_Reynolds, 1.0);
	/* DESCRIPTION: Factor for converting the grid to meters */
	AddScalarOption("CONVERT_TO_METER", Conversion_Factor, 1.0);
	/* DESCRIPTION: Write a new mesh converted to meters */
	AddSpecialOption("WRITE_CONVERTED_MESH", Write_Converted_Mesh, SetBoolOption, false);
	/* DESCRIPTION: Value of the Bulk Modulus  */
	AddScalarOption("BULK_MODULUS", Bulk_Modulus, 2.15E9);
	/* DESCRIPTION: Artifical compressibility factor  */
	AddScalarOption("ARTCOMP_FACTOR", ArtComp_Factor, 1.0);
	/* DESCRIPTION:  */
	AddScalarOption("CHARGE_COEFF", ChargeCoeff, -1.0);
	/* DESCRIPTION:  */
	AddListOption("CFL_MS", nSpeciesCFL, CFL_FineGrid_Species);
	/* DESCRIPTION:  */
	AddListOption("CFL_RATE_MS", nSpeciesCFL, CFL_Rate_Species);
	/* DESCRIPTION:  */
	AddListOption("CFL_ITER_MS", nSpeciesCFL, CFL_Iter_Species);
	/* DESCRIPTION:  */
	AddListOption("CFL_MAX_MS", nSpeciesCFL, CFL_Max_Species);
	/* DESCRIPTION:  */
	AddListOption("FREESTREAM_SPECIES_TEMPERATURE", nTemp, Species_Temperature_FreeStream);
	/* DESCRIPTION:  */
	AddListOption("INLET_SPECIES_TEMPERATURE", nTemp, Species_Temperature_Inlet);
	/* DESCRIPTION:  */
	AddListOption("OUTLET_SPECIES_TEMPERATURE", nTemp, Species_Temperature_Outlet);
	/* DESCRIPTION:  */
	AddListOption("INLET_SPECIES_PRESSURE", nTemp, Species_Pressure_Inlet);
	/* DESCRIPTION:  */
	AddListOption("OUTLET_SPECIES_PRESSURE", nTemp, Species_Pressure_Outlet);
	/* DESCRIPTION:  */
	AddListOption("INLET_SPECIES_VELOCITY", nTemp, Species_Velocity_Inlet);
	/* DESCRIPTION:  */
	AddListOption("OUTLET_SPECIES_VELOCITY", nTemp, Species_Velocity_Outlet);
	/* DESCRIPTION:  */
	AddListOption("GAS_COMPOSITION", nTemp, Gas_Composition);
	/* DESCRIPTION:  */
	AddListOption("PARTICLE_REFERENCE_TEMPERATURE", nRef_Temperature, Species_Ref_Temperature);
	/* DESCRIPTION:  */
	AddListOption("PARTICLE_REFERENCE_VISCOSITY", nRef_Viscosity, Species_Ref_Viscosity);

	/*--- options related to free surface simulation ---*/
	/* CONFIG_CATEGORY: Free surface simulation */

	/* DESCRIPTION: Ratio of density for two phase problems */
	AddScalarOption("RATIO_DENSITY", RatioDensity, 0.1);
	/* DESCRIPTION: Ratio of viscosity for two phase problems */
	AddScalarOption("RATIO_VISCOSITY", RatioViscosity, 0.1);
	/* DESCRIPTION: Location of the freesurface (y or z coordinate) */
	AddScalarOption("FREESURFACE_ZERO", FreeSurface_Zero, 0.0);
	/* DESCRIPTION: Free surface depth surface (x or y coordinate) */
	AddScalarOption("FREESURFACE_DEPTH", FreeSurface_Depth, 1.0);
	/* DESCRIPTION: Thickness of the interface in a free surface problem */
	AddScalarOption("FREESURFACE_THICKNESS", FreeSurface_Thickness, 0.1);
	/* DESCRIPTION: Free surface damping coefficient */
	AddScalarOption("FREESURFACE_DAMPING_COEFF", FreeSurface_Damping_Coeff, 0.0);
	/* DESCRIPTION: Free surface damping length (times the baseline wave) */
	AddScalarOption("FREESURFACE_DAMPING_LENGTH", FreeSurface_Damping_Length, 1.0);
	/* DESCRIPTION: Location of the free surface outlet surface (x or y coordinate) */
	AddScalarOption("FREESURFACE_OUTLET", FreeSurface_Outlet, 0.0);
	/* DESCRIPTION: Writing convergence history frequency for the dual time */
	AddScalarOption("FREESURFACE_REEVALUATION",  FreeSurface_Reevaluation, 1);

	/*--- options related to the grid deformation ---*/
	// these options share nDV as their size in the option references; not a good idea
	/* CONFIG_CATEGORY: Grid deformation */

	/* DESCRIPTION: Kind of deformation */
	AddEnumListOption("DV_KIND", nDV, Design_Variable, Param_Map);
	/* DESCRIPTION: Marker of the surface to which we are going apply the shape deformation */
	AddMarkerOption("DV_MARKER", nMarker_Moving, Marker_Moving);
	/* DESCRIPTION: Old value of the deformation for incremental deformations */
	AddListOption("DV_VALUE_OLD", nDV, DV_Value_Old);
	/* DESCRIPTION: New value of the shape deformation */
	AddListOption("DV_VALUE_NEW", nDV, DV_Value_New);
	/* DESCRIPTION: Parameters of the shape deformation
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
    	- FFD_CAMBER ( FFDBox ID, i_Ind, j_Ind )
    	- FFD_THICKNESS ( FFDBox ID, i_Ind, j_Ind )
    	- FFD_VOLUME ( FFDBox ID, i_Ind, j_Ind ) */
	AddDVParamOption("DV_PARAM", nDV, ParamDV, Design_Variable);
	/* DESCRIPTION: Hold the grid fixed in a region */
	AddSpecialOption("HOLD_GRID_FIXED", Hold_GridFixed, SetBoolOption, false);
	default_vec_6d[0] = -1E15; default_vec_6d[1] = -1E15; default_vec_6d[2] = -1E15;
	default_vec_6d[3] =  1E15; default_vec_6d[4] =  1E15; default_vec_6d[5] =  1E15;
	/* DESCRIPTION: Coordinates of the box where the grid will be deformed (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax) */
	AddArrayOption("HOLD_GRID_FIXED_COORD", 6, Hold_GridFixed_Coord, default_vec_6d);
	/* DESCRIPTION: Grid deformation technique */
	AddEnumOption("GRID_DEFORM_METHOD", Kind_GridDef_Method, Deform_Map, "SPRING");
	/* DESCRIPTION: Visualize the deformation */
	AddSpecialOption("VISUALIZE_DEFORMATION", Visualize_Deformation, SetBoolOption, false);
	/* DESCRIPTION: Number of iterations for FEA mesh deformation (surface deformation increments) */
	AddScalarOption("FEA_ITER", FEA_Iter, 1);
  
	/*--- option related to rotorcraft problems ---*/
	/* CONFIG_CATEGORY: Rotorcraft problem */

	AddScalarOption("CYCLIC_PITCH", Cyclic_Pitch, 0.0);
	AddScalarOption("COLLECTIVE_PITCH", Collective_Pitch, 0.0);


	/*--- options related to the FEA solver ---*/
	/* CONFIG_CATEGORY: FEA solver */

	/* DESCRIPTION: Modulus of elasticity */
	AddScalarOption("ELASTICITY_MODULUS", ElasticyMod, 2E11);
	/* DESCRIPTION: Poisson ratio */
	AddScalarOption("POISSON_RATIO", PoissonRatio, 0.30);
	/* DESCRIPTION: Material density */
	AddScalarOption("MATERIAL_DENSITY", MaterialDensity, 7854); 

	/*--- options related to the wave solver ---*/
	/* CONFIG_CATEGORY: Wave solver */

	/* DESCRIPTION: Constant wave speed */
	AddScalarOption("WAVE_SPEED", Wave_Speed, 331.79);

	/*--- options related to the heat solver ---*/
	/* CONFIG_CATEGORY: Heat solver */

	/* DESCRIPTION: Thermal diffusivity constant */
	AddScalarOption("THERMAL_DIFFUSIVITY", Thermal_Diffusivity, 1.172E-5);


	/* END_CONFIG_OPTIONS */

	/*--- Read the configuration file ---*/
	case_file.open(case_filename, ios::in);

	if (case_file.fail()) {
		cout << "There is no configuration file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get(); exit(1);
	}

	/*--- Parse the configuration file and set the options ---*/
	while (getline (case_file,text_line)) {
		if (TokenizeString(text_line, option_name, option_value)) {
			map<string, CAnyOptionRef*>::iterator it;
			it = param.find(option_name);
			if (it != param.end()) {
				param[option_name]->SetValue(option_value);
			} else {
				if ( !GetPython_Option(option_name) && (rank == MASTER_NODE) )
					cout << "WARNING: unrecognized option in the config. file: " << option_name << "." << endl;
			}
		}					
	}	

	case_file.close();

	/*--- Configuration file postprocessing ---*/
	SetPostprocessing(val_software, val_iZone);

	/*--- Configuration file boundaries/markers seting ---*/
	SetMarkers(val_software, val_iZone);

	/*--- Configuration file output ---*/
	if ((rank == MASTER_NODE) && (verb_level == VERB_HIGH) && (val_iZone != 1))
		SetOutput(val_software, val_iZone);

}

CConfig::CConfig(char case_filename[200]) {
	string text_line, text2find, option_name, keyword;
	ifstream case_file;
	vector<string> option_value;	

	/*--- Mesh information ---*/
	AddEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, "SU2");	
	AddScalarOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));

	/*--- Information about type of simulation and time-spectral instances ---*/
	AddEnumOption("UNSTEADY_SIMULATION", Unsteady_Simulation, Unsteady_Map, "NO");
	AddScalarOption("TIME_INSTANCES", nTimeInstances, 1);

	/*--- Read the configuration file ---*/
	case_file.open(case_filename, ios::in);

	if (case_file.fail()) {
		cout << "There is no configuration file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get(); exit(1);
	}

	/*--- Parse the configuration file and set the options ---*/
	while (getline (case_file,text_line)) {
		if (TokenizeString(text_line, option_name, option_value)) {
			map<string, CAnyOptionRef*>::iterator it;
			it = param.find(option_name);
			if (it != param.end()) {
				param[option_name]->SetValue(option_value);
			}
		}
	}	

	case_file.close();

}


void CConfig::SetPostprocessing(unsigned short val_software, unsigned short val_izone) {

	Kind_SU2 = val_software;

  /*--- Divide grid if runnning SU2_MDC ---*/
  if (Kind_SU2 == SU2_MDC) Divide_Element = true;
   
	/*--- Identification of free-surface problem, this problems are always unsteady and incompressible. ---*/
	if (FreeSurface) {
		Incompressible = true;
		if (Kind_Solver == EULER) Kind_Solver = FREE_SURFACE_EULER;
		if (Kind_Solver == NAVIER_STOKES) Kind_Solver = FREE_SURFACE_NAVIER_STOKES;
		if (Unsteady_Simulation != DT_STEPPING_2ND) Unsteady_Simulation = DT_STEPPING_1ST;
	}

	/*--- Decide whether we should be writing unsteady solution files. ---*/
	if (Unsteady_Simulation == STEADY ||
			Unsteady_Simulation == TIME_STEPPING ||
			Unsteady_Simulation == TIME_SPECTRAL ||
			Kind_Solver == FREE_SURFACE_EULER ||
			Kind_Solver == FREE_SURFACE_NAVIER_STOKES)
		Wrt_Unsteady = false;
	else
		Wrt_Unsteady = true;

	/*--- If we're solving a steady problem, make sure that
   there is no grid motion ---*/
	if (Kind_SU2 == SU2_CFD && Unsteady_Simulation == STEADY)
		Grid_Movement = false;

	/*--- If it is not specified, set the mesh motion mach number
   equal to the freestream value. ---*/
	if ((Rotating_Frame || Grid_Movement) && Mach_Motion == 0.0)
		Mach_Motion = Mach;

	/*--- In case rigid-motion parameters have not been declared in the config file, set them equal to zero. ---*/
	if ( Grid_Movement && (Kind_GridMovement[ZONE_0] == RIGID_MOTION) ) {

		unsigned short iZone;

		/*--- Motion Origin: ---*/
		if (Motion_Origin_X == NULL) {
			Motion_Origin_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Motion_Origin_X[iZone] = 0.0;
		}

		if (Motion_Origin_Y == NULL) {
			Motion_Origin_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Motion_Origin_Y[iZone] = 0.0;
		}

		if (Motion_Origin_Z == NULL) {
			Motion_Origin_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Motion_Origin_Z[iZone] = 0.0;
		}

		/*--- Translation: ---*/
		if (Translation_Rate_X == NULL) {
			Translation_Rate_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Translation_Rate_X[iZone] = 0.0;
		}

		if (Translation_Rate_Y == NULL) {
			Translation_Rate_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Translation_Rate_Y[iZone] = 0.0;
		}

		if (Translation_Rate_Z == NULL) {
			Translation_Rate_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Translation_Rate_Z[iZone] = 0.0;
		}

		/*--- Rotation: ---*/
		if (Rotation_Rate_X == NULL) {
			Rotation_Rate_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Rotation_Rate_X[iZone] = 0.0;
		}

		if (Rotation_Rate_Y == NULL) {
			Rotation_Rate_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Rotation_Rate_Y[iZone] = 0.0;
		}

		if (Rotation_Rate_Z == NULL) {
			Rotation_Rate_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Rotation_Rate_Z[iZone] = 0.0;
		}

		/*--- Pitching: ---*/
		if (Pitching_Omega_X == NULL) {
			Pitching_Omega_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Omega_X[iZone] = 0.0;
		}

		if (Pitching_Omega_Y == NULL) {
			Pitching_Omega_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Omega_Y[iZone] = 0.0;
		}

		if (Pitching_Omega_Z == NULL) {
			Pitching_Omega_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Omega_Z[iZone] = 0.0;
		}

		/*--- Pitching Amplitude: ---*/
		if (Pitching_Ampl_X == NULL) {
			Pitching_Ampl_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Ampl_X[iZone] = 0.0;
		}

		if (Pitching_Ampl_Y == NULL) {
			Pitching_Ampl_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Ampl_Y[iZone] = 0.0;
		}

		if (Pitching_Ampl_Z == NULL) {
			Pitching_Ampl_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Ampl_Z[iZone] = 0.0;
		}

		/*--- Pitching Phase: ---*/
		if (Pitching_Phase_X == NULL) {
			Pitching_Phase_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Phase_X[iZone] = 0.0;
		}

		if (Pitching_Phase_Y == NULL) {
			Pitching_Phase_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Phase_Y[iZone] = 0.0;
		}

		if (Pitching_Phase_Z == NULL) {
			Pitching_Phase_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Pitching_Phase_Z[iZone] = 0.0;
		}

		/*--- Plunging: ---*/
		if (Plunging_Omega_X == NULL) {
			Plunging_Omega_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Plunging_Omega_X[iZone] = 0.0;
		}

		if (Plunging_Omega_Y == NULL) {
			Plunging_Omega_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Plunging_Omega_Y[iZone] = 0.0;
		}

		if (Plunging_Omega_Z == NULL) {
			Plunging_Omega_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Plunging_Omega_Z[iZone] = 0.0;
		}

		/*--- Plunging Amplitude: ---*/
		if (Plunging_Ampl_X == NULL) {
			Plunging_Ampl_X = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Plunging_Ampl_X[iZone] = 0.0;
		}

		if (Plunging_Ampl_Y == NULL) {
			Plunging_Ampl_Y = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Plunging_Ampl_Y[iZone] = 0.0;
		}

		if (Plunging_Ampl_Z == NULL) {
			Plunging_Ampl_Z = new double[nZone];
			for (iZone = 0; iZone < nZone; iZone++ )
				Plunging_Ampl_Z[iZone] = 0.0;
		}

	}

	/*--- Use the various rigid-motion input frequencies to determine the period to be used with time-spectral cases.
		  There are THREE types of motion to consider, namely: rotation, pitching, and plunging.
		  The largest period of motion is the one to be used for time-spectral calculations. ---*/
	if (Unsteady_Simulation == TIME_SPECTRAL) {

		unsigned short N_MOTION_TYPES = 3;
		double periods[N_MOTION_TYPES];

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

	}

	/*--- Allocating memory for previous time step solutions of Aeroelastic problem and Intializing variables. ---*/
	if (Grid_Movement && (Kind_GridMovement[ZONE_0] == AEROELASTIC)) {
		Aeroelastic_pitch = 0.0;
		Aeroelastic_plunge = 0.0;
		Aeroelastic_np1 = new double[4];
		Aeroelastic_n   = new double[4];
		Aeroelastic_n1  = new double[4];
		for (int i=0; i<4; i++) {
			Aeroelastic_np1[i] = 0.0;
			Aeroelastic_n[i]   = 0.0;
			Aeroelastic_n1[i]  = 0.0;
		}
	}

	/*--- Set a flag for sliding interfaces so that a search
   and interpolation is performed after each time step. ---*/
	if (nMarker_Sliding > 0)
		Relative_Motion = true;
	else
		Relative_Motion = false;

	if (Kind_Solver == FLUID_STRUCTURE_EULER || Kind_Solver == FLUID_STRUCTURE_NAVIER_STOKES) {
		if (Kind_Turb_Model == SA || Kind_Turb_Model == SST)
			Kind_Solver = FLUID_STRUCTURE_RANS;
		Grid_Movement = true;
	}

	if (FullMG) FinestMesh = nMultiLevel;
	else FinestMesh = MESH_0;

	if ((Kind_Solver == NAVIER_STOKES) &&
			(Kind_Turb_Model == SA || Kind_Turb_Model == SST))
		Kind_Solver = RANS;

	if ((Kind_Solver == FREE_SURFACE_NAVIER_STOKES) &&
			(Kind_Turb_Model == SA || Kind_Turb_Model == SST))
		Kind_Solver = FREE_SURFACE_RANS;

	if (Kind_Solver == FREE_SURFACE_EULER ||
			Kind_Solver == FREE_SURFACE_NAVIER_STOKES ||
			Kind_Solver == FREE_SURFACE_RANS) GravityForce = true;

	Kappa_1st_Flow = Kappa_Flow[0];
	Kappa_2nd_Flow = Kappa_Flow[1];
	Kappa_4th_Flow = Kappa_Flow[2];   
	Kappa_1st_AdjFlow = Kappa_AdjFlow[0];
	Kappa_2nd_AdjFlow = Kappa_AdjFlow[1];
	Kappa_4th_AdjFlow = Kappa_AdjFlow[2];
	Kappa_1st_LinFlow = Kappa_AdjFlow[0];
	Kappa_4th_LinFlow = Kappa_AdjFlow[1];

	Kappa_1st_Plasma = Kappa_Plasma[0];
	Kappa_2nd_Plasma = Kappa_Plasma[1];
	Kappa_4th_Plasma = Kappa_Plasma[2];   
	Kappa_1st_AdjPlasma = Kappa_AdjPlasma[0];
	Kappa_2nd_AdjPlasma = Kappa_AdjPlasma[1];
	Kappa_4th_AdjPlasma = Kappa_AdjPlasma[2];

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

	// override MG Smooth parameters
	if (nMG_PreSmooth != 0)
		MG_PreSmooth[MESH_0] = 1;
	if (nMG_PostSmooth != 0) {
		MG_PostSmooth[MESH_0] = 0;
		MG_PostSmooth[nMultiLevel] = 0;
	}
	if (nMG_CorrecSmooth != 0)
		MG_CorrecSmooth[nMultiLevel] = 0;

	if (Restart) FullMG = false;

	if (Adjoint) {
		if (Kind_Solver == EULER) Kind_Solver = ADJ_EULER;
		if (Kind_Solver == FREE_SURFACE_EULER) Kind_Solver = ADJ_FREE_SURFACE_EULER;
		if (Kind_Solver == AEROACOUSTIC_EULER) Kind_Solver = ADJ_AEROACOUSTIC_EULER;
		if (Kind_Solver == PLASMA_EULER) Kind_Solver = ADJ_PLASMA_EULER;
		if (Kind_Solver == PLASMA_NAVIER_STOKES) Kind_Solver = ADJ_PLASMA_NAVIER_STOKES;
		if (Kind_Solver == NAVIER_STOKES) Kind_Solver = ADJ_NAVIER_STOKES;
		if (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) Kind_Solver = ADJ_FREE_SURFACE_NAVIER_STOKES;
		if (Kind_Solver == FREE_SURFACE_RANS) Kind_Solver = ADJ_FREE_SURFACE_RANS;
		if (Kind_Solver == RANS) Kind_Solver = ADJ_RANS;
	}

	if (Linearized) {
		if (Kind_Solver == EULER) Kind_Solver = LIN_EULER;
	}

	if (Unsteady_Simulation == TIME_STEPPING) {
		nMultiLevel = 0;
		MGCycle = 0;
	}

	if (Unsteady_Simulation == TIME_SPECTRAL) {
		Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
		Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
	}

	if (Rotating_Frame == YES) {
		Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
		Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
	}

	if (Axisymmetric == YES) {
		Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
		Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
	}

	if (GravityForce == YES) {
		Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
		Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
	}
  
	nCFL = nMultiLevel+1; 
	CFL = new double[nCFL];
	CFL[0] = CFLFineGrid;
	if (Adjoint) CFL[0] = CFL[0] * Adj_CFLRedCoeff;
	for (unsigned short iCFL = 1; iCFL < nCFL; iCFL++) 
		CFL[iCFL] = CFL[iCFL-1]*MG_CFLRedCoeff;

	if (nRKStep == 0) {
		RK_Alpha_Step = new double[1]; RK_Alpha_Step[0] = 1.0;
	}


	if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == ADJ_FREE_SURFACE_RANS))
			&& (Kind_ViscNumScheme_Flow == NONE)) {
		cout << "You must define a viscous numerical method for the flow equations!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}

	if (((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == ADJ_FREE_SURFACE_RANS)) && (Kind_ViscNumScheme_AdjFlow == NONE)) {
		cout << "You must define a viscous numerical method for the adjoint Navier-Stokes equations!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}	

	if (((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == ADJ_FREE_SURFACE_RANS)) && (Kind_SourNumScheme_AdjFlow == NONE)) {
		cout << "You must define a source numerical method for the adjoint Navier-Stokes equations!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}

    /*--- Set a flag for viscous simulations ---*/
    Viscous = ((Kind_Solver == NAVIER_STOKES) ||
               (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) ||
               (Kind_Solver == PLASMA_NAVIER_STOKES) ||
               (Kind_Solver == ADJ_NAVIER_STOKES) ||
               (Kind_Solver == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
               (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES) ||
               (Kind_Solver == RANS) ||
               (Kind_Solver == FREE_SURFACE_RANS) ||
               (Kind_Solver == ADJ_RANS) ||
               (Kind_Solver == ADJ_FREE_SURFACE_RANS));
    
	if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_EULER) ||
			(Kind_Solver == PLASMA_NAVIER_STOKES) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
		unsigned short iSpecies, jSpecies, iReaction, ii;
		double sum, conversionFact;
		double GammaMonatomic, GammaDiatomic;
		if (val_izone == ZONE_1 ) {
			Divide_Element = true;
			Restart_FlowFileName = "restart_phi.dat";
			SurfFlowCoeff_FileName = "surface_phi.dat";
		}

		GammaMonatomic = 5.0/3.0;
		GammaDiatomic  = 7.0/5.0;

		switch (Kind_GasModel) {

		case ARGON:
			/*--- Species definitions ---*/
			// Species ordering: Ar, Ar+, e-
			nMonatomics = 3;
			nDiatomics = 0;
			nSpecies = nMonatomics + nDiatomics;
			nReactions = 1;

			/*-- Allocation ---*/
			Molar_Mass         = new double[nSpecies];
			Particle_Mass      = new double[nSpecies];
			Species_Gas_Constant = new double[nSpecies];
			Species_Gamma = new double[nSpecies];
			//	Gas_Composition    = new double[nSpecies];
			Charge_Number      = new int[nSpecies];
			Enthalpy_Formation = new double[nSpecies];
			Reactions = new int**[nReactions];
			Molecular_Diameter = new double [nSpecies];
			for (iReaction = 0; iReaction < nReactions; iReaction++) {
				Reactions[iReaction] = new int*[2];
				for (ii = 0; ii < 2; ii++)
					Reactions[iReaction][ii] = new int[6];
			}
			CharVibTemp          = new double[nSpecies];

			/*--- Molecular properties of constituent species ---*/
			Particle_Mass[0] = 6.63053168E-26;					// [kg/kmol] Ar
			Particle_Mass[1] = 6.6304405861812E-026;		// [kg/kmol] Ar+
			Particle_Mass[2] = 9.10938188E-31;					// [kg/kmol] e-

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Molar_Mass[iSpecies] = Particle_Mass[iSpecies] * AVOGAD_CONSTANT;

			/* uncomment me to run naca_plasma test run for debug
				 cout << "Warning - Running the Debug, NACA 0012 case " << endl;
				 cout << "Warning - Running the Debug, NACA 0012 case " << endl;
				 cout << "Warning - Running the Debug, NACA 0012 case " << endl;
				 Molar_Mass[0] = 4.7960842E-26*AVOGAD_CONSTANT;						// [kg/kmol] air
				 Molar_Mass[1] = 4.7960842E-26*AVOGAD_CONSTANT;			// [kg/kmol] air
				 Molar_Mass[2] = 4.7960842E-26*AVOGAD_CONSTANT;					// [kg/kmol] air
				 GammaMonatomic  = 1.4;
				 cout << "Warning - Running the Debug, NACA 0012 case " << endl;
				 cout << "Warning - Running the Debug, NACA 0012 case " << endl;
				 cout << "Warning - Running the Debug, NACA 0012 case " << endl;
			 */

			Charge_Number[0] = 0;
			Charge_Number[1] = 1;
			Charge_Number[2] = -1;
			Molecular_Diameter[0] = 4E-10;
			Molecular_Diameter[1] = 4E-10;
			Molecular_Diameter[2] = 0.0;
			//  JANAF VALUES [KJ/Kmol]
			Enthalpy_Formation[0] = 0.0;
			Enthalpy_Formation[1] = 0.0;
			Enthalpy_Formation[2] = 0.0;

			/*--- Set initial fraction of number density (Ns / Ntotal) ---*/
			/*		Gas_Composition[0] = 0.98;
				 Gas_Composition[1] = 0.01;
				 Gas_Composition[2] = 0.01;
			 */
			/*--- Set reaction maps ---*/
			Reactions[0][0][0]=0;	Reactions[0][0][1]=nSpecies;	Reactions[0][0][2] =nSpecies;		Reactions[0][1][0] =1;	Reactions[0][1][1]=2;	Reactions[0][1][2] =nSpecies;

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				CharVibTemp[iSpecies] = 0.0;

			break;

		case AIR7:


			/*--- Species definitions ---*/
			// Species ordering: N2, O2, NO, NO+, N, O, e-
			nMonatomics = 7;
			nDiatomics = 0;
			nSpecies = nMonatomics + nDiatomics;
			nReactions = 24;

			/*-- Allocation ---*/
			Species_Gas_Constant = new double[nSpecies];
			Species_Gamma = new double[nSpecies];
			Molar_Mass         = new double[nSpecies];
			Particle_Mass      = new double[nSpecies];
			Molecular_Diameter = new double[nSpecies];
			Gas_Composition    = new double[nSpecies];
			Charge_Number      = new int[nSpecies];
			Enthalpy_Formation = new double[nSpecies];
			Reactions = new int**[nReactions];
			for (iReaction = 0; iReaction < nReactions; iReaction++) {
				Reactions[iReaction] = new int*[2];
				for (ii = 0; ii < 2; ii++)
					Reactions[iReaction][ii] = new int[6];
			}
			ArrheniusCoefficient    = new double[nReactions];
			ArrheniusEta   = new double[nReactions];
			ArrheniusTheta = new double[nReactions];
			CharVibTemp = new double[nSpecies];

			/*--- Molecular properties of constituent species ---*/
			Molar_Mass[0] = 2.0*14.0067;									// [kg/kmol] N2
			Molar_Mass[1] = 2.0*15.9994;									// [kg/kmol] O2
			Molar_Mass[2] = (14.0067+15.9994);						// [kg/kmol] NO
			Molar_Mass[3] = (14.0067+15.9994+5.4858E-4);	// [kg/kmol] NO+
			Molar_Mass[4] = 14.0067;											// [kg/kmol] N
			Molar_Mass[5] = 15.9994;											// [kg/kmol] O
			Molar_Mass[6] = 5.4858E-4;										// [kg/kmol] e-

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Particle_Mass[iSpecies] = Molar_Mass[iSpecies] / AVOGAD_CONSTANT;

			Molecular_Diameter[0] = 1.0E-10;
			Molecular_Diameter[1] = 1.0E-10;
			Molecular_Diameter[2] = 1.0E-10;
			Molecular_Diameter[3] = 1.0E-10;
			Molecular_Diameter[4] = 1.0E-10;
			Molecular_Diameter[5] = 1.0E-10;
			Molecular_Diameter[6] = 2.8179402894E-15;
			Charge_Number[0] = 0;
			Charge_Number[1] = 0;
			Charge_Number[2] = 0;
			Charge_Number[3] = 1;
			Charge_Number[4] = 0;
			Charge_Number[5] = 0;
			Charge_Number[6] = -1;
			//  JANAF VALUES [KJ/Kmol]
			Enthalpy_Formation[0] = 0.0;					//N2
			Enthalpy_Formation[1] = 0.0;					//O2
			Enthalpy_Formation[2] = 90.291E3;			//NO
			Enthalpy_Formation[3] = 990.185E3;		//NO+
			Enthalpy_Formation[4] = 472.683E3;		//N
			Enthalpy_Formation[5] = 249.173E3;		//O
			Enthalpy_Formation[6] = 0.0;					//e-

			/*--- Set initial fraction of number density (Ns / Ntotal) ---*/
			/*				Gas_Composition[0] = 0.78;
				 Gas_Composition[1] = 0.21;
				 Gas_Composition[2] = 0.002;
				 Gas_Composition[3] = 0.002;
				 Gas_Composition[4] = 0.002;
				 Gas_Composition[5] = 0.002;
				 Gas_Composition[6] = 0.002;*/
			Gas_Composition[0] = 0.35;
			Gas_Composition[1] = 0.35;
			Gas_Composition[2] = 0.1;
			Gas_Composition[3] = 0.05;
			Gas_Composition[4] = 0.05;
			Gas_Composition[5] = 0.05;
			Gas_Composition[6] = 0.05;

			/*--- Set reaction maps ---*/
			Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=4;		Reactions[0][1][1]=4;		Reactions[0][1][2] =0;
			Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=4;		Reactions[1][1][1]=4;		Reactions[1][1][2] =1;
			Reactions[2][0][0]=0;		Reactions[2][0][1]=2;		Reactions[2][0][2]=nSpecies;		Reactions[2][1][0]=4;		Reactions[2][1][1]=4;		Reactions[2][1][2] =2;
			Reactions[3][0][0]=0;		Reactions[3][0][1]=3;		Reactions[3][0][2]=nSpecies;		Reactions[3][1][0]=4;		Reactions[3][1][1]=4;		Reactions[3][1][2] =3;
			Reactions[4][0][0]=0;		Reactions[4][0][1]=4;		Reactions[4][0][2]=nSpecies;		Reactions[4][1][0]=4;		Reactions[4][1][1]=4;		Reactions[4][1][2] =4;
			Reactions[5][0][0]=0;		Reactions[5][0][1]=5;		Reactions[5][0][2]=nSpecies;		Reactions[5][1][0]=4;		Reactions[5][1][1]=4;		Reactions[5][1][2] =5;
			Reactions[6][0][0]=0;		Reactions[6][0][1]=6;		Reactions[6][0][2]=nSpecies;		Reactions[6][1][0]=4;		Reactions[6][1][1]=4;		Reactions[6][1][2] =6;
			Reactions[7][0][0] =1;	Reactions[7][0][1] =0;	Reactions[7][0][2] =nSpecies;		Reactions[7][1][0] =5;	Reactions[7][1][1] =5;	Reactions[7][1][2] =0;
			Reactions[8][0][0] =1;	Reactions[8][0][1] =1;	Reactions[8][0][2] =nSpecies;		Reactions[8][1][0] =5;	Reactions[8][1][1] =5;	Reactions[8][1][2] =1;
			Reactions[9][0][0] =1;	Reactions[9][0][1] =2;	Reactions[9][0][2] =nSpecies;		Reactions[9][1][0] =5;	Reactions[9][1][1] =5;	Reactions[9][1][2] =2;
			Reactions[10][0][0]=1;	Reactions[10][0][1]=3;	Reactions[10][0][2]=nSpecies;		Reactions[10][1][0]=5;	Reactions[10][1][1]=5;	Reactions[10][1][2]=3;
			Reactions[11][0][0]=1;	Reactions[11][0][1]=4;	Reactions[11][0][2]=nSpecies;		Reactions[11][1][0]=5;	Reactions[11][1][1]=5;	Reactions[11][1][2]=4;
			Reactions[12][0][0]=1;	Reactions[12][0][1]=5;	Reactions[12][0][2]=nSpecies;		Reactions[12][1][0]=5;	Reactions[12][1][1]=5;	Reactions[12][1][2]=5;
			Reactions[13][0][0]=1;	Reactions[13][0][1]=6;	Reactions[13][0][2]=nSpecies;		Reactions[13][1][0]=5;	Reactions[13][1][1]=5;	Reactions[13][1][2]=6;
			Reactions[14][0][0]=2;	Reactions[14][0][1]=0;	Reactions[14][0][2]=nSpecies;		Reactions[14][1][0]=4;	Reactions[14][1][1]=5;	Reactions[14][1][2]=0;
			Reactions[15][0][0]=2;	Reactions[15][0][1]=1;	Reactions[15][0][2]=nSpecies;		Reactions[15][1][0]=4;	Reactions[15][1][1]=5;	Reactions[15][1][2]=1;
			Reactions[16][0][0]=2;	Reactions[16][0][1]=2;	Reactions[16][0][2]=nSpecies;		Reactions[16][1][0]=4;	Reactions[16][1][1]=5;	Reactions[16][1][2]=2;
			Reactions[17][0][0]=2;	Reactions[17][0][1]=3;	Reactions[17][0][2]=nSpecies;		Reactions[17][1][0]=4;	Reactions[17][1][1]=5;	Reactions[17][1][2]=3;
			Reactions[18][0][0]=2;	Reactions[18][0][1]=4;	Reactions[18][0][2]=nSpecies;		Reactions[18][1][0]=4;	Reactions[18][1][1]=5;	Reactions[18][1][2]=4;
			Reactions[19][0][0]=2;	Reactions[19][0][1]=5;	Reactions[19][0][2]=nSpecies;		Reactions[19][1][0]=4;	Reactions[19][1][1]=5;	Reactions[19][1][2]=5;
			Reactions[20][0][0]=2;	Reactions[20][0][1]=6;	Reactions[20][0][2]=nSpecies;		Reactions[20][1][0]=4;	Reactions[20][1][1]=5;	Reactions[20][1][2]=6;
			Reactions[21][0][0]=0;	Reactions[21][0][1]=5;	Reactions[21][0][2]=nSpecies;		Reactions[21][1][0]=2;	Reactions[21][1][1]=4;	Reactions[21][1][2]=nSpecies;
			Reactions[22][0][0]=2;	Reactions[22][0][1]=5;	Reactions[22][0][2]=nSpecies;		Reactions[22][1][0]=1;	Reactions[22][1][1]=4;	Reactions[22][1][2]=nSpecies;
			Reactions[23][0][0]=4;	Reactions[23][0][1]=5;	Reactions[23][0][2]=nSpecies;		Reactions[23][1][0]=3;	Reactions[23][1][1]=6;	Reactions[23][1][2]=nSpecies;

			/*--- Set Arrhenius coefficients for reactions ---*/
			// Pre-exponential factor
			ArrheniusCoefficient[0]  = 7.0E21;
			ArrheniusCoefficient[1]  = 7.0E21;
			ArrheniusCoefficient[2]  = 7.0E21;
			ArrheniusCoefficient[3]  = 7.0E21;
			ArrheniusCoefficient[4]  = 3.0E22;
			ArrheniusCoefficient[5]  = 3.0E22;
			ArrheniusCoefficient[6]  = 0.0;
			ArrheniusCoefficient[7]  = 2.0E21;
			ArrheniusCoefficient[8]  = 2.0E21;
			ArrheniusCoefficient[9]  = 2.0E21;
			ArrheniusCoefficient[10] = 2.0E21;
			ArrheniusCoefficient[11] = 1.0E22;
			ArrheniusCoefficient[12] = 1.0E22;
			ArrheniusCoefficient[13] = 0.0;
			ArrheniusCoefficient[14] = 5.0E15;
			ArrheniusCoefficient[15] = 5.0E15;
			ArrheniusCoefficient[16] = 5.0E15;
			ArrheniusCoefficient[17] = 5.0E15;
			ArrheniusCoefficient[18] = 1.1E17;
			ArrheniusCoefficient[19] = 1.1E17;
			ArrheniusCoefficient[20] = 0.0;
			ArrheniusCoefficient[21] = 6.4E17;
			//		ArrheniusCoefficient[22] = 8.4E12;
			ArrheniusCoefficient[22] = 0.0;
			ArrheniusCoefficient[23] = 5.3E12;
			// Rate-controlling temperature exponent
			ArrheniusEta[0]  = -1.60;
			ArrheniusEta[1]  = -1.60;
			ArrheniusEta[2]  = -1.60;
			ArrheniusEta[3]  = -1.60;
			ArrheniusEta[4]  = -1.60;
			ArrheniusEta[5]  = -1.60;
			ArrheniusEta[6]  = -1.60;
			ArrheniusEta[7]  = -1.50;
			ArrheniusEta[8]  = -1.50;
			ArrheniusEta[9]  = -1.50;
			ArrheniusEta[10] = -1.50;
			ArrheniusEta[11] = -1.50;
			ArrheniusEta[12] = -1.50;
			ArrheniusEta[13] = -1.50;
			ArrheniusEta[14] = 0.0;
			ArrheniusEta[15] = 0.0;
			ArrheniusEta[16] = 0.0;
			ArrheniusEta[17] = 0.0;
			ArrheniusEta[18] = 0.0;
			ArrheniusEta[19] = 0.0;
			ArrheniusEta[20] = 0.0;
			ArrheniusEta[21] = -1.00;
			ArrheniusEta[22] = 0.0;
			ArrheniusEta[23] = 0.0;
			// Characteristic temperature
			ArrheniusTheta[0] = 113200.0;
			ArrheniusTheta[1] = 113200.0;
			ArrheniusTheta[2] = 113200.0;
			ArrheniusTheta[3] = 113200.0;
			ArrheniusTheta[4] = 113200.0;
			ArrheniusTheta[5] = 113200.0;
			ArrheniusTheta[6] = 113200.0;
			ArrheniusTheta[7] = 59500.0;
			ArrheniusTheta[8] = 59500.0;
			ArrheniusTheta[9] = 59500.0;
			ArrheniusTheta[10] = 59500.0;
			ArrheniusTheta[11] = 59500.0;
			ArrheniusTheta[12] = 59500.0;
			ArrheniusTheta[13] = 59500.0;
			ArrheniusTheta[14] = 75500.0;
			ArrheniusTheta[15] = 75500.0;
			ArrheniusTheta[16] = 75500.0;
			ArrheniusTheta[17] = 75500.0;
			ArrheniusTheta[18] = 75500.0;
			ArrheniusTheta[19] = 75500.0;
			ArrheniusTheta[20] = 75500.0;
			ArrheniusTheta[21] = 38400.0;
			ArrheniusTheta[22] = 19450.0;
			ArrheniusTheta[23] = 31900.0;

			//Characteristic vibrational temperatures for calculating e_vib (K)
			CharVibTemp[0] = 3395.0;
			CharVibTemp[1] = 2239.0;
			CharVibTemp[2] = 2817.0;
			CharVibTemp[3] = 2817.0;

			break;

		case O2:
			/*--- Species definitions ---*/
			// Species ordering: N2, O2, NO, NO+, N, O, e-
			nMonatomics = 1;
			nDiatomics = 1;
			nSpecies = nMonatomics + nDiatomics;
			nReactions = 2;

			/*-- Allocation ---*/
			Species_Gas_Constant = new double[nSpecies];
			Species_Gamma = new double[nSpecies];
			Molar_Mass         = new double[nSpecies];
			Particle_Mass      = new double[nSpecies];
			Molecular_Diameter = new double[nSpecies];
			Gas_Composition    = new double[nSpecies];
			Charge_Number      = new int[nSpecies];
			Enthalpy_Formation = new double[nSpecies];
			Reactions = new int**[nReactions];
			for (iReaction = 0; iReaction < nReactions; iReaction++) {
				Reactions[iReaction] = new int*[2];
				for (ii = 0; ii < 2; ii++)
					Reactions[iReaction][ii] = new int[6];
			}
			ArrheniusCoefficient = new double[nReactions];
			ArrheniusEta         = new double[nReactions];
			ArrheniusTheta       = new double[nReactions];
			CharVibTemp          = new double[nSpecies];

			/*--- Molecular properties of constituent species ---*/
			Molar_Mass[0] = 2.0*15.9994;									// [kg/kmol] O2
			Molar_Mass[1] = 15.9994;											// [kg/kmol] O

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Particle_Mass[iSpecies] = Molar_Mass[iSpecies] / AVOGAD_CONSTANT;

			Molecular_Diameter[0] = 1.0E-10;
			Molecular_Diameter[1] = 1.0E-10;

			Charge_Number[0] = 0;
			Charge_Number[1] = 0;

			//  JANAF VALUES [KJ/Kmol]
			Enthalpy_Formation[0] = 0.0;					//O2
			Enthalpy_Formation[1] = 249.173E3;		//O

			/*--- Set initial fraction of number density (Ns / Ntotal) ---*/
			Gas_Composition[0] = 0.95; //O2
			Gas_Composition[1] = 0.05;  //O

			/*--- Set reaction maps ---*/
			Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=1;		Reactions[0][1][1]=1;		Reactions[0][1][2] =0;
			Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=1;		Reactions[1][1][1]=1;		Reactions[1][1][2] =1;

			/*--- Set Arrhenius coefficients for reactions ---*/
			// Pre-exponential factor
			ArrheniusCoefficient[0]  = 2.0E21;
			ArrheniusCoefficient[1] = 1.0E22;
			// Rate-controlling temperature exponent
			ArrheniusEta[0]  = -1.50;
			ArrheniusEta[1] = -1.50;
			// Characteristic temperature
			ArrheniusTheta[0] = 59500.0;
			ArrheniusTheta[1] = 59500.0;

			//Characteristic vibrational temperatures for calculating e_vib (K)
			CharVibTemp[0] = 2239.0;
			CharVibTemp[1] = 0.0;

			break;

		case N2:
			/*--- Species definitions ---*/
			// Species ordering: N2, N
			nMonatomics = 1;
			nDiatomics = 1;
			nSpecies = nMonatomics + nDiatomics;
			nReactions = 2;

			/*-- Allocation ---*/
			Species_Gas_Constant = new double[nSpecies];
			Species_Gamma = new double[nSpecies];
			Molar_Mass         = new double[nSpecies];
			Particle_Mass      = new double[nSpecies];
			Molecular_Diameter = new double[nSpecies];
			Gas_Composition    = new double[nSpecies];
			Charge_Number      = new int[nSpecies];
			Enthalpy_Formation = new double[nSpecies];
			Reactions = new int**[nReactions];
			for (iReaction = 0; iReaction < nReactions; iReaction++) {
				Reactions[iReaction] = new int*[2];
				for (ii = 0; ii < 2; ii++)
					Reactions[iReaction][ii] = new int[6];
			}
			ArrheniusCoefficient    = new double[nReactions];
			ArrheniusEta   = new double[nReactions];
			ArrheniusTheta = new double[nReactions];
			CharVibTemp = new double[nSpecies];

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

			Blottner = new double*[nSpecies];
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Blottner[iSpecies] = new double[3];

			/*--- Molecular properties of constituent species ---*/
			Molar_Mass[0] = 2.0*14.0067;									// [kg/kmol] N2
			Molar_Mass[1] = 14.0067;											// [kg/kmol] N

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Particle_Mass[iSpecies] = Molar_Mass[iSpecies] / AVOGAD_CONSTANT;

			Molecular_Diameter[0] = 1.0E-10;
			Molecular_Diameter[1] = 1.0E-10;

			Charge_Number[0] = 0;
			Charge_Number[1] = 0;

			/*--- Define formation enthalpy from JANAF tables [KJ/Kmol] ---*/
			Enthalpy_Formation[0] = 0.0;					//N2
			Enthalpy_Formation[1] = 472.683E3;		//N

			/*--- Set initial fraction of number density (Ns / Ntotal) ---*/
			Gas_Composition[0] = 0.99;//0.90; //N2
			Gas_Composition[1] = 0.01;  //N

			/*--- Set reaction maps ---*/
			Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=1;		Reactions[0][1][1]=1;		Reactions[0][1][2] =0;
			Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=1;		Reactions[1][1][1]=1;		Reactions[1][1][2] =1;

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

			//Characteristic vibrational temperatures for calculating e_vib (K)
			CharVibTemp[0] = 3395.0;
			CharVibTemp[1] = 0.0;

			/*--- Collision integral data ---*/
			Omega00[0][0][0] = -6.0614558E-03;  Omega00[0][0][1] = 1.2689102E-01;   Omega00[0][0][2] = -1.0616948E+00;  Omega00[0][0][3] = 8.0955466E+02;
			Omega00[0][1][0] = -1.0796249E-02;  Omega00[0][1][1] = 2.2656509E-01;   Omega00[0][1][2] = -1.7910602E+00;  Omega00[0][1][3] = 4.0455218E+03;
			Omega00[1][0][0] = -1.0796249E-02;  Omega00[1][0][1] = 2.2656509E-01;   Omega00[1][0][2] = -1.7910602E+00;  Omega00[1][0][3] = 4.0455218E+03;
			Omega00[1][1][0] = -9.6083779E-03;  Omega00[1][1][1] = 2.0938971E-01;   Omega00[1][1][2] = -1.7386904E+00;  Omega00[1][1][3] = 3.3587983E+03;

			Omega11[0][0][0] = -7.6303990E-03;  Omega11[0][0][1] = 1.6878089E-01;   Omega11[0][0][2] = -1.4004234E+00;  Omega11[0][0][3] = 2.1427708E+03;
			Omega11[0][1][0] = -8.3493693E-03;  Omega11[0][1][1] = 1.7808911E-01;   Omega11[0][1][2] = -1.4466155E+00;  Omega11[0][1][3] = 1.9324210E+03;
			Omega11[1][0][0] = -8.3493693E-03;  Omega11[1][0][1] = 1.7808911E-01;   Omega11[1][0][2] = -1.4466155E+00;  Omega11[1][0][3] = 1.9324210E+03;
			Omega11[1][1][0] = -7.7439615E-03;  Omega11[1][1][1] = 1.7129007E-01;   Omega11[1][1][2] = -1.4809088E+00;  Omega11[1][1][3] = 2.1284951E+03;

			/*--- Viscosity coefficients for the Blottner et. al. (1971) model ---*/
			Blottner[0][0] = 0.0268142;
			Blottner[0][1] = 0.3177838;
			Blottner[0][2] = -11.3155513;

			Blottner[1][0] = 0.0115572;
			Blottner[1][1] = 0.6031679;
			Blottner[1][2] = -12.4327495;

			break;

		case ARGON_SID:
			/*--- Species definitions ---*/
			// Species ordering: Ar, Ar+, e-
			nMonatomics = 3;
			nDiatomics = 0;
			nSpecies = nMonatomics + nDiatomics;
			nReactions = 3;

			/*-- Allocation ---*/
			Species_Gas_Constant = new double[nSpecies];
			Species_Gamma = new double[nSpecies];
			Molar_Mass         = new double[nSpecies];
			Particle_Mass      = new double[nSpecies];
			Molecular_Diameter = new double[nSpecies];
			Gas_Composition    = new double[nSpecies];
			Charge_Number      = new int[nSpecies];
			Enthalpy_Formation = new double[nSpecies];
			Reactions = new int**[nReactions];
			for (iReaction = 0; iReaction < nReactions; iReaction++) {
				Reactions[iReaction] = new int*[2];
				for (ii = 0; ii < 2; ii++)
					Reactions[iReaction][ii] = new int[6];
			}
			ArrheniusCoefficient    = new double[nReactions];
			ArrheniusEta   = new double[nReactions];
			ArrheniusTheta = new double[nReactions];
			CharVibTemp = new double[nSpecies];

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

			Blottner = new double*[nSpecies];
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Blottner[iSpecies] = new double[3];

			/*--- Molecular properties of constituent species ---*/
			Molar_Mass[0] = 6.63053168E-26*AVOGAD_CONSTANT;									// [kg/kmol] N2
			Molar_Mass[1] = 6.6304405861812E-026*AVOGAD_CONSTANT;		// [kg/kmol] Ar+
			Molar_Mass[2] = 9.10938188E-31*AVOGAD_CONSTANT;													// [kg/kmol] N

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Particle_Mass[iSpecies] = Molar_Mass[iSpecies] / AVOGAD_CONSTANT;

			Molecular_Diameter[0] = 1.0E-10;
			Molecular_Diameter[1] = 1.0E-10;
			Molecular_Diameter[2] = 1.0E-20;

			Charge_Number[0] = 0;
			Charge_Number[1] = 1;
			Charge_Number[2] =-1;

			/*--- Define formation enthalpy from JANAF tables [KJ/Kmol] ---*/
			Enthalpy_Formation[0] = 0.0;					//N2
			Enthalpy_Formation[1] = 472.683E3;		//N
			Enthalpy_Formation[2] = 0.0;

			/*--- Set initial fraction of number density (Ns / Ntotal) ---*/
			Gas_Composition[0] = 0.98;//0.90; //N2
			Gas_Composition[1] = 0.01;  //N
			Gas_Composition[2] = 0.01;

			/*--- Set reaction maps ---*/
			Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=1;		Reactions[0][1][1]=2;		Reactions[0][1][2] =0;
			Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=1;		Reactions[1][1][1]=2;		Reactions[1][1][2] =1;
			Reactions[2][0][0]=0;		Reactions[2][0][1]=2;		Reactions[2][0][2]=nSpecies;		Reactions[2][1][0]=1;		Reactions[2][1][1]=2;		Reactions[2][1][2] =2;

			/*--- Set Arrhenius coefficients for chemical reactions ---*/
			// Pre-exponential factor
			ArrheniusCoefficient[0]  = 7.0E21;
			ArrheniusCoefficient[1]  = 3.0E22;
			ArrheniusCoefficient[2]  = 0.0;

			// Rate-controlling temperature exponent
			ArrheniusEta[0]  = -1.60;
			ArrheniusEta[1]  = -1.60;
			ArrheniusEta[2]  =  0.0;

			// Characteristic temperature
			ArrheniusTheta[0] = 113200.0;
			ArrheniusTheta[1] = 113200.0;
			ArrheniusTheta[2] = 113200.0;

			//Characteristic vibrational temperatures for calculating e_vib (K)
			CharVibTemp[0] = 0.0;
			CharVibTemp[1] = 0.0;
			CharVibTemp[2] = 0.0;

			/*--- Collision integral data ---*/
			Omega00[0][0][0] = -6.0614558E-03;  Omega00[0][0][1] = 1.2689102E-01;   Omega00[0][0][2] = -1.0616948E+00;  Omega00[0][0][3] = 8.0955466E+02;
			Omega00[0][1][0] = -1.0796249E-02;  Omega00[0][1][1] = 2.2656509E-01;   Omega00[0][1][2] = -1.7910602E+00;  Omega00[0][1][3] = 4.0455218E+03;
			Omega00[0][2][0] = -1.0796249E-02;  Omega00[0][2][1] = 2.2656509E-01;   Omega00[0][2][2] = -1.7910602E+00;  Omega00[0][2][3] = 4.0455218E+03;
			Omega00[1][0][0] = -1.0796249E-02;  Omega00[1][0][1] = 2.2656509E-01;   Omega00[1][0][2] = -1.7910602E+00;  Omega00[1][0][3] = 4.0455218E+03;
			Omega00[1][1][0] = -9.6083779E-03;  Omega00[1][1][1] = 2.0938971E-01;   Omega00[1][1][2] = -1.7386904E+00;  Omega00[1][1][3] = 3.3587983E+03;
			Omega00[1][2][0] = -9.6083779E-03;  Omega00[1][2][1] = 2.0938971E-01;   Omega00[1][2][2] = -1.7386904E+00;  Omega00[1][2][3] = 3.3587983E+03;
			Omega00[2][0][0] = -6.0614558E-03;  Omega00[2][0][1] = 1.2689102E-01;   Omega00[2][0][2] = -1.0616948E+00;  Omega00[2][0][3] = 8.0955466E+02;
			Omega00[2][1][0] = -1.0796249E-02;  Omega00[2][1][1] = 2.2656509E-01;   Omega00[2][1][2] = -1.7910602E+00;  Omega00[2][1][3] = 4.0455218E+03;
			Omega00[2][2][0] = -1.0796249E-02;  Omega00[2][2][1] = 2.2656509E-01;   Omega00[2][2][2] = -1.7910602E+00;  Omega00[2][2][3] = 4.0455218E+03;


			Omega11[0][0][0] = -7.6303990E-03;  Omega11[0][0][1] = 1.6878089E-01;   Omega11[0][0][2] = -1.4004234E+00;  Omega11[0][0][3] = 2.1427708E+03;
			Omega11[0][1][0] = -8.3493693E-03;  Omega11[0][1][1] = 1.7808911E-01;   Omega11[0][1][2] = -1.4466155E+00;  Omega11[0][1][3] = 1.9324210E+03;
			Omega11[0][2][0] = -8.3493693E-03;  Omega11[0][2][1] = 1.7808911E-01;   Omega11[0][2][2] = -1.4466155E+00;  Omega11[0][2][3] = 1.9324210E+03;
			Omega11[1][0][0] = -8.3493693E-03;  Omega11[1][0][1] = 1.7808911E-01;   Omega11[1][0][2] = -1.4466155E+00;  Omega11[1][0][3] = 1.9324210E+03;
			Omega11[1][1][0] = -7.7439615E-03;  Omega11[1][1][1] = 1.7129007E-01;   Omega11[1][1][2] = -1.4809088E+00;  Omega11[1][1][3] = 2.1284951E+03;
			Omega11[1][2][0] = -7.7439615E-03;  Omega11[1][2][1] = 1.7129007E-01;   Omega11[1][2][2] = -1.4809088E+00;  Omega11[1][2][3] = 2.1284951E+03;
			Omega11[2][0][0] = -7.6303990E-03;  Omega11[2][0][1] = 1.6878089E-01;   Omega11[2][0][2] = -1.4004234E+00;  Omega11[2][0][3] = 2.1427708E+03;
			Omega11[2][1][0] = -8.3493693E-03;  Omega11[2][1][1] = 1.7808911E-01;   Omega11[2][1][2] = -1.4466155E+00;  Omega11[2][1][3] = 1.9324210E+03;
			Omega11[2][2][0] = -8.3493693E-03;  Omega11[2][2][1] = 1.7808911E-01;   Omega11[2][2][2] = -1.4466155E+00;  Omega11[2][2][3] = 1.9324210E+03;


			/*--- Viscosity coefficients for the Blottner et. al. (1971) model ---*/
			Blottner[0][0] = 0.0268142;
			Blottner[0][1] = 0.3177838;
			Blottner[0][2] = -11.3155513;

			Blottner[1][0] = 0.0115572;
			Blottner[1][1] = 0.6031679;
			Blottner[1][2] = -12.4327495;

			Blottner[2][0] = 1.0;
			Blottner[2][1] = 1.0;
			Blottner[2][2] = 1.0;

			break;

		case AIR5:


			/*--- Species definitions ---*/
			// Species ordering: N2, O2, NO, NO+, N, O, e-
			nMonatomics = 2;
			nDiatomics = 3;
			nSpecies = nMonatomics + nDiatomics;
			nReactions = 17;

			/*-- Allocation ---*/
			Species_Gas_Constant = new double[nSpecies];
			Species_Gamma = new double[nSpecies];
			Molar_Mass         = new double[nSpecies];
			Particle_Mass      = new double[nSpecies];
			Molecular_Diameter = new double[nSpecies];
			Gas_Composition    = new double[nSpecies];
			Charge_Number      = new int[nSpecies];
			Enthalpy_Formation = new double[nSpecies];
			Reactions = new int**[nReactions];
			for (iReaction = 0; iReaction < nReactions; iReaction++) {
				Reactions[iReaction] = new int*[2];
				for (ii = 0; ii < 2; ii++)
					Reactions[iReaction][ii] = new int[6];
			}
			ArrheniusCoefficient    = new double[nReactions];
			ArrheniusEta   = new double[nReactions];
			ArrheniusTheta = new double[nReactions];
			CharVibTemp = new double[nSpecies];

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

			/*--- Molecular properties of constituent species ---*/
			Molar_Mass[0] = 2.0*14.0067;									// [kg/kmol] N2
			Molar_Mass[1] = 2.0*15.9994;									// [kg/kmol] O2
			Molar_Mass[2] = 14.0067 + 15.9994;						// [kg/kmol] NO
			Molar_Mass[3] = 14.0067;											// [kg/kmol] N
			Molar_Mass[4] = 15.9994;											// [kg/kmol] O
			// [kg/kmol] O
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Particle_Mass[iSpecies] = Molar_Mass[iSpecies] / AVOGAD_CONSTANT;

			Molecular_Diameter[0] = 1.0E-10;
			Molecular_Diameter[1] = 1.0E-10;
			Molecular_Diameter[2] = 1.0E-10;
			Molecular_Diameter[3] = 1.0E-10;
			Molecular_Diameter[4] = 1.0E-10;

			Charge_Number[0] = 0;
			Charge_Number[1] = 0;
			Charge_Number[2] = 0;
			Charge_Number[3] = 0;
			Charge_Number[4] = 0;

			//  JANAF VALUES [KJ/Kmol]
			Enthalpy_Formation[0] = 0.0;					//N2
			Enthalpy_Formation[1] = 0.0;					//O2
			Enthalpy_Formation[2] = 90.291E3;			//NO
			Enthalpy_Formation[3] = 472.683E3;		//N
			Enthalpy_Formation[4] = 249.173E3;		//O


			/*--- Set initial fraction of number density (Ns / Ntotal) ---*/
			Gas_Composition[0] = 0.788;	//N2
			Gas_Composition[1] = 0.199;	//O2
			Gas_Composition[2] = 0.011; //NO
			Gas_Composition[3] = 0.001; //N
			Gas_Composition[4] = 0.001;  //O


			/*--- Set reaction maps ---*/
			Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=3;		Reactions[0][1][1]=3;		Reactions[0][1][2] =0;
			Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=3;		Reactions[1][1][1]=3;		Reactions[1][1][2] =1;
			Reactions[2][0][0]=0;		Reactions[2][0][1]=2;		Reactions[2][0][2]=nSpecies;		Reactions[2][1][0]=3;		Reactions[2][1][1]=3;		Reactions[2][1][2] =2;
			Reactions[3][0][0]=0;		Reactions[3][0][1]=3;		Reactions[3][0][2]=nSpecies;		Reactions[3][1][0]=3;		Reactions[3][1][1]=3;		Reactions[3][1][2] =3;
			Reactions[4][0][0]=0;		Reactions[4][0][1]=4;		Reactions[4][0][2]=nSpecies;		Reactions[4][1][0]=3;		Reactions[4][1][1]=3;		Reactions[4][1][2] =4;

			Reactions[5][0][0]=1;		Reactions[5][0][1]=0;		Reactions[5][0][2]=nSpecies;		Reactions[5][1][0]=4;		Reactions[5][1][1]=4;		Reactions[5][1][2] =0;
			Reactions[6][0][0]=1;		Reactions[6][0][1]=1;		Reactions[6][0][2]=nSpecies;		Reactions[6][1][0]=4;		Reactions[6][1][1]=4;		Reactions[6][1][2] =1;
			Reactions[7][0][0]=1;		Reactions[7][0][1]=2;		Reactions[7][0][2]=nSpecies;		Reactions[7][1][0]=4;		Reactions[7][1][1]=4;		Reactions[7][1][2] =2;
			Reactions[8][0][0]=1;		Reactions[8][0][1]=3;		Reactions[8][0][2]=nSpecies;		Reactions[8][1][0]=4;		Reactions[8][1][1]=4;		Reactions[8][1][2] =3;
			Reactions[9][0][0]=1;		Reactions[9][0][1]=4;		Reactions[9][0][2]=nSpecies;		Reactions[9][1][0]=4;		Reactions[9][1][1]=4;		Reactions[9][1][2] =4;

			Reactions[10][0][0]=2;		Reactions[10][0][1]=0;		Reactions[10][0][2]=nSpecies;		Reactions[10][1][0]=3;		Reactions[10][1][1]=4;		Reactions[10][1][2] =0;
			Reactions[11][0][0]=2;		Reactions[11][0][1]=1;		Reactions[11][0][2]=nSpecies;		Reactions[11][1][0]=3;		Reactions[11][1][1]=4;		Reactions[11][1][2] =1;
			Reactions[12][0][0]=2;		Reactions[12][0][1]=2;		Reactions[12][0][2]=nSpecies;		Reactions[12][1][0]=3;		Reactions[12][1][1]=4;		Reactions[12][1][2] =2;
			Reactions[13][0][0]=2;		Reactions[13][0][1]=3;		Reactions[13][0][2]=nSpecies;		Reactions[13][1][0]=3;		Reactions[13][1][1]=4;		Reactions[13][1][2] =3;
			Reactions[14][0][0]=2;		Reactions[14][0][1]=4;		Reactions[14][0][2]=nSpecies;		Reactions[14][1][0]=3;		Reactions[14][1][1]=4;		Reactions[14][1][2] =4;

			Reactions[15][0][0]=0;		Reactions[15][0][1]=3;		Reactions[15][0][2]=nSpecies;		Reactions[15][1][0]=2;		Reactions[15][1][1]=4;		Reactions[15][1][2]= nSpecies;
			Reactions[16][0][0]=2;		Reactions[16][0][1]=3;		Reactions[16][0][2]=nSpecies;		Reactions[16][1][0]=1;		Reactions[16][1][1]=4;		Reactions[16][1][2]= nSpecies;

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
			ArrheniusCoefficient[8] = 1.0E22;
			ArrheniusCoefficient[9] = 1.0E22;
			ArrheniusCoefficient[10] = 5.0E15;
			ArrheniusCoefficient[11] = 5.0E15;
			ArrheniusCoefficient[12] = 5.0E15;
			ArrheniusCoefficient[13] = 1.1E17;
			ArrheniusCoefficient[14] = 1.1E17;
			ArrheniusCoefficient[15] = 5.7E12;
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
			ArrheniusEta[8] = -1.50;
			ArrheniusEta[9] = -1.50;
			ArrheniusEta[10] = 0.0;
			ArrheniusEta[11] = 0.0;
			ArrheniusEta[12] = 0.0;
			ArrheniusEta[13] = 0.0;
			ArrheniusEta[14] = 0.0;
			ArrheniusEta[15] = 0.42;
			ArrheniusEta[16] = 0.00;

			// Characteristic temperature
			ArrheniusTheta[0] = 113200.0;
			ArrheniusTheta[1] = 113200.0;
			ArrheniusTheta[2] = 113200.0;
			ArrheniusTheta[3] = 113200.0;
			ArrheniusTheta[4] = 113200.0;
			ArrheniusTheta[5] = 59500.0;
			ArrheniusTheta[6] = 59500.0;
			ArrheniusTheta[7] = 59500.0;
			ArrheniusTheta[8] = 59500.0;
			ArrheniusTheta[9] = 59500.0;
			ArrheniusTheta[10] = 75500.0;
			ArrheniusTheta[11] = 75500.0;
			ArrheniusTheta[12] = 75500.0;
			ArrheniusTheta[13] = 75500.0;
			ArrheniusTheta[14] = 75500.0;
			ArrheniusTheta[15] = 42938.0;
			ArrheniusTheta[16] = 19400.0;

			//Characteristic vibrational temperatures for calculating e_vib (K)
			CharVibTemp[0] = 3395.0;
			CharVibTemp[1] = 2239.0;
			CharVibTemp[2] = 2817.0;

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


		/*--- Check for consistency ---*/
		sum = 0;
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
			sum += Gas_Composition[iSpecies];
		if ( fabs(sum - 1.0) > 1E-10)
			cout << "WARNING: Sum of ratios of initial number density fractions != 1.0, sum = " << sum << endl;


		Mixture_Molar_mass = 0.0;
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
			Mixture_Molar_mass += Molar_Mass[iSpecies]*Gas_Composition[iSpecies];


		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			/*--- Convert formation enthalpy from KJ/Kmol -> J/kg ---*/
			conversionFact = 1000 / Molar_Mass[iSpecies];
			Enthalpy_Formation[iSpecies] = Enthalpy_Formation[iSpecies] * conversionFact;

			/*--- Array of gas constants [J/kmol] / [kg/kmol] = [J/kg] ---*/
			Species_Gas_Constant[iSpecies] = UNIVERSAL_GAS_CONSTANT / Molar_Mass[iSpecies];

			/*--- Array of ratios of specific heat ---*/
			if (iSpecies < nDiatomics)
				Species_Gamma[iSpecies] = GammaDiatomic;
			else
				Species_Gamma[iSpecies] = GammaMonatomic;
		}    

		if (Species_Temperature_FreeStream == NULL) {
#ifdef NO_MPI
			cout << "WARNING: No species temperature specified, using mean flow freestream value for all species." << endl;
#else
			if (MPI::COMM_WORLD.Get_rank() == MASTER_NODE)
				cout << "WARNING: No species temperature specified, using mean flow freestream value for all species." << endl;
#endif
			Species_Temperature_FreeStream = new double[nSpecies];
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
				Species_Temperature_FreeStream[iSpecies] = Temperature_FreeStream;
		}

		Inlet_Outlet_Defined = false;
		if (Species_Temperature_Inlet != NULL || Species_Pressure_Inlet != NULL || Species_Pressure_Outlet != NULL )
			Inlet_Outlet_Defined = true;

		if (PlasmaMultiTimeSteps) {

			nCFL = nMultiLevel+1;
			CFL_MS = new double*[nSpecies];

			if (CFL_FineGrid_Species == NULL) {
#ifdef NO_MPI
				cout << "WARNING: No species CFL numbers specified for plasma multi-timestepping, using mean flow CFL parameters" << endl;
#else
				if (MPI::COMM_WORLD.Get_rank() == MASTER_NODE)
					cout << "WARNING: No species CFL numbers specified for plasma multi-timestepping, using mean flow CFL parameters" << endl;
#endif
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					CFL_MS[iSpecies] = new double[nCFL];
					CFL_MS[iSpecies][0] = CFL[0];
					for (unsigned short iCFL = 1; iCFL < nCFL; iCFL++)
						CFL_MS[iSpecies][iCFL] = CFL_MS[iSpecies][iCFL-1]*MG_CFLRedCoeff;
				}
			} else {
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					CFL_MS[iSpecies] = new double[nCFL];
					CFL_MS[iSpecies][0] = CFL_FineGrid_Species[iSpecies];
					for (unsigned short iCFL = 1; iCFL < nCFL; iCFL++)
						CFL_MS[iSpecies][iCFL] = CFL_MS[iSpecies][iCFL-1]*MG_CFLRedCoeff;
				}
			}

			if (CFL_Iter_Species == NULL) {
#ifdef NO_MPI
				cout << "WARNING: No species CFL ramp iteration specified for plasma multi-timestepping, using mean flow CFL parameters" << endl;
#else
				if (MPI::COMM_WORLD.Get_rank() == MASTER_NODE)
					cout << "WARNING: No species CFL ramp iteration specified for plasma multi-timestepping, using mean flow CFL parameters" << endl;
#endif
				CFL_Iter_Species = new unsigned short[nSpecies];
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					CFL_Iter_Species[iSpecies] = CFLRamp[1];
				}
			}

			if (CFL_Max_Species == NULL) {
#ifdef NO_MPI
				cout << "WARNING: No maximum CFL specified for plasma multi-timestepping, using mean flow CFL parameters" << endl;
#else
				if (MPI::COMM_WORLD.Get_rank() == MASTER_NODE)
					cout << "WARNING: No maximum CFL specified for plasma multi-timestepping, using mean flow CFL parameters" << endl;
#endif
				CFL_Max_Species = new double[nSpecies];
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					CFL_Max_Species[iSpecies] = CFLRamp[2];
				}
			}
		}

	}


	/*--- Some discrete adjoint requirements ---*/
	if ((GetAdjoint() && (GetKind_Adjoint() == DISCRETE)) && (Kind_Solver==ADJ_EULER||ADJ_NAVIER_STOKES||ADJ_RANS)) {
		SetnExtIter(1);
		SetMGLevels(0);
	}

	delete [] tmp_smooth;


}

void CConfig::SetMarkers(unsigned short val_software, unsigned short val_izone) {

#ifndef NO_MPI
	if ((val_software != SU2_DDC) && (val_software != SU2_MAC) && (val_software != SU2_GDC))
		nDomain = MPI::COMM_WORLD.Get_size();
#endif

	/*--- Boundary (marker) treatment ---*/
	nMarker_All = nMarker_Euler + nMarker_FarField + nMarker_SymWall + nMarker_PerBound + nMarker_NearFieldBound + nMarker_Supersonic_Inlet
			+ nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Inlet + nMarker_Outlet + nMarker_Isothermal + nMarker_HeatFlux
			+ nMarker_NacelleInflow + nMarker_NacelleExhaust + nMarker_Dirichlet_Elec + nMarker_Displacement + nMarker_Load
			+ nMarker_FlowLoad + nMarker_FWH + nMarker_Observer + nMarker_Custom + nMarker_Sliding + 2*nDomain;

	Marker_All_Tag = new string[nMarker_All+2];									// Store the tag that correspond with each marker.
	Marker_All_SendRecv = new short[nMarker_All+2];							// +#domain (send), -#domain (receive) o 0 (no send neither receive).
	Marker_All_Boundary = new unsigned short[nMarker_All+2];		// Store the kind of boundary condition
	Marker_All_Monitoring = new unsigned short[nMarker_All+2];	// Store if the boundary should be monitored.
	Marker_All_Designing = new unsigned short[nMarker_All+2];   // Store if the boundary should be designed.
	Marker_All_Plotting = new unsigned short[nMarker_All+2];		// Store if the boundary should be plotted.
	Marker_All_Moving = new unsigned short[nMarker_All+2];			// Store if the boundary should be moved.
	Marker_All_PerBound = new short[nMarker_All+2];							// Store if the boundary belongs to a periodic boundary.
	Marker_All_Sliding = new unsigned short[nMarker_All+2];			// Store if the boundary belongs to a sliding interface.

	unsigned short iMarker_All, iMarker_Config, iMarker_Euler, iMarker_Custom, iMarker_FarField,
	iMarker_SymWall, iMarker_PerBound, iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Dirichlet,
	iMarker_Inlet, iMarker_Outlet, iMarker_Isothermal, iMarker_HeatFlux, iMarker_NacelleInflow, iMarker_NacelleExhaust, iMarker_Displacement, iMarker_Load,
	iMarker_FlowLoad, iMarker_FWH, iMarker_Observer, iMarker_Neumann, iMarker_Monitoring, iMarker_Designing, iMarker_Plotting, iMarker_Moving,
	iMarker_Supersonic_Inlet, iMarker_Sliding;

	for (iMarker_All = 0; iMarker_All < nMarker_All; iMarker_All++) {
		Marker_All_Tag[iMarker_All] = "NONE";
		Marker_All_SendRecv[iMarker_All] = 0;
		Marker_All_Boundary[iMarker_All] = 0;
		Marker_All_Monitoring[iMarker_All] = 0;
		Marker_All_Designing[iMarker_All] = 0;
		Marker_All_Plotting[iMarker_All] = 0;
		Marker_All_Moving[iMarker_All] = 0;
		Marker_All_PerBound[iMarker_All] = 0;
		Marker_All_Sliding[iMarker_All] = 0;
	}

	nMarker_Config = nMarker_Euler + nMarker_FarField + nMarker_SymWall + nMarker_PerBound + nMarker_NearFieldBound
			+ nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Inlet + nMarker_Outlet + nMarker_Isothermal + nMarker_HeatFlux + nMarker_NacelleInflow + nMarker_NacelleExhaust + nMarker_Supersonic_Inlet + nMarker_Displacement + nMarker_Load + nMarker_FlowLoad + nMarker_FWH + nMarker_Observer + nMarker_Custom + nMarker_Sliding;

	Marker_Config_Tag = new string[nMarker_Config];
	Marker_Config_Boundary = new unsigned short[nMarker_Config];
	Marker_Config_Monitoring = new unsigned short[nMarker_Config];
	Marker_Config_Plotting = new unsigned short[nMarker_Config];
	Marker_Config_Moving = new unsigned short[nMarker_Config];
	Marker_Config_Designing = new unsigned short[nMarker_Config];
	Marker_Config_PerBound = new unsigned short[nMarker_Config];
	Marker_Config_Sliding = new unsigned short[nMarker_Config];
    
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Tag[iMarker_Config] = "NONE";
		Marker_Config_Boundary[iMarker_Config] = 0;
		Marker_Config_Monitoring[iMarker_Config] = 0;
		Marker_Config_Designing[iMarker_Config] = 0;
		Marker_Config_Plotting[iMarker_Config] = 0;
		Marker_Config_Moving[iMarker_Config] = 0;
		Marker_Config_PerBound[iMarker_Config] = 0;
		Marker_Config_Sliding[iMarker_Config] = 0;
	}

	iMarker_Config = 0;
	for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Euler[iMarker_Euler];
		Marker_Config_Boundary[iMarker_Config] = EULER_WALL;
		iMarker_Config++;
	}

	for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
		Marker_Config_Tag[iMarker_Config] = Marker_FarField[iMarker_FarField];
		Marker_Config_Boundary[iMarker_Config] = FAR_FIELD;
		iMarker_Config++;
	}

	for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
		Marker_Config_Tag[iMarker_Config] = Marker_SymWall[iMarker_SymWall];
		Marker_Config_Boundary[iMarker_Config] = SYMMETRY_PLANE;
		iMarker_Config++;
	}

	for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
		Marker_Config_Tag[iMarker_Config] = Marker_PerBound[iMarker_PerBound];
		Marker_Config_Boundary[iMarker_Config] = PERIODIC_BOUNDARY;
		Marker_Config_PerBound[iMarker_Config] = iMarker_PerBound + 1;
		iMarker_Config++;
	}

	for (iMarker_Sliding = 0; iMarker_Sliding < nMarker_Sliding; iMarker_Sliding++) {
		Marker_Config_Tag[iMarker_Config] = Marker_SlideBound[iMarker_Sliding];
		Marker_Config_Boundary[iMarker_Config] = SLIDING_INTERFACE;
		Marker_Config_Sliding[iMarker_Config]  = iMarker_Sliding + 1;
		iMarker_Config++;
	}

	for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
		Marker_Config_Tag[iMarker_Config] = Marker_NearFieldBound[iMarker_NearFieldBound];
		Marker_Config_Boundary[iMarker_Config] = NEARFIELD_BOUNDARY;
		iMarker_Config++;
	}

	for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
		Marker_Config_Tag[iMarker_Config] = Marker_InterfaceBound[iMarker_InterfaceBound];
		Marker_Config_Boundary[iMarker_Config] = INTERFACE_BOUNDARY;
		iMarker_Config++;
	}

	for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Dirichlet[iMarker_Dirichlet];
		Marker_Config_Boundary[iMarker_Config] = DIRICHLET;
		iMarker_Config++;
	}

	for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Inlet[iMarker_Inlet];
		Marker_Config_Boundary[iMarker_Config] = INLET_FLOW;
		iMarker_Config++;
	}

    FanFace_Mach = new double[nMarker_NacelleInflow];
    FanFace_Pressure = new double[nMarker_NacelleInflow];

	for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
		Marker_Config_Tag[iMarker_Config] = Marker_NacelleInflow[iMarker_NacelleInflow];
		Marker_Config_Boundary[iMarker_Config] = NACELLE_INFLOW;
        FanFace_Mach[iMarker_NacelleInflow] = 0.0;
        FanFace_Pressure[iMarker_NacelleInflow] = 0.0;
		iMarker_Config++;
	}

	for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
		Marker_Config_Tag[iMarker_Config] = Marker_NacelleExhaust[iMarker_NacelleExhaust];
		Marker_Config_Boundary[iMarker_Config] = NACELLE_EXHAUST;
		iMarker_Config++;
	}

	for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
		Marker_Config_Boundary[iMarker_Config] = SUPERSONIC_INLET;
		iMarker_Config++;
	}

	for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Neumann[iMarker_Neumann];
		Marker_Config_Boundary[iMarker_Config] = NEUMANN;
		iMarker_Config++;
	}

	for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Custom[iMarker_Custom];
		Marker_Config_Boundary[iMarker_Config] = CUSTOM_BOUNDARY;
		iMarker_Config++;
	}

	for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Outlet[iMarker_Outlet];
		Marker_Config_Boundary[iMarker_Config] = OUTLET_FLOW;
		iMarker_Config++;
	}

	for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Isothermal[iMarker_Isothermal];
		Marker_Config_Boundary[iMarker_Config] = ISOTHERMAL;
		iMarker_Config++;
	}

	for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
		Marker_Config_Tag[iMarker_Config] = Marker_HeatFlux[iMarker_HeatFlux];
		Marker_Config_Boundary[iMarker_Config] = HEAT_FLUX;
		iMarker_Config++;
	}

	for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Displacement[iMarker_Displacement];
		Marker_Config_Boundary[iMarker_Config] = DISPLACEMENT_BOUNDARY;
		iMarker_Config++;
	}

	for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Load[iMarker_Load];
		Marker_Config_Boundary[iMarker_Config] = LOAD_BOUNDARY;
		iMarker_Config++;
	}

	for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
		Marker_Config_Tag[iMarker_Config] = Marker_FlowLoad[iMarker_FlowLoad];
		Marker_Config_Boundary[iMarker_Config] = FLOWLOAD_BOUNDARY;
		iMarker_Config++;
	}

	for (iMarker_FWH = 0; iMarker_FWH < nMarker_FWH; iMarker_FWH++) {
		Marker_Config_Tag[iMarker_Config] = Marker_FWH[iMarker_FWH];
		Marker_Config_Boundary[iMarker_Config] = FWH_SURFACE;
		iMarker_Config++;
	}

	for (iMarker_Observer = 0; iMarker_Observer < nMarker_Observer; iMarker_Observer++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Observer[iMarker_Observer];
		Marker_Config_Boundary[iMarker_Config] = WAVE_OBSERVER;
		iMarker_Config++;
	}

	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Monitoring[iMarker_Config] = NO;
		for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++)
			if (Marker_Config_Tag[iMarker_Config] == Marker_Monitoring[iMarker_Monitoring])
				Marker_Config_Monitoring[iMarker_Config] = YES;
	}
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Designing[iMarker_Config] = NO;
		for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++)
			if (Marker_Config_Tag[iMarker_Config] == Marker_Designing[iMarker_Designing])
				Marker_Config_Designing[iMarker_Config] = YES;
	}

	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Plotting[iMarker_Config] = NO;
		for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++)
			if (Marker_Config_Tag[iMarker_Config] == Marker_Plotting[iMarker_Plotting])
				Marker_Config_Plotting[iMarker_Config] = YES;
	}

	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Moving[iMarker_Config] = NO;
		for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
			if (Marker_Config_Tag[iMarker_Config] == Marker_Moving[iMarker_Moving])
				Marker_Config_Moving[iMarker_Config] = YES;
	}

}

void CConfig::SetOutput(unsigned short val_software, unsigned short val_izone) {
	unsigned short iMarker_Euler, iMarker_Custom, iMarker_FarField,
	iMarker_SymWall, iMarker_PerBound, iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Dirichlet,
	iMarker_Inlet, iMarker_Outlet, iMarker_Isothermal, iMarker_HeatFlux, iMarker_NacelleInflow, iMarker_NacelleExhaust, iMarker_Displacement, iMarker_Load, iMarker_FlowLoad, iMarker_FWH, iMarker_Observer,
	iMarker_Neumann, iMarker_Monitoring, iMarker_Designing, iMarker_Plotting, iMarker_Moving, iMarker_Supersonic_Inlet;

	cout << endl <<"-------------------------------------------------------------------------" << endl;
	cout <<"|     _____   _    _   ___                                              |" << endl;
	cout <<"|    / ____| | |  | | |__ \\       Web: su2.stanford.edu                 |" << endl;
	cout <<"|   | (___   | |  | |    ) |      Twitter: @su2code                     |" << endl;
	cout <<"|    \\___ \\  | |  | |   / /       Forum: www.cfd-online.com/Forums/su2/ |" << endl;
	cout <<"|    ____) | | |__| |  / /_                                             |" << endl;
	switch (val_software) {
	case SU2_CFD: cout << "|   |_____/   \\____/  |____|  Suite (Computational Fluid Dyn. Code)     |" << endl; break;
	case SU2_MDC: cout << "|   |_____/   \\____/  |____|  Suite (Mesh Deformation Code)             |" << endl; break;
	case SU2_GPC: cout << "|   |_____/   \\____/  |____|  Suite (Gradient Projection Code)          |" << endl; break;
	case SU2_DDC: cout << "|   |_____/   \\____/  |____|  Suite (Domain Decomposition Code)         |" << endl; break;
	case SU2_MAC: cout << "|   |_____/   \\____/  |____|  Suite (Mesh Adaptation Code)              |" << endl; break;
	case SU2_GDC: cout << "|   |_____/   \\____/  |____|  Suite (Geometry Design Code)              |" << endl; break;
	case SU2_PBC: cout << "|   |_____/   \\____/  |____|  Suite (Periodic Boundary Code)            |" << endl; break;
	case SU2_SMC: cout << "|   |_____/   \\____/  |____|  Suite (Sliding Mesh Code)                 |" << endl; break;
	case SU2_SOL: cout << "|   |_____/   \\____/  |____|  Suite (Solution Exporting Code)           |" << endl; break;
	}

	cout << "|                             Release 2.0.6                             |" << endl;
	cout <<"-------------------------------------------------------------------------" << endl;

	cout << endl <<"------------------------ Physical case definition -----------------------" << endl;
	if (val_software == SU2_CFD) {
		switch (Kind_Solver) {
		case EULER:
			if (Incompressible) cout << "Incompressible Euler equations." << endl;
			else cout << "Compressible Euler equations." << endl; break;
		case NAVIER_STOKES:
			if (Incompressible) cout << "Incompressible Laminar Navier-Stokes equations." << endl;
			else cout << "Compressible Laminar Navier-Stokes' equations." << endl; break;
		case RANS:
			if (Incompressible) cout << "Incompressible RANS equations." << endl;
			else cout << "Compressible RANS equations." << endl;
			cout << "Turbulence model: ";
			switch (Kind_Turb_Model) {
			case SA:  cout << "Spalart Allmaras" << endl; break;
			case SST: cout << "Menter's SST"     << endl; break;
			}
			break;
			case PLASMA_EULER:
				cout << "Plasma equations (without viscosity)." << endl;
				if (Kind_GasModel == ARGON) cout << "Using 3 species Argon gas model." << endl;
				if (Kind_GasModel == AIR7) cout << "Using 7 species Air gas model." << endl;
				if (Kind_GasModel == AIR5) cout << "Using 5 species Air gas model." << endl;
				if (Kind_GasModel == O2) cout << "Using 2 species Oxygen gas model." << endl;
				if (Kind_GasModel == N2) cout << "Using 2 species Nitrogen gas model." << endl;
				if (Kind_GasModel == ARGON_SID) cout << "Using 2 species Sid gas model." << endl;
				break;
			case ADJ_PLASMA_EULER:
				cout << "Plasma adjoint equations (without viscosity)." << endl;
				if (Kind_GasModel == ARGON) cout << "Using 3 species Argon gas model." << endl;
				if (Kind_GasModel == AIR7) cout << "Using 7 species Air gas model." << endl;
				if (Kind_GasModel == AIR5) cout << "Using 5 species Air gas model." << endl;
				if (Kind_GasModel == O2) cout << "Using 2 species Oxygen gas model." << endl;
				if (Kind_GasModel == N2) cout << "Using 2 species Nitrogen gas model." << endl;
				if (Kind_GasModel == ARGON_SID) cout << "Using 2 species Sid gas model." << endl;
				break;
			case PLASMA_NAVIER_STOKES:
				cout << "Plasma equations (with viscosity)." << endl;
				if (Kind_GasModel == ARGON) cout << "Using 3 species Argon gas model." << endl;
				if (Kind_GasModel == AIR7) cout << "Using 7 species Air gas model." << endl;
				if (Kind_GasModel == AIR5) cout << "Using 5 species Air gas model." << endl;
				if (Kind_GasModel == O2) cout << "Using 2 species Oxygen gas model." << endl;
				if (Kind_GasModel == N2) cout << "Using 2 species Nitrogen gas model." << endl;
				if (Kind_GasModel == ARGON_SID) cout << "Using 2 species Sid gas model." << endl;
				break;
			case ADJ_PLASMA_NAVIER_STOKES:
				cout << "Plasma continuous adjoint equations (with viscosity)." << endl;
				if (Kind_GasModel == ARGON) cout << "Using 3 species Argon gas model." << endl;
				if (Kind_GasModel == AIR7) cout << "Using 7 species Air gas model." << endl;
				if (Kind_GasModel == AIR5) cout << "Using 5 species Air gas model." << endl;
				if (Kind_GasModel == O2) cout << "Using 2 species Oxygen gas model." << endl;
				if (Kind_GasModel == N2) cout << "Using 2 species Nitrogen gas model." << endl;
				if (Kind_GasModel == ARGON_SID) cout << "Using 2 species Sid gas model." << endl;
				break;
			case ELECTRIC_POTENTIAL: cout << "Electric potential equation." << endl; break;
			case WAVE_EQUATION: cout << "Wave equation." << endl; break;
			case LINEAR_ELASTICITY: cout << "Finite Element Analysis." << endl; break;
			case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: cout << "Fluid-structure interaction." << endl; break;
			case ADJ_EULER: cout << "Continuous Euler adjoint equations." << endl; break;
			case ADJ_NAVIER_STOKES:
				if (Frozen_Visc)
					cout << "Continuous Navier-Stokes adjoint equations with frozen (laminar) viscosity." << endl;
				else
					cout << "Continuous Navier-Stokes adjoint equations." << endl;
				break;
			case ADJ_RANS:
				if (Kind_Adjoint == CONTINUOUS) {
					if (Frozen_Visc)
						cout << "Continuous RANS adjoint equations with frozen (laminar and eddy) viscosity." << endl;
					else
						cout << "Continuous RANS adjoint equations." << endl;
				}

				break;
			case LIN_EULER: cout << "Linearized Euler equations." << endl; break;
			case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
				if (Kind_Solver == FREE_SURFACE_EULER)
					cout << "Free surface flow equation. Density ratio: " << RatioDensity << "." << endl;
				if (Kind_Solver == FREE_SURFACE_NAVIER_STOKES)
					cout << "Free surface flow equation. Density ratio: " << RatioDensity <<". Viscosity ratio: "<< RatioViscosity << "." << endl;
				if (Kind_Solver == FREE_SURFACE_RANS) {
					cout << "Free surface flow equation. Density ratio: " << RatioDensity <<". Viscosity ratio: "<< RatioViscosity << "." << endl;
					cout << "Turbulence model: ";
					switch (Kind_Turb_Model) {
					case SA:  cout << "Spalart Allmaras" << endl; break;
					case SST: cout << "Menter's SST"     << endl; break;
					}
				}
				cout << "The free surface is located at: " << FreeSurface_Zero <<", and its thickness is: " << FreeSurface_Thickness << "." << endl;
				break;

		}

		if (!Incompressible) {
			cout << "Mach number: " << Mach <<"."<< endl;
			cout << "Angle of attack (AoA): " << AoA <<" deg, and angle of sideslip (AoS): " << AoS <<" deg."<< endl;
			if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
					(Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS))
				cout << "Reynolds number: " << Reynolds <<"."<< endl;
		}

		if (EquivArea) {
			cout <<"The equivalent area is going to be evaluated on the near-field."<< endl;
			cout <<"The lower integration limit is "<<EA_IntLimit[0]<<", and the upper is "<<EA_IntLimit[1]<<"."<<endl;
			cout <<"The near-field is situated at "<<EA_IntLimit[2]<<"."<<endl;
		}

		if (Grid_Movement) {
			cout << "Performing a dynamic mesh simulation." << endl;
		}

		if (Rotating_Frame) {
			cout << "Performing a simulation in a rotating frame." << endl;
			cout << "Origin of the rotation axis: " << RotAxisOrigin[0] <<","<< RotAxisOrigin[1] <<","<< RotAxisOrigin[2] << ". ";
			cout << "Angular velocity vector: " << Omega[0] <<","<< Omega[1] <<","<< Omega[2] << "." << endl;
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

		if (RefAreaCoeff == 0) cout << "The reference length/area will be computed using y(2D) or z(3D) projection." <<endl;
		else cout << "The reference length/area (force coefficient) is " << RefAreaCoeff << "." <<endl;
		cout << "The reference length (moment computation) is " << RefLengthMoment << "." <<endl;
		cout << "Reference origin (moment computation) is (" <<RefOriginMoment[0]<<", "<<RefOriginMoment[1]<<", "<<RefOriginMoment[2]<<")."<< endl;

    cout << "Surface(s) where the force coefficients are evaluated: ";
		for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
			cout << Marker_Monitoring[iMarker_Monitoring];
			if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
			else cout <<"."<<endl;
		}
    
    cout << "Surface(s) where the objective function is evaluated: ";
		for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++) {
			cout << Marker_Designing[iMarker_Designing];
			if (iMarker_Designing < nMarker_Designing-1) cout << ", ";
			else cout <<"."<<endl;
		}
    
    cout << "Surface(s) plotted in the output file: ";
		for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++) {
			cout << Marker_Plotting[iMarker_Plotting];
			if (iMarker_Plotting < nMarker_Plotting-1) cout << ", ";
			else cout <<"."<<endl;
		}
    
    cout << "Surface(s) afected by a grid movement: ";
		for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
			cout << Marker_Moving[iMarker_Moving];
			if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
			else cout <<"."<<endl;
		}
    
	}

	cout << "Input mesh file name: " << Mesh_FileName << endl;

	if (Divide_Element) cout << "Divide grid elements into triangles and tetrahedra." << endl;

	if (val_software == SU2_GPC) {
		cout << "Input sensitivity file name: " << SurfAdjCoeff_FileName << "." << endl;
	}

	if (val_software == SU2_MAC) {
		switch (Kind_Adaptation) {
		case FULL: case WAKE: case TWOPHASE: case FULL_FLOW: case FULL_ADJOINT: case FULL_LINEAR: case SMOOTHING: case SUPERSONIC_SHOCK:
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

	if (val_software == SU2_MDC) {
		cout << endl <<"---------------------- Grid deformation parameters ----------------------" << endl;
		switch (Kind_GridDef_Method) {
      case SPRING: cout << "Grid deformation using a classical spring method." << endl; break;
      case FEA: cout << "Grid deformation using a linear elasticity method." << endl; break;
		}

		if (Design_Variable[0] != NO_DEFORMATION && Design_Variable[0] != SURFACE_FILE) {
			if (Hold_GridFixed == YES) cout << "Hold some regions of the mesh fixed (hardcode implementation)." <<endl;
			cout << "Geo. design var. definition (markers <-> old def., new def. <-> param):" <<endl;
			for (unsigned short iDV = 0; iDV < nDV; iDV++) {
				switch (Design_Variable[iDV]) {
				case NO_DEFORMATION: cout << "There isn't any deformation." ; break;
				case HICKS_HENNE: cout << "Hicks Henne <-> " ; break;
				case COSINE_BUMP: cout << "Cosine bump <-> " ; break;
				case FOURIER: cout << "Fourier <-> " ; break;
				case SPHERICAL: cout << "Spherical design <-> " ; break;
				case MACH_NUMBER: cout << "Mach number <-> " ; break;
				case DISPLACEMENT: cout << "Displacement design variable."; break;
				case NACA_4DIGITS: cout << "NACA four digits <-> "; break;
				case PARABOLIC: cout << "Parabolic <-> "; break;
				case OBSTACLE: cout << "Obstacle <-> "; break;
				case STRETCH: cout << "Stretch <-> "; break;
				case ROTATION: cout << "Rotation <-> "; break;
				case FFD_CONTROL_POINT: cout << "FFD (control point) <-> "; break;
				case FFD_DIHEDRAL_ANGLE: cout << "FFD (dihedral angle) <-> "; break;
				case FFD_TWIST_ANGLE: cout << "FFD (twist angle) <-> "; break;
				case FFD_ROTATION: cout << "FFD (rotation) <-> "; break;
				case FFD_CAMBER: cout << "FFD (camber) <-> "; break;
				case FFD_THICKNESS: cout << "FFD (thickness) <-> "; break;
				case FFD_VOLUME: cout << "FFD (volume) <-> "; break;
				}
				for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
					cout << Marker_Moving[iMarker_Moving];
					if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
					else cout << " <-> ";
				}
				cout << DV_Value_Old[iDV] <<", "<< DV_Value_New[iDV] << " <-> ";

				if (Design_Variable[iDV] == NO_DEFORMATION) nParamDV = 0;
				if (Design_Variable[iDV] == HICKS_HENNE) nParamDV = 2;
				if (Design_Variable[iDV] == SPHERICAL) nParamDV = 3;
				if (Design_Variable[iDV] == COSINE_BUMP) nParamDV = 3;
				if (Design_Variable[iDV] == FOURIER) nParamDV = 3;
				if (Design_Variable[iDV] == DISPLACEMENT) nParamDV = 3;
				if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
				if (Design_Variable[iDV] == NACA_4DIGITS) nParamDV = 3;
				if (Design_Variable[iDV] == PARABOLIC) nParamDV = 2;
				if (Design_Variable[iDV] == OBSTACLE) nParamDV = 2;
				if (Design_Variable[iDV] == STRETCH) nParamDV = 2;
				if (Design_Variable[iDV] == FFD_CONTROL_POINT) nParamDV = 7;
				if (Design_Variable[iDV] == FFD_DIHEDRAL_ANGLE) nParamDV = 7;
				if (Design_Variable[iDV] == FFD_TWIST_ANGLE) nParamDV = 7;
				if (Design_Variable[iDV] == FFD_ROTATION) nParamDV = 7;
				if (Design_Variable[iDV] == FFD_CAMBER) nParamDV = 3;
				if (Design_Variable[iDV] == FFD_THICKNESS) nParamDV = 3;
				if (Design_Variable[iDV] == FFD_VOLUME) nParamDV = 3;

				for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {
					if (iParamDV == 0) cout << "( ";
					cout << ParamDV[iDV][iParamDV];
					if (iParamDV < nParamDV-1) cout << ", ";
					else cout <<" )"<<endl;
				}
			}
		}
	}

	if (((val_software == SU2_CFD) && ( Linearized )) || (val_software == SU2_GPC)) {
		cout << endl <<"-------------------- Surface deformation parameters ---------------------" << endl;
		cout << "Geo. design var. definition (markers <-> old def., new def. <-> param):" <<endl;
		for (unsigned short iDV = 0; iDV < nDV; iDV++) {
			switch (Design_Variable[iDV]) {
			case NO_DEFORMATION: cout << "There isn't any deformation." ; break;
			case HICKS_HENNE: cout << "Hicks Henne <-> " ; break;
			case COSINE_BUMP: cout << "Cosine bump <-> " ; break;
			case FOURIER: cout << "Fourier <-> " ; break;
			case SPHERICAL: cout << "Spherical design <-> " ; break;
			case MACH_NUMBER: cout << "Mach number <-> " ; break;
			case DISPLACEMENT: cout << "Displacement design variable."; break;
			case NACA_4DIGITS: cout << "NACA four digits <-> "; break;
			case PARABOLIC: cout << "Parabolic <-> "; break;
			case OBSTACLE: cout << "Obstacle <-> "; break;
			case STRETCH: cout << "Stretch <-> "; break;
			case ROTATION: cout << "Rotation <-> "; break;
			case FFD_CONTROL_POINT: cout << "FFD (control point) <-> "; break;
			case FFD_DIHEDRAL_ANGLE: cout << "FFD (dihedral angle) <-> "; break;
			case FFD_TWIST_ANGLE: cout << "FFD (twist angle) <-> "; break;
			case FFD_ROTATION: cout << "FFD (rotation) <-> "; break;
			case FFD_CAMBER: cout << "FFD (camber) <-> "; break;
			case FFD_THICKNESS: cout << "FFD (thickness) <-> "; break;
			case FFD_VOLUME: cout << "FFD (volume) <-> "; break;
			}
			for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
				cout << Marker_Moving[iMarker_Moving];
				if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
				else cout << " <-> ";
			}
			cout << DV_Value_Old[iDV] <<", "<< DV_Value_New[iDV] << " <-> ";

			if (Design_Variable[iDV] == NO_DEFORMATION) nParamDV = 0;
			if (Design_Variable[iDV] == HICKS_HENNE) nParamDV = 2;
			if (Design_Variable[iDV] == COSINE_BUMP) nParamDV = 3;
			if (Design_Variable[iDV] == FOURIER) nParamDV = 3;
			if (Design_Variable[iDV] == SPHERICAL) nParamDV = 3;
			if (Design_Variable[iDV] == DISPLACEMENT) nParamDV = 3;
			if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
			if (Design_Variable[iDV] == NACA_4DIGITS) nParamDV = 3;
			if (Design_Variable[iDV] == PARABOLIC) nParamDV = 2;
			if (Design_Variable[iDV] == OBSTACLE) nParamDV = 2;
			if (Design_Variable[iDV] == STRETCH) nParamDV = 2;
			if (Design_Variable[iDV] == FFD_CONTROL_POINT) nParamDV = 7;
			if (Design_Variable[iDV] == FFD_DIHEDRAL_ANGLE) nParamDV = 7;
			if (Design_Variable[iDV] == FFD_TWIST_ANGLE) nParamDV = 7;
			if (Design_Variable[iDV] == FFD_ROTATION) nParamDV = 7;
			if (Design_Variable[iDV] == FFD_CAMBER) nParamDV = 3;
			if (Design_Variable[iDV] == FFD_THICKNESS) nParamDV = 3;
			if (Design_Variable[iDV] == FFD_VOLUME) nParamDV = 3;

			for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {
				if (iParamDV == 0) cout << "( ";
				cout << ParamDV[iDV][iParamDV];
				if (iParamDV < nParamDV-1) cout << ", ";
				else cout <<" )"<<endl;
			}
		}
	}

	if (((val_software == SU2_CFD) && ( Adjoint || OneShot )) || (val_software == SU2_GPC)) {

		cout << endl <<"----------------------- Design problem definition -----------------------" << endl;
		switch (Kind_ObjFunc) {
		case DRAG_COEFFICIENT: cout << "Drag objective function." << endl; break;
		case LIFT_COEFFICIENT: cout << "Lift objective function." << endl; break;
		case SIDEFORCE_COEFFICIENT: cout << "Side force objective function." << endl; break;
		case MOMENT_X_COEFFICIENT: cout << "Mx objective function." << endl; break;
		case MOMENT_Y_COEFFICIENT: cout << "My objective function." << endl; break;
		case MOMENT_Z_COEFFICIENT: cout << "Mz objective function." << endl; break;
		case EFFICIENCY: cout << "Efficiency objective function." << endl; break;
		case PRESSURE_COEFFICIENT: cout << "Pressure objective function." << endl; break;
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
    case HEAT_LOAD: cout << "Integrated heat flux objective function." << endl; break;
		case FIGURE_OF_MERIT: cout << "Rotor Figure of Merit objective function." << endl; break;
		case FREE_SURFACE: cout << "Free-Surface objective function." << endl; break;
		case NOISE: cout << "Noise objective function." << endl; break;
		}

	}

	if (val_software == SU2_CFD) {
		cout << endl <<"---------------------- Space numerical integration ----------------------" << endl;

		if (SmoothNumGrid) cout << "There are some smoothing iterations on the grid coordinates." <<endl;

		if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)
				|| (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES)
				|| (Kind_Solver == FREE_SURFACE_RANS)) {
			if ((Kind_ConvNumScheme_Flow == SPACE_CENTERED) && (Kind_Centered_Flow == JST)) {
				cout << "Jameson-Schmidt-Turkel scheme for the flow inviscid terms."<< endl;
				cout << "JST viscous coefficients (1st, 2nd & 4th): " << Kappa_1st_Flow
						<< ", " << Kappa_2nd_Flow << ", " << Kappa_4th_Flow <<"."<< endl;
				cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
			}
			if ((Kind_ConvNumScheme_Flow == SPACE_CENTERED) && (Kind_Centered_Flow == LAX))
				cout << "Lax-Friedrich scheme for the flow inviscid terms."<< endl;
			if (Kind_ConvNumScheme_Flow == SPACE_UPWIND) {
				if (Kind_Upwind_Flow == ROE_1ST) cout << "1st order Roe solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == ROE_TURKEL_1ST) cout << "1st order Roe-Turkel solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == AUSM_1ST)	cout << "1st order AUSM solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == HLLC_1ST)	cout << "1st order HLLC solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == SW_1ST)	cout << "1st order Steger-Warming solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == MSW_1ST)	cout << "1st order Modified Steger-Warming solver for the flow inviscid terms."<< endl;
			}
			if ((Kind_ConvNumScheme_Flow == SPACE_UPWIND) &&
					((Kind_Upwind_Flow == ROE_2ND) || (Kind_Upwind_Flow == AUSM_2ND) || (Kind_Upwind_Flow == HLLC_2ND)
							|| (Kind_Upwind_Flow == SW_2ND) || (Kind_Upwind_Flow == MSW_2ND) || (Kind_Upwind_Flow == ROE_TURKEL_2ND))) {
				if (Kind_Upwind_Flow == ROE_2ND) cout << "2nd order Roe solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == ROE_TURKEL_2ND) cout << "2nd order Roe-Turkel solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == AUSM_2ND) cout << "2nd order AUSM solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == HLLC_2ND) cout << "2nd order HLLC solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == SW_2ND) cout << "2nd order Steger-Warming solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Flow == MSW_2ND) cout << "2nd order Modified Steger-Warming solver for the flow inviscid terms."<< endl;
				switch (Kind_SlopeLimit_Flow) {
				case NONE: cout << "Without slope-limiting method." << endl; break;
				case VENKATAKRISHNAN:
					cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
					cout << "The reference element size is: " << RefElemLength <<". "<< endl;
					break;
				case MINMOD:
					cout << "Minmod slope-limiting method." << endl;
					break;
				}
			}
		}

		if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
			if (Kind_ConvNumScheme_Plasma == SPACE_CENTERED) {
				if (Kind_Centered_Plasma == JST)
					cout << "Jameson-Schmidt-Turkel scheme for the plasma inviscid terms."<< endl;
				if (Kind_Centered_Plasma == LAX)
					cout << "Lax-Friedrich scheme for the plasma inviscid terms."<< endl;
			}

			if (Kind_ConvNumScheme_Plasma == SPACE_UPWIND) {
				if (Kind_Upwind_Plasma == ROE_1ST)	cout << "1st order Roe solver for the plasma inviscid terms."<< endl;
				if (Kind_Upwind_Plasma == ROE_TURKEL_1ST) cout << "1st order Roe-Turkel solver for the flow inviscid terms."<< endl;
				if (Kind_Upwind_Plasma == ROE_TURKEL_2ND) cout << "2nd order Roe-Turkel solver for the flow inviscid terms."<< endl;

				if (Kind_Upwind_Plasma == ROE_2ND) {
					cout << "2nd order Roe solver for the plasma inviscid terms."<< endl;
					switch (Kind_SlopeLimit_Plasma) {
					case NONE: cout << "Without slope-limiting method." << endl; break;
					case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method." << endl; break;
					case MINMOD: cout << "Minmod slope-limiting method." << endl; break;
					}
				}
				if (Kind_Upwind_Plasma == SW_1ST) {
					cout << "1st order Steger-Warming solver for the flow inviscid terms."<< endl;
				}
				if (Kind_Upwind_Plasma == SW_2ND) {
					cout << "2nd order Steger-Warming solver for the plasma inviscid terms."<< endl;
					switch (Kind_SlopeLimit_Plasma) {
					case NONE: cout << "Without slope-limiting method." << endl; break;
					case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method." << endl; break;
					case MINMOD: cout << "Minmod slope-limiting method." << endl; break;
					}
				}
				if (Kind_Upwind_Plasma == MSW_1ST) {
					cout << "1st order Modified Steger-Warming solver for the flow inviscid terms."<< endl;
				}
				if (Kind_Upwind_Plasma == MSW_2ND) {
					cout << "2nd order Modified Steger-Warming solver for the plasma inviscid terms."<< endl;
					switch (Kind_SlopeLimit_Plasma) {
					case NONE: cout << "Without slope-limiting method." << endl; break;
					case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method." << endl; break;
					case MINMOD: cout << "Minmod slope-limiting method." << endl; break;
					}
				}
			}
		}

		if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
			if ((Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED) && (Kind_Centered_AdjFlow == JST)) {
				cout << "Jameson-Schmidt-Turkel scheme for the adjoint inviscid terms."<< endl;
				cout << "JST viscous coefficients (1st, 2nd, & 4th): " << Kappa_1st_AdjFlow
						<< ", " << Kappa_2nd_AdjFlow << ", " << Kappa_4th_AdjFlow <<"."<< endl;
				cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
			}
			if ((Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED) && (Kind_Centered_AdjFlow == LAX))
				cout << "Lax-Friedrich scheme for the adjoint inviscid terms."<< endl;
			if ((Kind_ConvNumScheme_AdjFlow == SPACE_UPWIND) && (Kind_Upwind_AdjFlow == ROE_1ST))
				cout << "1st order Roe solver for the adjoint inviscid terms."<< endl;
			if ((Kind_ConvNumScheme_AdjFlow == SPACE_UPWIND) && (Kind_Upwind_AdjFlow == ROE_2ND)) {
				cout << "2nd order Roe solver for the adjoint inviscid terms."<< endl;
				switch (Kind_SlopeLimit_Flow) {
				case NONE: cout << "Without slope-limiting method." << endl; break;
				case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method." << endl; break;
				}
			}
		}

		if (Kind_Solver == LIN_EULER) {
			if ((Kind_ConvNumScheme_LinFlow == SPACE_CENTERED) && (Kind_Centered_LinFlow == JST)) {
				cout << "Jameson-Schmidt-Turkel scheme for the linearized inviscid terms."<< endl;
				cout << "JST viscous coefficients (1st, & 4th): " << Kappa_1st_LinFlow
						<< ", " << Kappa_4th_LinFlow <<"."<< endl;
				cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
			}
			if ((Kind_ConvNumScheme_LinFlow == SPACE_CENTERED) && (Kind_Centered_LinFlow == LAX))
				cout << "Lax-Friedrich scheme for the linearized inviscid terms."<< endl;
		}

		if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
			switch (Kind_ViscNumScheme_Flow) {
			case AVG_GRAD: cout << "Average of gradients (viscous flow terms)." << endl; break;
			case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (viscous flow terms)." << endl; break;
			}
		}
		if (Kind_Solver == PLASMA_NAVIER_STOKES) {
			switch (Kind_ViscNumScheme_Plasma) {
			case AVG_GRAD: cout << "Average of gradients (1st order) for computation of viscous flow terms." << endl; break;
			case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (2nd order) for computation of viscous flow terms." << endl; break;
			}
		}
		if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
			if (Kind_SourNumScheme_Flow == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the flow source terms." << endl;
		}

		if (Kind_Solver == ADJ_EULER) {
			if (Kind_SourNumScheme_AdjFlow == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the adjoint source terms." << endl;
		}

		if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES) ){
			if (Kind_SourNumScheme_Plasma == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the plasma source terms." << endl;
		}

		if (Kind_Solver == RANS) {
			if ((Kind_ConvNumScheme_Turb == SPACE_UPWIND) && (Kind_Upwind_Turb == SCALAR_UPWIND_1ST))
				cout << "Scalar upwind solver (first order) for the turbulence model."<< endl;
			if ((Kind_ConvNumScheme_Turb == SPACE_UPWIND) && (Kind_Upwind_Turb == SCALAR_UPWIND_2ND))
				cout << "Scalar upwind solver (second order) for the turbulence model."<< endl;
		}

		if ((Kind_Solver == ADJ_RANS) && (!Frozen_Visc)) {
			if ((Kind_ConvNumScheme_AdjTurb == SPACE_UPWIND) && (Kind_Upwind_AdjTurb == SCALAR_UPWIND_1ST))
				cout << "Adjoint turbulent eq - Scalar upwind solver (first order)"<< endl;
			if ((Kind_ConvNumScheme_AdjTurb == SPACE_UPWIND) && (Kind_Upwind_AdjTurb == SCALAR_UPWIND_2ND))
				cout << "Adjoint turbulent eq - Scalar upwind solver (second order)"<< endl;
		}

		if ((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
			switch (Kind_ViscNumScheme_AdjFlow) {
			case AVG_GRAD: cout << "Average of gradients (viscous adjoint terms)." << endl; break;
			case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (viscous adjoint terms)." << endl; break;
			}
		}		

		if (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES) {
			switch (Kind_ViscNumScheme_AdjPlasma) {
			case AVG_GRAD: cout << "Average of gradients (1st order) for computation of viscous adjoint terms." << endl; break;
			case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (2nd order) for computation of viscous adjoint terms." << endl; break;
			}
		}

		if (Kind_Solver == RANS) {
			if (Kind_ViscNumScheme_Turb == AVG_GRAD) cout << "Average of gradients (viscous turbulence terms)." << endl;
			if (Kind_ViscNumScheme_Turb == AVG_GRAD_CORRECTED) cout << "Average of gradients with correction (viscous turbulence terms)." << endl;
			if (Kind_SourNumScheme_Turb == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the turbulence model source terms." << endl;
		}

		if ((Kind_Solver == ELECTRIC_POTENTIAL) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
			if (Kind_ViscNumScheme_Elec == GALERKIN) cout << "Galerkin method for viscous terms computation of the electric potential equation." << endl;
		}

		if ((Kind_Solver == ADJ_RANS) && (!Frozen_Visc)) {
			if (Kind_ViscNumScheme_AdjTurb == AVG_GRAD) cout << "Average of gradients (1st order) for computation of adjoint viscous turbulence terms." << endl;
			if (Kind_ViscNumScheme_AdjTurb == AVG_GRAD_CORRECTED) cout << "Average of gradients with correction (2nd order) for computation of adjoint viscous turbulence terms." << endl;
			if (Kind_SourNumScheme_AdjTurb == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the turbulence adjoint model source terms." << endl;
			if (Kind_TimeIntScheme_AdjTurb == EULER_IMPLICIT) cout << "Euler implicit method for the turbulent adjoint equation." << endl;
		}

		if ((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
			if (Kind_SourNumScheme_AdjFlow == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the Navier-Stokes eq. source terms." << endl;
		}

		if (Kind_Solver == ELECTRIC_POTENTIAL) {
			if (Kind_SourNumScheme_Elec == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the electric potential source terms." << endl;
		}

		switch (Kind_Gradient_Method) {
		case GREEN_GAUSS: cout << "Gradient computation using Green-Gauss theorem." << endl; break;
		case WEIGHTED_LEAST_SQUARES: cout << "Gradient Computation using weighted Least-Squares method." << endl; break;
		}

		if (Incompressible) {
			cout << "Artificial compressibility factor: " << ArtComp_Factor << "." <<endl;
		}

		cout << endl <<"---------------------- Time numerical integration -----------------------" << endl;
		switch (Unsteady_Simulation) {
		case NO:
			cout << "Local time stepping (steady state simulation)." << endl; break;
		case TIME_STEPPING:
			cout << "Unsteady simulation using a time stepping strategy."<< endl;
			if (Unst_CFL != 0.0) cout << "Time step computed by the code. Unsteady CFL number: " << Unst_CFL <<"."<<endl;
			else cout << "Unsteady time step provided by the user (s): "<< Delta_UnstTime << "." << endl;
			break;
		case DT_STEPPING_1ST: case DT_STEPPING_2ND:
			if (Unsteady_Simulation == DT_STEPPING_1ST) cout << "Unsteady simulation, dual time stepping strategy (first order in time)."<< endl;
			if (Unsteady_Simulation == DT_STEPPING_2ND) cout << "Unsteady simulation, dual time stepping strategy (second order in time)."<< endl;
			if (Unst_CFL != 0.0) cout << "Time step computed by the code. Unsteady CFL number: " << Unst_CFL <<"."<<endl;
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
					cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<<endl;
					cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<<endl;
					cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<<endl;
					break;
				case FGMRES:
					cout << "FGMRES is used for solving the linear system." << endl;
					cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<<endl;
					cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<<endl;
					cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<<endl;
					break;
				}
				break;
			}
		}

		if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
			switch (Kind_TimeIntScheme_Plasma) {
			case RUNGE_KUTTA_EXPLICIT:
				cout << "Runge-Kutta explicit method for the plasma equations." << endl;
				cout << "Number of steps: " << nRKStep << endl;
				cout << "Alpha coefficients: ";
				for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
					cout << "\t" << RK_Alpha_Step[iRKStep];
				}
				cout << endl;
				break;
			case EULER_EXPLICIT: cout << "Euler explicit method for the plasma equations." << endl; break;
			case EULER_IMPLICIT: cout << "Euler implicit method for the plasma equations." << endl; break;
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

		if (FullMG) {
			cout << "Full Multigrid." << endl;
		}

		if (nMultiLevel !=0) {
			if (nStartUpIter != 0) cout << "A total of " << nStartUpIter << " start up iterations on the fine grid."<< endl;
			if (MGCycle == 0) cout << "V Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;
			if (MGCycle == 1) cout << "W Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;

			cout << "Reduction of the CFL coefficient in the coarse levels: " << MG_CFLRedCoeff <<"."<<endl;
			cout << "Max. number of children in the agglomeration stage: " << MaxChildren <<"."<<endl;
			cout << "Max. length of an agglom. elem. (compared with the domain): " << MaxDimension <<"."<<endl;
			cout << "Damping factor for the residual restriction: " << Damp_Res_Restric <<"."<<endl;
			cout << "Damping factor for the correction prolongation: " << Damp_Correc_Prolong <<"."<<endl;
		}
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
		for (unsigned short iCFL = 0; iCFL < nMultiLevel+1; iCFL++) {
			cout.width(6); cout << CFL[iCFL];
		}
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

		if (Kind_Solver == RANS)
			if (Kind_TimeIntScheme_Turb == EULER_IMPLICIT)
				cout << "Euler implicit time integration for the turbulence model." << endl;
	}

	if (val_software == SU2_CFD)
		cout << endl <<"------------------------- Convergence criteria --------------------------" << endl;

	if (val_software == SU2_CFD) {
		cout << "Maximum number of iterations: " << nExtIter <<"."<<endl;

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
			if (OneShot) cout << "Cauchy criteria for one shot method " << Cauchy_Eps_OneShot<< "."<< endl;
			if (FullMG) cout << "Cauchy criteria for full multigrid " << Cauchy_Eps_FullMG<< "."<< endl;
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

	if (val_software == SU2_MAC) {
		cout << endl <<"----------------------- Grid adaptation strategy ------------------------" << endl;

		switch (Kind_Adaptation) {
		case NONE: break;
		case FULL: cout << "Grid adaptation using a complete refinement." << endl; break;
		case WAKE: cout << "Grid adaptation of the wake." << endl; break;
		case TWOPHASE: cout << "Grid adaptation of the interphase of a free surface flow." << endl; break;
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
			cout << "Power of the dual volume in the adaptation sensor: " << DualVol_Power <<endl;
			cout << "Percentage of new elements in the adaptation process: " << New_Elem_Adapt << "."<<endl;
			break;
		}

		if (Analytical_Surface != NONE)
			cout << "Use analytical definition for including points in the surfaces." <<endl;

	}

	cout << endl <<"-------------------------- Output information ---------------------------" << endl;

	if ((Output_FileFormat==STL) && (val_software != SU2_MDC)){
		cerr << "Error: STL output file format only valid for SU2_MDC" << endl; throw(-1);
	}

	if (val_software == SU2_CFD) {

		cout << "Writing a flow solution every " << Wrt_Sol_Freq <<" iterations."<<endl;
		cout << "Writing the convergence history every " << Wrt_Con_Freq <<" iterations."<<endl;
		if ((Unsteady_Simulation == DT_STEPPING_1ST) || (Unsteady_Simulation == DT_STEPPING_2ND))  {
			cout << "Writing the dual time flow solution every " << Wrt_Sol_Freq_DualTime <<" iterations."<<endl;
			cout << "Writing the dual time convergence history every " << Wrt_Con_Freq_DualTime <<" iterations."<<endl;
		}

		switch (Output_FileFormat) {
		case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
		case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
		case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
		}

		cout << "Convergence history file name: " << Conv_FileName << "." << endl;

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

		if (Adjoint || OneShot) {
			cout << "Adjoint solution file name: " << Solution_AdjFileName << "." << endl;
			cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
			cout << "Adjoint variables file name: " << Adj_FileName << "." << endl;
			cout << "Surface adjoint coefficients file name: " << SurfAdjCoeff_FileName << "." << endl;
		}
	}

	if (val_software == SU2_SOL) {
		switch (Output_FileFormat) {
		case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
		case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
		case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
		}
		cout << "Flow variables file name: " << Flow_FileName << "." << endl;
	}

	if (val_software == SU2_MDC) {
		cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
		if (Visualize_Deformation) cout << "A file will be created to visualize the deformation." << endl;
		else cout << "No file for visualizing the deformation." << endl;
	}

	if (val_software == SU2_PBC) {
		cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
	}

	if (val_software == SU2_SMC) {
		cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
	}

	if (val_software == SU2_GPC) {
		cout << "Output gradient file name: " << ObjFunc_Grad_FileName << ". " << endl;
	}

	if (val_software == SU2_MAC) {
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

	if (val_software == SU2_DDC) {
		if (Visualize_Partition) cout << "Visualize the partitions. " << endl;
		else cout << "Don't visualize the partitions. " << endl;
	}

	cout << endl <<"------------------- Config file boundary information --------------------" << endl;

	if (nMarker_Euler != 0) {
		cout << "Euler wall boundary marker(s): ";
		for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
			cout << Marker_Euler[iMarker_Euler];
			if (iMarker_Euler < nMarker_Euler-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_FarField != 0) {
		cout << "Far-field boundary marker(s): ";
		for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
			cout << Marker_FarField[iMarker_FarField];
			if (iMarker_FarField < nMarker_FarField-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_SymWall != 0) {
		cout << "Symmetry plane boundary marker(s): ";
		for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
			cout << Marker_SymWall[iMarker_SymWall];
			if (iMarker_SymWall < nMarker_SymWall-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_PerBound != 0) {
		cout << "Periodic boundary marker(s): ";
		for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
			cout << Marker_PerBound[iMarker_PerBound];
			if (iMarker_PerBound < nMarker_PerBound-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_NearFieldBound != 0) {
		cout << "Near-field boundary marker(s): ";
		for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
			cout << Marker_NearFieldBound[iMarker_NearFieldBound];
			if (iMarker_NearFieldBound < nMarker_NearFieldBound-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_InterfaceBound != 0) {
		cout << "Interface boundary marker(s): ";
		for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
			cout << Marker_InterfaceBound[iMarker_InterfaceBound];
			if (iMarker_InterfaceBound < nMarker_InterfaceBound-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Dirichlet != 0) {
		cout << "Dirichlet boundary marker(s): ";
		for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
			cout << Marker_Dirichlet[iMarker_Dirichlet];
			if (iMarker_Dirichlet < nMarker_Dirichlet-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_FlowLoad != 0) {
		cout << "Flow Load boundary marker(s): ";
		for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
			cout << Marker_FlowLoad[iMarker_FlowLoad];
			if (iMarker_FlowLoad < nMarker_FlowLoad-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_FWH != 0) {
		cout << "FW-H boundary marker(s): ";
		for (iMarker_FWH = 0; iMarker_FWH < nMarker_FWH; iMarker_FWH++) {
			cout << Marker_FWH[iMarker_FWH];
			if (iMarker_FWH < nMarker_FWH-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Observer != 0) {
		cout << "Wave observer boundary marker(s): ";
		for (iMarker_Observer = 0; iMarker_Observer < nMarker_Observer; iMarker_Observer++) {
			cout << Marker_Observer[iMarker_Observer];
			if (iMarker_Observer < nMarker_Observer-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Neumann != 0) {
		cout << "Neumann boundary marker(s): ";
		for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
			cout << Marker_Neumann[iMarker_Neumann];
			if (iMarker_Neumann < nMarker_Neumann-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Inlet != 0) {
		cout << "Inlet boundary marker(s): ";
		for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
			cout << Marker_Inlet[iMarker_Inlet];
			if (iMarker_Inlet < nMarker_Inlet-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_NacelleInflow != 0) {
		cout << "Nacelle inflow boundary marker(s): ";
		for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
			cout << Marker_NacelleInflow[iMarker_NacelleInflow];
			if (iMarker_NacelleInflow < nMarker_NacelleInflow-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_NacelleExhaust != 0) {
		cout << "Nacelle exhaust boundary marker(s): ";
		for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
			cout << Marker_NacelleExhaust[iMarker_NacelleExhaust];
			if (iMarker_NacelleExhaust < nMarker_NacelleExhaust-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Supersonic_Inlet != 0) {
		cout << "Supersonic inlet boundary marker(s): ";
		for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
			cout << Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
			if (iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Outlet != 0) {
		cout << "Outlet boundary marker(s): ";
		for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
			cout << Marker_Outlet[iMarker_Outlet];
			if (iMarker_Outlet < nMarker_Outlet-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Isothermal != 0) {
		cout << "Isothermal wall boundary marker(s): ";
		for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
			cout << Marker_Isothermal[iMarker_Isothermal];
			if (iMarker_Isothermal < nMarker_Isothermal-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_HeatFlux != 0) {
		cout << "Constant heat flux wall boundary marker(s): ";
		for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
			cout << Marker_HeatFlux[iMarker_HeatFlux];
			if (iMarker_HeatFlux < nMarker_HeatFlux-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Displacement != 0) {
		cout << "Displacement boundary marker(s): ";
		for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
			cout << Marker_Displacement[iMarker_Displacement];
			if (iMarker_Displacement < nMarker_Displacement-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Load != 0) {
		cout << "Load boundary marker(s): ";
		for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
			cout << Marker_Load[iMarker_Load];
			if (iMarker_Load < nMarker_Load-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Neumann != 0) {
		cout << "Neumann boundary marker(s): ";
		for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
			cout << Marker_Neumann[iMarker_Neumann];
			if (iMarker_Neumann < nMarker_Neumann-1) cout << ", ";
			else cout <<"."<<endl;
		}
	}

	if (nMarker_Custom != 0) {
		cout << "Custom boundary marker(s): ";
		for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
			cout << Marker_Custom[iMarker_Custom];
			if (iMarker_Custom < nMarker_Custom-1) cout << ", ";
			else cout <<"."<<endl;
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

void CConfig::AddMarkerOption(const string & name, unsigned short & num_marker, string* & marker) {
	//cout << "Adding Marker option " << name << endl;
	num_marker = 0;
	CAnyOptionRef* option_ref = new CMarkerOptionRef(marker, num_marker);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddConvectOption(const string & name, unsigned short & space, unsigned short & centered,
		unsigned short & upwind) {
	//cout << "Adding Convect option " << name << endl;
	centered = NO_CENTERED;
	upwind = NO_UPWIND;
	space = SPACE_CENTERED;
	CAnyOptionRef* option_ref = new CConvOptionRef(space, centered, upwind);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMathProblem(const string & name, bool & Adjoint, const bool & Adjoint_default,
		bool & OneShot, const bool & OneShot_default,
		bool & Linearized, const bool & Linearized_default,
		bool & Restart_Flow, const bool & Restart_Flow_default) {
	//cout << "Adding Math Problem option " << name << endl;
	Adjoint = Adjoint_default;
	OneShot = OneShot_default;
	Linearized = Linearized_default;
	Restart_Flow = Restart_Flow_default;
	CAnyOptionRef* option_ref = new CMathProblemRef(Adjoint, OneShot, Linearized,
			Restart_Flow);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddDVParamOption(const string & name, unsigned short & nDV, double** & ParamDV,
		unsigned short* & Design_Variable) {
	//cout << "Adding DV Param option " << name << endl;
	CAnyOptionRef* option_ref = new CDVParamOptionRef(nDV, ParamDV, Design_Variable);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerPeriodic(const string & name, unsigned short & nMarker_PerBound,
		string* & Marker_PerBound, string* & Marker_PerDonor,
		double** & RotCenter, double** & RotAngles, double** & Translation) {
	//cout << "Adding Marker Periodic option " << name << endl;
	nMarker_PerBound = 0;
	CAnyOptionRef* option_ref = new CMarkerPeriodicRef(nMarker_PerBound, Marker_PerBound,
			Marker_PerDonor, RotCenter,
			RotAngles, Translation);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerSliding(const string & name, unsigned short & nMarker_Sliding,
		string* & Marker_SlideBound, string* & Marker_SlideDonor,
		unsigned short* & SlideBound_Zone, unsigned short* & SlideDonor_Zone) {
	nMarker_Sliding = 0;
	CAnyOptionRef* option_ref = new CMarkerSlidingRef(nMarker_Sliding, Marker_SlideBound,
			Marker_SlideDonor, SlideBound_Zone, SlideDonor_Zone);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerInlet(const string & name, unsigned short & nMarker_Inlet,
		string* & Marker_Inlet, double* & Ttotal, double* & Ptotal,
		double** & FlowDir) {
	nMarker_Inlet = 0;
	CAnyOptionRef* option_ref = new CMarkerInletRef(nMarker_Inlet, Marker_Inlet,
			Ttotal, Ptotal, FlowDir);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerInlet(const string & name, unsigned short & nMarker_Inlet,
		string* & Marker_Inlet, double* & Ttotal, double* & Ptotal) {
	nMarker_Inlet = 0;
	CAnyOptionRef* option_ref = new CMarkerInletRef_(nMarker_Inlet, Marker_Inlet,
			Ttotal, Ptotal);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerDirichlet(const string & name, unsigned short & nMarker_Dirichlet_Elec,
		string* & Marker_Dirichlet_Elec, double* & Dirichlet_Value) {
	nMarker_Dirichlet_Elec = 0;
	CAnyOptionRef* option_ref = new CMarkerDirichletRef(nMarker_Dirichlet_Elec, Marker_Dirichlet_Elec,
			Dirichlet_Value);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );

}

void CConfig::AddMarkerOutlet(const string & name, unsigned short & nMarker_Outlet, 
		string* & Marker_Outlet, double* & Pressure) {
	nMarker_Outlet = 0;
	CAnyOptionRef* option_ref = new CMarkerOutletRef(nMarker_Outlet, Marker_Outlet,
			Pressure);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerDisplacement(const string & name, unsigned short & nMarker_Displacement, 
		string* & Marker_Displacement, double* & Displ) {
	nMarker_Displacement = 0;
	CAnyOptionRef* option_ref = new CMarkerDisplacementRef(nMarker_Displacement, Marker_Displacement, Displ);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerLoad(const string & name, unsigned short & nMarker_Load, 
		string* & Marker_Load, double* & Force) {
	nMarker_Load = 0;
	CAnyOptionRef* option_ref = new CMarkerLoadRef(nMarker_Load, Marker_Load, Force);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerFlowLoad(const string & name, unsigned short & nMarker_FlowLoad, 
		string* & Marker_FlowLoad, double* & FlowForce) {
	nMarker_FlowLoad = 0;
	CAnyOptionRef* option_ref = new CMarkerLoadRef(nMarker_FlowLoad, Marker_FlowLoad, FlowForce);
	param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::SetBoolOption(bool* ref, const vector<string> & value) {
	if ( (value[0] != "YES") && (value[0] != "NO") ) {
		cerr << "Error in CConfig::SetBoolOption(): "
				<< "option value provided must be \"YES\" or \"NO\";"
				<< "value given is " << value[0] << endl;
		throw(-1);
	}
	if (value[0] == "YES") {
		*ref = true;
	} else {
		*ref = false;
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

bool CConfig::GetPython_Option(string & option_name) {

	bool isPython_Option = false;

	/*--- Check option name against all known Python options
	 for a match. These are the design options that are
	 never read by the SU2 C++ codes, and we would like
	 to ignore them while processing the config file. ---*/
	if (option_name == "OBJFUNC")       isPython_Option = true;
	if (option_name == "OBJFUNC_SCALE")    isPython_Option = true;
	if (option_name == "CONST_IEQ")      isPython_Option = true;
	if (option_name == "CONST_IEQ_SCALE")      isPython_Option = true;
	if (option_name == "CONST_IEQ_SIGN")     isPython_Option = true;
	if (option_name == "CONST_IEQ_VALUE") isPython_Option = true;
	if (option_name == "CONST_EQ")      isPython_Option = true;
	if (option_name == "CONST_EQ_SCALE")      isPython_Option = true;
	if (option_name == "CONST_EQ_SIGN")     isPython_Option = true;
	if (option_name == "CONST_EQ_VALUE") isPython_Option = true;
	if (option_name == "DEFINITION_DV") isPython_Option = true;
	if (option_name == "TASKS") isPython_Option = true;
	if (option_name == "OPT_OBJECTIVE") isPython_Option = true;
	if (option_name == "OPT_CONSTRAINT") isPython_Option = true;
	if (option_name == "GRADIENTS") isPython_Option = true;
	if (option_name == "FIN_DIFF_STEP") isPython_Option = true;
	if (option_name == "ADAPT_CYCLES") isPython_Option = true;
	if (option_name == "CONSOLE") isPython_Option = true;
	if (option_name == "DECOMPOSED") isPython_Option = true;

	return isPython_Option;
}

unsigned short CConfig::GetMarker_Config_Tag(string val_marker) {

	unsigned short iMarker_Config;

	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker)
			return iMarker_Config;

	cout <<"The configuration file doesn't have any definition for marker "<< val_marker <<"!!" << endl;
	cout <<"Press any key to exit..." << endl;
	cin.get();
	exit(1);
}

unsigned short CConfig::GetMarker_Config_Boundary(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_Boundary[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Monitoring(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_Monitoring[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Designing(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_Designing[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Plotting(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_Plotting[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Moving(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_Moving[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_PerBound(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_PerBound[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Sliding(string val_marker) {
	unsigned short iMarker_Config;
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
		if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
	return Marker_Config_Sliding[iMarker_Config];
}

CConfig::~CConfig(void)
{
	delete [] RK_Alpha_Step;
	delete [] MG_PreSmooth;
	delete [] MG_PostSmooth;
	delete [] U_FreeStreamND;

	/*--- If allocated, delete arrays for Plasma solver ---*/
	if (Species_Gas_Constant != NULL)
		delete [] Species_Gas_Constant;
	if (Species_Gamma != NULL)
		delete [] Species_Gamma;
	if (Molar_Mass != NULL)
		delete[] Molar_Mass;
	if (Particle_Mass != NULL)
		delete [] Particle_Mass;
	if (Molecular_Diameter != NULL)
		delete [] Molecular_Diameter;
	if (Gas_Composition != NULL)
		delete [] Gas_Composition;
	if (Charge_Number != NULL)
		delete [] Charge_Number;
	if (Enthalpy_Formation != NULL)
		delete [] Enthalpy_Formation;
	if (ArrheniusCoefficient != NULL)
		delete [] ArrheniusCoefficient;
	if (ArrheniusEta != NULL)
		delete [] ArrheniusEta;
	if (ArrheniusTheta != NULL)
		delete [] ArrheniusTheta;
	if (CharVibTemp != NULL)
		delete [] CharVibTemp;
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
	if (Grid_Movement && (Kind_GridMovement[ZONE_0] == AEROELASTIC)) {
		delete[] Aeroelastic_np1;
		delete[] Aeroelastic_n;
		delete[] Aeroelastic_n1;
	}

	/*--- Free memory for unspecified Rigid Motion parameters ---*/

	/*--- motion origin: ---*/
	if (Motion_Origin_X != NULL)
		delete [] Motion_Origin_X;
	if (Motion_Origin_Y != NULL)
		delete [] Motion_Origin_Y;
	if (Motion_Origin_Z != NULL)
		delete [] Motion_Origin_Z;

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

	if (Unsteady_Farfield && Unsteady_Simulation) {
		delete [] Density_FreeStreamND_Time;
		delete [] Pressure_FreeStreamND_Time;
		delete [] Mach_Inf_Time;
		delete [] Energy_FreeStreamND_Time;
		delete [] Density_FreeStreamND_Time;
		for (unsigned long iter = 0; iter < nExtIter; iter ++)
			delete []Velocity_FreeStreamND_Time[iter];
		delete []Velocity_FreeStreamND_Time;

	}
}

void CConfig::SetFileNameDomain(unsigned short val_domain) {

#ifndef NO_MPI

	string old_name;
	char buffer[10]; 

	/*--- Standard surface output ---*/
	old_name = SurfFlowCoeff_FileName;
	if (MPI::COMM_WORLD.Get_size() > 1) {
		sprintf (buffer, "_%d", int(val_domain)); 
		SurfFlowCoeff_FileName = old_name + buffer;	
	}

	old_name = SurfAdjCoeff_FileName;
	if (MPI::COMM_WORLD.Get_size() > 1) {
		sprintf (buffer, "_%d", int(val_domain)); 
		SurfAdjCoeff_FileName = old_name + buffer;
	}

	if (MPI::COMM_WORLD.Get_size() > 1) {

		/*--- Standard flow and adjoint output ---*/
		sprintf (buffer, "_%d", int(val_domain));
		old_name = Flow_FileName;
		Flow_FileName = old_name + buffer;

		sprintf (buffer, "_%d", int(val_domain));
		old_name = Structure_FileName;
		Structure_FileName = old_name + buffer;

		sprintf (buffer, "_%d", int(val_domain));
		old_name = Adj_FileName;
		Adj_FileName = old_name + buffer;

		/*--- Mesh files ---*/	
		sprintf (buffer, "_%d.su2", int(val_domain));
		old_name = Mesh_FileName;
		old_name.erase (old_name.end()-4, old_name.end());
		Mesh_FileName = old_name + buffer;

	}
#endif
}

unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {

	switch (val_eqsystem) {
	case RUNTIME_POT_SYS: return FLOW_SOL;
	case RUNTIME_PLASMA_SYS: return PLASMA_SOL;
	case RUNTIME_FLOW_SYS: return FLOW_SOL;
	case RUNTIME_TURB_SYS: return TURB_SOL;
	case RUNTIME_TRANS_SYS: return TRANS_SOL;
	case RUNTIME_ELEC_SYS: return ELEC_SOL;
	case RUNTIME_WAVE_SYS: return WAVE_SOL;
	case RUNTIME_FEA_SYS: return FEA_SOL;
	case RUNTIME_ADJPOT_SYS: return ADJFLOW_SOL;
	case RUNTIME_ADJFLOW_SYS: return ADJFLOW_SOL;
	case RUNTIME_ADJTURB_SYS: return ADJTURB_SOL;
	case RUNTIME_ADJPLASMA_SYS: return ADJPLASMA_SOL;
	case RUNTIME_ADJLEVELSET_SYS: return ADJLEVELSET_SOL;
	case RUNTIME_LINPOT_SYS: return LINFLOW_SOL;
	case RUNTIME_LINFLOW_SYS: return LINFLOW_SOL;
	case RUNTIME_LEVELSET_SYS: return LEVELSET_SOL;
	case RUNTIME_MULTIGRID_SYS: return 0;
	}
	return 0;
}

void CConfig::SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme, 
		unsigned short val_kind_centered, unsigned short val_kind_upwind,
		unsigned short val_kind_slopelimit) {

	Kind_ConvNumScheme = val_kind_convnumscheme; 
	Kind_Centered = val_kind_centered;
	Kind_Upwind = val_kind_upwind;
	Kind_SlopeLimit = val_kind_slopelimit;
}

void CConfig::UpdateCFL(unsigned long val_iter) {
	double coeff;
	bool change;
	unsigned short iCFL;

	if ((val_iter % int(CFLRamp[1]) == 0 ) && (val_iter != 0)) {
		change = false;
		for (iCFL = 0; iCFL <= nMultiLevel; iCFL++) {
			coeff = pow(MG_CFLRedCoeff, double(iCFL));
			if (Adjoint) coeff = coeff * Adj_CFLRedCoeff;

			if (CFL[iCFL]*CFLRamp[0] < CFLRamp[2]*coeff) {
				CFL[iCFL] = CFL[iCFL]*CFLRamp[0];
				change = true;
			}
		}

#ifdef NO_MPI
		if (change) {
			cout <<"\n New value of the CFL number: ";
			for (iCFL = 0; iCFL < nMultiLevel; iCFL++)
				cout << CFL[iCFL] <<", ";
			cout << CFL[nMultiLevel] <<".\n"<< endl;
		}
#else
		int rank = MPI::COMM_WORLD.Get_rank();
		if ((change) && (rank == MASTER_NODE)) {
			cout <<"\n New value of the CFL number: ";
			for (iCFL = 0; iCFL < nMultiLevel; iCFL++)
				cout << CFL[iCFL] <<", ";
			cout << CFL[nMultiLevel] <<".\n"<< endl;
		}
#endif
	}

	if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_EULER) ||
			(Kind_Solver == PLASMA_NAVIER_STOKES) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
		if (CFL_FineGrid_Species != NULL && PlasmaMultiTimeSteps) {

			for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ((val_iter % CFL_Iter_Species[iSpecies] == 0 ) && (val_iter != 0)) {
					change = false;
					for (iCFL = 0; iCFL <= nMultiLevel; iCFL++) {
						coeff = pow(MG_CFLRedCoeff, double(iCFL));
						if (CFL_MS[iSpecies][iCFL]*CFL_Rate_Species[iSpecies] < CFL_Max_Species[iSpecies]*coeff) {
							CFL_MS[iSpecies][iCFL] = CFL_MS[iSpecies][iCFL]*CFL_Rate_Species[iSpecies];
							change = true;
						}
					}

#ifdef NO_MPI
					if (change) {
						cout <<"\n New value of the CFL number for Species: " << iSpecies << " = ";
						for (iCFL = 0; iCFL < nMultiLevel; iCFL++)
							cout << CFL_MS[iSpecies][iCFL] <<", ";
						cout << CFL_MS[iSpecies][nMultiLevel] <<".\n"<< endl;
					}
#else
					int rank = MPI::COMM_WORLD.Get_rank();
					if ((change) && (rank == MASTER_NODE)) {
						cout <<"\n New value of the CFL number for Species: " << iSpecies << " = ";
						for (iCFL = 0; iCFL < nMultiLevel; iCFL++)
							cout << CFL_MS[iSpecies][iCFL] <<", ";
						cout << CFL_MS[iSpecies][nMultiLevel] <<".\n"<< endl;
					}
#endif
				}
			}
		}
	}

}

void CConfig::SetGlobalParam(unsigned short val_solver, unsigned short val_system, unsigned long val_extiter) {

	/*--- Set the simulation global time ---*/
	Current_UnstTime = static_cast<double>(val_extiter)*Delta_UnstTime;
	Current_UnstTimeND = static_cast<double>(val_extiter)*Delta_UnstTimeND;

	/*--- Set the solver methods ---*/
	switch (val_solver) {
	case EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		break;
	case NAVIER_STOKES:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		break;
	case RANS:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_TURB_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centered_Turb(),
					GetKind_Upwind_Turb(), GetKind_SlopeLimit_Turb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
		}
		if (val_system == RUNTIME_TRANS_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centered_Turb(),
					GetKind_Upwind_Turb(), GetKind_SlopeLimit_Turb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
		}
		break;
	case PLASMA_EULER:
		if (val_system == RUNTIME_PLASMA_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Plasma(), GetKind_Centered_Plasma(),
					GetKind_Upwind_Plasma(), GetKind_SlopeLimit_Plasma());
			SetKind_ViscNumScheme(NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Plasma());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Plasma());
		}
		if (val_system == RUNTIME_ELEC_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
			SetKind_TimeIntScheme(NONE);
		}
		break;
	case PLASMA_NAVIER_STOKES:
		if (val_system == RUNTIME_PLASMA_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Plasma(), GetKind_Centered_Plasma(),
					GetKind_Upwind_Plasma(), GetKind_SlopeLimit_Plasma());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Plasma());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Plasma());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Plasma());
		}
		if (val_system == RUNTIME_ELEC_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
			SetKind_TimeIntScheme(NONE);
		}
		break;
	case ADJ_PLASMA_EULER:
		if (val_system == RUNTIME_PLASMA_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Plasma(), GetKind_Centered_Plasma(),
					GetKind_Upwind_Plasma(), GetKind_SlopeLimit_Plasma());
			SetKind_ViscNumScheme(NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Plasma());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Plasma());
		}
		if (val_system == RUNTIME_ADJPLASMA_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjPlasma(), GetKind_Centered_AdjPlasma(),
					GetKind_Upwind_AdjPlasma(), GetKind_SlopeLimit_AdjPlasma());
			SetKind_ViscNumScheme(NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjPlasma());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjPlasma());
		}
		break;
	case ADJ_PLASMA_NAVIER_STOKES:
		if (val_system == RUNTIME_PLASMA_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Plasma(), GetKind_Centered_Plasma(),
					GetKind_Upwind_Plasma(), GetKind_SlopeLimit_Plasma());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Plasma());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Plasma());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Plasma());
		}
		if (val_system == RUNTIME_ADJPLASMA_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjPlasma(), GetKind_Centered_AdjPlasma(),
					GetKind_Upwind_AdjPlasma(), GetKind_SlopeLimit_AdjPlasma());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjPlasma());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjPlasma());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjPlasma());
		}
		break;
	case FREE_SURFACE_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_LEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centered_LevelSet(),
					GetKind_Upwind_LevelSet(), GetKind_SlopeLimit_LevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_LevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
		}
		break;
	case FREE_SURFACE_NAVIER_STOKES:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_LEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centered_LevelSet(),
					GetKind_Upwind_LevelSet(), GetKind_SlopeLimit_LevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_LevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
		}
		break;
	case FREE_SURFACE_RANS:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_TURB_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centered_Turb(),
					GetKind_Upwind_Turb(), GetKind_SlopeLimit_Turb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
		}
		if (val_system == RUNTIME_LEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centered_LevelSet(),
					GetKind_Upwind_LevelSet(), GetKind_SlopeLimit_LevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_LevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
		}
		break;
	case ADJ_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		break;
	case ADJ_NAVIER_STOKES:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(NONE);
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		break;
	case ADJ_RANS:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_SourNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		if (val_system == RUNTIME_TURB_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centered_Turb(),
					GetKind_Upwind_Turb(), GetKind_SlopeLimit_Turb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
		}
		if (val_system == RUNTIME_ADJTURB_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjTurb(), GetKind_Centered_AdjTurb(),
					GetKind_Upwind_AdjTurb(), GetKind_SlopeLimit_AdjTurb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjTurb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjTurb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjTurb());
		}
		break;
	case ADJ_FREE_SURFACE_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_LEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centered_LevelSet(),
					GetKind_Upwind_LevelSet(), GetKind_SlopeLimit_LevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_LevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		if (val_system == RUNTIME_ADJLEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjLevelSet(), GetKind_Centered_AdjLevelSet(),
					GetKind_Upwind_AdjLevelSet(), GetKind_SlopeLimit_AdjLevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjLevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjLevelSet());
		}
		break;
	case ADJ_FREE_SURFACE_NAVIER_STOKES:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_LEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centered_LevelSet(),
					GetKind_Upwind_LevelSet(), GetKind_SlopeLimit_LevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_LevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		if (val_system == RUNTIME_ADJLEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjLevelSet(), GetKind_Centered_AdjLevelSet(),
					GetKind_Upwind_AdjLevelSet(), GetKind_SlopeLimit_AdjLevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjLevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjLevelSet());
		}
		break;
	case ADJ_FREE_SURFACE_RANS:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_TURB_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centered_Turb(),
					GetKind_Upwind_Turb(), GetKind_SlopeLimit_Turb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
		}
		if (val_system == RUNTIME_LEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centered_LevelSet(),
					GetKind_Upwind_LevelSet(), GetKind_SlopeLimit_LevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_LevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		if (val_system == RUNTIME_ADJTURB_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjTurb(), GetKind_Centered_AdjTurb(),
					GetKind_Upwind_AdjTurb(), GetKind_SlopeLimit_AdjTurb());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjTurb());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjTurb());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjTurb());
		}
		if (val_system == RUNTIME_ADJLEVELSET_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjLevelSet(), GetKind_Centered_AdjLevelSet(),
					GetKind_Upwind_AdjLevelSet(), GetKind_SlopeLimit_AdjLevelSet());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjLevelSet());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjLevelSet());
		}
		break;
	case LIN_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(NONE); SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_LINFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_LinFlow(), GetKind_Centered_LinFlow(),
					GetKind_Upwind_LinFlow(), NONE);
			SetKind_ViscNumScheme(NONE); SetKind_SourNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_LinFlow());
		}
		break;
	case ELECTRIC_POTENTIAL:
		if (val_system == RUNTIME_ELEC_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
			SetKind_TimeIntScheme(NONE);
		}
		break;
	case WAVE_EQUATION:
		if (val_system == RUNTIME_WAVE_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Wave());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Wave());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Wave());
		}
		break;
	case LINEAR_ELASTICITY:
		if (val_system == RUNTIME_FEA_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_FEA());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_FEA());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_FEA());
		}
		break;
	case FLUID_STRUCTURE_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_FEA_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_FEA());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_FEA());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_FEA());
		}
		break;
	case AEROACOUSTIC_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_WAVE_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Wave());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Wave());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Wave());
		}
		break;
	case ADJ_AEROACOUSTIC_EULER:
		if (val_system == RUNTIME_FLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
					GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
		}
		if (val_system == RUNTIME_WAVE_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Wave());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Wave());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Wave());
		}
		if (val_system == RUNTIME_ADJFLOW_SYS) {
			SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centered_AdjFlow(),
					GetKind_Upwind_AdjFlow(), GetKind_SlopeLimit_AdjFlow());
			SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
			SetKind_ViscNumScheme(NONE);
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
		}
		if (val_system == RUNTIME_WAVE_SYS) {
			SetKind_ConvNumScheme(NONE, NONE, NONE, NONE);
			SetKind_SourNumScheme(GetKind_SourNumScheme_Wave());
			SetKind_ViscNumScheme(GetKind_ViscNumScheme_Wave());
			SetKind_TimeIntScheme(GetKind_TimeIntScheme_Wave());
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
		if (Marker_PerBound[jMarker_PerBound] == Marker_All_Tag[kMarker_All]) break;

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

string CConfig::GetMarker_Sliding_Donor(string val_marker) {
	unsigned short iMarker_SlideBound;

	/*--- Find the marker for this sliding boundary. ---*/
	for (iMarker_SlideBound = 0; iMarker_SlideBound < nMarker_Sliding; iMarker_SlideBound++)
		if (Marker_SlideBound[iMarker_SlideBound] == val_marker) break;

	/*--- Return the tag for the sliding donor boundary. ---*/
	return Marker_SlideDonor[iMarker_SlideBound];
}

unsigned short CConfig::GetSlideDonor_Zone(string val_marker) {
	unsigned short iMarker_SlideBound;

	/*--- Find the marker for this sliding boundary. ---*/
	for (iMarker_SlideBound = 0; iMarker_SlideBound < nMarker_Sliding; iMarker_SlideBound++)
		if (Marker_SlideBound[iMarker_SlideBound] == val_marker) break;

	return SlideDonor_Zone[iMarker_SlideBound];
}

unsigned short CConfig::GetSlideBound_Zone(string val_marker) {
	unsigned short iMarker_SlideBound;

	/*--- Find the marker for this sliding boundary. ---*/
	for (iMarker_SlideBound = 0; iMarker_SlideBound < nMarker_Sliding; iMarker_SlideBound++)
		if (Marker_SlideBound[iMarker_SlideBound] == val_marker) break;

	return SlideBound_Zone[iMarker_SlideBound];
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
	unsigned short iMarker_NacelleExhaust;
	for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++)
		if (Marker_NacelleExhaust[iMarker_NacelleExhaust] == val_marker) break;
	return Nozzle_Ttotal[iMarker_NacelleExhaust];
}

double CConfig::GetNozzle_Ptotal(string val_marker) {
	unsigned short iMarker_NacelleExhaust;
	for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++)
		if (Marker_NacelleExhaust[iMarker_NacelleExhaust] == val_marker) break;
	return Nozzle_Ptotal[iMarker_NacelleExhaust];
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

double CConfig::GetIsothermal_Temperature(string val_marker) {
	unsigned short iMarker_Isothermal;
	for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++)
		if (Marker_Isothermal[iMarker_Isothermal] == val_marker) break;
	return Isothermal_Temperature[iMarker_Isothermal];
}

double CConfig::GetWall_HeatFlux(string val_marker) {
	unsigned short iMarker_HeatFlux;
	for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++)
		if (Marker_HeatFlux[iMarker_HeatFlux] == val_marker) break;
	return Heat_Flux[iMarker_HeatFlux];
}

double CConfig::GetFanFace_Mach_Target(string val_marker) {
	unsigned short iMarker_NacelleInflow;
	for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++)
		if (Marker_NacelleInflow[iMarker_NacelleInflow] == val_marker) break;
	return FanFace_Mach_Target[iMarker_NacelleInflow];
}

double CConfig::GetFanFace_Pressure(string val_marker) {
	unsigned short iMarker_NacelleInflow;
	for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++)
		if (Marker_NacelleInflow[iMarker_NacelleInflow] == val_marker) break;
	return FanFace_Pressure[iMarker_NacelleInflow];
}

double CConfig::GetFanFace_Mach(string val_marker) {
	unsigned short iMarker_NacelleInflow;
	for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++)
		if (Marker_NacelleInflow[iMarker_NacelleInflow] == val_marker) break;
	return FanFace_Mach[iMarker_NacelleInflow];
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

double CConfig::GetMotion_Ramp(unsigned long val_iter) {
	double coeff;

	if ((val_iter < int(Motion_Ramp[0])) ) {
		coeff = 0.0;
	} else if (((val_iter >= int(Motion_Ramp[0])) && (int(Motion_Ramp[1]) == 0)) ||
			(val_iter >= int(Motion_Ramp[0])+int(Motion_Ramp[1]))) {
		coeff = 1.0;
	} else if ((val_iter >= int(Motion_Ramp[0])) &&
			(val_iter <  int(Motion_Ramp[0])+int(Motion_Ramp[1]))) {
		coeff = ((double)val_iter - Motion_Ramp[0])/Motion_Ramp[1];
	} else
		coeff = 1.0;

	return coeff;
}

void CConfig::SetNondimensionalization(unsigned short val_nDim, unsigned short val_iZone) {
	double Mach2Vel_FreeStream, ModVel_FreeStream, Energy_FreeStream = 0.0, ModVel_FreeStreamND;
	double Velocity_Reynolds;
	unsigned short iDim;
	int rank = MASTER_NODE;
  
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	Velocity_FreeStreamND = new double[val_nDim];
  
	/*--- Local variables and memory allocation ---*/
	double Alpha = AoA*PI_NUMBER/180.0;
	double Beta  = AoS*PI_NUMBER/180.0;
	double Gamma_Minus_One = Gamma - 1.0;
	bool Compressible = (!Incompressible);
	bool Unsteady = (Unsteady_Simulation != NO);
	bool turbulent = (Kind_Solver == RANS);
  
	if (Unsteady_Farfield && Unsteady) {
		double time, time_ext;
		unsigned long iter, n;
		string text_line;
		double temp_Temperature, temp_Mach, temp_Pressure, Derivative_Begin, Derivative_End;
		vector<double> Time_Spline, Temperature_Spline, Mach_Spline, Pressure_Spline;
		vector<double> Temperature2_Spline, Mach2_Spline, Pressure2_Spline;
    
		ifstream farfield_file;
		string filename = GetFarfield_FileName();
		farfield_file.open(filename.data(), ios::in);
    
		/*--- In case there is no restart file ---*/
		if (farfield_file.fail()) {
			cout << "There is no farfield bounadry data file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}
    
		/*--- The first line is the header ---*/
		getline (farfield_file, text_line);
    
		while (getline (farfield_file,text_line) ) {
			istringstream point_line(text_line);
      
			/*--- First value is the point index, then the conservative vars. ---*/
			point_line >> time;
			point_line >> temp_Temperature;// Temperature_FreeStream;
			point_line >> temp_Mach; // Mach
			point_line >> temp_Pressure;// Pressure_FreeStream;
      
			/*--- Spline interpolation, dummy function ---*/
			Time_Spline.push_back(time);
			Temperature_Spline.push_back(temp_Temperature);
			Mach_Spline.push_back(temp_Mach);
			Pressure_Spline.push_back(temp_Pressure);
			cout << " pressure = " << temp_Pressure << endl;
		}
    
		/*--- Close the restart file ---*/
		farfield_file.close();
    
		/*--- Create the spline ---*/
		n = Time_Spline.size();
		Temperature2_Spline.resize(n);
		Mach2_Spline.resize(n);
		Pressure2_Spline.resize(n);
    
		Derivative_Begin = (Temperature_Spline[1]-Temperature_Spline[0])/(Time_Spline[1]-Time_Spline[0]);
		Derivative_End = (Temperature_Spline[n-1]-Temperature_Spline[n-2])/(Time_Spline[n-1]-Time_Spline[n-2]);
		SetSpline(Time_Spline, Temperature_Spline, n, Derivative_Begin, Derivative_End, Temperature2_Spline);
    
    
		Derivative_Begin = (Mach_Spline[1]-Mach_Spline[0])/(Time_Spline[1]-Time_Spline[0]);
		Derivative_End = (Mach_Spline[n-1]-Mach_Spline[n-2])/(Time_Spline[n-1]-Time_Spline[n-2]);
		SetSpline(Time_Spline, Mach_Spline, n, Derivative_Begin, Derivative_End, Mach2_Spline);
    
    
		Derivative_Begin = (Pressure_Spline[1]-Pressure_Spline[0])/(Time_Spline[1]-Time_Spline[0]);
		Derivative_End = (Pressure_Spline[n-1]-Pressure_Spline[n-2])/(Time_Spline[n-1]-Time_Spline[n-2]);
		SetSpline(Time_Spline, Pressure_Spline, n, Derivative_Begin, Derivative_End, Pressure2_Spline);
    
		/*--- Allocate arrays to store time dependent farfield information ---*/
		Density_FreeStreamND_Time   = new double  [nExtIter];
		Pressure_FreeStreamND_Time  = new double  [nExtIter];
		Mach_Inf_Time 				= new double  [nExtIter];
		Energy_FreeStreamND_Time    = new double  [nExtIter];
		Velocity_FreeStreamND_Time  = new double* [nExtIter];
		for (iter = 0; iter < nExtIter; iter ++) {
			Velocity_FreeStreamND_Time[iter]  = new double [val_nDim];
			for (iDim = 0; iDim < val_nDim; iDim ++)
				Velocity_FreeStreamND_Time[iter][iDim] = 0.0;
		}
    
		/*--- Loop over all external time instances to store the time dependent values ---*/
    
		for (iter = 0; iter < nExtIter; iter ++) {
      
			time_ext =  iter*Delta_UnstTime;
      
			/*--- Get the interpolated values at current time step ---*/
			Temperature_FreeStream = GetSpline(Time_Spline, Temperature_Spline, Temperature2_Spline, n, time_ext);
			Mach = GetSpline(Time_Spline, Mach_Spline, Mach2_Spline, n, time_ext);
			Pressure_FreeStream = GetSpline(Time_Spline, Pressure_Spline, Pressure2_Spline, n, time_ext);
      
      
			Mach2Vel_FreeStream = sqrt(Gamma*Gas_Constant*Temperature_FreeStream);
      
			/*--- Compute the Free Stream velocity, using the Mach number ---*/
			if (val_nDim == 2) {
				Velocity_FreeStream[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
				Velocity_FreeStream[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
			}
			if (val_nDim == 3) {
				Velocity_FreeStream[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
				Velocity_FreeStream[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
				Velocity_FreeStream[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
			}
      
			Velocity_Ref = sqrt(Pressure_Ref/Density_Ref);
      
			for (iDim = 0; iDim < val_nDim; iDim++)
				Velocity_FreeStreamND_Time[iter][iDim] = Velocity_FreeStream[iDim]/Velocity_Ref;
      
			Density_FreeStream  				= Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
			Density_FreeStreamND_Time[iter]  	= Density_FreeStream/Density_Ref;
			Pressure_FreeStreamND_Time[iter] 	= Pressure_FreeStream/Pressure_Ref;
			Mach_Inf_Time[iter]				    = Mach;
      
			ModVel_FreeStreamND = 0;
			for (iDim = 0; iDim < val_nDim; iDim++)
				ModVel_FreeStreamND += Velocity_FreeStreamND_Time[iter][iDim]*Velocity_FreeStreamND_Time[iter][iDim];
			ModVel_FreeStreamND    	 = sqrt(ModVel_FreeStreamND);
			Energy_FreeStreamND_Time[iter]    = Pressure_FreeStreamND_Time[iter]/(Density_FreeStreamND_Time[iter]*Gamma_Minus_One)+0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
		}
    
	}
  
	if (Compressible) {
    
		Mach2Vel_FreeStream = sqrt(Gamma*Gas_Constant*Temperature_FreeStream);
    
		/*--- Compute the Free Stream velocity, using the Mach number ---*/
		if (val_nDim == 2) {
			Velocity_FreeStream[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
			Velocity_FreeStream[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
		}
		if (val_nDim == 3) {
			Velocity_FreeStream[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
			Velocity_FreeStream[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
			Velocity_FreeStream[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
		}
    
		/*--- Compute the modulus of the free stream velocity ---*/
		ModVel_FreeStream = 0;
		for (iDim = 0; iDim < val_nDim; iDim++)
			ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
		ModVel_FreeStream = sqrt(ModVel_FreeStream);
    
		if (Viscous) {
      
			/*--- First, check if there is mesh motion. If yes, use the Mach
       number relative to the body to initialize the flow. ---*/
			if (Rotating_Frame || Grid_Movement)
				Velocity_Reynolds = Mach_Motion*Mach2Vel_FreeStream;
			else
				Velocity_Reynolds = ModVel_FreeStream;
      
			/*--- For viscous flows, pressure will be computed from a density
       that is found from the Reynolds number. The viscosity is computed
       from the dimensional version of Sutherland's law ---*/
			Viscosity_FreeStream = 1.853E-5*(pow(Temperature_FreeStream/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_FreeStream+110.3));
			Density_FreeStream   = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*Length_Reynolds);
			Pressure_FreeStream  = Density_FreeStream*Gas_Constant*Temperature_FreeStream;
		} else {
			/*--- For inviscid flow, density is calculated from the specified
			 total temperature and pressure using the gas law. ---*/
			Density_FreeStream  = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
		}
		/*-- Compute the freestream energy. ---*/
		Energy_FreeStream = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+0.5*ModVel_FreeStream*ModVel_FreeStream;
    
		/*--- Additional reference values defined by Pref, Tref, RHOref. By definition,
		 Lref is one because we have converted the grid to meters.---*/
		Length_Ref         = 1.0;
		Velocity_Ref      = sqrt(Pressure_Ref/Density_Ref);
		Time_Ref          = Length_Ref/Velocity_Ref;
		Omega_Ref         = Velocity_Ref/Length_Ref;
		Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;
		Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/Temperature_Ref;
		Viscosity_Ref     = Density_Ref*Velocity_Ref*Length_Ref;
		Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);
    
	}
  
	else {
    
		/*--- Reference length = 1 (by default)
     Reference density = liquid density or freestream
     Reference viscosity = liquid viscosity or freestream
     Reference velocity = liquid velocity or freestream
     Reference pressure = Reference density * Reference velocity * Reference velocity
     Reynolds number based on the liquid or reference viscosity ---*/
    
		Pressure_FreeStream = 0.0;
		Length_Ref = 1.0;
		Density_Ref = Density_FreeStream;
		ModVel_FreeStream = 0;
		for (iDim = 0; iDim < val_nDim; iDim++)
			ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
		ModVel_FreeStream = sqrt(ModVel_FreeStream);
		Velocity_Ref = ModVel_FreeStream;
		Pressure_Ref = Density_Ref*(Velocity_Ref*Velocity_Ref);
    
		if (Viscous) {
			Reynolds = Density_Ref*Velocity_Ref*Length_Ref / Viscosity_FreeStream;
			Viscosity_Ref = Viscosity_FreeStream * Reynolds;
		}
    
		/*--- Compute mach number ---*/
		Mach = ModVel_FreeStream / sqrt(Bulk_Modulus/Density_FreeStream);
		if (val_nDim == 2) AoA = atan(Velocity_FreeStream[1]/Velocity_FreeStream[0])*180.0/PI_NUMBER;
		else AoA = atan(Velocity_FreeStream[2]/Velocity_FreeStream[0])*180.0/PI_NUMBER;
		if (val_nDim == 2) AoS = 0.0;
		else AoS = asin(Velocity_FreeStream[1]/ModVel_FreeStream)*180.0/PI_NUMBER;
    
		Froude = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);
    
		Time_Ref = Length_Ref/Velocity_Ref;
    
	}
  
	/*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
	Pressure_FreeStreamND = Pressure_FreeStream/Pressure_Ref;
	Density_FreeStreamND  = Density_FreeStream/Density_Ref;
  
	for (iDim = 0; iDim < val_nDim; iDim++)
		Velocity_FreeStreamND[iDim] = Velocity_FreeStream[iDim]/Velocity_Ref;
	Temperature_FreeStreamND = Temperature_FreeStream/Temperature_Ref;
  
	Gas_ConstantND = Gas_Constant/Gas_Constant_Ref;
  
	/*--- Perform non-dim. for rotating terms ---*/
	Omega_Mag = 0.0;
	for (iDim = 0; iDim < 3; iDim++) {
		Omega_FreeStreamND[iDim] = Omega[iDim]/Omega_Ref;
		Omega_Mag += Omega[iDim]*Omega[iDim];
	}
	Omega_Mag = sqrt(Omega_Mag)/Omega_Ref;
  
	ModVel_FreeStreamND = 0;
	for (iDim = 0; iDim < val_nDim; iDim++)
		ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
	ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND);
	Energy_FreeStreamND    = Pressure_FreeStreamND/(Density_FreeStreamND*Gamma_Minus_One)+0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
	Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;
	Total_UnstTimeND = Total_UnstTime / Time_Ref;
	Delta_UnstTimeND = Delta_UnstTime / Time_Ref;
  
	double kine_Inf  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*TurbulenceIntensity_FreeStream*TurbulenceIntensity_FreeStream);
	double omega_Inf = Density_FreeStreamND*kine_Inf/(Viscosity_FreeStreamND*Turb2LamViscRatio_FreeStream);
  
	/*--- Write output to the console if this is the master node and first domain ---*/
	if ((rank == MASTER_NODE) && (val_iZone == 0)) {
    
		cout << endl <<"---------------- Flow & Non-dimensionalization information ---------------" << endl;
    
		cout.precision(6);
    
		if (Compressible) {
			if (Viscous) {
				cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
				cout << "based on the freestream temperature and a density computed" << endl;
				cout << "from the Reynolds number." << endl;
			} else {
				cout << "Inviscid flow: Computing density based on freestream" << endl;
				cout << "temperature and pressure using the ideal gas law." << endl;
			}
		}
		else {
			cout << "Viscous and Inviscid flow: rho_ref, and vel_ref" << endl;
			cout << "are based on the freestream values, p_ref = rho_ref*vel_ref^2." << endl;
			cout << "The freestream value of the pressure is 0." << endl;
			cout << "Mach number: "<< Mach << ", computed using the Bulk modulus." << endl;
			cout << "Angle of attack (deg): "<< AoA << ", computed using the the free-stream velocity." << endl;
			cout << "Side slip angle (deg): "<< AoS << ", computed using the the free-stream velocity." << endl;
			if (Viscous) cout << "Reynolds number: " << Reynolds << ", computed using free-stream values."<<endl;
			cout << "Only dimensional computation, the grid should be dimensional." << endl;
		}
    
		cout <<"--Input conditions:"<< endl;
		cout << "Grid conversion factor to meters: " << Conversion_Factor << endl;
    
		if (Compressible) {
			cout << "Ratio of specific heats: " << Gamma           << endl;
			cout << "Specific gas constant (J/(kg.K)): "   << Gas_Constant  << endl;
		}
		else {
			cout << "Bulk modulus (N/m^2): "						<< Bulk_Modulus    << endl;
			cout << "Artificial compressibility factor (N/m^2): "						<< ArtComp_Factor    << endl;
		}
    
		cout << "Freestream pressure (N/m^2): "          << Pressure_FreeStream    << endl;
		if (Compressible)
			cout << "Freestream temperature (K): "       << Temperature_FreeStream << endl;
		cout << "Freestream density (kg/m^3): "					 << Density_FreeStream << endl;
		if (val_nDim == 2) {
			cout << "Freestream velocity (m/s): (" << Velocity_FreeStream[0] << ",";
			cout << Velocity_FreeStream[1] << ")";
		} else if (val_nDim == 3) {
			cout << "Freestream velocity (m/s): (" << Velocity_FreeStream[0] << ",";
			cout << Velocity_FreeStream[1] << "," << Velocity_FreeStream[2] << ")";
		}
    
		cout << " -> Modulus: "					 << ModVel_FreeStream << endl;
    
		if (Compressible)
			cout << "Freestream energy (kg.m/s^2): "					 << Energy_FreeStream << endl;
    
		if (Viscous)
			cout << "Freestream viscosity (N.s/m^2): "				 << Viscosity_FreeStream << endl;
    
		if (Rotating_Frame) {
			cout << "Freestream rotation (rad/s): (" << Omega[0];
			cout << "," << Omega[1] << "," << Omega[2] << ")" << endl;
		}
		if (Unsteady) {
			cout << "Total time (s): " << Total_UnstTime << ". Time step (s): " << Delta_UnstTime << endl;
		}
    
		/*--- Print out reference values. ---*/
		cout <<"--Reference values:"<< endl;
		cout << "Reference pressure (N/m^2): "      << Pressure_Ref    << endl;
    
		if (Compressible) {
			cout << "Reference temperature (K): "   << Temperature_Ref << endl;
			cout << "Reference energy (kg.m/s^2): "       << Energy_FreeStream/Energy_FreeStreamND     << endl;
		}
		else {
			cout << "Reference length (m): 1.0" << endl;
		}
		cout << "Reference density (kg/m^3): "       << Density_Ref     << endl;
		cout << "Reference velocity (m/s): "       << Velocity_Ref     << endl;
    
		if (Viscous)
			cout << "Reference viscosity (N.s/m^2): "       << Viscosity_Ref     << endl;
    
		if (Unsteady)
			cout << "Reference time (s): "        << Time_Ref      << endl;
    
		/*--- Print out resulting non-dim values here. ---*/
		cout << "--Resulting non-dimensional state:" << endl;
		cout << "Mach number (non-dimensional): " << Mach << endl;
		if (Viscous) {
			cout << "Reynolds number (non-dimensional): " << Reynolds << endl;
			cout << "Reynolds length (m): "       << Length_Reynolds     << endl;
		}
		cout << "Froude number (non-dimensional): " << Froude << endl;
		cout << "Lenght of the baseline wave (non-dimensional): " << 2.0*PI_NUMBER*Froude*Froude << endl;
    
		if (Compressible) {
			cout << "Negative pressure, temperature or density is not allowed!" << endl;
			cout << "Specific gas constant (non-dimensional): "   << Gas_Constant << endl;
			cout << "Freestream temperature (non-dimensional): "  << Temperature_FreeStreamND << endl;
		}
		cout << "Freestream pressure (non-dimensional): "     << Pressure_FreeStreamND    << endl;
		cout << "Freestream density (non-dimensional): "      << Density_FreeStreamND     << endl;
		if (val_nDim == 2) {
			cout << "Freestream velocity (non-dimensional): (" << Velocity_FreeStreamND[0] << ",";
			cout << Velocity_FreeStreamND[1] << ")";
		} else if (val_nDim == 3) {
			cout << "Freestream velocity (non-dimensional): (" << Velocity_FreeStreamND[0] << ",";
			cout << Velocity_FreeStreamND[1] << "," << Velocity_FreeStreamND[2] << ")";
		}
		cout << " -> Modulus: "					 << ModVel_FreeStreamND << endl;
    
		if (turbulent){
			cout << "Free-stream turb. kinetic energy (non-dimensional): " << kine_Inf << endl;
			cout << "Free-stream specific dissipation (non-dimensional): " << omega_Inf << endl;
		}
    
		if (Compressible)
			cout << "Freestream energy (non-dimensional): "					 << Energy_FreeStreamND << endl;
    
		if (Viscous)
			cout << "Freestream viscosity (non-dimensional): " << Viscosity_FreeStreamND << endl;
    
		if (Rotating_Frame) {
			cout << "Freestream rotation (non-dimensional): (" << Omega_FreeStreamND[0];
			cout << "," << Omega_FreeStreamND[1] << "," << Omega_FreeStreamND[2] << ")" << endl;
		}
		if (Unsteady) {
			cout << "Total time (non-dimensional): "				 << Total_UnstTimeND << endl;
			cout << "Time step (non-dimensional): "				 << Delta_UnstTimeND << endl;
		}
		if (Mach <= 0.0) cout << "Force coefficients computed using reference values." << endl;
		else cout << "Force coefficients computed using freestream values." << endl;
	}
  
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
