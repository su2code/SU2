/*!
 * \file config_structure.cpp
 * \brief Main file for managing the config file
 * \author F. Palacios, T. Economon, B. Tracey, H. Kline
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = MASTER_NODE;
#endif

  /*--- Initialize pointers to Null---*/

  SetPointersNull();

  /*--- Reading config options  ---*/

  SetConfig_Options(val_iZone, val_nZone);

  /*--- Parsing the config file  ---*/

  SetConfig_Parsing(case_filename);

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, val_iZone, val_nDim);

  /*--- Configuration file boundaries/markers setting ---*/

  SetMarkers(val_software);

  /*--- Configuration file output ---*/

  if ((rank == MASTER_NODE) && (verb_level == VERB_HIGH) && (val_iZone == 0))
    SetOutput(val_software, val_iZone);

}

CConfig::CConfig(char case_filename[MAX_STRING_SIZE], unsigned short val_software) {

  /*--- Initialize pointers to Null---*/

  SetPointersNull();

  /*--- Reading config options  ---*/

  SetConfig_Options(0, 1);

  /*--- Parsing the config file  ---*/

  SetConfig_Parsing(case_filename);

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, 0, 1);

}

CConfig::CConfig(char case_filename[MAX_STRING_SIZE], CConfig *config) {

  bool runtime_file = false;

  /*--- Initialize pointers to Null---*/

  SetPointersNull();

  /*--- Reading config options  ---*/

  SetRunTime_Options();

  /*--- Parsing the config file  ---*/

  runtime_file = SetRunTime_Parsing(case_filename);

  /*--- Update original config file ---*/

  if (runtime_file) {
    config->SetnExtIter(nExtIter);
  }

}

SU2_Comm CConfig::GetMPICommunicator() {

  return SU2_Communicator;

}

void CConfig::SetMPICommunicator(SU2_Comm Communicator) {

  SU2_Communicator = Communicator;

}

unsigned short CConfig::GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nZone = 1; // Default value
  unsigned short iLine, nLine = 10;
  char cstr[200];
  string::size_type position;

  /*--- Search the mesh file for the 'NZONE' keyword. ---*/

  switch (val_format) {
    case SU2:

      /*--- Open grid file ---*/

      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      if (mesh_file.fail()) {
        cout << "cstr=" << cstr << endl;
        cout << "There is no geometry file (GetnZone))!" << endl;

#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }

      /*--- Read the SU2 mesh file ---*/

      for (iLine = 0; iLine < nLine ; iLine++) {

        getline (mesh_file, text_line);

        /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/

        position = text_line.find ("NZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nZone = atoi(text_line.c_str());
        }
      }

      break;

  }

  /*--- For harmonic balance integration, nZones = nTimeInstances. ---*/

  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE && (config->GetKind_SU2() != SU2_DEF)   ) {
  	nZone = config->GetnTimeInstances();
  }

  return (unsigned short) nZone;
}

unsigned short CConfig::GetnDim(string val_mesh_filename, unsigned short val_format) {

  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nDim = 3;
  unsigned short iLine, nLine = 10;
  char cstr[200];
  string::size_type position;

  /*--- Open grid file ---*/

  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);

  switch (val_format) {
  case SU2:

    /*--- Read SU2 mesh file ---*/

    for (iLine = 0; iLine < nLine ; iLine++) {

      getline (mesh_file, text_line);

      /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/

      position = text_line.find ("NDIME=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nDim = atoi(text_line.c_str());
      }
    }
    break;

  case CGNS:

#ifdef HAVE_CGNS

    /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/

    int fn, nbases = 0, nzones = 0, file_type;
    int cell_dim = 0, phys_dim = 0;
    char basename[CGNS_STRING_SIZE];

    /*--- Check whether the supplied file is truly a CGNS file. ---*/

    if ( cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK ) {
      printf( "\n\n   !!! Error !!!\n" );
      printf( " %s is not a CGNS file.\n", val_mesh_filename.c_str());
      printf( " Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }

    /*--- Open the CGNS file for reading. The value of fn returned
       is the specific index number for this file and will be
       repeatedly used in the function calls. ---*/

    if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn)) cg_error_exit();

    /*--- Get the number of databases. This is the highest node
       in the CGNS heirarchy. ---*/

    if (cg_nbases(fn, &nbases)) cg_error_exit();

    /*--- Check if there is more than one database. Throw an
       error if there is because this reader can currently
       only handle one database. ---*/

    if ( nbases > 1 ) {
      printf("\n\n   !!! Error !!!\n" );
      printf("CGNS reader currently incapable of handling more than 1 database.");
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }

    /*--- Read the databases. Note that the indexing starts at 1. ---*/

    for ( int i = 1; i <= nbases; i++ ) {

      if (cg_base_read(fn, i, basename, &cell_dim, &phys_dim)) cg_error_exit();

      /*--- Get the number of zones for this base. ---*/

      if (cg_nzones(fn, i, &nzones)) cg_error_exit();

    }

    /*--- Set the problem dimension as read from the CGNS file ---*/

    nDim = cell_dim;

#endif

    break;

  }

  mesh_file.close();

  return (unsigned short) nDim;
}

void CConfig::SetPointersNull(void) {
  
  Marker_CfgFile_Out_1D       = NULL;   Marker_All_Out_1D        = NULL;
  Marker_CfgFile_GeoEval      = NULL;   Marker_All_GeoEval       = NULL;
  Marker_CfgFile_Monitoring   = NULL;   Marker_All_Monitoring    = NULL;
  Marker_CfgFile_Designing    = NULL;   Marker_All_Designing     = NULL;
  Marker_CfgFile_Plotting     = NULL;   Marker_All_Plotting      = NULL;
  Marker_CfgFile_Analyze      = NULL;   Marker_All_Analyze       = NULL;
  Marker_CfgFile_DV           = NULL;   Marker_All_DV            = NULL;
  Marker_CfgFile_Moving       = NULL;   Marker_All_Moving        = NULL;
  Marker_CfgFile_PerBound     = NULL;   Marker_All_PerBound      = NULL;    Marker_PerBound   = NULL;
  Marker_CfgFile_FSIinterface = NULL;
  
  Marker_DV                   = NULL;   Marker_Moving            = NULL;    Marker_Monitoring = NULL;
  Marker_Designing            = NULL;   Marker_GeoEval           = NULL;    Marker_Plotting   = NULL;
  Marker_Analyze              = NULL;
  Marker_CfgFile_KindBC       = NULL;   Marker_All_KindBC        = NULL;
  
  /*--- Marker Pointers ---*/

  Marker_Euler                = NULL;    Marker_FarField         = NULL;    Marker_Custom         = NULL;
  Marker_SymWall              = NULL;    Marker_Pressure         = NULL;    Marker_PerBound       = NULL;
  Marker_PerDonor             = NULL;    Marker_NearFieldBound   = NULL;    Marker_InterfaceBound = NULL;
  Marker_Dirichlet            = NULL;    Marker_Inlet            = NULL;    
  Marker_Supersonic_Inlet     = NULL;    Marker_Outlet           = NULL;    Marker_Out_1D         = NULL;
  Marker_Isothermal           = NULL;    Marker_HeatFlux         = NULL;    Marker_EngineInflow   = NULL;
  Marker_Supersonic_Outlet    = NULL;    Marker_Load             = NULL;
  Marker_EngineExhaust        = NULL;    Marker_Displacement     = NULL;    Marker_Load           = NULL;
  Marker_Load_Dir             = NULL;    Marker_Load_Sine        = NULL;    Marker_Clamped        = NULL;
  Marker_FlowLoad             = NULL;    Marker_Neumann          = NULL;    Marker_Internal       = NULL;
  Marker_All_TagBound         = NULL;    Marker_CfgFile_TagBound = NULL;    Marker_All_KindBC     = NULL;
  Marker_CfgFile_KindBC       = NULL;    Marker_All_SendRecv     = NULL;    Marker_All_PerBound   = NULL;
  Marker_FSIinterface         = NULL;    Marker_All_FSIinterface = NULL;    Marker_Riemann        = NULL;
  Marker_Fluid_InterfaceBound = NULL;

  
  /*--- Boundary Condition settings ---*/

  Dirichlet_Value = NULL;    Isothermal_Temperature = NULL;
  Heat_Flux       = NULL;    Displ_Value            = NULL;    Load_Value = NULL;
  FlowLoad_Value  = NULL;
  
  /*--- Inlet Outlet Boundary Condition settings ---*/

  Inlet_Ttotal    = NULL;    Inlet_Ptotal      = NULL;
  Inlet_FlowDir   = NULL;    Inlet_Temperature = NULL;    Inlet_Pressure = NULL;
  Inlet_Velocity  = NULL;
  Outlet_Pressure = NULL;
  
  /*--- Engine Boundary Condition settings ---*/
  
  Inflow_Pressure      = NULL;    Inflow_MassFlow    = NULL;    Inflow_ReverseMassFlow  = NULL;
  Inflow_TotalPressure = NULL;    Inflow_Temperature = NULL;    Inflow_TotalTemperature = NULL;
  Inflow_RamDrag       = NULL;    Inflow_Force       = NULL;    Inflow_Power            = NULL;
  Inflow_Mach          = NULL;
  
  Exhaust_Pressure        = NULL;   Exhaust_Temperature        = NULL;    Exhaust_MassFlow = NULL;
  Exhaust_TotalPressure   = NULL;   Exhaust_TotalTemperature   = NULL;
  Exhaust_GrossThrust     = NULL;   Exhaust_Force              = NULL;
  Exhaust_Power           = NULL;   Exhaust_Temperature_Target = NULL;
  Exhaust_Pressure_Target = NULL;
  
  Engine_Mach  = NULL;    Engine_Force        = NULL;
  Engine_Power = NULL;    Engine_NetThrust    = NULL;    Engine_GrossThrust = NULL;
  Engine_Area  = NULL;    EngineInflow_Target = NULL;
  
  Periodic_Translate   = NULL;   Periodic_Rotation  = NULL;   Periodic_Center    = NULL;
  Periodic_Translation = NULL;   Periodic_RotAngles = NULL;   Periodic_RotCenter = NULL;

  Dirichlet_Value           = NULL;     Exhaust_Temperature_Target	= NULL;	    Exhaust_Temperature   = NULL;
  Exhaust_Pressure_Target   = NULL;		Inlet_Ttotal                = NULL;	    Inlet_Ptotal          = NULL;
  Inlet_FlowDir             = NULL;     Inlet_Temperature           = NULL;     Inlet_Pressure        = NULL;
  Inlet_Velocity            = NULL;     Inflow_Mach                 = NULL;     Inflow_Pressure       = NULL;
  Exhaust_Pressure          = NULL;     Outlet_Pressure             = NULL;     Isothermal_Temperature= NULL;
  Heat_Flux                 = NULL;     Displ_Value                 = NULL;     Load_Value            = NULL;
  FlowLoad_Value            = NULL;     Periodic_RotCenter          = NULL;     Periodic_RotAngles    = NULL;
  Periodic_Translation      = NULL;     Periodic_Center             = NULL;     Periodic_Rotation     = NULL;
  Periodic_Translate        = NULL;

  Load_Dir            = NULL;    Load_Dir_Value      = NULL;    Load_Dir_Multiplier = NULL;
  Load_Sine_Dir       = NULL;    Load_Sine_Amplitude = NULL;    Load_Sine_Frequency = NULL;

  /*--- Actuator Disk Boundary Condition settings ---*/
  
  ActDiskInlet_Pressure         = NULL;    ActDiskInlet_TotalPressure = NULL;    ActDiskInlet_Temperature = NULL;
  ActDiskInlet_TotalTemperature = NULL;    ActDiskInlet_MassFlow      = NULL;    ActDiskInlet_RamDrag     = NULL;
  ActDiskInlet_Force            = NULL;    ActDiskInlet_Power         = NULL;

  ActDiskOutlet_Pressure      = NULL;
  ActDiskOutlet_TotalPressure = NULL;   ActDiskOutlet_GrossThrust = NULL;  ActDiskOutlet_Force            = NULL;
  ActDiskOutlet_Power         = NULL;   ActDiskOutlet_Temperature = NULL;  ActDiskOutlet_TotalTemperature = NULL;
  ActDiskOutlet_MassFlow      = NULL;
  
  ActDisk_DeltaPress      = NULL;    ActDisk_DeltaTemp      = NULL;
  ActDisk_TotalPressRatio = NULL;    ActDisk_TotalTempRatio = NULL;    ActDisk_StaticPressRatio = NULL;
  ActDisk_StaticTempRatio = NULL;    ActDisk_NetThrust      = NULL;    ActDisk_GrossThrust      = NULL;
  ActDisk_Power           = NULL;    ActDisk_MassFlow       = NULL;    ActDisk_Area             = NULL;
  ActDisk_ReverseMassFlow = NULL;    Surface_MassFlow       = NULL;    Surface_DC60             = NULL;    Surface_IDC = NULL;
  Surface_IDC_Mach        = NULL;    Surface_IDR            = NULL;    ActDisk_Mach             = NULL;
  ActDisk_Force           = NULL;    ActDisk_BCThrust       = NULL;    ActDisk_BCThrust_Old     = NULL;
  
  /*--- Miscellaneous/unsorted ---*/

  Aeroelastic_plunge  = NULL;
  Aeroelastic_pitch   = NULL;
  MassFrac_FreeStream = NULL;
  Velocity_FreeStream = NULL;

  RefOriginMoment     = NULL;
  CFL_AdaptParam      = NULL;            
  CFL                 = NULL;
  HTP_Axis = NULL;
  PlaneTag            = NULL;
  Kappa_Flow	      = NULL;    
  Kappa_AdjFlow       = NULL;
  Section_WingBounds    = NULL;
  ParamDV             = NULL;     
  DV_Value            = NULL;    
  Design_Variable     = NULL;

  Hold_GridFixed_Coord= NULL;
  SubsonicEngine_Cyl  = NULL;
  EA_IntLimit         = NULL;
  RK_Alpha_Step       = NULL;
  MG_CorrecSmooth     = NULL;
  MG_PreSmooth        = NULL;
  MG_PostSmooth       = NULL;
  Int_Coeffs          = NULL;

  Kind_ObjFunc   = NULL;

  Weight_ObjFunc = NULL;

  /*--- Moving mesh pointers ---*/

  Kind_GridMovement	  = NULL;    LocationStations	  = NULL;
  Motion_Origin_X     = NULL;    Motion_Origin_Y     = NULL;    Motion_Origin_Z	    = NULL;
  Translation_Rate_X  = NULL;    Translation_Rate_Y  = NULL;    Translation_Rate_Z  = NULL;
  Rotation_Rate_X     = NULL;    Rotation_Rate_Y     = NULL;    Rotation_Rate_Z     = NULL;
  Pitching_Omega_X    = NULL;    Pitching_Omega_Y    = NULL;    Pitching_Omega_Z    = NULL;
  Pitching_Ampl_X     = NULL;    Pitching_Ampl_Y     = NULL;    Pitching_Ampl_Z     = NULL;
  Pitching_Phase_X    = NULL;    Pitching_Phase_Y    = NULL;    Pitching_Phase_Z    = NULL;
  Plunging_Omega_X    = NULL;    Plunging_Omega_Y    = NULL;    Plunging_Omega_Z    = NULL;
  Plunging_Ampl_X     = NULL;    Plunging_Ampl_Y     = NULL;    Plunging_Ampl_Z     = NULL;
  RefOriginMoment_X   = NULL;    RefOriginMoment_Y   = NULL;    RefOriginMoment_Z   = NULL;
  MoveMotion_Origin   = NULL;
  Periodic_Translate  = NULL;    Periodic_Rotation 	 = NULL;    Periodic_Center	    = NULL;
  Periodic_Translation= NULL;    Periodic_RotAngles	 = NULL;    Periodic_RotCenter  = NULL;


  /* Harmonic Balance Frequency pointer */
  Omega_HB = NULL;
    
  /*--- Initialize some default arrays to NULL. ---*/
  
  default_vel_inf       = NULL;
  default_ffd_axis      = NULL;
  default_eng_cyl       = NULL;
  default_eng_val       = NULL;
  default_heat_capacity = NULL;
  default_cfl_adapt     = NULL;
  default_ad_coeff_flow = NULL;
  default_ad_coeff_adj  = NULL;
  default_obj_coeff     = NULL;
  default_geo_loc       = NULL;
  default_distortion    = NULL;
  default_ea_lim        = NULL;
  default_grid_fix      = NULL;
  default_inc_crit      = NULL;
  default_htp_axis      = NULL;
  default_body_force    = NULL;

  Riemann_FlowDir= NULL;
  NRBC_FlowDir = NULL;
  CoordFFDBox= NULL;
  DegreeFFDBox= NULL;
  FFDTag = NULL;
  nDV_Value = NULL;
  TagFFDBox = NULL;
 
  Kind_Data_Riemann      = NULL;
  Riemann_Var1           = NULL;
  Riemann_Var2           = NULL;
  Kind_Data_NRBC         = NULL;
  NRBC_Var1              = NULL;
  NRBC_Var2              = NULL;
  Marker_TurboBoundIn    = NULL;
  Marker_TurboBoundOut   = NULL;

  Kind_TurboPerformance  = NULL;

  Marker_NRBC            = NULL;
  
  /*--- Variable initialization ---*/
  
  ExtIter    = 0;
  IntIter    = 0;
  nIntCoeffs = 0;
  FSIIter    = 0;
  
  AoA_Offset = 0;
  AoS_Offset = 0;

  nMarker_PerBound = 0;
  nPeriodic_Index  = 0;

  Grid_Movement = false;
  Aeroelastic_Simulation = false;
  
}

void CConfig::SetRunTime_Options(void) {
  
  /* DESCRIPTION: Number of external iterations */
  
  addUnsignedLongOption("EXT_ITER", nExtIter, 999999);

}

void CConfig::SetConfig_Options(unsigned short val_iZone, unsigned short val_nZone) {
  
  nZone = val_nZone;
  iZone = val_iZone;

  /*--- Allocate some default arrays needed for lists of doubles. ---*/
  
  default_vel_inf       = new su2double[3];
  default_ffd_axis      = new su2double[3];
  default_eng_cyl       = new su2double[7];
  default_eng_val       = new su2double[5];
  default_heat_capacity = new su2double[5];
  default_cfl_adapt     = new su2double[4];
  default_ad_coeff_flow = new su2double[3];
  default_ad_coeff_adj  = new su2double[3];
  default_obj_coeff     = new su2double[5];
  default_geo_loc       = new su2double[2];
  default_distortion    = new su2double[2];
  default_ea_lim        = new su2double[3];
  default_grid_fix      = new su2double[6];
  default_inc_crit      = new su2double[3];
  default_htp_axis      = new su2double[2];
  default_body_force    = new su2double[3];

  // This config file is parsed by a number of programs to make it easy to write SU2
  // wrapper scripts (in python, go, etc.) so please do
  // the best you can to follow the established format. It's very hard to parse c++ code
  // and none of us that write the parsers want to write a full c++ interpreter. Please
  // play nice with the existing format so that you don't break the existing scripts.

  /* BEGIN_CONFIG_OPTIONS */

  /*!\par CONFIG_CATEGORY: Problem Definition \ingroup Config */
  /*--- Options related to problem definition and partitioning ---*/

  /*!\brief REGIME_TYPE \n  DESCRIPTION: Regime type \n OPTIONS: see \link Regime_Map \endlink \ingroup Config*/
  addEnumOption("REGIME_TYPE", Kind_Regime, Regime_Map, COMPRESSIBLE);
  
  /*!\brief PHYSICAL_PROBLEM \n DESCRIPTION: Physical governing equations \n Options: see \link Solver_Map \endlink \n DEFAULT: NO_SOLVER \ingroup Config*/
  addEnumOption("PHYSICAL_PROBLEM", Kind_Solver, Solver_Map, NO_SOLVER);
  /*!\brief MATH_PROBLEM  \n DESCRIPTION: Mathematical problem \n  Options: DIRECT, ADJOINT \ingroup Config*/
  addMathProblemOption("MATH_PROBLEM", ContinuousAdjoint, false, DiscreteAdjoint, false, Restart_Flow, false);
  /*!\brief KIND_TURB_MODEL \n DESCRIPTION: Specify turbulence model \n Options: see \link Turb_Model_Map \endlink \n DEFAULT: NO_TURB_MODEL \ingroup Config*/
  addEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, NO_TURB_MODEL);
  /*!\brief KIND_2PHASE_MODEL \n DESCRIPTION: Specify 2phase model \n Options: DEFAULT: NO_2PHASE_MODEL \ingroup Config*/
  addEnumOption("KIND_2PHASE_MODEL", Kind_2phase_Model, Two_phase_Model_Map, NO_2PHASE_MODEL);

  /*!\brief KIND_NUCLEATION_MODEL \n DESCRIPTION: Specify nucleation model \n Options: DEFAULT: classical theory \ingroup Config*/
  addEnumOption("KIND_NUCLEATION_MODEL", Kind_Nucleation_Model, Nucleation_Model_Map, CLASSICAL_THEORY);

  /*!\brief KIND_LIQUID_MODEL \n DESCRIPTION: Specify liquid model \n Options: DEFAULT: WATER \ingroup Config*/
  addEnumOption("KIND_LIQUID_MODEL", Kind_Liquid_Model, Liquid_Model_Map, WATER);

  /*!\brief KIND_TRANS_MODEL \n DESCRIPTION: Specify transition model OPTIONS: see \link Trans_Model_Map \endlink \n DEFAULT: NO_TRANS_MODEL \ingroup Config*/
  addEnumOption("KIND_TRANS_MODEL", Kind_Trans_Model, Trans_Model_Map, NO_TRANS_MODEL);
  /*\brief AXISYMMETRIC \n DESCRIPTION: Axisymmetric simulation \n DEFAULT: false \ingroup Config */
  addBoolOption("AXISYMMETRIC", Axisymmetric, false);
  /* DESCRIPTION: Add the gravity force */
  addBoolOption("GRAVITY_FORCE", GravityForce, false);
  /* DESCRIPTION: Apply a body force as a source term (NO, YES) */
  addBoolOption("BODY_FORCE", Body_Force, false);
  default_body_force[0] = 0.0; default_body_force[1] = 0.0; default_body_force[2] = 0.0;
  /* DESCRIPTION: Vector of body force values (BodyForce_X, BodyForce_Y, BodyForce_Z) */
  addDoubleArrayOption("BODY_FORCE_VECTOR", 3, Body_Force_Vector, default_body_force);
  /* DESCRIPTION: Perform a low fidelity simulation */
  addBoolOption("LOW_FIDELITY_SIMULATION", LowFidelitySim, false);
  /*!\brief RESTART_SOL \n DESCRIPTION: Restart solution from native solution file \n Options: NO, YES \ingroup Config */
  addBoolOption("RESTART_SOL", Restart, false);
  /*!\brief UPDATE_RESTART_PARAMS \n DESCRIPTION: Update some parameters from a metadata file when restarting \n Options: NO, YES \ingroup Config */
  addBoolOption("UPDATE_RESTART_PARAMS", Update_Restart_Params, false);
  /*!\brief BINARY_RESTART \n DESCRIPTION: Read / write binary SU2 native restart files. \n Options: YES, NO \ingroup Config */
  addBoolOption("WRT_BINARY_RESTART", Wrt_Binary_Restart, true);
  /*!\brief BINARY_RESTART \n DESCRIPTION: Read / write binary SU2 native restart files. \n Options: YES, NO \ingroup Config */
  addBoolOption("READ_BINARY_RESTART", Read_Binary_Restart, true);
  /*!\brief SYSTEM_MEASUREMENTS \n DESCRIPTION: System of measurements \n OPTIONS: see \link Measurements_Map \endlink \n DEFAULT: SI \ingroup Config*/
  addEnumOption("SYSTEM_MEASUREMENTS", SystemMeasurements, Measurements_Map, SI);

  /*!\par CONFIG_CATEGORY: FluidModel \ingroup Config*/
  /*!\brief FLUID_MODEL \n DESCRIPTION: Fluid model \n OPTIONS: See \link FluidModel_Map \endlink \n DEFAULT: STANDARD_AIR \ingroup Config*/
  addEnumOption("FLUID_MODEL", Kind_FluidModel, FluidModel_Map, STANDARD_AIR);

  /*!\par CONFIG_CATEGORY: Freestream Conditions \ingroup Config*/
  /*--- Options related to freestream specification ---*/

  /*!\brief GAS_CONSTANT \n DESCRIPTION: Specific gas constant (287.058 J/kg*K (air), only for compressible flows) \ingroup Config*/
  addDoubleOption("GAS_CONSTANT", Gas_Constant, 287.058);
  /*!\brief GAMMA_VALUE  \n DESCRIPTION: Ratio of specific heats (1.4 (air), only for compressible flows) \ingroup Config*/
  addDoubleOption("GAMMA_VALUE", Gamma, 1.4);

  /*!\brief CONSTANT_GAMMA  \n DESCRIPTION: if yes the gas is polytropic (constant gamma) \ingroup Config*/
  addBoolOption("CONSTANT_GAMMA", Constant_Gamma, true);

  /*!\brief AUTORESET_NEGATIVE_SOL  \n DESCRIPTION: if true two phase solutions are forced to be positive \ingroup Config*/
  addBoolOption("AUTORESET_NEGATIVE_SOL", AutoReset_NegativeSol, true);

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
  /*!\brief VISCOSITY_MODEL \n DESCRIPTION: model of the viscosity \n OPTIONS: See \link ViscosityModel_Map \endlink \n DEFAULT: SUTHERLAND \ingroup Config*/
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

  /*!\brief REYNOLDS_NUMBER \n DESCRIPTION: Reynolds number (non-dimensional, based on the free-stream values). Needed for viscous solvers. For incompressible solvers the Reynolds length will always be 1.0 \n DEFAULT: 0.0 \ingroup Config */
  addDoubleOption("REYNOLDS_NUMBER", Reynolds, 0.0);
  /*!\brief REYNOLDS_LENGTH \n DESCRIPTION: Reynolds length (1 m by default). Used for compressible solver: incompressible solver will use 1.0. \ingroup Config */
  addDoubleOption("REYNOLDS_LENGTH", Length_Reynolds, 1.0);
  /*!\brief PRANDTL_LAM \n DESCRIPTION: Laminar Prandtl number (0.72 (air), only for compressible flows) \n DEFAULT: 0.72 \ingroup Config*/
  addDoubleOption("PRANDTL_LAM", Prandtl_Lam, 0.72);
  /*!\brief PRANDTL_TURB \n DESCRIPTION: Turbulent Prandtl number (0.9 (air), only for compressible flows) \n DEFAULT 0.90 \ingroup Config*/
  addDoubleOption("PRANDTL_TURB", Prandtl_Turb, 0.90);
  /*!\brief BULK_MODULUS \n DESCRIPTION: Value of the Bulk Modulus  \n DEFAULT 1.42E5 \ingroup Config*/
  addDoubleOption("BULK_MODULUS", Bulk_Modulus, 1.42E5);
  /* DESCRIPTION: Artifical compressibility factor  */
  addDoubleOption("ARTCOMP_FACTOR", ArtComp_Factor, 1.0);
  /*!\brief MACH_NUMBER  \n DESCRIPTION:  Mach number (non-dimensional, based on the free-stream values). 0.0 by default \ingroup Config*/
  addDoubleOption("MACH_NUMBER", Mach, 0.0);
  /*!\brief INIT_OPTION \n DESCRIPTION: Init option to choose between Reynolds or thermodynamics quantities for initializing the solution \n OPTIONS: see \link InitOption_Map \endlink \n DEFAULT REYNOLDS \ingroup Config*/
  addEnumOption("INIT_OPTION", Kind_InitOption, InitOption_Map, REYNOLDS);
  /* DESCRIPTION: Free-stream option to choose between density and temperature for initializing the solution */
  addEnumOption("FREESTREAM_OPTION", Kind_FreeStreamOption, FreeStreamOption_Map, TEMPERATURE_FS);
  /*!\brief FREESTREAM_PRESSURE\n DESCRIPTION: Free-stream pressure (101325.0 N/m^2 by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_PRESSURE", Pressure_FreeStream, 101325.0);
  /*!\brief FREESTREAM_DENSITY\n DESCRIPTION: Free-stream density (1.2886 Kg/m^3 (air), 998.2 Kg/m^3 (water)) \n DEFAULT -1.0 (calculated from others) \ingroup Config*/
  addDoubleOption("FREESTREAM_DENSITY", Density_FreeStream, -1.0);
  /*!\brief FREESTREAM_TEMPERATURE\n DESCRIPTION: Free-stream temperature (288.15 K by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_TEMPERATURE", Temperature_FreeStream, 288.15);
  /*!\brief FREESTREAM_TEMPERATURE_VE\n DESCRIPTION: Free-stream vibrational-electronic temperature (288.15 K by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_TEMPERATURE_VE", Temperature_ve_FreeStream, 288.15);
  default_vel_inf[0] = 1.0; default_vel_inf[1] = 0.0; default_vel_inf[2] = 0.0;
  /*!\brief FREESTREAM_VELOCITY\n DESCRIPTION: Free-stream velocity (m/s) */
  addDoubleArrayOption("FREESTREAM_VELOCITY", 3, Velocity_FreeStream, default_vel_inf);
  /* DESCRIPTION: Free-stream viscosity (1.853E-5 Ns/m^2 (air), 0.798E-3 Ns/m^2 (water)) */
  addDoubleOption("FREESTREAM_VISCOSITY", Viscosity_FreeStream, -1.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_INTERMITTENCY", Intermittency_FreeStream, 1.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_TURBULENCEINTENSITY", TurbulenceIntensity_FreeStream, 0.05);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_NU_FACTOR", NuFactor_FreeStream, 3.0);
  /* DESCRIPTION:  */
  addDoubleOption("ENGINE_NU_FACTOR", NuFactor_Engine, 3.0);
  /* DESCRIPTION:  */
  addDoubleOption("ACTDISK_SECONDARY_FLOW", SecondaryFlow_ActDisk, 0.0);
  /* DESCRIPTION:  */
  addDoubleOption("INITIAL_BCTHRUST", Initial_BCThrust, 4000.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_TURB2LAMVISCRATIO", Turb2LamViscRatio_FreeStream, 10.0);
  /* DESCRIPTION: Side-slip angle (degrees, only for compressible flows) */
  addDoubleOption("SIDESLIP_ANGLE", AoS, 0.0);
  /*!\brief AOA  \n DESCRIPTION: Angle of attack (degrees, only for compressible flows) \ingroup Config*/
  addDoubleOption("AOA", AoA, 0.0);
  /* DESCRIPTION: Activate fixed CL mode (specify a CL instead of AoA). */
  addBoolOption("FIXED_CL_MODE", Fixed_CL_Mode, false);
  /* DESCRIPTION: Activate fixed CM mode (specify a CM instead of iH). */
  addBoolOption("FIXED_CM_MODE", Fixed_CM_Mode, false);
  /* DESCRIPTION: Evaluate the dCD_dCL or dCD_dCMy during run time. */
  addBoolOption("EVAL_DCD_DCX", Eval_dCD_dCX, true);
  /* DESCRIPTION: DIscard the angle of attack in the solution and the increment in the geometry files. */
  addBoolOption("DISCARD_INFILES", Discard_InFiles, false);
  /* DESCRIPTION: Specify a fixed coefficient of lift instead of AoA (only for compressible flows) */
  addDoubleOption("TARGET_CL", Target_CL, 0.0);
  /* DESCRIPTION: Specify a fixed coefficient of lift instead of AoA (only for compressible flows) */
  addDoubleOption("TARGET_CM", Target_CM, 0.0);
  /* DESCRIPTION: Damping factor for fixed CL mode. */
  addDoubleOption("DCL_DALPHA", dCL_dAlpha, 0.2);
  /* DESCRIPTION: Damping factor for fixed CL mode. */
  addDoubleOption("DCM_DIH", dCM_diH, 0.05);
  /* DESCRIPTION: Number of times Alpha is updated in a fix CL problem. */
  addUnsignedLongOption("UPDATE_ALPHA", Update_Alpha, 5);
  /* DESCRIPTION: Number of times Alpha is updated in a fix CL problem. */
  addUnsignedLongOption("UPDATE_IH", Update_iH, 5);
  /* DESCRIPTION: Damping factor for fixed CL mode. */
  addDoubleOption("DNETTHRUST_DBCTHRUST", dNetThrust_dBCThrust, 2.0);
  /* DESCRIPTION: Number of times Alpha is updated in a fix CL problem. */
  addUnsignedLongOption("UPDATE_BCTHRUST", Update_BCThrust, 5);


  /*!\par CONFIG_CATEGORY: Reference Conditions \ingroup Config*/
  /*--- Options related to reference values for nondimensionalization ---*/

  Length_Ref = 1.0; //<---- NOTE: this should be given an option or set as a const

  /*!\brief REF_ORIGIN_MOMENT_X\n DESCRIPTION: X Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_X", nRefOriginMoment_X, RefOriginMoment_X);
  /*!\brief REF_ORIGIN_MOMENT_Y\n DESCRIPTION: Y Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_Y", nRefOriginMoment_Y, RefOriginMoment_Y);
  /*!\brief REF_ORIGIN_MOMENT_Z\n DESCRIPTION: Z Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_Z", nRefOriginMoment_Z, RefOriginMoment_Z);
  /*!\brief REF_AREA\n DESCRIPTION: Reference area for force coefficients (0 implies automatic calculation) \ingroup Config*/
  addDoubleOption("REF_AREA", RefAreaCoeff, 1.0);
  /*!\brief SEMI_SPAN\n DESCRIPTION: Wing semi-span (1 by deafult) \ingroup Config*/
  addDoubleOption("SEMI_SPAN", SemiSpan, 1.0);
  /*!\brief REF_LENGTH_MOMENT\n DESCRIPTION: Reference length for pitching, rolling, and yawing non-dimensional moment \ingroup Config*/
  addDoubleOption("REF_LENGTH_MOMENT", RefLengthMoment, 1.0);
  /*!\brief REF_ELEM_LENGTH\n DESCRIPTION: Reference element length for computing the slope limiter epsilon \ingroup Config*/
  addDoubleOption("REF_ELEM_LENGTH", RefElemLength, 0.1);
  /*!\brief REF_SHARP_EDGES\n DESCRIPTION: Reference coefficient for detecting sharp edges \ingroup Config*/
  addDoubleOption("REF_SHARP_EDGES", RefSharpEdges, 3.0);
	/*!\brief REF_VELOCITY\n DESCRIPTION: Reference velocity (incompressible only)  \ingroup Config*/
  addDoubleOption("REF_VELOCITY", Velocity_Ref, -1.0);
	/* !\brief REF_VISCOSITY  \n DESCRIPTION: Reference viscosity (incompressible only)  \ingroup Config*/
  addDoubleOption("REF_VISCOSITY", Viscosity_Ref, -1.0);
  /* DESCRIPTION: Type of mesh motion */
  addEnumOption("REF_DIMENSIONALIZATION", Ref_NonDim, NonDim_Map, DIMENSIONAL);

  /*!\par CONFIG_CATEGORY: Boundary Markers \ingroup Config*/
  /*--- Options related to various boundary markers ---*/

  /*!\brief HTP_AXIS\n DESCRIPTION: Location of the HTP axis*/
  default_htp_axis[0] = 0.0; default_htp_axis[1] = 0.0;
  addDoubleArrayOption("HTP_AXIS", 2, HTP_Axis, default_htp_axis);
  /*!\brief MARKER_PLOTTING\n DESCRIPTION: Marker(s) of the surface in the surface flow solution file  \ingroup Config*/
  addStringListOption("MARKER_PLOTTING", nMarker_Plotting, Marker_Plotting);
  /*!\brief MARKER_MONITORING\n DESCRIPTION: Marker(s) of the surface where evaluate the non-dimensional coefficients \ingroup Config*/
  addStringListOption("MARKER_MONITORING", nMarker_Monitoring, Marker_Monitoring);
  /*!\brief MARKER_CONTROL_VOLUME\n DESCRIPTION: Marker(s) of the surface in the surface flow solution file  \ingroup Config*/
  addStringListOption("MARKER_ANALYZE", nMarker_Analyze, Marker_Analyze);
  /*!\brief MARKER_DESIGNING\n DESCRIPTION: Marker(s) of the surface where objective function (design problem) will be evaluated \ingroup Config*/
  addStringListOption("MARKER_DESIGNING", nMarker_Designing, Marker_Designing);
  /*!\brief MARKER_OUT_1D \n DESCRIPTION: Outlet boundary marker(s) over which to calculate 1-D flow properties
   Format: ( outlet marker) \ingroup Config*/
  addStringListOption("MARKER_OUT_1D", nMarker_Out_1D, Marker_Out_1D);
  /*!\brief GEO_MARKER\n DESCRIPTION: Marker(s) of the surface where evaluate the geometrical functions \ingroup Config*/
  addStringListOption("GEO_MARKER", nMarker_GeoEval, Marker_GeoEval);
  /*!\brief MARKER_EULER\n DESCRIPTION: Euler wall boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_EULER", nMarker_Euler, Marker_Euler);
  /*!\brief MARKER_FAR\n DESCRIPTION: Far-field boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_FAR", nMarker_FarField, Marker_FarField);
  /*!\brief MARKER_SYM\n DESCRIPTION: Symmetry boundary condition \ingroup Config*/
  addStringListOption("MARKER_SYM", nMarker_SymWall, Marker_SymWall);
  /*!\brief MARKER_PRESSURE\n DESCRIPTION: Symmetry boundary condition \ingroup Config*/
  addStringListOption("MARKER_PRESSURE", nMarker_Pressure, Marker_Pressure);
  /*!\brief MARKER_NEARFIELD\n DESCRIPTION: Near-Field boundary condition \ingroup Config*/
  addStringListOption("MARKER_NEARFIELD", nMarker_NearFieldBound, Marker_NearFieldBound);
  /*!\brief MARKER_FLUID_INTERFACE\n DESCRIPTION: Fluid interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_FLUID_INTERFACE", nMarker_Fluid_InterfaceBound, Marker_Fluid_InterfaceBound);
  /*!\brief MARKER_INTERFACE\n DESCRIPTION: Zone interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_INTERFACE", nMarker_InterfaceBound, Marker_InterfaceBound);
  /*!\brief MARKER_FSI_INTERFACE \n DESCRIPTION: FSI interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_FSI_INTERFACE", nMarker_FSIinterface, Marker_FSIinterface);
  /*!\brief MARKER_DIRICHLET  \n DESCRIPTION: Dirichlet boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_DIRICHLET", nMarker_Dirichlet, Marker_Dirichlet);
  /* DESCRIPTION: Neumann boundary marker(s) */
  addStringListOption("MARKER_NEUMANN", nMarker_Neumann, Marker_Neumann);
  /* DESCRIPTION: Neumann boundary marker(s) */
  addStringListOption("MARKER_INTERNAL", nMarker_Internal, Marker_Internal);
  /* DESCRIPTION: Custom boundary marker(s) */
  addStringListOption("MARKER_CUSTOM", nMarker_Custom, Marker_Custom);
  /* DESCRIPTION: Periodic boundary marker(s) for use with SU2_MSH
   Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
   rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
   rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) */
  addPeriodicOption("MARKER_PERIODIC", nMarker_PerBound, Marker_PerBound, Marker_PerDonor,
                    Periodic_RotCenter, Periodic_RotAngles, Periodic_Translation);

  /*!\brief ACTDISK_TYPE  \n DESCRIPTION: Actuator Disk boundary type \n OPTIONS: see \link ActDisk_Map \endlink \n Default: VARIABLES_JUMP \ingroup Config*/
  addEnumOption("ACTDISK_TYPE", Kind_ActDisk, ActDisk_Map, VARIABLES_JUMP);

  /*!\brief MARKER_ACTDISK\n DESCRIPTION: Periodic boundary marker(s) for use with SU2_MSH
   Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
   rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
   rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) \ingroup Config*/
  addActDiskOption("MARKER_ACTDISK",
                   nMarker_ActDiskInlet, nMarker_ActDiskOutlet,  Marker_ActDiskInlet, Marker_ActDiskOutlet,
                   ActDisk_PressJump, ActDisk_TempJump, ActDisk_Omega);

  /*!\brief INLET_TYPE  \n DESCRIPTION: Inlet boundary type \n OPTIONS: see \link Inlet_Map \endlink \n DEFAULT: TOTAL_CONDITIONS \ingroup Config*/
  addEnumOption("INLET_TYPE", Kind_Inlet, Inlet_Map, TOTAL_CONDITIONS);

  /*!\brief MARKER_INLET  \n DESCRIPTION: Inlet boundary marker(s) with the following formats,
   Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector.
   Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector. \ingroup Config*/
  addInletOption("MARKER_INLET", nMarker_Inlet, Marker_Inlet, Inlet_Ttotal, Inlet_Ptotal, Inlet_FlowDir);

  /*!\brief MARKER_RIEMANN \n DESCRIPTION: Riemann boundary marker(s) with the following formats, a unit vector.
   * \n OPTIONS: See \link Riemann_Map \endlink. The variables indicated by the option and the flow direction unit vector must be specified. \ingroup Config*/
  addRiemannOption("MARKER_RIEMANN", nMarker_Riemann, Marker_Riemann, Kind_Data_Riemann, Riemann_Map, Riemann_Var1, Riemann_Var2, Riemann_FlowDir);
  /*!\brief MARKER_NRBC \n DESCRIPTION: Riemann boundary marker(s) with the following formats, a unit vector. \ingroup Config*/
  addNRBCOption("MARKER_NRBC", nMarker_NRBC, Marker_NRBC, Kind_Data_NRBC, NRBC_Map, NRBC_Var1, NRBC_Var2, NRBC_FlowDir);
  /*!\brief MIXING_PROCESS_TYPE \n DESCRIPTION: types of mixing process for averaging quantities at the boundaries.
    \n OPTIONS: see \link MixingProcess_Map \endlink \n DEFAULT: AREA_AVERAGE \ingroup Config*/
  addEnumOption("MIXING_PROCESS_TYPE", Kind_MixingProcess, MixingProcess_Map, AREA_AVERAGE);
  /*!\brief MARKER_MIXINGPLANE \n DESCRIPTION: Identify the boundaries in which the mixing plane is applied. \ingroup Config*/
  addMixingPlaneOption("MARKER_MIXINGPLANE", nMarker_MixBound, Marker_MixBound, Marker_MixDonor);
  /*!\brief MARKER_MIXINGPLANE \n DESCRIPTION: Identify the boundaries in which the mixing plane is applied. \ingroup Config*/
  addTurboPerfOption("MARKER_TURBO_PERFORMANCE", nMarker_TurboPerf, Marker_TurboBoundIn, Marker_TurboBoundOut, Kind_TurboPerformance, TurboPerformance_Map);
  /*!\brief MARKER_SUPERSONIC_INLET  \n DESCRIPTION: Supersonic inlet boundary marker(s)
   * \n   Format: (inlet marker, temperature, static pressure, velocity_x,   velocity_y, velocity_z, ... ), i.e. primitive variables specified. \ingroup Config*/
  addInletOption("MARKER_SUPERSONIC_INLET", nMarker_Supersonic_Inlet, Marker_Supersonic_Inlet, Inlet_Temperature, Inlet_Pressure, Inlet_Velocity);
  /*!\brief MARKER_SUPERSONIC_OUTLET \n DESCRIPTION: Supersonic outlet boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_SUPERSONIC_OUTLET", nMarker_Supersonic_Outlet, Marker_Supersonic_Outlet);
  /*!\brief MARKER_OUTLET  \n DESCRIPTION: Outlet boundary marker(s)\n
   Format: ( outlet marker, back pressure (static), ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_OUTLET", nMarker_Outlet, Marker_Outlet, Outlet_Pressure);
  /*!\brief MARKER_ISOTHERMAL DESCRIPTION: Isothermal wall boundary marker(s)\n
   * Format: ( isothermal marker, wall temperature (static), ... ) \ingroup Config  */
  addStringDoubleListOption("MARKER_ISOTHERMAL", nMarker_Isothermal, Marker_Isothermal, Isothermal_Temperature);
  /*!\brief MARKER_HEATFLUX  \n DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_HEATFLUX", nMarker_HeatFlux, Marker_HeatFlux, Heat_Flux);
  /*!\brief MARKER_ENGINE_INFLOW  \n DESCRIPTION: Engine inflow boundary marker(s)
   Format: ( nacelle inflow marker, fan face Mach, ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_ENGINE_INFLOW", nMarker_EngineInflow, Marker_EngineInflow, EngineInflow_Target);
  /* DESCRIPTION: Highlite area */
  addDoubleOption("HIGHLITE_AREA", Highlite_Area, 1.0);
  /* DESCRIPTION: Fan poly efficiency */
  addDoubleOption("FAN_POLY_EFF", Fan_Poly_Eff, 1.0);
  /*!\brief SUBSONIC_ENGINE\n DESCRIPTION: Engine subsonic intake region \ingroup Config*/
  addBoolOption("SUBSONIC_ENGINE", SubsonicEngine, false);
  /* DESCRIPTION: Actuator disk double surface */
  addBoolOption("ACTDISK_DOUBLE_SURFACE", ActDisk_DoubleSurface, false);
  /* DESCRIPTION: Only half engine is in the computational grid */
  addBoolOption("ENGINE_HALF_MODEL", Engine_HalfModel, false);
  /* DESCRIPTION: Actuator disk double surface */
  addBoolOption("ACTDISK_SU2_DEF", ActDisk_SU2_DEF, false);
  /* DESCRIPTION: Definition of the distortion rack (radial number of proves / circumferential density (degree) */
  default_distortion[0] =  5.0; default_distortion[1] =  15.0;
  addDoubleArrayOption("DISTORTION_RACK", 2, DistortionRack, default_distortion);
  /* DESCRIPTION: Values of the box to impose a subsonic nacellle (mach, Pressure, Temperature) */
  default_eng_val[0]=0.0; default_eng_val[1]=0.0; default_eng_val[2]=0.0;
  default_eng_val[3]=0.0;  default_eng_val[4]=0.0;
  addDoubleArrayOption("SUBSONIC_ENGINE_VALUES", 5, SubsonicEngine_Values, default_eng_val);

  default_heat_capacity[0]= 0.0; default_heat_capacity[1]=0.0; default_heat_capacity[2]=0.0;
  default_heat_capacity[3]=0.0;  default_heat_capacity[4]=0.0;
  addDoubleArrayOption("HEAT_CAPACITY_MODEL", 5, Coeff_HeatCapacity, default_heat_capacity);

  addDoubleOption("MOLECULAR_MASS", Molecular_Mass, 18.0);

  /* DESCRIPTION: Coordinates of the box to impose a subsonic nacellle cylinder (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax, Radius) */
  default_eng_cyl[0] = 0.0; default_eng_cyl[1] = 0.0; default_eng_cyl[2] = 0.0;
  default_eng_cyl[3] =  1E15; default_eng_cyl[4] =  1E15; default_eng_cyl[5] =  1E15; default_eng_cyl[6] =  1E15;
  addDoubleArrayOption("SUBSONIC_ENGINE_CYL", 7, SubsonicEngine_Cyl, default_eng_cyl);
  /* DESCRIPTION: Engine exhaust boundary marker(s)
   Format: (nacelle exhaust marker, total nozzle temp, total nozzle pressure, ... )*/
  addExhaustOption("MARKER_ENGINE_EXHAUST", nMarker_EngineExhaust, Marker_EngineExhaust, Exhaust_Temperature_Target, Exhaust_Pressure_Target);
  /* DESCRIPTION: Clamped boundary marker(s) */
  addStringListOption("MARKER_CLAMPED", nMarker_Clamped, Marker_Clamped);
  /* DESCRIPTION: Displacement boundary marker(s) */
  addStringDoubleListOption("MARKER_NORMAL_DISPL", nMarker_Displacement, Marker_Displacement, Displ_Value);
  /* DESCRIPTION: Load boundary marker(s) */
  addStringDoubleListOption("MARKER_NORMAL_LOAD", nMarker_Load, Marker_Load, Load_Value);
  /* DESCRIPTION: Load boundary marker(s)
   Format: (inlet marker, load, multiplier, dir_x, dir_y, dir_z, ... ), i.e. primitive variables specified. */
  addInletOption("MARKER_LOAD", nMarker_Load_Dir, Marker_Load_Dir, Load_Dir_Value, Load_Dir_Multiplier, Load_Dir);
  /* DESCRIPTION: Sine load boundary marker(s)
   Format: (inlet marker, load, multiplier, dir_x, dir_y, dir_z, ... ), i.e. primitive variables specified. */
  addInletOption("MARKER_SINE_LOAD", nMarker_Load_Sine, Marker_Load_Sine, Load_Sine_Amplitude, Load_Sine_Frequency, Load_Sine_Dir);

  /* DESCRIPTION: Flow load boundary marker(s) */
  addStringDoubleListOption("MARKER_FLOWLOAD", nMarker_FlowLoad, Marker_FlowLoad, FlowLoad_Value);
  /* DESCRIPTION: Damping factor for engine inlet condition */
  addDoubleOption("DAMP_ENGINE_INFLOW", Damp_Engine_Inflow, 0.95);
  /* DESCRIPTION: Damping factor for engine exhaust condition */
  addDoubleOption("DAMP_ENGINE_EXHAUST", Damp_Engine_Exhaust, 0.95);
  /*!\brief ENGINE_INFLOW_TYPE  \n DESCRIPTION: Inlet boundary type \n OPTIONS: see \link Engine_Inflow_Map \endlink \n Default: FAN_FACE_MACH \ingroup Config*/
  addEnumOption("ENGINE_INFLOW_TYPE", Kind_Engine_Inflow, Engine_Inflow_Map, FAN_FACE_MACH);
  /* DESCRIPTION: Evaluate a problem with engines */
  addBoolOption("ENGINE", Engine, false);


  /*!\par CONFIG_CATEGORY: Time-marching \ingroup Config*/
  /*--- Options related to time-marching ---*/

  /* DESCRIPTION: Unsteady simulation  */
  addEnumOption("UNSTEADY_SIMULATION", Unsteady_Simulation, Unsteady_Map, STEADY);
  /* DESCRIPTION:  Courant-Friedrichs-Lewy condition of the finest grid */
  addDoubleOption("CFL_NUMBER", CFLFineGrid, 1.25);
  /* DESCRIPTION:  Max time step in local time stepping simulations */
  addDoubleOption("MAX_DELTA_TIME", Max_DeltaTime, 1000000);
  /* DESCRIPTION: Activate The adaptive CFL number. */
  addBoolOption("CFL_ADAPT", CFL_Adapt, false);
  /* !\brief CFL_ADAPT_PARAM
   * DESCRIPTION: Parameters of the adaptive CFL number (factor down, factor up, CFL limit (min and max) )
   * Factor down generally >1.0, factor up generally < 1.0 to cause the CFL to increase when residual is decreasing,
   * and decrease when the residual is increasing or stalled. \ingroup Config*/
  default_cfl_adapt[0] = 0.0; default_cfl_adapt[1] = 0.0; default_cfl_adapt[2] = 1.0; default_cfl_adapt[3] = 100.0;
  addDoubleArrayOption("CFL_ADAPT_PARAM", 4, CFL_AdaptParam, default_cfl_adapt);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the adjoint problem */
  addDoubleOption("CFL_REDUCTION_ADJFLOW", CFLRedCoeff_AdjFlow, 0.8);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
  addDoubleOption("CFL_REDUCTION_TURB", CFLRedCoeff_Turb, 1.0);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
  addDoubleOption("CFL_REDUCTION_2PHASE", CFLRedCoeff_2phase, 1.0);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the turbulent adjoint problem */
  addDoubleOption("CFL_REDUCTION_ADJTURB", CFLRedCoeff_AdjTurb, 1.0);
  /* DESCRIPTION: Number of total iterations */
  addUnsignedLongOption("EXT_ITER", nExtIter, 999999);
  /* DESCRIPTION: External iteration offset due to restart */
  addUnsignedLongOption("EXT_ITER_OFFSET", ExtIter_OffSet, 0);
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
  /* DESCRIPTION: Integer number of periodic time instances for Harmonic Balance */
  addUnsignedShortOption("TIME_INSTANCES", nTimeInstances, 1);
  /* DESCRIPTION: Time period for Harmonic Balance wihtout moving meshes */
  addDoubleOption("HB_PERIOD", HarmonicBalance_Period, -1.0);
  /* DESCRIPTION: Iteration number to begin unsteady restarts (dual time method) */
  addLongOption("UNST_RESTART_ITER", Unst_RestartIter, 0);
  /* DESCRIPTION: Starting direct solver iteration for the unsteady adjoint */
  addLongOption("UNST_ADJOINT_ITER", Unst_AdjointIter, 0);
  /* DESCRIPTION: Number of iterations to average the objective */
  addLongOption("ITER_AVERAGE_OBJ", Iter_Avg_Objective , 0);
  /* DESCRIPTION: Iteration number to begin unsteady restarts (structural analysis) */
  addLongOption("DYN_RESTART_ITER", Dyn_RestartIter, 0);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FLOW", Kind_TimeIntScheme_Flow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJFLOW", Kind_TimeIntScheme_AdjFlow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_TURB", Kind_TimeIntScheme_Turb, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_2PHASE", Kind_TimeIntScheme_2phase, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJTURB", Kind_TimeIntScheme_AdjTurb, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_WAVE", Kind_TimeIntScheme_Wave, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FEA", Kind_TimeIntScheme_FEA, Time_Int_Map_FEA, NEWMARK_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_HEAT", Kind_TimeIntScheme_Heat, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_POISSON", Kind_TimeIntScheme_Poisson, Time_Int_Map, EULER_IMPLICIT);

  /*!\par CONFIG_CATEGORY: Linear solver definition \ingroup Config*/
  /*--- Options related to the linear solvers ---*/

  /*!\brief LINEAR_SOLVER
   *  \n DESCRIPTION: Linear solver for the implicit, mesh deformation, or discrete adjoint systems \n OPTIONS: see \link Linear_Solver_Map \endlink \n DEFAULT: FGMRES \ingroup Config*/
  addEnumOption("LINEAR_SOLVER", Kind_Linear_Solver, Linear_Solver_Map, FGMRES);
  /*!\brief LINEAR_SOLVER_PREC
   *  \n DESCRIPTION: Preconditioner for the Krylov linear solvers \n OPTIONS: see \link Linear_Solver_Prec_Map \endlink \n DEFAULT: LU_SGS \ingroup Config*/
  addEnumOption("LINEAR_SOLVER_PREC", Kind_Linear_Solver_Prec, Linear_Solver_Prec_Map, LU_SGS);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("LINEAR_SOLVER_ERROR", Linear_Solver_Error, 1E-5);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("LINEAR_SOLVER_ITER", Linear_Solver_Iter, 10);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("LINEAR_SOLVER_RESTART_FREQUENCY", Linear_Solver_Restart_Frequency, 10);
  /* DESCRIPTION: Relaxation of the flow equations solver for the implicit formulation */
  addDoubleOption("RELAXATION_FACTOR_FLOW", Relaxation_Factor_Flow, 1.0);
  /* DESCRIPTION: Relaxation of the turb equations solver for the implicit formulation */
  addDoubleOption("RELAXATION_FACTOR_TURB", Relaxation_Factor_Turb, 1.0);
  /* DESCRIPTION: Relaxation of the 2phase equations solver for the implicit formulation */
    addDoubleOption("RELAXATION_FACTOR_2PHASE", Relaxation_Factor_2phase, 1.0);

  /* DESCRIPTION: Relaxation of the adjoint flow equations solver for the implicit formulation */
  addDoubleOption("RELAXATION_FACTOR_ADJFLOW", Relaxation_Factor_AdjFlow, 1.0);
  /* DESCRIPTION: Roe coefficient */
  addDoubleOption("ROE_KAPPA", Roe_Kappa, 0.5);
  /* DESCRIPTION: Roe-Turkel preconditioning for low Mach number flows */
  addBoolOption("ROE_TURKEL_PREC", Low_Mach_Precon, false);
  /* DESCRIPTION: Post-reconstruction correction for low Mach number flows */
  addBoolOption("LOW_MACH_CORR", Low_Mach_Corr, false);
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
  addDoubleOption("ENTROPY_FIX_COEFF", EntropyFix_Coeff, 0.001);
  /* DESCRIPTION: Linear solver for the discete adjoint systems */
  addEnumOption("DISCADJ_LIN_SOLVER", Kind_DiscAdj_Linear_Solver, Linear_Solver_Map, FGMRES);
  /* DESCRIPTION: Preconditioner for the discrete adjoint Krylov linear solvers */
  addEnumOption("DISCADJ_LIN_PREC", Kind_DiscAdj_Linear_Prec, Linear_Solver_Prec_Map, ILU);
  
  /*!\par CONFIG_CATEGORY: Convergence\ingroup Config*/
  /*--- Options related to convergence ---*/
  
  /*!\brief CONV_CRITERIA
   *  \n DESCRIPTION: Convergence criteria \n OPTIONS: see \link Converge_Crit_Map \endlink \n DEFAULT: RESIDUAL \ingroup Config*/
  addEnumOption("CONV_CRITERIA", ConvCriteria, Converge_Crit_Map, RESIDUAL);
  /*!\brief RESIDUAL_REDUCTION \n DESCRIPTION: Residual reduction (order of magnitude with respect to the initial value)\n DEFAULT: 3.0 \ingroup Config*/
  addDoubleOption("RESIDUAL_REDUCTION", OrderMagResidual, 3.0);
  /*!\brief RESIDUAL_MINVAL\n DESCRIPTION: Min value of the residual (log10 of the residual)\n DEFAULT: -8.0 \ingroup Config*/
  addDoubleOption("RESIDUAL_MINVAL", MinLogResidual, -8.0);
  /* DESCRIPTION: Residual reduction (order of magnitude with respect to the initial value) */
  addDoubleOption("RESIDUAL_REDUCTION_FSI", OrderMagResidualFSI, 3.0);
  /* DESCRIPTION: Min value of the residual (log10 of the residual) */
  addDoubleOption("RESIDUAL_MINVAL_FSI", MinLogResidualFSI, -5.0);
  /* DESCRIPTION: FEM: UTOL = norm(Delta_U(k)) / norm(U(k)) */
  addDoubleOption("RESIDUAL_FEM_UTOL", Res_FEM_UTOL, -9.0);
  /* DESCRIPTION: FEM: RTOL = norm(Residual(k)) / norm(Residual(0)) */
  addDoubleOption("RESIDUAL_FEM_RTOL", Res_FEM_RTOL, -9.0);
  /* DESCRIPTION: FEM: ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) */
  addDoubleOption("RESIDUAL_FEM_ETOL", Res_FEM_ETOL, -9.0);
  /*!\brief RESIDUAL_FUNC_FLOW\n DESCRIPTION: Flow functional for the Residual criteria\n OPTIONS: See \link Residual_Map \endlink \n DEFAULT: RHO_RESIDUAL \ingroup Config*/
  addEnumOption("RESIDUAL_FUNC_FLOW", Residual_Func_Flow, Residual_Map, RHO_RESIDUAL);
  /*!\brief STARTCONV_ITER\n DESCRIPTION: Iteration number to begin convergence monitoring\n DEFAULT: 5 \ingroup Config*/
  addUnsignedLongOption("STARTCONV_ITER", StartConv_Iter, 5);
  /*!\brief CAUCHY_ELEMS\n DESCRIPTION: Number of elements to apply the criteria. \n DEFAULT 100 \ingroup Config*/
  addUnsignedShortOption("CAUCHY_ELEMS", Cauchy_Elems, 100);
  /*!\brief CAUCHY_EPS\n DESCRIPTION: Epsilon to control the series convergence \n DEFAULT: 1e-10 \ingroup Config*/
  addDoubleOption("CAUCHY_EPS", Cauchy_Eps, 1E-10);
  /*!\brief CAUCHY_FUNC_FLOW
   *  \n DESCRIPTION: Flow functional for the Cauchy criteria \n OPTIONS: see \link Objective_Map \endlink \n DEFAULT: DRAG_COEFFICIENT \ingroup Config*/
  addEnumOption("CAUCHY_FUNC_FLOW", Cauchy_Func_Flow, Objective_Map, DRAG_COEFFICIENT);
  /*!\brief CAUCHY_FUNC_ADJFLOW\n DESCRIPTION: Adjoint functional for the Cauchy criteria.\n OPTIONS: See \link Sens_Map \endlink. \n DEFAULT: SENS_GEOMETRY \ingroup Config*/
  addEnumOption("CAUCHY_FUNC_ADJFLOW", Cauchy_Func_AdjFlow, Sens_Map, SENS_GEOMETRY);

  /*!\par CONFIG_CATEGORY: Multi-grid \ingroup Config*/
  /*--- Options related to Multi-grid ---*/

  /*!\brief START_UP_ITER \n DESCRIPTION: Start up iterations using the fine grid only. DEFAULT: 0 \ingroup Config*/
  addUnsignedShortOption("START_UP_ITER", nStartUpIter, 0);
  /*!\brief MGLEVEL\n DESCRIPTION: Multi-grid Levels. DEFAULT: 0 \ingroup Config*/
  addUnsignedShortOption("MGLEVEL", nMGLevels, 0);
  /*!\brief MGCYCLE\n DESCRIPTION: Multi-grid cycle. OPTIONS: See \link MG_Cycle_Map \endlink. Defualt V_CYCLE \ingroup Config*/
  addEnumOption("MGCYCLE", MGCycle, MG_Cycle_Map, V_CYCLE);
  /*!\brief MG_PRE_SMOOTH\n DESCRIPTION: Multi-grid pre-smoothing level \ingroup Config*/
  addUShortListOption("MG_PRE_SMOOTH", nMG_PreSmooth, MG_PreSmooth);
  /*!\brief MG_POST_SMOOTH\n DESCRIPTION: Multi-grid post-smoothing level \ingroup Config*/
  addUShortListOption("MG_POST_SMOOTH", nMG_PostSmooth, MG_PostSmooth);
  /*!\brief MG_CORRECTION_SMOOTH\n DESCRIPTION: Jacobi implicit smoothing of the correction \ingroup Config*/
  addUShortListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth, MG_CorrecSmooth);
  /*!\brief MG_DAMP_RESTRICTION\n DESCRIPTION: Damping factor for the residual restriction. DEFAULT: 0.75 \ingroup Config*/
  addDoubleOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.75);
  /*!\brief MG_DAMP_PROLONGATION\n DESCRIPTION: Damping factor for the correction prolongation. DEFAULT 0.75 \ingroup Config*/
  addDoubleOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.75);

  /*!\par CONFIG_CATEGORY: Spatial Discretization \ingroup Config*/
  /*--- Options related to the spatial discretization ---*/

  /*!\brief NUM_METHOD_GRAD
   *  \n DESCRIPTION: Numerical method for spatial gradients \n OPTIONS: See \link Gradient_Map \endlink. \n DEFAULT: WEIGHTED_LEAST_SQUARES. \ingroup Config*/
  addEnumOption("NUM_METHOD_GRAD", Kind_Gradient_Method, Gradient_Map, WEIGHTED_LEAST_SQUARES);
  /*!\brief LIMITER_COEFF
   *  \n DESCRIPTION: Coefficient for the limiter. DEFAULT value 0.5. Larger values decrease the extent of limiting, values approaching zero cause lower-order approximation to the solution. \ingroup Config */
  addDoubleOption("LIMITER_COEFF", LimiterCoeff, 0.1);
  /*!\brief LIMITER_COEFF_2PHASE
   *  \n DESCRIPTION: Coefficient for the limiter. DEFAULT value 0.5. Larger values decrease the extent of limiting, values approaching zero cause lower-order approximation to the solution. \ingroup Config */
  addDoubleOption("LIMITER_COEFF_2PHASE", LimiterCoeff_2phase, 0.1);
  /*!\brief LIMITER_ITER
   *  \n DESCRIPTION: Freeze the value of the limiter after a number of iterations. DEFAULT value 999999. \ingroup Config*/
  addUnsignedLongOption("LIMITER_ITER", LimiterIter, 999999);
  /*!\brief SHARP_EDGES_COEFF
   *  \n DESCRIPTION: Coefficient for detecting the limit of the sharp edges. DEFAULT value 3.0.  Use with sharp edges limiter. \ingroup Config*/
  addDoubleOption("SHARP_EDGES_COEFF", SharpEdgesCoeff, 3.0);

  /*!\brief CONV_NUM_METHOD_FLOW
   *  \n DESCRIPTION: Convective numerical method \n OPTIONS: See \link Upwind_Map \endlink , \link Centered_Map \endlink. \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centered_Flow, Kind_Upwind_Flow);
  /*!\brief SPATIAL_ORDER_FLOW
   *  \n DESCRIPTION: Spatial numerical order integration \n OPTIONS: See \link SpatialOrder_Map \endlink \n DEFAULT: SECOND_ORDER \ingroup Config*/
  addEnumOption("SPATIAL_ORDER_FLOW", SpatialOrder_Flow, SpatialOrder_Map, SECOND_ORDER);
  /*!\brief SLOPE_LIMITER_FLOW
   * DESCRIPTION: Slope limiter for the direct solution. \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_FLOW", Kind_SlopeLimit_Flow, Limiter_Map, VENKATAKRISHNAN);
  default_ad_coeff_flow[0] = 0.15; default_ad_coeff_flow[1] = 0.5; default_ad_coeff_flow[2] = 0.02;
  /*!\brief AD_COEFF_FLOW \n DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients \ingroup Config*/
  addDoubleArrayOption("AD_COEFF_FLOW", 3, Kappa_Flow, default_ad_coeff_flow);

  /*!\brief CONV_NUM_METHOD_ADJFLOW
   *  \n DESCRIPTION: Convective numerical method for the adjoint solver.
   *  \n OPTIONS:  See \link Upwind_Map \endlink , \link Centered_Map \endlink. Note: not all methods are guaranteed to be implemented for the adjoint solver. \ingroup Config */
  addConvectOption("CONV_NUM_METHOD_ADJFLOW", Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow, Kind_Upwind_AdjFlow);
  /*!\brief SPATIAL_ORDER_ADJFLOW
   *  \n DESCRIPTION: Spatial numerical order integration \n OPTIONS: See \link SpatialOrder_Map \endlink \n DEFAULT: SECOND_ORDER \ingroup Config*/
  addEnumOption("SPATIAL_ORDER_ADJFLOW", SpatialOrder_AdjFlow, SpatialOrder_Map, SECOND_ORDER);
  /*!\brief SLOPE_LIMITER_ADJFLOW
     * DESCRIPTION: Slope limiter for the adjoint solution. \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_ADJFLOW", Kind_SlopeLimit_AdjFlow, Limiter_Map, VENKATAKRISHNAN);
  default_ad_coeff_adj[0] = 0.15; default_ad_coeff_adj[1] = 0.5; default_ad_coeff_adj[2] = 0.02;
  /*!\brief AD_COEFF_ADJFLOW
   *  \n DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients for the adjoint solver.
   *  \n FORMAT and default values: AD_COEFF_ADJFLOW = (0.15, 0.5, 0.02) \ingroup Config*/
  addDoubleArrayOption("AD_COEFF_ADJFLOW", 3, Kappa_AdjFlow, default_ad_coeff_adj);

  /*!\brief SPATIAL_ORDER_TURB
   *  \n DESCRIPTION: Spatial numerical order integration.\n OPTIONS: See \link SpatialOrder_Map \endlink \n DEFAULT: FIRST_ORDER \ingroup Config*/
  addEnumOption("SPATIAL_ORDER_TURB", SpatialOrder_Turb, SpatialOrder_Map, FIRST_ORDER);
  /*!\brief SLOPE_LIMITER_TURB
   *  \n DESCRIPTION: Slope limiter  \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_TURB", Kind_SlopeLimit_Turb, Limiter_Map, VENKATAKRISHNAN);
  /*!\brief CONV_NUM_METHOD_TURB
   *  \n DESCRIPTION: Convective numerical method \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_TURB", Kind_ConvNumScheme_Turb, Kind_Centered_Turb, Kind_Upwind_Turb);
  /*!\brief SPATIAL_ORDER_2PHASE
   *  \n DESCRIPTION: Spatial numerical order integration.\n OPTIONS: See \link SpatialOrder_Map \endlink \n DEFAULT: FIRST_ORDER \ingroup Config*/
  addEnumOption("SPATIAL_ORDER_2PHASE", SpatialOrder_2phase, SpatialOrder_Map, FIRST_ORDER);
  /*!\brief SLOPE_LIMITER_2PHASE
   *  \n DESCRIPTION: Slope limiter  \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_2PHASE", Kind_SlopeLimit_2phase, Limiter_Map, VENKATAKRISHNAN);
  /*!\brief CONV_NUM_METHOD_2PHASE
   *  \n DESCRIPTION: Convective numerical method \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_2PHASE", Kind_ConvNumScheme_2phase, Kind_Centered_2phase, Kind_Upwind_2phase);
  /*!\brief SPATIAL_ORDER_ADJTURB
   *  \n DESCRIPTION: Spatial numerical order integration \n OPTIONS: See \link SpatialOrder_Map \endlink \n DEFAULT: FIRST_ORDER \ingroup Config*/
  addEnumOption("SPATIAL_ORDER_ADJTURB", SpatialOrder_AdjTurb, SpatialOrder_Map, FIRST_ORDER);
  /*!\brief SLOPE_LIMITER_ADJTURB
   *  \n DESCRIPTION: Slope limiter \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config */
  addEnumOption("SLOPE_LIMITER_ADJTURB", Kind_SlopeLimit_AdjTurb, Limiter_Map, VENKATAKRISHNAN);
  /*!\brief CONV_NUM_METHOD_ADJTURB\n DESCRIPTION: Convective numerical method for the adjoint/turbulent problem \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_ADJTURB", Kind_ConvNumScheme_AdjTurb, Kind_Centered_AdjTurb, Kind_Upwind_AdjTurb);

  /* DESCRIPTION: Viscous limiter mean flow equations */
  addBoolOption("VISCOUS_LIMITER_FLOW", Viscous_Limiter_Flow, false);
  /* DESCRIPTION: Viscous limiter turbulent equations */
  addBoolOption("VISCOUS_LIMITER_TURB", Viscous_Limiter_Turb, false);
  
  /*!\par CONFIG_CATEGORY: Adjoint and Gradient \ingroup Config*/
  /*--- Options related to the adjoint and gradient ---*/

  /*!\brief LIMIT_ADJFLOW \n DESCRIPTION: Limit value for the adjoint variable.\n DEFAULT: 1E6. \ingroup Config*/
  addDoubleOption("LIMIT_ADJFLOW", AdjointLimit, 1E6);
  /*!\brief MG_ADJFLOW\n DESCRIPTION: Multigrid with the adjoint problem. \n Defualt: YES \ingroup Config*/
  addBoolOption("MG_ADJFLOW", MG_AdjointFlow, true);

  /*!\brief OBJECTIVE_WEIGHT  \n DESCRIPTION: Adjoint problem boundary condition weights. Applies scaling factor to objective(s) \ingroup Config*/
  addDoubleListOption("OBJECTIVE_WEIGHT", nObjW, Weight_ObjFunc);
  /*!\brief OBJECTIVE_FUNCTION
   *  \n DESCRIPTION: Adjoint problem boundary condition \n OPTIONS: see \link Objective_Map \endlink \n DEFAULT: DRAG_COEFFICIENT \ingroup Config*/
  addEnumListOption("OBJECTIVE_FUNCTION", nObj, Kind_ObjFunc, Objective_Map);

  /* DESCRIPTION: parameter for the definition of a complex objective function */
  addDoubleOption("DCD_DCL_VALUE", dCD_dCL, 0.0);
  /* DESCRIPTION: parameter for the definition of a complex objective function */
  addDoubleOption("DCD_DCM_VALUE", dCD_dCM, 0.0);

  default_obj_coeff[0]=0.0; default_obj_coeff[1]=0.0; default_obj_coeff[2]=0.0;
  default_obj_coeff[3]=0.0;  default_obj_coeff[4]=0.0;
  /*!\brief OBJ_CHAIN_RULE_COEFF
  * \n DESCRIPTION: Coefficients defining the objective function gradient using the chain rule
  * with area-averaged outlet primitive variables. This is used with the genereralized outflow
  * objective.  \ingroup Config   */
  addDoubleArrayOption("OBJ_CHAIN_RULE_COEFF",5,Obj_ChainRuleCoeff,default_obj_coeff);

  default_geo_loc[0] = 0.0; default_geo_loc[1] = 1.0;
  /* DESCRIPTION: Definition of the airfoil section */
  addDoubleArrayOption("GEO_WING_BOUNDS", 2, Section_WingBounds, default_geo_loc);
  /* DESCRIPTION: Identify the axis of the section */
  addEnumOption("GEO_AXIS_STATIONS", Axis_Stations, Axis_Stations_Map, Y_AXIS);
  /* DESCRIPTION: Number of section cuts to make when calculating internal volume */
  addUnsignedShortOption("GEO_WING_STATIONS", nWingStations, 101);
  /* DESCRIPTION: Definition of the airfoil sections */
  addDoubleListOption("GEO_LOCATION_STATIONS", nLocationStations, LocationStations);
  /* DESCRIPTION: Output sectional forces for specified markers. */
  addBoolOption("GEO_PLOT_STATIONS", Plot_Section_Forces, false);
  /* DESCRIPTION: Mode of the GDC code (analysis, or gradient) */
  addEnumOption("GEO_MODE", GeometryMode, GeometryMode_Map, FUNCTION);

  /* DESCRIPTION: Drag weight in sonic boom Objective Function (from 0.0 to 1.0) */
  addDoubleOption("DRAG_IN_SONICBOOM", WeightCd, 0.0);
  /* DESCRIPTION: Sensitivity smoothing  */
  addEnumOption("SENS_SMOOTHING", Kind_SensSmooth, Sens_Smoothing_Map, NO_SMOOTH);
  /* DESCRIPTION: Adjoint frozen viscosity */
  addBoolOption("FROZEN_VISC", Frozen_Visc, true);
   /* DESCRIPTION:  */
  addDoubleOption("FIX_AZIMUTHAL_LINE", FixAzimuthalLine, 90.0);
  /*!\brief SENS_REMOVE_SHARP
   * \n DESCRIPTION: Remove sharp edges from the sensitivity evaluation  \n Format: SENS_REMOVE_SHARP = YES \n DEFAULT: NO \ingroup Config*/
  addBoolOption("SENS_REMOVE_SHARP", Sens_Remove_Sharp, false);

  /*!\par CONFIG_CATEGORY: Input/output files and formats \ingroup Config */
  /*--- Options related to input/output files and formats ---*/

  /*!\brief OUTPUT_FORMAT \n DESCRIPTION: I/O format for output plots. \n OPTIONS: see \link Output_Map \endlink \n DEFAULT: TECPLOT \ingroup Config */
  addEnumOption("OUTPUT_FORMAT", Output_FileFormat, Output_Map, TECPLOT);
  /*!\brief ACTDISK_JUMP \n DESCRIPTION: The jump is given by the difference in values or a ratio */
  addEnumOption("ACTDISK_JUMP", ActDisk_Jump, Jump_Map, DIFFERENCE);
  /*!\brief MESH_FORMAT \n DESCRIPTION: Mesh input file format \n OPTIONS: see \link Input_Map \endlink \n DEFAULT: SU2 \ingroup Config*/
  addEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, SU2);
  /* DESCRIPTION:  Mesh input file */
  addStringOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));
  /*!\brief MESH_OUT_FILENAME \n DESCRIPTION: Mesh output file name. Used when converting, scaling, or deforming a mesh. \n DEFAULT: mesh_out.su2 \ingroup Config*/
  addStringOption("MESH_OUT_FILENAME", Mesh_Out_FileName, string("mesh_out.su2"));

  /*!\brief CONV_FILENAME \n DESCRIPTION: Output file convergence history (w/o extension) \n DEFAULT: history \ingroup Config*/
  addStringOption("CONV_FILENAME", Conv_FileName, string("history"));
  /*!\brief BREAKDOWN_FILENAME \n DESCRIPTION: Output file forces breakdown \ingroup Config*/
  addStringOption("BREAKDOWN_FILENAME", Breakdown_FileName, string("forces_breakdown.dat"));
  /*!\brief CONV_FILENAME \n DESCRIPTION: Output file convergence history (w/o extension) \n DEFAULT: history \ingroup Config*/
  addStringOption("CONV_FILENAME_FSI", Conv_FileName_FSI, string("historyFSI.csv"));
  /* DESCRIPTION: Viscous limiter turbulent equations */
  addBoolOption("WRITE_CONV_FILENAME_FSI", Write_Conv_FSI, false);
  /*!\brief SOLUTION_FLOW_FILENAME \n DESCRIPTION: Restart flow input file (the file output under the filename set by RESTART_FLOW_FILENAME) \n DEFAULT: solution_flow.dat \ingroup Config */
  addStringOption("SOLUTION_FLOW_FILENAME", Solution_FlowFileName, string("solution_flow.dat"));

  addStringOption("SOLUTION_2PHASE_FILENAME", Solution_2phaseFileName, string("solution_2phase.dat"));
  /*!\brief SOLUTION_ADJ_FILENAME\n DESCRIPTION: Restart adjoint input file. Objective function abbreviation is expected. \ingroup Config*/
  addStringOption("SOLUTION_ADJ_FILENAME", Solution_AdjFileName, string("solution_adj.dat"));
  /*!\brief SOLUTION_FLOW_FILENAME \n DESCRIPTION: Restart structure input file (the file output under the filename set by RESTART_FLOW_FILENAME) \n Default: solution_flow.dat \ingroup Config */
  addStringOption("SOLUTION_STRUCTURE_FILENAME", Solution_FEMFileName, string("solution_structure.dat"));
  /*!\brief RESTART_FLOW_FILENAME \n DESCRIPTION: Output file restart flow \ingroup Config*/
  addStringOption("RESTART_FLOW_FILENAME", Restart_FlowFileName, string("restart_flow.dat"));

  addStringOption("RESTART_2PHASE_FILENAME", Restart_2phaseFileName, string("restart_2phase.dat"));
  /*!\brief RESTART_ADJ_FILENAME  \n DESCRIPTION: Output file restart adjoint. Objective function abbreviation will be appended. \ingroup Config*/
  addStringOption("RESTART_ADJ_FILENAME", Restart_AdjFileName, string("restart_adj.dat"));
  /*!\brief RESTART_WAVE_FILENAME \n DESCRIPTION: Output file restart wave \ingroup Config*/
  addStringOption("RESTART_WAVE_FILENAME", Restart_WaveFileName, string("restart_wave.dat"));
  /*!\brief RESTART_FLOW_FILENAME \n DESCRIPTION: Output file restart structure \ingroup Config*/
  addStringOption("RESTART_STRUCTURE_FILENAME", Restart_FEMFileName, string("restart_structure.dat"));
  /*!\brief VOLUME_FLOW_FILENAME  \n DESCRIPTION: Output file flow (w/o extension) variables \ingroup Config */
  addStringOption("VOLUME_FLOW_FILENAME", Flow_FileName, string("flow"));

  addStringOption("VOLUME_2PHASE_FILENAME", Two_phase_FileName, string("two_phase"));
  /*!\brief VOLUME_STRUCTURE_FILENAME
   * \n  DESCRIPTION: Output file structure (w/o extension) variables \ingroup Config*/
  addStringOption("VOLUME_STRUCTURE_FILENAME", Structure_FileName, string("structure"));
  /*!\brief SURFACE_STRUCTURE_FILENAME
   *  \n DESCRIPTION: Output file structure (w/o extension) variables \ingroup Config*/
  addStringOption("SURFACE_STRUCTURE_FILENAME", SurfStructure_FileName, string("surface_structure"));
  /*!\brief SURFACE_WAVE_FILENAME
   *  \n DESCRIPTION: Output file structure (w/o extension) variables \ingroup Config*/
  addStringOption("SURFACE_WAVE_FILENAME", SurfWave_FileName, string("surface_wave"));
  /*!\brief SURFACE_HEAT_FILENAME
   *  \n DESCRIPTION: Output file structure (w/o extension) variables \ingroup Config */
  addStringOption("SURFACE_HEAT_FILENAME", SurfHeat_FileName, string("surface_heat"));
  /*!\brief VOLUME_WAVE_FILENAME
   *  \n DESCRIPTION: Output file wave (w/o extension) variables  \ingroup Config*/
  addStringOption("VOLUME_WAVE_FILENAME", Wave_FileName, string("wave"));
  /*!\brief VOLUME_HEAT_FILENAME
   *  \n DESCRIPTION: Output file wave (w/o extension) variables  \ingroup Config*/
  addStringOption("VOLUME_HEAT_FILENAME", Heat_FileName, string("heat"));
  /*!\brief VOLUME_ADJWAVE_FILENAME
   *  \n DESCRIPTION: Output file adj. wave (w/o extension) variables  \ingroup Config*/
  addStringOption("VOLUME_ADJWAVE_FILENAME", AdjWave_FileName, string("adjoint_wave"));
  /*!\brief VOLUME_ADJ_FILENAME
   *  \n DESCRIPTION: Output file adjoint (w/o extension) variables  \ingroup Config*/
  addStringOption("VOLUME_ADJ_FILENAME", Adj_FileName, string("adjoint"));
  /*!\brief GRAD_OBJFUNC_FILENAME
   *  \n DESCRIPTION: Output objective function gradient  \ingroup Config*/
  addStringOption("GRAD_OBJFUNC_FILENAME", ObjFunc_Grad_FileName, string("of_grad.dat"));
  /*!\brief VALUE_OBJFUNC_FILENAME
   *  \n DESCRIPTION: Output objective function  \ingroup Config*/
  addStringOption("VALUE_OBJFUNC_FILENAME", ObjFunc_Value_FileName, string("of_func.dat"));
  /*!\brief SURFACE_FLOW_FILENAME
   *  \n DESCRIPTION: Output file surface flow coefficient (w/o extension)  \ingroup Config*/
  addStringOption("SURFACE_FLOW_FILENAME", SurfFlowCoeff_FileName, string("surface_flow"));
  /*!\brief SURFACE_ADJ_FILENAME
   *  \n DESCRIPTION: Output file surface adjoint coefficient (w/o extension)  \ingroup Config*/
  addStringOption("SURFACE_ADJ_FILENAME", SurfAdjCoeff_FileName, string("surface_adjoint"));
  /*!\brief SURFACE_SENS_FILENAME_FILENAME
   *  \n DESCRIPTION: Output file surface sensitivity (discrete adjoint) (w/o extension)  \ingroup Config*/
  addStringOption("SURFACE_SENS_FILENAME", SurfSens_FileName, string("surface_sens"));
  /*!\brief VOLUME_SENS_FILENAME
   *  \n DESCRIPTION: Output file volume sensitivity (discrete adjoint))  \ingroup Config*/
  addStringOption("VOLUME_SENS_FILENAME", VolSens_FileName, string("volume_sens"));
  /*!\brief WRT_SOL_FREQ
   *  \n DESCRIPTION: Writing solution file frequency  \ingroup Config*/
  addUnsignedLongOption("WRT_SOL_FREQ", Wrt_Sol_Freq, 1000);
  /*!\brief WRT_SOL_FREQ_DUALTIME
   *  \n DESCRIPTION: Writing solution file frequency for dual time  \ingroup Config*/
  addUnsignedLongOption("WRT_SOL_FREQ_DUALTIME", Wrt_Sol_Freq_DualTime, 1);
  /*!\brief WRT_CON_FREQ
   *  \n DESCRIPTION: Writing convergence history frequency  \ingroup Config*/
  addUnsignedLongOption("WRT_CON_FREQ",  Wrt_Con_Freq, 1);
  /*!\brief WRT_CON_FREQ_DUALTIME
   *  \n DESCRIPTION: Writing convergence history frequency for the dual time  \ingroup Config*/
  addUnsignedLongOption("WRT_CON_FREQ_DUALTIME",  Wrt_Con_Freq_DualTime, 10);
  /*!\brief LOW_MEMORY_OUTPUT
   *  \n DESCRIPTION: Output less information for lower memory use.  \ingroup Config*/
  addBoolOption("LOW_MEMORY_OUTPUT", Low_MemoryOutput, false);
  /*!\brief WRT_VOL_SOL
   *  \n DESCRIPTION: Write a volume solution file  \ingroup Config*/
  addBoolOption("WRT_VOL_SOL", Wrt_Vol_Sol, true);
  /*!\brief WRT_SRF_SOL
   *  \n DESCRIPTION: Write a surface solution file  \ingroup Config*/
  addBoolOption("WRT_SRF_SOL", Wrt_Srf_Sol, true);
  /*!\brief WRT_CSV_SOL
   *  \n DESCRIPTION: Write a surface CSV solution file  \ingroup Config*/
  addBoolOption("WRT_CSV_SOL", Wrt_Csv_Sol, true);
  /*!\brief WRT_RESIDUALS
   *  \n DESCRIPTION: Output residual info to solution/restart file  \ingroup Config*/
  addBoolOption("WRT_RESIDUALS", Wrt_Residuals, false);

  /*!\brief WRT_LIQUID_PROPS
   *  \n DESCRIPTION: Output residual info to solution/restart file  \ingroup Config*/
  addBoolOption("WRT_LIQUID_PROPS", Wrt_LiquidProps, false);

  /*!\brief WRT_LIMITERS
   *  \n DESCRIPTION: Output limiter value information to solution/restart file  \ingroup Config*/
  addBoolOption("WRT_LIMITERS", Wrt_Limiters, false);
  /*!\brief WRT_SHARPEDGES
   *  \n DESCRIPTION: Output sharp edge limiter information to solution/restart file  \ingroup Config*/
  addBoolOption("WRT_SHARPEDGES", Wrt_SharpEdges, false);
  /* DESCRIPTION: Output the rind layers in the solution files  \ingroup Config*/
  addBoolOption("WRT_HALO", Wrt_Halo, false);
  /*!\brief ONE_D_OUTPUT
   *  \n DESCRIPTION: Output averaged outlet flow values on specified exit marker. \n Use with MARKER_OUT_1D. \ingroup Config*/
  addBoolOption("ONE_D_OUTPUT", Wrt_1D_Output, false);
  /*!\brief CONSOLE_OUTPUT_VERBOSITY
   *  \n DESCRIPTION: Verbosity level for console output  \ingroup Config*/
  addEnumOption("CONSOLE_OUTPUT_VERBOSITY", Console_Output_Verb, Verb_Map, VERB_HIGH);


  /*!\par CONFIG_CATEGORY: Dynamic mesh definition \ingroup Config*/
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

  /*!\par CONFIG_CATEGORY: Grid adaptation \ingroup Config*/
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

  /*!\par CONFIG_CATEGORY: Aeroelastic Simulation (Typical Section Model) \ingroup Config*/
  /*--- Options related to aeroelastic simulations using the Typical Section Model) ---*/
  /* DESCRIPTION: The flutter speed index (modifies the freestream condition) */
  addDoubleOption("FLUTTER_SPEED_INDEX", FlutterSpeedIndex, 0.6);
  /* DESCRIPTION: Natural frequency of the spring in the plunging direction (rad/s). */
  addDoubleOption("PLUNGE_NATURAL_FREQUENCY", PlungeNaturalFrequency, 100);
  /* DESCRIPTION: Natural frequency of the spring in the pitching direction (rad/s). */
  addDoubleOption("PITCH_NATURAL_FREQUENCY", PitchNaturalFrequency, 100);
  /* DESCRIPTION: The airfoil mass ratio. */
  addDoubleOption("AIRFOIL_MASS_RATIO", AirfoilMassRatio, 60);
  /* DESCRIPTION: Distance in semichords by which the center of gravity lies behind the elastic axis. */
  addDoubleOption("CG_LOCATION", CG_Location, 1.8);
  /* DESCRIPTION: The radius of gyration squared (expressed in semichords) of the typical section about the elastic axis. */
  addDoubleOption("RADIUS_GYRATION_SQUARED", RadiusGyrationSquared, 3.48);
  /* DESCRIPTION: Solve the aeroelastic equations every given number of internal iterations. */
  addUnsignedShortOption("AEROELASTIC_ITER", AeroelasticIter, 3);
  
  /*!\par CONFIG_CATEGORY: Optimization Problem*/
  
  /* DESCRIPTION: Setup for design variables (upper bound) */
  addDoubleOption("OPT_BOUND_UPPER", DVBound_Upper, 1E6);
  /* DESCRIPTION: Setup for design variables (lower bound) */
  addDoubleOption("OPT_BOUND_LOWER", DVBound_Lower, -1E6);

  /*!\par CONFIG_CATEGORY: Wind Gust \ingroup Config*/
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

  /* Harmonic Balance config */
  /* DESCRIPTION: Omega_HB = 2*PI*frequency - frequencies for Harmonic Balance method */
  addDoubleListOption("OMEGA_HB", nOmega_HB, Omega_HB);

  /*!\par CONFIG_CATEGORY: Equivalent Area \ingroup Config*/
  /*--- Options related to the equivalent area ---*/

  /* DESCRIPTION: Evaluate equivalent area on the Near-Field  */
  addBoolOption("EQUIV_AREA", EquivArea, false);
  default_ea_lim[0] = 0.0; default_ea_lim[1] = 1.0; default_ea_lim[2] = 1.0;
  /* DESCRIPTION: Integration limits of the equivalent area ( xmin, xmax, Dist_NearField ) */
  addDoubleArrayOption("EA_INT_LIMIT", 3, EA_IntLimit, default_ea_lim);
  /* DESCRIPTION: Equivalent area scaling factor */
  addDoubleOption("EA_SCALE_FACTOR", EA_ScaleFactor, 1.0);

	// these options share nDV as their size in the option references; not a good idea
	/*!\par CONFIG_CATEGORY: Grid deformation \ingroup Config*/
  /*--- Options related to the grid deformation ---*/

	/* DESCRIPTION: Kind of deformation */
	addEnumListOption("DV_KIND", nDV, Design_Variable, Param_Map);
	/* DESCRIPTION: Marker of the surface to which we are going apply the shape deformation */
  addStringListOption("DV_MARKER", nMarker_DV, Marker_DV);
	/* DESCRIPTION: Parameters of the shape deformation
   - FFD_CONTROL_POINT_2D ( FFDBox ID, i_Ind, j_Ind, x_Disp, y_Disp )
   - FFD_RADIUS_2D ( FFDBox ID )
   - FFD_CAMBER_2D ( FFDBox ID, i_Ind )
   - FFD_THICKNESS_2D ( FFDBox ID, i_Ind )
   - HICKS_HENNE ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc )
   - SURFACE_BUMP ( x_start, x_end, x_Loc )
   - CST ( Lower Surface (0)/Upper Surface (1), Kulfan parameter number, Total number of Kulfan parameters for surface )
   - NACA_4DIGITS ( 1st digit, 2nd digit, 3rd and 4th digit )
   - PARABOLIC ( Center, Thickness )
   - TRANSLATION ( x_Disp, y_Disp, z_Disp )
   - ROTATION ( x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - OBSTACLE ( Center, Bump size )
   - SPHERICAL ( ControlPoint_Index, Theta_Disp, R_Disp )
   - FFD_CONTROL_POINT ( FFDBox ID, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )
   - FFD_TWIST_ANGLE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_ROTATION ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CONTROL_SURFACE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CAMBER ( FFDBox ID, i_Ind, j_Ind )
   - FFD_THICKNESS ( FFDBox ID, i_Ind, j_Ind ) */
	addDVParamOption("DV_PARAM", nDV, ParamDV, FFDTag, Design_Variable);
  /* DESCRIPTION: New value of the shape deformation */
  addDVValueOption("DV_VALUE", nDV_Value, DV_Value, nDV, ParamDV, Design_Variable);
	/* DESCRIPTION: Hold the grid fixed in a region */
  addBoolOption("HOLD_GRID_FIXED", Hold_GridFixed, false);
	default_grid_fix[0] = -1E15; default_grid_fix[1] = -1E15; default_grid_fix[2] = -1E15;
	default_grid_fix[3] =  1E15; default_grid_fix[4] =  1E15; default_grid_fix[5] =  1E15;
	/* DESCRIPTION: Coordinates of the box where the grid will be deformed (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax) */
  addDoubleArrayOption("HOLD_GRID_FIXED_COORD", 6, Hold_GridFixed_Coord, default_grid_fix);
	/* DESCRIPTION: Visualize the deformation */
  addBoolOption("VISUALIZE_DEFORMATION", Visualize_Deformation, false);
  /* DESCRIPTION: Print the residuals during mesh deformation to the console */
  addBoolOption("DEFORM_CONSOLE_OUTPUT", Deform_Output, true);
  /* DESCRIPTION: Number of nonlinear deformation iterations (surface deformation increments) */
  addUnsignedLongOption("DEFORM_NONLINEAR_ITER", GridDef_Nonlinear_Iter, 1);
  /* DESCRIPTION: Number of smoothing iterations for FEA mesh deformation */
  addUnsignedLongOption("DEFORM_LINEAR_ITER", GridDef_Linear_Iter, 1000);
  /* DESCRIPTION: Factor to multiply smallest volume for deform tolerance (0.001 default) */
  addDoubleOption("DEFORM_TOL_FACTOR", Deform_Tol_Factor, 1E-6);
  /* DESCRIPTION: Deform coefficient (-1.0 to 0.5) */
  addDoubleOption("DEFORM_COEFF", Deform_Coeff, 1E6);
  /* DESCRIPTION: Type of element stiffness imposed for FEA mesh deformation (INVERSE_VOLUME, WALL_DISTANCE, CONSTANT_STIFFNESS) */
  addEnumOption("DEFORM_STIFFNESS_TYPE", Deform_Stiffness_Type, Deform_Stiffness_Map, WALL_DISTANCE);
  /* DESCRIPTION: Poisson's ratio for constant stiffness FEA method of grid deformation*/
  addDoubleOption("DEFORM_ELASTICITY_MODULUS", Deform_ElasticityMod, 2E11);
  /* DESCRIPTION: Young's modulus and Poisson's ratio for constant stiffness FEA method of grid deformation*/
  addDoubleOption("DEFORM_POISSONS_RATIO", Deform_PoissonRatio, 0.3);
  /*  DESCRIPTION: Linear solver for the mesh deformation\n OPTIONS: see \link Linear_Solver_Map \endlink \n DEFAULT: FGMRES \ingroup Config*/
  addEnumOption("DEFORM_LINEAR_SOLVER", Kind_Deform_Linear_Solver, Linear_Solver_Map, FGMRES);
  /*  \n DESCRIPTION: Preconditioner for the Krylov linear solvers \n OPTIONS: see \link Linear_Solver_Prec_Map \endlink \n DEFAULT: LU_SGS \ingroup Config*/
  addEnumOption("DEFORM_LINEAR_SOLVER_PREC", Kind_Deform_Linear_Solver_Prec, Linear_Solver_Prec_Map, LU_SGS);

  /*!\par CONFIG_CATEGORY: Rotorcraft problem \ingroup Config*/
  /*--- option related to rotorcraft problems ---*/

  /* DESCRIPTION: MISSING ---*/
  addDoubleOption("CYCLIC_PITCH", Cyclic_Pitch, 0.0);
  /* DESCRIPTION: MISSING ---*/
  addDoubleOption("COLLECTIVE_PITCH", Collective_Pitch, 0.0);


  /*!\par CONFIG_CATEGORY: FEA solver \ingroup Config*/
  /*--- Options related to the FEA solver ---*/

  /* DESCRIPTION: Modulus of elasticity */
  addDoubleOption("ELASTICITY_MODULUS", ElasticyMod, 2E11);
  /* DESCRIPTION: Poisson ratio */
  addDoubleOption("POISSON_RATIO", PoissonRatio, 0.30);
  /* DESCRIPTION: Material density */
  addDoubleOption("MATERIAL_DENSITY", MaterialDensity, 7854);
  /*!\brief BULK_MODULUS_STRUCT \n DESCRIPTION: Value of the Bulk Modulus for a structural problem \n DEFAULT 160E9 */
  /* This is a temporal definition */
  addDoubleOption("BULK_MODULUS_STRUCT", Bulk_Modulus_Struct, 160E9);

  /*!\brief REGIME_TYPE \n  DESCRIPTION: Geometric condition \n OPTIONS: see \link Struct_Map \endlink \ingroup Config*/
  addEnumOption("GEOMETRIC_CONDITIONS", Kind_Struct_Solver, Struct_Map, SMALL_DEFORMATIONS);
  /*!\brief REGIME_TYPE \n  DESCRIPTION: Material model \n OPTIONS: see \link Material_Map \endlink \ingroup Config*/
  addEnumOption("MATERIAL_MODEL", Kind_Material, Material_Map, LINEAR_ELASTIC);
  /*!\brief REGIME_TYPE \n  DESCRIPTION: Compressibility of the material \n OPTIONS: see \link MatComp_Map \endlink \ingroup Config*/
  addEnumOption("MATERIAL_COMPRESSIBILITY", Kind_Material_Compress, MatComp_Map, COMPRESSIBLE_MAT);

  /*  DESCRIPTION: Consider a prestretch in the structural domain
  *  Options: NO, YES \ingroup Config */
  addBoolOption("PRESTRETCH", Prestretch, false);
  /*!\brief PRESTRETCH_FILENAME \n DESCRIPTION: Filename to input for prestretching membranes \n Default: prestretch_file.dat \ingroup Config */
  addStringOption("PRESTRETCH_FILENAME", Prestretch_FEMFileName, string("prestretch_file.dat"));

  /* DESCRIPTION: Iterative method for non-linear structural analysis */
  addEnumOption("NONLINEAR_FEM_SOLUTION_METHOD", Kind_SpaceIteScheme_FEA, Space_Ite_Map_FEA, NEWTON_RAPHSON);
  /* DESCRIPTION: Number of internal iterations for Newton-Raphson Method in nonlinear structural applications */
  addUnsignedLongOption("NONLINEAR_FEM_INT_ITER", Dyn_nIntIter, 10);

  /* DESCRIPTION: Formulation for bidimensional elasticity solver */
  addEnumOption("FORMULATION_ELASTICITY_2D", Kind_2DElasForm, ElasForm_2D, PLANE_STRAIN);
  /*  DESCRIPTION: Apply dead loads
  *  Options: NO, YES \ingroup Config */
  addBoolOption("DEAD_LOAD", DeadLoad, false);
  /* DESCRIPTION: Dynamic or static structural analysis */
  addEnumOption("DYNAMIC_ANALYSIS", Dynamic_Analysis, Dynamic_Map, STATIC);
  /* DESCRIPTION: Time Step for dynamic analysis (s) */
  addDoubleOption("DYN_TIMESTEP", Delta_DynTime, 0.0);
  /* DESCRIPTION: Total Physical Time for dual time stepping simulations (s) */
  addDoubleOption("DYN_TIME", Total_DynTime, 1.0);
  /* DESCRIPTION: Parameter alpha for Newmark scheme (s) */
  addDoubleOption("NEWMARK_ALPHA", Newmark_alpha, 0.25);
  /* DESCRIPTION: Parameter delta for Newmark scheme (s) */
  addDoubleOption("NEWMARK_DELTA", Newmark_delta, 0.5);
  /* DESCRIPTION: Apply the load slowly or suddenly */
  addBoolOption("SIGMOID_LOADING", Sigmoid_Load, false);
  /* DESCRIPTION: Apply the load as a ramp */
  addBoolOption("RAMP_LOADING", Ramp_Load, false);
  /* DESCRIPTION: Time while the load is to be increased linearly */
  addDoubleOption("RAMP_TIME", Ramp_Time, 1.0);
  /* DESCRIPTION: Time while the load is to be increased linearly */
  addDoubleOption("SIGMOID_TIME", Sigmoid_Time, 1.0);
  /* DESCRIPTION: Constant of steepness of the sigmoid */
  addDoubleOption("SIGMOID_K", Sigmoid_K, 10.0);

  /* DESCRIPTION: Newmark - Generalized alpha - coefficients */
  addDoubleListOption("TIME_INT_STRUCT_COEFFS", nIntCoeffs, Int_Coeffs);

  /*  DESCRIPTION: Apply dead loads. Options: NO, YES \ingroup Config */
  addBoolOption("INCREMENTAL_LOAD", IncrementalLoad, false);
  /* DESCRIPTION: Maximum number of increments of the  */
  addUnsignedLongOption("NUMBER_INCREMENTS", IncLoad_Nincrements, 10);

  default_inc_crit[0] = 0.0; default_inc_crit[1] = 0.0; default_inc_crit[2] = 0.0;
  /* DESCRIPTION: Definition of the  UTOL RTOL ETOL*/
  addDoubleArrayOption("INCREMENTAL_CRITERIA", 3, IncLoad_Criteria, default_inc_crit);

  /* DESCRIPTION: Time while the structure is static */
  addDoubleOption("STATIC_TIME", Static_Time, 0.0);

  /* DESCRIPTION: Order of the predictor */
  addUnsignedShortOption("PREDICTOR_ORDER", Pred_Order, 0);

  /* DESCRIPTION: Transfer method used for multiphysics problems */
  addEnumOption("MULTIPHYSICS_TRANSFER_METHOD", Kind_TransferMethod, Transfer_Method_Map, BROADCAST_DATA);


  /* CONFIG_CATEGORY: FSI solver */
  /*--- Options related to the FSI solver ---*/

  /*!\brief PHYSICAL_PROBLEM_FLUID_FSI
   *  DESCRIPTION: Physical governing equations \n
   *  Options: NONE (default),EULER, NAVIER_STOKES, RANS,
   *  \ingroup Config*/
  addEnumOption("FSI_FLUID_PROBLEM", Kind_Solver_Fluid_FSI, FSI_Fluid_Solver_Map, NO_SOLVER_FFSI);

  /*!\brief PHYSICAL_PROBLEM_STRUCTURAL_FSI
   *  DESCRIPTION: Physical governing equations \n
   *  Options: NONE (default), FEM_ELASTICITY
   *  \ingroup Config*/
  addEnumOption("FSI_STRUCTURAL_PROBLEM", Kind_Solver_Struc_FSI, FSI_Struc_Solver_Map, NO_SOLVER_SFSI);

  /* DESCRIPTION: Linear solver for the structural side on FSI problems */
  addEnumOption("FSI_LINEAR_SOLVER_STRUC", Kind_Linear_Solver_FSI_Struc, Linear_Solver_Map, FGMRES);
  /* DESCRIPTION: Preconditioner for the Krylov linear solvers */
  addEnumOption("FSI_LINEAR_SOLVER_PREC_STRUC", Kind_Linear_Solver_Prec_FSI_Struc, Linear_Solver_Prec_Map, LU_SGS);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("FSI_LINEAR_SOLVER_ITER_STRUC", Linear_Solver_Iter_FSI_Struc, 500);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("FSI_LINEAR_SOLVER_ERROR_STRUC", Linear_Solver_Error_FSI_Struc, 1E-6);

  /* DESCRIPTION: Restart from a steady state (sets grid velocities to 0 when loading the restart). */
  addBoolOption("RESTART_STEADY_STATE", SteadyRestart, false);

  /*  DESCRIPTION: Apply dead loads
  *  Options: NO, YES \ingroup Config */
  addBoolOption("MATCHING_MESH", MatchingMesh, true);

  /*!\par KIND_INTERPOLATION \n
   * DESCRIPTION: Type of interpolation to use for multi-zone problems. \n OPTIONS: see \link Interpolator_Map \endlink
   * Sets Kind_Interpolation \ingroup Config
   */
  addEnumOption("KIND_INTERPOLATION", Kind_Interpolation, Interpolator_Map, NEAREST_NEIGHBOR);

  /* DESCRIPTION: Maximum number of FSI iterations */
  addUnsignedShortOption("FSI_ITER", nIterFSI, 1);
  /* DESCRIPTION: Aitken's static relaxation factor */
  addDoubleOption("STAT_RELAX_PARAMETER", AitkenStatRelax, 0.4);
  /* DESCRIPTION: Aitken's dynamic maximum relaxation factor for the first iteration */
  addDoubleOption("AITKEN_DYN_MAX_INITIAL", AitkenDynMaxInit, 0.5);
  /* DESCRIPTION: Aitken's dynamic minimum relaxation factor for the first iteration */
  addDoubleOption("AITKEN_DYN_MIN_INITIAL", AitkenDynMinInit, 0.5);
  /* DESCRIPTION: Type of gust */
  addEnumOption("BGS_RELAXATION", Kind_BGS_RelaxMethod, AitkenForm_Map, NO_RELAXATION);


  /*!\par CONFIG_CATEGORY: Wave solver \ingroup Config*/
  /*--- options related to the wave solver ---*/

  /* DESCRIPTION: Constant wave speed */
  addDoubleOption("WAVE_SPEED", Wave_Speed, 331.79);

  /*!\par CONFIG_CATEGORY: Heat solver \ingroup Config*/
  /*--- options related to the heat solver ---*/

  /* DESCRIPTION: Thermal diffusivity constant */
  addDoubleOption("THERMAL_DIFFUSIVITY", Thermal_Diffusivity, 1.172E-5);

  /*!\par CONFIG_CATEGORY: Visualize Control Volumes \ingroup Config*/
  /*--- options related to visualizing control volumes ---*/

  /* DESCRIPTION: Node number for the CV to be visualized */
  addLongOption("VISUALIZE_CV", Visualize_CV, -1);

  /*!\par CONFIG_CATEGORY: Inverse design problem \ingroup Config*/
  /*--- options related to inverse design problem ---*/

  /* DESCRIPTION: Evaluate inverse design on the surface  */
  addBoolOption("INV_DESIGN_CP", InvDesign_Cp, false);

  /* DESCRIPTION: Evaluate inverse design on the surface  */
  addBoolOption("INV_DESIGN_HEATFLUX", InvDesign_HeatFlux, false);

  /*!\par CONFIG_CATEGORY: Unsupported options \ingroup Config*/
  /*--- Options that are experimental and not intended for general use ---*/

  /* DESCRIPTION: Write extra output */
  addBoolOption("EXTRA_OUTPUT", ExtraOutput, false);

  /*--- options related to the FFD problem ---*/
  /*!\par CONFIG_CATEGORY:FFD point inversion \ingroup Config*/
  
  /* DESCRIPTION: Fix I plane */
  addShortListOption("FFD_FIX_I", nFFD_Fix_IDir, FFD_Fix_IDir);
  
  /* DESCRIPTION: Fix J plane */
  addShortListOption("FFD_FIX_J", nFFD_Fix_JDir, FFD_Fix_JDir);
  
  /* DESCRIPTION: Fix K plane */
  addShortListOption("FFD_FIX_K", nFFD_Fix_KDir, FFD_Fix_KDir);
  
  /* DESCRIPTION: FFD symmetry plane (j=0) */
  addBoolOption("FFD_SYMMETRY_PLANE", FFD_Symmetry_Plane, false);

  /* DESCRIPTION: Define different coordinates systems for the FFD */
  addEnumOption("FFD_COORD_SYSTEM", FFD_CoordSystem, CoordSystem_Map, CARTESIAN);

  /* DESCRIPTION: Axis information for the spherical and cylindrical coord system */
  default_ffd_axis[0] = 0.0; default_ffd_axis[1] = 0.0; default_ffd_axis[2] =0.0;
  addDoubleArrayOption("FFD_AXIS", 3, FFD_Axis, default_ffd_axis);

  /* DESCRIPTION: Number of total iterations in the FFD point inversion */
  addUnsignedShortOption("FFD_ITERATIONS", nFFD_Iter, 500);

  /* DESCRIPTION: Free surface damping coefficient */
	addDoubleOption("FFD_TOLERANCE", FFD_Tol, 1E-10);
  
  /* DESCRIPTION: Free surface damping coefficient */
  addDoubleOption("FFD_SCALE", FFD_Scale, 1.0);

  /* DESCRIPTION: Definition of the FFD boxes */
  addFFDDefOption("FFD_DEFINITION", nFFDBox, CoordFFDBox, TagFFDBox);
  
  /* DESCRIPTION: Definition of the FFD boxes */
  addFFDDegreeOption("FFD_DEGREE", nFFDBox, DegreeFFDBox);
  
  /* DESCRIPTION: Surface continuity at the intersection with the FFD */
  addEnumOption("FFD_CONTINUITY", FFD_Continuity, Continuity_Map, DERIVATIVE_2ND);

  /* DESCRIPTION: Kind of blending for the FFD definition */
  addEnumOption("FFD_BLENDING", FFD_Blending, Blending_Map, BEZIER );

  /* DESCRIPTION: Order of the BSplines for BSpline Blending function */
  default_ad_coeff_flow[0] = 2; default_ad_coeff_flow[1] = 2; default_ad_coeff_flow[2] = 2;
  addDoubleArrayOption("FFD_BSPLINE_ORDER", 3, FFD_BSpline_Order,default_ad_coeff_flow);

  /*--- Options for the automatic differentiation methods ---*/
  /*!\par CONFIG_CATEGORY: Automatic Differentation options\ingroup Config*/

  /* DESCRIPTION: Direct differentiation mode (forward) */
  addEnumOption("DIRECT_DIFF", DirectDiff, DirectDiff_Var_Map, NO_DERIVATIVE);

  /* DESCRIPTION: Automatic differentiation mode (reverse) */
  addBoolOption("AUTO_DIFF", AD_Mode, NO);

  /*--- options that are used in the python optimization scripts. These have no effect on the c++ toolsuite ---*/
  /*!\par CONFIG_CATEGORY:Python Options\ingroup Config*/

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
  
  /*!\brief OPT_COMBINE_OBJECTIVE
   *  \n DESCRIPTION: Flag specifying whether to internally combine a multi-objective function or treat separately */
  addPythonOption("OPT_COMBINE_OBJECTIVE");

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

  /* DESCRIPTION: Activate ParMETIS mode for testing */
  addBoolOption("PARMETIS", ParMETIS, false);
  
  /* END_CONFIG_OPTIONS */

}

void CConfig::SetConfig_Parsing(char case_filename[MAX_STRING_SIZE]) {
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

  map<string, bool> included_options;

  /*--- Parse the configuration file and set the options ---*/
  
  while (getline (case_file, text_line)) {
    
    if (err_count >= max_err_count) {
      errorString.append("too many errors. Stopping parse");

      cout << errorString << endl;
      throw(1);
    }
    
    if (TokenizeString(text_line, option_name, option_value)) {
      
      /*--- See if it's a python option ---*/

      if (option_map.find(option_name) == option_map.end()) {
          string newString;
          newString.append(option_name);
          newString.append(": invalid option name");
          newString.append(". Check current SU2 options in config_template.cfg.");
          newString.append("\n");
          errorString.append(newString);
          err_count++;
        continue;
      }

      /*--- Option exists, check if the option has already been in the config file ---*/
      
      if (included_options.find(option_name) != included_options.end()) {
        string newString;
        newString.append(option_name);
        newString.append(": option appears twice");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }


      /*--- New found option. Add it to the map, and delete from all options ---*/
      
      included_options.insert(pair<string, bool>(option_name, true));
      all_options.erase(option_name);

      /*--- Set the value and check error ---*/
      
      string out = option_map[option_name]->SetValue(option_value);
      if (out.compare("") != 0) {
        errorString.append(out);
        errorString.append("\n");
        err_count++;
      }
    }
  }

  /*--- See if there were any errors parsing the config file ---*/
      
  if (errorString.size() != 0) {
    if (rank == MASTER_NODE) cout << errorString << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Set the default values for all of the options that weren't set ---*/
      
  for (map<string, bool>::iterator iter = all_options.begin(); iter != all_options.end(); ++iter) {
    option_map[iter->first]->SetDefault();
  }

  case_file.close();
  
}

bool CConfig::SetRunTime_Parsing(char case_filename[MAX_STRING_SIZE]) {
  string text_line, option_name;
  ifstream case_file;
  vector<string> option_value;
  int rank = MASTER_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Read the configuration file ---*/
  
  case_file.open(case_filename, ios::in);
  
  if (case_file.fail()) { return false; }
  
  string errorString;
  
  int err_count = 0;  // How many errors have we found in the config file
  int max_err_count = 30; // Maximum number of errors to print before stopping
  
  map<string, bool> included_options;
  
  /*--- Parse the configuration file and set the options ---*/
  
  while (getline (case_file, text_line)) {
    
    if (err_count >= max_err_count) {
      errorString.append("too many errors. Stopping parse");
      
      cout << errorString << endl;
      throw(1);
    }
    
    if (TokenizeString(text_line, option_name, option_value)) {
      
      if (option_map.find(option_name) == option_map.end()) {
        
        /*--- See if it's a python option ---*/
        
        string newString;
        newString.append(option_name);
        newString.append(": invalid option name");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }
      
      /*--- Option exists, check if the option has already been in the config file ---*/
      
      if (included_options.find(option_name) != included_options.end()) {
        string newString;
        newString.append(option_name);
        newString.append(": option appears twice");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }
      
      /*--- New found option. Add it to the map, and delete from all options ---*/
      
      included_options.insert(pair<string, bool>(option_name, true));
      all_options.erase(option_name);
      
      /*--- Set the value and check error ---*/
      
      string out = option_map[option_name]->SetValue(option_value);
      if (out.compare("") != 0) {
        errorString.append(out);
        errorString.append("\n");
        err_count++;
      }
      
    }
  }
  
  /*--- See if there were any errors parsing the runtime file ---*/
  
  if (errorString.size() != 0) {
    if (rank == MASTER_NODE) cout << errorString << endl;
    exit(EXIT_FAILURE);
  }
  
  case_file.close();
  
  return true;
  
}

void CConfig::SetPostprocessing(unsigned short val_software, unsigned short val_izone, unsigned short val_nDim) {
  
  unsigned short iZone, iCFL, iMarker;
  bool ideal_gas       = (Kind_FluidModel == STANDARD_AIR || Kind_FluidModel == IDEAL_GAS );
  bool standard_air       = (Kind_FluidModel == STANDARD_AIR);
  
#ifdef HAVE_MPI
  int size = SINGLE_NODE;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
#ifndef HAVE_TECIO
  if (Output_FileFormat == TECPLOT_BINARY) {
    cout << "Tecplot binary file requested but SU2 was built without TecIO support." << "\n";
    Output_FileFormat = TECPLOT;
  }
#endif
  
  /*--- Fixed CM mode requires a static movement of the grid ---*/
  
  if (Fixed_CM_Mode) {
    Grid_Movement= true;
  	 nGridMovement = 1;
  	 Kind_GridMovement = new unsigned short[nGridMovement];
  	 Kind_GridMovement[0] = MOVING_HTP;
  }

  /*--- Store the SU2 module that we are executing. ---*/
  
  Kind_SU2 = val_software;

  /*--- Set the default for thrust in ActDisk ---*/
  
  if ((Kind_ActDisk == NET_THRUST) || (Kind_ActDisk == BC_THRUST)
      || (Kind_ActDisk == DRAG_MINUS_THRUST) || (Kind_ActDisk == MASSFLOW)
      || (Kind_ActDisk == POWER))
    ActDisk_Jump = RATIO;

  /*--- If Kind_Obj has not been specified, these arrays need to take a default --*/

  if (Weight_ObjFunc == NULL and Kind_ObjFunc == NULL) {
    Kind_ObjFunc = new unsigned short[1];
    Kind_ObjFunc[0] = DRAG_COEFFICIENT;
    Weight_ObjFunc = new su2double[1];
    Weight_ObjFunc[0] = 1.0;
    nObj=1;
    nObjW=1;
  }
  /*-- Correct for case where Weight_ObjFunc has not been provided or has length < kind_objfunc---*/
  
  if (nObjW<nObj) {
    if (Weight_ObjFunc!= NULL) {
      cout <<"The option OBJECTIVE_WEIGHT must either have the same length as OBJECTIVE_FUNCTION,\n"<<
          "or be deleted from the config file (equal weights will be applied)."<< endl;
      exit(EXIT_FAILURE);
    }
    Weight_ObjFunc = new su2double[nObj];
    for (unsigned short iObj=0; iObj<nObj; iObj++)
      Weight_ObjFunc[iObj] = 1.0;
  }
  /*--- Ignore weights if only one objective provided ---*/
  
  if (nObj == 1 )
      Weight_ObjFunc[0] = 1.0;

  /*--- Maker sure that nMarker = nObj ---*/

  if (nObj>0) {
    if (nMarker_Monitoring!=nObj and Marker_Monitoring!= NULL) {
      if (nMarker_Monitoring==1) {
        /*-- If only one marker was listed with multiple objectives, set that marker as the marker for each objective ---*/
        nMarker_Monitoring = nObj;
        string marker = Marker_Monitoring[0];
        delete[] Marker_Monitoring;
        Marker_Monitoring = new string[nMarker_Monitoring];
        for (iMarker=0; iMarker<nMarker_Monitoring; iMarker++)
          Marker_Monitoring[iMarker] = marker;
      }
      else if(nObj>1) {
        cout <<"When using more than one OBJECTIVE_FUNCTION, MARKER_MONTIOR must be the same length or length 1. \n "<<
            "For multiple surfaces per objective, list the objective multiple times. \n"<<
            "For multiple objectives per marker either use one marker overall or list the marker multiple times."<<endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  /*--- Low memory only for ASCII Tecplot ---*/

  if (Output_FileFormat != TECPLOT) Low_MemoryOutput = NO;
  
  /*--- Deactivate the multigrid in the adjoint problem ---*/
  
  if ((ContinuousAdjoint && !MG_AdjointFlow) ||
      (Unsteady_Simulation == TIME_STEPPING)) { nMGLevels = 0; }

  /*--- If Fluid Structure Interaction, set the solver for each zone.
   *--- ZONE_0 is the zone of the fluid.
   *--- All the other zones are structure.
   *--- This will allow us to define multiple physics structural problems */

  if (Kind_Solver == FLUID_STRUCTURE_INTERACTION) {
	  if (val_izone == 0) {	Kind_Solver = Kind_Solver_Fluid_FSI; 		FSI_Problem = true;}

	  else {			 	Kind_Solver = Kind_Solver_Struc_FSI;	  	FSI_Problem = true;
	  	  	  	  	  	  	Kind_Linear_Solver = Kind_Linear_Solver_FSI_Struc;
	  	  	  	  	  	  	Kind_Linear_Solver_Prec = Kind_Linear_Solver_Prec_FSI_Struc;
	  	  	  	  	  	  	Linear_Solver_Error = Linear_Solver_Error_FSI_Struc;
	  	  	  	  	  	  	Linear_Solver_Iter = Linear_Solver_Iter_FSI_Struc;}
  }
  else { FSI_Problem = false; }

  if ((rank==MASTER_NODE) && ContinuousAdjoint && (Ref_NonDim == DIMENSIONAL) && (Kind_SU2 == SU2_CFD)) {
    cout << "WARNING: The adjoint solver should use a non-dimensional flow solution." << endl;
  }
  
  /*--- Initialize non-physical points/reconstructions to zero ---*/
  
  Nonphys_Points   = 0;
  Nonphys_Reconstr = 0;
  
  /*--- Apply a bound to the deformation if there is any design variable ---*/
  
  if (Design_Variable != NULL) {
    for (unsigned short iDV = 0; iDV < nDV; iDV++) {
      if (DV_Value != NULL)
        DV_Value[iDV][0] = max(min(DV_Value[iDV][0], DVBound_Upper), DVBound_Lower);
    }
  }
  
  if (Kind_Solver == POISSON_EQUATION) {
    Unsteady_Simulation = STEADY;
  }
  
  /*--- Set the number of external iterations to 1 for the steady state problem ---*/

  if ((Kind_Solver == HEAT_EQUATION) ||
      (Kind_Solver == WAVE_EQUATION) || (Kind_Solver == POISSON_EQUATION)) {
    nMGLevels = 0;
    if (Unsteady_Simulation == STEADY) nExtIter = 1;
    else Unst_nIntIter = 2;
  }
  
  if (Kind_Solver == FEM_ELASTICITY) {
    nMGLevels = 0;
    if (Dynamic_Analysis == STATIC)
	nExtIter = 1;
  }

  /*--- Decide whether we should be writing unsteady solution files. ---*/
  
  if (Unsteady_Simulation == STEADY ||
      Unsteady_Simulation == HARMONIC_BALANCE)
 { Wrt_Unsteady = false; }
  else { Wrt_Unsteady = true; }

  if (Kind_Solver == FEM_ELASTICITY) {

	  if (Dynamic_Analysis == STATIC) { Wrt_Dynamic = false; }
	  else { Wrt_Dynamic = true; }

  } else {
    Wrt_Dynamic = false;
  }

  /*--- Check for unsupported features. ---*/

  if ((Kind_Regime == INCOMPRESSIBLE) && (Unsteady_Simulation == HARMONIC_BALANCE)){
    cout << "Harmonic Balance not yet implemented for the incompressible solver." << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Check for Fluid model consistency ---*/

  if (standard_air) {
	if (Gamma != 1.4 || Gas_Constant != 287.058) {
		Gamma = 1.4;
		Gas_Constant = 287.058;
        }
  }
  /*--- Check for Measurement System ---*/
  
  if (SystemMeasurements == US && !standard_air) {
    cout << "Only STANDARD_AIR fluid model can be used with US Measurement System" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Check for Convective scheme available for NICFD ---*/
  
  if (!ideal_gas) {
    if (Kind_ConvNumScheme_Flow != SPACE_UPWIND) {
      cout << "Only ROE Upwind and HLLC Upwind scheme can be used for Non-Ideal Compressible Fluids" << endl;
      exit(EXIT_FAILURE);
    }
    else {
      if (Kind_Upwind_Flow != ROE && Kind_Upwind_Flow != HLLC) {
        cout << "Only ROE Upwind and HLLC Upwind scheme can be used for Non-Ideal Compressible Fluids" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  
  /*--- Check for Boundary condition available for NICFD ---*/
  
  if (!ideal_gas) {
    if (nMarker_Inlet != 0) {
      cout << "Riemann Boundary conditions or NRBC must be used for inlet and outlet with Not Ideal Compressible Fluids " << endl;
      exit(EXIT_FAILURE);
    }
    if (nMarker_Outlet != 0) {
      cout << "Riemann Boundary conditions or NRBC must be used outlet with Not Ideal Compressible Fluids " << endl;
      exit(EXIT_FAILURE);
    }
    
    if (nMarker_FarField != 0) {
      cout << "Riemann Boundary conditions or NRBC must be used outlet with Not Ideal Compressible Fluids " << endl;
      exit(EXIT_FAILURE);
    }
    
  }
  
  /*--- Check for Boundary condition available for NICF ---*/
  
  if (ideal_gas) {
    if (SystemMeasurements == US && standard_air) {
      if (Kind_ViscosityModel != SUTHERLAND) {
        cout << "Only SUTHERLAND viscosity model can be used with US Measurement  " << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (Kind_ConductivityModel != CONSTANT_PRANDTL ) {
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
       (Kind_GridMovement[ZONE_0] != ROTATING_FRAME) &&
       (Kind_GridMovement[ZONE_0] != STEADY_TRANSLATION) &&
       (Kind_GridMovement[ZONE_0] != FLUID_STRUCTURE)))
    Grid_Movement = false;
  
  if ((Kind_SU2 == SU2_CFD || Kind_SU2 == SU2_SOL) &&
      (Unsteady_Simulation == STEADY) &&
      ((Kind_GridMovement[ZONE_0] == MOVING_HTP)))
    Grid_Movement = true;

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
  
  if (Grid_Movement &&
      (Kind_GridMovement[ZONE_0] != RIGID_MOTION) &&
      (Kind_GridMovement[ZONE_0] != ROTATING_FRAME) &&
      (Kind_GridMovement[ZONE_0] != MOVING_HTP) &&
      (Kind_GridMovement[ZONE_0] != STEADY_TRANSLATION) &&
      (Kind_GridMovement[ZONE_0] != FLUID_STRUCTURE) &&
      (Kind_GridMovement[ZONE_0] != GUST) &&
      (nGridMovement != nMarker_Moving)) {
    cout << "Number of GRID_MOVEMENT_KIND must match number of MARKER_MOVING!!" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Make sure that there aren't more than one rigid motion or
   rotating frame specified in GRID_MOVEMENT_KIND. ---*/
 /* 
  if (Grid_Movement && (Kind_GridMovement[ZONE_0] == RIGID_MOTION) &&
      (nGridMovement > 2)) {
    cout << "Can not support more than 2 type of rigid motion in GRID_MOVEMENT_KIND!!" << endl;
    exit(EXIT_FAILURE);
  }
 */ 
  /*--- In case the grid movement parameters have not been declared in the
   config file, set them equal to zero for safety. Also check to make sure
   that for each option, a value has been declared for each moving marker. ---*/
  
  unsigned short nMoving;
  if (nGridMovement > nZone) nMoving = nGridMovement;
  else nMoving = nZone;

  /*--- Motion Origin: ---*/
  
  if (Motion_Origin_X == NULL) {
    Motion_Origin_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Motion_Origin_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nMotion_Origin_X != nGridMovement)) {
      cout << "Length of MOTION_ORIGIN_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Motion_Origin_Y == NULL) {
    Motion_Origin_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Motion_Origin_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nMotion_Origin_Y != nGridMovement)) {
      cout << "Length of MOTION_ORIGIN_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Motion_Origin_Z == NULL) {
    Motion_Origin_Z = new su2double[nMoving];
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
    Translation_Rate_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Translation_Rate_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nTranslation_Rate_X != nGridMovement)) {
      cout << "Length of TRANSLATION_RATE_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Translation_Rate_Y == NULL) {
    Translation_Rate_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Translation_Rate_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nTranslation_Rate_Y != nGridMovement)) {
      cout << "Length of TRANSLATION_RATE_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Translation_Rate_Z == NULL) {
    Translation_Rate_Z = new su2double[nMoving];
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
    Rotation_Rate_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Rotation_Rate_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nRotation_Rate_X != nGridMovement)) {
      cout << "Length of ROTATION_RATE_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Rotation_Rate_Y == NULL) {
    Rotation_Rate_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Rotation_Rate_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nRotation_Rate_Y != nGridMovement)) {
      cout << "Length of ROTATION_RATE_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Rotation_Rate_Z == NULL) {
    Rotation_Rate_Z = new su2double[nMoving];
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
    Pitching_Omega_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Omega_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Omega_X != nGridMovement)) {
      cout << "Length of PITCHING_OMEGA_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Pitching_Omega_Y == NULL) {
    Pitching_Omega_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Omega_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Omega_Y != nGridMovement)) {
      cout << "Length of PITCHING_OMEGA_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Pitching_Omega_Z == NULL) {
    Pitching_Omega_Z = new su2double[nMoving];
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
    Pitching_Ampl_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Ampl_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Ampl_X != nGridMovement)) {
      cout << "Length of PITCHING_AMPL_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Pitching_Ampl_Y == NULL) {
    Pitching_Ampl_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Ampl_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Ampl_Y != nGridMovement)) {
      cout << "Length of PITCHING_AMPL_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Pitching_Ampl_Z == NULL) {
    Pitching_Ampl_Z = new su2double[nMoving];
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
    Pitching_Phase_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Phase_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Phase_X != nGridMovement)) {
      cout << "Length of PITCHING_PHASE_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Pitching_Phase_Y == NULL) {
    Pitching_Phase_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Phase_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Phase_Y != nGridMovement)) {
      cout << "Length of PITCHING_PHASE_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Pitching_Phase_Z == NULL) {
    Pitching_Phase_Z = new su2double[nMoving];
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
    Plunging_Omega_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Omega_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Omega_X != nGridMovement)) {
      cout << "Length of PLUNGING_OMEGA_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Plunging_Omega_Y == NULL) {
    Plunging_Omega_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Omega_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Omega_Y != nGridMovement)) {
      cout << "Length of PLUNGING_OMEGA_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Plunging_Omega_Z == NULL) {
    Plunging_Omega_Z = new su2double[nMoving];
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
    Plunging_Ampl_X = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Ampl_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Ampl_X != nGridMovement)) {
      cout << "Length of PLUNGING_AMPL_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Plunging_Ampl_Y == NULL) {
    Plunging_Ampl_Y = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Ampl_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Ampl_Y != nGridMovement)) {
      cout << "Length of PLUNGING_AMPL_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if (Plunging_Ampl_Z == NULL) {
    Plunging_Ampl_Z = new su2double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Ampl_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Ampl_Z != nGridMovement)) {
      cout << "Length of PLUNGING_AMPL_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  /*-- Setting Harmonic Balance period from the config file */

  if (Unsteady_Simulation == HARMONIC_BALANCE) {
  	HarmonicBalance_Period = GetHarmonicBalance_Period();
  	if (HarmonicBalance_Period < 0)  {
  		cout << "Not a valid value for time period!!" << endl;
  		exit(EXIT_FAILURE);
  	}
  	/* Initialize the Harmonic balance Frequency pointer */
  	if (Omega_HB == NULL) {
  		Omega_HB = new su2double[nOmega_HB];
  		for (iZone = 0; iZone < nOmega_HB; iZone++ )
  			Omega_HB[iZone] = 0.0;
  	}else {
  		if (nOmega_HB != nTimeInstances) {
  			cout << "Length of omega_HB  must match the number TIME_INSTANCES!!" << endl;
  			exit(EXIT_FAILURE);
  		}
  	}
  }

    /*--- Use the various rigid-motion input frequencies to determine the period to be used with harmonic balance cases.
     There are THREE types of motion to consider, namely: rotation, pitching, and plunging.
     The largest period of motion is the one to be used for harmonic balance  calculations. ---*/
    
  /*if (Unsteady_Simulation == HARMONIC_BALANCE) {
      if (!(GetGrid_Movement())) {
          // No grid movement - Time period from config file //
          HarmonicBalance_Period = GetHarmonicBalance_Period();
      }
      
      else {
          unsigned short N_MOTION_TYPES = 3;
          su2double *periods;
          periods = new su2double[N_MOTION_TYPES];
          
          //--- rotation: ---//
          
          su2double Omega_mag_rot = sqrt(pow(Rotation_Rate_X[ZONE_0],2)+pow(Rotation_Rate_Y[ZONE_0],2)+pow(Rotation_Rate_Z[ZONE_0],2));
          if (Omega_mag_rot > 0)
              periods[0] = 2*PI_NUMBER/Omega_mag_rot;
          else
              periods[0] = 0.0;
          
          //--- pitching: ---//
          
          su2double Omega_mag_pitch = sqrt(pow(Pitching_Omega_X[ZONE_0],2)+pow(Pitching_Omega_Y[ZONE_0],2)+pow(Pitching_Omega_Z[ZONE_0],2));
          if (Omega_mag_pitch > 0)
              periods[1] = 2*PI_NUMBER/Omega_mag_pitch;
          else
              periods[1] = 0.0;
          
          //--- plunging: ---//
          
          su2double Omega_mag_plunge = sqrt(pow(Plunging_Omega_X[ZONE_0],2)+pow(Plunging_Omega_Y[ZONE_0],2)+pow(Plunging_Omega_Z[ZONE_0],2));
          if (Omega_mag_plunge > 0)
              periods[2] = 2*PI_NUMBER/Omega_mag_plunge;
          else
              periods[2] = 0.0;
          
          //--- determine which period is largest ---//
          
          unsigned short iVar;
          HarmonicBalance_Period = 0.0;
          for (iVar = 0; iVar < N_MOTION_TYPES; iVar++) {
              if (periods[iVar] > HarmonicBalance_Period)
                  HarmonicBalance_Period = periods[iVar];
          }
          
          delete periods;
      }
    
  }*/
  

  
    
  /*--- Initialize the RefOriginMoment Pointer ---*/
  
  RefOriginMoment = NULL;
  RefOriginMoment = new su2double[3];
  RefOriginMoment[0] = 0.0; RefOriginMoment[1] = 0.0; RefOriginMoment[2] = 0.0;
  
  /*--- In case the moment origin coordinates have not been declared in the
   config file, set them equal to zero for safety. Also check to make sure
   that for each marker, a value has been declared for the moment origin.
   Unless only one value was specified, then set this value for all the markers
   being monitored. ---*/
  
  
  if ((nRefOriginMoment_X != nRefOriginMoment_Y) || (nRefOriginMoment_X != nRefOriginMoment_Z) ) {
    cout << "ERROR: Length of REF_ORIGIN_MOMENT_X, REF_ORIGIN_MOMENT_Y and REF_ORIGIN_MOMENT_Z must be the same!!" << endl;
    exit(EXIT_FAILURE);
  }
  
  if (RefOriginMoment_X == NULL) {
    RefOriginMoment_X = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_X[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_X == 1) {
      
      su2double aux_RefOriginMoment_X = RefOriginMoment_X[0];
      delete [] RefOriginMoment_X;
      RefOriginMoment_X = new su2double[nMarker_Monitoring];
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
    RefOriginMoment_Y = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_Y[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_Y == 1) {
      
      su2double aux_RefOriginMoment_Y = RefOriginMoment_Y[0];
      delete [] RefOriginMoment_Y;
      RefOriginMoment_Y = new su2double[nMarker_Monitoring];
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
    RefOriginMoment_Z = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_Z[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_Z == 1) {
      
      su2double aux_RefOriginMoment_Z = RefOriginMoment_Z[0];
      delete [] RefOriginMoment_Z;
      RefOriginMoment_Z = new su2double[nMarker_Monitoring];
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
  
  if (Grid_Movement && (Kind_GridMovement[ZONE_0] == AEROELASTIC || Kind_GridMovement[ZONE_0] == AEROELASTIC_RIGID_MOTION)) Aeroelastic_Simulation = true;
  else Aeroelastic_Simulation = false;
  
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
    Aeroelastic_pitch = new su2double[nMarker_Monitoring];
    Aeroelastic_plunge = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ ) {
      Aeroelastic_pitch[iMarker] = 0.0;
      Aeroelastic_plunge[iMarker] = 0.0;
    }
  }

  /*--- Fluid-Structure Interaction problems ---*/

  if (FSI_Problem) {
	  Kind_GridMovement[val_izone] = FLUID_STRUCTURE;
	  Grid_Movement = true;
  }
  
  if (MGCycle == FULLMG_CYCLE) FinestMesh = nMGLevels;
  else FinestMesh = MESH_0;
  
  if ((Kind_Solver == NAVIER_STOKES) &&
      (Kind_Turb_Model != NONE))
    Kind_Solver = RANS;

  if ((Kind_Solver == TWO_PHASE_NAVIER_STOKES) &&
      (Kind_Turb_Model != NONE))
    Kind_Solver = TWO_PHASE_RANS;
    
  Kappa_1st_Flow = Kappa_Flow[0];
  Kappa_2nd_Flow = Kappa_Flow[1];
  Kappa_4th_Flow = Kappa_Flow[2];
  Kappa_1st_AdjFlow = Kappa_AdjFlow[0];
  Kappa_2nd_AdjFlow = Kappa_AdjFlow[1];
  Kappa_4th_AdjFlow = Kappa_AdjFlow[2];
  
  /*--- Make the MG_PreSmooth, MG_PostSmooth, and MG_CorrecSmooth
   arrays consistent with nMGLevels ---*/
  
  unsigned short * tmp_smooth = new unsigned short[nMGLevels+1];
  
  if ((nMG_PreSmooth != nMGLevels+1) && (nMG_PreSmooth != 0)) {
    if (nMG_PreSmooth > nMGLevels+1) {
      
      /*--- Truncate by removing unnecessary elements at the end ---*/
      
      for (unsigned int i = 0; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PreSmooth[i];
      delete [] MG_PreSmooth;
      MG_PreSmooth=NULL;
    } else {
      
      /*--- Add additional elements equal to last element ---*/
      
      for (unsigned int i = 0; i < nMG_PreSmooth; i++)
        tmp_smooth[i] = MG_PreSmooth[i];
      for (unsigned int i = nMG_PreSmooth; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PreSmooth[nMG_PreSmooth-1];
      delete [] MG_PreSmooth;
      MG_PreSmooth=NULL;
    }
    
    nMG_PreSmooth = nMGLevels+1;
    MG_PreSmooth = new unsigned short[nMG_PreSmooth];
    for (unsigned int i = 0; i < nMG_PreSmooth; i++)
      MG_PreSmooth[i] = tmp_smooth[i];
  }
  if ((nMGLevels != 0) && (nMG_PreSmooth == 0)) {
    delete [] MG_PreSmooth;
    nMG_PreSmooth = nMGLevels+1;
    MG_PreSmooth = new unsigned short[nMG_PreSmooth];
    for (unsigned int i = 0; i < nMG_PreSmooth; i++)
      MG_PreSmooth[i] = i+1;
  }
  
  if ((nMG_PostSmooth != nMGLevels+1) && (nMG_PostSmooth != 0)) {
    if (nMG_PostSmooth > nMGLevels+1) {
      
      /*--- Truncate by removing unnecessary elements at the end ---*/
      
      for (unsigned int i = 0; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PostSmooth[i];
      delete [] MG_PostSmooth;
      MG_PostSmooth=NULL;
    } else {
      
      /*--- Add additional elements equal to last element ---*/
       
      for (unsigned int i = 0; i < nMG_PostSmooth; i++)
        tmp_smooth[i] = MG_PostSmooth[i];
      for (unsigned int i = nMG_PostSmooth; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PostSmooth[nMG_PostSmooth-1];
      delete [] MG_PostSmooth;
      MG_PostSmooth=NULL;
    }
    
    nMG_PostSmooth = nMGLevels+1;
    MG_PostSmooth = new unsigned short[nMG_PostSmooth];
    for (unsigned int i = 0; i < nMG_PostSmooth; i++)
      MG_PostSmooth[i] = tmp_smooth[i];
    
  }
  
  if ((nMGLevels != 0) && (nMG_PostSmooth == 0)) {
    delete [] MG_PostSmooth;
    nMG_PostSmooth = nMGLevels+1;
    MG_PostSmooth = new unsigned short[nMG_PostSmooth];
    for (unsigned int i = 0; i < nMG_PostSmooth; i++)
      MG_PostSmooth[i] = 0;
  }
  
  if ((nMG_CorrecSmooth != nMGLevels+1) && (nMG_CorrecSmooth != 0)) {
    if (nMG_CorrecSmooth > nMGLevels+1) {
      
      /*--- Truncate by removing unnecessary elements at the end ---*/
      
      for (unsigned int i = 0; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_CorrecSmooth[i];
      delete [] MG_CorrecSmooth;
      MG_CorrecSmooth = NULL;
    } else {
      
      /*--- Add additional elements equal to last element ---*/
      
      for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
        tmp_smooth[i] = MG_CorrecSmooth[i];
      for (unsigned int i = nMG_CorrecSmooth; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_CorrecSmooth[nMG_CorrecSmooth-1];
      delete [] MG_CorrecSmooth;
      MG_CorrecSmooth = NULL;
    }
    nMG_CorrecSmooth = nMGLevels+1;
    MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
    for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
      MG_CorrecSmooth[i] = tmp_smooth[i];
  }
  
  if ((nMGLevels != 0) && (nMG_CorrecSmooth == 0)) {
    delete [] MG_CorrecSmooth;
    nMG_CorrecSmooth = nMGLevels+1;
    MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
    for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
      MG_CorrecSmooth[i] = 0;
  }
  
  /*--- Override MG Smooth parameters ---*/
  
  if (nMG_PreSmooth != 0) MG_PreSmooth[MESH_0] = 1;
  if (nMG_PostSmooth != 0) {
    MG_PostSmooth[MESH_0] = 0;
    MG_PostSmooth[nMGLevels] = 0;
  }
  if (nMG_CorrecSmooth != 0) MG_CorrecSmooth[nMGLevels] = 0;
  
  if (Restart) MGCycle = V_CYCLE;
  
  if (ContinuousAdjoint) {
    if (Kind_Solver == EULER) Kind_Solver = ADJ_EULER;
    if (Kind_Solver == NAVIER_STOKES) Kind_Solver = ADJ_NAVIER_STOKES;
    if (Kind_Solver == RANS) Kind_Solver = ADJ_RANS;
  }
  
  nCFL = nMGLevels+1;
  CFL = new su2double[nCFL];
  CFL[0] = CFLFineGrid;
  
  /*--- Evaluate when the Cl should be evaluated ---*/
  
  Iter_Fixed_CL        = SU2_TYPE::Int(nExtIter / (su2double(Update_Alpha)+5.0));
  Iter_Fixed_CM        = SU2_TYPE::Int(nExtIter / (su2double(Update_iH)+5.0));
  Iter_Fixed_NetThrust = SU2_TYPE::Int(nExtIter / (su2double(Update_BCThrust)+5.0));

  if (ContinuousAdjoint) {
    CFL[0] = CFL[0] * CFLRedCoeff_AdjFlow;
    CFL_AdaptParam[2] *= CFLRedCoeff_AdjFlow;
    CFL_AdaptParam[3] *= CFLRedCoeff_AdjFlow;
    Iter_Fixed_CL = SU2_TYPE::Int(su2double (Iter_Fixed_CL) / CFLRedCoeff_AdjFlow);
    Iter_Fixed_CM = SU2_TYPE::Int(su2double (Iter_Fixed_CM) / CFLRedCoeff_AdjFlow);
    Iter_Fixed_NetThrust = SU2_TYPE::Int(su2double (Iter_Fixed_NetThrust) / CFLRedCoeff_AdjFlow);
  }
  
  if (Iter_Fixed_CL == 0) { Iter_Fixed_CL = nExtIter+1; Update_Alpha = 0; }
  if (Iter_Fixed_CM == 0) { Iter_Fixed_CM = nExtIter+1; Update_iH = 0; }
  if (Iter_Fixed_NetThrust == 0) { Iter_Fixed_NetThrust = nExtIter+1; Update_BCThrust = 0; }

  for (iCFL = 1; iCFL < nCFL; iCFL++)
    CFL[iCFL] = CFL[iCFL-1];
  
  if (nRKStep == 0) {
    nRKStep = 1;
    RK_Alpha_Step = new su2double[1]; RK_Alpha_Step[0] = 1.0;
  }
  
  if (nIntCoeffs == 0) {
	nIntCoeffs = 2;
	Int_Coeffs = new su2double[2]; Int_Coeffs[0] = 0.25; Int_Coeffs[1] = 0.5;
  }

  if ((Kind_SU2 == SU2_CFD) && (Kind_Solver == NO_SOLVER)) {
    cout << "PHYSICAL_PROBLEM must be set in the configuration file" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Set a flag for viscous simulations ---*/
  
  Viscous = (( Kind_Solver == NAVIER_STOKES          ) ||
             ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
			 ( Kind_Solver == TWO_PHASE_NAVIER_STOKES) ||
             ( Kind_Solver == RANS                   ) ||
             ( Kind_Solver == ADJ_RANS               ) ||
			 ( Kind_Solver == TWO_PHASE_RANS         ));
  
  /*--- To avoid boundary intersections, let's add a small constant to the planes. ---*/

  Section_WingBounds[0] += EPS;
  Section_WingBounds[1] += EPS;

  for (unsigned short iSections = 0; iSections < nLocationStations; iSections++) {
    LocationStations[iSections] += EPS;
  }

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
    if ((val_nDim == 2) && (!Axisymmetric)) RefAreaCoeff = RefAreaCoeff/12.0;
    else RefAreaCoeff = RefAreaCoeff/144.0;
    Length_Reynolds = Length_Reynolds/12.0;
    RefElemLength = RefElemLength/12.0;
    Highlite_Area = Highlite_Area/144.0;
    SemiSpan = SemiSpan/12.0;

    EA_IntLimit[0] = EA_IntLimit[0]/12.0;
    EA_IntLimit[1] = EA_IntLimit[1]/12.0;
    EA_IntLimit[2] = EA_IntLimit[2]/12.0;
    
    for (unsigned short iSections = 0; iSections < nLocationStations; iSections++) {
      LocationStations[iSections] = LocationStations[iSections]/12.0;
    }

    Section_WingBounds[0] = Section_WingBounds[0]/12.0;
    Section_WingBounds[1] = Section_WingBounds[1]/12.0;
    
    SubsonicEngine_Cyl[0] = SubsonicEngine_Cyl[0]/12.0;
    SubsonicEngine_Cyl[1] = SubsonicEngine_Cyl[1]/12.0;
    SubsonicEngine_Cyl[2] = SubsonicEngine_Cyl[2]/12.0;
    SubsonicEngine_Cyl[3] = SubsonicEngine_Cyl[3]/12.0;
    SubsonicEngine_Cyl[4] = SubsonicEngine_Cyl[4]/12.0;
    SubsonicEngine_Cyl[5] = SubsonicEngine_Cyl[5]/12.0;
    SubsonicEngine_Cyl[6] = SubsonicEngine_Cyl[6]/12.0;
    
  }
  
  /*--- Check for constant lift mode. Initialize the update flag for
   the AoA with each iteration to false  ---*/
  
  if (Fixed_CL_Mode) Update_AoA = false;
  if (Fixed_CM_Mode) Update_HTPIncidence = false;

  if (DirectDiff != NO_DERIVATIVE) {
#if !defined COMPLEX_TYPE && !defined ADOLC_FORWARD_TYPE && !defined CODI_FORWARD_TYPE
      if (Kind_SU2 == SU2_CFD) {
        cout << "SU2_CFD: Config option DIRECT_DIFF= YES requires AD or complex support!" << endl;
        cout << "Please use SU2_CFD_DIRECTDIFF (configuration/compilation is done using the preconfigure.py script)." << endl;
        exit(EXIT_FAILURE);
      }
#endif
    /*--- Initialize the derivative values ---*/
    switch (DirectDiff) {
      case D_MACH:
        SU2_TYPE::SetDerivative(Mach, 1.0);
        break;
      case D_AOA:
        SU2_TYPE::SetDerivative(AoA, 1.0);
        break;
      case D_SIDESLIP:
        SU2_TYPE::SetDerivative(AoS, 1.0);
        break;
      case D_REYNOLDS:
        SU2_TYPE::SetDerivative(Reynolds, 1.0);
        break;
      case D_TURB2LAM:
       SU2_TYPE::SetDerivative(Turb2LamViscRatio_FreeStream, 1.0);
        break;
      default:
        /*--- All other cases are handled in the specific solver ---*/
        break;
      }
  }

#if defined CODI_REVERSE_TYPE
  AD_Mode = YES;
#else
  if (AD_Mode == YES) {
    cout << "AUTO_DIFF=YES requires Automatic Differentiation support." << endl;
    cout << "Please use correct executables (configuration/compilation is done using the preconfigure.py script)." << endl;
  }
#endif

  if (DiscreteAdjoint) {
#if !defined ADOLC_REVERSE_TYPE && !defined CODI_REVERSE_TYPE
    if (Kind_SU2 == SU2_CFD) {
      cout << "SU2_CFD: Config option MATH_PROBLEM= DISCRETE_ADJOINT requires AD support!" << endl;
      cout << "Please use SU2_CFD_AD (configuration/compilation is done using the preconfigure.py script)." << endl;
      exit(EXIT_FAILURE);
    }
#endif

    /*--- Disable writing of limiters if enabled ---*/
    Wrt_Limiters = false;

    if (Unsteady_Simulation) {

      Restart_Flow = false;

      if (Grid_Movement) {
        cout << "Dynamic mesh movement currently not supported for the discrete adjoint solver." << endl;
        exit(EXIT_FAILURE);
      }

      /*--- If the averaging interval is not set, we average over all time-steps ---*/

      if (Iter_Avg_Objective == 0.0) {
        Iter_Avg_Objective = nExtIter;
      }
    }

    switch(Kind_Solver) {
      case EULER:
        Kind_Solver = DISC_ADJ_EULER;
        break;
      case RANS:
        Kind_Solver = DISC_ADJ_RANS;
        Frozen_Visc = false;
        break;
      case NAVIER_STOKES:
        Kind_Solver = DISC_ADJ_NAVIER_STOKES;
        break;
      default:
        break;
    }
  }

  /*--- Check for 2nd order w/ limiting for JST and correct ---*/
  
  if ((Kind_ConvNumScheme_Flow == SPACE_CENTERED) && (Kind_Centered_Flow == JST) && (SpatialOrder_Flow == SECOND_ORDER_LIMITER))
    SpatialOrder_Flow = SECOND_ORDER;
  
  if ((Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED) && (Kind_Centered_AdjFlow == JST) && (SpatialOrder_AdjFlow == SECOND_ORDER_LIMITER))
    SpatialOrder_AdjFlow = SECOND_ORDER;


  delete [] tmp_smooth;
  
}

void CConfig::SetMarkers(unsigned short val_software) {

  unsigned short iMarker_All, iMarker_CfgFile, iMarker_Euler, iMarker_Custom,
  iMarker_FarField, iMarker_SymWall, iMarker_Pressure, iMarker_PerBound,
  iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Fluid_InterfaceBound, iMarker_Dirichlet,
  iMarker_Inlet, iMarker_Riemann, iMarker_NRBC, iMarker_Outlet, iMarker_Isothermal,
  iMarker_HeatFlux, iMarker_EngineInflow, iMarker_EngineExhaust,
  iMarker_Displacement, iMarker_Load, iMarker_FlowLoad, iMarker_Neumann, iMarker_Internal,
  iMarker_Monitoring, iMarker_Designing, iMarker_GeoEval, iMarker_Plotting, iMarker_Analyze,
  iMarker_DV, iMarker_Moving, iMarker_Supersonic_Inlet, iMarker_Supersonic_Outlet,
  iMarker_Clamped, iMarker_FSIinterface, iMarker_Load_Dir, iMarker_Load_Sine,
  iMarker_ActDiskInlet, iMarker_ActDiskOutlet, iMarker_Out_1D;

  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  if (val_software != SU2_MSH)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /*--- Compute the total number of markers in the config file ---*/
  
  nMarker_CfgFile = nMarker_Euler + nMarker_FarField + nMarker_SymWall +
  nMarker_Pressure + nMarker_PerBound + nMarker_NearFieldBound + nMarker_Fluid_InterfaceBound +
  nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Inlet + nMarker_Riemann +
  nMarker_NRBC + nMarker_Outlet + nMarker_Isothermal + nMarker_HeatFlux +
  nMarker_EngineInflow + nMarker_EngineExhaust + nMarker_Internal +
  nMarker_Supersonic_Inlet + nMarker_Supersonic_Outlet + nMarker_Displacement + nMarker_Load +
  nMarker_FlowLoad + nMarker_Custom +
  nMarker_Clamped + nMarker_Load_Sine + nMarker_Load_Dir +
  nMarker_ActDiskInlet + nMarker_ActDiskOutlet + nMarker_Out_1D;
  
  /*--- Add the possible send/receive domains ---*/

  nMarker_Max = nMarker_CfgFile + OVERHEAD*size;
  
  /*--- Basic dimensionalization of the markers (worst scenario) ---*/

  nMarker_All = nMarker_Max;

  /*--- Allocate the memory (markers in each domain) ---*/
  
  Marker_All_TagBound       = new string[nMarker_All];			    // Store the tag that correspond with each marker.
  Marker_All_SendRecv       = new short[nMarker_All];						// +#domain (send), -#domain (receive).
  Marker_All_KindBC         = new unsigned short[nMarker_All];	// Store the kind of boundary condition.
  Marker_All_Monitoring     = new unsigned short[nMarker_All];	// Store whether the boundary should be monitored.
  Marker_All_Designing      = new unsigned short[nMarker_All];  // Store whether the boundary should be designed.
  Marker_All_Plotting       = new unsigned short[nMarker_All];	// Store whether the boundary should be plotted.
  Marker_All_Analyze  = new unsigned short[nMarker_All];	// Store whether the boundary should be plotted.
  Marker_All_FSIinterface   = new unsigned short[nMarker_All];	// Store whether the boundary is in the FSI interface.
  Marker_All_GeoEval        = new unsigned short[nMarker_All];	// Store whether the boundary should be geometry evaluation.
  Marker_All_DV             = new unsigned short[nMarker_All];	// Store whether the boundary should be affected by design variables.
  Marker_All_Moving         = new unsigned short[nMarker_All];	// Store whether the boundary should be in motion.
  Marker_All_PerBound       = new short[nMarker_All];						// Store whether the boundary belongs to a periodic boundary.
  Marker_All_Out_1D         = new unsigned short[nMarker_All];  // Store whether the boundary belongs to a 1-d output boundary.

  for (iMarker_All = 0; iMarker_All < nMarker_All; iMarker_All++) {
    Marker_All_TagBound[iMarker_All]      = "SEND_RECEIVE";
    Marker_All_SendRecv[iMarker_All]      = 0;
    Marker_All_KindBC[iMarker_All]        = 0;
    Marker_All_Monitoring[iMarker_All]    = 0;
    Marker_All_GeoEval[iMarker_All]       = 0;
    Marker_All_Designing[iMarker_All]     = 0;
    Marker_All_Plotting[iMarker_All]      = 0;
    Marker_All_Analyze[iMarker_All] = 0;
    Marker_All_FSIinterface[iMarker_All]  = 0;
    Marker_All_DV[iMarker_All]            = 0;
    Marker_All_Moving[iMarker_All]        = 0;
    Marker_All_PerBound[iMarker_All]      = 0;
    Marker_All_Out_1D[iMarker_All]        = 0;
  }

  /*--- Allocate the memory (markers in the config file) ---*/

  Marker_CfgFile_TagBound      = new string[nMarker_CfgFile];
  Marker_CfgFile_KindBC        = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_Monitoring    = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_Designing     = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_Plotting      = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_Analyze       = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_GeoEval       = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_FSIinterface	 = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_DV            = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_Moving        = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_PerBound      = new unsigned short[nMarker_CfgFile];
  Marker_CfgFile_Out_1D        = new unsigned short[nMarker_CfgFile];

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile]      = "SEND_RECEIVE";
    Marker_CfgFile_KindBC[iMarker_CfgFile]        = 0;
    Marker_CfgFile_Monitoring[iMarker_CfgFile]    = 0;
    Marker_CfgFile_GeoEval[iMarker_CfgFile]       = 0;
    Marker_CfgFile_Designing[iMarker_CfgFile]     = 0;
    Marker_CfgFile_Plotting[iMarker_CfgFile]      = 0;
    Marker_CfgFile_Analyze[iMarker_CfgFile] = 0;
    Marker_CfgFile_FSIinterface[iMarker_CfgFile]  = 0;
    Marker_CfgFile_DV[iMarker_CfgFile]            = 0;
    Marker_CfgFile_Moving[iMarker_CfgFile]        = 0;
    Marker_CfgFile_PerBound[iMarker_CfgFile]      = 0;
    Marker_CfgFile_Out_1D[iMarker_CfgFile]        = 0;
  }

  /*--- Allocate memory to store surface information (Analyze BC) ---*/

  Surface_MassFlow = new su2double[nMarker_Analyze];
  Surface_DC60 = new su2double[nMarker_Analyze];
  Surface_IDC = new su2double[nMarker_Analyze];
  Surface_IDC_Mach = new su2double[nMarker_Analyze];
  Surface_IDR = new su2double[nMarker_Analyze];
  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
     Surface_MassFlow[iMarker_Analyze] = 0.0;
   	 Surface_DC60[iMarker_Analyze] = 0.0;
     Surface_IDC[iMarker_Analyze] = 0.0;
     Surface_IDC_Mach[iMarker_Analyze] = 0.0;
     Surface_IDR[iMarker_Analyze] = 0.0;
   }

  /*--- Populate the marker information in the config file (all domains) ---*/

  iMarker_CfgFile = 0;
  for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Euler[iMarker_Euler];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = EULER_WALL;
    iMarker_CfgFile++;
  }

  for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_FarField[iMarker_FarField];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = FAR_FIELD;
    iMarker_CfgFile++;
  }

  for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_SymWall[iMarker_SymWall];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SYMMETRY_PLANE;
    iMarker_CfgFile++;
  }

  for (iMarker_Pressure = 0; iMarker_Pressure < nMarker_Pressure; iMarker_Pressure++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Pressure[iMarker_Pressure];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = PRESSURE_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_PerBound[iMarker_PerBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = PERIODIC_BOUNDARY;
    Marker_CfgFile_PerBound[iMarker_CfgFile] = iMarker_PerBound + 1;
    iMarker_CfgFile++;
  }

  ActDisk_DeltaPress = new su2double[nMarker_ActDiskInlet];
  ActDisk_DeltaTemp = new su2double[nMarker_ActDiskInlet];
  ActDisk_TotalPressRatio = new su2double[nMarker_ActDiskInlet];
  ActDisk_TotalTempRatio = new su2double[nMarker_ActDiskInlet];
  ActDisk_StaticPressRatio = new su2double[nMarker_ActDiskInlet];
  ActDisk_StaticTempRatio = new su2double[nMarker_ActDiskInlet];
  ActDisk_Power = new su2double[nMarker_ActDiskInlet];
  ActDisk_MassFlow = new su2double[nMarker_ActDiskInlet];
  ActDisk_Mach = new su2double[nMarker_ActDiskInlet];
  ActDisk_Force = new su2double[nMarker_ActDiskInlet];
  ActDisk_NetThrust = new su2double[nMarker_ActDiskInlet];
  ActDisk_BCThrust = new su2double[nMarker_ActDiskInlet];
  ActDisk_BCThrust_Old = new su2double[nMarker_ActDiskInlet];
  ActDisk_GrossThrust = new su2double[nMarker_ActDiskInlet];
  ActDisk_Area = new su2double[nMarker_ActDiskInlet];
  ActDisk_ReverseMassFlow = new su2double[nMarker_ActDiskInlet];
  
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
    ActDisk_DeltaPress[iMarker_ActDiskInlet] = 0.0;
    ActDisk_DeltaTemp[iMarker_ActDiskInlet] = 0.0;
    ActDisk_TotalPressRatio[iMarker_ActDiskInlet] = 0.0;
    ActDisk_TotalTempRatio[iMarker_ActDiskInlet] = 0.0;
    ActDisk_StaticPressRatio[iMarker_ActDiskInlet] = 0.0;
    ActDisk_StaticTempRatio[iMarker_ActDiskInlet] = 0.0;
    ActDisk_Power[iMarker_ActDiskInlet] = 0.0;
    ActDisk_MassFlow[iMarker_ActDiskInlet] = 0.0;
    ActDisk_Mach[iMarker_ActDiskInlet] = 0.0;
    ActDisk_Force[iMarker_ActDiskInlet] = 0.0;
    ActDisk_NetThrust[iMarker_ActDiskInlet] = 0.0;
    ActDisk_BCThrust[iMarker_ActDiskInlet] = 0.0;
    ActDisk_BCThrust_Old[iMarker_ActDiskInlet] = 0.0;
    ActDisk_GrossThrust[iMarker_ActDiskInlet] = 0.0;
    ActDisk_Area[iMarker_ActDiskInlet] = 0.0;
    ActDisk_ReverseMassFlow[iMarker_ActDiskInlet] = 0.0;
  }
  
  
  ActDiskInlet_MassFlow = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_Temperature = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_TotalTemperature = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_Pressure = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_TotalPressure = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_RamDrag = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_Force = new su2double[nMarker_ActDiskInlet];
  ActDiskInlet_Power = new su2double[nMarker_ActDiskInlet];
  
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_ActDiskInlet[iMarker_ActDiskInlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ACTDISK_INLET;
    ActDiskInlet_MassFlow[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_Temperature[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_TotalTemperature[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_Pressure[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_TotalPressure[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_RamDrag[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_Force[iMarker_ActDiskInlet] = 0.0;
    ActDiskInlet_Power[iMarker_ActDiskInlet] = 0.0;
    iMarker_CfgFile++;
  }
  
  ActDiskOutlet_MassFlow = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_Temperature = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_TotalTemperature = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_Pressure = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_TotalPressure = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_GrossThrust = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_Force = new su2double[nMarker_ActDiskOutlet];
  ActDiskOutlet_Power = new su2double[nMarker_ActDiskOutlet];
  
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_ActDiskOutlet[iMarker_ActDiskOutlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ACTDISK_OUTLET;
    ActDiskOutlet_MassFlow[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_Temperature[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_TotalTemperature[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_Pressure[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_TotalPressure[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_GrossThrust[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_Force[iMarker_ActDiskOutlet] = 0.0;
    ActDiskOutlet_Power[iMarker_ActDiskOutlet] = 0.0;
    iMarker_CfgFile++;
  }

  for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_NearFieldBound[iMarker_NearFieldBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = NEARFIELD_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_InterfaceBound[iMarker_InterfaceBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = INTERFACE_BOUNDARY;
    iMarker_CfgFile++;
  }
  
  for (iMarker_Fluid_InterfaceBound = 0; iMarker_Fluid_InterfaceBound < nMarker_Fluid_InterfaceBound; iMarker_Fluid_InterfaceBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Fluid_InterfaceBound[iMarker_Fluid_InterfaceBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = FLUID_INTERFACE;
    iMarker_CfgFile++;
  }

  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Dirichlet[iMarker_Dirichlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = DIRICHLET;
    iMarker_CfgFile++;
  }

  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Inlet[iMarker_Inlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = INLET_FLOW;
    iMarker_CfgFile++;
  }

  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Riemann[iMarker_Riemann];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = RIEMANN_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_NRBC = 0; iMarker_NRBC < nMarker_NRBC; iMarker_NRBC++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_NRBC[iMarker_NRBC];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = NRBC_BOUNDARY;
    iMarker_CfgFile++;
  }

  Engine_Power       = new su2double[nMarker_EngineInflow];
  Engine_Mach        = new su2double[nMarker_EngineInflow];
  Engine_Force       = new su2double[nMarker_EngineInflow];
  Engine_NetThrust   = new su2double[nMarker_EngineInflow];
  Engine_GrossThrust = new su2double[nMarker_EngineInflow];
  Engine_Area        = new su2double[nMarker_EngineInflow];
  
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    Engine_Power[iMarker_EngineInflow] = 0.0;
    Engine_Mach[iMarker_EngineInflow] = 0.0;
    Engine_Force[iMarker_EngineInflow] = 0.0;
    Engine_NetThrust[iMarker_EngineInflow] = 0.0;
    Engine_GrossThrust[iMarker_EngineInflow] = 0.0;
    Engine_Area[iMarker_EngineInflow] = 0.0;
  }
  
  Inflow_Mach = new su2double[nMarker_EngineInflow];
  Inflow_Pressure = new su2double[nMarker_EngineInflow];
  Inflow_MassFlow = new su2double[nMarker_EngineInflow];
  Inflow_ReverseMassFlow = new su2double[nMarker_EngineInflow];
  Inflow_TotalPressure = new su2double[nMarker_EngineInflow];
  Inflow_Temperature = new su2double[nMarker_EngineInflow];
  Inflow_TotalTemperature = new su2double[nMarker_EngineInflow];
  Inflow_RamDrag = new su2double[nMarker_EngineInflow];
  Inflow_Force = new su2double[nMarker_EngineInflow];
  Inflow_Power = new su2double[nMarker_EngineInflow];
  
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_EngineInflow[iMarker_EngineInflow];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ENGINE_INFLOW;
    Inflow_Mach[iMarker_EngineInflow] = 0.0;
    Inflow_Pressure[iMarker_EngineInflow] = 0.0;
    Inflow_MassFlow[iMarker_EngineInflow] = 0.0;
    Inflow_ReverseMassFlow[iMarker_EngineInflow] = 0.0;
    Inflow_TotalPressure[iMarker_EngineInflow] = 0.0;
    Inflow_Temperature[iMarker_EngineInflow] = 0.0;
    Inflow_TotalTemperature[iMarker_EngineInflow] = 0.0;
    Inflow_RamDrag[iMarker_EngineInflow] = 0.0;
    Inflow_Force[iMarker_EngineInflow] = 0.0;
    Inflow_Power[iMarker_EngineInflow] = 0.0;
    iMarker_CfgFile++;
  }
  
  Exhaust_Pressure = new su2double[nMarker_EngineExhaust];
  Exhaust_Temperature = new su2double[nMarker_EngineExhaust];
  Exhaust_MassFlow = new su2double[nMarker_EngineExhaust];
  Exhaust_TotalPressure = new su2double[nMarker_EngineExhaust];
  Exhaust_TotalTemperature = new su2double[nMarker_EngineExhaust];
  Exhaust_GrossThrust = new su2double[nMarker_EngineExhaust];
  Exhaust_Force = new su2double[nMarker_EngineExhaust];
  Exhaust_Power = new su2double[nMarker_EngineExhaust];
  
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_EngineExhaust[iMarker_EngineExhaust];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ENGINE_EXHAUST;
    Exhaust_Pressure[iMarker_EngineExhaust] = 0.0;
    Exhaust_Temperature[iMarker_EngineExhaust] = 0.0;
    Exhaust_MassFlow[iMarker_EngineExhaust] = 0.0;
    Exhaust_TotalPressure[iMarker_EngineExhaust] = 0.0;
    Exhaust_TotalTemperature[iMarker_EngineExhaust] = 0.0;
    Exhaust_GrossThrust[iMarker_EngineExhaust] = 0.0;
    Exhaust_Force[iMarker_EngineExhaust] = 0.0;
    Exhaust_Power[iMarker_EngineExhaust] = 0.0;
    iMarker_CfgFile++;
  }
  
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SUPERSONIC_INLET;
    iMarker_CfgFile++;
  }
  
  for (iMarker_Supersonic_Outlet = 0; iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet; iMarker_Supersonic_Outlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Supersonic_Outlet[iMarker_Supersonic_Outlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SUPERSONIC_OUTLET;
    iMarker_CfgFile++;
  }

  for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Neumann[iMarker_Neumann];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = NEUMANN;
    iMarker_CfgFile++;
  }
  
  for (iMarker_Internal = 0; iMarker_Internal < nMarker_Internal; iMarker_Internal++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Internal[iMarker_Internal];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = INTERNAL_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Custom[iMarker_Custom];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = CUSTOM_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Outlet[iMarker_Outlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = OUTLET_FLOW;
    iMarker_CfgFile++;
  }

  for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Isothermal[iMarker_Isothermal];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ISOTHERMAL;
    iMarker_CfgFile++;
  }

  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_HeatFlux[iMarker_HeatFlux];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = HEAT_FLUX;
    iMarker_CfgFile++;
  }

  for (iMarker_Clamped = 0; iMarker_Clamped < nMarker_Clamped; iMarker_Clamped++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Clamped[iMarker_Clamped];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = CLAMPED_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Displacement[iMarker_Displacement];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = DISPLACEMENT_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Load[iMarker_Load];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = LOAD_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Load_Dir[iMarker_Load_Dir];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = LOAD_DIR_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Load_Sine = 0; iMarker_Load_Sine < nMarker_Load_Sine; iMarker_Load_Sine++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Load_Sine[iMarker_Load_Sine];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = LOAD_SINE_BOUNDARY;
    iMarker_CfgFile++;
  }


  for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_FlowLoad[iMarker_FlowLoad];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = FLOWLOAD_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Monitoring[iMarker_CfgFile] = NO;
    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Monitoring[iMarker_Monitoring])
        Marker_CfgFile_Monitoring[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_GeoEval[iMarker_CfgFile] = NO;
    for (iMarker_GeoEval = 0; iMarker_GeoEval < nMarker_GeoEval; iMarker_GeoEval++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_GeoEval[iMarker_GeoEval])
        Marker_CfgFile_GeoEval[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Designing[iMarker_CfgFile] = NO;
    for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Designing[iMarker_Designing])
        Marker_CfgFile_Designing[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Plotting[iMarker_CfgFile] = NO;
    for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Plotting[iMarker_Plotting])
        Marker_CfgFile_Plotting[iMarker_CfgFile] = YES;
  }
  
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Analyze[iMarker_CfgFile] = NO;
    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Analyze[iMarker_Analyze])
        Marker_CfgFile_Analyze[iMarker_CfgFile] = YES;
  }

  /*--- Identification of Fluid-Structure interface markers ---*/

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
	unsigned short indexMarker = 0;
    Marker_CfgFile_FSIinterface[iMarker_CfgFile] = NO;
    for (iMarker_FSIinterface = 0; iMarker_FSIinterface < nMarker_FSIinterface; iMarker_FSIinterface++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_FSIinterface[iMarker_FSIinterface])
			indexMarker = (int)(iMarker_FSIinterface/2+1);
	  Marker_CfgFile_FSIinterface[iMarker_CfgFile] = indexMarker;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_DV[iMarker_CfgFile] = NO;
    for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_DV[iMarker_DV])
        Marker_CfgFile_DV[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Moving[iMarker_CfgFile] = NO;
    for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Moving[iMarker_Moving])
        Marker_CfgFile_Moving[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Out_1D[iMarker_CfgFile] = NO;
    for (iMarker_Out_1D = 0; iMarker_Out_1D < nMarker_Out_1D; iMarker_Out_1D++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Out_1D[iMarker_Out_1D])
        Marker_CfgFile_Out_1D[iMarker_CfgFile] = YES;
  }

}

void CConfig::SetOutput(unsigned short val_software, unsigned short val_izone) {

  unsigned short iMarker_Euler, iMarker_Custom, iMarker_FarField,
  iMarker_SymWall, iMarker_PerBound, iMarker_Pressure, iMarker_NearFieldBound,
  iMarker_InterfaceBound, iMarker_Fluid_InterfaceBound, iMarker_Dirichlet, iMarker_Inlet, iMarker_Riemann,
  iMarker_NRBC, iMarker_MixBound, iMarker_Outlet, iMarker_Isothermal, iMarker_HeatFlux,
  iMarker_EngineInflow, iMarker_EngineExhaust, iMarker_Displacement,
  iMarker_Load, iMarker_FlowLoad, iMarker_Neumann, iMarker_Internal, iMarker_Monitoring,
  iMarker_Designing, iMarker_GeoEval, iMarker_Plotting, iMarker_Analyze, iMarker_DV, iDV_Value,
  iMarker_FSIinterface, iMarker_Load_Dir, iMarker_Load_Sine, iMarker_Clamped,
  iMarker_Moving, iMarker_Supersonic_Inlet, iMarker_Supersonic_Outlet, iMarker_ActDiskInlet,
  iMarker_ActDiskOutlet;
  
  
  /*--- WARNING: when compiling on Windows, ctime() is not available. Comment out
   the two lines below that use the dt variable. ---*/
  //time_t now = time(0);
  //string dt = ctime(&now); dt[24] = '.';

  cout << endl << "-------------------------------------------------------------------------" << endl;
  cout << "|    ___ _   _ ___                                                      |" << endl;
  cout << "|   / __| | | |_  )   Release 5.0.0  \"Raven\"                            |" << endl;
  cout << "|   \\__ \\ |_| |/ /                                                      |" << endl;
  switch (val_software) {
    case SU2_CFD: cout << "|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |" << endl; break;
    case SU2_DEF: cout << "|   |___/\\___//___|   Suite (Mesh Deformation Code)                     |" << endl; break;
    case SU2_DOT: cout << "|   |___/\\___//___|   Suite (Gradient Projection Code)                  |" << endl; break;
    case SU2_MSH: cout << "|   |___/\\___//___|   Suite (Mesh Adaptation Code)                      |" << endl; break;
    case SU2_GEO: cout << "|   |___/\\___//___|   Suite (Geometry Definition Code)                  |" << endl; break;
    case SU2_SOL: cout << "|   |___/\\___//___|   Suite (Solution Exporting Code)                   |" << endl; break;
  }

  cout << "|                                                                       |" << endl;
  //cout << "|   Local date and time: " << dt << "                      |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << "| SU2 Lead Dev.: Dr. Francisco Palacios, Francisco.D.Palacios@boeing.com|" << endl;
  cout << "|                Dr. Thomas D. Economon, economon@stanford.edu          |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << "| SU2 Developers:                                                       |" << endl;
  cout << "| - Prof. Juan J. Alonso's group at Stanford University.                |" << endl;
  cout << "| - Prof. Piero Colonna's group at Delft University of Technology.      |" << endl;
  cout << "| - Prof. Nicolas R. Gauger's group at Kaiserslautern U. of Technology. |" << endl;
  cout << "| - Prof. Alberto Guardone's group at Polytechnic University of Milan.  |" << endl;
  cout << "| - Prof. Rafael Palacios' group at Imperial College London.            |" << endl;
  cout << "| - Prof. Edwin van der Weide's group at the University of Twente.      |" << endl;
  cout << "| - Prof. Vincent Terrapon's group at the University of Liege.          |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << "| Copyright (C) 2012-2017 SU2, the open-source CFD code.                |" << endl;
  cout << "|                                                                       |" << endl;
  cout << "| SU2 is free software; you can redistribute it and/or                  |" << endl;
  cout << "| modify it under the terms of the GNU Lesser General Public            |" << endl;
  cout << "| License as published by the Free Software Foundation; either          |" << endl;
  cout << "| version 2.1 of the License, or (at your option) any later version.    |" << endl;
  cout << "|                                                                       |" << endl;
  cout << "| SU2 is distributed in the hope that it will be useful,                |" << endl;
  cout << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |" << endl;
  cout << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |" << endl;
  cout << "| Lesser General Public License for more details.                       |" << endl;
  cout << "|                                                                       |" << endl;
  cout << "| You should have received a copy of the GNU Lesser General Public      |" << endl;
  cout << "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;

  cout << endl <<"------------------------ Physical Case Definition -----------------------" << endl;
  if (val_software == SU2_CFD) {
	if (FSI_Problem) {
	   cout << "Fluid-Structure Interaction." << endl;
	}

  if (DiscreteAdjoint) {
     cout <<"Discrete Adjoint equations using Algorithmic Differentiation " << endl;
     cout <<"based on the physical case: ";
  }
    switch (Kind_Solver) {
      case EULER: case DISC_ADJ_EULER:
        if (Kind_Regime == COMPRESSIBLE) cout << "Compressible Euler equations." << endl;
        if (Kind_Regime == INCOMPRESSIBLE) cout << "Incompressible Euler equations." << endl;
        switch (Kind_2phase_Model) {
                  case NONE:     cout << "Single phase simulation" << endl; break;
                  case HILL_RUS:     cout << "Hill's method of moments, Rusanov flux" << endl; break;
                }
        break;
      case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES:
        if (Kind_Regime == COMPRESSIBLE) cout << "Compressible Laminar Navier-Stokes' equations." << endl;
        if (Kind_Regime == INCOMPRESSIBLE) cout << "Incompressible Laminar Navier-Stokes' equations." << endl;
        switch (Kind_2phase_Model) {
                  case NONE:     cout << "Single phase simulation" << endl; break;
                  case HILL_RUS:     cout << "Hill's method of moments, Rusanov flux" << endl; break;
                }
        break;
      case RANS: case DISC_ADJ_RANS:
        if (Kind_Regime == COMPRESSIBLE) cout << "Compressible RANS equations." << endl;
        if (Kind_Regime == INCOMPRESSIBLE) cout << "Incompressible RANS equations." << endl;
        cout << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case SA:     cout << "Spalart Allmaras" << endl; break;
          case SA_NEG: cout << "Negative Spalart Allmaras" << endl; break;
          case SST:    cout << "Menter's SST"     << endl; break;
        }
        switch (Kind_2phase_Model) {
                  case NONE:     cout << "Single phase simulation" << endl; break;
                  case HILL_RUS:     cout << "Hill's method of moments, Rusanov flux" << endl; break;
                }
        break;
      case POISSON_EQUATION: cout << "Poisson equation." << endl; break;
      case WAVE_EQUATION: cout << "Wave equation." << endl; break;
      case HEAT_EQUATION: cout << "Heat equation." << endl; break;
      case FEM_ELASTICITY:
    	  if (Kind_Struct_Solver == SMALL_DEFORMATIONS) cout << "Geometrically linear elasticity solver." << endl;
    	  if (Kind_Struct_Solver == LARGE_DEFORMATIONS) cout << "Geometrically non-linear elasticity solver." << endl;
    	  if (Kind_Material == LINEAR_ELASTIC) cout << "Linear elastic material." << endl;
    	  if (Kind_Material == NEO_HOOKEAN) {
    		  if (Kind_Material_Compress == COMPRESSIBLE_MAT) cout << "Compressible Neo-Hookean material model." << endl;
    		  if (Kind_Material_Compress == INCOMPRESSIBLE_MAT) cout << "Incompressible Neo-Hookean material model (mean dilatation method)." << endl;
    	  }
    	  break;
      case ADJ_EULER: cout << "Continuous Euler adjoint equations." << endl; break;
      case ADJ_NAVIER_STOKES:
        if (Frozen_Visc)
          cout << "Continuous Navier-Stokes adjoint equations with frozen (laminar) viscosity." << endl;
        else
          cout << "Continuous Navier-Stokes adjoint equations." << endl;
        break;
      case ADJ_RANS:
        if (Frozen_Visc)
          cout << "Continuous RANS adjoint equations with frozen (laminar and eddy) viscosity." << endl;
        else
          cout << "Continuous RANS adjoint equations." << endl;

        break;

    }

    if ((Kind_Regime == COMPRESSIBLE) && (Kind_Solver != FEM_ELASTICITY) &&
        (Kind_Solver != HEAT_EQUATION) && (Kind_Solver != WAVE_EQUATION)) {
      cout << "Mach number: " << Mach <<"."<< endl;
      cout << "Angle of attack (AoA): " << AoA <<" deg, and angle of sideslip (AoS): " << AoS <<" deg."<< endl;
      if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
          (Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS))
        cout << "Reynolds number: " << Reynolds <<". Reference length "  << Length_Reynolds << "." << endl;
      if (Fixed_CL_Mode) cout << "Fixed CL mode, target value: " << Target_CL << "." << endl;
      if (Fixed_CM_Mode) {
      		cout << "Fixed CM mode, target value:  " << Target_CM << "." << endl;
        cout << "HTP rotation axis (X,Z): ("<< HTP_Axis[0] <<", "<< HTP_Axis[1] <<")."<< endl;
      }
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
        case MOVING_HTP:      cout << "HTP moving." << endl; break;
        case ROTATING_FRAME:  cout << "rotating reference frame." << endl; break;
        case AEROELASTIC:     cout << "aeroelastic motion." << endl; break;
        case FLUID_STRUCTURE: cout << "fluid-structure motion." << endl; break;
        case EXTERNAL:        cout << "externally prescribed motion." << endl; break;
        case AEROELASTIC_RIGID_MOTION:  cout << "rigid mesh motion plus aeroelastic motion." << endl; break;
      }
    }

    if (Restart) {
      if (Read_Binary_Restart) cout << "Reading and writing binary SU2 native restart files." << endl;
      else cout << "Reading and writing ASCII SU2 native restart files." << endl;
      if (!ContinuousAdjoint && Kind_Solver != FEM_ELASTICITY) cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
      if (ContinuousAdjoint) cout << "Read adjoint solution from: " << Solution_AdjFileName << "." << endl;
      if (Kind_Solver == FEM_ELASTICITY) cout << "Read structural solution from: " << Solution_FEMFileName << "." << endl;
    }
    else {
      cout << "No restart solution, use the values at infinity (freestream)." << endl;
    }

    if (ContinuousAdjoint)
      cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;

    
    if (Ref_NonDim == DIMENSIONAL) { cout << "Dimensional simulation." << endl; }
    else if (Ref_NonDim == FREESTREAM_PRESS_EQ_ONE) { cout << "Non-Dimensional simulation (P=1.0, Rho=1.0, T=1.0 at the farfield)." << endl; }
    else if (Ref_NonDim == FREESTREAM_VEL_EQ_MACH) { cout << "Non-Dimensional simulation (V=Mach, Rho=1.0, T=1.0 at the farfield)." << endl; }
    else if (Ref_NonDim == FREESTREAM_VEL_EQ_ONE) { cout << "Non-Dimensional simulation (V=1.0, Rho=1.0, T=1.0 at the farfield)." << endl; }
    
    if (RefAreaCoeff == 0) cout << "The reference length/area will be computed using y(2D) or z(3D) projection." << endl;
    else cout << "The reference length/area (force coefficient) is " << RefAreaCoeff << "." << endl;
    cout << "The reference length (moment computation) is " << RefLengthMoment << "." << endl;

    if ((nRefOriginMoment_X > 1) || (nRefOriginMoment_Y > 1) || (nRefOriginMoment_Z > 1)) {
      cout << "Surface(s) where the force coefficients are evaluated and \n";
      cout << "their reference origin for moment computation: \n";

      for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
        cout << "   - " << Marker_Monitoring[iMarker_Monitoring] << " (" << RefOriginMoment_X[iMarker_Monitoring] <<", "<<RefOriginMoment_Y[iMarker_Monitoring] <<", "<< RefOriginMoment_Z[iMarker_Monitoring] << ")";
        if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ".\n";
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
    
    if (nMarker_Plotting != 0) {
      cout << "Surface(s) plotted in the output file: ";
      for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++) {
        cout << Marker_Plotting[iMarker_Plotting];
        if (iMarker_Plotting < nMarker_Plotting-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }
    
    if (nMarker_Analyze != 0) {
      cout << "Surface(s) to be analyzed in detail: ";
      for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
        cout << Marker_Analyze[iMarker_Analyze];
        if (iMarker_Analyze < nMarker_Analyze-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }
    
    if (nMarker_FSIinterface != 0) {
      cout << "Surface(s) belonging to the Fluid-Structure Interaction problem: ";
      for (iMarker_FSIinterface = 0; iMarker_FSIinterface < nMarker_FSIinterface; iMarker_FSIinterface++) {
        cout << Marker_FSIinterface[iMarker_FSIinterface];
        if (iMarker_FSIinterface < nMarker_FSIinterface-1) cout << ", ";
        else cout <<".";
      }
      cout<<endl;
    }
    
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
    if (DiscreteAdjoint) {
      cout << "Input sensitivity file name: " << GetObjFunc_Extension(Solution_AdjFileName) << "." << endl;
    }else {
		cout << "Input sensitivity file name: " << SurfAdjCoeff_FileName << "." << endl;
	}
  }

	if (val_software == SU2_MSH) {
		switch (Kind_Adaptation) {
		case FULL: case WAKE: case FULL_FLOW: case FULL_ADJOINT: case SMOOTHING: case SUPERSONIC_SHOCK:
			break;
		case GRAD_FLOW:
			cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
			break;
		case GRAD_ADJOINT:
			cout << "Read adjoint flow solution from: " << Solution_AdjFileName << "." << endl;
			break;
		case GRAD_FLOW_ADJ: case COMPUTABLE: case REMAINING:
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

  if (((val_software == SU2_DEF) || (val_software == SU2_DOT)) && (Design_Variable[0] != NONE)) {

    for (unsigned short iDV = 0; iDV < nDV; iDV++) {

      
      if ((Design_Variable[iDV] != NO_DEFORMATION) &&
          (Design_Variable[iDV] != FFD_SETTING) &&
          (Design_Variable[iDV] != SURFACE_FILE)) {
        
        if (iDV == 0)
          cout << "Design variables definition (markers <-> value <-> param):" << endl;
        
        switch (Design_Variable[iDV]) {
          case FFD_CONTROL_POINT_2D:  cout << "FFD 2D (control point) <-> "; break;
          case FFD_CAMBER_2D:         cout << "FFD 2D (camber) <-> "; break;
          case FFD_THICKNESS_2D:      cout << "FFD 2D (thickness) <-> "; break;
          case FFD_TWIST_2D:          cout << "FFD 2D (twist) <-> "; break;
          case HICKS_HENNE:           cout << "Hicks Henne <-> " ; break;
          case SURFACE_BUMP:          cout << "Surface bump <-> " ; break;
          case ANGLE_OF_ATTACK:       cout << "Angle of attack <-> " ; break;
	        case CST:           	      cout << "Kulfan parameter number (CST) <-> " ; break;
          case TRANSLATION:           cout << "Translation design variable."; break;
          case SCALE:                 cout << "Scale design variable."; break;
          case NACA_4DIGITS:          cout << "NACA four digits <-> "; break;
          case PARABOLIC:             cout << "Parabolic <-> "; break;
          case AIRFOIL:               cout << "Airfoil <-> "; break;
          case ROTATION:              cout << "Rotation <-> "; break;
          case FFD_CONTROL_POINT:     cout << "FFD (control point) <-> "; break;
          case FFD_NACELLE:           cout << "FFD (nacelle) <-> "; break;
          case FFD_GULL:              cout << "FFD (gull) <-> "; break;
          case FFD_TWIST:             cout << "FFD (twist) <-> "; break;
          case FFD_ROTATION:          cout << "FFD (rotation) <-> "; break;
          case FFD_CONTROL_SURFACE:   cout << "FFD (control surface) <-> "; break;
          case FFD_CAMBER:            cout << "FFD (camber) <-> "; break;
          case FFD_THICKNESS:         cout << "FFD (thickness) <-> "; break;
          case FFD_ANGLE_OF_ATTACK:   cout << "FFD (angle of attack) <-> "; break;
          case CUSTOM:                cout << "Custom DV <-> "; break;
        }
        
        for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
          cout << Marker_DV[iMarker_DV];
          if (iMarker_DV < nMarker_DV-1) cout << ", ";
          else cout << " <-> ";
        }

        for (iDV_Value = 0; iDV_Value < nDV_Value[iDV]; iDV_Value++) {
          cout << DV_Value[iDV][iDV_Value];
          if (iDV_Value != nDV_Value[iDV]-1) cout << ", ";
        }
        cout << " <-> ";

        if ((Design_Variable[iDV] == NO_DEFORMATION) ||
            (Design_Variable[iDV] == FFD_SETTING) ||
            (Design_Variable[iDV] == SCALE) ) nParamDV = 0;
        if (Design_Variable[iDV] == ANGLE_OF_ATTACK) nParamDV = 1;
        if ((Design_Variable[iDV] == FFD_CAMBER_2D) ||
            (Design_Variable[iDV] == FFD_THICKNESS_2D) ||
            (Design_Variable[iDV] == HICKS_HENNE) ||
            (Design_Variable[iDV] == PARABOLIC) ||
            (Design_Variable[iDV] == AIRFOIL) ||
            (Design_Variable[iDV] == FFD_GULL) ||
            (Design_Variable[iDV] == FFD_ANGLE_OF_ATTACK) ) nParamDV = 2;
        if ((Design_Variable[iDV] ==  TRANSLATION) ||
            (Design_Variable[iDV] ==  NACA_4DIGITS) ||
            (Design_Variable[iDV] ==  CST) ||
            (Design_Variable[iDV] ==  SURFACE_BUMP) ||
            (Design_Variable[iDV] ==  FFD_CAMBER) ||
            (Design_Variable[iDV] ==  FFD_TWIST_2D) ||
            (Design_Variable[iDV] ==  FFD_THICKNESS) ) nParamDV = 3;
        if (Design_Variable[iDV] == FFD_CONTROL_POINT_2D) nParamDV = 5;
        if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
        if ((Design_Variable[iDV] ==  FFD_CONTROL_POINT) ||
            (Design_Variable[iDV] ==  FFD_ROTATION) ||
            (Design_Variable[iDV] ==  FFD_CONTROL_SURFACE) ) nParamDV = 7;
        if (Design_Variable[iDV] ==  CUSTOM) nParamDV = 1;
        if (Design_Variable[iDV] == FFD_TWIST) nParamDV = 8;

        for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {

          if (iParamDV == 0) cout << "( ";

          if ((iParamDV == 0) &&
              ((Design_Variable[iDV] == NO_DEFORMATION) ||
               (Design_Variable[iDV] == FFD_SETTING) ||
               (Design_Variable[iDV] == FFD_ANGLE_OF_ATTACK) ||
               (Design_Variable[iDV] == FFD_CONTROL_POINT_2D) ||
               (Design_Variable[iDV] == FFD_CAMBER_2D) ||
               (Design_Variable[iDV] == FFD_THICKNESS_2D) ||
               (Design_Variable[iDV] == FFD_TWIST_2D) ||
               (Design_Variable[iDV] == FFD_CONTROL_POINT) ||
               (Design_Variable[iDV] == FFD_NACELLE) ||
               (Design_Variable[iDV] == FFD_GULL) ||
               (Design_Variable[iDV] == FFD_TWIST) ||
               (Design_Variable[iDV] == FFD_ROTATION) ||
               (Design_Variable[iDV] == FFD_CONTROL_SURFACE) ||
               (Design_Variable[iDV] == FFD_CAMBER) ||
               (Design_Variable[iDV] == FFD_THICKNESS))) cout << FFDTag[iDV];
          else cout << ParamDV[iDV][iParamDV];

          if (iParamDV < nParamDV-1) cout << ", ";
          else cout <<" )"<< endl;
          
        }

      }
      
      else if (Design_Variable[iDV] == FFD_SETTING) {
        
        cout << "Setting the FFD box structure." << endl;
        cout << "FFD boxes definition (FFD tag <-> degree <-> coord):" << endl;
        
        for (unsigned short iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {
          
          cout << TagFFDBox[iFFDBox] << " <-> ";
          
          for (unsigned short iDegreeFFD = 0; iDegreeFFD < 3; iDegreeFFD++) {
            if (iDegreeFFD == 0) cout << "( ";
            cout << DegreeFFDBox[iFFDBox][iDegreeFFD];
            if (iDegreeFFD < 2) cout << ", ";
            else cout <<" )";
          }
          
          cout << " <-> ";

          for (unsigned short iCoordFFD = 0; iCoordFFD < 24; iCoordFFD++) {
            if (iCoordFFD == 0) cout << "( ";
            cout << CoordFFDBox[iFFDBox][iCoordFFD];
            if (iCoordFFD < 23) cout << ", ";
            else cout <<" )"<< endl;
          }
          
        }
        
      }
      
      else cout << endl;

		}
	}

	if (((val_software == SU2_CFD) && ( ContinuousAdjoint || DiscreteAdjoint)) || (val_software == SU2_DOT)) {

		cout << endl <<"----------------------- Design problem definition -----------------------" << endl;
		if (nObj==1) {
      switch (Kind_ObjFunc[0]) {
        case DRAG_COEFFICIENT:
          cout << "CD objective function." << endl;
          if (Fixed_CL_Mode) cout << "dCD/dCL = " << dCD_dCL << "." << endl;
          if (Fixed_CM_Mode) cout << "dCD/dCM = " << dCD_dCM << "." << endl;
          break;
        case LIFT_COEFFICIENT:        cout << "CL objective function." << endl; break;
        case MOMENT_X_COEFFICIENT:    cout << "CMx objective function." << endl; break;
        case MOMENT_Y_COEFFICIENT:    cout << "CMy objective function." << endl; break;
        case MOMENT_Z_COEFFICIENT:    cout << "CMz objective function." << endl; break;
        case INVERSE_DESIGN_PRESSURE: cout << "Inverse design (Cp) objective function." << endl; break;
        case INVERSE_DESIGN_HEATFLUX: cout << "Inverse design (Heat Flux) objective function." << endl; break;
        case SIDEFORCE_COEFFICIENT:   cout << "Side force objective function." << endl; break;
        case EFFICIENCY:              cout << "CL/CD objective function." << endl; break;
        case EQUIVALENT_AREA:         cout << "Equivalent area objective function. CD weight: " << WeightCd <<"."<< endl;  break;
        case NEARFIELD_PRESSURE:      cout << "Nearfield pressure objective function. CD weight: " << WeightCd <<"."<< endl;  break;
        case FORCE_X_COEFFICIENT:     cout << "X-force objective function." << endl; break;
        case FORCE_Y_COEFFICIENT:     cout << "Y-force objective function." << endl; break;
        case FORCE_Z_COEFFICIENT:     cout << "Z-force objective function." << endl; break;
        case THRUST_COEFFICIENT:      cout << "Thrust objective function." << endl; break;
        case TORQUE_COEFFICIENT:      cout << "Torque efficiency objective function." << endl; break;
        case TOTAL_HEATFLUX:          cout << "Total heat flux objective function." << endl; break;
        case MAXIMUM_HEATFLUX:        cout << "Maximum heat flux objective function." << endl; break;
        case FIGURE_OF_MERIT:         cout << "Rotor Figure of Merit objective function." << endl; break;
        case AVG_TOTAL_PRESSURE:      cout << "Average total objective pressure." << endl; break;
        case AVG_OUTLET_PRESSURE:     cout << "Average static objective pressure." << endl; break;
        case MASS_FLOW_RATE:          cout << "Mass flow rate objective function." << endl; break;
        case OUTFLOW_GENERALIZED:     cout << "Generalized outflow objective function." << endl; break;
        case AERO_DRAG_COEFFICIENT:   cout << "Aero CD objective function." << endl; break;
        case RADIAL_DISTORTION:       cout << "Radial distortion objective function." << endl; break;
        case CIRCUMFERENTIAL_DISTORTION:   cout << "Circumferential distortion objective function." << endl; break;

      }
		}
		else {
		  cout << "Weighted sum objective function." << endl;
		}

	}

	if (val_software == SU2_CFD) {
		cout << endl <<"---------------------- Space Numerical Integration ----------------------" << endl;

		if (SmoothNumGrid) cout << "There are some smoothing iterations on the grid coordinates." << endl;

    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
         (Kind_Solver == DISC_ADJ_EULER) || (Kind_Solver == DISC_ADJ_NAVIER_STOKES) || (Kind_Solver == DISC_ADJ_RANS) ) {

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
				if (Kind_Upwind_Flow == ROE) cout << "Roe (with entropy fix) solver for the flow inviscid terms."<< endl;
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

    if ((Kind_Solver == RANS) || (Kind_Solver == DISC_ADJ_RANS)) {
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


    if ((Kind_2phase_Model != NONE)) {
      if (Kind_ConvNumScheme_2phase == SPACE_UPWIND) {
        if (Kind_Upwind_2phase == SCALAR_UPWIND) cout << "Scalar upwind solver (first order) for the 2phase model."<< endl;
        switch (SpatialOrder_2phase) {
          case FIRST_ORDER: cout << "First order integration." << endl; break;
          case SECOND_ORDER: cout << "Second order integration." << endl; break;
          case SECOND_ORDER_LIMITER: cout << "Second order integration with slope limiter." << endl;
            switch (Kind_SlopeLimit_2phase) {
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
        if (Kind_Upwind_AdjFlow == ROE) cout << "Roe (with entropy fix) solver for the adjoint inviscid terms."<< endl;
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
              case SOLID_WALL_DISTANCE:
                cout << "Wall distance slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                cout << "The reference wall distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;
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
              case SOLID_WALL_DISTANCE:
                cout << "Wall distance slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                cout << "The reference wall distance is: " << SharpEdgesCoeff*RefElemLength*LimiterCoeff <<". "<< endl;
                break;
              case BARTH_JESPERSEN:
                cout << "Barth-Jespersen slope-limiting method." << endl;
                break;
            }
            break;
        }
      }
    }

    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) || (Kind_Solver == DISC_ADJ_RANS)) {
        cout << "Average of gradients with correction (viscous flow terms)." << endl;
    }

    if ((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
      cout << "Average of gradients with correction (viscous adjoint terms)." << endl;
    }

    if ((Kind_Solver == RANS) || (Kind_Solver == DISC_ADJ_RANS)) {
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

    if (Kind_Regime == INCOMPRESSIBLE) {
      cout << "Artificial compressibility factor: " << ArtComp_Factor << "." << endl;
    }

    cout << endl <<"---------------------- Time Numerical Integration -----------------------" << endl;

    if (Kind_Solver != FEM_ELASTICITY) {
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
    }
	else {
		switch (Dynamic_Analysis) {
		  case NO:
			cout << "Static structural analysis." << endl; break;
		  case YES:
			cout << "Dynamic structural analysis."<< endl;
			cout << "Time step provided by the user for the dynamic analysis(s): "<< Delta_DynTime << "." << endl;
			break;
		}
	}

    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
        (Kind_Solver == DISC_ADJ_EULER) || (Kind_Solver == DISC_ADJ_NAVIER_STOKES) || (Kind_Solver == DISC_ADJ_RANS)) {
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
              break;
            case FGMRES || RESTARTED_FGMRES:
              cout << "FGMRES is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              break;
            case SMOOTHER_JACOBI:
              cout << "A Jacobi method is used for smoothing the linear system." << endl;
              break;
            case SMOOTHER_ILU:
              cout << "A ILU0 method is used for smoothing the linear system." << endl;
              break;
            case SMOOTHER_LUSGS:
              cout << "A LU-SGS method is used for smoothing the linear system." << endl;
              break;
            case SMOOTHER_LINELET:
              cout << "A Linelet method is used for smoothing the linear system." << endl;
              break;
          }
          break;
        case CLASSICAL_RK4_EXPLICIT:
          cout << "Classical RK4 explicit method for the flow equations." << endl;
          cout << "Number of steps: " << 4 << endl;
          cout << "Time coefficients: {0.5, 0.5, 1, 1}" << endl;
          cout << "Function coefficients: {1/6, 1/3, 1/3, 1/6}" << endl;
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

    if (nMGLevels !=0) {
      
      if (nStartUpIter != 0) cout << "A total of " << nStartUpIter << " start up iterations on the fine grid."<< endl;
      if (MGCycle == V_CYCLE) cout << "V Multigrid Cycle, with " << nMGLevels << " multigrid levels."<< endl;
      if (MGCycle == W_CYCLE) cout << "W Multigrid Cycle, with " << nMGLevels << " multigrid levels."<< endl;
      if (MGCycle == FULLMG_CYCLE) cout << "Full Multigrid Cycle, with " << nMGLevels << " multigrid levels."<< endl;

      cout << "Damping factor for the residual restriction: " << Damp_Res_Restric <<"."<< endl;
      cout << "Damping factor for the correction prolongation: " << Damp_Correc_Prolong <<"."<< endl;
    }

    if ((Kind_Solver != FEM_ELASTICITY) && (Kind_Solver != HEAT_EQUATION) && (Kind_Solver != WAVE_EQUATION)) {

      if (!CFL_Adapt) cout << "No CFL adaptation." << endl;
      else cout << "CFL adaptation. Factor down: "<< CFL_AdaptParam[0] <<", factor up: "<< CFL_AdaptParam[1]
        <<",\n                lower limit: "<< CFL_AdaptParam[2] <<", upper limit: " << CFL_AdaptParam[3] <<"."<< endl;

      if (nMGLevels !=0) {
        cout << "Multigrid Level:                  ";
        for (unsigned short iLevel = 0; iLevel < nMGLevels+1; iLevel++) {
          cout.width(6); cout << iLevel;
        }
        cout << endl;
      }

			if (Unsteady_Simulation != TIME_STEPPING) {
				cout << "Courant-Friedrichs-Lewy number:   ";
				cout.precision(3);
				cout.width(6); cout << CFL[0];
				cout << endl;
			}
			

      if (nMGLevels !=0) {
        cout.precision(3);
        cout << "MG PreSmooth coefficients:        ";
        for (unsigned short iMG_PreSmooth = 0; iMG_PreSmooth < nMGLevels+1; iMG_PreSmooth++) {
          cout.width(6); cout << MG_PreSmooth[iMG_PreSmooth];
        }
        cout << endl;
      }

      if (nMGLevels !=0) {
        cout.precision(3);
        cout << "MG PostSmooth coefficients:       ";
        for (unsigned short iMG_PostSmooth = 0; iMG_PostSmooth < nMGLevels+1; iMG_PostSmooth++) {
          cout.width(6); cout << MG_PostSmooth[iMG_PostSmooth];
        }
        cout << endl;
      }

      if (nMGLevels !=0) {
        cout.precision(3);
        cout << "MG CorrecSmooth coefficients:     ";
        for (unsigned short iMG_CorrecSmooth = 0; iMG_CorrecSmooth < nMGLevels+1; iMG_CorrecSmooth++) {
          cout.width(6); cout << MG_CorrecSmooth[iMG_CorrecSmooth];
        }
        cout << endl;
      }

    }

    if ((Kind_Solver == RANS) || (Kind_Solver == DISC_ADJ_RANS))
      if (Kind_TimeIntScheme_Turb == EULER_IMPLICIT)
        cout << "Euler implicit time integration for the turbulence model." << endl;

    if (Kind_2phase_Model != NONE)
      if (Kind_TimeIntScheme_2phase == EULER_IMPLICIT)
        cout << "Euler implicit time integration for the 2-phase model." << endl;
  }

  if (val_software == SU2_CFD) {

    cout << endl <<"------------------------- Convergence Criteria --------------------------" << endl;

    cout << "Maximum number of iterations: " << nExtIter <<"."<< endl;

    if (ConvCriteria == CAUCHY) {
      if (!ContinuousAdjoint && !DiscreteAdjoint)
        switch (Cauchy_Func_Flow) {
          case LIFT_COEFFICIENT: cout << "Cauchy criteria for Lift using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case DRAG_COEFFICIENT: cout << "Cauchy criteria for Drag using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }

      if (ContinuousAdjoint || DiscreteAdjoint)
        switch (Cauchy_Func_AdjFlow) {
          case SENS_GEOMETRY: cout << "Cauchy criteria for geo. sensitivity using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case SENS_MACH: cout << "Cauchy criteria for Mach number sensitivity using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }

      cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
      
    }


    if (ConvCriteria == RESIDUAL) {
      if (!ContinuousAdjoint && !DiscreteAdjoint) {
        cout << "Reduce the density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum bound for the density residual is 10^(" << MinLogResidual<< ")."<< endl;
        cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
      }

      if (ContinuousAdjoint || DiscreteAdjoint) {
        cout << "Reduce the adjoint density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum value for the adjoint density residual is 10^(" << MinLogResidual<< ")."<< endl;
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
      case GRAD_FLOW: cout << "Grid adaptation using gradient based strategy (density)." << endl; break;
      case GRAD_ADJOINT: cout << "Grid adaptation using gradient based strategy (adjoint density)." << endl; break;
      case GRAD_FLOW_ADJ: cout << "Grid adaptation using gradient based strategy (density and adjoint density)." << endl; break;
      case COMPUTABLE: cout << "Grid adaptation using computable correction."<< endl; break;
      case REMAINING: cout << "Grid adaptation using remaining error."<< endl; break;
      case SMOOTHING: cout << "Grid smoothing using an implicit method."<< endl; break;
      case SUPERSONIC_SHOCK: cout << "Grid adaptation for a supersonic shock at Mach: " << Mach <<"."<< endl; break;
    }

    switch (Kind_Adaptation) {
      case GRAD_FLOW: case GRAD_ADJOINT: case GRAD_FLOW_ADJ: case COMPUTABLE: case REMAINING:
        cout << "Power of the dual volume in the adaptation sensor: " << DualVol_Power << endl;
        cout << "Percentage of new elements in the adaptation process: " << New_Elem_Adapt << "."<< endl;
        break;
    }

    if (Analytical_Surface != NONE)
      cout << "Use analytical definition for including points in the surfaces." << endl;

  }

  cout << endl <<"-------------------------- Output Information ---------------------------" << endl;

  if (val_software == SU2_CFD) {

    if (Low_MemoryOutput) cout << "Writing output files with low memory RAM requirements."<< endl;
    cout << "Writing a flow solution every " << Wrt_Sol_Freq <<" iterations."<< endl;
    cout << "Writing the convergence history every " << Wrt_Con_Freq <<" iterations."<< endl;
    if ((Unsteady_Simulation == DT_STEPPING_1ST) || (Unsteady_Simulation == DT_STEPPING_2ND)) {
      cout << "Writing the dual time flow solution every " << Wrt_Sol_Freq_DualTime <<" iterations."<< endl;
      cout << "Writing the dual time convergence history every " << Wrt_Con_Freq_DualTime <<" iterations."<< endl;
    }

    switch (Output_FileFormat) {
      case PARAVIEW: cout << "The output file format is Paraview ASCII (.vtk)." << endl; break;
      case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
      case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
      case FIELDVIEW: cout << "The output file format is FieldView ASCII (.uns)." << endl; break;
      case FIELDVIEW_BINARY: cout << "The output file format is FieldView binary (.uns)." << endl; break;
      case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
    }

    cout << "Convergence history file name: " << Conv_FileName << "." << endl;

    cout << "Forces breakdown file name: " << Breakdown_FileName << "." << endl;

    if ((Kind_Solver != FEM_ELASTICITY) && (Kind_Solver != HEAT_EQUATION) && (Kind_Solver != WAVE_EQUATION)) {
      if (!ContinuousAdjoint && !DiscreteAdjoint) {
        cout << "Surface flow coefficients file name: " << SurfFlowCoeff_FileName << "." << endl;
        cout << "Flow variables file name: " << Flow_FileName << "." << endl;
        cout << "Restart flow file name: " << Restart_FlowFileName << "." << endl;
      }

      if (ContinuousAdjoint || DiscreteAdjoint) {
        cout << "Adjoint solution file name: " << Solution_AdjFileName << "." << endl;
        cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
        cout << "Adjoint variables file name: " << Adj_FileName << "." << endl;
        cout << "Surface adjoint coefficients file name: " << SurfAdjCoeff_FileName << "." << endl;
      }
    }
    else {
      cout << "Surface structure coefficients file name: " << SurfStructure_FileName << "." << endl;
      cout << "Structure variables file name: " << Structure_FileName << "." << endl;
      cout << "Restart structure file name: " << Restart_FEMFileName << "." << endl;
    }

  }

  if (val_software == SU2_SOL) {
    if (Low_MemoryOutput) cout << "Writing output files with low memory RAM requirements."<< endl;
    switch (Output_FileFormat) {
      case PARAVIEW: cout << "The output file format is Paraview ASCII (.vtk)." << endl; break;
      case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
      case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
      case FIELDVIEW: cout << "The output file format is FieldView ASCII (.uns)." << endl; break;
      case FIELDVIEW_BINARY: cout << "The output file format is FieldView binary (.uns)." << endl; break;
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
    if (DiscreteAdjoint) {
      cout << "Output Volume Sensitivity file name: " << VolSens_FileName << ". " << endl;
      cout << "Output Surface Sensitivity file name: " << SurfSens_FileName << ". " << endl;
    }
    cout << "Output gradient file name: " << ObjFunc_Grad_FileName << ". " << endl;
  }

  if (val_software == SU2_MSH) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
    cout << "Restart flow file name: " << Restart_FlowFileName << "." << endl;
    if ((Kind_Adaptation == FULL_ADJOINT) || (Kind_Adaptation == GRAD_ADJOINT) || (Kind_Adaptation == GRAD_FLOW_ADJ) ||
        (Kind_Adaptation == COMPUTABLE) || (Kind_Adaptation == REMAINING)) {
      if (Kind_ObjFunc[0] == DRAG_COEFFICIENT) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc[0] == EQUIVALENT_AREA) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc[0] == NEARFIELD_PRESSURE) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc[0] == LIFT_COEFFICIENT) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
    }
  }

  cout << endl <<"------------------- Config File Boundary Information --------------------" << endl;

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
  
  if (nMarker_Fluid_InterfaceBound != 0) {
    cout << "Fluid interface boundary marker(s): ";
    for (iMarker_Fluid_InterfaceBound = 0; iMarker_Fluid_InterfaceBound < nMarker_Fluid_InterfaceBound; iMarker_Fluid_InterfaceBound++) {
      cout << Marker_Fluid_InterfaceBound[iMarker_Fluid_InterfaceBound];
      if (iMarker_Fluid_InterfaceBound < nMarker_Fluid_InterfaceBound-1) cout << ", ";
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
  
  if (nMarker_Internal != 0) {
    cout << "Internal boundary marker(s): ";
    for (iMarker_Internal = 0; iMarker_Internal < nMarker_Internal; iMarker_Internal++) {
      cout << Marker_Internal[iMarker_Internal];
      if (iMarker_Internal < nMarker_Internal-1) cout << ", ";
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
  
  if (nMarker_NRBC != 0) {
      cout << "NRBC boundary marker(s): ";
      for (iMarker_NRBC = 0; iMarker_NRBC < nMarker_NRBC; iMarker_NRBC++) {
        cout << Marker_NRBC[iMarker_NRBC];
        if (iMarker_NRBC < nMarker_NRBC-1) cout << ", ";
        else cout <<"."<< endl;
    }
  }

  if (nMarker_MixBound != 0) {
      cout << "MixingPlane boundary marker(s): ";
      for (iMarker_MixBound = 0; iMarker_MixBound < nMarker_MixBound; iMarker_MixBound++) {
        cout << Marker_MixBound[iMarker_MixBound];
        if (iMarker_MixBound < nMarker_MixBound-1) cout << ", ";
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
  
  if (nMarker_Supersonic_Outlet != 0) {
    cout << "Supersonic outlet boundary marker(s): ";
    for (iMarker_Supersonic_Outlet = 0; iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet; iMarker_Supersonic_Outlet++) {
      cout << Marker_Supersonic_Outlet[iMarker_Supersonic_Outlet];
      if (iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet-1) cout << ", ";
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

  if (nMarker_HeatFlux != 0) {
    cout << "Constant heat flux wall boundary marker(s): ";
    for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
      cout << Marker_HeatFlux[iMarker_HeatFlux];
      if (iMarker_HeatFlux < nMarker_HeatFlux-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Clamped != 0) {
    cout << "Clamped boundary marker(s): ";
    for (iMarker_Clamped = 0; iMarker_Clamped < nMarker_Clamped; iMarker_Clamped++) {
      cout << Marker_Clamped[iMarker_Clamped];
      if (iMarker_Clamped < nMarker_Clamped-1) cout << ", ";
      else cout <<"."<<endl;
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
    cout << "Normal load boundary marker(s): ";
    for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
      cout << Marker_Load[iMarker_Load];
      if (iMarker_Load < nMarker_Load-1) cout << ", ";
      else cout <<"."<< endl;
    }
  }

  if (nMarker_Load_Dir != 0) {
    cout << "Load boundary marker(s) in cartesian coordinates: ";
    for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++) {
      cout << Marker_Load_Dir[iMarker_Load_Dir];
      if (iMarker_Load_Dir < nMarker_Load_Dir-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }

  if (nMarker_Load_Sine != 0) {
    cout << "Sine-Wave Load boundary marker(s): ";
    for (iMarker_Load_Sine = 0; iMarker_Load_Sine < nMarker_Load_Sine; iMarker_Load_Sine++) {
      cout << Marker_Load_Sine[iMarker_Load_Sine];
      if (iMarker_Load_Sine < nMarker_Load_Sine-1) cout << ", ";
      else cout <<"."<<endl;
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

  if (nMarker_ActDiskInlet != 0) {
		cout << "Actuator disk (inlet) boundary marker(s): ";
		for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
			cout << Marker_ActDiskInlet[iMarker_ActDiskInlet];
			if (iMarker_ActDiskInlet < nMarker_ActDiskInlet-1) cout << ", ";
			else cout <<"."<< endl;
		}
	}

  if (nMarker_ActDiskOutlet != 0) {
		cout << "Actuator disk (outlet) boundary marker(s): ";
		for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
			cout << Marker_ActDiskOutlet[iMarker_ActDiskOutlet];
			if (iMarker_ActDiskOutlet < nMarker_ActDiskOutlet-1) cout << ", ";
			else cout <<"."<< endl;
		}
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
  value_part = str.substr(pos+1, string::npos);
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
    cerr << "Error in TokenizeString(): "
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
      string after_semi= it->substr(pos+1, string::npos);
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

  unsigned short iMarker_CfgFile;

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker)
      return iMarker_CfgFile;

  cout <<"The configuration file doesn't have any definition for marker "<< val_marker <<"!!" << endl;
  exit(EXIT_FAILURE);
  
}

string CConfig::GetMarker_CfgFile_TagBound(unsigned short val_marker) {
  return Marker_CfgFile_TagBound[val_marker];
}

unsigned short CConfig::GetMarker_CfgFile_KindBC(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_KindBC[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Monitoring(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Monitoring[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_GeoEval(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_GeoEval[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Designing(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Designing[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Plotting(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Plotting[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Analyze(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Analyze[iMarker_CfgFile];
}


unsigned short CConfig::GetMarker_CfgFile_FSIinterface(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_FSIinterface[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Out_1D(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Out_1D[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_DV(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_DV[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Moving(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Moving[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_PerBound(string val_marker) {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_PerBound[iMarker_CfgFile];
}

int CConfig::GetMarker_FSIinterface(string val_marker) {	
	  unsigned short iMarker_CfgFile;
	  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    
		  if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker)
				return  Marker_CfgFile_FSIinterface[iMarker_CfgFile];
    return 0;
}


CConfig::~CConfig(void) {
	
  unsigned long iDV, iMarker, iPeriodic, iFFD;

  /*--- Delete all of the option objects in the global option map ---*/
    
  for(map<string, COptionBase*>::iterator itr = option_map.begin(); itr != option_map.end(); itr++) {
    delete itr->second;
  }
 
  if (RK_Alpha_Step != NULL) delete [] RK_Alpha_Step;
  if (MG_PreSmooth  != NULL) delete [] MG_PreSmooth;
  if (MG_PostSmooth != NULL) delete [] MG_PostSmooth;
  
  /*--- Free memory for Aeroelastic problems. ---*/

  if (Grid_Movement && Aeroelastic_Simulation) {
    if (Aeroelastic_pitch  != NULL) delete[] Aeroelastic_pitch;
    if (Aeroelastic_plunge != NULL) delete[] Aeroelastic_plunge;
  }

  /*--- Free memory for unspecified grid motion parameters ---*/

 if (Kind_GridMovement != NULL) delete [] Kind_GridMovement;

 /*--- Free memory for airfoil sections ---*/

 if (LocationStations   != NULL) delete [] LocationStations;

  /*--- motion origin: ---*/
  
  if (Motion_Origin_X   != NULL) delete [] Motion_Origin_X;
  if (Motion_Origin_Y   != NULL) delete [] Motion_Origin_Y;
  if (Motion_Origin_Z   != NULL) delete [] Motion_Origin_Z;
  if (MoveMotion_Origin != NULL) delete [] MoveMotion_Origin;

  /*--- translation: ---*/
  
  if (Translation_Rate_X != NULL) delete [] Translation_Rate_X;
  if (Translation_Rate_Y != NULL) delete [] Translation_Rate_Y;
  if (Translation_Rate_Z != NULL) delete [] Translation_Rate_Z;

  /*--- rotation: ---*/
  
  if (Rotation_Rate_X != NULL) delete [] Rotation_Rate_X;
  if (Rotation_Rate_Y != NULL) delete [] Rotation_Rate_Y;
  if (Rotation_Rate_Z != NULL) delete [] Rotation_Rate_Z;

  /*--- pitching: ---*/
  
  if (Pitching_Omega_X != NULL) delete [] Pitching_Omega_X;
  if (Pitching_Omega_Y != NULL) delete [] Pitching_Omega_Y;
  if (Pitching_Omega_Z != NULL) delete [] Pitching_Omega_Z;

  /*--- pitching amplitude: ---*/
  
  if (Pitching_Ampl_X != NULL) delete [] Pitching_Ampl_X;
  if (Pitching_Ampl_Y != NULL) delete [] Pitching_Ampl_Y;
  if (Pitching_Ampl_Z != NULL) delete [] Pitching_Ampl_Z;

  /*--- pitching phase: ---*/
  
  if (Pitching_Phase_X != NULL) delete [] Pitching_Phase_X;
  if (Pitching_Phase_Y != NULL) delete [] Pitching_Phase_Y;
  if (Pitching_Phase_Z != NULL) delete [] Pitching_Phase_Z;

  /*--- plunging: ---*/
  
  if (Plunging_Omega_X != NULL) delete [] Plunging_Omega_X;
  if (Plunging_Omega_Y != NULL) delete [] Plunging_Omega_Y;
  if (Plunging_Omega_Z != NULL) delete [] Plunging_Omega_Z;

  /*--- plunging amplitude: ---*/
  
  if (Plunging_Ampl_X != NULL) delete [] Plunging_Ampl_X;
  if (Plunging_Ampl_Y != NULL) delete [] Plunging_Ampl_Y;
  if (Plunging_Ampl_Z != NULL) delete [] Plunging_Ampl_Z;

  /*--- reference origin for moments ---*/
  
  if (RefOriginMoment   != NULL) delete [] RefOriginMoment;
  if (RefOriginMoment_X != NULL) delete [] RefOriginMoment_X;
  if (RefOriginMoment_Y != NULL) delete [] RefOriginMoment_Y;
  if (RefOriginMoment_Z != NULL) delete [] RefOriginMoment_Z;

  /*--- Free memory for Harmonic Blance Frequency  pointer ---*/
    
  if (Omega_HB != NULL) delete [] Omega_HB;
    
  /*--- Marker pointers ---*/
  
  if (Marker_CfgFile_Out_1D != NULL) delete[] Marker_CfgFile_Out_1D;
  if (Marker_All_Out_1D     != NULL) delete[] Marker_All_Out_1D;
  
  if (Marker_CfgFile_GeoEval != NULL) delete[] Marker_CfgFile_GeoEval;
  if (Marker_All_GeoEval     != NULL) delete[] Marker_All_GeoEval;
  
  if (Marker_CfgFile_TagBound != NULL) delete[] Marker_CfgFile_TagBound;
  if (Marker_All_TagBound     != NULL) delete[] Marker_All_TagBound;
  
  if (Marker_CfgFile_KindBC != NULL) delete[] Marker_CfgFile_KindBC;
  if (Marker_All_KindBC     != NULL) delete[] Marker_All_KindBC;
  
  if (Marker_CfgFile_Monitoring != NULL) delete[] Marker_CfgFile_Monitoring;
  if (Marker_All_Monitoring     != NULL) delete[] Marker_All_Monitoring;
  
  if (Marker_CfgFile_Designing != NULL) delete[] Marker_CfgFile_Designing;
  if (Marker_All_Designing     != NULL) delete[] Marker_All_Designing;
  
  if (Marker_CfgFile_Plotting != NULL) delete[] Marker_CfgFile_Plotting;
  if (Marker_All_Plotting     != NULL) delete[] Marker_All_Plotting;
  
  if (Marker_CfgFile_Analyze != NULL) delete[] Marker_CfgFile_Analyze;
  if (Marker_All_Analyze  != NULL) delete[] Marker_All_Analyze;

  if (Marker_CfgFile_FSIinterface != NULL) delete[] Marker_CfgFile_FSIinterface;
  if (Marker_All_FSIinterface     != NULL) delete[] Marker_All_FSIinterface;
  
  if (Marker_CfgFile_DV != NULL) delete[] Marker_CfgFile_DV;
  if (Marker_All_DV     != NULL) delete[] Marker_All_DV;
  
  if (Marker_CfgFile_Moving != NULL) delete[] Marker_CfgFile_Moving;
  if (Marker_All_Moving     != NULL) delete[] Marker_All_Moving;
  
  if (Marker_CfgFile_PerBound != NULL) delete[] Marker_CfgFile_PerBound;
  if (Marker_All_PerBound     != NULL) delete[] Marker_All_PerBound;

  if (Marker_DV!= NULL)               delete[] Marker_DV;
  if (Marker_Moving != NULL)           delete[] Marker_Moving;
  if (Marker_Monitoring != NULL)      delete[] Marker_Monitoring;
  if (Marker_Designing != NULL)       delete[] Marker_Designing;
  if (Marker_GeoEval != NULL)         delete[] Marker_GeoEval;
  if (Marker_Plotting != NULL)        delete[] Marker_Plotting;
  if (Marker_Analyze != NULL)        delete[] Marker_Analyze;
  if (Marker_FSIinterface != NULL)        delete[] Marker_FSIinterface;
  if (Marker_All_SendRecv != NULL)    delete[] Marker_All_SendRecv;

  if (Kind_ObjFunc != NULL)      delete[] Kind_ObjFunc;
  if (Weight_ObjFunc != NULL)      delete[] Weight_ObjFunc;

  if (DV_Value != NULL) {
    for (iDV = 0; iDV < nDV; iDV++) delete[] DV_Value[iDV];
    delete [] DV_Value;
  }
  
  if (ParamDV != NULL) {
    for (iDV = 0; iDV < nDV; iDV++) delete[] ParamDV[iDV];
    delete [] ParamDV;
  }
  
  if (CoordFFDBox != NULL) {
    for (iFFD = 0; iFFD < nFFDBox; iFFD++) delete[] CoordFFDBox[iFFD];
    delete [] CoordFFDBox;
  }
  
  if (DegreeFFDBox != NULL) {
    for (iFFD = 0; iFFD < nFFDBox; iFFD++) delete[] DegreeFFDBox[iFFD];
    delete [] DegreeFFDBox;
  }
  
  if (Design_Variable != NULL)    delete[] Design_Variable;
  if (Dirichlet_Value != NULL)    delete[] Dirichlet_Value;
  
  if (Exhaust_Temperature_Target != NULL)    delete[]  Exhaust_Temperature_Target;
  if (Exhaust_Pressure_Target != NULL)    delete[]  Exhaust_Pressure_Target;
  if (Exhaust_Pressure != NULL)    delete[] Exhaust_Pressure;
  if (Exhaust_Temperature != NULL)    delete[] Exhaust_Temperature;
  if (Exhaust_MassFlow != NULL)    delete[] Exhaust_MassFlow;
  if (Exhaust_TotalPressure != NULL)    delete[] Exhaust_TotalPressure;
  if (Exhaust_TotalTemperature != NULL)    delete[] Exhaust_TotalTemperature;
  if (Exhaust_GrossThrust != NULL)    delete[] Exhaust_GrossThrust;
  if (Exhaust_Force != NULL)    delete[] Exhaust_Force;
  if (Exhaust_Power != NULL)    delete[] Exhaust_Power;

  if (Inflow_Mach != NULL)    delete[]  Inflow_Mach;
  if (Inflow_Pressure != NULL)    delete[] Inflow_Pressure;
  if (Inflow_MassFlow != NULL)    delete[] Inflow_MassFlow;
  if (Inflow_ReverseMassFlow != NULL)    delete[] Inflow_ReverseMassFlow;
  if (Inflow_TotalPressure != NULL)    delete[] Inflow_TotalPressure;
  if (Inflow_Temperature != NULL)    delete[] Inflow_Temperature;
  if (Inflow_TotalTemperature != NULL)    delete[] Inflow_TotalTemperature;
  if (Inflow_RamDrag != NULL)    delete[] Inflow_RamDrag;
  if (Inflow_Force != NULL)    delete[]  Inflow_Force;
  if (Inflow_Power != NULL)    delete[] Inflow_Power;

  if (Engine_Power != NULL)    delete[]  Engine_Power;
  if (Engine_Mach != NULL)    delete[]  Engine_Mach;
  if (Engine_Force != NULL)    delete[]  Engine_Force;
  if (Engine_NetThrust != NULL)    delete[]  Engine_NetThrust;
  if (Engine_GrossThrust != NULL)    delete[]  Engine_GrossThrust;
  if (Engine_Area != NULL)    delete[]  Engine_Area;
  if (EngineInflow_Target != NULL)    delete[] EngineInflow_Target;

  if (ActDiskInlet_MassFlow != NULL)    delete[]  ActDiskInlet_MassFlow;
  if (ActDiskInlet_Temperature != NULL)    delete[]  ActDiskInlet_Temperature;
  if (ActDiskInlet_TotalTemperature != NULL)    delete[]  ActDiskInlet_TotalTemperature;
  if (ActDiskInlet_Pressure != NULL)    delete[]  ActDiskInlet_Pressure;
  if (ActDiskInlet_TotalPressure != NULL)    delete[]  ActDiskInlet_TotalPressure;
  if (ActDiskInlet_RamDrag != NULL)    delete[]  ActDiskInlet_RamDrag;
  if (ActDiskInlet_Force != NULL)    delete[]  ActDiskInlet_Force;
  if (ActDiskInlet_Power != NULL)    delete[]  ActDiskInlet_Power;

  if (ActDiskOutlet_MassFlow != NULL)    delete[]  ActDiskOutlet_MassFlow;
  if (ActDiskOutlet_Temperature != NULL)    delete[]  ActDiskOutlet_Temperature;
  if (ActDiskOutlet_TotalTemperature != NULL)    delete[]  ActDiskOutlet_TotalTemperature;
  if (ActDiskOutlet_Pressure != NULL)    delete[]  ActDiskOutlet_Pressure;
  if (ActDiskOutlet_TotalPressure != NULL)    delete[]  ActDiskOutlet_TotalPressure;
  if (ActDiskOutlet_GrossThrust != NULL)    delete[]  ActDiskOutlet_GrossThrust;
  if (ActDiskOutlet_Force != NULL)    delete[]  ActDiskOutlet_Force;
  if (ActDiskOutlet_Power != NULL)    delete[]  ActDiskOutlet_Power;

  if (ActDisk_DeltaPress != NULL)    delete[]  ActDisk_DeltaPress;
  if (ActDisk_DeltaTemp != NULL)    delete[]  ActDisk_DeltaTemp;
  if (ActDisk_TotalPressRatio != NULL)    delete[]  ActDisk_TotalPressRatio;
  if (ActDisk_TotalTempRatio != NULL)    delete[]  ActDisk_TotalTempRatio;
  if (ActDisk_StaticPressRatio != NULL)    delete[]  ActDisk_StaticPressRatio;
  if (ActDisk_StaticTempRatio != NULL)    delete[]  ActDisk_StaticTempRatio;
  if (ActDisk_Power != NULL)    delete[]  ActDisk_Power;
  if (ActDisk_MassFlow != NULL)    delete[]  ActDisk_MassFlow;
  if (ActDisk_Mach != NULL)    delete[]  ActDisk_Mach;
  if (ActDisk_Force != NULL)    delete[]  ActDisk_Force;
  if (ActDisk_NetThrust != NULL)    delete[]  ActDisk_NetThrust;
  if (ActDisk_BCThrust != NULL)    delete[]  ActDisk_BCThrust;
  if (ActDisk_BCThrust_Old != NULL)    delete[]  ActDisk_BCThrust_Old;
  if (ActDisk_GrossThrust != NULL)    delete[]  ActDisk_GrossThrust;
  if (ActDisk_Area != NULL)    delete[]  ActDisk_Area;
  if (ActDisk_ReverseMassFlow != NULL)    delete[]  ActDisk_ReverseMassFlow;
  
  if (Surface_MassFlow != NULL)    delete[]  Surface_MassFlow;
  if (Surface_DC60 != NULL)    delete[]  Surface_DC60;
  if (Surface_IDC != NULL)    delete[]  Surface_IDC;
  if (Surface_IDC_Mach != NULL)    delete[]  Surface_IDC_Mach;
  if (Surface_IDR != NULL)    delete[]  Surface_IDR;

  if (Inlet_Ttotal != NULL) delete[]  Inlet_Ttotal;
  if (Inlet_Ptotal != NULL) delete[]  Inlet_Ptotal;
  if (Inlet_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker_Inlet; iMarker++)
      delete [] Inlet_FlowDir[iMarker];
    delete [] Inlet_FlowDir;
  }
  
  if (Inlet_Velocity != NULL) {
    for (iMarker = 0; iMarker < nMarker_Supersonic_Inlet; iMarker++)
      delete [] Inlet_Velocity[iMarker];
    delete [] Inlet_Velocity;
  }
  
  if (Riemann_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker_Riemann; iMarker++)
      delete [] Riemann_FlowDir[iMarker];
    delete [] Riemann_FlowDir;
  }
  
  if (NRBC_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker_NRBC; iMarker++)
      delete [] NRBC_FlowDir[iMarker];
    delete [] NRBC_FlowDir;
  }
  
  if (Load_Sine_Dir != NULL) {
    for (iMarker = 0; iMarker < nMarker_Load_Sine; iMarker++)
      delete [] Load_Sine_Dir[iMarker];
    delete [] Load_Sine_Dir;
  }
  
  if (Load_Dir != NULL) {
    for (iMarker = 0; iMarker < nMarker_Load_Dir; iMarker++)
      delete [] Load_Dir[iMarker];
    delete [] Load_Dir;
  }
  
  if (Inlet_Temperature != NULL)    delete[] Inlet_Temperature;
  if (Inlet_Pressure != NULL)    delete[] Inlet_Pressure;
  if (Outlet_Pressure != NULL)    delete[] Outlet_Pressure;
  if (Isothermal_Temperature != NULL)    delete[] Isothermal_Temperature;
  if (Heat_Flux != NULL)    delete[] Heat_Flux;
  if (Displ_Value != NULL)    delete[] Displ_Value;
  if (Load_Value != NULL)    delete[] Load_Value;
  if (Load_Dir_Multiplier != NULL)    delete[] Load_Dir_Multiplier;
  if (Load_Dir_Value != NULL)    delete[] Load_Dir_Value;
  if (Load_Sine_Amplitude != NULL)    delete[] Load_Sine_Amplitude;
  if (Load_Sine_Frequency != NULL)    delete[] Load_Sine_Frequency;
  if (FlowLoad_Value != NULL)    delete[] FlowLoad_Value;

  /*--- related to periodic boundary conditions ---*/
  
  for (iMarker = 0; iMarker < nMarker_PerBound; iMarker++) {
    if (Periodic_RotCenter   != NULL) delete [] Periodic_RotCenter[iMarker];
    if (Periodic_RotAngles   != NULL) delete [] Periodic_RotAngles[iMarker];
    if (Periodic_Translation != NULL) delete [] Periodic_Translation[iMarker];
  }
  if (Periodic_RotCenter   != NULL) delete[] Periodic_RotCenter;
  if (Periodic_RotAngles   != NULL) delete[] Periodic_RotAngles;
  if (Periodic_Translation != NULL) delete[] Periodic_Translation;

  for (iPeriodic = 0; iPeriodic < nPeriodic_Index; iPeriodic++) {
    if (Periodic_Center    != NULL) delete [] Periodic_Center[iPeriodic];
    if (Periodic_Rotation  != NULL) delete [] Periodic_Rotation[iPeriodic];
    if (Periodic_Translate != NULL) delete [] Periodic_Translate[iPeriodic];
  }
  if (Periodic_Center      != NULL) delete[] Periodic_Center;
  if (Periodic_Rotation    != NULL) delete[] Periodic_Rotation;
  if (Periodic_Translate   != NULL) delete[] Periodic_Translate;
  
  if (MG_CorrecSmooth != NULL) delete[] MG_CorrecSmooth;
  if (PlaneTag != NULL)        delete[] PlaneTag;
  if (CFL != NULL)             delete[] CFL;

  /*--- String markers ---*/
  
  if (Marker_Euler != NULL )              delete[] Marker_Euler;
  if (Marker_FarField != NULL )           delete[] Marker_FarField;
  if (Marker_Custom != NULL )             delete[] Marker_Custom;
  if (Marker_SymWall != NULL )            delete[] Marker_SymWall;
  if (Marker_Pressure != NULL )           delete[] Marker_Pressure;
  if (Marker_PerBound != NULL )           delete[] Marker_PerBound;
  if (Marker_PerDonor != NULL )           delete[] Marker_PerDonor;
  if (Marker_NearFieldBound != NULL )     delete[] Marker_NearFieldBound;
  if (Marker_InterfaceBound != NULL )     delete[] Marker_InterfaceBound;
  if (Marker_Fluid_InterfaceBound != NULL )     delete[] Marker_Fluid_InterfaceBound;
  if (Marker_Dirichlet != NULL )          delete[] Marker_Dirichlet;
  if (Marker_Inlet != NULL )              delete[] Marker_Inlet;
  if (Marker_Supersonic_Inlet != NULL )   delete[] Marker_Supersonic_Inlet;
  if (Marker_Supersonic_Outlet != NULL )   delete[] Marker_Supersonic_Outlet;
  if (Marker_Outlet != NULL )             delete[] Marker_Outlet;
  if (Marker_Out_1D != NULL )             delete[] Marker_Out_1D;
  if (Marker_Isothermal != NULL )         delete[] Marker_Isothermal;
  if (Marker_EngineInflow != NULL )      delete[] Marker_EngineInflow;
  if (Marker_EngineExhaust != NULL )     delete[] Marker_EngineExhaust;
  if (Marker_Displacement != NULL )       delete[] Marker_Displacement;
  if (Marker_Load != NULL )               delete[] Marker_Load;
  if (Marker_Load_Dir != NULL )               delete[] Marker_Load_Dir;
  if (Marker_Load_Sine != NULL )               delete[] Marker_Load_Sine;
  if (Marker_FlowLoad != NULL )           delete[] Marker_FlowLoad;
  if (Marker_Neumann != NULL )            delete[] Marker_Neumann;
  if (Marker_Internal != NULL )            delete[] Marker_Internal;
  if (Marker_HeatFlux != NULL )               delete[] Marker_HeatFlux;

  if (Int_Coeffs != NULL) delete [] Int_Coeffs;
  
  /*--- Delete some arrays needed just for initializing options. ---*/
  
  if (default_vel_inf       != NULL) delete [] default_vel_inf;
  if (default_ffd_axis      != NULL) delete [] default_ffd_axis;
  if (default_eng_cyl       != NULL) delete [] default_eng_cyl;
  if (default_eng_val       != NULL) delete [] default_eng_val;
  if (default_heat_capacity != NULL) delete [] default_heat_capacity;
  if (default_cfl_adapt     != NULL) delete [] default_cfl_adapt;
  if (default_ad_coeff_flow != NULL) delete [] default_ad_coeff_flow;
  if (default_ad_coeff_adj  != NULL) delete [] default_ad_coeff_adj;
  if (default_obj_coeff     != NULL) delete [] default_obj_coeff;
  if (default_geo_loc       != NULL) delete [] default_geo_loc;
  if (default_distortion    != NULL) delete [] default_distortion;
  if (default_ea_lim        != NULL) delete [] default_ea_lim;
  if (default_grid_fix      != NULL) delete [] default_grid_fix;
  if (default_inc_crit      != NULL) delete [] default_inc_crit;
  if (default_htp_axis      != NULL) delete [] default_htp_axis;
  if (default_body_force    != NULL) delete [] default_body_force;

  if (FFDTag != NULL) delete [] FFDTag;
  if (nDV_Value != NULL) delete [] nDV_Value;
  if (TagFFDBox != NULL) delete [] TagFFDBox;
  
  if (Kind_Data_Riemann != NULL) delete [] Kind_Data_Riemann;
  if (Riemann_Var1 != NULL) delete [] Riemann_Var1;
  if (Riemann_Var2 != NULL) delete [] Riemann_Var2;
  if (Kind_Data_NRBC != NULL) delete [] Kind_Data_NRBC;
  if (NRBC_Var1 != NULL) delete [] NRBC_Var1;
  if (NRBC_Var2 != NULL) delete [] NRBC_Var2;
  if (Marker_TurboBoundIn != NULL) delete [] Marker_TurboBoundIn;
  if (Marker_TurboBoundOut != NULL) delete [] Marker_TurboBoundOut;
  if (Kind_TurboPerformance != NULL) delete [] Kind_TurboPerformance;

  if (Marker_Riemann != NULL) delete [] Marker_Riemann;
  if (Marker_NRBC != NULL) delete [] Marker_NRBC;
 
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

  if ((Wrt_Unsteady) || (Wrt_Dynamic)) {
    unsigned short lastindex = UnstFilename.find_last_of(".");
    UnstFilename = UnstFilename.substr(0, lastindex);
    if ((val_iter >= 0)    && (val_iter < 10))    SPRINTF (buffer, "_0000%d.dat", val_iter);
    if ((val_iter >= 10)   && (val_iter < 100))   SPRINTF (buffer, "_000%d.dat",  val_iter);
    if ((val_iter >= 100)  && (val_iter < 1000))  SPRINTF (buffer, "_00%d.dat",   val_iter);
    if ((val_iter >= 1000) && (val_iter < 10000)) SPRINTF (buffer, "_0%d.dat",    val_iter);
    if (val_iter >= 10000) SPRINTF (buffer, "_%d.dat", val_iter);
    string UnstExt = string(buffer);
    UnstFilename.append(UnstExt);
  }

  return UnstFilename;
}

string CConfig::GetMultizone_FileName(string val_filename, int val_iZone) {

    string multizone_filename = val_filename;
    char buffer[50];
    
    if (GetnZone() > 1 ) {
        unsigned short lastindex = multizone_filename.find_last_of(".");
        multizone_filename = multizone_filename.substr(0, lastindex);
        SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
        multizone_filename.append(string(buffer));
    }
    
    return multizone_filename;
}

string CConfig::GetObjFunc_Extension(string val_filename) {

  string AdjExt, Filename = val_filename;

  if (ContinuousAdjoint || DiscreteAdjoint) {

    /*--- Remove filename extension (.dat) ---*/
    unsigned short lastindex = Filename.find_last_of(".");
    Filename = Filename.substr(0, lastindex);
    if (nObj==1) {
      switch (Kind_ObjFunc[0]) {
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
      case AVG_TOTAL_PRESSURE:      AdjExt = "_pt";       break;
      case AVG_OUTLET_PRESSURE:     AdjExt = "_pe";       break;
      case MASS_FLOW_RATE:          AdjExt = "_mfr";      break;
      case OUTFLOW_GENERALIZED:     AdjExt = "_chn";      break;
      case AERO_DRAG_COEFFICIENT:   AdjExt = "_acd";       break;
      case RADIAL_DISTORTION:           AdjExt = "_rdis";      break;
      case CIRCUMFERENTIAL_DISTORTION:  AdjExt = "_cdis";      break;
      }
    }
    else{
      AdjExt = "_combo";
    }
    Filename.append(AdjExt);

    /*--- Lastly, add the .dat extension ---*/
    Filename.append(".dat");

  }

  return Filename;
}

unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {

  switch (val_eqsystem) {
    case RUNTIME_FLOW_SYS:        return FLOW_SOL;
    case RUNTIME_TURB_SYS:        return TURB_SOL;
    case RUNTIME_2PHASE_SYS:      return TWO_PHASE_SOL;
    case RUNTIME_TRANS_SYS:       return TRANS_SOL;
    case RUNTIME_POISSON_SYS:     return POISSON_SOL;
    case RUNTIME_WAVE_SYS:        return WAVE_SOL;
    case RUNTIME_HEAT_SYS:        return HEAT_SOL;
    case RUNTIME_FEA_SYS:         return FEA_SOL;
    case RUNTIME_ADJPOT_SYS:      return ADJFLOW_SOL;
    case RUNTIME_ADJFLOW_SYS:     return ADJFLOW_SOL;
    case RUNTIME_ADJTURB_SYS:     return ADJTURB_SOL;
    case RUNTIME_MULTIGRID_SYS:   return 0;
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

void CConfig::SetGlobalParam(unsigned short val_solver,
                             unsigned short val_system,
                             unsigned long val_extiter) {

  /*--- Set the simulation global time ---*/
  Current_UnstTime = static_cast<su2double>(val_extiter)*Delta_UnstTime;
  Current_UnstTimeND = static_cast<su2double>(val_extiter)*Delta_UnstTimeND;

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
    case TWO_PHASE_EULER:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_2PHASE_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_2phase, Kind_Centered_2phase,
                              Kind_Upwind_2phase, Kind_SlopeLimit_2phase,
                              SpatialOrder_2phase);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_2phase);
      }
      break;
    case TWO_PHASE_NAVIER_STOKES:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              SpatialOrder_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
      }
      if (val_system == RUNTIME_2PHASE_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_2phase, Kind_Centered_2phase,
                              Kind_Upwind_2phase, Kind_SlopeLimit_2phase,
                              SpatialOrder_2phase);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_2phase);
      }
      break;
    case TWO_PHASE_RANS:
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
      if (val_system == RUNTIME_2PHASE_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_2phase, Kind_Centered_2phase,
                              Kind_Upwind_2phase, Kind_SlopeLimit_2phase,
                              SpatialOrder_2phase);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_2phase);
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
    case FEM_ELASTICITY:

      Current_DynTime = static_cast<su2double>(val_extiter)*Delta_DynTime;

      if (val_system == RUNTIME_FEA_SYS) {
        SetKind_ConvNumScheme(NONE, NONE, NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_FEA);
      }
      break;
  }
}

su2double* CConfig::GetPeriodicRotCenter(string val_marker) {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_RotCenter[iMarker_PerBound];
}

su2double* CConfig::GetPeriodicRotAngles(string val_marker) {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_RotAngles[iMarker_PerBound];
}

su2double* CConfig::GetPeriodicTranslation(string val_marker) {
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
  for (kMarker_All = 0; kMarker_All < nMarker_CfgFile; kMarker_All++)
    if (Marker_PerBound[jMarker_PerBound] == Marker_All_TagBound[kMarker_All]) break;

  return kMarker_All;
}

su2double CConfig::GetActDisk_NetThrust(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_NetThrust[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_Power(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Power[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_MassFlow(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_MassFlow[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_Mach(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Mach[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_Force(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Force[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_BCThrust(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_BCThrust[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_BCThrust_Old(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_BCThrust_Old[iMarker_ActDisk];
}

void CConfig::SetActDisk_BCThrust(string val_marker, su2double val_actdisk_bcthrust) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  ActDisk_BCThrust[iMarker_ActDisk] = val_actdisk_bcthrust;
}

void CConfig::SetActDisk_BCThrust_Old(string val_marker, su2double val_actdisk_bcthrust_old) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  ActDisk_BCThrust_Old[iMarker_ActDisk] = val_actdisk_bcthrust_old;
}

su2double CConfig::GetActDisk_Area(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Area[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_ReverseMassFlow(string val_marker) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_ReverseMassFlow[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_PressJump(string val_marker, unsigned short val_value) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_PressJump[iMarker_ActDisk][val_value];
}

su2double CConfig::GetActDisk_TempJump(string val_marker, unsigned short val_value) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_TempJump[iMarker_ActDisk][val_value];;
}

su2double CConfig::GetActDisk_Omega(string val_marker, unsigned short val_value) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Omega[iMarker_ActDisk][val_value];;
}

unsigned short CConfig::GetMarker_CfgFile_ActDiskOutlet(string val_marker) {
  unsigned short iMarker_ActDisk, kMarker_All;
  
  /*--- Find the marker for this actuator disk inlet. ---*/
  
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if (Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) break;
  
  /*--- Find and return global marker index for the actuator disk outlet. ---*/
  
  for (kMarker_All = 0; kMarker_All < nMarker_CfgFile; kMarker_All++)
    if (Marker_ActDiskOutlet[iMarker_ActDisk] == Marker_CfgFile_TagBound[kMarker_All]) break;
  
  return kMarker_All;
}

unsigned short CConfig::GetMarker_CfgFile_EngineExhaust(string val_marker) {
  unsigned short iMarker_Engine, kMarker_All;
  
  /*--- Find the marker for this engine inflow. ---*/
  
  for (iMarker_Engine = 0; iMarker_Engine < nMarker_EngineInflow; iMarker_Engine++)
    if (Marker_EngineInflow[iMarker_Engine] == val_marker) break;
  
  /*--- Find and return global marker index for the engine exhaust. ---*/
  
  for (kMarker_All = 0; kMarker_All < nMarker_CfgFile; kMarker_All++)
    if (Marker_EngineExhaust[iMarker_Engine] == Marker_CfgFile_TagBound[kMarker_All]) break;
  
  return kMarker_All;
}

void CConfig::SetnPeriodicIndex(unsigned short val_index) {

  /*--- Store total number of transformations. ---*/
  nPeriodic_Index = val_index;

  /*--- Allocate memory for centers, angles, translations. ---*/
  Periodic_Center    = new su2double*[nPeriodic_Index];
  Periodic_Rotation  = new su2double*[nPeriodic_Index];
  Periodic_Translate = new su2double*[nPeriodic_Index];
  
  for (unsigned long i = 0; i < nPeriodic_Index; i++) {
    Periodic_Center[i]    = new su2double[3];
    Periodic_Rotation[i]  = new su2double[3];
    Periodic_Translate[i] = new su2double[3];
  }
  
}

unsigned short CConfig::GetMarker_Moving(string val_marker) {
  unsigned short iMarker_Moving;

  /*--- Find the marker for this moving boundary. ---*/
  for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
    if (Marker_Moving[iMarker_Moving] == val_marker) break;

  return iMarker_Moving;
}

su2double CConfig::GetDirichlet_Value(string val_marker) {
  unsigned short iMarker_Dirichlet;
  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++)
    if (Marker_Dirichlet[iMarker_Dirichlet] == val_marker) break;
  return Dirichlet_Value[iMarker_Dirichlet];
}

bool CConfig::GetDirichlet_Boundary(string val_marker) {
  unsigned short iMarker_Dirichlet;
  bool Dirichlet = false;
  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++)
    if (Marker_Dirichlet[iMarker_Dirichlet] == val_marker) {
      Dirichlet = true;
      break;
    }
  return Dirichlet;
}

su2double CConfig::GetExhaust_Temperature_Target(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Temperature_Target[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Pressure_Target(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Pressure_Target[iMarker_EngineExhaust];
}

su2double CConfig::GetInlet_Ttotal(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ttotal[iMarker_Inlet];
}

su2double CConfig::GetInlet_Ptotal(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ptotal[iMarker_Inlet];
}

su2double* CConfig::GetInlet_FlowDir(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_FlowDir[iMarker_Inlet];
}

su2double CConfig::GetInlet_Temperature(string val_marker) {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Temperature[iMarker_Supersonic_Inlet];
}

su2double CConfig::GetInlet_Pressure(string val_marker) {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Pressure[iMarker_Supersonic_Inlet];
}

su2double* CConfig::GetInlet_Velocity(string val_marker) {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Velocity[iMarker_Supersonic_Inlet];
}

su2double CConfig::GetOutlet_Pressure(string val_marker) {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
  return Outlet_Pressure[iMarker_Outlet];
}

su2double CConfig::GetRiemann_Var1(string val_marker) {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_Var1[iMarker_Riemann];
}

su2double CConfig::GetRiemann_Var2(string val_marker) {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_Var2[iMarker_Riemann];
}

su2double* CConfig::GetRiemann_FlowDir(string val_marker) {
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


su2double CConfig::GetNRBC_Var1(string val_marker) {
  unsigned short iMarker_NRBC;
  for (iMarker_NRBC = 0; iMarker_NRBC < nMarker_NRBC; iMarker_NRBC++)
    if (Marker_NRBC[iMarker_NRBC] == val_marker) break;
  return NRBC_Var1[iMarker_NRBC];
}

su2double CConfig::GetNRBC_Var2(string val_marker) {
  unsigned short iMarker_NRBC;
  for (iMarker_NRBC = 0; iMarker_NRBC < nMarker_NRBC; iMarker_NRBC++)
    if (Marker_NRBC[iMarker_NRBC] == val_marker) break;
  return NRBC_Var2[iMarker_NRBC];
}

su2double* CConfig::GetNRBC_FlowDir(string val_marker) {
  unsigned short iMarker_NRBC;
  for (iMarker_NRBC = 0; iMarker_NRBC < nMarker_NRBC; iMarker_NRBC++)
    if (Marker_NRBC[iMarker_NRBC] == val_marker) break;
  return NRBC_FlowDir[iMarker_NRBC];
}

unsigned short CConfig::GetKind_Data_NRBC(string val_marker) {
  unsigned short iMarker_NRBC;
  for (iMarker_NRBC = 0; iMarker_NRBC < nMarker_NRBC; iMarker_NRBC++)
    if (Marker_NRBC[iMarker_NRBC] == val_marker) break;
  return Kind_Data_NRBC[iMarker_NRBC];
}


su2double CConfig::GetIsothermal_Temperature(string val_marker) {

  unsigned short iMarker_Isothermal = 0;

  if (nMarker_Isothermal > 0) {
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++)
      if (Marker_Isothermal[iMarker_Isothermal] == val_marker) break;
  }

  return Isothermal_Temperature[iMarker_Isothermal];
}

su2double CConfig::GetWall_HeatFlux(string val_marker) {
  unsigned short iMarker_HeatFlux = 0;

  if (nMarker_HeatFlux > 0) {
  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++)
    if (Marker_HeatFlux[iMarker_HeatFlux] == val_marker) break;
  }

  return Heat_Flux[iMarker_HeatFlux];
}

su2double CConfig::GetEngineInflow_Target(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return EngineInflow_Target[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Pressure(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Pressure[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_MassFlow(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_MassFlow[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_ReverseMassFlow(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_ReverseMassFlow[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_TotalPressure(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_TotalPressure[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Temperature(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Temperature[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_TotalTemperature(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_TotalTemperature[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_RamDrag(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_RamDrag[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Force(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Force[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Power(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Power[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Mach(string val_marker) {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Mach[iMarker_EngineInflow];
}

su2double CConfig::GetExhaust_Pressure(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Pressure[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Temperature(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Temperature[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_MassFlow(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_MassFlow[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_TotalPressure(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_TotalPressure[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_TotalTemperature(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_TotalTemperature[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_GrossThrust(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_GrossThrust[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Force(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Force[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Power(string val_marker) {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Power[iMarker_EngineExhaust];
}

su2double CConfig::GetActDiskInlet_Pressure(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Pressure[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_TotalPressure(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_TotalPressure[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_RamDrag(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_RamDrag[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_Force(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Force[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_Power(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Power[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskOutlet_Pressure(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Pressure[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_TotalPressure(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_TotalPressure[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_GrossThrust(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_GrossThrust[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_Force(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Force[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_Power(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Power[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskInlet_Temperature(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Temperature[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_TotalTemperature(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_TotalTemperature[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskOutlet_Temperature(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Temperature[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_TotalTemperature(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_TotalTemperature[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskInlet_MassFlow(string val_marker) {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_MassFlow[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskOutlet_MassFlow(string val_marker) {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_MassFlow[iMarker_ActDiskOutlet];
}

su2double CConfig::GetDispl_Value(string val_marker) {
  unsigned short iMarker_Displacement;
  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++)
    if (Marker_Displacement[iMarker_Displacement] == val_marker) break;
  return Displ_Value[iMarker_Displacement];
}

su2double CConfig::GetLoad_Value(string val_marker) {
  unsigned short iMarker_Load;
  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++)
    if (Marker_Load[iMarker_Load] == val_marker) break;
  return Load_Value[iMarker_Load];
}

su2double CConfig::GetLoad_Dir_Value(string val_marker) {
  unsigned short iMarker_Load_Dir;
  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++)
    if (Marker_Load_Dir[iMarker_Load_Dir] == val_marker) break;
  return Load_Dir_Value[iMarker_Load_Dir];
}

su2double CConfig::GetLoad_Dir_Multiplier(string val_marker) {
  unsigned short iMarker_Load_Dir;
  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++)
    if (Marker_Load_Dir[iMarker_Load_Dir] == val_marker) break;
  return Load_Dir_Multiplier[iMarker_Load_Dir];
}

su2double* CConfig::GetLoad_Dir(string val_marker) {
  unsigned short iMarker_Load_Dir;
  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++)
    if (Marker_Load_Dir[iMarker_Load_Dir] == val_marker) break;
  return Load_Dir[iMarker_Load_Dir];
}


su2double CConfig::GetLoad_Sine_Amplitude(string val_marker) {
  unsigned short iMarker_Load_Sine;
  for (iMarker_Load_Sine = 0; iMarker_Load_Sine < nMarker_Load_Sine; iMarker_Load_Sine++)
    if (Marker_Load_Sine[iMarker_Load_Sine] == val_marker) break;
  return Load_Sine_Amplitude[iMarker_Load_Sine];
}

su2double CConfig::GetLoad_Sine_Frequency(string val_marker) {
  unsigned short iMarker_Load_Sine;
  for (iMarker_Load_Sine = 0; iMarker_Load_Sine < nMarker_Load_Sine; iMarker_Load_Sine++)
    if (Marker_Load_Sine[iMarker_Load_Sine] == val_marker) break;
  return Load_Sine_Frequency[iMarker_Load_Sine];
}

su2double* CConfig::GetLoad_Sine_Dir(string val_marker) {
  unsigned short iMarker_Load_Sine;
  for (iMarker_Load_Sine = 0; iMarker_Load_Sine < nMarker_Load_Sine; iMarker_Load_Sine++)
    if (Marker_Load_Sine[iMarker_Load_Sine] == val_marker) break;
  return Load_Sine_Dir[iMarker_Load_Sine];
}

su2double CConfig::GetFlowLoad_Value(string val_marker) {
  unsigned short iMarker_FlowLoad;
  for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++)
    if (Marker_FlowLoad[iMarker_FlowLoad] == val_marker) break;
  return FlowLoad_Value[iMarker_FlowLoad];
}

void CConfig::SetSpline(vector<su2double> &x, vector<su2double> &y, unsigned long n, su2double yp1, su2double ypn, vector<su2double> &y2) {
  unsigned long i, k;
  su2double p, qn, sig, un, *u;

  u = new su2double [n];

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

su2double CConfig::GetSpline(vector<su2double>&xa, vector<su2double>&ya, vector<su2double>&y2a, unsigned long n, su2double x) {
  unsigned long klo, khi, k;
  su2double h, b, a, y;

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
