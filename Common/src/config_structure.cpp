/*!
 * \file config_structure.cpp
 * \brief Main file for reading the config file.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

CConfig::CConfig(char case_filename[200], unsigned short val_software) {
  
  /*--- Store the component of SU2. ---*/
  Kind_SU2 = val_software;
  
  string text_line, text2find, option_name, keyword;
	ifstream case_file;
	vector<string> option_value;
	
	double default_vec_3d[3];
	string string_array[2];
	
	/*--- Definition of default values ---*/
	
	/*--- options related to problem definition and partitioning ---*/
	AddEnumOption("PHYSICAL_PROBLEM", Kind_Solver, Solver_Map, "NONE");
	AddEnumOption("GAS_MODEL", Kind_GasModel, GasModel_Map, "ARGON");        
	AddMathProblem("MATH_PROBLEM" , ContAdj, false , OneShot, false, 
		       Linerized, false, Restart_Flow, false);
	AddSpecialOption("INCOMPRESSIBLE_FORMULATION", Incompressible, SetBoolOption, false);
	AddSpecialOption("RESTART_SOL", Restart, SetBoolOption, false);
	AddScalarOption("NUMBER_PART", nDomain, 0);
	AddSpecialOption("VISUALIZE_PART", Visualize_Partition, SetBoolOption, false);
	AddSpecialOption("DIVIDE_ELEMENTS", Divide_Element, SetBoolOption, false);
	AddEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, "NONE");
	AddSpecialOption("GRAVITY_FORCE", GravityForce, SetBoolOption, false);
	AddSpecialOption("SHOCKTUBE_PROBLEM", ShockTube, SetBoolOption, false);
	AddScalarOption("NUMBER_OF_SPECIES_IN_PLASMA", nSpecies, 0);
	AddScalarOption("NUMBER_OF_FLUIDS_IN_PLASMA", nFluids, 0);
	
	/*--- options related to various boundary markers ---*/
	AddMarkerOption("MARKER_MONITORING", nMarker_Monitoring, Marker_Monitoring);
	AddMarkerOption("MARKER_PLOTTING", nMarker_Plotting, Marker_Plotting);
	AddMarkerOption("DV_MARKER", nMarker_Moving, Marker_Moving);
	AddMarkerOption("MARKER_EULER", nMarker_Euler, Marker_Euler);
	AddMarkerOption("MARKER_NS", nMarker_NS, Marker_NS);
	AddMarkerOption("MARKER_FAR", nMarker_FarField, Marker_FarField);
	AddMarkerOption("MARKER_SYM", nMarker_SymWall, Marker_SymWall);
	AddMarkerOption("MARKER_NEARFIELD", nMarker_NearFieldBound, Marker_NearFieldBound);
	AddMarkerOption("MARKER_INTERFACE", nMarker_InterfaceBound, Marker_InterfaceBound);
	AddMarkerOption("MARKER_DIRICHLET", nMarker_Dirichlet, Marker_Dirichlet);
	AddMarkerOption("MARKER_CUSTOM", nMarker_Custom, Marker_Custom);
	AddMarkerPeriodic("MARKER_PERIODIC", nMarker_PerBound, Marker_PerBound, Marker_PerDonor,
										Periodic_RotCenter, Periodic_RotAngles, Periodic_Translation);
	AddMarkerInlet("MARKER_INLET", nMarker_Inlet, Marker_Inlet, Inlet_Ttotal, Inlet_Ptotal,
								 Inlet_FlowDir);
	AddMarkerOutlet("MARKER_OUTLET", nMarker_Outlet, Marker_Outlet, Outlet_Pressure);
	
	/*--- options related to grid adaptation ---*/
	AddEnumOption("KIND_ADAPT", Kind_Adaptation, Adapt_Map, "NONE");
	AddScalarOption("NEW_ELEMS", New_Elem_Adapt, -1.0);
	AddScalarOption("DUALVOL_POWER", DualVol_Power, 0.5);        
	AddEnumOption("ANALYTICAL_SURFDEF", Analytical_Surface, Geo_Analytic_Map, "NONE");        
	AddSpecialOption("SMOOTH_GEOMETRY", SmoothNumGrid, SetBoolOption, false);
	AddSpecialOption("VISUALIZE_DEFORMATION", Visualize_Deformation, SetBoolOption, false);
	
	/*--- options related to time-marching ---*/
	AddEnumOption("UNSTEADY_SIMULATION", Unsteady_Simulation, Unsteady_Map, "NO");
	AddScalarOption("CFL_NUMBER", CFLFineGrid, 1.25);
	default_vec_3d[0] = 1.0; default_vec_3d[1] = 100.0; default_vec_3d[2] = 1.0;
	AddArrayOption("CFL_RAMP", 3, CFLRamp, default_vec_3d);
	AddScalarOption("EXT_ITER", nExtIter, 999999);
	// these options share nRKStep as their size, which is not a good idea in general
	AddListOption("RK_ALPHA_COEFF", nRKStep, RK_Alpha_Step);
	AddListOption("RK_BETA_COEFF", nRKStep, RK_Beta_Step);
	nIntIter = 1; // NOTE: should this be a formal option?
	AddScalarOption("UNST_TIMESTEP", Unst_TimeStep, 0.1);
	AddScalarOption("UNST_TIME", Unst_Time, 1.0);
	AddEnumOption("TIME_DISCRE_FLOW", Kind_TimeIntScheme_Flow, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	AddEnumOption("TIME_DISCRE_COMB", Kind_TimeIntScheme_Combustion, Time_Int_Map, "EULER_EXPLICIT");
	AddEnumOption("TIME_DISCRE_LEVELSET", Kind_TimeIntScheme_LevelSet, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	AddEnumOption("TIME_DISCRE_PLASMA", Kind_TimeIntScheme_Plasma, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	AddEnumOption("TIME_DISCRE_ADJ", Kind_TimeIntScheme_Adj, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	AddEnumOption("TIME_DISCRE_LIN", Kind_TimeIntScheme_Lin, Time_Int_Map, "RUNGE-KUTTA_EXPLICIT");
	AddEnumOption("TIME_DISCRE_TURB", Kind_TimeIntScheme_Turb, Time_Int_Map, "EULER_IMPLICIT");
	AddEnumOption("TIME_DISCRE_ADJTURB", Kind_TimeIntScheme_AdjTurb, Time_Int_Map, "EULER_IMPLICIT");
	
	/*--- options related to convergence ---*/
	AddEnumOption("CONV_CRITERIA", ConvCriteria, Converge_Crit_Map, "RESIDUAL");
	AddScalarOption("CAUCHY_EPS", Cauchy_Eps, 1E-10);
	AddScalarOption("RESIDUAL_REDUCTION", OrderMagResidual, 3.0);
	AddScalarOption("RESIDUAL_MINVAL", MinLogResidual, -8.0);
	AddScalarOption("STARTCONV_ITER", StartConv_Iter, 5);
	AddScalarOption("CAUCHY_ELEMS", Cauchy_Elems, 100);
	AddEnumOption("CAUCHY_FUNC_FLOW", Cauchy_Func_Flow, Objective_Map, "DRAG");
	AddEnumOption("CAUCHY_FUNC_ADJ", Cauchy_Func_AdjFlow, Sens_Map, "SENS_GEOMETRY");
	AddEnumOption("CAUCHY_FUNC_LIN", Cauchy_Func_LinFlow, Linear_Obj_Map, "DELTA_DRAG");
	AddScalarOption("ONESHOT_CAUCHY_EPS", Cauchy_Eps_OneShot, 1E-5);
	AddScalarOption("FULLMG_CAUCHY_EPS", Cauchy_Eps_FullMG, 1E-4);
	
	/*--- options related to Multi-grid ---*/
	AddSpecialOption("FULLMG", FullMG, SetBoolOption, false);
	AddScalarOption("START_UP_ITER", nStartUpIter, 0);
	AddScalarOption("MGLEVEL", nMultiLevel, 0);
	AddScalarOption("MGCYCLE", MGCycle, 0);
	AddListOption("MG_PRE_SMOOTH", nMG_PreSmooth, MG_PreSmooth);
	AddListOption("MG_POST_SMOOTH", nMG_PostSmooth, MG_PostSmooth);
	AddListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth, MG_CorrecSmooth);       
	AddScalarOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.75);
	AddScalarOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.75);
	AddSpecialOption("MG_RESTART_CYCLE", ReStartMGCycle, SetBoolOption, true);
	AddScalarOption("MG_CFL_REDUCTION", MG_CFLRedCoeff, 0.8);        
	AddScalarOption("MAX_CHILDREN", MaxChildren, 50);
	AddScalarOption("MAX_DIMENSION", MaxDimension, 0.04);
	AddScalarOption("RES_SMOOTHING_COEFF", SmoothCoeff, 0.5);
	AddScalarOption("RES_SMOOTHING_ITER", nSmooth, 0);
	
	/*--- options related to the discretization ---*/
	AddEnumOption("SLOPE_LIMITER_FLOW", Kind_SlopeLimit_Flow, Limiter_Map, "NONE");
	AddEnumOption("SLOPE_LIMITER_TURB", Kind_SlopeLimit_Turb, Limiter_Map, "NONE");
	AddEnumOption("SLOPE_LIMITER_ADJTURB", Kind_SlopeLimit_AdjTurb, Limiter_Map, "NONE");
	AddEnumOption("SLOPE_LIMITER_LEVELSET", Kind_SlopeLimit_LevelSet, Limiter_Map, "NONE");
	AddScalarOption("LIMITER_COEFF", LimiterCoeff, 0.1);
	AddConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centred_Flow, Kind_Upwind_Flow);
	AddConvectOption("CONV_NUM_METHOD_COMB", Kind_ConvNumScheme_Combustion, Kind_Centred_Combustion, Kind_Upwind_Combustion);
	AddConvectOption("CONV_NUM_METHOD_LEVELSET", Kind_ConvNumScheme_LevelSet, Kind_Centred_LevelSet, Kind_Upwind_LevelSet);
	AddConvectOption("CONV_NUM_METHOD_PLASMA", Kind_ConvNumScheme_Plasma, Kind_Centred_Plasma, Kind_Upwind_Plasma);
	AddConvectOption("CONV_NUM_METHOD_ADJ", Kind_ConvNumScheme_Adj, Kind_Centred_Adj, Kind_Upwind_Adj);
	AddConvectOption("CONV_NUM_METHOD_LIN", Kind_ConvNumScheme_Lin, Kind_Centred_Lin, Kind_Upwind_Lin);
	AddConvectOption("CONV_NUM_METHOD_TURB", Kind_ConvNumScheme_Turb, Kind_Centred_Turb, Kind_Upwind_Turb);
	AddConvectOption("CONV_NUM_METHOD_ADJTURB", Kind_ConvNumScheme_AdjTurb, Kind_Centred_AdjTurb, Kind_Upwind_AdjTurb);
	AddEnumOption("NUM_METHOD_GRAD", Kind_Gradient_Method, Gradient_Map, "GREEN_GAUSS");
	AddEnumOption("VISC_NUM_METHOD_FLOW", Kind_ViscNumScheme_Flow, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_FLOW", Kind_SourNumScheme_Flow, Source_Map, "NONE");
	AddEnumOption("VISC_NUM_METHOD_ADJ", Kind_ViscNumScheme_Adj, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_ADJ", Kind_SourNumScheme_Adj, Source_Map, "NONE");
	AddEnumOption("VISC_NUM_METHOD_LIN", Kind_ViscNumScheme_Lin, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_LIN", Kind_SourNumScheme_Lin, Source_Map, "NONE");
	AddEnumOption("VISC_NUM_METHOD_PLASMA", Kind_ViscNumScheme_Plasma, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_PLASMA", Kind_SourNumScheme_Plasma, Source_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_COMB", Kind_SourNumScheme_Combustion, Source_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_TEMPLATE", Kind_SourNumScheme_Plasma, Source_Map, "NONE");
	AddEnumOption("VISC_NUM_METHOD_TURB", Kind_ViscNumScheme_Turb, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_TURB", Kind_SourNumScheme_Turb, Source_Map, "NONE");
	AddEnumOption("VISC_NUM_METHOD_ELEC", Kind_ViscNumScheme_Elec, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_ELEC", Kind_SourNumScheme_Elec, Source_Map, "NONE");
	AddEnumOption("VISC_NUM_METHOD_ADJTURB", Kind_ViscNumScheme_AdjTurb, Viscous_Map, "NONE");
	AddEnumOption("SOUR_NUM_METHOD_ADJTURB", Kind_SourNumScheme_AdjTurb, Source_Map, "NONE");
	
	/*--- options related to the artificial dissipation ---*/
	default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
	AddArrayOption("AD_COEFF_FLOW", 3, Kappa_Flow, default_vec_3d);
	default_vec_3d[0] = 0.5; default_vec_3d[1] = 0.02;
	AddArrayOption("AD_COEFF_ADJ", 2, Kappa_Adj, default_vec_3d);
	default_vec_3d[0] = 0.5; default_vec_3d[1] = 0.02;
	AddArrayOption("AD_COEFF_LIN", 2, Kappa_Lin, default_vec_3d);
	AddSpecialOption("AD_STRETCHING_FLOW", StretchingFactor_Flow, SetBoolOption, true);
	AddSpecialOption("AD_STRETCHING_ADJ", StretchingFactor_Adj, SetBoolOption, true);
	AddSpecialOption("AD_STRETCHING_LIN", StretchingFactor_Lin, SetBoolOption, true);
	
	/*--- options related to the adjoint and gradient ---*/
	AddScalarOption("ADJ_CFL_REDUCTION", Adj_CFLRedCoeff, 0.8);
	AddScalarOption("PRIMGRAD_THRESHOLD", PrimGrad_Threshold, 1E14);
	AddEnumOption("CADJ_OBJFUNC", Kind_ObjFunc, Objective_Map, "DRAG");
	AddScalarOption("SENS_LIMIT", Sens_Limit, 1000);        
	AddSpecialOption("FROZEN_VISC", Frozen_Visc, SetBoolOption, true);
	AddSpecialOption("BIGRID_FILTER", BiGrid, SetBoolOption, false);
	
	/*--- options related to input/output files and formats ---*/
	AddEnumOption("OUTPUT_FORMAT", Output_FileFormat, Output_Map, "PARAVIEW");
	AddEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, "SU2");
	AddSpecialOption("CGNS_TO_SU2", CGNS_To_SU2, SetBoolOption, false);
	AddScalarOption("NEW_SU2_FILENAME", New_SU2_FileName, string("new_mesh.su2"));
	AddScalarOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));
	AddScalarOption("MESH_OUT_FILENAME", Mesh_Out_FileName, string("mesh_out.su2"));
	AddScalarOption("CONV_FILENAME", Conv_FileName, string("history"));
	AddScalarOption("SOLUTION_FLOW_FILENAME", Solution_FlowFileName, string("solution_flow.dat"));
	AddScalarOption("SOLUTION_LIN_FILENAME", Solution_LinFileName, string("solution_lin.dat"));
	AddScalarOption("SOLUTION_ADJ_FILENAME", Solution_AdjFileName, string("solution_adj.dat"));
	AddScalarOption("RESTART_FLOW_FILENAME", ReStart_FlowFileName, string("restart_flow.dat"));
	AddScalarOption("RESTART_LIN_FILENAME",ReStart_LinFileName, string("restart_lin.dat"));
	AddScalarOption("RESTART_ADJ_FILENAME", ReStart_AdjFileName, string("restart_adj.dat"));
	AddScalarOption("VOLUME_FLOW_FILENAME", Flow_FileName, string("flow"));
	AddScalarOption("VOLUME_ADJ_FILENAME", Adj_FileName, string("adjoint"));
	AddScalarOption("VOLUME_LIN_FILENAME", Lin_FileName, string("linearized"));
	AddScalarOption("GRAD_OBJFUNC_FILENAME", ObjFunc_Grad_FileName, string("of_grad.dat"));
	AddScalarOption("SURFACE_FLOW_FILENAME", SurfFlowCoeff_FileName, string("surface_flow"));
	AddScalarOption("SURFACE_ADJ_FILENAME", SurfAdjCoeff_FileName, string("surface_adjoint"));
	AddScalarOption("SURFACE_LIN_FILENAME", SurfLinCoeff_FileName, string("surface_linear"));
	AddScalarOption("WRT_SOL_FREQ", Wrt_Sol_Freq, 1000);
	AddScalarOption("WRT_CON_FREQ",  Wrt_Con_Freq, 1);
	
	/*--- options related to the equivalent area ---*/
	AddSpecialOption("EQUIV_AREA", EquivArea, SetBoolOption, false);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 1.0; default_vec_3d[2] = 1.0;
	AddArrayOption("EA_INT_LIMIT", 3, EA_IntLimit, default_vec_3d);
	AddScalarOption("HORIZONTAL_PLANE_POSITION", Position_Plane, -1.0);
	AddScalarOption("DRAG_IN_SONICBOOM", WeightCd, 0.0);
	string_array[0] = "97"; string_array[1] = "98";
	AddArrayOption("HORIZONTAL_PLANE_MARKER", 2, PlaneTag, string_array);
	
	/*--- options related to rotational problems ---*/
	AddSpecialOption("ROTATING_FRAME", Rotating_Frame, SetBoolOption, false);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	AddArrayOption("ROTATIONAL_ORIGIN", 3, RotAxisOrigin, default_vec_3d);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	AddArrayOption("ROTATION_RATE", 3, Omega, default_vec_3d);
  AddScalarOption("ROT_VEL_REF", RotVel_Ref, 1.0);
	AddScalarOption("RED_FREC", Reduced_Frequency, 0.0);
	AddScalarOption("MAX_PITCH", Pitching_Amplitude, 0.0);        
	
	/*--- options related to reference and freestream quantities ---*/
	Length_Ref = 1.0; //<---- NOTE: this should be given an option or set as a const
	AddScalarOption("GAMMA_VALUE", Gamma, 1.4);
	AddScalarOption("GAMMA_MONATOMIC_VALUE", GammaMonatomic, 5.0/3.0);
	AddScalarOption("GAMMA_DIATOMIC_VALUE", GammaDiatomic, 7.0/5.0);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	AddArrayOption("FREESTREAM_VELOCITY", 3, Velocity_FreeStream, default_vec_3d);
	default_vec_3d[0] = 0.0; default_vec_3d[1] = 0.0; default_vec_3d[2] = 0.0;
	AddArrayOption("REF_ORIGIN_MOMENT", 3, RefOriginMoment, default_vec_3d);
	AddScalarOption("REF_AREA", RefAreaCoeff, 1.0);
	AddScalarOption("REF_LENGTH_MOMENT", RefLengthMoment, 1.0);
	AddScalarOption("SIDESLIP_ANGLE", AoS, 0.0);
	AddScalarOption("AOA", AoA, 0.0);
	AddScalarOption("GAS_CONSTANT", Gas_Constant, 287.87);
	AddScalarOption("FREESTREAM_PRESSURE", Pressure_FreeStream, 101325.0);
	AddScalarOption("FREESTREAM_DENSITY", Density_FreeStream, -1.0); 
	AddScalarOption("FREESTREAM_TEMPERATURE", Temperature_FreeStream, 273.15);
	AddScalarOption("FREESTREAM_VISCOSITY", Viscosity_FreeStream, -1.0);
	AddScalarOption("PRANDTL_LAM", Prandtl_Lam, 0.72);
	AddScalarOption("PRANDTL_TURB", Prandtl_Turb, 0.90);
	AddScalarOption("REF_PRESSURE", Pressure_Ref, 1.0);
	AddScalarOption("REF_TEMPERATURE", Temperature_Ref, 1.0);
	AddScalarOption("REF_DENSITY", Density_Ref, 1.0);
	AddScalarOption("REF_VELOCITY", Velocity_Ref, -1.0);
	AddScalarOption("REF_VISCOSITY", Viscosity_Ref, -1.0);
	AddScalarOption("MACH_NUMBER", Mach, 0.0);
	AddScalarOption("REYNOLDS_NUMBER", Reynolds, 0.0);
	AddScalarOption("REYNOLDS_LENGTH", Length_Reynolds, 1.0);
	AddScalarOption("CONVERT_TO_METER", Conversion_Factor, 1.0);
	AddScalarOption("BULK_MODULUS", Bulk_Modulus, 2.15E9);
	AddScalarOption("ARTCOMP_FACTOR", ArtComp_Factor, 1.0);
	AddScalarOption("ARTCOMP_MIN", ArtComp_Min, 0.3);
	AddScalarOption("CHARGE_COEFF", ChargeCoeff, -1.0);
	AddListOption("PARTICLE_MASS", nMass, Particle_Mass);
	
	/*--- options related to two phase simluation ---*/
	AddScalarOption("LEVELSET_ZERO", LevelSet_Zero, 1.0);
	AddScalarOption("RATIO_DENSITY", RatioDensity, 0.1);
	AddScalarOption("RATIO_VISCOSITY", RatioViscosity, 0.1);
	AddScalarOption("INTERFASE_THICKNESS", InterfaseThickness, 0.1);
	AddSpecialOption("FREESURFACE_MOVEMENT", FreeSurface_Movement, SetBoolOption, false);
	
	/*--- options related to the grid deformation ---*/
	// these options share nDV as their size in the option references; not a good idea
	AddSpecialOption("GRID_MOVEMENT", Grid_Movement, SetBoolOption, false);
	AddListOption("DV_VALUE_OLD", nDV, DV_Value_Old); 
	AddListOption("DV_VALUE_NEW", nDV, DV_Value_New);
	AddEnumListOption("DV_KIND", nDV, Design_Variable, Param_Map);
	AddDVParamOption("DV_PARAM", nDV, ParamDV, Design_Variable); 
	AddSpecialOption("GRID_DEFORM_BOX", Move_FFDBox, SetBoolOption, false);
	AddEnumOption("GRID_DEFORM_METHOD", Kind_GridDef_Method, Deform_Map, "SPRING");
	AddEnumOption("GRID_DEFORM_SOLVER", Kind_GridDef_Solver, Deform_Solver_Map, "SYM_GAUSS_SEIDEL");
	AddScalarOption("GRID_DEFORM_ERROR", GridDef_Error, 1E-14);
	
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
				//cout << option_name << ": value = "; param[option_name]->WriteValue();
			} else {
				cout << "WARNING: unrecognized option in the config. file: " << option_name << "." << endl; 
			}
		}					
	}	
	
	case_file.close();
	
	/*--- Identify crossed definitions in the config file (a parameter depend 
	 on other parameter) ---*/
	
	if (FullMG) FinestMesh = nMultiLevel;
	else FinestMesh = MESH_0;
	if (Kind_Turb_Model == SA || Kind_Turb_Model == SST) Kind_Solver = RANS;
	
	Kappa_1st_Flow = Kappa_Flow[0];
	Kappa_2nd_Flow = Kappa_Flow[1];
	Kappa_4th_Flow = Kappa_Flow[2];   
	Kappa_1st_Adj = Kappa_Adj[0];
	Kappa_4th_Adj = Kappa_Adj[1];
	Kappa_1st_Lin = Kappa_Adj[0];
	Kappa_4th_Lin = Kappa_Adj[1];

        // make the MG_PreSmooth, MG_PostSmooth, and MG_CorrecSmooth arrays consistent with nMultiLevel
        unsigned short * tmp_smooth = new unsigned short[nMultiLevel+1];
	if (nMG_PreSmooth != nMultiLevel+1) {
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

	if (nMG_PostSmooth != nMultiLevel+1) {
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

	if (nMG_CorrecSmooth != nMultiLevel+1) {
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
	MG_PreSmooth[MESH_0] = 1;
	MG_PostSmooth[MESH_0] = 0;
	MG_PostSmooth[nMultiLevel] = 0;
	MG_CorrecSmooth[nMultiLevel] = 0;
	
	if (Restart) FullMG = false;
	
	if (!ContAdj && !Linerized) PrimGrad_Threshold = 1E14;
	
	if (ContAdj) {
		if (Kind_Solver == EULER) Kind_Solver = ADJ_EULER;
		if (Kind_Solver == NAVIER_STOKES) Kind_Solver = ADJ_NAVIER_STOKES;
		if ((Kind_Solver == RANS) && (Frozen_Visc)) Kind_Solver = ADJ_NAVIER_STOKES;
		if ((Kind_Solver == RANS) && (!Frozen_Visc)) Kind_Solver = ADJ_RANS;
		if (Kind_Solver == ELECTRIC_POTENTIAL) Kind_Solver = ADJ_ELECTRIC_POTENTIAL;
	}
	if (OneShot) {
		if (Kind_Solver == EULER) Kind_Solver = ONE_SHOT_ADJ_EULER;
		if (Kind_Solver == NAVIER_STOKES) Kind_Solver = ONE_SHOT_ADJ_NAVIER_STOKES;
	}
	if (Linerized) {
		if (Kind_Solver == EULER) Kind_Solver = LIN_EULER;
		if (Kind_Solver == ELECTRIC_POTENTIAL) Kind_Solver = LIN_ELECTRIC_POTENTIAL;
	}
	
	Kind_Solver_MPI = Kind_Solver;
	
	if (Unsteady_Simulation == TIME_STEPPING) {
		nMultiLevel = 0;
		MGCycle = 0;
	}
	
	if (Rotating_Frame == YES)
		Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
	
	if (Unsteady_Simulation == DUAL_TIME_STEPPING) {
		Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
		nIntIter = nExtIter;
		Wrt_Sol_Freq = 1;
		Wrt_Con_Freq = 1;
	}
	
	nCFL = nMultiLevel+1; 
	CFL = new double[nCFL];
	CFL[0] = CFLFineGrid;
	if (ContAdj) CFL[0] = CFL[0] * Adj_CFLRedCoeff;
	for (unsigned short iCFL = 1; iCFL < nCFL; iCFL++) 
		CFL[iCFL] = CFL[iCFL-1]*MG_CFLRedCoeff;
	
	if (nRKStep == 0) {
		RK_Alpha_Step = new double[1]; RK_Alpha_Step[0] = 1.0;
		RK_Beta_Step = new double[1]; RK_Beta_Step[0] = 1.0;
	}
	
	
	if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) && (Kind_ViscNumScheme_Flow == NONE)) {
		cout << "You must define a viscous numerical method for the flow equations!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}
	
	/*	if ((Kind_Solver == MULTI_SPECIES_NAVIER_STOKES) && (Kind_ViscNumScheme_Plasma == NONE)) {
	 cout << "You must define a viscous numerical method for the plasma equations!!" << endl;
	 cout << "Press any key to exit..." << endl;
	 cin.get();
	 exit(1);
	 }
	 */
	
	if (((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) && (Kind_ViscNumScheme_Adj == NONE)) {
		cout << "You must define a viscous numerical method for the adjoint equations!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}
	
	/*--- Set the number of domains in the MPI computation ---*/
#ifndef NO_MPI
	if (val_software == SU2_CFD)
		nDomain = MPI::COMM_WORLD.Get_size()-1;
#endif
	
	/*--- The multigrid do not work with periodic bc ---*/
#ifndef NO_MPI
	/*	if ((nDomain != 0) && (nMarker_PerBound != 0) && (nMultiLevel != 0)) {
	 cout << "The multigrid do not work with periodic BC" << endl;
	 cout << "Stopping the computation" << endl;
	 MPI::COMM_WORLD.Abort(1);
	 MPI::Finalize();
	 }*/
#endif
	
	/*--- Boundary (marker) treatment ---*/
	nMarker_All = nMarker_Euler + nMarker_FarField + nMarker_SymWall + nMarker_PerBound + nMarker_NearFieldBound 
	+ nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Inlet + nMarker_Outlet + nMarker_NS + nMarker_Custom 
	+ 2*nDomain;
	
	Marker_All_Tag = new string[nMarker_All+2];									// Store the tag that correspond with each marker.
	Marker_All_SendRecv = new short[nMarker_All+2];							// +#domain (send), -#domain (receive) o 0 (no send neither receive).
	Marker_All_Boundary = new unsigned short[nMarker_All+2];		// Store the kind of boundary condition
	Marker_All_Monitoring = new unsigned short[nMarker_All+2];	// Store if the boundary should be monitored.
	Marker_All_Plotting = new unsigned short[nMarker_All+1];		// Store if the boundary should be plotted.
	Marker_All_Moving = new unsigned short[nMarker_All+1];			// Store if the boundary should be moved.
	Marker_All_PerBound = new short[nMarker_All+2];							// Store if the boundary belong to a periodic boundary.
	
	unsigned short iMarker_All, iMarker_Config, iMarker_Euler, iMarker_NS, iMarker_Custom, iMarker_FarField, 
	iMarker_SymWall, iMarker_PerBound, iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Dirichlet,
	iMarker_Inlet, iMarker_Outlet, iMarker_Monitoring, iMarker_Plotting, iMarker_Moving;
	
	for (iMarker_All = 0; iMarker_All < nMarker_All; iMarker_All++) {
		Marker_All_Tag[iMarker_All] = "NONE";
		Marker_All_SendRecv[iMarker_All] = 0;
		Marker_All_Boundary[iMarker_All] = 0;
		Marker_All_Monitoring[iMarker_All] = 0;
		Marker_All_Plotting[iMarker_All] = 0;
		Marker_All_Moving[iMarker_All] = 0;
		Marker_All_PerBound[iMarker_All] = 0;
	}
	
	nMarker_Config = nMarker_Euler + nMarker_FarField + nMarker_SymWall + nMarker_PerBound + nMarker_NearFieldBound 
	+ nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Inlet + nMarker_Outlet + nMarker_NS + nMarker_Custom;
	
	Marker_Config_Tag = new string[nMarker_Config];
	Marker_Config_Boundary = new unsigned short[nMarker_Config];
	Marker_Config_Monitoring = new unsigned short[nMarker_Config];
	Marker_Config_Plotting = new unsigned short[nMarker_Config];
	Marker_Config_Moving = new unsigned short[nMarker_Config];
	Marker_Config_PerBound = new unsigned short[nMarker_Config];
	
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Tag[iMarker_Config] = "NONE";
		Marker_Config_Boundary[iMarker_Config] = 0;
		Marker_Config_Monitoring[iMarker_Config] = 0;
		Marker_Config_Plotting[iMarker_Config] = 0;
		Marker_Config_Moving[iMarker_Config] = 0;
		Marker_Config_PerBound[iMarker_Config] = 0;
	}
	
	iMarker_Config = 0;
	for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Euler[iMarker_Euler];
		Marker_Config_Boundary[iMarker_Config] = EULER_WALL;
		iMarker_Config++;
	}
	
	for (iMarker_NS = 0; iMarker_NS < nMarker_NS; iMarker_NS++) {
		Marker_Config_Tag[iMarker_Config] = Marker_NS[iMarker_NS];
		Marker_Config_Boundary[iMarker_Config] = NO_SLIP_WALL;
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
	
	for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Outlet[iMarker_Outlet];
		Marker_Config_Boundary[iMarker_Config] = OUTLET_FLOW;
		iMarker_Config++;
	}
	
	for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
		Marker_Config_Tag[iMarker_Config] = Marker_Custom[iMarker_Custom];
		Marker_Config_Boundary[iMarker_Config] = CUSTOM_BOUNDARY;
		iMarker_Config++;
	}
	
	for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
		Marker_Config_Monitoring[iMarker_Config] = NO;
		for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++)
			if (Marker_Config_Tag[iMarker_Config] == Marker_Monitoring[iMarker_Monitoring])
				Marker_Config_Monitoring[iMarker_Config] = YES;
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
	
	/*--- Begin screen output ---*/
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
	
	/*--- The screen output is done only for the master node ---*/
	if (rank == MASTER_NODE) {
		
		cout << endl <<"-------------------------------------------------------------------------" << endl;
		switch (val_software) {
			case SU2_CFD: cout << "|             SU2 Suite (Computational Fluid Dynamics Code)             |" << endl; break;
			case SU2_MDC: cout << "|                   SU2 Suite (Mesh Deformation Code)                   |" << endl; break;
			case SU2_GPC: cout << "|                 SU2 Suite (Gradient Projection Code)                  |" << endl; break;
			case SU2_DDC: cout << "|                SU2 Suite (Domain Decomposition Code)                  |" << endl; break;
			case SU2_MAC: cout << "|                   SU2 Suite (Mesh Adaptation Code)                    |" << endl; break;
			case SU2_GDC: cout << "|                   SU2 Suite (Geometry Design Code)                    |" << endl; break;
			case SU2_PBC: cout << "|                  SU2 Suite (Periodic Boundary Code)                   |" << endl; break;
		}
		cout <<"-------------------------------------------------------------------------" << endl;
		
		
		cout << endl <<"------------------------ Physical case definition -----------------------" << endl;
		if (val_software == SU2_CFD) {
			switch (Kind_Solver) {
				case POTENTIAL_FLOW: cout << "Potential flow equation." << endl; break;
				case EULER:
					if (Incompressible) cout << "Incompressible Euler equations (only for 2D flows!!!)." << endl;
					else cout << "Compressible Euler equations." << endl; break;
				case NAVIER_STOKES: cout << "Laminar Navier-Stokes' equations." << endl; break;
				case RANS: cout << "RANS equations with:"; break;
				case MULTI_SPECIES_NAVIER_STOKES: 
					cout << "Plasma equations (flow + plasma + electrical potential)." << endl; 
					if (Kind_GasModel == ARGON)
						cout << "Using 3 species Argon gas model." << endl;
					if (Kind_GasModel == AIR7)
						cout << "Using 7 species Air gas model." << endl;
					break;
				case ELECTRIC_POTENTIAL: cout << "Electric potential equation." << endl; break;
				case ADJ_EULER: cout << "Continuous Euler adjoint equations." << endl; break;
				case ADJ_NAVIER_STOKES: cout << "Continuous Navier-Stokes adjoint equations with frozen viscosity." << endl; break;
				case ADJ_RANS: cout << "Continuous RANS adjoint equations." << endl; break;
				case ADJ_ELECTRIC_POTENTIAL: cout << "Continuous electric potential adjoint equation." << endl; break;
				case ONE_SHOT_ADJ_EULER: cout << "One shot method for the Euler equations." << endl; break;
				case ONE_SHOT_ADJ_NAVIER_STOKES: cout << "One shot method for the Navier-Stokes equations." << endl; break;
				case LIN_EULER: cout << "Linearized Euler equations." << endl; break;
				case LIN_ELECTRIC_POTENTIAL: cout << "Linearized electric potential equation." << endl; break;
				case TWO_PHASE_FLOW:
					cout << "Two-Phase flow equation. Density ratio: " << RatioDensity <<". Viscosity ratio: "<< RatioViscosity << "." << endl;
					break;
			}
			if (GravityForce) cout << "Include gravity force in the formulation (only for incompressible flows!!!)." << endl;
			
			if (Kind_Solver == RANS){
				switch (Kind_Turb_Model){
					case NONE: cout << " No turbulence model selected." << endl; break;
					case SA:   cout << " Spalart-Allmaras turbulence model." << endl; break;
					case SST:  cout << " Menter SST turbulence model." << endl; break;
				}
			}
			if (!Incompressible) {
				cout << "Mach number: " << Mach <<"."<< endl;
				cout << "Angle of attack (AoA): " << AoA <<" deg, and angle of sideslip (AoS): " << AoS <<" deg."<< endl;
				if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
						(Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS) || 
						(Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES))
					cout << "Reynolds number: " << Reynolds <<"."<< endl;
			}
			if (EquivArea) {
				cout <<"The equivalent area is going to be evaluated on the near-field."<< endl;
				cout <<"The lower integration limit is "<<EA_IntLimit[0]<<", and the upper is "<<EA_IntLimit[1]<<"."<<endl; 
				cout <<"The near-field is situated at "<<EA_IntLimit[2]<<"."<<endl; 
			}
		}
		
		cout << "Input mesh file name: " << Mesh_FileName << endl;
		
		if (Divide_Element) cout << "Divide elements into triangles and tetrahedrons." << endl;
		
		if (val_software == SU2_GPC) {
			cout << "Input sensitivity file name: " << SurfAdjCoeff_FileName << "." << endl;
		}
		
		if (val_software == SU2_CFD) {
			switch (Grid_Movement) {
				case YES:
					cout << "Grid movement (fluttering)." << endl;
					cout << "Reduced frequency: "<< Reduced_Frequency << "." << endl;
					cout << "Maximum amplitude of pitching: "<< Pitching_Amplitude << "." << endl;
					break;
			}
			
			switch (Rotating_Frame) {
			case YES:
					cout << "Performing a simulation in a rotating frame." << endl;
					cout << "Origin of the rotation axis: " << RotAxisOrigin[0] <<","<< RotAxisOrigin[1] <<","<< RotAxisOrigin[2] << ". ";
					cout << "Angular velocity vector: " << Omega[0] <<","<< Omega[1] <<","<< Omega[2] << "." << endl;
					break;
			}
			
			switch (Restart) {
				case YES:
					if (!ContAdj && !Linerized) cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
					if (ContAdj) cout << "Read adjoint solution from: " << Solution_AdjFileName << "." << endl;
					if (Linerized) cout << "Read linearized solution from: " << Solution_LinFileName << "." << endl;
					
					break;
				case NO: cout << "No restart solution, use the values at infinity (freestream)." << endl; break;
			}
			
			if (ContAdj || Linerized)
				cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl;
			
			cout << "Surface(s) where the force coefficients are to be evaluated: ";
			for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
				cout << Marker_Monitoring[iMarker_Monitoring];
				if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
				else cout <<"."<<endl;
			}
			
			if (RefAreaCoeff == 0) cout << "The reference length/area will be computed using y(2D) or z(3D) projection." <<endl;
			else cout << "The reference length/area (force coefficient) is " << RefAreaCoeff << "." <<endl;
			cout << "The reference length (moment computation) is " << RefLengthMoment << "." <<endl;
			cout << "Reference origin (moment computation) is (" <<RefOriginMoment[0]<<", "<<RefOriginMoment[1]<<", "<<RefOriginMoment[1]<<")."<< endl;		
			
		}
		
		if (val_software == SU2_MAC) {
			switch (Kind_Adaptation) {
				case FULL: case WAKE: case TWOPHASE: case FULL_FLOW: case FULL_ADJOINT: case FULL_LINEAR: case HORIZONTAL_PLANE: case SMOOTHING: case SUPERSONIC_SHOCK:
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
				case TORSIONAL_SPRING: cout << "Grid deformation using a torsional spring method." << endl; break;
				case ALGEBRAIC: cout << "Grid deformation using an algebraic method." << endl; break;
			}
			switch (Kind_GridDef_Solver) {
				case SYM_GAUSS_SEIDEL: cout << "A symmetric Gauss-Seidel method is used for solving the linear system." << endl; break;
				case CONJUGATE_GRADIENT: cout << "A conjugate gradient method is used for solving the linear system." << endl; break;
				case PREC_CONJUGATE_GRADIENT: cout << "A precond. conjugate gradient is used for solving the linear system." << endl; break;
			}
			cout << "Convergence criteria of the linear solver: "<<GridDef_Error<<"."<<endl;			
			
			if (Design_Variable[0] != NO_DEFORMATION) {
				/*				if ((Design_Variable == FFD_CONTROL_POINT) || (Design_Variable == FFD_DIHEDRAL_ANGLE) ||
				 (Design_Variable == FFD_TWIST_ANGLE) || (Design_Variable == FFD_ROTATION)) {
				 if (Move_FFDBox == YES) cout << "The 3D volumetric grid movement is limited to the FFD Box." <<endl;
				 else cout << "The 3D volumetric grid movement is not limited." <<endl;
				 }*/
				cout << "Geo. design var. definition (markers <-> old def., new def. <-> param):" <<endl;
				for (unsigned short iDV = 0; iDV < nDV; iDV++) {
					switch (Design_Variable[iDV]) {
						case NO_DEFORMATION: cout << "There isn't any deformation." ; break;
						case HICKS_HENNE: cout << "Hicks Henne <-> " ; break;
						case MACH_NUMBER: cout << "Mach number <-> " ; break;
						case DISPLACEMENT: cout << "Displacement design variable."; break;
						case NACA_4DIGITS: cout << "NACA four digits <-> "; break;
						case PARABOLIC: cout << "Parabolic <-> "; break;
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
					if (Design_Variable[iDV] == DISPLACEMENT) nParamDV = 3;
					if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
					if (Design_Variable[iDV] == NACA_4DIGITS) nParamDV = 3;
					if (Design_Variable[iDV] == PARABOLIC) nParamDV = 2;
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
		
		if (((val_software == SU2_CFD) && ( Linerized )) || (val_software == SU2_GPC)) {
			cout << endl <<"-------------------- Surface deformation parameters ---------------------" << endl;
			switch (Design_Variable[0]) {
				case NO_DEFORMATION: cout << "No deformation." << endl; break;
				case HICKS_HENNE: cout << "Hicks Henne design variables." << endl; break;
				case HICKS_HENNE_SHOCK: cout << "Shock Hicks Henne design variables." << endl; break;
				case HICKS_HENNE_NORMAL: cout << "Shock Hicks Henne design variables." << endl; break;
				case MACH_NUMBER: cout << "Mach number design variable." << endl; break;
				case DISPLACEMENT: cout << "Displacement design variable."<< endl; break;
				case ROTATION: cout << "Rotation around an axis design variable."<< endl; break;
				case FFD_CONTROL_POINT: cout << "Free form deformation (control point change)."<< endl; break;
				case FFD_DIHEDRAL_ANGLE: cout << "Free form deformation (dihedral angle change)."<< endl; break;
				case FFD_TWIST_ANGLE: cout << "Free form deformation (twist angle change)."<< endl; break;
				case FFD_ROTATION: cout << "Free form deformation (rotation change)."<< endl; break;
			}
			if (Design_Variable[0] != NO_DEFORMATION) {
				cout << "Previous value of the design variable with respect to the non-deformed geometry: "<< DV_Value_Old[0] << "." << endl;
				cout << "New value of the design variable: "<< DV_Value_New[0] << "." << endl;
				cout << "Parameters of the design variable: ";
				for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {
					cout << ParamDV[0][iParamDV];
					if (iParamDV < nParamDV-1) cout << ", ";
					else cout <<"."<<endl;
				} 
			}
			cout << "Surface(s) to be deformed: ";
			for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
				cout << Marker_Moving[iMarker_Moving];
				if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
				else cout <<"."<<endl;
			}	
		}
		
		if (((val_software == SU2_CFD) && ( ContAdj || OneShot )) || (val_software == SU2_GPC)) {
			
			cout << endl <<"----------------------- Design problem definition -----------------------" << endl;
			switch (Kind_ObjFunc) {
				case DRAG_COEFFICIENT: cout << "Drag objective function." << endl; break;
				case LIFT_COEFFICIENT: cout << "Lift objective function." << endl; break;
				case SIDEFORCE_COEFFICIENT: cout << "Side force objective function." << endl; break;
				case MOMENT_X_COEFFICIENT: cout << "Pitching moment objective function." << endl; break;
				case MOMENT_Y_COEFFICIENT: cout << "Rolling moment objective function." << endl; break;
				case MOMENT_Z_COEFFICIENT: cout << "Yawing moment objective function." << endl; break;
				case EFFICIENCY: cout << "Efficiency objective function." << endl; break;
				case PRESSURE_COEFFICIENT: cout << "Pressure objective function." << endl; break;
				case ELECTRIC_CHARGE: cout << "Electric charge objective function." << endl; break;
				case EQUIVALENT_AREA:
					cout << "Equivalent area objective function." << endl;
					cout << "Drag coefficient weight in the objective function: " << WeightCd <<"."<< endl;  break;
				case NEARFIELD_PRESSURE:
					cout << "Nearfield pressure objective function." << endl;
					cout << "Drag coefficient weight in the objective function: " << WeightCd <<"."<< endl;  break;
					
					break;
			}
			
			cout << "Sensitivity limit: "<< Sens_Limit << "." << endl;
			cout << "Primitive variables gradient threshold: "<< PrimGrad_Threshold << "." << endl;
			
		}
		
		if (val_software == SU2_CFD) {
			cout << endl <<"---------------------- Space numerical integration ----------------------" << endl;
			
			if (SmoothNumGrid) cout << "There are some smoothing iterations on the grid coordinates." <<endl;
			
			if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) || (Kind_Solver == ONE_SHOT_ADJ_EULER) 
					|| (Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES)) {
				if ((Kind_ConvNumScheme_Flow == SPACE_CENTRED) && (Kind_Centred_Flow == JST)) {
					cout << "Jameson-Schmidt-Turkel scheme for the flow inviscid terms."<< endl;
					cout << "JST viscous coefficients (1st, 2nd & 4th): " << Kappa_1st_Flow 
					<< ", " << Kappa_2nd_Flow << ", " << Kappa_4th_Flow <<"."<< endl;
					if (StretchingFactor_Flow) cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
					else cout << "The method doesn't include a grid stretching correction."<< endl;
				}
				if ((Kind_ConvNumScheme_Flow == SPACE_CENTRED) && (Kind_Centred_Flow == LAX)) 
					cout << "Lax-Friedrich scheme for the flow inviscid terms."<< endl;
				if (Kind_ConvNumScheme_Flow == SPACE_UPWIND) {
					if (Kind_Upwind_Flow == ROE_1ST) cout << "1st order Roe solver for the flow inviscid terms."<< endl;
					if (Kind_Upwind_Flow == AUSM_1ST)	cout << "1st order AUSM solver for the flow inviscid terms."<< endl;
					if (Kind_Upwind_Flow == HLLC_1ST)	cout << "1st order HLLC solver for the flow inviscid terms."<< endl;
				}
				if ((Kind_ConvNumScheme_Flow == SPACE_UPWIND) && 
						((Kind_Upwind_Flow == ROE_2ND) || (Kind_Upwind_Flow == AUSM_2ND) || (Kind_Upwind_Flow == HLLC_2ND))) {
					if (Kind_Upwind_Flow == ROE_2ND) cout << "2nd order Roe solver for the flow inviscid terms."<< endl;
					if (Kind_Upwind_Flow == AUSM_2ND) cout << "2nd order AUSM solver for the flow inviscid terms."<< endl;
					if (Kind_Upwind_Flow == HLLC_2ND) cout << "2nd order HLLC solver for the flow inviscid terms."<< endl;
					switch (Kind_SlopeLimit_Flow) {
						case NONE: cout << "Without slope-limiting method." << endl; break;
						case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl; break;
					}
				}
			}
			
			if (Kind_Solver == MULTI_SPECIES_NAVIER_STOKES) {
				if ((Kind_ConvNumScheme_Plasma == SPACE_CENTRED) && (Kind_Centred_Plasma == JST)) {
					cout << "Jameson-Schmidt-Turkel scheme for the flow inviscid terms."<< endl;
					//					cout << "JST Viscous coefficients (1st, 2nd & 4th): " << Kappa_1st_Flow 
					//					<< ", " << Kappa_2nd_Flow << ", " << Kappa_4th_Flow <<"."<< endl;
					//					if (StretchingFactor_Flow) cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
					//					else cout << "The method doesn't include a grid stretching correction."<< endl;
				}
				if ((Kind_ConvNumScheme_Plasma == SPACE_CENTRED) && (Kind_Centred_Plasma == LAX)) 
					cout << "Lax-Friedrich scheme for the plasma inviscid terms."<< endl;
				if ((Kind_ConvNumScheme_Plasma == SPACE_UPWIND) && (Kind_Upwind_Plasma == ROE_1ST)) 
					cout << "1st order Roe solver for the plasma inviscid terms."<< endl;
				if ((Kind_ConvNumScheme_Flow == SPACE_UPWIND) && (Kind_Upwind_Plasma == ROE_2ND)) {
					cout << "2nd order Roe solver for the plasma inviscid terms."<< endl;
					switch (Kind_SlopeLimit_Plasma) {
						case NONE: cout << "Without slope-limiting method." << endl; break;
						case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method." << endl; break;
					}
				}
			}
			
			
			if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) 
					|| (Kind_Solver == ONE_SHOT_ADJ_EULER) || (Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES)) {
				if ((Kind_ConvNumScheme_Adj == SPACE_CENTRED) && (Kind_Centred_Adj == JST)) {
					cout << "Jameson-Schmidt-Turkel scheme for the adjoint inviscid terms."<< endl;
					cout << "JST viscous coefficients (1st, & 4th): " << Kappa_1st_Adj 
					<< ", " << Kappa_4th_Adj <<"."<< endl;
					if (StretchingFactor_Adj) cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
					else cout << "The method doesn't include a grid stretching correction."<< endl;
				}
				if ((Kind_ConvNumScheme_Adj == SPACE_CENTRED) && (Kind_Centred_Adj == LAX)) 
					cout << "Lax-Friedrich scheme for the adjoint inviscid terms."<< endl;
				if ((Kind_ConvNumScheme_Adj == SPACE_UPWIND) && (Kind_Upwind_Adj == ROE_1ST)) 
					cout << "1st order Roe solver for the adjoint inviscid terms."<< endl;
				if ((Kind_ConvNumScheme_Adj == SPACE_UPWIND) && (Kind_Upwind_Adj == ROE_2ND)) {
					cout << "2nd order Roe solver for the adjoint inviscid terms."<< endl;
					switch (Kind_SlopeLimit_Flow) {
						case NONE: cout << "Without slope-limiting method." << endl; break;
						case VENKATAKRISHNAN: cout << "Venkatakrishnan slope-limiting method." << endl; break;
					}			
				}
			}
			
			if (Kind_Solver == LIN_EULER) {
				if ((Kind_ConvNumScheme_Lin == SPACE_CENTRED) && (Kind_Centred_Lin == JST)) {
					cout << "Jameson-Schmidt-Turkel scheme for the linearized inviscid terms."<< endl;
					cout << "JST viscous coefficients (1st, & 4th): " << Kappa_1st_Lin
					<< ", " << Kappa_4th_Lin <<"."<< endl;
					if (StretchingFactor_Lin) cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
					else cout << "The method doesn't include a grid stretching correction."<< endl;
				}
				if ((Kind_ConvNumScheme_Lin == SPACE_CENTRED) && (Kind_Centred_Lin == LAX)) 
					cout << "Lax-Friedrich scheme for the linearized inviscid terms."<< endl;
			}
			
			if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) || (Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES) || 
					(Kind_Solver == MULTI_SPECIES_NAVIER_STOKES)) {
				switch (Kind_ViscNumScheme_Flow) {
					case AVG_GRAD: cout << "Average of gradients (1st order) for computation of viscous flow terms." << endl; break;
					case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (2nd order) for computation of viscous flow terms." << endl; break;
				}
			}
			
			if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
				if (Kind_SourNumScheme_Flow == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the flow source terms." << endl;
			}
			
			if (Kind_Solver == MULTI_SPECIES_NAVIER_STOKES) {
				if (Kind_SourNumScheme_Flow == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the plasma source terms." << endl;
			}
			
			if (Kind_Solver == RANS) {
				if ((Kind_ConvNumScheme_Turb == SPACE_UPWIND) && (Kind_Upwind_Turb == SCALAR_UPWIND_1ST)) 
					cout << "Scalar upwind solver (first order) for the turbulence model."<< endl;
				if ((Kind_ConvNumScheme_Turb == SPACE_UPWIND) && (Kind_Upwind_Turb == SCALAR_UPWIND_2ND)) 
					cout << "Scalar upwind solver (second order) for the turbulence model."<< endl;
			}
			
			if (Kind_Solver == ADJ_RANS) {
				if ((Kind_ConvNumScheme_AdjTurb == SPACE_UPWIND) && (Kind_Upwind_AdjTurb == SCALAR_UPWIND_1ST)) 
					cout << "Adjoint turbulent eq - Scalar upwind solver (first order)"<< endl;
				if ((Kind_ConvNumScheme_AdjTurb == SPACE_UPWIND) && (Kind_Upwind_AdjTurb == SCALAR_UPWIND_2ND)) 
					cout << "Adjoint turbulent eq - Scalar upwind solver (second order)"<< endl;
			}
			
			if ((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES)) {	
				switch (Kind_ViscNumScheme_Adj) {
					case AVG_GRAD: cout << "Average of gradients (1st order) for computation of viscous adjoint terms." << endl; break;
					case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (2nd order) for computation of viscous adjoint terms." << endl; break;
				}
			}
			
			if (Kind_Solver == RANS) {
				if (Kind_ViscNumScheme_Turb == AVG_GRAD) cout << "Average of gradients (1st order) for computation of viscous turbulence terms." << endl;
				if (Kind_ViscNumScheme_Turb == AVG_GRAD_CORRECTED) cout << "Average of gradients with correction (2nd order) for computation of viscous turbulence terms." << endl;
				if (Kind_SourNumScheme_Turb == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the turbulence model source terms." << endl;
			}
			
			if ((Kind_Solver == ELECTRIC_POTENTIAL) || (Kind_Solver == MULTI_SPECIES_NAVIER_STOKES)) {
				if (Kind_ViscNumScheme_Elec == GALERKIN) cout << "Galerkin method for viscous terms computation of the electric potential equation." << endl;
			}
			
			if (Kind_Solver == ADJ_RANS) {
				if (Kind_ViscNumScheme_AdjTurb == AVG_GRAD) cout << "Average of gradients (1st order) for computation of adjoint viscous turbulence terms." << endl;
				if (Kind_ViscNumScheme_AdjTurb == AVG_GRAD_CORRECTED) cout << "Average of gradients with correction (2nd order) for computation of adjoint viscous turbulence terms." << endl;
				if (Kind_SourNumScheme_AdjTurb == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the turbulence adjoint model source terms." << endl;
				if (Kind_TimeIntScheme_AdjTurb == EULER_IMPLICIT) cout << "Euler implicit method for the turbulent adjoint equation." << endl;
			}
			
			if ((Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
				if (Kind_SourNumScheme_Adj == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the Navier-Stokes eq. source terms." << endl;
			}
			
			if ((Kind_Solver == ELECTRIC_POTENTIAL) || (Kind_Solver == MULTI_SPECIES_NAVIER_STOKES)) {
				if (Kind_SourNumScheme_Elec == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the electric potential source terms." << endl;
			}
			
			switch (Kind_Gradient_Method) {
				case GREEN_GAUSS: cout << "Gradient computation by Green-Gauss theorem." << endl; break;
				case LEAST_SQUARES: cout << "Gradient Computation by unweighted Least-Squares method." << endl; break;
				case WEIGHTED_LEAST_SQUARES: cout << "Gradient Computation by weighted Least-Squares method." << endl; break;
			}
			
			if (Incompressible) {
				cout << "Artificial compressibility factor: " << ArtComp_Factor << "." <<endl;
				cout << "Min value for the artificial compressibility: " << ArtComp_Factor*ArtComp_Min << "." <<endl;
			}
			
			cout << endl <<"---------------------- Time numerical integration -----------------------" << endl;
			switch (Unsteady_Simulation) {
				case NO: cout << "Local time stepping (steady state simulation)." << endl; break;
				case TIME_STEPPING:
					cout << "Unsteady simulation using a time stepping strategy. ";
					if (Unst_TimeStep == 0.0) cout << "Time step computed by the code." << endl;
					else cout << "Time step: "<< Unst_TimeStep << "." << endl;
					break;
				case DUAL_TIME_STEPPING:
					cout << "Unsteady simulation using a dual time stepping strategy. ";
					if (Unst_TimeStep == 0.0) cout << "Time step computed by the code." << endl;
					else cout << "Time step: "<< Unst_TimeStep << "." << endl;
					break;
			}
			
			if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) || (Kind_Solver == ONE_SHOT_ADJ_EULER) 
					|| (Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES)) {
				switch (Kind_TimeIntScheme_Flow) {
					case RUNGE_KUTTA_EXPLICIT:
						cout << "Runge-Kutta explicit method for the flow equations." << endl;
						cout << "Number of steps: " << nRKStep << endl;
						cout << "Alpha coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Alpha_Step[iRKStep];
						}
						cout << endl;
						cout << "Beta coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Beta_Step[iRKStep];
						}
						cout << endl;
						break;
					case EULER_EXPLICIT: cout << "Euler explicit method for the flow equations." << endl; break;
					case EULER_IMPLICIT: cout << "Euler implicit method for the flow equations." << endl; break;
				}
			}
			
			if (Kind_Solver == MULTI_SPECIES_NAVIER_STOKES) {
				switch (Kind_TimeIntScheme_Plasma) {
					case RUNGE_KUTTA_EXPLICIT:
						cout << "Runge-Kutta explicit method for the plasma equations." << endl;
						cout << "Number of steps: " << nRKStep << endl;
						cout << "Alpha coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Alpha_Step[iRKStep];
						}
						cout << endl;
						cout << "Beta coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Beta_Step[iRKStep];
						}
						cout << endl;
						break;
					case EULER_EXPLICIT: cout << "Euler explicit method for the plasma equations." << endl; break;
					case EULER_IMPLICIT: cout << "Euler implicit method for the plasma equations." << endl; break;
				}
			}
			
			if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) 
					|| (Kind_Solver == ONE_SHOT_ADJ_EULER) || (Kind_Solver == ONE_SHOT_ADJ_NAVIER_STOKES)) {
				switch (Kind_TimeIntScheme_Adj) {
					case RUNGE_KUTTA_EXPLICIT:
						cout << "Runge-Kutta explicit method for the adjoint equations." << endl;
						cout << "Number of steps: " << nRKStep << endl;
						cout << "Alpha coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Alpha_Step[iRKStep];
						}
						cout << endl;
						cout << "Beta coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Beta_Step[iRKStep];
						}
						cout << endl;
						break;
					case EULER_EXPLICIT: cout << "Euler explicit method for the adjoint equations." << endl; break;
					case EULER_IMPLICIT: cout << "Euler implicit method for the adjoint equations." << endl; break;
				}
			}
			
			if (Kind_Solver == LIN_EULER) {
				switch (Kind_TimeIntScheme_Lin) {
					case RUNGE_KUTTA_EXPLICIT:
						cout << "Runge-Kutta explicit explicit method for the linearized equations." << endl;
						cout << "Number of steps: " << nRKStep << endl;
						cout << "Alpha coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Alpha_Step[iRKStep];
						}
						cout << endl;
						cout << "Beta coefficients: ";
						for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
							cout << "\t" << RK_Beta_Step[iRKStep];
						}
						cout << endl;
						break;
					case EULER_EXPLICIT: cout << "Euler explicit method for the linearized equations." << endl; break;
				}
			}
			
			switch (FullMG) {
				case YES: cout << "Full Multigrid." << endl; break;
			}
			
			if (nMultiLevel !=0) {
				if (nStartUpIter != 0) cout << "A total of " << nStartUpIter << " start up iterations on the fine grid."<< endl;
				if (MGCycle == 0) cout << "V Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;
				if (MGCycle == 1) cout << "W Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;	
				
				cout << "Reduction of the CFL coefficient in the coarse levels: " << MG_CFLRedCoeff <<"."<<endl;
				cout << "Max. number of children in the agglomeration stage: " << MaxChildren <<"."<<endl;
				cout << "Max. length of an agglom. elem. (compared with the domain): " << MaxDimension <<"."<<endl;
				if (ReStartMGCycle) cout <<"The coarse level solution is restarted at each MG cycle."<<endl;
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
			
			if (nSmooth != 0) cout << "Residual smoothing strategy with "<<nSmooth<<
				" iterations, and a factor of "<< SmoothCoeff <<"."<<endl;
			
			if (Kind_Solver == RANS)
				if (Kind_TimeIntScheme_Turb == EULER_IMPLICIT) 
					cout << "Euler implicit time integration for the turbulence model." << endl;
		}
		
		if (val_software == SU2_CFD)
			cout << endl <<"------------------------- Convergence criteria --------------------------" << endl;
		
		if (val_software == SU2_CFD) {
			cout << "Maximum number of iterations: " << nExtIter <<"."<<endl;
			
			if (ConvCriteria == CAUCHY) {
				if (!ContAdj && !Linerized)
					switch (Cauchy_Func_Flow) {
						case LIFT_COEFFICIENT: cout << "Cauchy criteria for Lift using "
							<< Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
						case DRAG_COEFFICIENT: cout << "Cauchy criteria for Drag using "
							<< Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
					}
				
				if (ContAdj)
					switch (Cauchy_Func_AdjFlow) {
						case SENS_GEOMETRY: cout << "Cauchy criteria for geo. sensitivity using "
							<< Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
						case SENS_MACH: cout << "Cauchy criteria for Mach number sensitivity using "
							<< Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
					}
				
				if (Linerized)
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
				if (!ContAdj && !Linerized) {
					cout << "Reduce the density residual " << OrderMagResidual << " orders of magnitude."<< endl;
					cout << "The minimum bound for the density residual is 10^(" << MinLogResidual<< ")."<< endl;
					cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
				}
				
				if (ContAdj) {
					cout << "Reduce the adjoint density residual " << OrderMagResidual << " orders of magnitude."<< endl;
					cout << "The minimum value for the adjoint density residual is 10^(" << MinLogResidual<< ")."<< endl;
				}
				
				if (Linerized) {
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
				case TWOPHASE: cout << "Grid adaptation of the interphase of a two phase flow." << endl; break;
				case HORIZONTAL_PLANE:
					cout << "Grid adaptation of a horizontal plane." << endl;
					cout << "The horizontal plane is situated at (y 2D, or z 3D): " << Position_Plane <<"."<< endl;
					cout << "Markers of the new horizontal plane: " << PlaneTag[0] <<" (upper), and "<< PlaneTag[1] <<" (lower)."<< endl;
					break;
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
				case CURVE_SURFACE: cout << "Create a curve surface using the current structured grid." << endl; break;
			}
			
			switch (Kind_Adaptation) {
				case GRAD_FLOW: case GRAD_ADJOINT: case GRAD_FLOW_ADJ: case ROBUST: case COMPUTABLE: case COMPUTABLE_ROBUST: case REMAINING:
					cout << "Power of the dual volume in the adaptation sensor: " << DualVol_Power <<endl;
					cout << "Percentage of new elements in the adaptation process: " << New_Elem_Adapt << "."<<endl;
			}
			
			if (Analytical_Surface != NONE)
				cout << "Use analytical definition for including points in the surfaces." <<endl;
			
		}
		
		cout << endl <<"-------------------------- Output information ---------------------------" << endl;
		if (val_software == SU2_CFD) {
			
			cout << "Writing a flow solution every " << Wrt_Sol_Freq <<" iterations."<<endl;
			cout << "Writing the convergence history every " << Wrt_Con_Freq <<" iterations."<<endl;
			
			switch (Output_FileFormat) {
				case PARAVIEW: cout << "The output file format is Paraview (.vtk)." << endl; break;
				case TECPLOT: cout << "The output file format is Tecplot (.plt)." << endl; break;
			}
			
			cout << "Convergence history file name: " << Conv_FileName << "." << endl;
			
			if (!Linerized && !ContAdj) {
				cout << "Surface flow coefficients file name: " << SurfFlowCoeff_FileName << "." << endl;
				cout << "Flow variables file name: " << Flow_FileName << "." << endl;
				cout << "Restart flow file name: " << ReStart_FlowFileName << "." << endl;
			}
			
			if (Linerized) {
				cout << "Linearized flow solution file name: " << Solution_LinFileName << "." << endl;
				cout << "Restart linearized flow file name: " << ReStart_LinFileName << "." << endl;
				cout << "Linearized variables file name: " << Lin_FileName << "." << endl;
				cout << "Surface linearized coefficients file name: " << SurfLinCoeff_FileName << "." << endl;
			}
			
			if (ContAdj || OneShot) {
				cout << "Adjoint solution file name: " << Solution_AdjFileName << "." << endl;
				cout << "Restart adjoint file name: " << ReStart_AdjFileName << "." << endl;
				cout << "Adjoint variables file name: " << Adj_FileName << "." << endl;
				cout << "Surface adjoint coefficients file name: " << SurfAdjCoeff_FileName << "." << endl;
			}
			cout << "Surface(s) to be plotted: ";
			for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++) {
				cout << Marker_Plotting[iMarker_Plotting];
				if (iMarker_Plotting < nMarker_Plotting-1) cout << ", ";
				else cout <<"."<<endl;
			}
		}
		
		if (val_software == SU2_MDC) {
			cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
			if (Visualize_Deformation) cout << "A paraview file will be created to visualize the deformation." << endl;
			else cout << "No paraview file to visualize the deformation." << endl;

		}
		
		if (val_software == SU2_PBC) {
			cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
		}
		
		if (val_software == SU2_GPC) {
			cout << "Output gradient file name: " << ObjFunc_Grad_FileName << ". " << endl;
		}
		if (val_software == SU2_DDC) {
			if (Visualize_Partition) cout << "A paraview file will be created for each domain." << endl;
			else cout << "No paraview file to visualize the domain partitioning." << endl;
		}
		
		if (val_software == SU2_MAC) {
			cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
			cout << "Restart flow file name: " << ReStart_FlowFileName << "." << endl;
			if ((Kind_Adaptation == FULL_ADJOINT) || (Kind_Adaptation == GRAD_ADJOINT) || (Kind_Adaptation == GRAD_FLOW_ADJ) ||
					(Kind_Adaptation == ROBUST) || (Kind_Adaptation == COMPUTABLE_ROBUST) || (Kind_Adaptation == COMPUTABLE) ||
					(Kind_Adaptation == REMAINING)) {
				if (Kind_ObjFunc == DRAG_COEFFICIENT) cout << "Restart adjoint file name: " << ReStart_AdjFileName << "." << endl; 
				if (Kind_ObjFunc == EQUIVALENT_AREA) cout << "Restart adjoint file name: " << ReStart_AdjFileName << "." << endl; 
				if (Kind_ObjFunc == NEARFIELD_PRESSURE) cout << "Restart adjoint file name: " << ReStart_AdjFileName << "." << endl; 
				if (Kind_ObjFunc == LIFT_COEFFICIENT) cout << "Restart adjoint file name: " << ReStart_AdjFileName << "." << endl; 
			}
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
		
		if (nMarker_NS != 0) {
			cout << "Navier-Stokes wall boundary marker(s): ";
			for (iMarker_NS = 0; iMarker_NS < nMarker_NS; iMarker_NS++) {
				cout << Marker_NS[iMarker_NS];
				if (iMarker_NS < nMarker_NS-1) cout << ", ";
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
		
		if (nMarker_Inlet != 0) {
			cout << "Inlet boundary marker(s): ";
			for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
				cout << Marker_Inlet[iMarker_Inlet];
				if (iMarker_Inlet < nMarker_Inlet-1) cout << ", ";
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
		
		switch (BiGrid) {
			case YES: cout << "BiGrid filtering in the gradient computation" << endl; break;
		}
	}
	
}

void CConfig::AddMarkerOption(const string & name, unsigned short & num_marker, string* & marker) {
  //cout << "Adding Marker option " << name << endl;
  num_marker = 0;
  CAnyOptionRef* option_ref = new CMarkerOptionRef(marker, num_marker);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddConvectOption(const string & name, unsigned short & space, unsigned short & centred,
                               unsigned short & upwind) {
  //cout << "Adding Convect option " << name << endl;
  centred = NO_CENTRED;
  upwind = NO_UPWIND;
  space = SPACE_CENTRED;
  CAnyOptionRef* option_ref = new CConvOptionRef(space, centred, upwind);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMathProblem(const string & name, bool & ContAdj, const bool & ContAdj_default,
                             bool & OneShot, const bool & OneShot_default,
                             bool & Linearized, const bool & Linearized_default,
                             bool & Restart_Flow, const bool & Restart_Flow_default) {
  //cout << "Adding Math Problem option " << name << endl;
  ContAdj = ContAdj_default;
  OneShot = OneShot_default;
  Linearized = Linearized_default;
  Restart_Flow = Restart_Flow_default;
  CAnyOptionRef* option_ref = new CMathProblemRef(ContAdj, OneShot, Linearized,
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

void CConfig::AddMarkerInlet(const string & name, unsigned short & nMarker_Inlet,
                             string* & Marker_Inlet, double* & Ttotal, double* & Ptotal,
                             double** & FlowDir) {
  nMarker_Inlet = 0;
  CAnyOptionRef* option_ref = new CMarkerInletRef(nMarker_Inlet, Marker_Inlet,
                                                  Ttotal, Ptotal, FlowDir);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerOutlet(const string & name, unsigned short & nMarker_Outlet, 
                              string* & Marker_Outlet, double* & Pressure) {
  nMarker_Outlet = 0;
  CAnyOptionRef* option_ref = new CMarkerOutletRef(nMarker_Outlet, Marker_Outlet,
                                                   Pressure);
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
      cout << "before_semi = " << before_semi << endl;
      cout << "after_semi = " << after_semi << endl;
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

CConfig::~CConfig(void)
{
	delete [] RK_Alpha_Step;
	delete [] RK_Beta_Step;
	delete [] MG_PreSmooth;
	delete [] MG_PostSmooth;
	delete [] U_FreeStreamND;
}

void CConfig::SetFileNameDomain(unsigned short val_domain)
{
	string old_name;
	char buffer[10]; 
	
	/*--- Standart flow output ---*/
	sprintf (buffer, "_%d", int(val_domain));
	old_name = SurfFlowCoeff_FileName;
	SurfFlowCoeff_FileName = old_name + buffer;
	
	sprintf (buffer, "_%d", int(val_domain));
	old_name = Flow_FileName;
	Flow_FileName = old_name + buffer;
	
	/*--- Standart flow adjoint output ---*/
	sprintf (buffer, "_%d", int(val_domain));
	old_name = Adj_FileName;
	Adj_FileName = old_name + buffer;
	
	sprintf (buffer, "_%d", int(val_domain));
	old_name = SurfAdjCoeff_FileName;
	SurfAdjCoeff_FileName = old_name + buffer;
	
	/*--- Mesh files ---*/	
	sprintf (buffer, "_%d.su2", int(val_domain));
	old_name = Mesh_FileName;
	old_name.erase (old_name.end()-4, old_name.end());
	Mesh_FileName = old_name + buffer;
	
	/*--- Flow solution and restart files ---*/	
	sprintf (buffer, "_%d.dat", int(val_domain));
	old_name = Solution_FlowFileName;
	old_name.erase (old_name.end()-4, old_name.end());
	Solution_FlowFileName = old_name + buffer;
	
	sprintf (buffer, "_%d.dat", int(val_domain));
	old_name = ReStart_FlowFileName;
	old_name.erase (old_name.end()-4, old_name.end());
	ReStart_FlowFileName = old_name + buffer;
	
	/*--- Adjoint flow solution and restart files ---*/
	if (ContAdj || OneShot) {
		sprintf (buffer, "_%d.dat", int(val_domain));
		old_name = Solution_AdjFileName;
		old_name.erase (old_name.end()-4, old_name.end());
		Solution_AdjFileName = old_name + buffer;
		
		sprintf (buffer, "_%d.dat", int(val_domain));
		old_name = ReStart_AdjFileName;
		old_name.erase (old_name.end()-4, old_name.end());
		ReStart_AdjFileName = old_name + buffer;
	}
	
}

unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {
	
	switch (val_eqsystem) {
		case RUNTIME_POT_SYS: return FLOW_SOL;
		case RUNTIME_PLASMA_SYS: return PLASMA_SOL;
		case RUNTIME_FLOW_SYS: return FLOW_SOL;
		case RUNTIME_TURB_SYS: return TURB_SOL;
		case RUNTIME_ELEC_SYS: return ELEC_SOL;
		case RUNTIME_ADJPOT_SYS: return ADJFLOW_SOL;
		case RUNTIME_ADJFLOW_SYS: return ADJFLOW_SOL;
		case RUNTIME_ADJTURB_SYS: return ADJTURB_SOL;
		case RUNTIME_ADJELEC_SYS: return ADJELEC_SOL;
		case RUNTIME_LINPOT_SYS: return LINFLOW_SOL;
		case RUNTIME_LINFLOW_SYS: return LINFLOW_SOL;
		case RUNTIME_LINELEC_SYS: return LINELEC_SOL;
		case RUNTIME_LEVELSET_SYS: return LEVELSET_SOL;
		case RUNTIME_COMBUSTION_SYS: return COMBUSTION_SOL;
		case RUNTIME_MULTIGRID_SYS: return 0;
	}
	
	return 0;
}

void CConfig::SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme, 
																		unsigned short val_kind_centred, unsigned short val_kind_upwind) {
	
	Kind_ConvNumScheme = val_kind_convnumscheme; 
	Kind_Centred = val_kind_centred;
	Kind_Upwind = val_kind_upwind;
}

void CConfig::UpdateCFL(unsigned long val_iter) {
	
	double coeff;
	bool change;
	unsigned short iCFL;
	
	if ((val_iter % int(CFLRamp[1]) == 0 ) && (val_iter != 0)) {
		change = false;
		for (iCFL = 0; iCFL <= nMultiLevel; iCFL++) {
			coeff = pow(MG_CFLRedCoeff, double(iCFL));
			if (ContAdj) coeff = coeff * Adj_CFLRedCoeff;
			
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
}

void CConfig::SetGlobalParam(unsigned short val_solver, unsigned short val_system) {
	switch (val_solver) {
		case POTENTIAL_FLOW:
			if (val_system == RUNTIME_POT_SYS) {
				SetKind_ConvNumScheme(NONE, NONE, NONE);
				
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(NONE);
			}
			break;
		case EULER:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
				SetKind_ViscNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			break;
		case NAVIER_STOKES:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			break;
		case RANS:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE);
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_TURB_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centred_Turb(),
															GetKind_Upwind_Turb());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
			}
			break;
		case MULTI_SPECIES_NAVIER_STOKES:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE);
				SetKind_ViscNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_PLASMA_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Plasma(), GetKind_Centred_Plasma(),
															GetKind_Upwind_Plasma());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Plasma());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Plasma());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Plasma());
			}
			if (val_system == RUNTIME_ELEC_SYS) {
				SetKind_ConvNumScheme(NONE, NONE, NONE);
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
				SetKind_TimeIntScheme(NONE);
			}
			break;
		case TWO_PHASE_FLOW:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_LEVELSET_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_LevelSet(), GetKind_Centred_LevelSet(),
															GetKind_Upwind_LevelSet());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_LevelSet());
			}
			break;
		case COMBUSTION:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_COMBUSTION_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Combustion(), GetKind_Centred_Combustion(),
															GetKind_Upwind_Combustion());
				SetKind_SourNumScheme(GetKind_SourNumScheme_Combustion());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Combustion());
			}
			break;
		case ADJ_EULER:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE); SetKind_ViscNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_ADJFLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centred_AdjFlow(),
															GetKind_Upwind_AdjFlow());
				SetKind_ViscNumScheme(NONE); SetKind_SourNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
			}
			break;
		case ADJ_NAVIER_STOKES:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE);
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_ADJFLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centred_AdjFlow(),
															GetKind_Upwind_AdjFlow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
			}
			break;
		case ADJ_RANS:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_SourNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_ADJFLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centred_AdjFlow(),
															GetKind_Upwind_AdjFlow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
			}
			if (val_system == RUNTIME_ADJTURB_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjTurb(), GetKind_Centred_AdjTurb(),
															GetKind_Upwind_AdjTurb());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjTurb());
				SetKind_SourNumScheme(GetKind_SourNumScheme_AdjTurb());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjTurb());
			}
			break;
		case LIN_EULER:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE); SetKind_ViscNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			if (val_system == RUNTIME_LINFLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_LinFlow(), GetKind_Centred_LinFlow(),
															GetKind_Upwind_LinFlow());
				SetKind_ViscNumScheme(NONE); SetKind_SourNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_LinFlow());
			}
			break;
		case ELECTRIC_POTENTIAL:
			if (val_system == RUNTIME_ELEC_SYS) {
				SetKind_ConvNumScheme(NONE, NONE, NONE);
				SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
				SetKind_TimeIntScheme(NONE);
			}
			break;
		case ADJ_ELECTRIC_POTENTIAL:
			if (val_system == RUNTIME_ADJELEC_SYS) {
				SetKind_ConvNumScheme(NONE, NONE, NONE);
				SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
				SetKind_TimeIntScheme(NONE);
			}
			break;
		case LIN_ELECTRIC_POTENTIAL:
			if (val_system == RUNTIME_LINELEC_SYS) {
				SetKind_ConvNumScheme(NONE, NONE, NONE);
				SetKind_SourNumScheme(GetKind_SourNumScheme_Elec());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Elec());
				SetKind_TimeIntScheme(NONE);
			}
			break;
		case ONE_SHOT_ADJ_EULER:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE);
				SetKind_ViscNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			break;
			if (val_system == RUNTIME_FLOW_ADJFLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centred_AdjFlow(),
															GetKind_Upwind_AdjFlow());
				SetKind_ViscNumScheme(NONE);
				SetKind_SourNumScheme(NONE);
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
			}
			break;
		case ONE_SHOT_ADJ_NAVIER_STOKES:
			if (val_system == RUNTIME_FLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centred_Flow(),
															GetKind_Upwind_Flow());
				SetKind_SourNumScheme(NONE);
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
			}
			break;
			if (val_system == RUNTIME_FLOW_ADJFLOW_SYS) {
				SetKind_ConvNumScheme(GetKind_ConvNumScheme_AdjFlow(), GetKind_Centred_AdjFlow(),
															GetKind_Upwind_AdjFlow());
				SetKind_ViscNumScheme(GetKind_ViscNumScheme_AdjFlow());
				SetKind_SourNumScheme(GetKind_SourNumScheme_AdjFlow());
				SetKind_TimeIntScheme(GetKind_TimeIntScheme_AdjFlow());
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

double CConfig::GetOutlet_Pressure(string val_marker) {
	unsigned short iMarker_Outlet;
	for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
		if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
	return Outlet_Pressure[iMarker_Outlet];
}

void CConfig::SetNondimensionalization(unsigned short val_nDim, int val_rank) { 
	double Mach2Vel_FreeStream, ModVel_FreeStream, Energy_FreeStream = 0.0, ModVel_FreeStreamND;
	unsigned short iDim;
	Velocity_FreeStreamND = new double[val_nDim];
	
	/*--- Local variables and memory allocation ---*/
	double Alpha = AoA*PI_NUMBER/180.0;
	double Beta  = AoS*PI_NUMBER/180.0;
	double Gamma_Minus_One = Gamma - 1.0;
	bool Viscous = (Kind_Solver == RANS || Kind_Solver == NAVIER_STOKES);
//	bool Viscous = (Kind_Solver == RANS || Kind_Solver == NAVIER_STOKES || Kind_Solver == TWO_PHASE_FLOW);
	bool Compressible = (!Incompressible);
	
 	if (val_rank == MASTER_NODE) {
    cout << endl <<"---------------- Flow & Non-dimensionalization information ---------------" << endl;
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
	}
	else {
		if (val_rank == MASTER_NODE)
			cout << "The incompressible solver is non-dimensional." << endl;
		if (Density_FreeStream == 0.0) {
			Density_FreeStream = 1.2886;
			if (val_rank == MASTER_NODE) {
				cout << "Free-stream density is not defined, default value will be used!." << endl;
				cout << "Reference density is not properly defined, default value will be used!." << endl;
			}
			Density_Ref = Density_FreeStream;
		}
		
		Mach2Vel_FreeStream = sqrt(Bulk_Modulus/Density_FreeStream);
		
		/*--- Compute the Free Stream Mach number and the angles ---*/
		ModVel_FreeStream = 0; 
		for (iDim = 0; iDim < val_nDim; iDim++) 
			ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
		ModVel_FreeStream = sqrt(ModVel_FreeStream);
		
		if (ModVel_FreeStream == 0.0) {
			/*--- If the velocity is not provided, the code will use the Mach and AoA ---*/
			cout << "Free-stream velocity is not defined, Mach and AoA will be used!." << endl;
			if (val_nDim == 2) {
				Velocity_FreeStream[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;	
				Velocity_FreeStream[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
			}
			if (val_nDim == 3) {
				Velocity_FreeStream[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
				Velocity_FreeStream[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
				Velocity_FreeStream[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
			}
			ModVel_FreeStream = 0; 
			for (iDim = 0; iDim < val_nDim; iDim++) 
				ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
			ModVel_FreeStream = sqrt(ModVel_FreeStream);
		}
		else {
			/*--- It is necessary to compute the angle of attack, from the velocity ---*/
			Alpha = 0.0; Beta  = 0.0;
			Mach = ModVel_FreeStream / Mach2Vel_FreeStream;
		}
		
	}
	
	if (Compressible) {
		if (Viscous) {
			/*--- For viscous flows, pressure will be computed from a density which
			 is found from the Reynolds number. The viscosity is computed from the Sutherland's law ---*/
			Viscosity_FreeStream = 1.853E-5*(pow(Temperature_FreeStream/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_FreeStream+110.3));
			Density_FreeStream   = Reynolds*Viscosity_FreeStream/(ModVel_FreeStream*Length_Reynolds);
			Pressure_FreeStream  = Density_FreeStream*Gas_Constant*Temperature_FreeStream;
		} else {
			/*--- For inviscid flow, density is calculated from the specified
			 total temperature and pressure using the gas law. ---*/
			Density_FreeStream  = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);    
		}	
		/*-- Compute the freestream energy. ---*/
		Energy_FreeStream = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+0.5*ModVel_FreeStream*ModVel_FreeStream;
	}
	else {
		/*--- The incompressible simulation used the provided desnsity, viscosity, and pressure at the free-stream.
		 The Reynolds number for the incompressible simulation is computed ---*/
		Viscosity_FreeStream = Density_FreeStream*ModVel_FreeStream*Length_Reynolds/Reynolds;
	}
	
	
	/*--- Additional reference values defined by Pref, Tref, RHOref, note 
	 that the incompressible formulation uses the reference velocity that is in the 
	 input file. By definition, Lref is one because we have converted the
   grid to meters.---*/
  Length_Ref        = 1.0;
	if (Compressible) Velocity_Ref = sqrt(Pressure_Ref/Density_Ref);
	else {
		if (Velocity_Ref == 0) {
			cout << "Reference velocity is not defined, the free-stream velocity will be used!." << endl;
			Velocity_Ref = ModVel_FreeStream;
		}
		Pressure_Ref = Density_Ref*(Velocity_Ref*Velocity_Ref);
	}
	
  Time_Ref          = Length_Ref/Velocity_Ref;
  Omega_Ref         = Velocity_Ref/Length_Ref;
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/Temperature_Ref;
	Viscosity_Ref     = Density_Ref*Velocity_Ref*Length_Ref;
	
	
  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  Pressure_FreeStreamND = Pressure_FreeStream/Pressure_Ref;
  Density_FreeStreamND  = Density_FreeStream/Density_Ref;
	
  for (iDim = 0; iDim < val_nDim; iDim++) 
    Velocity_FreeStreamND[iDim] = Velocity_FreeStream[iDim]/Velocity_Ref;
  Temperature_FreeStreamND = Temperature_FreeStream/Temperature_Ref;
	
	Gas_Constant = Gas_Constant/Gas_Constant_Ref;
	if (Rotating_Frame) {
		Omega_FreeStreamND = new double[3];
		for (iDim = 0; iDim < 3; iDim++)
			Omega_FreeStreamND[iDim] = Omega[iDim]/Omega_Ref;
	}
	ModVel_FreeStreamND = 0; 
	for (iDim = 0; iDim < val_nDim; iDim++)
		ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
	ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND);
	Energy_FreeStreamND    = Pressure_FreeStreamND/(Density_FreeStreamND*Gamma_Minus_One)+0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
	Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;
	
  /*--- Write output to the console if this is the master node ---*/
  if (val_rank == MASTER_NODE) {
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
			cout << "Viscous and Inviscid flow: Density, pressure, and " << endl;
			cout << "velocity are based on the freestream values, the " << endl;
			cout << "viscosity is computed from the Reynolds number, and " << endl;
			cout << "Mach number: "<< Mach << ". Computed using the Bulk modulus." << endl;
			if (Viscous) cout << "Reynolds number: " << Reynolds << "."<<endl;
		}
		
		cout <<"--------------------------------------------------------------------------" << endl;
    cout <<"  Input conditions:"<< endl;
		cout << "Grid conversion factor to meters: " << Conversion_Factor << endl;
		
		if (Compressible) {
			cout << "Ratio of specific heats: " << Gamma           << endl;
			cout << "Specific gas constant (J/(kg.K)): "   << Gas_Constant*Gas_Constant_Ref  << endl;
		}
		else
			cout << "Bulk modulus (N/m^2): "						<< Bulk_Modulus    << endl;
		
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
		
		if (Viscous) {
      cout << "Reynolds number (non-dimensional): " << Reynolds << endl;
			cout << "Reynolds length (m): "       << Length_Reynolds     << endl;
    }
		cout << "Reference length (m): "        << Length_Ref      << endl;
		cout << "Reference pressure (N/m^2): "      << Pressure_Ref    << endl;
		if (Compressible)
			cout << "Reference temperature (K): "   << Temperature_Ref << endl;
		cout << "Reference density (kg/m^3): "       << Density_Ref     << endl;
		cout << "Reference velocity (m/s): "       << Velocity_Ref     << endl;
    if (Compressible) {
      cout << "Reference energy (kg.m/s^2): "       << Energy_FreeStream/Energy_FreeStreamND     << endl;
    }
		if (Viscous)
			cout << "Reference viscosity (N.s/m^2): "       << Viscosity_Ref     << endl;
      
    /*--- Print out resulting non-dim values here. ---*/
    cout << "  Resulting non-dimensional state:" << endl;
		if (Compressible) {
			cout << "Specific gas constant (non-dimensional): "   << Gas_Constant << endl;
		}
		cout << "Freestream pressure (non-dimensional): "     << Pressure_FreeStreamND    << endl;
		if (Compressible)
			cout << "Freestream temperature (non-dimensional): "  << Temperature_FreeStreamND << endl;
		cout << "Freestream density (non-dimensional): "      << Density_FreeStreamND     << endl;
		if (val_nDim == 2) {
			cout << "Freestream velocity (non-dimensional): (" << Velocity_FreeStreamND[0] << ",";
			cout << Velocity_FreeStreamND[1] << ")";
		} else if (val_nDim == 3) {
			cout << "Freestream velocity (non-dimensional): (" << Velocity_FreeStreamND[0] << ",";
			cout << Velocity_FreeStreamND[1] << "," << Velocity_FreeStreamND[2] << ")";
		}
		cout << " -> Modulus: "					 << ModVel_FreeStreamND << endl;
		
		if (Compressible)
			cout << "Freestream energy (non-dimensional): "					 << Energy_FreeStreamND << endl;
		
		if (Viscous)
			cout << "Freestream viscosity (non-dimensional): " << Viscosity_FreeStreamND << endl;
		
		if (Rotating_Frame) { 
			cout << "Freestream rotation (non-dimensional): (" << Omega_FreeStreamND[0];
			cout << "," << Omega_FreeStreamND[1] << "," << Omega_FreeStreamND[2] << ")" << endl;
		}
		
		if (Mach <= 0.0) cout << "Force coefficients computed using reference values." << endl;
		else cout << "Force coefficients computed using freestream values." << endl;
	}
}
