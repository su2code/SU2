/*!
 * \file solution_direct_plasma.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

/*!
 * \class  CPlasmaSolver
 * \brief Initialization of the plasma solution class
 * \author A. Lonkar
 */
CPlasmaSolver::CPlasmaSolver(void) : CSolver() { }

CPlasmaSolver::CPlasmaSolver(CGeometry *geometry, CConfig *config) : CSolver() {
	unsigned long iPoint, index;
	unsigned short iVar, iDim, iSpecies, iMarker, nPrimVar;
	double Vel2 = 0.0;

	restart = (config->GetRestart() || config->GetRestart_Flow());
	implicit = config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT;
	centered_scheme = config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED;
	axisymmetric = config->GetAxisymmetric();
	weighted_LS = config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES;
	roe_turkel = ((config->GetKind_Upwind_Plasma() == ROE_TURKEL_1ST) || (config->GetKind_Upwind_Plasma() == ROE_TURKEL_2ND) );
	magnet = (config->GetMagnetic_Force() == YES);

#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#endif

	double Density_Inf_mean, NumDensity_Inf_mean, Pressure_Inf_mean, Temperature_Inf_mean, Mach_Inf_mean;
	double Gas_Constant_mean, SoundSpeed_mean, Gamma, Molar_Mass_mean, CharVibTemp;
	double Gas_Constant;
	double NumberDensity, NumberDensity_ratio;

	double Pressure_Inlet, Pressure_Outlet, Temperature_Inlet, Temperature_Outlet;

	double AoA = (config->GetAoA()*PI_NUMBER) / 180.0;
	double AoS = (config->GetAoS()*PI_NUMBER) / 180.0;

	/*--- Get relevant constants ---*/
	nDim        = geometry->GetnDim();
	nMonatomics = config->GetnMonatomics();
	nDiatomics  = config->GetnDiatomics();
	nSpecies    = config->GetnSpecies();
	nVar        = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);
	node        = new CVariable*[nPoint];
	nPrimVar    = nDim+3;
	nPoint      = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();

	/*--- Define auxiliary vectors for residuals at nodes i & j ---*/
	Residual = new double[nVar];    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]      = 0.0;
	Residual_Max = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]      = 0.0;
	Point_Max = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]      = 0;
	Residual_i = new double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]      = 0.0;
	Residual_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]      = 0.0;
	Residual_Chemistry = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_Chemistry[iVar]      = 0.0;
	Residual_MomentumExch = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_MomentumExch[iVar]      = 0.0;
	Residual_ElecForce = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_ElecForce[iVar]      = 0.0;
	Residual_EnergyExch = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_EnergyExch[iVar]      = 0.0;
	Res_Conv = new double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
	Res_Visc = new double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
	Res_Sour = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;

	if (axisymmetric) {
		Residual_Axisymmetric = new double[nVar];
	}

	Species_Delta = new double [nSpecies];

	Min_Delta_Time = new double[nSpecies];
	Max_Delta_Time = new double[nSpecies];

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_RMS(iVar, 0.0);

	/*--- Define auxiliary vectors for the solution at nodes i & j ---*/
	Solution = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]      = 0.0;
	Solution_i = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]      = 0.0;
	Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]      = 0.0;

	/*--- Define auxiliary vectors for the geometry ---*/
	Vector = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim]   = 0.0;
	Vector_j = new double[nDim];for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim]   = 0.0;
	PrimVar_i = new double [nDim+6]; for (iVar = 0; iVar < nDim+6; iVar++) PrimVar_i[iVar]      = 0.0;
	PrimVar_j = new double [nDim+6]; for (iVar = 0; iVar < nDim+6; iVar++) PrimVar_j[iVar]      = 0.0;

	PrimVar_max = new double**[nPoint];
	PrimVar_min = new double**[nPoint];
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		PrimVar_max[iPoint] = new double*[nSpecies];
		PrimVar_min[iPoint] = new double*[nSpecies];
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			PrimVar_max[iPoint][iSpecies] = new double[nPrimVar];
			PrimVar_min[iPoint][iSpecies] = new double[nPrimVar];
		}
	}

	PrimGrad_i = new double**[nSpecies];
	PrimGrad_j = new double**[nSpecies];
	PrimGradLimiter_i = new double*[nSpecies];
	PrimGradLimiter_j = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		PrimGrad_i[iSpecies] = new double*[nPrimVar];
		PrimGrad_j[iSpecies] = new double*[nPrimVar];
		PrimGradLimiter_i[iSpecies] = new double[nPrimVar];
		PrimGradLimiter_j[iSpecies] = new double[nPrimVar];
		for (iVar = 0; iVar < nPrimVar; iVar++){
			PrimGrad_i[iSpecies][iVar] = new double[nDim];
			PrimGrad_j[iSpecies][iVar] = new double[nDim];
		}
	}


	/*--- Define some auxiliary vector related with the undivided lapalacian computation ---*/
	if (centered_scheme) {
		p1_Und_Lapl = new double * [nSpecies];
		p2_Und_Lapl = new double * [nSpecies];
		for (iSpecies =0; iSpecies < nSpecies; iSpecies ++ ) {
			p1_Und_Lapl[iSpecies] = new double [nPoint];
			p2_Und_Lapl[iSpecies] = new double [nPoint];
		}
	}

	if (roe_turkel) {
		Precon_Mat_inv = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar ++)
			Precon_Mat_inv[iVar] = new double[nVar];
	}

	if (magnet) Mag_Force = new double [nDim];

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (implicit) {
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar]; Jacobian_j = new double* [nVar];
		Jacobian_Chemistry    = new double* [nVar];
		Jacobian_ElecForce    = new double* [nVar];
		Jacobian_MomentumExch = new double* [nVar];
		Jacobian_EnergyExch   = new double* [nVar];

		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; Jacobian_j[iVar] = new double [nVar];
			Jacobian_Chemistry[iVar]    = new double [nVar];
			Jacobian_ElecForce[iVar]    = new double [nVar];
			Jacobian_MomentumExch[iVar] = new double [nVar];
			Jacobian_EnergyExch[iVar]   = new double [nVar];
		}

		if (axisymmetric) {
			Jacobian_Axisymmetric = new double*[nVar];
			for (iVar = 0; iVar < nVar; iVar++)
				Jacobian_Axisymmetric[iVar] = new double[nVar];
		}

		/*--- Initialization of the structure of the whole Jacobian ---*/
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    
	}

	/*--- Computation of gradients by least squares ---*/
	if (weighted_LS) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++)
			cvector[iVar] = new double [nDim];
	}

	/*--- Forces definition and coefficient in all the markers ---*/
	nMarker = config->GetnMarker_All();
	CPressure = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CPressure[iMarker] = new double [geometry->nVertex[iMarker]];

	CSkinFriction = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];

	//	CHeatTransfer = new double* [nMarker];
	//	for (iMarker = 0; iMarker < nMarker; iMarker++)
	//		CHeatTransfer[iMarker] = new double [geometry->nVertex[iMarker]];


	CHeatTransfer   = new double** [nMarker];
	CViscForce      = new double*** [nMarker];
	CPressForce     = new double*** [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		CHeatTransfer[iMarker]  = new double* [nSpecies];
		CViscForce[iMarker]     = new double** [nSpecies];
		CPressForce[iMarker]    = new double** [nSpecies];
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			CViscForce[iMarker][iSpecies]   = new double*[nDim];
			CPressForce[iMarker][iSpecies]  = new double*[nDim];
			CHeatTransfer[iMarker][iSpecies]        = new double[geometry->nVertex[iMarker]];
			for (iDim = 0; iDim < nDim; iDim++) {
				CViscForce[iMarker][iSpecies][iDim]     = new double[geometry->nVertex[iMarker]];
				CPressForce[iMarker][iSpecies][iDim]    = new double[geometry->nVertex[iMarker]];
			}
		}
	}

	ForceInviscid  = new double[nDim];
	MomentInviscid = new double[3];
	CDrag_Inv      = new double[nMarker];
	CLift_Inv      = new double[nMarker];
	CSideForce_Inv = new double[nMarker];
	CMx_Inv        = new double[nMarker];
	CMy_Inv        = new double[nMarker];
	CMz_Inv        = new double[nMarker];
	CEff_Inv       = new double[nMarker];
	CEquivArea_Inv = new double[nMarker];
	CNearFieldOF_Inv = new double[nMarker];
	CFx_Inv        = new double[nMarker];
	CFy_Inv        = new double[nMarker];
	CFz_Inv        = new double[nMarker];
	CMerit_Inv    = new double[nMarker];
	CT_Inv      = new double[nMarker];
	CQ_Inv      = new double[nMarker];
	Total_CDrag = 0.0; Total_CLift = 0.0; Total_CSideForce = 0.0;
	Total_CMx = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
	Total_CEff = 0.0; Total_CEquivArea = 0.0; Total_CNearFieldOF = 0.0;
	Total_CFx = 0.0; Total_CFy = 0.0; Total_CFz = 0.0;
	Total_CT = 0.0; Total_CQ = 0.0; Total_CMerit = 0.0;

	ForceViscous = new double[nDim];
	MomentViscous = new double[nDim];
	CDrag_Visc = new double[config->GetnMarker_All()];
	CLift_Visc = new double[config->GetnMarker_All()];
	CMx_Visc = new double[config->GetnMarker_All()];
	CMy_Visc = new double[config->GetnMarker_All()];
	CMz_Visc = new double[config->GetnMarker_All()];
	CEff_Visc = new double[config->GetnMarker_All()];
	CFx_Visc = new double[config->GetnMarker_All()];
	CFy_Visc = new double[config->GetnMarker_All()];
	CFz_Visc = new double[config->GetnMarker_All()];
	CMerit_Visc = new double[config->GetnMarker_All()];
	CT_Visc = new double[config->GetnMarker_All()];
	CQ_Visc = new double[config->GetnMarker_All()];
  Q_Visc = new double[config->GetnMarker_All()];


	Prandtl_Lam   = config->GetPrandtl_Lam();



	/*--- Flow infinity array initialization ---*/
	Velocity_Inf_mean = new double[nDim];

	Density_Inf     = new double [nSpecies];
	Pressure_Inf    = new double [nSpecies];
	Mach_Inf        = new double [nSpecies];
	Temperature_Inf = new double [nSpecies];
	Energy_Inf      = new double [nSpecies];
	Energy_vib_Inf  = new double [nSpecies];

	Density_Inlet  = new double [nSpecies];
	Mach_Inlet 		 = new double [nSpecies];
	Energy_Inlet   = new double [nSpecies];

	Density_Outlet  = new double [nSpecies];
	Energy_Outlet   = new double [nSpecies];
	Mach_Outlet     = new double [nSpecies];

	Velocity_Inf   = new double*[nSpecies];
	Velocity_Inlet = new double*[nSpecies];
	Velocity_Outlet = new double*[nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Velocity_Inf[iSpecies] = new double [nDim];
		Velocity_Inlet[iSpecies] = new double [nDim];
		Velocity_Outlet[iSpecies] = new double [nDim];
	}

	/*--- Get Mean Flow Properties from the configuration file ---*/
	Pressure_Inf_mean    = config->GetPressure_FreeStream();
	Temperature_Inf_mean = config->GetTemperature_FreeStream();
	Mach_Inf_mean        = config->GetMach_FreeStreamND();
	//	Gas_Constant_mean    = config->GetGas_ConstantND();
	//	SoundSpeed_mean      = sqrt(config->GetGamma()*Gas_Constant_mean*Temperature_Inf_mean);
	Molar_Mass_mean      = config->GetMixtureMolar_Mass();

	/*--- Initialize flow conditions ---*/
	Density_Inf_mean = 0.0;
	Gas_Constant_mean = 0.0;
	SoundSpeed_mean = 0.0;
	NumDensity_Inf_mean = Pressure_Inf_mean / (BOLTZMANN_CONSTANT*Temperature_Inf_mean);
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		NumberDensity             = config->GetInitial_Gas_Composition(iSpecies) * NumDensity_Inf_mean;
		Molar_Mass                = config->GetMolar_Mass(iSpecies);
		Gas_Constant              = config->GetSpecies_Gas_Constant(iSpecies);
		Gamma                     = config->GetSpecies_Gamma(iSpecies);
		Temperature_Inf[iSpecies] = config->GetSpecies_Temperature(iSpecies);
		Density_Inf[iSpecies]     = NumberDensity*Molar_Mass/AVOGAD_CONSTANT;
		Pressure_Inf[iSpecies]    = Density_Inf[iSpecies]*Gas_Constant*Temperature_Inf[iSpecies];
		Density_Inf_mean         += Density_Inf[iSpecies];
		Gas_Constant_mean        += Density_Inf[iSpecies]*Gas_Constant;
		SoundSpeed_mean          += Density_Inf[iSpecies]*sqrt(Gamma*Gas_Constant*Temperature_Inf[iSpecies]);
	}
	Gas_Constant_mean = Gas_Constant_mean/Density_Inf_mean;
	SoundSpeed_mean   = SoundSpeed_mean/Density_Inf_mean;

	/*--- Calculate mean flow velocity ---*/
	if (nDim == 2) {
		Velocity_Inf_mean[0] = cos(AoA) * Mach_Inf_mean * SoundSpeed_mean;
		Velocity_Inf_mean[1] = sin(AoA) * Mach_Inf_mean * SoundSpeed_mean;
	}
	if (nDim == 3) {
		Velocity_Inf_mean[0] = cos(AoA) * cos(AoS) * Mach_Inf_mean * SoundSpeed_mean;
		Velocity_Inf_mean[1] = sin(AoS) * Mach_Inf_mean * SoundSpeed_mean;
		Velocity_Inf_mean[2] = sin(AoA) * cos(AoS) * Mach_Inf_mean * SoundSpeed_mean;
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Gas_Constant       = config->GetSpecies_Gas_Constant(iSpecies);
		Gamma              = config->GetSpecies_Gamma(iSpecies);
		Molar_Mass         = config->GetMolar_Mass(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Vel2         = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_Inf[iSpecies][iDim] = Velocity_Inf_mean[iDim];
			Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];
		}
		Mach_Inf[iSpecies] = sqrt(Vel2/(Gamma*Gas_Constant*Temperature_Inf[iSpecies]));
		Energy_vib_Inf[iSpecies] = 0.0;
		Energy_el_Inf = 0.0;
		if (iSpecies < nDiatomics) {
			CharVibTemp = config->GetCharVibTemp(iSpecies);
			Energy_vib_Inf[iSpecies] = UNIVERSAL_GAS_CONSTANT/Molar_Mass * CharVibTemp / (exp(CharVibTemp/Temperature_Inf[iSpecies]) - 1.0);
		}
		Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*(Gamma-1.0)) + 1.0/2.0*Vel2 + Enthalpy_formation + Energy_vib_Inf[iSpecies] + Energy_el_Inf;
	}

#ifdef NO_MPI
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		cout << " iSpecies = " << iSpecies << endl;
		cout << " Temperature_inf = " << Temperature_Inf[iSpecies] << endl;
		cout << " Pressure_Inf = " << Pressure_Inf[iSpecies] << endl;
		cout << " Density_Inf = " << Density_Inf[iSpecies] << endl;
		cout << " Mach_Inf = " << Mach_Inf[iSpecies] << endl;
		cout << " Velocity_Inf = " << Velocity_Inf[iSpecies][0] << endl;
		cout << " *********" << endl;
	}
#else
	if (rank == MASTER_NODE) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			cout << " iSpecies = " << iSpecies << endl;
			cout << " Temperature_inf = " << Temperature_Inf[iSpecies] << endl;
			cout << " Pressure_Inf = " << Pressure_Inf[iSpecies] << endl;
			cout << " Density_Inf = " << Density_Inf[iSpecies] << endl;
			cout << " Mach_Inf = " << Mach_Inf[iSpecies] << endl;
			cout << " Velocity_Inf = " << Velocity_Inf[iSpecies][0] << endl;
			cout << " *********" << endl;
		}
	}
#endif

	bool Inlet_Outlet_Defined = config->GetInletConditionsDefined();

	if(!Inlet_Outlet_Defined) {
#ifdef NO_MPI
		cout << " Now defining the inlet and outlet conditions " << endl;
#else
		if (rank == MASTER_NODE)
			cout << " Now defining the inlet and outlet conditions " << endl;
#endif

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

			Density_Inlet[iSpecies]  = Density_Inf[iSpecies];
			Density_Outlet[iSpecies] = Density_Inf[iSpecies];

			Energy_Inlet[iSpecies]   = Energy_Inf[iSpecies];
			Energy_Outlet[iSpecies]  = Energy_Inf[iSpecies];

			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity_Inlet[iSpecies][iDim] = Velocity_Inf_mean[iDim];
				Velocity_Outlet[iSpecies][iDim] = Velocity_Inf_mean[iDim];
			}
		}
	}
	else {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

			Pressure_Inlet 				= config->GetInlet_Species_Pressure(iSpecies);
			Pressure_Outlet  			= config->GetOutlet_Species_Pressure(iSpecies);

			Temperature_Inlet  			= config->GetInlet_Species_Temperature(iSpecies);
			Temperature_Outlet  		= config->GetOutlet_Species_Temperature(iSpecies);

			Gas_Constant                = config->GetSpecies_Gas_Constant(iSpecies);
			Gamma                       = config->GetSpecies_Gamma(iSpecies);

			Density_Inlet[iSpecies] 	= Pressure_Inlet/(Gas_Constant * Temperature_Inlet);
			Density_Outlet[iSpecies] 	= Pressure_Outlet/(Gas_Constant * Temperature_Outlet);

			Velocity_Inlet[iSpecies][0]  = config->GetInlet_Species_Velocity(iSpecies);
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Inlet[iSpecies][iDim] * Velocity_Inlet[iSpecies][iDim];

			Energy_Inlet[iSpecies]      = Pressure_Inlet/(Density_Inlet[iSpecies]*(Gamma-1.0)) + 0.5*Vel2;

			Velocity_Outlet[iSpecies][0] = config->GetOutlet_Species_Velocity(iSpecies);
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += Velocity_Outlet[iSpecies][iDim] * Velocity_Outlet[iSpecies][iDim];

			Energy_Outlet[iSpecies]      = Pressure_Outlet/(Density_Outlet[iSpecies]*(Gamma-1.0)) + 0.5*Vel2;


			for (iDim = 1; iDim < nDim; iDim ++) {
				Velocity_Inlet[iSpecies][iDim] = 0.0;
				Velocity_Outlet[iSpecies][iDim] = 0.0;
			}
		}
	}

	/*--- Initialize the solution from configuration file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CPlasmaVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, nDim, nVar, config);
		}

		/*--- Restart the solution from the specified restart file ---*/
	} else {
		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
		restart_file.open(filename.data(), ios::in);

		/*--- In case there is no restart file ---*/
		if (restart_file.fail()) {
			cout << "There is no plasma restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}

		/*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
		long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];

		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
			Global2Local[iPoint] = -1;

		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;

		/*--- The first line is the header ---*/
		getline (restart_file, text_line);

		while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);

			/*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
			iPoint_Local = Global2Local[iPoint_Global];
			if (iPoint_Local >= 0) {

				/*--- Restart from a previous plasma simulation ---*/
				if (!restart_from_Euler) {

					/*--- First value is the point index, then the conservative vars. ---*/
					point_line >> index;
					for (iVar = 0; iVar < nVar; iVar++)
						point_line >> Solution[iVar];

					node[iPoint_Local] = new CPlasmaVariable(Solution, nDim, nVar, config);

					/*--- Restart from a previous euler simulation ---*/
				} else {

					/*--- Load the line from the Euler restart file. ---*/
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];


					double Molecular_Mass_mean;
					Molecular_Mass_mean = 0.0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						Molar_Mass           = config->GetMolar_Mass(iSpecies);
						NumberDensity_ratio  = config->GetInitial_Gas_Composition(iSpecies);
						Molecular_Mass_mean += (Molar_Mass/AVOGAD_CONSTANT)*NumberDensity_ratio;
					}
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

						/*--- Initialize all species velocites equal to the mean flow velocities ---*/
						NumberDensity_ratio       = config->GetInitial_Gas_Composition(iSpecies);
						Molar_Mass                = config->GetMolar_Mass(iSpecies);
						Enthalpy_formation        = config->GetEnthalpy_Formation(iSpecies);
						Density_Inf[iSpecies]  	  = (Molar_Mass/AVOGAD_CONSTANT) * NumberDensity_ratio
								* Solution[0] / Molecular_Mass_mean;
						for (iDim = 0; iDim < nDim; iDim++) {
							Velocity_Inf[iSpecies][iDim] = Solution[iDim+1]/Solution[0];
						}
						Energy_vib_Inf[iSpecies] = 0.0;
						Energy_Inf[iSpecies]     = Solution[nDim+1]/Solution[0] + Enthalpy_formation;

					}
					node[iPoint_Local] = new CPlasmaVariable(Solution, nDim, nVar, config);
				}
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CPlasmaVariable(Solution, nDim, nVar, config);
		}

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
}

CPlasmaSolver::~CPlasmaSolver(void) {
	unsigned short iVar, iDim, iMarker;
	unsigned long iPoint;

	/*--- Delete variables stored in at each grid node ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		delete node[iPoint];
	delete [] node;

	delete [] Residual;		delete [] Residual_Max;
	delete [] Residual_i;	delete [] Residual_j;
	delete [] Residual_Chemistry;			delete [] Residual_ElecForce;
	delete [] Residual_MomentumExch;	delete[] Residual_EnergyExch;
	delete [] Res_Conv;		delete [] Res_Visc;		delete [] Res_Sour;
	if (axisymmetric)
		delete [] Residual_Axisymmetric;

	delete [] Solution;		delete [] Solution_i; delete [] Solution_j;
	delete [] PrimVar_i;  delete [] PrimVar_j;
	delete [] Vector;			delete [] Vector_i;		delete [] Vector_j;

	delete [] Species_Delta;

	for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		for (iVar = 0; iVar < nDim+3; iVar++){
			delete [] PrimGrad_i[iSpecies][iVar];
			delete [] PrimGrad_j[iSpecies][iVar];
		}
		delete [] PrimGradLimiter_i[iSpecies];
		delete [] PrimGradLimiter_j[iSpecies];
	}
	delete [] PrimGrad_i;
	delete [] PrimGrad_j;
	delete [] PrimGradLimiter_i;
	delete [] PrimGradLimiter_j;

	for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
		for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			delete [] PrimVar_max[iPoint][iSpecies];
			delete [] PrimVar_min[iPoint][iSpecies];
		}
		delete [] PrimVar_max[iPoint];
		delete [] PrimVar_min[iPoint];
	}
	delete [] PrimVar_max;
	delete [] PrimVar_min;

	delete [] Min_Delta_Time;
	delete [] Max_Delta_Time;

	if (centered_scheme) {
		for (unsigned short iSpecies =0; iSpecies < nSpecies; iSpecies ++) {
			delete [] p1_Und_Lapl[iSpecies];
			delete [] p2_Und_Lapl[iSpecies];
		}
		delete [] p1_Und_Lapl;
		delete [] p1_Und_Lapl;
	}

	if (roe_turkel) {
		for (iVar = 0; iVar < nVar; iVar ++)
			delete [] Precon_Mat_inv[iVar];
		delete [] Precon_Mat_inv;
	}

	if (magnet) delete Mag_Force;

	/*--- De-allocate Jacobians if running an implicit calculation ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			delete [] Jacobian_i[iVar];
			delete [] Jacobian_j[iVar];
			delete [] Jacobian_Chemistry[iVar];
			delete [] Jacobian_ElecForce[iVar];
			delete [] Jacobian_MomentumExch[iVar];
			delete [] Jacobian_EnergyExch[iVar];
			if (axisymmetric)
				delete [] Jacobian_Axisymmetric[iVar];
		}
		delete [] Jacobian_i;							delete [] Jacobian_j;
		delete [] Jacobian_Chemistry;			delete [] Jacobian_ElecForce;
		delete [] Jacobian_MomentumExch;	delete [] Jacobian_EnergyExch;
		if (axisymmetric)
			delete [] Jacobian_Axisymmetric;

	}

	/*--- De-allocate arrays associated with the weighted least-squares gradient calculation method ---*/
	if (weighted_LS) {
		for (iDim = 0; iDim < this->nDim; iDim++)
			delete [] Smatrix[iDim];
		delete [] Smatrix;

		for (iVar = 0; iVar < nDim+3; iVar++)
			delete [] cvector[iVar];
		delete [] cvector;
	}

	/*--- De-allocate forces and coefficients on all markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		delete [] CPressure[iMarker];
		delete [] CSkinFriction[iMarker];
		for (unsigned short iSpecies = 0; iVar < nSpecies; iVar ++) {
			delete [] CHeatTransfer[iMarker][iSpecies];
			for (unsigned short iDim = 0; iDim < nDim; iDim ++) {
				delete [] CViscForce[iMarker][iSpecies][iDim];
				delete [] CPressForce[iMarker][iSpecies][iDim];
			}
		}
	}


	delete [] CPressure; delete [] CSkinFriction;	delete [] CHeatTransfer;
	delete [] CViscForce; delete [] CPressForce;

	delete [] ForceInviscid;
	delete [] MomentInviscid;
	delete [] CDrag_Inv;
	delete [] CLift_Inv;
	delete [] CSideForce_Inv;
	delete [] CMx_Inv;
	delete [] CMy_Inv;
	delete [] CMz_Inv;
	delete [] CEff_Inv;
	delete [] CEquivArea_Inv;
	delete [] CNearFieldOF_Inv;
	delete [] CFx_Inv;
	delete [] CFy_Inv;
	delete [] CFz_Inv;
	delete [] CMerit_Inv;
	delete [] CT_Inv;
	delete [] CQ_Inv;

	delete [] ForceViscous;			delete [] MomentViscous;		delete [] CDrag_Visc;
	delete [] CMx_Visc;				delete [] CMy_Visc;				delete [] CMz_Visc;
	delete [] CFx_Visc;				delete [] CFy_Visc;				delete [] CFz_Visc;
	delete [] CEff_Visc;			delete [] CMerit_Visc;			delete [] CT_Visc;
	delete [] CQ_Visc;
  delete [] Q_Visc;

	/*--- De-allocate flow infinity condition arrays ---*/
	delete [] Velocity_Inf_mean;

	delete [] Density_Inf;
	delete [] Pressure_Inf;
	delete [] Mach_Inf;
	delete [] Temperature_Inf;
	delete [] Energy_Inf;
	delete [] Energy_vib_Inf;

	delete [] Density_Inlet;	delete [] Density_Outlet;
	delete [] Energy_Inlet;		delete [] Energy_Outlet;
	delete [] Mach_Inlet;			delete [] Mach_Outlet;

	for (iVar = 0; iVar < nSpecies; iVar ++) {
		delete [] Velocity_Inf[iVar];
		delete [] Velocity_Inlet[iVar];
		delete [] Velocity_Outlet[iVar];
	}
	delete[] Velocity_Inf; 	delete [] Velocity_Inlet; delete [] Velocity_Outlet;

}

void CPlasmaSolver::SetVel_Residual_Zero(unsigned long val_ipoint, unsigned short iSpecies) {
	unsigned short loc, iDim;
	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
  
	for (iDim = 0; iDim < nDim; iDim++)
		LinSysRes[val_ipoint*nVar+loc+iDim+1] = 0.0;
}

void CPlasmaSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {

	unsigned short iVar, iMarker, iPeriodic_Index, iSpecies, loc;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi,
	sinPsi, *newSolution = NULL, *Buffer_Receive_U = NULL;
	short SendRecv;
	int send_to, receive_from;

#ifndef NO_MPI
	double *Buffer_Send_U = NULL, *Buffer_Send_LaminarViscosity = NULL,
			*Buffer_Send_EddyViscosity = NULL, *Buffer_Send_VGrad = NULL, *Buffer_Send_UGrad = NULL, *Buffer_Send_Limit = NULL, *Buffer_Send_Undivided_Laplacian = NULL,
			*Buffer_Send_Sensor = NULL, *Buffer_Send_Lambda = NULL;
	unsigned short *Buffer_Send_Neighbor = NULL;
#endif

	newSolution = new double[nVar];

	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];

			nBuffer_Vector			= nVertex*nVar;

			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;

#ifndef NO_MPI

			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
				/*--- Allocate upwind variables ---*/
				Buffer_Send_U = new double[nBuffer_Vector];

				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Copy data ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						Buffer_Send_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
					}
				}

				/*--- Send the buffer, and deallocate information using MPI ---*/
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_Vector, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_U;
			}

#endif

			/*--- Receive information  ---*/
			if (SendRecv < 0) {

				/*--- Allocate upwind variables ---*/
				Buffer_Receive_U = new double [nBuffer_Vector];


#ifdef NO_MPI

				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Copy data ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						Buffer_Receive_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
					}
				}

#else
				/*--- Receive the information using MPI---*/
				MPI::COMM_WORLD.Recv(Buffer_Receive_U, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);

#endif

				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();

					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);

					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);

					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;

					/*--- Copy conserved variables before performing transformation. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						newSolution[iVar] = Buffer_Receive_U[iVar*nVertex+iVertex];


					/*--- Rotate the momentum components. ---*/
					if (nDim == 2) {
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
							if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
							newSolution[loc + 1] = rotMatrix[0][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex];
							newSolution[loc + 2] = rotMatrix[1][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex];
						}
					}
					else {
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
							if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
							newSolution[loc + 1] = rotMatrix[0][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_U[(loc + 3)*nVertex+iVertex];
							newSolution[loc + 2] = rotMatrix[1][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_U[(loc + 3)*nVertex+iVertex];
							newSolution[loc + 3] = rotMatrix[2][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_U[(loc + 3)*nVertex+iVertex];
						}
					}
					/*--- Copy transformed conserved variables back into buffer. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar*nVertex+iVertex] = newSolution[iVar];


					/*--- Upwind method. Store the received information ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
					}
				}
				delete [] Buffer_Receive_U;
			}
		}
	}
	delete [] newSolution;
}

void CPlasmaSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index;
	unsigned long iVertex, iPoint, nVertex, nBuffer_VectorGrad;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi,
	sinPsi, **newGradient = NULL, *Buffer_Receive_UGrad = NULL;
	short SendRecv;
	int send_to, receive_from;
    
#ifndef NO_MPI
    
    MPI::COMM_WORLD.Barrier();
	double *Buffer_Send_UGrad = NULL;
    
#endif
    
	newGradient = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		newGradient[iVar] = new double[3];
    
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_VectorGrad = nVertex*nVar*nDim;
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
            
#ifndef NO_MPI
            
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
                Buffer_Send_UGrad = new double[nBuffer_VectorGrad];
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex] = node[iPoint]->GetGradient(iVar,iDim);
				}
                MPI::COMM_WORLD.Bsend(Buffer_Send_UGrad, nBuffer_VectorGrad, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_UGrad;
			}
            
#endif
            
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
                Buffer_Receive_UGrad = new double [nBuffer_VectorGrad];
                
#ifdef NO_MPI
                
				/*--- Receive information without MPI ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex] = node[iPoint]->GetGradient(iVar,iDim);
				}
                
#else
                
                MPI::COMM_WORLD.Recv(Buffer_Receive_UGrad, nBuffer_VectorGrad, MPI::DOUBLE, receive_from, 0);
                
#endif
                
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
                    
					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();
                    
					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);
                    
					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
                    
					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;
                    
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            newGradient[iVar][iDim] = Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex];
                    
                    /*--- Need to rotate the gradients for all conserved variables. ---*/
                    for (iVar = 0; iVar < nVar; iVar++) {
                        if (nDim == 2) {
                            newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex];
                            newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex];
                        }
                        else {
                            newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_UGrad[2*nVar*nVertex+iVar*nVertex+iVertex];
                            newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_UGrad[2*nVar*nVertex+iVar*nVertex+iVertex];
                            newGradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_UGrad[2*nVar*nVertex+iVar*nVertex+iVertex];
                        }
                    }
                    
                    /*--- Copy transformed gradients back into buffer. ---*/
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex] = newGradient[iVar][iDim];
                    
                    
					/*--- Store the received information ---*/
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            node[iPoint]->SetGradient(iVar, iDim, Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex]);
				}
                delete [] Buffer_Receive_UGrad;
			}
		}
	}
    
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] newGradient[iVar];
	delete [] newGradient;
    
#ifndef NO_MPI
    
    MPI::COMM_WORLD.Barrier();
    
#endif
    
}


void CPlasmaSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	bool upwind_2nd = ((config->GetKind_Upwind() == ROE_2ND) || (config->GetKind_Upwind() == AUSM_2ND)
			|| (config->GetKind_Upwind() == HLLC_2ND) || (config->GetKind_Upwind() == ROE_TURKEL_2ND) || (config->GetKind_Upwind() == SW_2ND) || (config->GetKind_Upwind() == MSW_2ND));
	bool viscous = (config->GetKind_Solver() == PLASMA_NAVIER_STOKES);
	bool limiter = ((config->GetKind_SlopeLimit() != NONE));

	for (iPoint = 0; iPoint < nPoint; iPoint ++) {

		/*--- Compute primitive variables ---*/
		node[iPoint]->SetPrimVar(config, geometry->node[iPoint]->GetCoord());

		/*--- Initialize the convective and viscous residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);

	}

	/*--- Upwind second order flux reconstruction ---*/
	if ((upwind_2nd) && ((iMesh == MESH_0))) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

		/*--- Limiter computation ---*/
		if (limiter) SetSolution_Limiter(geometry, config);
	}

	if (viscous) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimVar_Gradient_GG(geometry, config);
		else if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimVar_Gradient_LS(geometry, config);

		/*--- Set the limiter for the primitive variables ---*/
		SetPrimVar_Limiter(geometry, config);
	}

	/*--- Initialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CPlasmaSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Area, dV, dij;
	double Mean_SoundSpeed, Mean_ProjVel, Mean_LaminarVisc, Mean_ThermalConductivity, Mean_ThermalConductivity_vib, Mean_Density;
	double CharVibTemp, Mean_Tvib, dTvdEv;
	double Lambda, Lambda_1, Lambda_2, Lambda_iSpecies, Lambda_Visc_iSpecies;
	double Local_Delta_Time, Local_Delta_Time_Visc, Global_Delta_Time = 1E6;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker, iSpecies, loc;
	double K_v, CFL;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Min_Delta_Time[iSpecies] = 1.E6;
		Max_Delta_Time[iSpecies] = 0.0;
	}

	K_v = 0.25;
	bool viscous = false;
	bool MultipleTimeSteps = (config->MultipleTimeSteps());
	bool centered = (config->GetKind_ConvNumScheme() == SPACE_CENTERED);
	if ((config->GetKind_Solver() == PLASMA_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES)) viscous = true;

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			node[iPoint]->SetMax_Lambda_Inv(0.0, iSpecies);
			node[iPoint]->SetMax_Lambda_Visc(0.0, iSpecies);
			if (centered) node[iPoint]->SetLambda(0.0,iSpecies);
		}
	}

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Identify points in the edge, the normal and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();

		Area = 0.0; dij = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Area += Normal[iDim]*Normal[iDim];
			dij += pow(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim), 2.0);
		}
		Area = sqrt(Area);
		dij  = sqrt(dij);

		/*--- Mean Values ---*/
		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal,iSpecies) + node[jPoint]->GetProjVel(Normal,iSpecies));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed(iSpecies) + node[jPoint]->GetSoundSpeed(iSpecies)) * Area;

			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
			if (centered) {
				if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda, iSpecies);
				if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda, iSpecies);
			}

			/*--- Viscous contribution ---*/
			if (viscous) {
				if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

				Gamma            = config->GetSpecies_Gamma(iSpecies);
				Gas_Constant     = config->GetSpecies_Gas_Constant(iSpecies);
				CharVibTemp      = config->GetCharVibTemp(iSpecies);
				Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity(iSpecies) + node[jPoint]->GetLaminarViscosity(iSpecies));
				Mean_ThermalConductivity     = 0.5*(node[iPoint]->GetThermalConductivity(iSpecies)
						+ node[jPoint]->GetThermalConductivity(iSpecies));
				Mean_ThermalConductivity_vib = 0.5*(node[iPoint]->GetThermalConductivity_vib(iSpecies)
						+ node[jPoint]->GetThermalConductivity_vib(iSpecies));
				Mean_Density     = 0.5*(node[iPoint]->GetSolution(loc+0) + node[jPoint]->GetSolution(loc+0));

				Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
				Lambda_2 = (Gamma-1.0) * Mean_ThermalConductivity / Gas_Constant;
				Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

				/*Mean_Tvib                    = 0.5*(node[iPoint]->GetPrimVar(iSpecies, nDim+1)
                                            + node[jPoint]->GetPrimVar(iSpecies, nDim+1));
				dTvdEv = Mean_Tvib*Mean_Tvib/(CharVibTemp*Mean_Density*Gas_Constant) * (exp(CharVibTemp/Mean_Tvib)-1.0)*(exp(CharVibTemp/Mean_Tvib)-1.0)/exp(CharVibTemp/Mean_Tvib);

				Lambda_1 = (4.0/3.0) * Mean_LaminarVisc * Area / (dij*Mean_Density);
				Lambda_2 = (Gamma-1.0) * Mean_ThermalConductivity * Area / (dij*Gas_Constant*Mean_Density);
				Lambda_3 = dTvdEv * Mean_ThermalConductivity_vib * Area / dij;*/

				/*--- Determine largest eigenvalue ---*/
				/*        Lambda = max(Lambda_1, Lambda_2);
        if (iSpecies < nDiatomics)
          Lambda = max(Lambda, Lambda_3);*/

				if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda, iSpecies);
				if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda, iSpecies);

			}
		}
	}

	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

			/*--- Compute closest normal neighbor ---*/
			jPoint = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

			Area = 0.0;
			dij = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Area += Normal[iDim]*Normal[iDim];
				dij  += pow(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim), 2.0);
			}
			Area = sqrt(Area);
			dij = sqrt(dij);

			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {

				/*--- Mean Values ---*/
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal,iSpecies);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed(iSpecies) * Area;

				Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
				if (geometry->node[iPoint]->GetDomain()) {

					/*--- Inviscid contribution ---*/
					node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
					if (centered) node[iPoint]->AddLambda(Lambda, iSpecies);

					/*--- Viscous contribution ---*/
					if (viscous) {
						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

						Gamma        = config->GetSpecies_Gamma(iSpecies);
						Gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
						CharVibTemp  = config->GetCharVibTemp(iSpecies);

						Mean_LaminarVisc             = 0.5*(node[iPoint]->GetLaminarViscosity(iSpecies)
								+ node[jPoint]->GetLaminarViscosity(iSpecies));
						Mean_ThermalConductivity     = 0.5*(node[iPoint]->GetThermalConductivity(iSpecies)
								+ node[jPoint]->GetThermalConductivity(iSpecies));
						Mean_ThermalConductivity_vib = 0.5*(node[iPoint]->GetThermalConductivity_vib(iSpecies)
								+ node[jPoint]->GetThermalConductivity_vib(iSpecies));
						Mean_Density                 = 0.5*(node[iPoint]->GetSolution(loc+0)
								+ node[jPoint]->GetSolution(loc+0));

						Mean_Tvib                    = 0.5*(node[iPoint]->GetPrimVar(iSpecies, nDim+1)
								+ node[jPoint]->GetPrimVar(iSpecies, nDim+1));
						dTvdEv = Mean_Tvib*Mean_Tvib/(CharVibTemp*Mean_Density*Gas_Constant) * (exp(CharVibTemp/Mean_Tvib)-1.0)*(exp(CharVibTemp/Mean_Tvib)-1.0)/exp(CharVibTemp/Mean_Tvib);

						Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
						//Lambda_2 = Gamma*Mean_LaminarVisc/Prandtl_Lam;
						Lambda_2 = (Gamma-1.0) * Mean_ThermalConductivity / Gas_Constant;
						Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

						/*Lambda_1 = (4.0/3.0) * Mean_LaminarVisc * Area / (dij*Mean_Density);
            Lambda_2 = (Gamma-1.0) * Mean_ThermalConductivity * Area / (dij*Gas_Constant*Mean_Density);
            Lambda_3 = dTvdEv * Mean_ThermalConductivity_vib * Area / dij;*/

						/*--- Determine largest eigenvalue ---*/
						/*            Lambda = max(Lambda_1, Lambda_2);
            if (iSpecies < nDiatomics)
              Lambda = max(Lambda, Lambda_3);*/

						node[iPoint]->AddMax_Lambda_Visc(Lambda, iSpecies);            

					}
				}
			}
		}
	}

	if (!MultipleTimeSteps) {

		/*--- Local time-stepping with all species advanced by the same time-step ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			dV = geometry->node[iPoint]->GetVolume();
			Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(0);
			for (iSpecies = 1; iSpecies < nSpecies; iSpecies++) {
				Lambda_iSpecies = max(Lambda_iSpecies, node[iPoint]->GetMax_Lambda_Inv(iSpecies));
			}
			Local_Delta_Time = config->GetCFL(iMesh)*dV / Lambda_iSpecies;

			/*--- Calculate viscous contribution ---*/
			if (viscous) {
				Lambda_Visc_iSpecies = node[iPoint]->GetMax_Lambda_Visc(0);
				for (iSpecies = 1; iSpecies < nSpecies; iSpecies++) {
					Lambda_Visc_iSpecies = max(Lambda_Visc_iSpecies, node[iPoint]->GetMax_Lambda_Visc(iSpecies));
				}
        Local_Delta_Time_Visc = config->GetCFL(iMesh)*dV / Lambda_Visc_iSpecies;
				//Local_Delta_Time_Visc =  config->GetCFL(iMesh)*K_v*dV*dV / Lambda_Visc_iSpecies;
				//Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
			}

			/*--- Set the time-step ---*/
			if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
				Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
			} else {
				node[iPoint]->SetDelta_Time(Local_Delta_Time);
			}

			/*--- Track maximum and minimum time-step values ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				Max_Delta_Time[iSpecies] = max(Max_Delta_Time[iSpecies], Local_Delta_Time);
				Min_Delta_Time[iSpecies] = min(Min_Delta_Time[iSpecies], Local_Delta_Time);
			}
		}
	}

	else {

		/*--- Local time-stepping with each species advanced at different speeds ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			dV = geometry->node[iPoint]->GetVolume();

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				CFL = config->GetCFL(iMesh, iSpecies);
				Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(iSpecies);
				Local_Delta_Time = CFL*dV / Lambda_iSpecies;

				/*--- Determine viscous contribution ---*/
				if (viscous) {
					Lambda_Visc_iSpecies = node[iPoint]->GetMax_Lambda_Visc(iSpecies);
					Local_Delta_Time_Visc = CFL*dV / Lambda_Visc_iSpecies;
					//Local_Delta_Time_Visc =  CFL*K_v*dV*dV / Lambda_Visc_iSpecies;
					//Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
				}

				/*--- Set the time step ---*/
				node[iPoint]->SetDelta_Time(Local_Delta_Time,iSpecies);

				/*--- Track maximum and minimum time-step values ---*/
				Max_Delta_Time[iSpecies] = max(Max_Delta_Time[iSpecies], Local_Delta_Time);
				Min_Delta_Time[iSpecies] = min(Min_Delta_Time[iSpecies], Local_Delta_Time);
			}
		}
	}

	/*--- For time-accurate solutions, use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}
}


void CPlasmaSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
		unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	double Volume_nM1, Volume_n, Volume_nP1, TimeStep;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

		/*--- Volume at time n-1 and n ---*/
		Volume_nM1 = geometry->node[iPoint]->GetVolume();
		Volume_n = geometry->node[iPoint]->GetVolume();
		Volume_nP1 = geometry->node[iPoint]->GetVolume();

		/*--- Time Step ---*/
		TimeStep = config->GetDelta_UnstTimeND();

		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
						+  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		}

		/*--- Add Residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);

		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;

				if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
					Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
				if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
					Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
			}
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
}



void CPlasmaSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) { 	unsigned long iEdge, iPoint, jPoint;

		bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
		bool high_order_diss = ((config->GetKind_Centered() == JST) && (iMesh == MESH_0));
		unsigned short iSpecies;

		/*--- Artificial dissipation preprocessing ---*/
		if (high_order_diss) {
			SetDissipation_Switch(geometry, config);
			SetUndivided_Laplacian(geometry, config);
		}

		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			/*--- Points in edge, set normal vectors, and number of neighbors ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
			numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());

			/*--- Set conservative variables w/o reconstruction ---*/
			numerics->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				numerics->SetPressure(node[iPoint]->GetPressure(iSpecies), node[jPoint]->GetPressure(iSpecies), iSpecies);
				numerics->SetSoundSpeed(node[iPoint]->GetSoundSpeed(iSpecies), node[jPoint]->GetSoundSpeed(iSpecies),iSpecies);
				numerics->SetEnthalpy(node[iPoint]->GetEnthalpy(iSpecies), node[jPoint]->GetEnthalpy(iSpecies),iSpecies);
				numerics->SetLambda(node[iPoint]->GetLambda(iSpecies), node[jPoint]->GetLambda(iSpecies),iSpecies);
			}


			if (high_order_diss) {
				numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					numerics->SetSensor(node[iPoint]->GetSensor(iSpecies), node[jPoint]->GetSensor(iSpecies), iSpecies);
				}
			}


			/*--- Compute residuals ---*/
			numerics->ComputeResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, config);

			/*--- Update convective and artificial dissipation residuals ---*/
			LinSysRes.AddBlock(iPoint, Res_Conv);
			LinSysRes.SubtractBlock(jPoint, Res_Conv);
      LinSysRes.AddBlock(iPoint, Res_Visc);
      LinSysRes.SubtractBlock(jPoint, Res_Visc);

			/*--- Set implicit stuff ---*/
			if (implicit) {
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
				Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
				Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);

			}
		}
}


void CPlasmaSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {

	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool high_order_diss = (( (config->GetKind_Upwind_Plasma() == ROE_2ND) || (config->GetKind_Upwind_Plasma() == SW_2ND) || (config->GetKind_Upwind_Plasma() == MSW_2ND) || (config->GetKind_Upwind_Plasma() == AUSM_2ND) || (config->GetKind_Upwind_Plasma() == HLLC_2ND)) && (iMesh == MESH_0));

	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Set conservative variables w/o reconstruction ---*/
		U_i = node[iPoint]->GetSolution();
		U_j = node[jPoint]->GetSolution();

		numerics->SetConservative(U_i, U_j);

		if (high_order_diss) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}

			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_j = node[jPoint]->GetLimiter(); Limiter_i = node[iPoint]->GetLimiter();
			}

			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (config->GetKind_SlopeLimit() == NONE) {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j;
				}
				else {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
			}
			numerics->SetConservative(Solution_i, Solution_j);
		}

		/*--- Compute the residual ---*/
		numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

		/*--- Update residual value ---*/
		LinSysRes.AddBlock(iPoint, Res_Conv);
		LinSysRes.SubtractBlock(jPoint, Res_Conv);

		/*--- Set implicit jacobians ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
		}
	}
}

void CPlasmaSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iVar, iEdge, iDim;
	unsigned short iSpecies, nPrimVar;
  
	nPrimVar = nDim+3;
  
	double min_value = 1.0;
  
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
  /*--- Compute gradients (not in the preprocessing, because we want to be sure that we have
   the gradient of the primitive varibles, not the conservative) ---*/
  for (iEdge = 0; iEdge <  geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (iVar = 0 ; iVar < nPrimVar; iVar++) {
        PrimGradLimiter_i[iSpecies][iVar] = node[iPoint]->GetLimiterPrimitive(iSpecies, iVar);
        PrimGradLimiter_j[iSpecies][iVar] = node[jPoint]->GetLimiterPrimitive(iSpecies, iVar);
        for (iDim = 0; iDim < nDim; iDim++) {
          PrimGrad_i[iSpecies][iVar][iDim] = node[iPoint]->GetGradient_Primitive(iSpecies, iVar, iDim);
          PrimGrad_j[iSpecies][iVar][iDim] = node[jPoint]->GetGradient_Primitive(iSpecies, iVar, iDim);
          
          /*--- Apply the limiter directly ---*/
          PrimGrad_i[iSpecies][iVar][iDim] = PrimGrad_i[iSpecies][iVar][iDim]*PrimGradLimiter_i[iSpecies][iVar];
          PrimGrad_j[iSpecies][iVar][iDim] = PrimGrad_j[iSpecies][iVar][iDim]*PrimGradLimiter_j[iSpecies][iVar];
          
          
          if (PrimGradLimiter_i[iSpecies][iVar] < min_value) min_value = PrimGradLimiter_i[iSpecies][iVar];
          if (PrimGradLimiter_j[iSpecies][iVar] < min_value) min_value = PrimGradLimiter_j[iSpecies][iVar];
          
        }
      }
    }
    
    /*--- Acquire primitive variables and gradients from the CPlasmaVariable class ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimVar_Plasma(), node[jPoint]->GetPrimVar_Plasma());
    numerics->SetPrimVarGradient(PrimGrad_i, PrimGrad_j);
    
    /*--- Pass transport coefficients from the variable to the numerics class ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
      numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(iSpecies), node[jPoint]->GetLaminarViscosity(iSpecies), iSpecies);
      numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(iSpecies), node[jPoint]->GetEddyViscosity(iSpecies), iSpecies);
      numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(iSpecies), node[jPoint]->GetThermalConductivity(iSpecies), iSpecies);
      numerics->SetThermalConductivity_vib(node[iPoint]->GetThermalConductivity_vib(iSpecies), node[jPoint]->GetThermalConductivity_vib(iSpecies), iSpecies);
    }
    
    /*--- Compute and update residual ---*/
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.SubtractBlock(iPoint, Res_Visc);
    LinSysRes.AddBlock(jPoint, Res_Visc);
    
    /*--- Update Jacobians for implicit calculations ---*/
    if (implicit) {
      Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
      Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
      Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
      Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
    }
  }
  
}

void CPlasmaSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
		CConfig *config, unsigned short iMesh) {
	unsigned short iVar, jVar, iSpecies, iDim;
	unsigned long iPoint;
	double vol;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool axisymmetric = config->GetAxisymmetric();
	if(magnet) {
		for (unsigned short iDim = 0; iDim < nDim; iDim ++)
			Mag_Force[iDim] = 0.0;
		vol = 1.0;
	}

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Initialize residual vectors ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			Residual[iVar] = 0.0;
			Residual_Chemistry[iVar] = 0.0;
			Residual_MomentumExch[iVar] = 0.0;
			Residual_ElecForce[iVar] = 0.0;
			Residual_EnergyExch[iVar] = 0.0;
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_Chemistry[iVar][jVar] = 0.0;
				Jacobian_MomentumExch[iVar][jVar] = 0.0;
				Jacobian_ElecForce[iVar][jVar] = 0.0;
				Jacobian_EnergyExch[iVar][jVar] = 0.0;
			}
		}
		if (axisymmetric) {
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual_Axisymmetric[iVar] = 0.0;
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_Axisymmetric[iVar][jVar] = 0.0;
			}
		}

		if (config->GetElectricSolver()) {
			numerics->SetElecField(node[iPoint]->GetElectricField());
		}

		/*--- Set y coordinate ---*/
		numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());

		/*--- Set solution  ---*/
		numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

		/*--- Set control volume ---*/
		numerics->SetVolume(geometry->node[iPoint]->GetVolume());

		//        double time = node[iPoint]->GetDelta_Time();
		//	numerics->SetTimeStep(node[iPoint]->GetDelta_Time());

		//        node[iPoint]->GetDelta_Time()
		//        [ELEC_SOL]->node[iPoint]->SetChargeDensity

		//     node[iPoint]->SetDelta_Time(node[iPoint]->GetDelta_Time());

		/*--- Acquire primitive variables and gradients from the CPlasmaVariable class ---*/
		numerics->SetPrimitive(node[iPoint]->GetPrimVar_Plasma(), node[iPoint]->GetPrimVar_Plasma());
		numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive_Plasma(),node[iPoint]->GetGradient_Primitive_Plasma());

		/*--- Compute Residual ---*/
		if (config->GetKind_GasModel() == ARGON) {
			numerics->ComputeResidual(Residual, Jacobian_i, config);

			if (magnet && nDim == 3) {
				node[iPoint]->SetMagneticField(numerics->GetMagneticField());
				vol = geometry->node[iPoint]->GetVolume();
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					for ( iDim = 0; iDim < nDim; iDim ++)
						Mag_Force[iDim] += numerics->GetMagneticForce(iSpecies, iDim)*vol;
				}
			}

			LinSysRes.SubtractBlock(iPoint, Residual);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
		}

		else if ((config->GetKind_GasModel() == ARGON_SID) || config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2
				|| config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {

			/*--- Axisymmetric source terms ---*/
			if (axisymmetric) {
				numerics->ComputeResidual_Axisymmetric(Residual_Axisymmetric, config);
				numerics->SetJacobian_Axisymmetric(Jacobian_Axisymmetric, config);
				LinSysRes.AddBlock(iPoint, Residual_Axisymmetric);
				if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_Axisymmetric);
			}

			/*--- Chemistry source terms ---*/
			numerics->ComputeResidual_Chemistry(Residual_Chemistry, config);
			numerics->SetJacobian_Chemistry(Jacobian_Chemistry, config);
			LinSysRes.SubtractBlock(iPoint, Residual_Chemistry);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_Chemistry);

			/*--- Electric force source terms ---*/
			/*		numerics->ComputeResidual_ElecForce(Residual_ElecForce, config);
			 numerics->SetJacobian_ElecForce(Jacobian_ElecForce, config);
			 LinSysRes.SubtractBlock(iPoint, Residual_ElecForce);
			 if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ElecForce);*/

			/*--- Momentum exchange source terms ---*/
			numerics->ComputeResidual_MomentumExch(Residual_MomentumExch, config);
			numerics->SetJacobian_MomentumExch(Jacobian_MomentumExch, config);
			LinSysRes.SubtractBlock(iPoint, Residual_MomentumExch);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_MomentumExch);

			/*--- Energy exchange source terms ---*/
			numerics->ComputeResidual_EnergyExch(Residual_EnergyExch, config);
			numerics->SetJacobian_EnergyExch(Jacobian_EnergyExch, config);
			LinSysRes.SubtractBlock(iPoint, Residual_EnergyExch);
			if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_EnergyExch);
		}
	}

}

void CPlasmaSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {

}

/*!
 * \method Copy_Zone_Solution
 * \brief Copy solution from solver 1 into solver 2
 * \author A. Lonkar
 */
void CPlasmaSolver::Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config,
		CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {
	unsigned long iPoint;

	double positive_charge, negative_charge;
	unsigned short nVar_Ion, nVar_Electron; //, iSpecies, iDim, loc, iVar;
	//	double **Gradient;

	nVar_Ion = 1*(nDim+2);
	nVar_Electron = 2*(nDim+2);
	for (iPoint = 0; iPoint < solver1_geometry[MESH_0]->GetnPointDomain(); iPoint++) {
		positive_charge = solver1_solution[MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(nVar_Ion);
		negative_charge = solver1_solution[MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(nVar_Electron);
		solver2_solution[MESH_0][ELEC_SOL]->node[iPoint]->SetChargeDensity(positive_charge, negative_charge);
		solver2_solution[MESH_0][ELEC_SOL]->node[iPoint]->SetPlasmaTimeStep(solver1_solution[MESH_0][PLASMA_SOL]->node[iPoint]->GetDelta_Time());

		/*	if (solver1_config->GetElectricSolver()) {
     Gradient = solver1_solution[MESH_0][PLASMA_SOL]->node[iPoint]->GetGradient();
     for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
     if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
     else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
     for (iDim = 0; iDim < nDim; iDim++) {
     iVar = loc + iDim + 1;
     solver2_solution[MESH_0][ELEC_SOL]->node[iPoint]->SetPlasmaRhoUGradient(iSpecies, Gradient[iVar][iDim], iDim);
     }
     }
     solver2_solution[MESH_0][ELEC_SOL]->node[iPoint]->SetPlasmaTimeStep(node[iPoint]->GetDelta_Time());
     } */
	}

}


void CPlasmaSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

	unsigned long iPoint, jPoint, iEdge, Point_Normal = 0, iVertex;
	double Pressure_i = 0, Pressure_j = 0, *Normal;
	unsigned short iVar, iMarker, iSpecies, loc, nVar_Species;
	double *Diff = new double[nVar];
	double *U_halo = new double[nVar];

	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetUnd_LaplZero();

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		/*--- Solution differences ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);

		/*--- Correction for compressible flows which use the enthalpy ---*/
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Pressure_i = node[iPoint]->GetPressure(iSpecies);
			Pressure_j = node[jPoint]->GetPressure(iSpecies);
			if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
			else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			nVar_Species = loc + (nDim+2);
			Diff[nVar_Species-1] = (node[iPoint]->GetSolution(nVar_Species-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar_Species-1) + Pressure_j);
		}

		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

	}

	/*--- Loop over all boundaries and include an extra contribution
   from a halo node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
				config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY )
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

					Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

					/*--- Interpolate & compute difference in the conserved variables ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						if ((config->GetMarker_All_Boundary(iMarker) == EULER_WALL) ||
								(config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
								(config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL)) {
							U_halo[iVar] = node[Point_Normal]->GetSolution(iVar);
						} else {
							U_halo[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
						}
						Diff[iVar]   = node[iPoint]->GetSolution(iVar) - U_halo[iVar];
					}
					/*--- Correction for compressible flows ---*/
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
						Pressure_i = node[iPoint]->GetPressure(iSpecies);
						Pressure_j = 2.0*node[iPoint]->GetPressure(iSpecies) - node[Point_Normal]->GetPressure(iSpecies);
						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
						nVar_Species = loc + (nDim+2);
						Diff[nVar_Species-1] = (node[iPoint]->GetSolution(nVar_Species-1) + Pressure_i) - (U_halo[nVar_Species-1] + Pressure_j);
					}
					/*--- Subtract contribution at the boundary node only ---*/
					if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
				}
			}

	delete [] Diff;
	delete [] U_halo;
}

void CPlasmaSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {

	unsigned long iVertex, iPoint;
	unsigned short iDim, iMarker, Boundary, Monitoring, iSpecies;
	double Pressure_allspecies,Pressure_allspecies_Inf, *Normal = NULL, dist[3], *Coord, Face_Area, PressInviscid;
	double factor, NFPressOF, RefVel2, RefDensity, RefPressure;

	double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
	double Beta            = config->GetAoS()*PI_NUMBER/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment = config->GetRefLengthMoment();
	double *Origin         = config->GetRefOriginMoment();

	RefVel2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		RefVel2  += Velocity_Inf_mean[iDim]*Velocity_Inf_mean[iDim];

	RefDensity = 0.0; RefPressure = 0.0; Pressure_allspecies_Inf = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		RefDensity   += Density_Inf[iSpecies];
		RefPressure  += Pressure_Inf[iSpecies];
		Pressure_allspecies_Inf += Pressure_Inf[iSpecies];
	}

	/*-- Initialization ---*/
	Total_CDrag = 0.0; Total_CLift = 0.0; Total_CSideForce = 0.0;
	Total_CMx = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
	Total_CFx = 0.0; Total_CFy = 0.0; Total_CFz = 0.0;
	Total_CEff = 0.0; Total_CMerit = 0.0; Total_CNearFieldOF = 0.0;
	Total_CT = 0.0; Total_CQ = 0.0; Total_Q = 0.0;
	AllBound_CDrag_Inv = 0.0; AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0;
	AllBound_CMx_Inv = 0.0; AllBound_CMy_Inv = 0.0; AllBound_CMz_Inv = 0.0;
	AllBound_CFx_Inv = 0.0; AllBound_CFy_Inv = 0.0; AllBound_CFz_Inv = 0.0;
	AllBound_CEff_Inv = 0.0; AllBound_CMerit_Inv = 0.0;
	AllBound_CNearFieldOF_Inv = 0.0;
	AllBound_CT_Inv = 0.0; AllBound_CQ_Inv = 0.0;
	PressureDrag = 0.0;  MagnetDrag = 0.0;
	/*--- Loop over the Euler and Navier-Stokes markers ---*/
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary   = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);

		if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {

			for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;

			MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;

			NFPressOF = 0.0; PressInviscid = 0.0;

			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				//                for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
				//                    for (iDim = 0; iDim < nDim; iDim++)
				//                        CPressForce[iMarker][iSpecies][iDim][iVertex] = 0.0;

				Pressure_allspecies = 0.0;
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
					Pressure_allspecies += node[iPoint]->GetPressure(iSpecies);

				CPressure[iMarker][iVertex] = (Pressure_allspecies - RefPressure)*factor*RefAreaCoeff;

				/*--- Note that the pressure coefficient is computed at the
				 halo cells (for visualization purposes), but not the forces ---*/
				if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Coord = geometry->node[iPoint]->GetCoord();

					Face_Area = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						/*--- Total force, distance computation, and face area ---*/
						ForceInviscid[iDim] -= Pressure_allspecies*Normal[iDim]*factor;
						dist[iDim] = Coord[iDim] - Origin[iDim];
						Face_Area += Normal[iDim]*Normal[iDim];
					}
					Face_Area = sqrt(Face_Area);

					for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
						for (iDim = 0; iDim < nDim; iDim++) {
							CPressForce[iMarker][iSpecies][iDim][iVertex] = -1.0*node[iPoint]->GetPressure(iSpecies)*Normal[iDim]*factor/(Face_Area);
						}
					}

					/*--- Quadratic objective function for the near field.
           This uses the infinity pressure regardless of Mach number. ---*/
					NFPressOF += 0.5*(Pressure_allspecies - Pressure_allspecies_Inf)*(Pressure_allspecies - Pressure_allspecies_Inf)*Normal[nDim-1];
					PressInviscid += CPressure[iMarker][iVertex]*Face_Area;

					/*--- Moment with respect to the reference axis ---*/
					if (iDim == 3) {
						MomentInviscid[0] -= Pressure_allspecies*(Normal[2]*dist[1]-Normal[1]*dist[2])*factor/RefLengthMoment;
						MomentInviscid[1] -= Pressure_allspecies*(Normal[0]*dist[2]-Normal[2]*dist[0])*factor/RefLengthMoment;
					}
					MomentInviscid[2] -= Pressure_allspecies*(Normal[1]*dist[0]-Normal[0]*dist[1])*factor/RefLengthMoment;
				}
			}

			/*--- Transform ForceInviscid and MomentInviscid into non-dimensional coefficient ---*/
			if  (Monitoring == YES) {
				if (nDim == 2) {
					if (Boundary != NEARFIELD_BOUNDARY) {
						CDrag_Inv[iMarker] =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
						CLift_Inv[iMarker] = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
						CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0;
						CMy_Inv[iMarker] = 0.0;
						CMz_Inv[iMarker] = MomentInviscid[2];
						CEff_Inv[iMarker] = CLift_Inv[iMarker]/(CDrag_Inv[iMarker]+config->GetCteViscDrag());
						CNearFieldOF_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = ForceInviscid[0];
						CFy_Inv[iMarker] = ForceInviscid[1];
						CFz_Inv[iMarker] = 0.0;
						CT_Inv[iMarker] = -CFx_Inv[iMarker];
						CQ_Inv[iMarker] = 0.0;
						CMerit_Inv[iMarker] = 0.0;
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = 0.0; CFy_Inv[iMarker] = 0.0; CFz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
						CNearFieldOF_Inv[iMarker] = NFPressOF;
						CT_Inv[iMarker] = 0.0;
						CQ_Inv[iMarker] = 0.0;
						CMerit_Inv[iMarker] = 0.0;
					}

				}
				if (nDim == 3) {
					if (Boundary != NEARFIELD_BOUNDARY) {
						CDrag_Inv[iMarker] =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
						CLift_Inv[iMarker] = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
						CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
						CMx_Inv[iMarker] = MomentInviscid[0];
						CMy_Inv[iMarker] = MomentInviscid[1];
						CMz_Inv[iMarker] = MomentInviscid[2];
						CEff_Inv[iMarker] = CLift_Inv[iMarker]/(CDrag_Inv[iMarker]+config->GetCteViscDrag());
						CNearFieldOF_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = ForceInviscid[0];
						CFy_Inv[iMarker] = ForceInviscid[1];
						CFz_Inv[iMarker] = ForceInviscid[2];
						CT_Inv[iMarker] = -CFx_Inv[iMarker];
						CQ_Inv[iMarker] = -CMx_Inv[iMarker];
						/*--- For now, define the FM as a simple rotor efficiency ---*/
						//CMerit_Inv[iMarker] = CT_Inv[iMarker]*sqrt(fabs(CT_Inv[iMarker]))/(sqrt(2.0)*CQ_Inv[iMarker]);
						CMerit_Inv[iMarker] = CT_Inv[iMarker]/CQ_Inv[iMarker];
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = 0.0; CFy_Inv[iMarker] = 0.0; CFz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
						CNearFieldOF_Inv[iMarker] = NFPressOF;
						CT_Inv[iMarker] = 0.0;
						CQ_Inv[iMarker] = 0.0;
						CMerit_Inv[iMarker] = 0.0;
					}
				}

				AllBound_CDrag_Inv += CDrag_Inv[iMarker];
				AllBound_CLift_Inv += CLift_Inv[iMarker];
				AllBound_CSideForce_Inv += CSideForce_Inv[iMarker];
				AllBound_CMx_Inv += CMx_Inv[iMarker];
				AllBound_CMy_Inv += CMy_Inv[iMarker];
				AllBound_CMz_Inv += CMz_Inv[iMarker];
				AllBound_CEff_Inv += CEff_Inv[iMarker];
				AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
				AllBound_CFx_Inv += CFx_Inv[iMarker];
				AllBound_CFy_Inv += CFy_Inv[iMarker];
				AllBound_CFz_Inv += CFz_Inv[iMarker];
				AllBound_CT_Inv += CT_Inv[iMarker];
				AllBound_CQ_Inv += CQ_Inv[iMarker];
				AllBound_CMerit_Inv += CMerit_Inv[iMarker];
			}
		}
	}

	/*-- Pressure Drag Calcualtion --- */
	PressureDrag = AllBound_CDrag_Inv;

	MagnetDrag = 0.0;
	if (magnet && nDim ==3) {
		double Cd_mag, Cl_mag, Csd_mag;
		Cd_mag =  factor*(Mag_Force[0]*cos(Alpha)*cos(Beta) + Mag_Force[1]*sin(Beta) + Mag_Force[2]*sin(Alpha)*cos(Beta));
		Cl_mag = factor*(-Mag_Force[0]*sin(Alpha) + Mag_Force[2]*cos(Alpha));
		Csd_mag = factor*(-Mag_Force[0]*sin(Beta)*cos(Alpha) + Mag_Force[1]*cos(Beta) - Mag_Force[2]*sin(Beta)*sin(Alpha));

		Total_CDrag -= Cd_mag;
		Total_CLift -= Cl_mag;
		Total_CSideForce -= Csd_mag;
		MagnetDrag  -= Cd_mag;
		//        cout << "Mag_Force[0] = " << Mag_Force[0] << "Mag_Force[1] = " << Mag_Force[1] << endl;
		//        cout << "PressureDrag = " << PressureDrag << endl;
	}

	Total_CDrag += AllBound_CDrag_Inv;
	Total_CLift += AllBound_CLift_Inv;
	Total_CSideForce += AllBound_CSideForce_Inv;
	Total_CMx += AllBound_CMx_Inv;
	Total_CMy += AllBound_CMy_Inv;
	Total_CMz += AllBound_CMz_Inv;
	Total_CEff += AllBound_CEff_Inv;
	Total_CNearFieldOF += AllBound_CNearFieldOF_Inv;
	Total_CFx += AllBound_CFx_Inv;
	Total_CFy += AllBound_CFy_Inv;
	Total_CFz += AllBound_CFz_Inv;
	Total_CT += AllBound_CT_Inv;
	Total_CQ += AllBound_CQ_Inv;
	Total_CMerit += AllBound_CMerit_Inv;




}

void CPlasmaSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint, iSpecies;
	unsigned short Boundary, Monitoring, iMarker, iDim, jDim;
	double **Tau, Delta, Viscosity, **Grad_PrimVar, div_vel, *Normal, *TauElem;
	double dist[3], *Coord, *UnitaryNormal, *TauTangent, Area, WallShearStress, TauNormal;
	double GradTemperature, GradTemperature_vib;
	double ThermalConductivity, ThermalConductivity_vib;
	double factor, RefVel2, RefDensity;
	double HeatLoad;

	double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double Beta         = config->GetAoS()*PI_NUMBER/180.0;
	double *Origin      = config->GetRefOriginMoment();

	//	double Gas_Constant = config->GetGas_ConstantND();
	//	cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
	RefVel2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		RefVel2  += Velocity_Inf_mean[iDim]*Velocity_Inf_mean[iDim];

	RefDensity = 0.0;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		RefDensity   += Density_Inf[iSpecies];
	}
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);


	/*-- Initialization --*/
	AllBound_CDrag_Visc = 0.0; AllBound_CLift_Visc = 0.0;
	AllBound_CMx_Visc = 0.0; AllBound_CMy_Visc = 0.0; AllBound_CMz_Visc = 0.0;
	AllBound_CFx_Visc = 0.0; AllBound_CFy_Visc = 0.0; AllBound_CFz_Visc = 0.0;
	AllBound_CEff_Visc = 0.0; AllBound_CMerit_Visc = 0.0;
	AllBound_CT_Visc = 0.0; AllBound_Q_Visc = 0.0; ViscDrag = 0.0;

	/*--- Vector and variables initialization ---*/
	UnitaryNormal      = new double [nDim];
	TauElem    = new double [nDim];
	TauTangent = new double [nDim];
	Tau        = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim]   = new double [nDim];

	/*--- Loop over the Navier-Stokes markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);

		if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {

			for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
			MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
			HeatLoad = 0.0;

			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {

					if (config->GetKind_GasModel() != ARGON) CHeatTransfer[iMarker][iSpecies][iVertex] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						CViscForce[iMarker][iSpecies][iDim][iVertex] = 0.0;
					}
				}

				CSkinFriction[iMarker][iVertex] = 0.0;


				if (geometry->node[iPoint]->GetDomain()) {

					/*-- Compute geometry parameters ---*/
					Coord = geometry->node[iPoint]->GetCoord();
					for (iDim = 0; iDim < nDim; iDim++) dist[iDim] = Coord[iDim] - Origin[iDim];
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
					for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = Normal[iDim]/Area;


					for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
						/*--- Initialize Tau ---*/
						for (iDim = 0; iDim < nDim; iDim++)
							for (jDim = 0 ; jDim < nDim; jDim++)
								Tau[iDim][jDim] = 0.0;

						/*--- Get the species primitive variables ---*/
						Grad_PrimVar            = node[iPoint]->GetGradient_Primitive(iSpecies);
						Viscosity               = node[iPoint]->GetLaminarViscosity(iSpecies);
						ThermalConductivity     = node[iPoint]->GetThermalConductivity(iSpecies);
						ThermalConductivity_vib = node[iPoint]->GetThermalConductivity_vib(iSpecies);

						div_vel = 0.0;
						for (iDim = 0; iDim < nDim; iDim++)
							div_vel += Grad_PrimVar[iDim+1][iDim];

						/*--- Compute tau for each species ---*/
						for (iDim = 0; iDim < nDim; iDim++)
							for (jDim = 0 ; jDim < nDim; jDim++) {
								Delta = 0.0; if (iDim == jDim) Delta = 1.0;
								Tau[iDim][jDim] += Viscosity*(Grad_PrimVar[jDim+1][iDim] +
										Grad_PrimVar[iDim+1][jDim]) - TWO3*Viscosity*div_vel*Delta;
							}

						/*--- Project tau in the normal direction --*/
						for (iDim = 0; iDim < nDim; iDim++) {
							TauElem[iDim] = 0.0;
							for (jDim = 0; jDim < nDim; jDim++)
								TauElem[iDim] += Tau[iDim][jDim]*UnitaryNormal[jDim] ;
						}

						/*--- Compute viscous forces, and moment ---*/
						if (Monitoring == YES) {
							for (iDim = 0; iDim < nDim; iDim++)
								ForceViscous[iDim] += TauElem[iDim]*Area*factor;
						}

						for (iDim = 0; iDim < nDim; iDim++)
							CViscForce[iMarker][iSpecies][iDim][iVertex] = TauElem[iDim]*factor;

						/*--- Calculate the heat flux at the wall ---*/
						GradTemperature = 0.0;
						GradTemperature_vib = 0.0;
						for (iDim = 0; iDim < nDim; iDim++){
							GradTemperature +=  Grad_PrimVar[0][iDim]*(-Normal[iDim]);
							GradTemperature_vib += Grad_PrimVar[nDim+1][iDim]*(-Normal[iDim]);
						}
						if (config->GetKind_GasModel() != ARGON) {
							CHeatTransfer[iMarker][iSpecies][iVertex] = (ThermalConductivity*GradTemperature + ThermalConductivity_vib*GradTemperature_vib)/(0.5*RefDensity*RefVel2);
              HeatLoad += CHeatTransfer[iMarker][iSpecies][iVertex];
						}
					}

					/*--- Compute wall shear stress, skin friction coefficient, and heat flux on the wall ---*/
					TauNormal = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						TauNormal += TauElem[iDim] * UnitaryNormal[iDim];

					for (iDim = 0; iDim < nDim; iDim++)
						TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitaryNormal[iDim];

					WallShearStress = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						WallShearStress += TauTangent[iDim]*TauTangent[iDim];
					WallShearStress = sqrt(WallShearStress);

					/*--- Note that the wall Shear Stress is just mu(delta u/delta y)---*/
					CSkinFriction[iMarker][iVertex] += WallShearStress / (0.5*RefDensity*RefVel2);

				}

			}

			/*--- Transform ForceInviscid into CLift and CDrag ---*/
			if  (Monitoring == YES) {
				if (nDim == 2) {
					CDrag_Visc[iMarker] =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
					CLift_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
					CMx_Visc[iMarker] = 0.0;
					CMy_Visc[iMarker] = 0.0;
					CMz_Visc[iMarker] = MomentViscous[2];
					CEff_Visc[iMarker] = CLift_Visc[iMarker]/CDrag_Visc[iMarker];
					CFx_Visc[iMarker] = ForceViscous[0];
					CFy_Visc[iMarker] = ForceViscous[1];
					CFz_Visc[iMarker] = 0.0;
					CQ_Visc[iMarker] = -CMz_Visc[iMarker];
          Q_Visc[iMarker] = HeatLoad;
				}
				if (nDim == 3) {
					CDrag_Visc[iMarker] =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
					CLift_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
					CMx_Visc[iMarker] = MomentViscous[0];
					CMy_Visc[iMarker] = MomentViscous[1];
					CMz_Visc[iMarker] = MomentViscous[2];
					CEff_Visc[iMarker] = CLift_Visc[iMarker]/CDrag_Visc[iMarker];
					CFx_Visc[iMarker] = ForceViscous[0];
					CFy_Visc[iMarker] = ForceViscous[1];
					CFz_Visc[iMarker] = ForceViscous[2];
					CQ_Visc[iMarker] = -CMz_Visc[iMarker];
          Q_Visc[iMarker] = HeatLoad;
				}

				AllBound_CDrag_Visc += CDrag_Visc[iMarker];
				AllBound_CLift_Visc += CLift_Visc[iMarker];
				AllBound_CMx_Visc += CMx_Visc[iMarker];
				AllBound_CMy_Visc += CMy_Visc[iMarker];
				AllBound_CMz_Visc += CMz_Visc[iMarker];
				AllBound_CEff_Visc += CEff_Visc[iMarker];
				AllBound_CFx_Visc += CFx_Visc[iMarker];
				AllBound_CFy_Visc += CFy_Visc[iMarker];
				AllBound_CFz_Visc += CFz_Visc[iMarker];
				AllBound_Q_Visc += Q_Visc[iMarker];
			}
		}
	}
	ViscDrag = AllBound_CDrag_Visc;
	Total_CDrag += AllBound_CDrag_Visc;
	Total_CLift += AllBound_CLift_Visc;
	Total_CMx += AllBound_CMx_Visc;
	Total_CMy += AllBound_CMy_Visc;
	Total_CMz += AllBound_CMz_Visc;
	Total_CEff = Total_CLift/Total_CDrag;
	Total_CFx += AllBound_CFx_Visc;
	Total_CFy += AllBound_CFy_Visc;
	Total_CFz += AllBound_CFz_Visc;
	Total_Q += AllBound_Q_Visc;

	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
	delete [] UnitaryNormal;
	delete [] TauTangent;
	delete [] TauElem;
}


void CPlasmaSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {

	unsigned long iEdge, iPoint, jPoint;
	double Pressure_i, Pressure_j;
	unsigned short iSpecies;
	/*--- Reset variables to store the undivided pressure ---*/
	for (iSpecies =0; iSpecies < nSpecies; iSpecies ++ ) {
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			p1_Und_Lapl[iSpecies][iPoint] = 0.0;
			p2_Und_Lapl[iSpecies][iPoint] = 0.0;
		}
	}
	/*--- Evaluate the pressure sensor ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {

			Pressure_i = node[iPoint]->GetPressure(iSpecies);
			Pressure_j = node[jPoint]->GetPressure(iSpecies);

			if (geometry->node[iPoint]->GetDomain()) p1_Und_Lapl[iSpecies][iPoint] += (Pressure_j - Pressure_i);
			if (geometry->node[jPoint]->GetDomain()) p1_Und_Lapl[iSpecies][jPoint] += (Pressure_i - Pressure_j);

			if (geometry->node[iPoint]->GetDomain()) p2_Und_Lapl[iSpecies][iPoint] += (Pressure_i + Pressure_j);
			if (geometry->node[jPoint]->GetDomain()) p2_Und_Lapl[iSpecies][jPoint] += (Pressure_i + Pressure_j);
		}
	}
	/*--- Set pressure switch for each point ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
			node[iPoint]->SetSensor(fabs(p1_Und_Lapl[iSpecies][iPoint]) / p2_Und_Lapl[iSpecies][iPoint], iSpecies);

		}
	}
}

void CPlasmaSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned short iVar, iSpecies, loc, nVar_Species;
	unsigned long iPoint, total_index, IterLinSol = 0;
	double Delta, *local_Res_TruncError, Vol;
	bool MultipleTimeSteps = (config->MultipleTimeSteps());
    
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
		SetRes_Max(iVar, 0.0, 0);
	}
    
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
		local_Res_TruncError = node[iPoint]->GetResTruncError();
		Vol = geometry->node[iPoint]->GetVolume();
        
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
			LinSysRes[total_index] = -(LinSysRes[total_index] + local_Res_TruncError[iVar]);
			LinSysSol[total_index] = 0.0;
			AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
			AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
		}
        
		if (roe_turkel) {
			SetPreconditioner(config, iPoint);
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				if (MultipleTimeSteps)  Delta = Vol/(node[iPoint]->GetDelta_Time(iSpecies));
				else Delta = Vol/(node[iPoint]->GetDelta_Time());
				if ( iSpecies < nDiatomics ) { loc = (nDim+3)*iSpecies; nVar_Species = (nDim+3);}
				else { loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics); nVar_Species = (nDim+2);}
                
				for (iVar = 0; iVar < nVar_Species; iVar ++ )
					for (unsigned short jVar = 0; jVar < nVar_Species; jVar ++ )
						Precon_Mat_inv[loc + iVar][loc + jVar] = Delta*Precon_Mat_inv[loc + iVar][loc + jVar];
			}
			Jacobian.AddBlock(iPoint, iPoint, Precon_Mat_inv);
		}
		else {
            
			if (!MultipleTimeSteps) {
				Delta = Vol/(node[iPoint]->GetDelta_Time());
				Jacobian.AddVal2Diag(iPoint,Delta);
			} else {
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					Species_Delta[iSpecies] = Vol/(node[iPoint]->GetDelta_Time(iSpecies));
				}
				//Jacobian.AddVal2Diag(iPoint, Species_Delta, nDim);
				Jacobian.AddVal2Diag(iPoint, Species_Delta, nDim, nDiatomics);
			}
		}
	}
    
	/*--- Initialize residual and solution at the ghost points ---*/
	for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar + iVar;
			LinSysRes[total_index] = 0.0;
			LinSysSol[total_index] = 0.0;
		}
	}
    
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    IterLinSol = system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                                config->GetLinear_Solver_Iter(), true);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    IterLinSol = system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                              config->GetLinear_Solver_Iter(), true);
  
  /*--- The the number of iterations of the linear solver ---*/
  SetIterLinSolver(IterLinSol);
  
  /*--- dealocate memory ---*/
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
		}
	}
    
	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

	/*--- Compute the root mean square residual ---*/
	SetResidual_RMS(geometry, config);
    
}


void CPlasmaSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iVar, iMarker, iSpecies;
	double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
	Partial_Gradient, Partial_Res, *Normal;

	unsigned short nPrimVar = nDim+3;

	PrimVar_i = new double [nPrimVar];
	PrimVar_j = new double [nPrimVar];
	PrimVar_Vertex = new double [nPrimVar];

	/*--- Set Gradient_Primitive to zero ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero(nPrimVar);

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

		/*--- Loop interior edges ---*/
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);

			for (iVar = 0; iVar < nPrimVar; iVar++) {
				PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iSpecies, iVar);
				PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iSpecies, iVar);
			}

			Normal = geometry->edge[iEdge]->GetNormal();
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
				for (iDim = 0; iDim < nDim; iDim++) {
					Partial_Res = PrimVar_Average*Normal[iDim];
					if (geometry->node[iPoint]->GetDomain())
						node[iPoint]->AddGradient_Primitive(iSpecies, iVar, iDim, Partial_Res);
					if (geometry->node[jPoint]->GetDomain())
						node[jPoint]->SubtractGradient_Primitive(iSpecies, iVar, iDim, Partial_Res);
				}
			}
		}

		/*--- Loop boundary edges ---*/
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {

					for (iVar = 0; iVar < nPrimVar; iVar++)
						PrimVar_Vertex[iVar] = node[iPoint]->GetPrimVar(iSpecies, iVar);

					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					for (iVar = 0; iVar < nPrimVar; iVar++) {
						for (iDim = 0; iDim < nDim; iDim++) {
							Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
							node[iPoint]->SubtractGradient_Primitive(iSpecies, iVar, iDim, Partial_Res);
						}
					}
				}
			}
		}

		/*--- Update gradient value ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				for (iDim = 0; iDim < nDim; iDim++) {
					Partial_Gradient = node[iPoint]->GetGradient_Primitive(iSpecies, iVar, iDim) / geometry->node[iPoint]->GetVolume();
					node[iPoint]->SetGradient_Primitive(iSpecies, iVar, iDim, Partial_Gradient);
				}
			}
		}
	}

	delete [] PrimVar_Vertex;
	delete [] PrimVar_i;
	delete [] PrimVar_j;

	Set_MPI_PrimVar_Gradient(geometry, config);
}

void CPlasmaSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config, unsigned long iVertex, unsigned short val_marker, double *val_PrimVar_i) {
	unsigned long iPoint;
	unsigned short iDim, iVar, iSpecies;
	double *PrimVar_Vertex, Partial_Gradient, Partial_Res, *Normal;

	unsigned short nPrimVar = nDim+3;

	PrimVar_Vertex = new double [nPrimVar];

	/*--- Set Gradient_Primitive to zero ---*/
	iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
	node[iPoint]->SetGradient_PrimitiveZero(nPrimVar);

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

		if (geometry->node[iPoint]->GetDomain()) {
			for (iVar = 0; iVar < nPrimVar; iVar++)
				PrimVar_Vertex[iVar] = node[iPoint]->GetPrimVar(iSpecies, iVar);

			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				for (iDim = 0; iDim < nDim; iDim++) {
					Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
					node[iPoint]->SubtractGradient_Primitive(iSpecies, iVar, iDim, Partial_Res);
				}
			}
		}

		/*--- Update gradient value ---*/
		for (iVar = 0; iVar < nPrimVar; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Partial_Gradient = node[iPoint]->GetGradient_Primitive(iSpecies, iVar, iDim) / geometry->node[iPoint]->GetVolume();
				node[iPoint]->SetGradient_Primitive(iSpecies, iVar, iDim, Partial_Gradient);
			}
		}
	}

	delete [] PrimVar_Vertex;
}

void CPlasmaSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
	unsigned short iSpecies, iVar, iDim, jDim, iNeigh;
	unsigned long iPoint, jPoint;
	double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, weight, product;

	unsigned short nPrimVar = nDim+3;

	/*--- Gradient primitive variables compressible (temp, vx, vy, vz, P)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
	PrimVar_i = new double [nPrimVar];
	PrimVar_j = new double [nPrimVar];

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Coord_i = geometry->node[iPoint]->GetCoord();

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

			for (iVar = 0; iVar < nPrimVar; iVar++)
				PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iSpecies, iVar);

			/*--- Inizialization of variables ---*/
			for (iVar = 0; iVar < nPrimVar; iVar++)
				for (iDim = 0; iDim < nDim; iDim++)
					cvector[iVar][iDim] = 0.0;

			r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
			r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

			for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
				jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
				Coord_j = geometry->node[jPoint]->GetCoord();

				for (iVar = 0; iVar < nPrimVar; iVar++)
					PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iSpecies, iVar);

				weight = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

				/*--- Sumations for entries of upper triangular matrix R ---*/
				r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
				r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
				r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
				if (nDim == 3) {
					r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
					r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
					r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
					r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
				}

				/*--- Entries of c:= transpose(A)*b ---*/
				for (iVar = 0; iVar < nPrimVar; iVar++)
					for (iDim = 0; iDim < nDim; iDim++)
						cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/(weight);
			}

			/*--- Entries of upper triangular matrix R ---*/
			r11 = sqrt(r11);
			r12 = r12/(r11);
			r22 = sqrt(r22-r12*r12);
			if (nDim == 3) {
				r13 = r13/(r11);
				r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
				r33 = sqrt(r33-r23*r23-r13*r13);
			}
			/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
			if (nDim == 2) {
				double detR2 = (r11*r22)*(r11*r22);
				Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
				Smatrix[0][1] = -r11*r12/(detR2);
				Smatrix[1][0] = Smatrix[0][1];
				Smatrix[1][1] = r11*r11/(detR2);
			}
			else {
				double detR2 = (r11*r22*r33)*(r11*r22*r33);
				double z11, z12, z13, z22, z23, z33;
				z11 = r22*r33;
				z12 = -r12*r33;
				z13 = r12*r23-r13*r22;
				z22 = r11*r33;
				z23 = -r11*r23;
				z33 = r11*r22;
				Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
				Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
				Smatrix[0][2] = (z13*z33)/(detR2);
				Smatrix[1][0] = Smatrix[0][1];
				Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
				Smatrix[1][2] = (z23*z33)/(detR2);
				Smatrix[2][0] = Smatrix[0][2];
				Smatrix[2][1] = Smatrix[1][2];
				Smatrix[2][2] = (z33*z33)/(detR2);
			}
			/*--- Computation of the gradient: S*c ---*/
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				for (iDim = 0; iDim < nDim; iDim++) {
					product = 0.0;
					for (jDim = 0; jDim < nDim; jDim++)
						product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
					node[iPoint]->SetGradient_Primitive(iSpecies, iVar, iDim, product);
				}
			}
		}
	}

	delete [] PrimVar_i;
	delete [] PrimVar_j;

	Set_MPI_PrimVar_Gradient(geometry, config);
}



void CPlasmaSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config, unsigned long iPoint, double *val_PrimVar_i) {
	unsigned short iSpecies, iVar, iDim, jDim, iNeigh;
	unsigned long jPoint;
	double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, weight, product;

	unsigned short nPrimVar = nDim+3;

	/*--- Gradient primitive variables compressible (temp, vx, vy, vz, P)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
	PrimVar_i = new double [nPrimVar];
	PrimVar_j = new double [nPrimVar];

	/*--- Calculate gradient at iPoint ---*/
	Coord_i = geometry->node[iPoint]->GetCoord();

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

		/*--- Overwrite PrimVar_i[0] with specified wall temperature ---*/
		for (iVar = 0; iVar < nPrimVar; iVar++)
			PrimVar_i[iVar] = val_PrimVar_i[iVar];

		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nPrimVar; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				cvector[iVar][iDim] = 0.0;
			}
		}

		r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();

			for (iVar = 0; iVar < nPrimVar; iVar++)
				PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iSpecies, iVar);

			weight = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

			/*--- Sumations for entries of upper triangular matrix R ---*/
			r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
			r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
			r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
			if (nDim == 3) {
				r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
			}

			/*--- Entries of c:= transpose(A)*b ---*/
			for (iVar = 0; iVar < nPrimVar; iVar++)
				for (iDim = 0; iDim < nDim; iDim++)
					cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/(weight);
		}

		/*--- Entries of upper triangular matrix R ---*/
		r11 = sqrt(r11);
		r12 = r12/(r11);
		r22 = sqrt(r22-r12*r12);
		if (nDim == 3) {
			r13 = r13/(r11);
			r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
			r33 = sqrt(r33-r23*r23-r13*r13);
		}
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		if (nDim == 2) {
			double detR2 = (r11*r22)*(r11*r22);
			Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
			Smatrix[0][1] = -r11*r12/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = r11*r11/(detR2);
		}
		else {
			double detR2 = (r11*r22*r33)*(r11*r22*r33);
			double z11, z12, z13, z22, z23, z33;
			z11 = r22*r33;
			z12 = -r12*r33;
			z13 = r12*r23-r13*r22;
			z22 = r11*r33;
			z23 = -r11*r23;
			z33 = r11*r22;
			Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
			Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
			Smatrix[0][2] = (z13*z33)/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
			Smatrix[1][2] = (z23*z33)/(detR2);
			Smatrix[2][0] = Smatrix[0][2];
			Smatrix[2][1] = Smatrix[1][2];
			Smatrix[2][2] = (z33*z33)/(detR2);
		}
		/*--- Computation of the gradient: S*c ---*/
		for (iVar = 0; iVar < nPrimVar; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				product = 0.0;
				for (jDim = 0; jDim < nDim; jDim++)
					product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
				node[iPoint]->SetGradient_Primitive(iSpecies, iVar, iDim, product);
			}
		}
	}

	delete [] PrimVar_i;
	delete [] PrimVar_j;
}

void CPlasmaSolver::Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iVar_Species, iDim, iSpecies, iMarker, iPeriodic_Index, nPrimVar, nSpeciesPrimVar;
	unsigned long iVertex, iPoint, nVertex, nBuffer_VectorGrad;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi,
	sinPsi, **newGradient = NULL, *Buffer_Receive_VGrad = NULL;
	short SendRecv;
	int send_to, receive_from;

#ifndef NO_MPI
	double *Buffer_Send_VGrad = NULL;
#endif

	nSpeciesPrimVar = nDim+3;
	nPrimVar = nSpecies*nSpeciesPrimVar;

	newGradient = new double* [nPrimVar];
	for (iVar = 0; iVar < nPrimVar; iVar++)
		newGradient[iVar] = new double[3];


	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_VectorGrad	= nVertex*nPrimVar*nDim;

			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;

#ifndef NO_MPI

			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {

				/*--- Allocate viscous variables ---*/
				Buffer_Send_VGrad = new double[nBuffer_VectorGrad];

				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Copy viscous data ---*/
					iVar = 0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						for (iVar_Species = 0; iVar_Species < nSpeciesPrimVar; iVar_Species++) {
							for (iDim = 0; iDim < nDim; iDim++) {
								Buffer_Send_VGrad[iDim*nPrimVar*nVertex+iVar*nVertex+iVertex] = node[iPoint]->GetGradient_Primitive(iSpecies,iVar_Species,iDim);
							}
							iVar++;
						}
					}
				}

				/*--- Send the buffer and deallocate information using MPI ---*/
				MPI::COMM_WORLD.Bsend(Buffer_Send_VGrad, nBuffer_VectorGrad, MPI::DOUBLE, send_to, 9); delete []  Buffer_Send_VGrad;

			}

#endif

			/*--- Receive information  ---*/
			if (SendRecv < 0) {

				/*--- Allocate viscous variables ---*/
				Buffer_Receive_VGrad = new double [nBuffer_VectorGrad];

#ifdef NO_MPI

				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Copy viscous data ---*/
					iVar = 0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						for (iVar_Species = 0; iVar_Species < nSpeciesPrimVar; iVar_Species++) {
							for (iDim = 0; iDim < nDim; iDim++) {
								Buffer_Receive_VGrad[iDim*nPrimVar*nVertex+iVar*nVertex+iVertex] = node[iPoint]->GetGradient_Primitive(iSpecies,iVar_Species,iDim);
							}
							iVar++;
						}
					}
				}

#else
				MPI::COMM_WORLD.Recv(Buffer_Receive_VGrad, nBuffer_VectorGrad, MPI::DOUBLE, receive_from, 9);
#endif

				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();

					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);

					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);

					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;

					for (iVar = 0; iVar < nPrimVar; iVar++)
						for (iDim = 0; iDim < nDim; iDim++)
							newGradient[iVar][iDim] = Buffer_Receive_VGrad[iDim*nPrimVar*nVertex+iVar*nVertex+iVertex];

					/*--- Need to rotate the gradients for all conserved variables. ---*/
					for (iVar = 0; iVar < nPrimVar; iVar ++) {
						if (nDim == 2) {
							newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_VGrad[0*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_VGrad[1*nPrimVar*nVertex+iVar*nVertex+iVertex];
							newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_VGrad[0*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_VGrad[1*nPrimVar*nVertex+iVar*nVertex+iVertex];
						}
						else {
							newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_VGrad[0*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_VGrad[1*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_VGrad[2*nPrimVar*nVertex+iVar*nVertex+iVertex];
							newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_VGrad[0*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_VGrad[1*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_VGrad[2*nPrimVar*nVertex+iVar*nVertex+iVertex];
							newGradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_VGrad[0*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_VGrad[1*nPrimVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_VGrad[2*nPrimVar*nVertex+iVar*nVertex+iVertex];
						}
					}

					/*--- Copy transformed gradients back into buffer. ---*/
					for (iVar = 0; iVar < nPrimVar; iVar++)
						for (iDim = 0; iDim < nDim; iDim++)
							Buffer_Receive_VGrad[iDim*nPrimVar*nVertex+iVar*nVertex+iVertex] = newGradient[iVar][iDim];

					/*--- Viscous method. Store the received information ---*/
					iVar = 0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						for (iVar_Species = 0; iVar_Species < nSpeciesPrimVar; iVar_Species++) {
							for (iDim = 0; iDim < nDim; iDim++)
								node[iPoint]->SetGradient_Primitive(iSpecies,iVar_Species, iDim, Buffer_Receive_VGrad[iDim*nPrimVar*nVertex+iVar*nVertex+iVertex]);
							iVar ++;
						}
					}
				}
				delete [] Buffer_Receive_VGrad;
			}
		}
	}
	for (iVar = 0; iVar < nPrimVar; iVar++)
		delete [] newGradient[iVar];
	delete [] newGradient;

}


void CPlasmaSolver::SetPrimVar_Limiter(CGeometry *geometry, CConfig *config) {

	unsigned long iEdge, iPoint, jPoint, nPoint;
	unsigned short iVar, iDim, iSpecies, nPrimVar;
	double ***Gradient_i, ***Gradient_j, *Coord_i, *Coord_j,
	dave, LimK, eps2, dm, dp, du, limiter;
	double **PrimVar_i, **PrimVar_j;
	// [T-tr, u, v, w, T-vib, P]
	nPrimVar = nDim+3;

	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				PrimVar_max[iPoint][iSpecies][iVar] = -EPS*EPS;
				PrimVar_min[iPoint][iSpecies][iVar] = EPS*EPS;
			}
		}
	}

	/*--- Establish bounds for monotonicity by finding max & min values of neighbor variables --*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		PrimVar_i = node[iPoint]->GetPrimVar_Plasma();
		PrimVar_j = node[jPoint]->GetPrimVar_Plasma();

		/*--- Compute the maximum, and minimum values for nodes i & j ---*/
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				du = PrimVar_j[iSpecies][iVar] - PrimVar_i[iSpecies][iVar];
				PrimVar_min[iPoint][iSpecies][iVar] = min(PrimVar_min[iPoint][iSpecies][iVar], du);
				PrimVar_max[iPoint][iSpecies][iVar] = max(PrimVar_max[iPoint][iSpecies][iVar], du);
				PrimVar_min[jPoint][iSpecies][iVar] = min(PrimVar_min[jPoint][iSpecies][iVar], -du);
				PrimVar_max[jPoint][iSpecies][iVar] = max(PrimVar_max[jPoint][iSpecies][iVar], -du);
			}
		}
	}

	/*--- Initialize the limiter --*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			for (iVar = 0; iVar < nPrimVar; iVar++) {
				node[iPoint]->SetLimiterPrimitive(iSpecies,iVar, 2.0);
			}
		}
	}

	switch (config->GetKind_SlopeLimit()) {

	/*--- Minmod (Roe 1984) limiter ---*/
	case MINMOD:

		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			iPoint     = geometry->edge[iEdge]->GetNode(0);
			jPoint     = geometry->edge[iEdge]->GetNode(1);
			Coord_i    = geometry->node[iPoint]->GetCoord();
			Coord_j    = geometry->node[jPoint]->GetCoord();
			Gradient_i = node[iPoint]->GetGradient_Primitive_Plasma();
			Gradient_j = node[jPoint]->GetGradient_Primitive_Plasma();

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				for (iVar = 0; iVar < nPrimVar; iVar++) {

					/*--- Calculate the interface left gradient, delta- (dm) ---*/
					dm = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iSpecies][iVar][iDim];

					/*--- Calculate the interface right gradient, delta+ (dp) ---*/
					if ( dm > 0.0 ) dp = PrimVar_max[iPoint][iSpecies][iVar];
					else dp = PrimVar_min[iPoint][iSpecies][iVar];

					limiter = max(0.0, min(1.0,dp/dm));

					if (limiter < node[iPoint]->GetLimiterPrimitive(iSpecies, iVar))
						if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiterPrimitive(iSpecies, iVar, limiter);

					/*-- Repeat for point j on the edge ---*/
					dm = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iSpecies][iVar][iDim];

					if ( dm > 0.0 ) dp = PrimVar_max[jPoint][iSpecies][iVar];
					else dp = PrimVar_min[jPoint][iSpecies][iVar];

					limiter = max(0.0, min(1.0,dp/dm));

					if (limiter < node[jPoint]->GetLimiterPrimitive(iSpecies,iVar))
						if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiterPrimitive(iSpecies, iVar, limiter);
				}
			}
		}
		break;

		/*--- Venkatakrishnan (Venkatakrishnan 1994) limiter ---*/
	case VENKATAKRISHNAN:

		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			iPoint     = geometry->edge[iEdge]->GetNode(0);
			jPoint     = geometry->edge[iEdge]->GetNode(1);
			Coord_i    = geometry->node[iPoint]->GetCoord();
			Coord_j    = geometry->node[jPoint]->GetCoord();
			Gradient_i = node[iPoint]->GetGradient_Primitive_Plasma();
			Gradient_j = node[jPoint]->GetGradient_Primitive_Plasma();

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				for (iVar = 0; iVar < nPrimVar; iVar++) {

					/*-- Get limiter parameters from the configuration file ---*/
					dave = config->GetRefElemLength();
					LimK = config->GetLimiterCoeff();
					eps2 = pow((LimK*dave), 3.0);

					/*--- Calculate the interface left gradient, delta- (dm) ---*/
					dm = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iSpecies][iVar][iDim];

					/*--- Calculate the interface right gradient, delta+ (dp) ---*/
					if ( dm > 0.0 ) dp = PrimVar_max[iPoint][iSpecies][iVar];
					else dp = PrimVar_min[iPoint][iSpecies][iVar];

					limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

					if (limiter < node[iPoint]->GetLimiterPrimitive(iSpecies, iVar))
						if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiterPrimitive(iSpecies, iVar, limiter);

					/*-- Repeat for point j on the edge ---*/
					dave = config->GetRefElemLength();
					LimK = config->GetLimiterCoeff();
					eps2 = pow((LimK*dave), 3.0);

					dm = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iSpecies][iVar][iDim];

					if ( dm > 0.0 ) dp = PrimVar_max[jPoint][iSpecies][iVar];
					else dp = PrimVar_min[jPoint][iSpecies][iVar];

					limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

					if (limiter < node[jPoint]->GetLimiterPrimitive(iSpecies, iVar)) {
						if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiterPrimitive(iSpecies, iVar, limiter);
					}
				}
			}
		}
		break;
	}

	SetPrimVar_Limiter_MPI(geometry, config);
}

void CPlasmaSolver::SetPrimVar_Limiter_MPI(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iVar_Species, iSpecies, iMarker, iPeriodic_Index, nPrimVar, nSpeciesPrimVar, loc;
	unsigned long iVertex, iPoint, nVertex, nBuffer_VectorGrad;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi,
	sinPsi, *newLimiter = NULL, *Buffer_Receive_VGrad = NULL;
	short SendRecv;
	int send_to, receive_from;

#ifndef NO_MPI
	double *Buffer_Send_VGrad = NULL;
#endif

	nSpeciesPrimVar = nDim+3;
	nPrimVar = nSpecies*nSpeciesPrimVar;

	newLimiter = new double [nPrimVar];

	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_VectorGrad	= nVertex*nPrimVar;

			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;

#ifndef NO_MPI

			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {

				/*--- Allocate viscous variables ---*/
				Buffer_Send_VGrad = new double[nBuffer_VectorGrad];

				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Copy viscous data ---*/
					iVar = 0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						for (iVar_Species = 0; iVar_Species < nSpeciesPrimVar; iVar_Species++) {

							Buffer_Send_VGrad[iVar*nVertex+iVertex] = node[iPoint]->GetLimiterPrimitive(iSpecies, iVar_Species);
							iVar++;
						}
					}
				}

				/*--- Send the buffer and deallocate information using MPI ---*/
				MPI::COMM_WORLD.Bsend(Buffer_Send_VGrad, nBuffer_VectorGrad, MPI::DOUBLE, send_to, 9);
				delete []  Buffer_Send_VGrad;

			}

#endif

			/*--- Receive information  ---*/
			if (SendRecv < 0) {

				/*--- Allocate viscous variables ---*/
				Buffer_Receive_VGrad = new double [nBuffer_VectorGrad];

#ifdef NO_MPI

				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Copy viscous data ---*/
					iVar = 0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						for (iVar_Species = 0; iVar_Species < nSpeciesPrimVar; iVar_Species++) {

							Buffer_Receive_VGrad[iVar*nVertex+iVertex] = node[iPoint]->GetLimiterPrimitive(iSpecies,iVar_Species);

							iVar++;
						}
					}
				}

#else
				MPI::COMM_WORLD.Recv(Buffer_Receive_VGrad, nBuffer_VectorGrad, MPI::DOUBLE, receive_from, 9);
#endif

				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {

					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();

					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);

					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);

					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;

					for (iVar = 0; iVar < nPrimVar; iVar++)
						newLimiter[iVar] = Buffer_Receive_VGrad[iVar*nVertex+iVertex];

					/*--- Need to rotate the gradients for all conserved variables. ---*/
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						loc = nSpeciesPrimVar*iSpecies;
						if (nDim == 2) {
							newLimiter[loc+1] = rotMatrix[0][0]*Buffer_Receive_VGrad[(loc+1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_VGrad[(loc+2)*nVertex+iVertex];
							newLimiter[loc+2] = rotMatrix[1][0]*Buffer_Receive_VGrad[(loc+1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_VGrad[(loc+2)*nVertex+iVertex];
						}
						else {
							newLimiter[loc+1] = rotMatrix[0][0]*Buffer_Receive_VGrad[(loc+1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_VGrad[(loc+2)*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_VGrad[(loc+3)*nVertex+iVertex];
							newLimiter[loc+2] = rotMatrix[1][0]*Buffer_Receive_VGrad[(loc+1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_VGrad[(loc+2)*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_VGrad[(loc+3)*nVertex+iVertex];
							newLimiter[loc+3] = rotMatrix[2][0]*Buffer_Receive_VGrad[(loc+1)*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_VGrad[(loc+2)*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_VGrad[(loc+3)*nVertex+iVertex];
						}
					}

					/*--- Copy transformed gradients back into buffer. ---*/
					for (iVar = 0; iVar < nPrimVar; iVar++)
						Buffer_Receive_VGrad[iVar*nVertex+iVertex] = newLimiter[iVar];

					/*--- Viscous method. Store the received information ---*/
					iVar = 0;
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						for (iVar_Species = 0; iVar_Species < nSpeciesPrimVar; iVar_Species++) {
							node[iPoint]->SetLimiterPrimitive(iSpecies, iVar_Species, Buffer_Receive_VGrad[iVar*nVertex+iVertex]);
							iVar++;
						}
					}

				}
				delete [] Buffer_Receive_VGrad;
			}
		}
	}
	delete [] newLimiter;
}


void CPlasmaSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
	unsigned short iVar, jVar;

	/*--- Initialise the preconditioning matrix to an identity matrix ---*/
	for (iVar = 0; iVar < nVar; iVar ++) {
		for (jVar = 0; jVar < nVar; jVar ++)
			Precon_Mat_inv[iVar][jVar] = 0.0;
		Precon_Mat_inv[iVar][iVar] = 1.0;
	}
  
#ifdef turkel_acceleration

  unsigned short iDim, jDim, electrons, loc, nVar_species;
  double Beta, local_Mach, Beta2;
	double rho, enthalpy, soundspeed, sq_vel;
	double *U_i;
	double Beta_min = config->GetminTurkelBeta();
  double Beta_max = config->GetmaxTurkelBeta();

	/*--- Variables to calculate the preconditioner parameter Beta ---*/
	electrons = nSpecies - 1;
	loc = (nDim+2)*electrons;
	nVar_species = (nDim+2);

	U_i 			= node[iPoint]->GetSolution();
	rho 			= U_i[loc + 0];
	enthalpy 		= node[iPoint]->GetEnthalpy(electrons);
	soundspeed 		= node[iPoint]->GetSoundSpeed(electrons);
	sq_vel 			= node[iPoint]->GetVelocity2(electrons);
	local_Mach		= sqrt(sq_vel)/node[iPoint]->GetSoundSpeed(electrons);
	Beta 		    = max(Beta_min,min(local_Mach,Beta_max));
	Beta2 		    = Beta*Beta;
	Gamma = config->GetSpecies_Gamma(nSpecies-1);

	/*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
	Precon_Mat_inv[loc + 0][loc + 0] = 0.5*sq_vel;
	Precon_Mat_inv[loc + 0][loc + nVar_species-1] = 1.0;
	for (iDim = 0; iDim < nDim; iDim ++)
		Precon_Mat_inv[loc + 0][loc + 1+iDim] = -1.0*U_i[loc + iDim+1]/rho;

	for (iDim = 0; iDim < nDim; iDim ++) {
		Precon_Mat_inv[loc + iDim+1][loc + 0] = 0.5*sq_vel*U_i[loc + iDim+1]/rho;
		Precon_Mat_inv[loc + iDim+1][loc + nVar_species-1] = U_i[loc + iDim+1]/rho;
		for (jDim = 0; jDim < nDim; jDim ++) {
			Precon_Mat_inv[loc + iDim+1][loc + 1+jDim] = -1.0*U_i[loc + jDim+1]/rho*U_i[loc + iDim+1]/rho;
		}
	}

	Precon_Mat_inv[loc + nVar_species-1][loc + 0] = 0.5*sq_vel*enthalpy;
	Precon_Mat_inv[loc + nVar_species-1][loc + nVar_species-1] = enthalpy;
	for (iDim = 0; iDim < nDim; iDim ++)
		Precon_Mat_inv[loc + nVar_species-1][loc + 1+iDim] = -1.0*U_i[loc + iDim+1]/rho*enthalpy;

	for (iVar = 0; iVar < nVar_species; iVar ++ ) {
		for (jVar = 0; jVar < nVar_species; jVar ++ ) {
			Precon_Mat_inv[loc + iVar][loc + jVar] = (1.0/(Beta2) - 1.0) * (Gamma-1.0)/(soundspeed*soundspeed)*Precon_Mat_inv[loc + iVar][loc + jVar];
			if (iVar == jVar)
				Precon_Mat_inv[loc + iVar][loc + iVar] += 1.0;
		}
	}
#endif
}


void CPlasmaSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iSpecies, iVar, loc;
	double Pressure, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Acquire geometry parameters ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			/*--- Initialize the residual ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual[iVar] = 0.0;
			}

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				Pressure = node[iPoint]->GetPressure(iSpecies);

				/*--- Apply the boundary condition ---*/
				Residual[loc+0] = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[loc + iDim+1] =  Pressure*UnitaryNormal[iDim]*Area;
				Residual[loc+nDim+1] = 0.0;
				if ( iSpecies < nDiatomics )
					Residual[loc+nDim+2] = 0.0;
			}

			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- If needed, determine the Jacobian ---*/
			if (implicit) {
				unsigned short jVar, jDim;
				double dPdrho, dPdE, dPdEv, EPorho, Gamma, Energy_el, Vel2, Density;

				/*--- Initialize the Jacobian ---*/
				for (iVar = 0; iVar < nVar; iVar++)
					for (jVar = 0; jVar < nVar; jVar++)
						Jacobian_i[iVar][jVar] = 0.0;

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

					Density = node[iPoint]->GetSolution(loc+0);
					Vel2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Vel2 += node[iPoint]->GetVelocity(iDim, iSpecies)*node[iPoint]->GetVelocity(iDim, iSpecies);
					Energy_el = 0.0;
					Gamma = config->GetSpecies_Gamma(iSpecies);
					dPdrho = (Gamma - 1.0) * (0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el);
					EPorho = (node[iPoint]->GetSolution(loc+nDim+1) + node[iPoint]->GetPressure(iSpecies))/Density;
					dPdE   = (Gamma - 1.0);
					dPdEv  = -(Gamma - 1.0);


					Jacobian_i[loc+0][loc+0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[loc+0][loc+iDim+1] = UnitaryNormal[iDim]*Area;
					Jacobian_i[loc+0][loc+nDim+1] = 0.0;

					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc+iDim+1][loc+0] = dPdrho * UnitaryNormal[iDim] * Area;
						for (jDim = 0; jDim < nDim; jDim++) {
							Jacobian_i[loc+iDim+1][loc+jDim+1] = (node[iPoint]->GetVelocity(iDim, iSpecies)*UnitaryNormal[jDim] -
									(Gamma-1.0) * node[iPoint]->GetVelocity(jDim, iSpecies)*UnitaryNormal[iDim])*Area;
						}
						Jacobian_i[loc+iDim+1][loc+nDim+1] = dPdE * UnitaryNormal[iDim] * Area;
					}

					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[loc+nDim+1][loc+iDim+1] = EPorho * UnitaryNormal[iDim] * Area;

					if ( iSpecies < nDiatomics) {
						for (iDim = 0; iDim < nDim; iDim++) {
							Jacobian_i[loc+iDim+1][loc+nDim+2] = dPdEv * UnitaryNormal[iDim] * Area;
							Jacobian_i[loc+nDim+2][loc+iDim+1] = node[iPoint]->GetSolution(loc+nDim+2) / Density * UnitaryNormal[iDim] * Area;
						}
					}
				}
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

void CPlasmaSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint,jPoint, total_index, Point_Normal = 0, iNeigh;
	unsigned short iVar, iDim, iSpecies, loc, jVar;
	double *Normal, *Coord_i, *Coord_j, *U_i, *U_j;
	double sq_vel,dist_ij, heat_flux_factor, cp, cpoR, theta, phi, Viscosity, gas_Constant;
	double factor, phi_rho, phi_p, rhoovisc, Twall, Pressure, Density, Temperature, Temperature_Gradient;

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool adiabatic = config->GetAdiabaticWall();
	bool isothermal = config->GetIsothermalWall();
	bool catalytic = config->GetCatalyticWall();

	Twall = 300.0;

	if (adiabatic) {
		for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {

				/*--- Vector --> Velocity_corrected ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;

				/*--- Set the residual, truncation error and velocity value ---*/
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					node[iPoint]->SetVelocity_Old(Vector, iSpecies);
					SetVel_Residual_Zero(iPoint, iSpecies);
					node[iPoint]->SetVel_ResTruncError_Zero(iSpecies);

					/*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal) ---*/
					if (implicit)
						for (iVar = loc + 1; iVar <= loc + nDim; iVar++) {
							total_index = iPoint*nVar+iVar;
							Jacobian.DeleteValsRowi(total_index);
						}
				}
			}
		}
	}
	if (isothermal) {
		for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {

				/*--- Compute the projected residual ---*/
				Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

				double Area = 0.0; double UnitaryNormal[3];
				for (iDim = 0; iDim < nDim; iDim++)
					Area += Normal[iDim]*Normal[iDim];
				Area = sqrt (Area);

				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim] = -Normal[iDim]/Area;
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
					CHeatTransfer[val_marker][iSpecies][iVertex] = 0.0;

				/*--- Compute closest normal neighbor ---*/
				double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;

				double *Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
				cos_max = -1.0;
				for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
					jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
					scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
					for(iDim = 0; iDim < nDim; iDim++) {
						diff_coord = geometry->node[jPoint]->GetCoord(iDim)-geometry->node[iPoint]->GetCoord(iDim);
						scalar_prod += diff_coord*Normal[iDim];
						norm_vect += diff_coord*diff_coord;
						norm_Normal += Normal[iDim]*Normal[iDim];
					}
					norm_vect = sqrt(norm_vect);
					norm_Normal = sqrt(norm_Normal);
					cos_alpha = scalar_prod/(norm_vect*norm_Normal);
					/*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
					if (cos_alpha >= cos_max) {
						Point_Normal = jPoint;
						cos_max = cos_alpha;
					}
				}
				// Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

				/*--- Vector --> Velocity_corrected ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;


				for (iVar = 0; iVar < nVar; iVar ++)
					Res_Visc[iVar] = 0.0;

				/*--- Set the residual, truncation error and velocity value ---*/
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					node[iPoint]->SetVelocity_Old(Vector, iSpecies);
					SetVel_Residual_Zero(iPoint, iSpecies);
					SetVel_Residual_Zero(iPoint, iSpecies);
					SetVel_Residual_Zero(iPoint, iSpecies);
					node[iPoint]->SetVel_ResTruncError_Zero(iSpecies);
				}

				Coord_i = geometry->node[iPoint]->GetCoord();
				Coord_j = geometry->node[Point_Normal]->GetCoord();

				dist_ij = 0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
				dist_ij = sqrt(dist_ij);

				U_i = node[iPoint]->GetSolution() ;
				U_j = node[Point_Normal]->GetSolution();


				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					CHeatTransfer[val_marker][iSpecies][iVertex] = 0.0;
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					Temperature = node[Point_Normal]->GetTemperature_tr(iSpecies);
					Temperature_Gradient = (Twall - Temperature)/dist_ij;
					Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
					gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
					Gamma = config->GetSpecies_Gamma(iSpecies);
					Gamma_Minus_One = Gamma - 1.0;
					cp = (Gamma / Gamma_Minus_One) * gas_Constant;
					heat_flux_factor = cp * Viscosity/PRANDTL;

					Res_Visc[loc + nDim+1] = heat_flux_factor * Temperature_Gradient*Area;

					CHeatTransfer[val_marker][iSpecies][iVertex] -= heat_flux_factor * Temperature_Gradient;
					//                    cout << " CHeatTransfer[iMarker][iSpecies][iVertex] = " <<CHeatTransfer[val_marker][iSpecies][iVertex] << endl;
				}
				LinSysRes.SubtractBlock(iPoint, Res_Visc);  // SIGN CHECK


				/*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal) ---*/
				if (implicit) {

					theta = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						theta += UnitaryNormal[iDim]*UnitaryNormal[iDim];

					for (iVar = 0; iVar < nVar; iVar ++)
						for (jVar = 0; jVar < nVar; jVar ++)
							Jacobian_i[iVar][jVar] = 0.0;

					for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
						sq_vel = 0.0;
						for (iDim = 0; iDim< nDim; iDim ++)
							sq_vel += (U_i[loc + iDim+1]/U_i[loc + 0])*(U_i[loc + iDim+1]/U_i[loc + 0]);

						Density = node[iPoint]->GetDensity(iSpecies);
						Pressure = node[iPoint]->GetPressure(iSpecies);
						Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
						gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
						Gamma = config->GetSpecies_Gamma(iSpecies);
						Gamma_Minus_One = Gamma - 1.0;
						cp =  (Gamma / Gamma_Minus_One)* gas_Constant;
						cpoR = (Gamma / Gamma_Minus_One);

						heat_flux_factor = cp * Viscosity/PRANDTL;
						phi = 0.5*(Gamma-1.0)*sq_vel;
						factor = Viscosity/(Density*dist_ij)*Area;
						phi_rho = -cpoR*heat_flux_factor*Pressure/(Density*Density);
						phi_p = cpoR*heat_flux_factor/Density;
						rhoovisc = Density/(Viscosity); // rho over viscosity


						Jacobian_i[loc + nDim+1][loc + 0] = -factor*rhoovisc*theta*(phi_rho+phi*phi_p);
						Jacobian_i[loc + nDim+1][loc + nDim+1] = -factor*rhoovisc*(Gamma-1)*theta*phi_p;
						for (iVar = loc + 1; iVar <= loc + nDim; iVar++) {
							total_index = iPoint*nVar+iVar;
							Jacobian.DeleteValsRowi(total_index);
						}
					}
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				}
			}
		}
	}

	if (catalytic) {

		double *Mass;
		double ion_flux, electron_flux, atom_flux, ion_density, atom_density, catalyticity_coeff, Kb, ionTemperature, ionmass;
		double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;

		Mass = new double[nSpecies];
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
			Mass[iSpecies] = config->GetParticle_Mass(iSpecies);

		catalyticity_coeff = 1.0;
		ionmass = Mass[1];
		Kb = BOLTZMANN_CONSTANT;


		for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {

				/*--- Compute the projected residual ---*/
				Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

				double Area = 0.0; double UnitaryNormal[3];
				for (iDim = 0; iDim < nDim; iDim++)
					Area += Normal[iDim]*Normal[iDim];
				Area = sqrt (Area);

				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim] = -Normal[iDim]/Area;

				/*--- Compute closest normal neighbor ---*/
				double *Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

				cos_max = -1.0;
				for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
					jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
					scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
					for(iDim = 0; iDim < nDim; iDim++) {
						diff_coord = geometry->node[jPoint]->GetCoord(iDim)-geometry->node[iPoint]->GetCoord(iDim);
						scalar_prod += diff_coord*Normal[iDim];
						norm_vect += diff_coord*diff_coord;
						norm_Normal += Normal[iDim]*Normal[iDim];
					}
					norm_vect = sqrt(norm_vect);
					norm_Normal = sqrt(norm_Normal);
					cos_alpha = scalar_prod/(norm_vect*norm_Normal);
					/*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
					if (cos_alpha >= cos_max) {
						Point_Normal = jPoint;
						cos_max = cos_alpha;
					}
				}
				//	Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

				/*--- Vector --> Velocity_corrected ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;


				for (iVar = 0; iVar < nVar; iVar ++)
					Res_Visc[iVar] = 0.0;

				/*--- Set the residual, truncation error and velocity value ---*/
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					node[iPoint]->SetVelocity_Old(Vector, iSpecies);
					SetVel_Residual_Zero(iPoint, iSpecies);
					SetVel_Residual_Zero(iPoint, iSpecies);
					node[iPoint]->SetVel_ResTruncError_Zero(iSpecies);
				}

				Coord_i = geometry->node[iPoint]->GetCoord();
				Coord_j = geometry->node[Point_Normal]->GetCoord();

				dist_ij = 0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
				dist_ij = sqrt(dist_ij);

				U_i = node[iPoint]->GetSolution() ;
				U_j = node[Point_Normal]->GetSolution();

				ion_density    = U_j[1*(nDim+2) + 0];
				atom_density   = U_j[0*(nDim+2) + 0];
				ionTemperature = node[Point_Normal]->GetTemperature_tr(1);
				ion_flux       = - 0.25 * ion_density * catalyticity_coeff * sqrt( 8*Kb*ionTemperature/ ( PI_NUMBER* ionmass));
				electron_flux  = ion_flux*Mass[2]/Mass[1];
				atom_flux      = - ( ion_flux + electron_flux);

				Res_Visc[0*(nDim+2) + 0] = atom_flux     * Area;
				Res_Visc[1*(nDim+2) + 0] = ion_flux      * Area;
				Res_Visc[2*(nDim+2) + 0] = electron_flux * Area;

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
					CHeatTransfer[val_marker][iSpecies][iVertex] = 0.0;

					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					Temperature = node[Point_Normal]->GetTemperature_tr(iSpecies);
					Temperature_Gradient = (Twall - Temperature)/dist_ij;
					Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
					gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
					Gamma = config->GetSpecies_Gamma(iSpecies);
					Gamma_Minus_One = Gamma - 1.0;
					cp = (Gamma / Gamma_Minus_One) * gas_Constant;
					heat_flux_factor = cp * Viscosity/PRANDTL;
					Res_Visc[loc + nDim+1] = heat_flux_factor * Temperature_Gradient*Area;
					CHeatTransfer[val_marker][iSpecies][iVertex] -= heat_flux_factor * Temperature_Gradient;
				}
				LinSysRes.SubtractBlock(iPoint, Res_Visc);  // SIGN CHECK


				/*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal) ---*/
				if (implicit) {

					theta = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						theta += UnitaryNormal[iDim]*UnitaryNormal[iDim];

					for (iVar = 0; iVar < nVar; iVar ++)
						for (jVar = 0; jVar < nVar; jVar ++)
							Jacobian_i[iVar][jVar] = 0.0;

					Jacobian_i[0*(nDim+2)+0][0*(nDim+2)+0] = atom_flux/atom_density*Mass[1]/Mass[0] * Area;
					Jacobian_i[1*(nDim+2)+0][1*(nDim+2)+0] = ion_flux/ion_density*Mass[1]/Mass[1] * Area;
					Jacobian_i[2*(nDim+2)+0][2*(nDim+2)+0] = electron_flux/ion_density*Mass[1]/Mass[2] * Area;

					for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
						sq_vel = 0.0;
						for (iDim = 0; iDim< nDim; iDim ++)
							sq_vel += (U_i[loc + iDim+1]/U_i[loc + 0])*(U_i[loc + iDim+1]/U_i[loc + 0]);

						Density = node[iPoint]->GetDensity(iSpecies);
						Pressure = node[iPoint]->GetPressure(iSpecies);
						Viscosity = node[iPoint]->GetLaminarViscosity(iSpecies);
						gas_Constant = config->GetSpecies_Gas_Constant(iSpecies);
						Gamma = config->GetSpecies_Gamma(iSpecies);
						Gamma_Minus_One = Gamma - 1.0;
						cp =  (Gamma / Gamma_Minus_One)* gas_Constant;
						cpoR = (Gamma / Gamma_Minus_One);

						heat_flux_factor = cp * Viscosity/PRANDTL;
						phi = 0.5*(Gamma-1.0)*sq_vel;
						factor = Viscosity/(Density*dist_ij)*Area;
						phi_rho = -cpoR*heat_flux_factor*Pressure/(Density*Density);
						phi_p = cpoR*heat_flux_factor/Density;
						rhoovisc = Density/(Viscosity); // rho over viscosity

						Jacobian_i[loc + nDim+1][loc + 0]= -factor*rhoovisc*theta*(phi_rho+phi*phi_p);
						Jacobian_i[loc + nDim+1][loc + nDim+1] = -factor*rhoovisc*(Gamma-1)*theta*phi_p;

						for (iVar = loc + 1; iVar <= loc + nDim; iVar++) {
							total_index = iPoint*nVar+iVar;
							Jacobian.DeleteValsRowi(total_index);
						}
					}
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				}
			}
		}
	}
}

//void CPlasmaSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
//	//Comment: This implementation allows for a specified wall heat flux (typically zero).
//
//	unsigned long iVertex, iPoint, total_index;
//	unsigned short iDim, iVar, iSpecies, loc;
//	double Wall_HeatFlux;
//
//	/*--- Identify the boundary ---*/
//	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
//
//	/*--- Get the specified wall heat flux ---*/
//	Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);
//
//	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//		if (geometry->node[iPoint]->GetDomain()) {
//
//			/*--- Vector --> Velocity_corrected ---*/
//			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
//
//			/*--- Set the residual, truncation error and velocity value ---*/
//			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
//				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
//				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
//				node[iPoint]->SetVelocity_Old(Vector, iSpecies);
//				SetVel_Residual_Zero(iPoint, iSpecies);
//				SetVel_Residual_Zero(iPoint, iSpecies);
//				SetVel_Residual_Zero(iPoint, iSpecies);
//				node[iPoint]->SetVel_ResTruncError_Zero(iSpecies);
//
//				/*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal) ---*/
//				if (implicit)
//					for (iVar = loc + 1; iVar <= loc + nDim; iVar++) {
//						total_index = iPoint*nVar+iVar;
//						Jacobian.DeleteValsRowi(total_index);
//					}
//			}
//		}
//	}
//}

void CPlasmaSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint, Point_Normal;
	unsigned short iSpecies, iVar, jVar, iDim, jDim, loc, nVarSpecies;
	double *Normal, *Coord_i, *Coord_j, Area, dist_ij;
	double UnitaryNormal[3];
	double Twall, Temperature_tr, Temperature_vib, dT_trdrho, dT_vibdrho, dT_vibdev;
	double Density, Pwall, Vel2, Energy_el;
	double Viscosity, ThermalConductivity, ThermalConductivity_vib, GasConstant, h_f, CharVibTemp;
	double Theta;
	double div_vel, **dudx;
	double **ViscFlux_Tensor;
	double ***GradPrimVar;
	unsigned short nSpecies_BC = nSpecies;

	//if (config->GetKind_GasModel() == ARGON) nSpecies_BC = nSpecies -1;


//	if (nDim == 2) {
//		cout << "2D Viscous Jacobian not yet implemented in BC_Isothermal_Wall!!" << endl;
//		cin.get();
//	}

	Point_Normal = 0;
	dudx = new double*[nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		dudx[iDim] = new double[nDim];
	}
	ViscFlux_Tensor = new double*[nDim+3];
	for (iVar = 0; iVar < nDim+3; iVar++)
		ViscFlux_Tensor[iVar] = new double[nDim];

	/*--- Identify the boundary ---*/
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);

	/*--- Retrieve the specified wall temperature ---*/
	Twall = config->GetIsothermal_Temperature(Marker_Tag);

	/*--- Loop over boundary points ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Store dummy vector with no-slip condition ---*/
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;

			/*--- Compute dual-grid area and boundary normal ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = 0.0;
			Theta = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			for (iDim = 0; iDim < nDim; iDim++) {
				UnitaryNormal[iDim] = -Normal[iDim]/Area;
				Theta += UnitaryNormal[iDim]*UnitaryNormal[iDim];
			}

			/*--- Compute closest normal neighbor ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

			/*--- Get cartesian coordinates of i & j and compute inter-node distance ---*/
			Coord_i = geometry->node[iPoint]->GetCoord();
			Coord_j = geometry->node[Point_Normal]->GetCoord();
			dist_ij = 0;
			for (iDim = 0; iDim < nDim; iDim++)
				dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
			dist_ij = sqrt(dist_ij);

			/*--- Load auxiliary vector with no-slip wall velocity ---*/
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;

			/*--- Initialize viscous residual (and Jacobian if implicit) to zero ---*/
			for (iVar = 0; iVar < nVar; iVar ++) {
				Res_Visc[iVar] = 0.0;
			}
			if (implicit) {
				for (iVar = 0; iVar < nVar; iVar ++) {
					for (jVar = 0; jVar < nVar; jVar ++) {
						Jacobian_i[iVar][jVar] = 0.0;
						Jacobian_j[iVar][jVar] = 0.0;
						Jacobian_ElecForce[iVar][jVar] = 0.0;
					}
				}
			}

			/*--- Acquire primitive variable gradients ---*/
			GradPrimVar = node[iPoint]->GetGradient_Primitive_Plasma();

			for (iSpecies = 0; iSpecies < nSpecies_BC; iSpecies++) {
				if ( iSpecies < nDiatomics ) {
					loc = (nDim+3)*iSpecies;
					nVarSpecies = nDim+3;
				}
				else {
					loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					nVarSpecies = nDim+2;
				}

				/*--- Retrieve wall and neighbor quantities ---*/
				Pwall           = node[iPoint]->GetPressure(iSpecies);
				Viscosity       = (node[iPoint]->GetLaminarViscosity(iSpecies) + node[Point_Normal]->GetLaminarViscosity(iSpecies))/2.0;
				Temperature_tr  = node[Point_Normal]->GetPrimVar(iSpecies, 0);
				Temperature_vib = node[Point_Normal]->GetPrimVar(iSpecies, nDim+1);
				ThermalConductivity     = node[iPoint]->GetThermalConductivity(iSpecies);
				ThermalConductivity_vib = node[iPoint]->GetThermalConductivity_vib(iSpecies);

				div_vel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					for (jDim = 0; jDim < nDim; jDim++)
						dudx[iDim][jDim] = GradPrimVar[iSpecies][iDim+1][jDim];
					div_vel += GradPrimVar[iSpecies][iDim+1][iDim];
				}

				/*--- Calculate wall gradients using finite differencing ---*/
				for (iVar = 0; iVar < nDim+3; iVar++)
					for (iDim = 0; iDim < nDim; iDim++)
						ViscFlux_Tensor[iVar][iDim] = 0.0;

				for (iDim = 0; iDim < nDim; iDim++) {
					for (jDim = 0; jDim < nDim; jDim++) {
						ViscFlux_Tensor[iDim+1][jDim] = Viscosity * (dudx[iDim][jDim]+dudx[jDim][iDim]);
					}
					ViscFlux_Tensor[iDim+1][iDim] -= 2.0/3.0 * Viscosity * div_vel;
				}
				/*for (jDim = 0; jDim < nDim; jDim++)
          ViscFlux_Tensor[nDim+1][jDim] = ThermalConductivity * (Twall-Temperature_tr)/dist_ij * UnitaryNormal[jDim];
        if (iSpecies < nDiatomics) {
          for (jDim = 0; jDim < nDim; jDim++) {
            ViscFlux_Tensor[nDim+1][jDim] += ThermalConductivity_vib * (Twall-Temperature_vib)/dist_ij * UnitaryNormal[jDim];
            ViscFlux_Tensor[nDim+2][jDim] = ThermalConductivity_vib * (Twall-Temperature_vib)/dist_ij * UnitaryNormal[jDim];
          }
        }*/

				/*--- Determine viscous contribution to boundary ---*/
				Res_Visc[loc+nDim+1] = ThermalConductivity * (Twall-Temperature_tr)/dist_ij * Area;
				if (iSpecies < nDiatomics) {
					Res_Visc[loc+nDim+1] += ThermalConductivity_vib * (Twall-Temperature_vib)/dist_ij * Area;
					Res_Visc[loc+nDim+2] = ThermalConductivity_vib * (Twall-Temperature_vib)/dist_ij * Area;
				}

				/*--- Apply the no-slip condition to the residual, truncation error and velocity value ---*/
				node[iPoint]->SetVelocity_Old(Vector, iSpecies);
				SetVel_Residual_Zero(iPoint, iSpecies);
				SetVel_Residual_Zero(iPoint, iSpecies);
				SetVel_Residual_Zero(iPoint, iSpecies);
				node[iPoint]->SetVel_ResTruncError_Zero(iSpecies);

				/*--- Calculate Jacobian for implicit time stepping ---*/
				if (implicit) {
					double dPdE, dPdEvib, Energy_tot, Energy_vib, Chi;
          
					/*--- Calculate useful quantities ---*/
					Gamma           = config->GetSpecies_Gamma(iSpecies);
					GasConstant     = config->GetSpecies_Gas_Constant(iSpecies);
					h_f             = config->GetEnthalpy_Formation(iSpecies);
					Density         = node[iPoint]->GetSolution(loc);
					Energy_tot      = node[iPoint]->GetSolution(loc+nDim+1);
					Energy_vib      = 0.0;
					if (iSpecies < nDiatomics)
						Energy_vib    = node[iPoint]->GetSolution(loc+nDim+2);
					Temperature_tr  = Twall;
					Temperature_vib = Twall;
					Vel2            = 0.0;
					Energy_el       = 0.0;
					dPdE            = Gamma - 1.0;
					dPdEvib         = -(Gamma - 1.0);
					dT_trdrho       = 1.0/Density * ( -Temperature_tr + (Gamma-1.0)/GasConstant*(Vel2/2.0 - h_f - Energy_el) );
					Chi             = (Gamma-1.0)/(GasConstant*Density*Density) * (-Energy_tot + Vel2 + Energy_vib);
          if (nDim == 3) {
            /*--- Jacobian of the convective terms (energy rows only) ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ElecForce[loc+nDim+1][loc+iDim+1] = (Energy_tot + Pwall) / Density * UnitaryNormal[iDim] * Area;
            if (iSpecies < nDiatomics) {
              for (iDim = 0; iDim < nDim; iDim++)
                Jacobian_ElecForce[loc+nDim+2][loc+iDim+1] = Energy_vib / Density * UnitaryNormal[iDim] * Area;
            }
            
            /*--- Jacobian of the viscous terms ---*/
            // More terms in loc+4 row follow after the diatomics check
            Jacobian_i[loc+4][loc+0] =  -Chi*ThermalConductivity*Theta/dist_ij * Area;
            Jacobian_i[loc+4][loc+nDim+1] = -(Gamma-1.0)/GasConstant * ThermalConductivity/Density * Theta/dist_ij * Area;
            
            if (iSpecies < nDiatomics) {
              CharVibTemp     = config->GetCharVibTemp(iSpecies);
              dT_vibdrho      =  -Temperature_vib*Temperature_vib/(CharVibTemp*Density)
              * ( 1.0 - 1.0/exp(CharVibTemp/Temperature_vib) );
              dT_vibdev       =  Temperature_vib*Temperature_vib/(CharVibTemp*CharVibTemp*Density*GasConstant)
              * ( (exp(CharVibTemp/Temperature_vib) - 1.0)*(exp(CharVibTemp/Temperature_vib) - 1.0)
                 / exp(CharVibTemp/Temperature_vib) );
              
              Jacobian_i[loc+4][loc+0] -= dT_vibdrho * ThermalConductivity_vib * Theta / dist_ij * Area;
              Jacobian_i[loc+5][loc+0] = -dT_vibdrho * ThermalConductivity_vib * Theta / dist_ij * Area;
              Jacobian_i[loc+0][loc+nDim+2] = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                Jacobian_i[loc+iDim+1][loc+nDim+2] = 0.0;
                Jacobian_i[loc+nDim+2][loc+iDim+1] = 0.0;
              }
              Jacobian_i[loc+nDim+2][loc+nDim+1] = 0.0;
              Jacobian_i[loc+nDim+1][loc+nDim+2] = ( (Gamma-1.0)/GasConstant * ThermalConductivity/Density * Theta/dist_ij * dist_ij - dT_vibdev * ThermalConductivity_vib * Theta/dist_ij ) * Area;
              Jacobian_i[loc+nDim+2][loc+nDim+2] = -dT_vibdev * ThermalConductivity_vib * Theta/dist_ij * Area;
            }
            
            /*--- Determine contribution from normal neighbor to i (store in elecforce Jacobian) ---*/
            for (iVar = 0; iVar < nVarSpecies; iVar++)
              for (jVar = 0; jVar < nVarSpecies; jVar++)
                Jacobian_j[loc+iVar][loc+jVar] = -Jacobian_i[loc+iVar][loc+jVar];
            
            for (iDim = 0; iDim < nDim; iDim++) {
              for (jDim = 0; jDim < nDim; jDim++) {
                Jacobian_i[loc+4][loc+iDim+1] += 1.0/(2.0*Density) * ViscFlux_Tensor[iDim+1][jDim] * UnitaryNormal[jDim] * Area;
                Jacobian_j[loc+4][loc+iDim+1] += 1.0/(2.0*Density) * ViscFlux_Tensor[iDim+1][jDim] * UnitaryNormal[jDim] * Area;
              }
            }
            /*--- Strong enforcement of the no-slip condition in the Jacobian ---*/
            for (iVar = loc+1; iVar < loc+nDim+1; iVar++) {
              unsigned long total_index = iPoint*nVar+iVar;
              Jacobian.DeleteValsRowi(total_index);
            }
            
          }//implicit
        }
			} //iSpecies

			/*--- Apply calculated residuals and Jacobians to the linear system ---*/
			LinSysRes.SubtractBlock(iPoint, Res_Visc);
//			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//			Jacobian.SubtractBlock(iPoint, Point_Normal, Jacobian_j);
//			Jacobian.AddBlock(iPoint, iPoint, Jacobian_ElecForce);
		}
	}

#ifdef Electrons_Invisc
	if (config->GetKind_GasModel() == ARGON) {
		iSpecies = nSpecies -1;

		/*--- Initialize the residual ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			Residual[iVar] = 0.0;
		}
		loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Pressure = node[iPoint]->GetPressure(iSpecies);

		/*--- Apply the boundary condition ---*/
		Residual[loc+0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Residual[loc + iDim+1] =  Pressure*UnitaryNormal[iDim]*Area;
		Residual[loc+nDim+1] = 0.0;


		/*--- Add value to the residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);

		/*--- If needed, determine the Jacobian ---*/
		if (implicit) {
			unsigned short jVar, jDim;
			double dPdrho, dPdE, dPdEv, EPorho, Gamma, Energy_el, Vel2, Density;

			/*--- Initialize the Jacobian ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;

			Density = node[iPoint]->GetSolution(loc+0);
			Vel2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Vel2 += node[iPoint]->GetVelocity(iDim, iSpecies)*node[iPoint]->GetVelocity(iDim, iSpecies);
			Energy_el = 0.0;
			Gamma = config->GetSpecies_Gamma(iSpecies);
			dPdrho = (Gamma - 1.0) * (0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el);
			EPorho = (node[iPoint]->GetSolution(loc+nDim+1) + node[iPoint]->GetPressure(iSpecies))/Density;
			dPdE   = (Gamma - 1.0);
			dPdEv  = -(Gamma - 1.0);


			Jacobian_i[loc+0][loc+0] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Jacobian_i[loc+0][loc+iDim+1] = UnitaryNormal[iDim]*Area;
			Jacobian_i[loc+0][loc+nDim+1] = 0.0;

			for (iDim = 0; iDim < nDim; iDim++) {
				Jacobian_i[loc+iDim+1][loc+0] = dPdrho * UnitaryNormal[iDim] * Area;
				for (jDim = 0; jDim < nDim; jDim++) {
					Jacobian_i[loc+iDim+1][loc+jDim+1] = (node[iPoint]->GetVelocity(iDim, iSpecies)*UnitaryNormal[jDim] -
							(Gamma-1.0) * node[iPoint]->GetVelocity(jDim, iSpecies)*UnitaryNormal[iDim])*Area;
				}
				Jacobian_i[loc+iDim+1][loc+nDim+1] = dPdE * UnitaryNormal[iDim] * Area;
			}

			for (iDim = 0; iDim < nDim; iDim++)
				Jacobian_i[loc+nDim+1][loc+iDim+1] = EPorho * UnitaryNormal[iDim] * Area;

			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
		}
	}
#endif

	for (iDim = 0; iDim < nDim; iDim++) {
		delete [] dudx[iDim];
	}
	delete [] dudx;
	for (iVar = 0; iVar < nDim+3; iVar++)
		delete [] ViscFlux_Tensor[iVar];
	delete [] ViscFlux_Tensor;
}


void CPlasmaSolver::BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

}

void CPlasmaSolver::BC_Dielectric(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

}


/*!
 * \method BC_Inlet
 * \brief Inlet Boundary Condition
 * \author A. Lonkar
 */

void CPlasmaSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iSpecies, loc = 0;
	unsigned short iVar, jVar, iDim;
	double dS, sq_vel,  **P_Matrix, **invP_Matrix, *Face_Normal,
	*kappa, *U_domain, *U_inlet, *U_update, *W_domain, *W_inlet, *W_update;
	double *vn = NULL, *rho = NULL, *rhoE = NULL, **velocity = NULL, *c = NULL, *pressure = NULL, *enthalpy = NULL;

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	kappa = new double[nDim];
	U_domain = new double[nVar]; U_inlet = new double[nVar]; U_update = new double[nVar];
	W_domain = new double[nVar]; W_inlet = new double[nVar]; W_update = new double[nVar];
	P_Matrix = new double* [nVar]; invP_Matrix = new double* [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		P_Matrix[iVar] = new double [nVar];
		invP_Matrix[iVar] = new double [nVar];
	}

	vn		 = new double [nSpecies];
	rho		 = new double [nSpecies];
	rhoE	 = new double [nSpecies];
	c		 = new double [nSpecies];
	velocity = new double*[nSpecies];
	pressure = new double [nSpecies];
	enthalpy = new double [nSpecies];


	/*--- Solution at the inlet ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		U_inlet[loc + 0] = GetDensity_Inlet(iSpecies);
		U_inlet[loc + 1] = GetDensity_Velocity_Inlet(0,iSpecies);
		U_inlet[loc + 2] = GetDensity_Velocity_Inlet(1,iSpecies);
		U_inlet[loc + 3] = GetDensity_Energy_Inlet(iSpecies);
		if (nDim == 3) {
			U_inlet[loc + 3] = GetDensity_Velocity_Inlet(2,iSpecies);
			U_inlet[loc + 4] = GetDensity_Energy_Inlet(iSpecies);
		}
	}

	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			Face_Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			dS = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				dS += Face_Normal[iDim]*Face_Normal[iDim];
			dS = sqrt (dS);

			for (iDim = 0; iDim < nDim; iDim++) {
				kappa[iDim] = -Face_Normal[iDim]/dS;
			}
			/*--- Computation of P and inverse P matrix using values at the infinity ---*/

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				sq_vel = 0.0;
				vn[iSpecies] = 0.0;
				velocity[iSpecies] = new double[nDim];
				loc = iSpecies * (nDim+2);

				rho[iSpecies] = (U_inlet[loc + 0] + U_domain[loc + 0])/2.0;
				rhoE[iSpecies] = (U_inlet[loc + nDim + 1] + U_domain[loc + nDim + 1])/2.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iSpecies][iDim] = (U_inlet[loc + iDim+1]/U_inlet[loc + 0] + U_domain[loc + iDim+1]/U_domain[loc + 0] )/2.0;
					sq_vel += velocity[iSpecies][iDim]*velocity[iSpecies][iDim];
					vn[iSpecies] += velocity[iSpecies][iDim]*kappa[iDim]*dS;
				}

				Gamma = config->GetSpecies_Gamma(iSpecies);
				Gamma_Minus_One = Gamma - 1.0;
				c[iSpecies] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iSpecies]/rho[iSpecies]-0.5*sq_vel)));
			}

			conv_numerics->GetPMatrix_inv(rho, velocity, c, kappa, invP_Matrix);
			conv_numerics->GetPMatrix(rho, velocity, c, kappa, P_Matrix);

			/*--- computation of characteristics variables at the infinity ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_inlet[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_inlet[iVar] +=invP_Matrix[iVar][jVar]*U_inlet[jVar];
			}

			/*--- computation of characteristics variables at the wall ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_domain[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_domain[iVar] +=invP_Matrix[iVar][jVar]*U_domain[jVar];
			}

			/*--- fix characteristics value ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				loc = iSpecies * (nDim+2);

				if (nDim == 2) {
					if(vn[iSpecies] > 0.0) {
						W_update[loc + 0] = W_inlet[loc + 0];
						W_update[loc + 1] = W_inlet[loc + 1];
					}
					else {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
					}

					if(vn[iSpecies]+c[iSpecies]*dS > 0.0) W_update[loc + 2] = W_inlet[loc + 2];
					else W_update[loc + 2] = W_domain[loc + 2];

					if(vn[iSpecies]-c[iSpecies]*dS > 0.0) W_update[loc + 3] = W_inlet[loc + 3];
					else W_update[loc + 3] = W_domain[loc + 3];
				}

				if (nDim == 3) {
					if(vn[iSpecies] > 0.0) {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
						W_update[loc + 2] = W_domain[loc + 2];
					}
					else {
						W_update[loc + 0] = W_inlet[loc + 0];
						W_update[loc + 1] = W_inlet[loc + 1];
						W_update[loc + 2] = W_inlet[loc + 2];
					}

					if(vn[iSpecies] + c[iSpecies]*dS > 0.0) W_update[loc + 3] = W_domain[loc + 3];
					else W_update[loc + 3] = W_inlet[loc + 3];

					if(vn[iSpecies]-c[iSpecies]*dS > 0.0) W_update[loc + 4] = W_domain[loc + 4];
					else W_update[loc + 4] = W_inlet[loc + 4];
				}

			}

			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				U_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					U_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}

			/*--- Residual computation ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->SetConservative(U_domain, U_update);
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);


		}
	}

	for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		delete [] velocity[iSpecies];
	}
	delete [] velocity;
	delete [] vn; delete [] rho; delete [] pressure;
	delete [] rhoE; delete [] c; delete [] enthalpy;
	delete [] kappa;
	delete [] U_domain; delete [] U_inlet; delete [] U_update;
	delete [] W_domain; delete [] W_inlet; delete [] W_update;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Matrix[iVar];
		delete [] invP_Matrix[iVar];
	}
	delete [] P_Matrix;
	delete [] invP_Matrix;
}

/*!
 * \method BC_Outlet
 * \brief Outlet Boundary Condition
 * \author A. Lonkar
 */

void CPlasmaSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint;
	unsigned short iSpecies, loc = 0;
	unsigned short iVar, jVar, iDim;
	double Area, sq_vel,  **P_Matrix, **invP_Matrix, *Normal,
	*UnitaryNormal, *U_domain, *U_outlet, *U_update, *W_domain, *W_outlet, *W_update;

	double *vn = NULL, *rho = NULL, *rhoE = NULL, **velocity = NULL, *c = NULL, *pressure = NULL, *enthalpy = NULL;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	UnitaryNormal = new double[nDim];
	U_domain = new double[nVar]; U_outlet = new double[nVar]; U_update = new double[nVar];
	W_domain = new double[nVar]; W_outlet = new double[nVar]; W_update = new double[nVar];
	P_Matrix = new double* [nVar]; invP_Matrix = new double* [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		P_Matrix[iVar] = new double [nVar];
		invP_Matrix[iVar] = new double [nVar];
	}
	vn		 = new double [nSpecies];
	rho		 = new double [nSpecies];
	rhoE	 = new double [nSpecies];
	c		 = new double [nSpecies];
	velocity = new double*[nSpecies];
	pressure = new double [nSpecies];
	enthalpy = new double [nSpecies];

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the Outlet ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				loc = iSpecies * (nDim+2);

				U_outlet[loc + 0] = GetDensity_Outlet(iSpecies);
				U_outlet[loc + 1] = GetDensity_Velocity_Outlet(0,iSpecies);
				U_outlet[loc + 2] = GetDensity_Velocity_Outlet(1,iSpecies);
				U_outlet[loc + 3] = GetDensity_Energy_Outlet(iSpecies);
				if (nDim == 3) {
					U_outlet[loc + 3] = GetDensity_Velocity_Outlet(2,iSpecies);
					U_outlet[loc + 4] = GetDensity_Energy_Outlet(iSpecies);
				}
			}

			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			/*--- Computation of P and inverse P matrix using values at the domain ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				sq_vel = 0.0;
				vn[iSpecies] = 0.0;
				velocity[iSpecies] = new double[nDim];
				loc = iSpecies * (nDim + 2);

				rho[iSpecies]  = (U_domain[loc + 0] + U_outlet[loc+0]) / 2.0;
				rhoE[iSpecies] = (U_domain[loc + nDim + 1] + U_outlet[loc + nDim + 1])/2.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iSpecies][iDim] = (U_domain[loc + iDim+1]/U_domain[loc + 0] + U_outlet[loc + iDim+1]/U_outlet[loc + 0] )/2.0;
					sq_vel += velocity[iSpecies][iDim]*velocity[iSpecies][iDim];
					vn[iSpecies] += velocity[iSpecies][iDim]*UnitaryNormal[iDim]*Area;
				}
				Gamma = config->GetSpecies_Gamma(iSpecies);
				Gamma_Minus_One = Gamma - 1.0;
				c[iSpecies] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iSpecies]/rho[iSpecies]-0.5*sq_vel)));
			}
			conv_numerics->GetPMatrix_inv(rho, velocity, c, UnitaryNormal, invP_Matrix);
			conv_numerics->GetPMatrix(rho, velocity, c, UnitaryNormal, P_Matrix);

			/*--- computation of characteristics variables at the wall ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_domain[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_domain[iVar] +=invP_Matrix[iVar][jVar]*U_domain[jVar];
			}

			/*--- computation of characteristics variables at the infinity ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_outlet[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_outlet[iVar] +=invP_Matrix[iVar][jVar]*U_outlet[jVar];
			}

			/*--- fix characteristics value ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				loc = iSpecies * (nDim +2);

				if (nDim == 2) {
					if(vn[iSpecies] > 0.0) {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
					}
					else {
						W_update[loc + 0] = W_outlet[loc + 0];
						W_update[loc + 1] = W_outlet[loc + 1];
					}

					if(vn[iSpecies]+c[iSpecies]*Area > 0.0) W_update[loc + 2] = W_domain[loc + 2];
					else W_update[loc + 2] = W_outlet[loc + 2];

					if(vn[iSpecies]-c[iSpecies]*Area > 0.0) W_update[loc + 3] = W_domain[loc + 3];
					else W_update[loc + 3] = W_outlet[loc + 3];

				}

				if (nDim == 3) {
					if(vn[iSpecies] > 0.0) {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
						W_update[loc + 2] = W_domain[loc + 2];
					}
					else {
						W_update[loc + 0] = W_outlet[loc + 0];
						W_update[loc + 1] = W_outlet[loc + 1];
						W_update[loc + 2] = W_outlet[loc + 2];
					}

					if(vn[iSpecies]+c[iSpecies]*Area > 0.0) W_update[loc + 3] = W_domain[loc + 3];
					else W_update[loc + 3] = W_outlet[loc + 3];


					if(vn[iSpecies]-c[iSpecies]*Area > 0.0) W_update[loc + 4] = W_domain[loc + 4];
					else W_update[loc + 4] = W_outlet[loc + 4];
				}
			}

			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				U_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					U_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}

			/*--- Residual computation compressible flow ---*/
			sq_vel = 0.0;
			for (iSpecies = 0; iSpecies <nSpecies; iSpecies ++) {
				loc = iSpecies * (nDim +2);
				rho[iSpecies] = U_update[loc + 0];
				rhoE[iSpecies] = U_update[loc + nDim +1];
				sq_vel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iSpecies][iDim] = U_update[loc + iDim+1]/rho[iSpecies];
					sq_vel +=velocity[iSpecies][iDim]*velocity[iSpecies][iDim];
				}
				Gamma = config->GetSpecies_Gamma(iSpecies);
				Gamma_Minus_One = Gamma - 1.0;
				c[iSpecies] = sqrt(Gamma*Gamma_Minus_One*(rhoE[iSpecies]/rho[iSpecies] - 0.5*sq_vel));
				pressure[iSpecies] = (c[iSpecies] * c[iSpecies] * rho[iSpecies]) / Gamma;
				enthalpy[iSpecies] = (rhoE[iSpecies] + pressure[iSpecies]) / rho[iSpecies];
			}

			conv_numerics->GetInviscidProjFlux(rho, velocity, pressure, enthalpy, UnitaryNormal, Residual);

			for	(iVar = 0; iVar < nVar; iVar++)
				Residual[iVar] = Residual[iVar]*Area;

			LinSysRes.AddBlock(iPoint, Residual);
			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
				conv_numerics->SetNormal(Vector);
				conv_numerics->SetConservative(U_domain, U_outlet);
				conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

			}
		}
	}

	for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
		delete [] velocity[iSpecies];
	}
	delete [] velocity;
	delete [] vn; delete [] rho; delete [] pressure;
	delete [] rhoE; delete [] c; delete [] enthalpy;
	delete [] UnitaryNormal;
	delete [] U_domain; delete [] U_outlet; delete [] U_update;
	delete [] W_domain; delete [] W_outlet; delete [] W_update;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Matrix[iVar];
		delete [] invP_Matrix[iVar];
	}
	delete [] P_Matrix;
	delete [] invP_Matrix;
}

/*!
 * \method BC_Neumann
 * \brief Neumann Boundary Condition
 * \author A. Lonkar
 */
void CPlasmaSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint, Point_Normal, total_index;
	unsigned short iVar;
	double *U_domain, *U_interior;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	U_domain = new double[nVar];
	U_interior = new double[nVar];

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			for (iVar = 0; iVar < nVar; iVar++) {
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);
				U_interior[iVar] = node[Point_Normal]->GetSolution(iVar);
			}
			node[iPoint]->SetSolution_Old(U_interior);
			node[iPoint]->SetSolution(U_interior);
			LinSysRes.SetBlock_Zero(iPoint);
			node[iPoint]->SetRes_TruncErrorZero();

			/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
			if (implicit)
				for (iVar = 0; iVar < nVar; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
		}
	}
	delete [] U_domain; delete [] U_interior;

}
/*!
 * \method BC_Far_Field
 * \brief Far field boundary condition
 * \author A. Lonkar
 */
void CPlasmaSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, Point_Normal, iSpecies, loc;
	unsigned short iVar, iDim;
	double *U_domain, *U_infty, **V_infty;
  double Gas_constant, Gamma, Vel2, energy_vib, energy_el, Char_temp_vib, Molar_mass;
  
  bool viscous = (config->GetKind_Solver() == PLASMA_NAVIER_STOKES);
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	U_domain = new double[nVar];
	U_infty = new double[nVar];
  
  V_infty = new double *[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    V_infty[iSpecies] = new double[nDim+6];

	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the infinity ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        
        Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
        Gamma = config->GetSpecies_Gamma(iSpecies);
        Char_temp_vib = config->GetCharVibTemp(iSpecies);
        Molar_mass = config->GetMolar_Mass(iSpecies);
        Vel2 = 0.0;
        energy_vib = 0.0;
        energy_el = 0.0;
        if (iSpecies < nDiatomics)
          energy_vib = GetDensity_Energy_vib_Inf(iSpecies);

				U_infty[loc + 0] = GetDensity_Inf(iSpecies);
        V_infty[iSpecies][nDim+3] = GetDensity_Inf(iSpecies);
				for (iDim = 0; iDim < nDim; iDim++) {
					U_infty[loc+iDim+1] = GetDensity_Velocity_Inf(iDim, iSpecies);
          V_infty[iSpecies][iDim+1] = U_infty[loc+iDim+1]/V_infty[iSpecies][nDim+3];
          Vel2 += V_infty[iSpecies][iDim+1];
        }
				U_infty[loc+nDim+1] = GetDensity_Energy_Inf(iSpecies);
        V_infty[iSpecies][0] = (Gamma-1.0)/Gas_constant * (U_infty[loc+nDim+1] - 0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies) - energy_vib - energy_el);
				if (iSpecies < nDiatomics)
					U_infty[loc+nDim+2] = energy_vib;
        V_infty[iSpecies][nDim+1] = Char_temp_vib / log(Char_temp_vib*UNIVERSAL_GAS_CONSTANT/(energy_vib*Molar_mass) + 1.0);
			}

			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->SetConservative(U_domain, U_infty);
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      if (viscous) {
        
        /*--- Index of the closest interior node ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
        
        /*--- Points, coordinates and normal vector in edge ---*/
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Vector);
        
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (iVar = 0 ; iVar < nPrimVar; iVar++) {
            for (iDim = 0; iDim < nDim; iDim++) {
              PrimGrad_i[iSpecies][iVar][iDim] = node[iPoint]->GetGradient_Primitive(iSpecies, iVar, iDim);
            }
          }
        }
        
        /*--- Acquire primitive variables and gradients from the CPlasmaVariable class ---*/
        visc_numerics->SetPrimitive(node[iPoint]->GetPrimVar_Plasma(), V_infty);
        visc_numerics->SetPrimVarGradient(PrimGrad_i, PrimGrad_i);
        
        /*--- Pass transport coefficients from the variable to the numerics class ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
          visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(iSpecies), node[iPoint]->GetLaminarViscosity(iSpecies), iSpecies);
          visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(iSpecies), node[iPoint]->GetEddyViscosity(iSpecies), iSpecies);
          visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(iSpecies), node[iPoint]->GetThermalConductivity(iSpecies), iSpecies);
          visc_numerics->SetThermalConductivity_vib(node[iPoint]->GetThermalConductivity_vib(iSpecies), node[iPoint]->GetThermalConductivity_vib(iSpecies), iSpecies);
        }
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Res_Visc);
        
        /*--- Update Jacobians for implicit calculations ---*/
        if (implicit) {
          Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
        }
      }
		}
	}
	delete [] U_domain;
	delete [] U_infty;
}

/*!
 * \method BC_Sym_Plane
 * \brief Symmetry Boundary Condition
 * \author A. Lonkar.  Modified by S. R. Copeland
 */

void CPlasmaSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iSpecies, iVar, loc;
	double Pressure, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Acquire geometry parameters ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			/*--- Initialize the residual ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual[iVar] = 0.0;
			}

			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				Pressure = node[iPoint]->GetPressure(iSpecies);

				/*--- Apply the boundary condition ---*/
				Residual[loc+0] = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[loc + iDim+1] =  Pressure*UnitaryNormal[iDim]*Area;
				Residual[loc+nDim+1] = 0.0;
				if ( iSpecies < nDiatomics )
					Residual[loc+nDim+2] = 0.0;
			}

			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- If needed, determine the Jacobian ---*/
			if (implicit) {
				unsigned short jVar, jDim;
				double dPdrho, dPdE, dPdEv, EPorho, Gamma, Energy_el, Vel2, Density;

				/*--- Initialize the Jacobian ---*/
				for (iVar = 0; iVar < nVar; iVar++)
					for (jVar = 0; jVar < nVar; jVar++)
						Jacobian_i[iVar][jVar] = 0.0;

				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

					Density = node[iPoint]->GetSolution(loc+0);
					Vel2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Vel2 += node[iPoint]->GetVelocity(iDim, iSpecies)*node[iPoint]->GetVelocity(iDim, iSpecies);
					Energy_el = 0.0;
					Gamma = config->GetSpecies_Gamma(iSpecies);
					dPdrho = (Gamma - 1.0) * (0.5*Vel2 - config->GetEnthalpy_Formation(iSpecies) - Energy_el);
					EPorho = (node[iPoint]->GetSolution(loc+nDim+1) + node[iPoint]->GetPressure(iSpecies))/Density;
					dPdE   = (Gamma - 1.0);
					dPdEv  = -(Gamma - 1.0);


					Jacobian_i[loc+0][loc+0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[loc+0][loc+iDim+1] = UnitaryNormal[iDim]*Area;
					Jacobian_i[loc+0][loc+nDim+1] = 0.0;

					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc+iDim+1][loc+0] = dPdrho * UnitaryNormal[iDim] * Area;
						for (jDim = 0; jDim < nDim; jDim++) {
							Jacobian_i[loc+iDim+1][loc+jDim+1] = (node[iPoint]->GetVelocity(iDim, iSpecies)*UnitaryNormal[jDim] -
									(Gamma-1.0) * node[iPoint]->GetVelocity(jDim, iSpecies)*UnitaryNormal[iDim])*Area;
						}
						Jacobian_i[loc+iDim+1][loc+nDim+1] = dPdE * UnitaryNormal[iDim] * Area;
					}

					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[loc+nDim+1][loc+iDim+1] = EPorho * UnitaryNormal[iDim] * Area;

					if ( iSpecies < nDiatomics) {
						for (iDim = 0; iDim < nDim; iDim++) {
							Jacobian_i[loc+iDim+1][loc+nDim+2] = dPdEv * UnitaryNormal[iDim] * Area;
							Jacobian_i[loc+nDim+2][loc+iDim+1] = node[iPoint]->GetSolution(loc+nDim+2) / Density * UnitaryNormal[iDim] * Area;
						}
					}
				}
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}
